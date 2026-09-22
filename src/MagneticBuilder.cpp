#include "mvb/MagneticBuilder.h"
#include <cstdio>
#include <cstdlib>
#include "mvb/StepExporter.h"
#include "mvb/Symmetry.h"
#include "mvb/Utils.h"
#include "mvb/shapes/ShapeBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/ConductorBuilder.h"
#include "mvb/BobbinBuilder.h"
#include "mvb/SpacerBuilder.h"   // ABT #1170 (WP1)
#include "mvb/PinBuilder.h"      // ABT #1171 (WP2)
#include "mvb/ShuntBuilder.h"    // ABT #1176 (WP7)
#include "mvb/BaseBuilder.h"     // ABT #1173 (WP4)
#include "mvb/FR4Builder.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/CorePiece.h"
#include "support/Utils.h"
#include <nlohmann/json.hpp>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepOffset_Mode.hxx>
#include <GeomAbs_JoinType.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <Standard_Failure.hxx>
#include <iostream>
#include <BRepBndLib.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <Bnd_Box.hxx>
#include "support/Settings.h"
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <numbers>

namespace mvb {

using json = nlohmann::json;

// The coating thickness to DRAW, resolved from the MAS the caller handed in.
//
// Only a coating the design actually declares is drawn. MKF's Core::get_coating_thickness()
// also invents a default jacket for an uncoated toroid -- correct for the dielectric path it
// exists for (a bare-ferrite toroid wound with bare wire is rare), but drawing that would put
// a shell on every toroid in every consumer, including designs whose author never mentioned a
// coating. So the PRESENCE of the layer is decided here, from functionalDescription.coating,
// while its THICKNESS is still resolved by MKF -- which is what turns the name-only form
// ("epoxy") into a datasheet number, and what nothing here should be re-deriving.
//
// MAGNETIC_EPOXY is not this kind of coating: it is the powder-loaded shield cap moulded over
// the winding of a semishielded drum, drawn separately (and translucently) by drawCoreShell.
// Same exclusion StrayCapacitance::resolve_core_jacket makes, for the same reason.
double MagneticBuilder::declaredCoreCoatingThickness(const MAS::MagneticCore& masCore) {
    const auto coating = masCore.get_functional_description().get_coating();
    if (!coating) return 0.0;
    if (std::holds_alternative<MAS::CoreCoating>(coating.value())) {
        const auto type = std::get<MAS::CoreCoating>(coating.value()).get_type();
        if (type && type.value() == MAS::CoatingType::MAGNETIC_EPOXY) return 0.0;
    }
    return OpenMagnetics::Core(masCore).get_coating_thickness();
}

// Build the core's insulating COATING (epoxy/parylene/etc., MAS CoreCoating) as a real conformal
// SHELL solid: offset the core surface outward by the coating thickness and subtract the core, so
// the result is a uniform-thickness layer wrapping the whole core (outer/inner/top/bottom). Returns
// a null shape if the offset/cut fails (e.g. a sharp-cornered core OCCT cannot offset) -- caller skips
// it (no fabrication). Geometry-agnostic: works on any core solid (toroid, PQ, ...).
TopoDS_Shape MagneticBuilder::buildCoreCoatingShell(const TopoDS_Shape& core, double thickness) {
    if (thickness <= 0.0 || core.IsNull()) return TopoDS_Shape();
    try {
        BRepOffsetAPI_MakeOffsetShape mk;
        mk.PerformByJoin(core, thickness, 1e-7, BRepOffset_Skin,
                         Standard_False, Standard_False, GeomAbs_Intersection);
        if (!mk.IsDone()) return TopoDS_Shape();
        TopoDS_Shape bigger = mk.Shape();
        if (bigger.IsNull()) return TopoDS_Shape();
        BRepAlgoAPI_Cut cut(bigger, core);            // shell = offset - core
        if (!cut.IsDone()) return TopoDS_Shape();
        TopoDS_Shape shell = cut.Shape();
        return shell.IsNull() ? TopoDS_Shape() : shell;
    } catch (const Standard_Failure&) {
        return TopoDS_Shape();
    }
}

static bool isCoreToroidal(const MAS::MagneticCore& core) {
    auto geo = core.get_geometrical_description();
    if (!geo) return false;
    for (const auto& piece : *geo) {
        if (piece.get_type() == MAS::CoreGeometricalDescriptionElementType::TOROIDAL) {
            return true;
        }
    }
    return false;
}

template<typename BobbinT, typename VariantT>
static MAS::CoreBobbinProcessedDescription getBobbinProcessedT(const VariantT& bobbinVar) {
    const BobbinT* bobbin = std::get_if<BobbinT>(&bobbinVar);
    if (bobbin) {
        auto pd = bobbin->get_processed_description();
        if (pd) return *pd;
    }
    return MAS::CoreBobbinProcessedDescription();
}

// ---- ABT #1248: toroid mounting -------------------------------------------------------------
template<typename BobbinT, typename VariantT>
static ConductorBuilder::ToroidMounting toroidMountingT(const VariantT& bobbinVar) {
    // A toroid BASE record states how the ring sits on it; that overrides the global setting.
    if (const BobbinT* bobbin = std::get_if<BobbinT>(&bobbinVar)) {
        const auto& fd = bobbin->get_functional_description();
        if (fd && fd->get_family() == MAS::BobbinFamily::T && fd->get_base()) {
            return fd->get_base()->get_mounting() == MAS::OrientationEnum::HORIZONTAL
                       ? ConductorBuilder::ToroidMounting::Horizontal
                       : ConductorBuilder::ToroidMounting::Vertical;
        }
    }
    return OpenMagnetics::Settings::GetInstance().get_toroid_mounting() == MAS::OrientationEnum::HORIZONTAL
               ? ConductorBuilder::ToroidMounting::Horizontal
               : ConductorBuilder::ToroidMounting::Vertical;
}

static bool isCoreToroidal(const MAS::MagneticCore& core);

// Depth of the core (with its coating shells) along the build-frame lead direction: rotate the
// solids so `down` is -Y and read the tight bounding box (AddOptimal: a ring's rim is exact).
static double coreDepthAlongDown(const std::vector<NamedShape>& coreShapes,
                                 const std::array<double, 3>& down) {
    const gp_Vec from(down[0], down[1], down[2]);
    const gp_Vec to(0.0, -1.0, 0.0);
    gp_Trsf align;
    const gp_Vec axis = from.Crossed(to);
    if (axis.Magnitude() > 1e-12) {
        align.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(axis)), from.AngleWithRef(to, axis));
    }
    else if (from.Dot(to) < 0.0) {
        align.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), std::numbers::pi);
    }
    Bnd_Box box;
    bool any = false;
    for (const auto& ns : coreShapes) {
        if (ns.role != Role::Core && ns.role != Role::CoreCoating) continue;
        if (ns.shape.IsNull()) continue;
        BRepBndLib::AddOptimal(BRepBuilderAPI_Transform(ns.shape, align, true).Shape(), box,
                               /*useTriangulation=*/false, /*useShapeTolerance=*/false);
        any = true;
    }
    if (!any || box.IsVoid())
        throw std::runtime_error("coreDepthAlongDown: the toroid has no core solid to measure");
    double xmin, ymin, zmin, xmax, ymax, zmax;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    return -ymin;
}

static MAS::CoreBobbinProcessedDescription getBobbinProcessed(const MAS::Coil& coil) {
    return getBobbinProcessedT<MAS::Bobbin>(coil.get_bobbin());
}

static MAS::CoreBobbinProcessedDescription getBobbinProcessed(const OpenMagnetics::Coil& coil) {
    return getBobbinProcessedT<OpenMagnetics::Bobbin>(coil.get_bobbin());
}

template<typename BobbinT, typename VariantT>
static std::string getBobbinNameT(const VariantT& bobbinVar, const std::string& fallback) {
    if (const BobbinT* b = std::get_if<BobbinT>(&bobbinVar)) {
        return b->get_name().value_or(fallback);
    }
    if (const std::string* s = std::get_if<std::string>(&bobbinVar)) {
        return s->empty() ? fallback : *s;
    }
    return fallback;
}

static void patchBobbinDimensions(MAS::CoreBobbinProcessedDescription& bobbinPd, const MAS::MagneticCore& core) {
    double colWidth = bobbinPd.get_column_width().value_or(0.0);

    // MKF sometimes leaves column_depth uninitialized and column_width empty.
    // Fall back to the core's central column dimensions, using MKF's convention:
    // column_width = core_column_width / 2 + wall_thickness
    // (this is the outer radius of the bobbin body tube).
    //
    // No silent default for wall_thickness: if MKF/MAS leaves it unset (NaN
    // or negative), we honour that as "no wall" and let the caller decide
    // — per the no-fallback policy, fabricating a 0.5 mm wall here would
    // hide whatever the upstream actually wanted.
    if (colWidth <= 0.0) {
        auto corePd = core.get_processed_description();
        if (corePd && !corePd->get_columns().empty()) {
            const auto& centralCol = corePd->get_columns()[0];
            double wallThickness = bobbinPd.get_wall_thickness();
            if (std::isnan(wallThickness) || wallThickness < 0.0) {
                wallThickness = 0.0;
            }
            if (centralCol.get_width() > 0.0) {
                bobbinPd.set_column_width(centralCol.get_width() / 2.0 + wallThickness);
            }
            if (centralCol.get_depth() > 0.0) {
                bobbinPd.set_column_depth(centralCol.get_depth() / 2.0 + wallThickness);
            }
            // Propagate the column shape too — FR4Builder picks its annulus /
            // rectangular / oblong branch off this. The default-constructed
            // ColumnShape is RECTANGULAR, which silently produces the wrong
            // board for round-column cores.
            bobbinPd.set_column_shape(centralCol.get_shape());
        }
    }
}

static std::string drawMagneticCommon(const std::vector<NamedShape>& named,
                                      const std::string& outputPath,
                                      const std::string& format,
                                      double scale,
                                      bool femReady) {
    // ABT #685: keep the NamedShapes WHOLE. Splitting them into parallel {shape, name}
    // vectors and calling the shapes+names overload silently dropped `partNames`, so every
    // STEP written through drawMagnetic (i.e. everything the CLI and the bindings produce)
    // came out with one unnamed multi-solid product per conductor and the viewer invented
    // its own numbering for the pieces. Only the tests, which call exportSTEP(NamedShape)
    // directly, ever saw the per-solid names.
    std::vector<NamedShape> scaled = named;
    std::vector<TopoDS_Shape> allShapes;
    allShapes.reserve(scaled.size());

    if (scale != 1.0) {
        gp_Trsf trsf;
        trsf.SetScale(gp_Pnt(0, 0, 0), scale);
        for (auto& ns : scaled)
            ns.shape = BRepBuilderAPI_Transform(ns.shape, trsf).Shape();
    }
    for (const auto& ns : scaled) allShapes.push_back(ns.shape);

    std::filesystem::path out = outputPath;
    if (format == "stl") {
        out /= "magnetic.stl";
        // Fuse all shapes into a single compound so exportSTL exports the
        // full assembly, not just the first shape.
        TopoDS_Compound compound;
        TopoDS_Builder b;
        b.MakeCompound(compound);
        bool any = false;
        for (const auto& s : allShapes) {
            if (s.IsNull()) continue;
            b.Add(compound, s);
            any = true;
        }
        exportSTL(any ? compound : TopoDS_Shape(), out.string());
        return out.string();
    }
    out /= "magnetic.step";
    StepExportOptions opts;
    opts.nurbsPeriodicSolids = femReady;
    exportSTEP(scaled, out.string(), opts);
    return out.string();
}

std::string MagneticBuilder::drawMagnetic(const MAS::Magnetic& magnetic,
                                           const std::string& outputPath,
                                           const std::string& format,
                                           bool includeBobbin,
                                           double scale,
                                           int symmetryPlanes,
                                           int wirePolygonSegments,
                                           int corePolygonSegments) const {
    DrawConfig cfg{format, includeBobbin, scale, symmetryPlanes,
                   wirePolygonSegments, corePolygonSegments};
    return drawMagnetic(magnetic, outputPath, cfg);
}

std::string MagneticBuilder::drawMagnetic(const MAS::Magnetic& magnetic,
                                           const std::string& outputPath,
                                           const DrawConfig& cfg) const {
    auto named = buildAllNamed(magnetic, cfg.includeBobbin, cfg.symmetryPlanes,
                               cfg.wirePolygonSegments, cfg.corePolygonSegments,
                               cfg.paintCoating, /*emitCoatingShells=*/false,
                               /*includeInsulation=*/false,
                               magnetic.get_core() ? declaredCoreCoatingThickness(magnetic.get_core().value()) : 0.0,
                               cfg.useRealWindingGeometry, cfg.femReady);
    return drawMagneticCommon(named, outputPath, cfg.format, cfg.scale, cfg.femReady);
}

std::string MagneticBuilder::drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                                           const std::string& outputPath,
                                           const std::string& format,
                                           bool includeBobbin,
                                           double scale,
                                           int symmetryPlanes,
                                           int wirePolygonSegments,
                                           int corePolygonSegments) const {
    DrawConfig cfg{format, includeBobbin, scale, symmetryPlanes,
                   wirePolygonSegments, corePolygonSegments};
    return drawMagnetic(magnetic, outputPath, cfg);
}

std::string MagneticBuilder::drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                                           const std::string& outputPath,
                                           const DrawConfig& cfg) const {
    auto named = buildAllNamed(magnetic, cfg.includeBobbin, cfg.symmetryPlanes,
                               cfg.wirePolygonSegments, cfg.corePolygonSegments,
                               cfg.paintCoating, /*emitCoatingShells=*/false,
                               /*includeInsulation=*/false,
                               declaredCoreCoatingThickness(magnetic.get_core()),
                               cfg.useRealWindingGeometry, cfg.femReady);
    return drawMagneticCommon(named, outputPath, cfg.format, cfg.scale, cfg.femReady);
}

namespace {

// Internal core-piece builder. Kept as a free helper so the only public API
// for emitting core shapes is buildCoreNamed(). Callers that don't need the
// piece names use this directly; callers that need names go through
// MagneticBuilder::buildCoreNamed.
std::vector<TopoDS_Shape> buildCoreShapes_impl(const MAS::MagneticCore& core,
                                               int corePolygonSegments) {
    std::vector<TopoDS_Shape> result;
    auto geoOpt = core.get_geometrical_description();
    if (!geoOpt) return result;

    for (const auto& piece : *geoOpt) {
        // CLOSED covers single-solid cores (UT, and DRUM per ABT #331): their family
        // builders emit the complete piece, so they go through the same path as
        // HALF_SET. PLATE (piece-and-plate closers, ABT #264) still needs a
        // plate-dimension build and is not handled yet.
        if (piece.get_type() != MAS::CoreGeometricalDescriptionElementType::HALF_SET
            && piece.get_type() != MAS::CoreGeometricalDescriptionElementType::TOROIDAL
            && piece.get_type() != MAS::CoreGeometricalDescriptionElementType::CLOSED) {
            continue;
        }

        auto shapeOpt = piece.get_shape();
        if (!shapeOpt) continue;
        const MAS::CoreShape* shapeData = std::get_if<MAS::CoreShape>(&*shapeOpt);
        if (!shapeData) continue;

        auto builder = shapes::createShapeBuilder(shapeData->get_family(), "", corePolygonSegments);
        if (!builder) continue;

        TopoDS_Shape shape = builder->buildPiece(*shapeData);
        if (shape.IsNull()) continue;

        auto dimsOpt = shapeData->get_dimensions();
        auto dims = dimsOpt ? flatten_dimensions(*dimsOpt) : std::map<std::string, double>{};

        // Apply rotation from geometrical description FIRST; Python MVB applies
        // machining in the post-rotation frame so the gap coordinates refer to
        // the already-flipped piece. Applying machining before rotation leaves
        // the gap tool positioned outside the flipped piece and the cut is a
        // silent no-op (observed on stacked E cores — example 18).
        auto rotOpt = piece.get_rotation();
        if (rotOpt && rotOpt->size() >= 3) {
            shape = rotate_shape(shape, (*rotOpt)[0], (*rotOpt)[1], (*rotOpt)[2]);
        }

        // Apply machining after rotation. EVERY cut must remove material. This used to say the
        // un-gapped shape was "a strictly-better fallback than no shape" -- it is the opposite:
        // a gap that silently is not there makes the FEM report the gapless inductance with no
        // sign anything is wrong (00_debug, 2026-09-05: +82% against MKF, chased for hours
        // through the mesher and the solver before the core turned out to be uncut). A boolean
        // that fails throws inside applyMachining; a boolean that "succeeds" without touching
        // the piece (tool placed off the column) is caught here by the volume.
        //
        // "Removed nothing" is measured against the piece as the knife FOUND it, and that answers
        // two different questions at once -- only one of which is a bug (ABT #1184):
        //   * the tool never overlapped the column, so the gap is missing from the geometry. That
        //     is the failure this guard exists for and it still throws;
        //   * the tool did overlap the column, but an EARLIER gap on the same column had already
        //     ground that material away, so the gap IS in the geometry, cut by its overlapping
        //     neighbour. MKF spaces a distributed gap by columnHeight / (numberGaps + 1) counting
        //     gap CENTRES only, never their lengths, so gaps as long as that spacing overlap and
        //     one can land wholly inside another: on RM 5/8 (3.8 mm central column) the 0.5 mm gap
        //     at y = 0 and the 2 mm gap at y = +0.95 mm put the latter's 0.05 mm bottom-half slab
        //     (y in [-0.05, 0]) entirely inside the former's (y in [-0.25, 0]).
        // So a cut that removes nothing is retried against the piece as it was before ANY gap, and
        // only a tool that removes nothing from THAT throws: a tool placed off the column finds no
        // material on the pristine piece either, so the original defect is still caught.
        const TopoDS_Shape pristinePiece = shape;
        auto machiningOpt = piece.get_machining();
        if (machiningOpt) {
            for (const auto& mach : *machiningOpt) {
                GProp_GProps before, after;
                BRepGProp::VolumeProperties(shape, before);
                shape = builder->applyMachining(shape, mach, dims);
                BRepGProp::VolumeProperties(shape, after);
                const double removed = before.Mass() - after.Mass();
                const auto& mc = mach.get_coordinates();
                if (std::getenv("MVB_CORE_DIAG"))
                    std::fprintf(stderr, "[machining] '%s' gap %.4g mm at (%.4g, %.4g) mm: removed %.5f mm3 (piece %.4f -> %.4f mm3)\n",
                                 shapeData->get_name().value_or(std::string("?")).c_str(), mach.get_length() * 1e3,
                                 mc.size() > 0 ? mc[0] * 1e3 : 0.0, mc.size() > 1 ? mc[1] * 1e3 : 0.0,
                                 removed * 1e9, before.Mass() * 1e9, after.Mass() * 1e9);
                bool alreadyCutByAnEarlierGap = false;
                if (std::abs(mach.get_length()) > 1e-12 && !(removed > 1e-6 * before.Mass())) {
                    // Same tool, pristine piece: did this gap ever have material to remove?
                    GProp_GProps pristineBefore, pristineAfter;
                    BRepGProp::VolumeProperties(pristinePiece, pristineBefore);
                    const TopoDS_Shape pristineCut = builder->applyMachining(pristinePiece, mach, dims);
                    BRepGProp::VolumeProperties(pristineCut, pristineAfter);
                    const double removedFromPristine = pristineBefore.Mass() - pristineAfter.Mass();
                    alreadyCutByAnEarlierGap = removedFromPristine > 1e-6 * pristineBefore.Mass();
                    if (std::getenv("MVB_CORE_DIAG"))
                        std::fprintf(stderr, "[machining] '%s' no-op cut re-tried on the pristine piece: removed %.5f mm3 -> %s\n",
                                     shapeData->get_name().value_or(std::string("?")).c_str(),
                                     removedFromPristine * 1e9,
                                     alreadyCutByAnEarlierGap ? "gap already cut by an overlapping neighbour" : "the tool misses the column");
                }
                if (std::abs(mach.get_length()) > 1e-12 && !(removed > 1e-6 * before.Mass())
                        && !alreadyCutByAnEarlierGap) {
                    throw std::runtime_error("core machining removed nothing: gap length " +
                        std::to_string(mach.get_length() * 1e3) + " mm at (" +
                        (mc.size() > 0 ? std::to_string(mc[0] * 1e3) : std::string("?")) + ", " +
                        (mc.size() > 1 ? std::to_string(mc[1] * 1e3) : std::string("?")) + ") mm on '" +
                        shapeData->get_name().value_or(std::string("?")) + "' -- the tool missed the column (piece volume " +
                        std::to_string(before.Mass() * 1e9) + " -> " + std::to_string(after.Mass() * 1e9) + " mm3)");
                }
            }
        }

        // Apply translation from geometrical description (after machining).
        // Toroidal pieces are rotated {π/2, π/2, 0} upstream so their hole axis
        // points along world Y. MKF still emits the stack offset in coords[2]
        // (MAS Z, the pre-rotation axial) for historical parity with
        // two-piece sets. Swap Y↔Z here so stacked rings displace along the
        // post-rotation axial direction and visually stack along the hole,
        // not perpendicular to it.
        auto coords = piece.get_coordinates();
        if (coords.size() >= 3) {
            if (piece.get_type() == MAS::CoreGeometricalDescriptionElementType::TOROIDAL) {
                shape = translate_shape(shape, coords[0], coords[2], coords[1]);
            } else {
                shape = translate_shape(shape, coords[0], coords[1], coords[2]);
            }
        }


        // Drop phantom sub-solids left by OCCT boolean fuse/cut artifacts
        // (observed on PQ) — but keep every solid that is a meaningful
        // fraction of the largest, so distributed-gap cores (e.g. two
        // central-column subtractive gaps on PQ5050) retain the floating
        // middle chunk between the two gaps.
        {
            std::vector<std::pair<double, TopoDS_Shape>> solids;
            double largestVol = 0.0;
            for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next()) {
                GProp_GProps props;
                BRepGProp::VolumeProperties(e.Current(), props);
                double v = std::abs(props.Mass());
                solids.emplace_back(v, e.Current());
                if (v > largestVol) largestVol = v;
            }
            if (solids.size() > 1) {
                // 1% of the largest volume is well below any real core
                // chunk but well above typical boolean residue.
                const double threshold = largestVol * 0.01;
                TopoDS_Compound compound;
                BRep_Builder builder;
                builder.MakeCompound(compound);
                int kept = 0;
                for (const auto& [v, s] : solids) {
                    if (v >= threshold) { builder.Add(compound, s); ++kept; }
                }
                if (kept >= 1) shape = compound;
            }
        }
        result.push_back(shape);
    }
    return result;
}

// ABT #1249: MKF's pin rails for the coil's bobbin (empty for a bobbin without pins). The MAS
// variant is wrapped without re-processing, so the rails are read from exactly the processed
// description (and pins) the rest of the build draws.
static std::vector<OpenMagnetics::Bobbin::PinRailBlock> getPinRails(const MAS::Coil& coil) {
    if (const auto* bobbin = std::get_if<MAS::Bobbin>(&coil.get_bobbin())) {
        return OpenMagnetics::Bobbin(*bobbin).get_pin_rails();
    }
    return {};
}

static std::vector<OpenMagnetics::Bobbin::PinRailBlock> getPinRails(const OpenMagnetics::Coil& coil) {
    if (const auto* bobbin = std::get_if<OpenMagnetics::Bobbin>(&coil.get_bobbin())) {
        return bobbin->get_pin_rails();
    }
    return {};
}

// ABT #1249 collision gate: a pin rail shares no volume with a core piece or a conductor. It runs
// BEFORE cut_bobbin, which would otherwise carve a rail that runs into the core away without a word.
static void checkPinRailCollisions(const std::vector<OpenMagnetics::Bobbin::PinRailBlock>& rails,
                                   const std::string& bobbinName,
                                   const std::vector<std::pair<std::string, TopoDS_Shape>>& obstacles) {
    if (rails.empty()) return;
    for (const auto& rail : BobbinBuilder::buildPinRailsNamed(rails, bobbinName)) {
        Bnd_Box railBox;
        BRepBndLib::Add(rail.shape, railBox);
        for (const auto& [name, obstacle] : obstacles) {
            if (obstacle.IsNull()) continue;
            Bnd_Box box;
            BRepBndLib::Add(obstacle, box);
            if (railBox.IsOut(box)) continue;
            BRepAlgoAPI_Common common(rail.shape, obstacle);
            if (!common.IsDone())
                throw std::runtime_error("checkPinRailCollisions: boolean common of '" + rail.name + "' and '" + name + "' failed");
            GProp_GProps props;
            BRepGProp::VolumeProperties(common.Shape(), props);
            if (props.Mass() > 1e-15)
                throw std::runtime_error("Pin rail '" + rail.name + "' overlaps '" + name + "' by " +
                                         std::to_string(props.Mass() * 1e9) + " mm^3: MKF's rail geometry and the " +
                                         "core/winding disagree, and cutting the rail away would hide it (ABT #1249)");
        }
    }
}

// Internal bobbin builder. Public surface goes through buildBobbinNamed().
template<typename CoilT>
TopoDS_Shape buildBobbinShape_impl(const CoilT& coil, const MAS::MagneticCore& core,
                                    int polygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) {
    // Toroidal cores have no bobbin: the winding is wound directly on the
    // core. MKF still emits a CoreBobbinProcessedDescription with whatever
    // column dims it can infer, so without this guard we'd build a phantom
    // rectangular bobbin around a ring core and overlap it with everything.
    if (isCoreToroidal(core)) return TopoDS_Shape();
    auto bobbinPd = getBobbinProcessed(coil);
    patchBobbinDimensions(bobbinPd, core);
    if (bobbinPd.get_column_width().value_or(0.0) <= 0.0) return TopoDS_Shape();
    double flangeThickness = bobbinPd.get_wall_thickness();
    if (flangeThickness < 0.0 || std::isnan(flangeThickness)) flangeThickness = 0.0;
    // If MAS/MKF reports both the column-wall thickness and the flange (wall)
    // thickness as zero, there is no bobbin material to render. This is what
    // we see on many MAS examples (PQ3230 buck, ER25/15 planar, E70/33
    // stacked) where MKF emits a CoreBobbinProcessedDescription that is
    // really just a column footprint, not a real bobbin. Building anyway
    // produces a degenerate solid that overlaps the core column (round path,
    // innerR == outerR) or an empty body (rect path, hole == outer).
    double columnThickness = bobbinPd.get_column_thickness();
    if (std::isnan(columnThickness) || columnThickness < 0.0) columnThickness = 0.0;
    if (flangeThickness == 0.0 && columnThickness == 0.0) return TopoDS_Shape();
    // NOTE: we do NOT post-process the bobbin to avoid overlap with the
    // core. The dimensions and position are taken verbatim from MAS/MKF;
    // if the resulting bobbin geometry overlaps a core leg, that means
    // either MAS describes a bobbin whose flange shape doesn't match the
    // actual winding-window cross-section, or BobbinBuilder is using the
    // wrong flange topology. Either way, the fix belongs upstream (in
    // BobbinBuilder or in MAS), not in a silent post-hoc cut here.
    return BobbinBuilder::buildBobbin(bobbinPd, flangeThickness, !isCoreToroidal(core), polygonSegments,
                                      getPinRails(coil));   // ABT #1249
}

// ---- Multi-column placement resolution -------------------------------------------------------
// Maps each turn/insulation layer to the core column it wraps, from the MAS placement fields:
// turn.section -> section.windingWindow (else section.group -> group.windingWindow) -> else
// winding.windingWindow -> else window 0; then the winding window's `column` edge -> column index
// (absent column edge = the main column, the schema default). Anything resolving to the MAIN
// column returns nullopt and takes the unchanged legacy path, so MAS without placement fields is
// byte-identical. Inconsistent placement data throws — no silent fallbacks.
class WoundColumnResolver {
public:
    template<typename CoilT>
    WoundColumnResolver(const CoilT& coil, const MAS::MagneticCore& core,
                        const MAS::CoreBobbinProcessedDescription& bobbinPd) {
        auto corePd = core.get_processed_description();
        if (corePd) {
            _columns = corePd->get_columns();
            for (std::size_t i = 0; i < _columns.size(); ++i) {
                if (_columns[i].get_type() == MAS::ColumnType::CENTRAL) {
                    _mainColumnIndex = i;
                    break;
                }
            }
        }
        // Same NaN/negative treatment as patchBobbinDimensions' wall_thickness: honour
        // "no wall" instead of fabricating one.
        _wall = bobbinPd.get_column_thickness();
        if (std::isnan(_wall) || _wall < 0.0) _wall = 0.0;

        // Winding window index -> column edge. The coil's governing bobbin wins; the core's
        // processed windows cover MAS where only the core carries the edges.
        const auto& bobbinWindows = bobbinPd.get_winding_windows();
        _windowCount = bobbinWindows.size();
        for (std::size_t i = 0; i < bobbinWindows.size(); ++i) {
            if (auto col = bobbinWindows[i].get_column()) {
                _windowColumn[static_cast<int64_t>(i)] = *col;
            }
        }
        if (corePd) {
            const auto coreWindows = corePd->get_winding_windows();
            _windowCount = std::max(_windowCount, coreWindows.size());
            for (std::size_t i = 0; i < coreWindows.size(); ++i) {
                auto col = coreWindows[i].get_column();
                if (col && !_windowColumn.count(static_cast<int64_t>(i))) {
                    _windowColumn[static_cast<int64_t>(i)] = *col;
                }
            }
        }

        for (const auto& winding : coil.get_functional_description()) {
            if (auto ww = winding.get_winding_window()) {
                _windingWindow[winding.get_name()] = *ww;
            }
        }
        std::map<std::string, int64_t> groupWindow;
        if (auto groupsOpt = coil.get_groups_description()) {
            for (const auto& group : *groupsOpt) {
                if (auto ww = group.get_winding_window()) {
                    groupWindow[group.get_name()] = *ww;
                }
            }
        }
        if (auto sectionsOpt = coil.get_sections_description()) {
            for (const auto& section : *sectionsOpt) {
                if (auto ww = section.get_winding_window()) {
                    _sectionWindow[section.get_name()] = *ww;
                } else if (auto grp = section.get_group()) {
                    auto it = groupWindow.find(*grp);
                    if (it != groupWindow.end()) {
                        _sectionWindow[section.get_name()] = it->second;
                    }
                }
            }
        }
    }

    std::optional<TurnBuilder::WoundColumnSpec>
    resolve(const std::optional<std::string>& sectionName, const std::string& windingName,
            const std::string& what) const {
        int64_t windowIndex = 0;
        bool found = false;
        if (sectionName) {
            auto it = _sectionWindow.find(*sectionName);
            if (it != _sectionWindow.end()) {
                windowIndex = it->second;
                found = true;
            }
        }
        if (!found) {
            auto it = _windingWindow.find(windingName);
            if (it != _windingWindow.end()) {
                windowIndex = it->second;
            }
        }
        if (windowIndex == 0 && _windowColumn.empty()) return std::nullopt;   // no placement anywhere
        if (windowIndex < 0 || static_cast<std::size_t>(windowIndex) >= std::max<std::size_t>(_windowCount, 1)) {
            throw std::runtime_error(
                "WoundColumnResolver: '" + what + "' resolves to winding window "
                + std::to_string(windowIndex) + " but the bobbin/core only defines "
                + std::to_string(_windowCount) + " winding windows — inconsistent MAS placement data");
        }
        auto colIt = _windowColumn.find(windowIndex);
        if (colIt == _windowColumn.end()) return std::nullopt;   // no column edge = main column
        int64_t columnIndex = colIt->second;
        if (columnIndex < 0 || static_cast<std::size_t>(columnIndex) >= _columns.size()) {
            throw std::runtime_error(
                "WoundColumnResolver: winding window " + std::to_string(windowIndex)
                + " references column " + std::to_string(columnIndex) + " but the core has "
                + std::to_string(_columns.size()) + " columns — inconsistent MAS placement data");
        }
        if (static_cast<std::size_t>(columnIndex) == _mainColumnIndex) return std::nullopt;

        const auto& column = _columns[static_cast<std::size_t>(columnIndex)];
        const auto& colCoords = column.get_coordinates();
        if (colCoords.empty()) {
            throw std::runtime_error("WoundColumnResolver: wound column "
                                     + std::to_string(columnIndex) + " has no coordinates");
        }
        if (colCoords.size() > 2 && std::abs(colCoords[2]) > 1e-12) {
            throw std::runtime_error(
                "WoundColumnResolver: wound column " + std::to_string(columnIndex)
                + " is depth-displaced (z = " + std::to_string(colCoords[2])
                + " m) — depth-displaced wound columns are not supported yet");
        }
        if (column.get_width() <= 0.0 || column.get_depth() <= 0.0) {
            throw std::runtime_error(
                "WoundColumnResolver: wound column " + std::to_string(columnIndex)
                + " has non-positive width/depth (" + std::to_string(column.get_width()) + "/"
                + std::to_string(column.get_depth()) + " m)");
        }
        TurnBuilder::WoundColumnSpec spec;
        spec.axisX = colCoords[0];
        spec.halfWidth = column.get_width() / 2.0 + _wall;
        spec.halfDepth = column.get_depth() / 2.0 + _wall;
        spec.shape = column.get_shape();
        return spec;
    }

private:
    std::vector<MAS::ColumnElement> _columns;
    std::size_t _mainColumnIndex = 0;
    std::size_t _windowCount = 0;
    double _wall = 0.0;
    std::map<int64_t, int64_t> _windowColumn;
    std::map<std::string, int64_t> _windingWindow;
    std::map<std::string, int64_t> _sectionWindow;
};

// ABT #871: the same placement resolution, keyed by SECTION, for the real-winding conductor
// builder. It works per (winding, parallel) rather than per turn, and needs the answer before
// it starts, so it takes the whole map up front; sections wrapping the main column are left out
// (the builder's own default), which makes the map empty for every single-window design.
template<typename CoilT>
std::map<std::string, TurnBuilder::WoundColumnSpec> resolveWoundColumnsPerSection(
    const CoilT& coil, const MAS::MagneticCore& core,
    const MAS::CoreBobbinProcessedDescription& bobbinPd) {
    std::map<std::string, TurnBuilder::WoundColumnSpec> out;
    auto sectionsOpt = coil.get_sections_description();
    if (!sectionsOpt) return out;
    const WoundColumnResolver resolver(coil, core, bobbinPd);
    const auto sections = *sectionsOpt;
    for (const auto& section : sections) {
        if (auto spec = resolver.resolve(section.get_name(), "", section.get_name())) {
            out[section.get_name()] = *spec;
        }
    }
    return out;
}

// Internal turns builder. Public surface goes through buildTurnsNamed() or
// buildTurnsNamedFromTurns(). The template covers both the MAS and the
// OpenMagnetics coil/wire variants used internally by buildAllNamed.
// emitCoatingShells: when true, emit TWO concentric solids per turn -- the bare COPPER core
// (named "<turn>") and the OUTER insulated footprint (named "<turn> coating"). They overlap (copper
// inside outer); the downstream mesher fragments them into a copper volume + the enamel annulus, so
// a thermal FEA gets the real low-k wire coating between turns. Reuses buildTurn(), so it works for
// every column shape (round/rect/oblong/toroidal) and wire type (round/rect/litz) unchanged.
template<typename CoilT, typename WireT>
std::vector<TopoDS_Shape> buildTurnsImpl(const CoilT& coil, const MAS::MagneticCore& core, std::vector<std::string>* outNames, int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS, bool paintCoating = true, bool emitCoatingShells = false) {
    std::vector<TopoDS_Shape> result;
    auto turnsOpt = coil.get_turns_description();
    if (!turnsOpt || turnsOpt->empty()) return result;

    TurnBuilder::clearCache();  // Reset per-call so different magnetics don't share cached shapes
    bool toroidal = isCoreToroidal(core);
    auto bobbinPd = getBobbinProcessed(coil);
    patchBobbinDimensions(bobbinPd, core);
    const WoundColumnResolver woundColumns(coil, core, bobbinPd);

    const auto& funcDesc = coil.get_functional_description();
    std::map<std::string, MAS::Wire> wireMap;
    for (const auto& winding : funcDesc) {
        const auto& wireVar = winding.get_wire();
        if (std::holds_alternative<std::string>(wireVar)) {
            wireMap[winding.get_name()] = OpenMagnetics::find_wire_by_name(std::get<std::string>(wireVar));
        } else {
            wireMap[winding.get_name()] = std::get<WireT>(wireVar);
        }
    }

    size_t turnIdx = 0;
    for (const auto& turn : *turnsOpt) {
        const std::string baseName = turn.get_name().empty()
                                        ? ("Turn_" + std::to_string(turnIdx)) : turn.get_name();
        auto it = wireMap.find(turn.get_winding());
        if (it == wireMap.end()) {
            throw std::runtime_error(
                "buildTurnsImpl: turn '" + baseName + "' references winding '" + turn.get_winding()
                + "' which has no wire in coil.functionalDescription — refusing to draw it with "
                  "an invented default wire");
        }
        const MAS::Wire& wire = it->second;
        const auto woundColumn = woundColumns.resolve(turn.get_section(), turn.get_winding(), baseName);
        auto emit = [&](bool coat, const std::string& suffix) {
            TopoDS_Shape t = TurnBuilder::buildTurn(turn, wire, bobbinPd, toroidal, wirePolygonSegments,
                                                    DEFAULT_WIRE_REVOLUTION_SEGMENTS, coat, woundColumn);
            if (!t.IsNull()) {
                result.push_back(t);
                if (outNames) outNames->push_back(baseName + suffix);
            }
        };
        if (emitCoatingShells) {
            emit(false, "");           // bare copper core
            emit(true, " coating");    // outer insulated footprint -> enamel shell after fragment
        } else {
            emit(paintCoating, "");
        }
        ++turnIdx;
    }
    return result;
}

// Build INSULATION-layer solids (inter-layer / inter-section tape). Each insulation layer is a thin
// ring spanning the layer footprint; we synthesize a RECTANGULAR "turn" from the layer's
// coordinates+dimensions and reuse buildTurn(), so every column shape (round/rect/oblong/toroidal)
// is covered. Zero-thickness layers (logical placeholders with no material/thickness) are SKIPPED
// (no silent fabrication). Names "insulation_layer_<i>" so a thermal FEA can tag them as the low-k
// conduction barrier between windings (unlike wire enamel, this sits on the winding<->winding
// conduction path with no convection film in series, so it can matter thermally).
template<typename CoilT>
std::vector<NamedShape> buildInsulationLayersImpl(const CoilT& coil, const MAS::MagneticCore& core,
                                                  int wirePolygonSegments) {
    std::vector<NamedShape> out;
    auto layersOpt = coil.get_layers_description();
    if (!layersOpt || layersOpt->empty()) return out;
    const bool toroidal = isCoreToroidal(core);
    auto bobbinPd = getBobbinProcessed(coil);
    patchBobbinDimensions(bobbinPd, core);
    const WoundColumnResolver woundColumns(coil, core, bobbinPd);
    int idx = 0;
    for (const auto& layer : *layersOpt) {
        const int i = idx++;
        if (layer.get_type() != MAS::ElectricalType::INSULATION) continue;
        const auto& dims = layer.get_dimensions();
        const auto& coords = layer.get_coordinates();
        if (dims.size() < 2 || coords.size() < 2) continue;
        const double w = dims[0], h = dims[1];
        if (w <= 1e-9 || h <= 1e-9) continue;   // zero-thickness placeholder -> nothing physical to build
        MAS::Turn turn;
        turn.set_name("insulation_layer_" + std::to_string(i));
        turn.set_coordinates(std::vector<double>{coords[0], coords[1]});
        turn.set_dimensions(std::vector<double>{w, h});
        turn.set_cross_sectional_shape(MAS::TurnCrossSectionalShape::RECTANGULAR);
        // Toroidal layers carry their outer-ring crossing in additionalCoordinates (written
        // by MKF); the synthesized turn must carry it too — TurnBuilder refuses to invent it.
        if (const auto& layerAdd = layer.get_additional_coordinates()) {
            turn.set_additional_coordinates(*layerAdd);
        }
        MAS::Wire wire; wire.set_type(MAS::WireType::RECTANGULAR);
        MAS::DimensionWithTolerance ww; ww.set_nominal(w);
        MAS::DimensionWithTolerance hh; hh.set_nominal(h);
        wire.set_outer_width(std::optional<MAS::DimensionWithTolerance>(ww));
        wire.set_outer_height(std::optional<MAS::DimensionWithTolerance>(hh));
        wire.set_conducting_width(std::optional<MAS::DimensionWithTolerance>(ww));
        wire.set_conducting_height(std::optional<MAS::DimensionWithTolerance>(hh));
        // Insulation layers in a lateral window wrap that window's column; the layer's
        // section carries the placement (insulation layers have no winding of their own).
        const auto woundColumn = woundColumns.resolve(layer.get_section(), "", layer.get_name());
        TopoDS_Shape s = TurnBuilder::buildTurn(turn, wire, bobbinPd, toroidal,
                                                wirePolygonSegments, DEFAULT_WIRE_REVOLUTION_SEGMENTS, true,
                                                woundColumn);
        if (!s.IsNull()) out.push_back({s, "insulation_layer_" + std::to_string(i), Role::Insulation});
    }
    return out;
}

// ABT #1169. A turn solid drawn at its OUTER insulated footprint is emitted with a " coating"
// suffix by buildTurnsImpl(emitCoatingShells=true); the copper body keeps the bare name. Role
// follows that one distinction, so a consumer never has to re-parse the suffix itself.
Role turn_role_for(const std::string& name) {
    static const std::string kCoating = " coating";
    const bool coated = name.size() > kCoating.size() &&
        name.compare(name.size() - kCoating.size(), kCoating.size(), kCoating) == 0;
    return coated ? Role::TurnCoating : Role::Turn;
}

} // anonymous namespace

// ---- Named-shape overloads ------------------------------------------------

std::vector<NamedShape> MagneticBuilder::buildCoreNamed(const MAS::MagneticCore& core,
                                                         int corePolygonSegments) const {
    auto shapes = buildCoreShapes_impl(core, corePolygonSegments);
    std::vector<NamedShape> out;
    out.reserve(shapes.size());
    const std::string base = core.get_name().value_or("Core");
    const bool single = shapes.size() == 1;
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        out.emplace_back(shapes[i],
                         single ? base : base + "_" + std::to_string(i), Role::Core);
    }
    return out;
}

std::vector<NamedShape> MagneticBuilder::buildTurnsNamed(const MAS::Coil& coil,
                                                         const MAS::MagneticCore& core,
                                                         int wirePolygonSegments,
                                                         bool paintCoating) const {
    std::vector<std::string> names;
    auto shapes = buildTurnsImpl<MAS::Coil, MAS::Wire>(coil, core, &names, wirePolygonSegments, paintCoating);
    std::vector<NamedShape> out;
    out.reserve(shapes.size());
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        const std::string n = (i < names.size() && !names[i].empty())
                                ? names[i]
                                : "Turn_" + std::to_string(i);
        out.emplace_back(shapes[i], n, turn_role_for(n));
    }
    return out;
}

std::vector<NamedShape> MagneticBuilder::buildTurnsNamed(const OpenMagnetics::Coil& coil,
                                                         const MAS::MagneticCore& core,
                                                         int wirePolygonSegments,
                                                         bool paintCoating) const {
    std::vector<std::string> names;
    auto shapes = buildTurnsImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(coil, core, &names, wirePolygonSegments, paintCoating);
    std::vector<NamedShape> out;
    out.reserve(shapes.size());
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        const std::string n = (i < names.size() && !names[i].empty())
                                ? names[i]
                                : "Turn_" + std::to_string(i);
        out.emplace_back(shapes[i], n, turn_role_for(n));
    }
    return out;
}

NamedShape MagneticBuilder::buildBobbinNamed(const MAS::Coil& coil,
                                             const MAS::MagneticCore& core,
                                             int corePolygonSegments) const {
    NamedShape ns;
    ns.shape = buildBobbinShape_impl(coil, core, corePolygonSegments);
    ns.name = getBobbinNameT<MAS::Bobbin>(coil.get_bobbin(), "Bobbin");
    ns.role = Role::Bobbin;
    return ns;
}

NamedShape MagneticBuilder::buildBobbinNamed(const OpenMagnetics::Coil& coil,
                                             const MAS::MagneticCore& core,
                                             int corePolygonSegments) const {
    NamedShape ns;
    ns.shape = buildBobbinShape_impl(coil, core, corePolygonSegments);
    ns.name = getBobbinNameT<OpenMagnetics::Bobbin>(coil.get_bobbin(), "Bobbin");
    ns.role = Role::Bobbin;
    return ns;
}

template<typename CoilT, typename WireT>
static std::vector<NamedShape> buildTurnsNamedImpl(const CoilT& coil, const MAS::MagneticCore& core,
                                                    int wirePolygonSegments) {
    std::vector<std::string> names;
    auto shapes = buildTurnsImpl<CoilT, WireT>(coil, core, &names, wirePolygonSegments);
    std::vector<NamedShape> out;
    out.reserve(shapes.size());
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        const std::string n = (i < names.size() && !names[i].empty())
                                ? names[i]
                                : "Turn_" + std::to_string(i);
        out.emplace_back(shapes[i], n, turn_role_for(n));
    }
    return out;
}

// Apply up to `numPlanes` symmetry cuts — defined in Symmetry.cpp.

std::vector<NamedShape> MagneticBuilder::buildAllNamed(const MAS::Magnetic& magnetic,
                                                         bool includeBobbin,
                                                         int symmetryPlanes,
                                                         int wirePolygonSegments,
                                                         int corePolygonSegments,
                                                         bool paintCoating,
                                                         bool emitCoatingShells,
                                                         bool includeInsulation,
                                                         double coreCoatingThickness,
                                                         bool useRealWindingGeometry,
                                                         bool femReady,
                                                         bool skipGeometryChecks) const {
    // MAS 1.x makes Magnetic.core / Magnetic.coil optional, but this builder
    // requires both present. The generated getters return the optional BY VALUE,
    // so bind COPIES (not references — a reference would dangle past the temporary).
    // Report WHICH one is missing. A bare .value() here throws std::bad_optional_access, whose
    // message names neither the field nor the builder, and the usual cause is a caller handing in a
    // full MAS file (masVersion/inputs/magnetic/outputs) instead of the nested magnetic object --
    // that parses cleanly to a Magnetic with both optionals empty.
    if (!magnetic.get_core().has_value() || !magnetic.get_coil().has_value())
        throw std::runtime_error(
            std::string("buildAllNamed: the magnetic has no ") +
            (!magnetic.get_core().has_value() ? "core" : "coil") +
            ". Both are required. If the input is a full MAS file, pass its 'magnetic' member.");
    const MAS::MagneticCore core = magnetic.get_core().value();
    const MAS::Coil coil = magnetic.get_coil().value();
    // If geometricalDescription is already present, skip expensive MKF
    // autocomplete and use the pre-enriched data directly.
    auto geoOpt = core.get_geometrical_description();
    if (geoOpt && geoOpt.has_value() && !geoOpt->empty()) {
        if (useRealWindingGeometry) {
            throw std::runtime_error(
                "buildAllNamed: useRealWindingGeometry requires MKF to (re-)wind the magnetic, "
                "but the input already carries a geometricalDescription — pass the raw "
                "functional design instead (caller-provided turn positions are never "
                "silently re-interpreted)");
        }
        // Build directly from MAS types — no MKF enrichment needed.
        auto all = buildCoreNamed(core, corePolygonSegments);

        if (coreCoatingThickness > 0.0) {   // conformal core-coating shells (offset core - core)
            std::vector<NamedShape> coatings;
            for (auto& ns : all) {
                auto shell = buildCoreCoatingShell(ns.shape, coreCoatingThickness);
                if (!shell.IsNull()) coatings.push_back({shell, ns.name + " coating", Role::CoreCoating});
            }
            for (auto& c : coatings) all.push_back(std::move(c));
        }

        std::vector<std::string> turnNames;
        auto turnShapes = buildTurnsImpl<MAS::Coil, MAS::Wire>(
            coil, core, &turnNames, wirePolygonSegments, paintCoating, emitCoatingShells);

        if (includeBobbin) {
            auto bobbin = buildBobbinNamed(coil, core, corePolygonSegments);
            if (!bobbin.shape.IsNull()) {
                std::vector<TopoDS_Shape> cutters;
                for (const auto& ns : all) cutters.push_back(ns.shape);
                cutters.insert(cutters.end(), turnShapes.begin(), turnShapes.end());
                {   // ABT #1249: rails vs core and copper, before the cut could hide an overlap.
                    std::vector<std::pair<std::string, TopoDS_Shape>> obstacles;
                    for (const auto& ns : all) obstacles.emplace_back(ns.name, ns.shape);
                    for (std::size_t i = 0; i < turnShapes.size(); ++i)
                        obstacles.emplace_back(i < turnNames.size() ? turnNames[i] : "Turn_" + std::to_string(i), turnShapes[i]);
                    checkPinRailCollisions(getPinRails(coil), bobbin.name, obstacles);
                }
                // ABT #1171: a pin owns the flange volume it passes through (see PinBuilder.h).
                for (auto& pinShape : PinBuilder::buildPins(getBobbinProcessed(coil)))
                    cutters.push_back(std::move(pinShape));
                // ABT #1176: a magnetic shunt owns the former volume it passes through.
                for (auto& shunt : ShuntBuilder::buildShuntsNamed(magnetic))
                    cutters.push_back(std::move(shunt.shape));
                bobbin.shape = cut_bobbin(bobbin.shape, cutters);
                if (!bobbin.shape.IsNull()) all.push_back(bobbin);
            }
        }

        appendAccessorySolids(all, magnetic,
                              AccessoryOptions{includeBobbin, wirePolygonSegments,
                                               corePolygonSegments, paintCoating,
                                               useRealWindingGeometry, femReady});

        for (std::size_t i = 0; i < turnShapes.size(); ++i) {
            const std::string n = (i < turnNames.size() && !turnNames[i].empty())
                                    ? turnNames[i]
                                    : "Turn_" + std::to_string(i);
            all.emplace_back(turnShapes[i], n, turn_role_for(n));
        }

        if (includeInsulation) {
            auto ins = buildInsulationLayersImpl<MAS::Coil>(coil, core, wirePolygonSegments);
            for (auto& ns : ins) all.push_back(std::move(ns));
        }

        // Planar (PCB) coils get an FR4 substrate board. Patch the bobbin
        // processed description first so column_shape/width/depth are
        // populated from the core when MKF left the bobbin variant empty.
        if (auto groupsOpt = coil.get_groups_description();
            groupsOpt && !groupsOpt->empty()) {
            auto bobbinPd = getBobbinProcessed(coil);
            patchBobbinDimensions(bobbinPd, core);
            auto fr4 = FR4Builder::buildFR4Board(*groupsOpt, bobbinPd);
            if (!fr4.IsNull()) {
                // ABT #1176: the board is one slab over the whole group height (the copper sits
                // inside it); a sheet laid in the stack takes its volume out of the board.
                NamedShape board{fr4, "FR4Board", Role::FR4};
                cutByShunts(board, all);
                all.push_back(std::move(board));
            }
        }

        checkShuntCollisions(all);   // ABT #1176: on the finished assembly, before symmetry
        // ABT #1248: a toroid is placed for its mounting, exactly as in the MKF-enriched overload.
        if (isCoreToroidal(core)) {
            const auto frame = toroidMountingFrameOf(magnetic);
            if (frame.mounting == ConductorBuilder::ToroidMounting::Vertical)
                for (auto& ns : all)
                    ns.shape = BRepBuilderAPI_Transform(ns.shape, frame.toExported).Shape();
        }
        return apply_symmetry(std::move(all), symmetryPlanes);
    }

    OpenMagnetics::Magnetic enriched = magnetic_autocomplete_safe(magnetic, useRealWindingGeometry);
    return buildAllNamed(enriched, includeBobbin, symmetryPlanes,
                         wirePolygonSegments, corePolygonSegments, paintCoating, emitCoatingShells,
                         includeInsulation, coreCoatingThickness, useRealWindingGeometry, femReady,
                         skipGeometryChecks);
}

std::vector<NamedShape> MagneticBuilder::buildAllNamed(const OpenMagnetics::Magnetic& magnetic,
                                                         bool includeBobbin,
                                                         int symmetryPlanes,
                                                         int wirePolygonSegments,
                                                         int corePolygonSegments,
                                                         bool paintCoating,
                                                         bool emitCoatingShells,
                                                         bool includeInsulation,
                                                         double coreCoatingThickness,
                                                         bool useRealWindingGeometry,
                                                         bool femReady,
                                                         bool skipGeometryChecks) const {
    auto all = buildCoreNamed(magnetic.get_core(), corePolygonSegments);

    if (coreCoatingThickness > 0.0) {   // conformal core-coating shells (offset core - core)
        std::vector<NamedShape> coatings;
        for (auto& ns : all) {
            auto shell = buildCoreCoatingShell(ns.shape, coreCoatingThickness);
            if (!shell.IsNull()) coatings.push_back({shell, ns.name + " coating", Role::CoreCoating});
        }
        for (auto& c : coatings) all.push_back(std::move(c));
    }

    // Build the winding solids ONCE, then reuse for both bobbin-cutting and the assembly.
    // Real winding: one continuous conductor per (winding, parallel) replaces the per-turn
    // closed loops; positions come verbatim from the MKF-wound coil.
    std::vector<std::string> turnNames;
    std::vector<TopoDS_Shape> turnShapes;
    // ABT #685: per-solid names, kept alongside the shape/name vectors this path splits things
    // into. Indexed like turnShapes/turnNames; empty for anything that has none.
    std::vector<std::vector<std::string>> turnPartNames;
    // ABT #1245: the role each conductor solid was built with. The real-winding builder emits
    // terminal caps (Role::Terminal) and foil solder bodies (Role::Solder) beside the copper;
    // re-deriving the role from the name at the append below turned all of them into Turn.
    std::vector<Role> turnRoles;
    double toroidTerminalPlane = std::numeric_limits<double>::quiet_NaN();   // ABT #1173
    // ABT #1173: a base this assembly cannot draw is refused before any copper is built.
    if (includeBobbin) {
        if (const auto base = BaseBuilder::baseOf<OpenMagnetics::Bobbin>(magnetic.get_coil().get_bobbin()))
            BaseBuilder::requireDrawableMounting(
                base.value(), getBobbinNameT<OpenMagnetics::Bobbin>(magnetic.get_coil().get_bobbin(), "Bobbin"));
    }
    if (useRealWindingGeometry) {
        for (auto& ns : buildRealWindingConductorsNamed(magnetic, all, wirePolygonSegments,
                                                        paintCoating, emitCoatingShells, femReady,
                                                        skipGeometryChecks,
                                                        /*cutterOnly=*/false,
                                                        &toroidTerminalPlane)) {
            turnShapes.push_back(ns.shape);
            turnNames.push_back(ns.name);
            turnPartNames.push_back(std::move(ns.partNames));
            turnRoles.push_back(ns.role);
        }
    } else {
        turnShapes = buildTurnsImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(
            magnetic.get_coil(), magnetic.get_core(), &turnNames, wirePolygonSegments, paintCoating, emitCoatingShells);
        // Per-turn loops carry no role of their own; the " coating" suffix is the one distinction.
        for (std::size_t i = 0; i < turnShapes.size(); ++i)
            turnRoles.push_back(turn_role_for(i < turnNames.size() ? turnNames[i] : std::string()));
    }

    if (includeBobbin) {
        auto bobbin = buildBobbinNamed(magnetic.get_coil(), magnetic.get_core(), corePolygonSegments);
        if (!bobbin.shape.IsNull()) {
            std::vector<TopoDS_Shape> cutters;
            for (const auto& ns : all) cutters.push_back(ns.shape);
            {   // ABT #1249: rails vs core and copper, before the cut could hide an overlap.
                std::vector<std::pair<std::string, TopoDS_Shape>> obstacles;
                for (const auto& ns : all) obstacles.emplace_back(ns.name, ns.shape);
                for (std::size_t i = 0; i < turnShapes.size(); ++i)
                    obstacles.emplace_back(i < turnNames.size() ? turnNames[i] : "Turn_" + std::to_string(i), turnShapes[i]);
                checkPinRailCollisions(getPinRails(magnetic.get_coil()), bobbin.name, obstacles);
            }
            // The bobbin is ALWAYS cut with the OUTER (coated) envelope. When the product
            // paints bare copper, cutting with the copper solids makes every lead
            // pass-through slot EXACTLY tangent to the conductor -- measured on 10_emi
            // (EP 13): a 0.043 mm^2 copper-on-bobbin contact patch whose fragment interface
            // degenerates ("overlapping facets", 12/17 volumes lost). Physically the enamel
            // occupies the corridor, so the slot clearance is the coating thickness
            // (ro - rc), never zero and never a made-up margin. The cutter build runs in
            // DRAWING mode (femReady=false): a cutting tool needs no welds, no mitre
            // certification and may overlap itself, so the cheap per-run compound is enough.
            if (!paintCoating && useRealWindingGeometry) {
                for (auto& ns : buildRealWindingConductorsNamed(
                         magnetic, all, wirePolygonSegments,
                         /*paintCoating=*/true, /*emitCoatingShells=*/false,
                         /*femReady=*/false, /*diagnosticSkipCollisionCheck=*/false,
                         /*cutterOnly=*/true)) {
                    cutters.push_back(ns.shape);
                }
            } else {
                cutters.insert(cutters.end(), turnShapes.begin(), turnShapes.end());
            }
            // ABT #1171: a pin owns the flange volume it passes through (see PinBuilder.h).
            for (auto& pinShape : PinBuilder::buildPins(getBobbinProcessed(magnetic.get_coil())))
                cutters.push_back(std::move(pinShape));
            // ABT #1176: a magnetic shunt owns the former volume it passes through.
            for (auto& shunt : ShuntBuilder::buildShuntsNamed(magnetic))
                cutters.push_back(std::move(shunt.shape));
            bobbin.shape = cut_bobbin(bobbin.shape, cutters);
            if (!bobbin.shape.IsNull()) all.push_back(bobbin);
        }
    }

    appendAccessorySolids(all, magnetic,
                          AccessoryOptions{includeBobbin, wirePolygonSegments,
                                           corePolygonSegments, paintCoating,
                                           useRealWindingGeometry, femReady,
                                           toroidTerminalPlane});

    // Append the already-built turns (no second build).
    for (std::size_t i = 0; i < turnShapes.size(); ++i) {
        const std::string n = (i < turnNames.size() && !turnNames[i].empty())
                                ? turnNames[i]
                                : "Turn_" + std::to_string(i);
        if (i < turnPartNames.size() && !turnPartNames[i].empty()) {
            all.emplace_back(turnShapes[i], n, turnPartNames[i], turnRoles.at(i));   // ABT #685
        }
        else {
            all.emplace_back(turnShapes[i], n, turnRoles.at(i));
        }
    }

    if (includeInsulation) {
        auto ins = buildInsulationLayersImpl<OpenMagnetics::Coil>(magnetic.get_coil(), magnetic.get_core(), wirePolygonSegments);
        for (auto& ns : ins) all.push_back(std::move(ns));
    }

    if (auto groupsOpt = magnetic.get_coil().get_groups_description();
        groupsOpt && !groupsOpt->empty()) {
        auto bobbinPd = getBobbinProcessed(magnetic.get_coil());
        patchBobbinDimensions(bobbinPd, magnetic.get_core());
        auto fr4 = FR4Builder::buildFR4Board(*groupsOpt, bobbinPd);
        if (!fr4.IsNull()) {
            // ABT #1176: see the MAS overload.
            NamedShape board{fr4, "FR4Board", Role::FR4};
            cutByShunts(board, all);
            all.push_back(std::move(board));
        }
    }

    checkShuntCollisions(all);   // ABT #1176: on the finished assembly, before rotation/symmetry
    checkBaseCollisions(all);    // ABT #1173: likewise for a toroid base

    // MKF's geometricalDescription rotates the toroid by {pi/2, pi/2, 0} (Core.cpp), tipping
    // the ring out of the MAS XY plane: it lands in XZ with the hole axis along world Y, and
    // the turns are placed to match. Counter-rotate the WHOLE assembled magnetic (core + turns
    // together, so they stay consistent) by -pi/2 about X, restoring the MAS convention: the
    // ring lies in XY with the hole axis along Z. 2D consumers can then project the toroid's
    // magnetic path in XY like any other core.
    const bool isToroidal = [&]{
        auto geo = magnetic.get_core().get_geometrical_description();
        if (!geo) return false;
        for (const auto& p : *geo)
            if (p.get_type() == MAS::CoreGeometricalDescriptionElementType::TOROIDAL) return true;
        return false;
    }();
    // ABT #1248: that counter-rotation is now the mounting. HORIZONTAL keeps MKF's frame (hole
    // axis Y, the ring flat, every terminal lead dropping in -Y); VERTICAL stands the ring on its
    // rim (hole axis Z) turned so the terminals are at the bottom and their leads run in -Y.
    if (isToroidal) {
        const auto frame = toroidMountingFrameOf(magnetic);
        if (frame.mounting == ConductorBuilder::ToroidMounting::Vertical)
            for (auto& ns : all) ns.shape = BRepBuilderAPI_Transform(ns.shape, frame.toExported).Shape();
    }

    // FEM product: the periodic-surface -> B-spline re-expression (ABT #490 class) is NOT done
    // here any more. It has to run in the frame the STEP is written in, so it lives in
    // exportSTEP (StepExportOptions::nurbsPeriodicSolids), which every FEM consumer goes
    // through. Converting here, in metres, and scaling x1000 at export left every
    // length-parametrised B-spline (facet strips, caps) with knots 1000x too compressed for its
    // geometry: BOPAlgo then reported sporadic self-intersections on faceted revolves (ABT #1111
    // -- 11 of 38 designs at --segments 12; the same solids are clean rescaled to metres).

    return apply_symmetry(std::move(all), symmetryPlanes);
}

// ---- ABT #1169 (WP0): accessory-solid hook --------------------------------
//
// Deliberately empty. See MagneticBuilder::AccessoryOptions in the header for what each
// later work package attaches here and the naming/role contract it must honour. The two
// overloads exist because buildAllNamed has a MAS-typed (pre-enriched geometricalDescription)
// path and an MKF-enriched path; an accessory that only MKF can size belongs in the second.
void MagneticBuilder::appendAccessorySolids(std::vector<NamedShape>& all,
                                            const MAS::Magnetic& magnetic,
                                            const AccessoryOptions& opts) const {
    // buildAllNamed has already refused a magnetic without a core or a coil, so both are engaged.
    appendSpacerSolids(all, magnetic.get_core().value());  // ABT #1170 (WP1)
    // ABT #1176 (WP7): magnetic shunts, after the spacers so a sheet held in an additive gap cuts
    // its shims. See ShuntBuilder.h for the MKF validation contract and the naming decision.
    appendShuntSolids(all, magnetic);

    // ABT #1171 (WP2): the bobbin's solder pins, exactly where MKF placed them. They belong to
    // the former, so they are drawn only with it. Absent `pins` means the catalogue record
    // states no footprint: nothing to draw and nothing to guess.
    if (opts.includeBobbin) {
        const MAS::Coil coil = magnetic.get_coil().value();
        // ABT #1173 (WP4): a toroid base sits under the real-winding terminal plane, which only the
        // MKF-enriched build knows; a magnetic arriving with its geometry already described has none.
        if (BaseBuilder::baseOf<MAS::Bobbin>(coil.get_bobbin()))
            throw std::runtime_error(
                "appendAccessorySolids: bobbin '" + getBobbinNameT<MAS::Bobbin>(coil.get_bobbin(), "Bobbin") +
                "' is a toroid base, which is drawn only from the MKF-enriched magnetic with real winding "
                "geometry (its top face is the terminal plane of the real-winding leads). Pass the functional "
                "design without geometricalDescription and set useRealWindingGeometry.");
        for (auto& pin : PinBuilder::buildPinsNamed(getBobbinProcessed(coil),
                                                    getBobbinNameT<MAS::Bobbin>(coil.get_bobbin(), "Bobbin")))
            all.push_back(std::move(pin));
    }
}

void MagneticBuilder::appendAccessorySolids(std::vector<NamedShape>& all,
                                            const OpenMagnetics::Magnetic& magnetic,
                                            const AccessoryOptions& opts) const {
    appendSpacerSolids(all, magnetic.get_core());          // ABT #1170 (WP1)
    appendShuntSolids(all, magnetic);                      // ABT #1176 (WP7): see the MAS overload

    // ABT #1171 (WP2): see the MAS overload.
    if (opts.includeBobbin) {
        const auto& coil = magnetic.get_coil();
        const std::string bobbinName = getBobbinNameT<OpenMagnetics::Bobbin>(coil.get_bobbin(), "Bobbin");
        // ABT #1173 (WP4): a toroid base and ITS pins (BaseBuilder.h: top face on the terminal plane,
        // MKF's pins moved along Y onto its bottom face), instead of the pins at MKF's bare-ring height.
        if (BaseBuilder::baseOf<OpenMagnetics::Bobbin>(coil.get_bobbin())) {
            if (!opts.useRealWindingGeometry || std::isnan(opts.toroidTerminalPlaneY))
                throw std::runtime_error(
                    "appendAccessorySolids: bobbin '" + bobbinName + "' is a toroid base; it is drawn only "
                    "with real winding geometry, whose terminal plane its top face sits on.");
            const auto functional = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin()).get_functional_description();
            for (auto& solid : BaseBuilder::buildBaseNamed(functional.value(), getBobbinProcessed(coil),
                                                           opts.toroidTerminalPlaneY, bobbinName))
                all.push_back(std::move(solid));
            return;
        }
        for (auto& pin : PinBuilder::buildPinsNamed(getBobbinProcessed(coil),
                                                    getBobbinNameT<OpenMagnetics::Bobbin>(coil.get_bobbin(), "Bobbin")))
            all.push_back(std::move(pin));
    }
}

// ---- Standalone builders for the unified bindings API ---------------------

NamedShape MagneticBuilder::buildCorePieceNamed(const MAS::CoreShape& shape,
                                                  int corePolygonSegments) const {
    // Validate / process the shape via MKF (computes effective parameters
    // and per-piece data; throws if the shape data is malformed).
    auto corePiece = OpenMagnetics::CorePiece::factory(shape, /*process=*/true);
    if (!corePiece) {
        throw std::runtime_error(
            "buildCorePieceNamed: OpenMagnetics::CorePiece::factory returned null");
    }

    auto family = shape.get_family();
    std::string subtype = shape.get_family_subtype().value_or("");
    auto builder = shapes::createShapeBuilder(family, subtype, corePolygonSegments);
    if (!builder) {
        throw std::runtime_error(
            "buildCorePieceNamed: no geometry builder for family '"
            + core_shape_family_to_string(family) + "'");
    }

    TopoDS_Shape geom = builder->buildPiece(shape);
    if (geom.IsNull()) {
        throw std::runtime_error(
            "buildCorePieceNamed: builder produced null shape for '"
            + shape.get_name().value_or(core_shape_family_to_string(family)) + "'");
    }

    std::string name = shape.get_name().value_or(core_shape_family_to_string(family));
    return NamedShape{geom, name, Role::Core};
}

NamedShape MagneticBuilder::buildBobbinNamedFromBobbin(const MAS::Bobbin& bobbin,
                                                       bool axisIsY,
                                                       int polygonSegments) const {
    auto pdOpt = bobbin.get_processed_description();
    if (!pdOpt) {
        throw std::runtime_error(
            "buildBobbinNamedFromBobbin: bobbin.processedDescription is required "
            "(MAS Bobbin is not enriched). Use OpenMagnetics::Bobbin::process_data() "
            "or feed a fully-populated MAS::Magnetic to drawMagnetic instead.");
    }
    if (pdOpt->get_column_width().value_or(0.0) <= 0.0) {
        throw std::runtime_error(
            "buildBobbinNamedFromBobbin: processedDescription.columnWidth must be > 0");
    }
    double flangeThickness = pdOpt->get_wall_thickness();
    if (flangeThickness < 0.0 || std::isnan(flangeThickness)) flangeThickness = 0.0;
    std::string name = bobbin.get_name().value_or("Bobbin");
    TopoDS_Shape s = BobbinBuilder::buildBobbin(*pdOpt, flangeThickness, axisIsY, polygonSegments);
    // "No bobbin to draw" — e.g. a planar / PCB winding whose column has zero wall
    // and column thickness, so the box body cut (outer minus same-size hole) leaves
    // an empty compound. That is NOT an error: normalise it to a null shape and let
    // the caller drop it, instead of emitting an empty 84-byte STL or (previously)
    // throwing, which surfaced as a spurious "[mvbpp] unknown C++ exception" on
    // every planar core. An OBLONG box-minus-box yields a non-null but solid-less
    // compound, so check for solids rather than just IsNull().
    bool hasSolid = false;
    if (!s.IsNull()) {
        TopExp_Explorer exp(s, TopAbs_SOLID);
        hasSolid = exp.More();
    }
    if (!hasSolid) {
        return NamedShape{TopoDS_Shape(), name, Role::Bobbin};
    }
    return NamedShape{s, name, Role::Bobbin};
}

// The conductor options every real-winding consumer builds from, so the emitted conductors and
// the lead-length measurement (ABT #1215) can never be planned with different inputs.
ConductorBuilder::Options MagneticBuilder::realWindingConductorOptions(
    const OpenMagnetics::Magnetic& magnetic, const std::vector<NamedShape>& coreShapes,
    int wirePolygonSegments, bool femReady, MAS::CoreBobbinProcessedDescription& bobbinPd,
    bool& toroidalCore) const {
    bobbinPd = getBobbinProcessed(magnetic.get_coil());
    patchBobbinDimensions(bobbinPd, magnetic.get_core());
    toroidalCore = isCoreToroidal(magnetic.get_core());
    ConductorBuilder::Options copts;
    if (toroidalCore) {   // ABT #1248
        copts.toroidMounting = toroidMountingFrameOf(magnetic);
        copts.toroidCoreDepthAlongDown = coreDepthAlongDown(coreShapes, copts.toroidMounting->down);
    }
    copts.wirePolygonSegments = wirePolygonSegments;
    copts.femReady = femReady;   // OM drawing -> fast compound; FEM export -> one-piece/conformal
    // Hand the CORE solids to the conductor builder so it can aim the terminal leads at
    // the true window opening (classified from the real geometry -- column metadata
    // under-describes cores like PQ whose plates wrap most of the perimeter; measured on
    // 03_buck_pq3230, whose lead tip landed on an oblique plate face at both column-derived
    // azimuths).
    if (!toroidalCore)
        for (const auto& ns : coreShapes) copts.coreObstacles.push_back(ns.shape);
    // ABT #871: which core column each section wraps. The conductor builder is handed the coil
    // and the bobbin only, and the leg geometry lives on the CORE — so the placement chain is
    // resolved here, with the same resolver the ideal turn path already uses.
    copts.woundColumnPerSection =
        resolveWoundColumnsPerSection(magnetic.get_coil(), magnetic.get_core(), bobbinPd);
    return copts;
}

std::map<std::string, ConductorBuilder::TerminalLeadLength>
MagneticBuilder::measureTerminalLeadLengths(const OpenMagnetics::Magnetic& magnetic,
                                            bool paintCoating, bool femReady,
                                            int wirePolygonSegments, int corePolygonSegments,
                                            double coreCoatingThickness) const {
    // The obstacles buildAllNamed hands the conductor builder: the core pieces, plus their
    // coating shells when a coating is drawn.
    auto obstacles = buildCoreNamed(magnetic.get_core(), corePolygonSegments);
    if (coreCoatingThickness > 0.0) {
        std::vector<NamedShape> coatings;
        for (auto& ns : obstacles) {
            auto shell = buildCoreCoatingShell(ns.shape, coreCoatingThickness);
            if (!shell.IsNull()) coatings.push_back({shell, ns.name + " coating", Role::CoreCoating});
        }
        for (auto& c : coatings) obstacles.push_back(std::move(c));
    }
    MAS::CoreBobbinProcessedDescription bobbinPd;
    bool toroidalCore = false;
    ConductorBuilder::Options copts = realWindingConductorOptions(
        magnetic, obstacles, wirePolygonSegments, femReady, bobbinPd, toroidalCore);
    copts.paintCoating = paintCoating;
    auto leads = ConductorBuilder::measureTerminalLeadLengths(magnetic.get_coil(), bobbinPd,
                                                              toroidalCore, copts);
    // ABT #1248: into the exported frame (lengths, radii and sweeps are invariant).
    if (toroidalCore) {
        const gp_Trsf& t = copts.toroidMounting->toExported;
        auto move = [&](std::array<double, 3>& p) {
            gp_Pnt q(p[0], p[1], p[2]);
            q.Transform(t);
            p = {q.X(), q.Y(), q.Z()};
        };
        for (auto& [winding, lead] : leads)
            for (auto& e : lead.per_end)
                for (auto& pc : e.pieces) {
                    move(pc.start);
                    move(pc.end);
                }
    }
    return leads;
}

nlohmann::json MagneticBuilder::terminalLeadLengthsToJson(
    const std::map<std::string, ConductorBuilder::TerminalLeadLength>& leads) {
    nlohmann::json j = nlohmann::json::object();
    for (const auto& [winding, t] : leads) {
        nlohmann::json ends = nlohmann::json::array();
        for (const auto& e : t.per_end) {
            nlohmann::json pieces = nlohmann::json::array();
            for (const auto& pc : e.pieces)
                pieces.push_back({{"label", pc.label},
                                  {"kind", pc.kind},
                                  {"length_m", pc.length_m},
                                  {"start_m", pc.start},
                                  {"end_m", pc.end},
                                  {"radius_m", pc.radius_m},
                                  {"sweep_rad", pc.sweep_rad}});
            ends.push_back({{"parallel", e.parallel},
                            {"end", e.end},
                            {"length_m", e.length_m},
                            {"pieces", std::move(pieces)}});
        }
        j["winding_" + winding] = {{"terminal_lead_length_m", t.total_m},
                                   {"parallels", t.parallels},
                                   {"ends", std::move(ends)}};
    }
    return j;
}

std::string MagneticBuilder::terminalLeadSidecarPath(const std::string& stepPath) {
    std::filesystem::path p(stepPath);
    p.replace_extension(".leads.json");
    return p.string();
}

std::string MagneticBuilder::writeTerminalLeadSidecar(
    const std::map<std::string, ConductorBuilder::TerminalLeadLength>& leads,
    const std::string& stepPath) {
    const std::string path = terminalLeadSidecarPath(stepPath);
    std::ofstream f(path);
    if (!f)
        throw std::runtime_error("writeTerminalLeadSidecar: cannot open '" + path + "' for writing");
    // max_digits10: the lengths round-trip bit-exactly, so a consumer's 1e-9 m comparison is
    // against the number MVB++ computed, not a printed approximation of it.
    f << std::setprecision(std::numeric_limits<double>::max_digits10)
      << terminalLeadLengthsToJson(leads).dump(2) << "\n";
    if (!f)
        throw std::runtime_error("writeTerminalLeadSidecar: write to '" + path + "' failed");
    return path;
}

// The ONE place real-winding conductors are emitted. buildAllNamed (whole assembly, one
// call) and buildRealWindingTurnsNamed (turns only, for a viewer that draws core, bobbin
// and turns as separate meshes) both come through here, so the two can never drift into
// drawing different copper for the same magnetic.
std::vector<NamedShape> MagneticBuilder::buildRealWindingConductorsNamed(
    const OpenMagnetics::Magnetic& magnetic,
    const std::vector<NamedShape>& coreShapes,
    int wirePolygonSegments,
    bool paintCoating,
    bool emitCoatingShells,
    bool femReady,
    bool diagnosticSkipCollisionCheck,
    bool cutterOnly,
    double* toroidTerminalPlaneOut) const {
    MAS::CoreBobbinProcessedDescription bobbinPd;
    bool toroidalCore = false;
    ConductorBuilder::Options copts = realWindingConductorOptions(
        magnetic, coreShapes, wirePolygonSegments, femReady, bobbinPd, toroidalCore);
    copts.diagnosticSkipCollisionCheck = diagnosticSkipCollisionCheck;
    copts.cutterOnly = cutterOnly;
    if (toroidTerminalPlaneOut) *toroidTerminalPlaneOut = std::numeric_limits<double>::quiet_NaN();

    std::vector<NamedShape> out;
    auto emitConductors = [&](bool coat, const std::string& suffix) {
        copts.paintCoating = coat;
        double plane = std::numeric_limits<double>::quiet_NaN();
        copts.toroidTerminalPlaneOut = &plane;
        auto built = ConductorBuilder::buildAll(magnetic.get_coil(), bobbinPd, toroidalCore, copts);
        copts.toroidTerminalPlaneOut = nullptr;
        if (toroidTerminalPlaneOut && !std::isnan(plane) &&
            (std::isnan(*toroidTerminalPlaneOut) || plane < *toroidTerminalPlaneOut))
            *toroidTerminalPlaneOut = plane;
        for (auto& ns : built) {
            // Carry the per-solid names through (ABT #685) — rebuilding the NamedShape from
            // {shape, name} alone silently dropped them, and the STEP went back to one unnamed
            // multi-solid product.
            // ABT #1169: the conductor builder already roled each solid (Turn / Terminal /
            // Solder); the coating pass re-draws the SAME conductors at their outer footprint,
            // so only a Turn becomes a TurnCoating — a terminal cap or a solder body does not.
            // ABT #1245: only the SHELL pass (emitCoatingShells, " coating" suffix) is a coating.
            // A conductor painted at its outer diameter (paintCoating alone, the web viewer) is
            // still THE conductor and stays a Turn, as on the per-turn path.
            const bool shellPass = !suffix.empty();
            const Role r = (shellPass && ns.role == Role::Turn) ? Role::TurnCoating : ns.role;
            out.push_back({ns.shape, ns.name + suffix, std::move(ns.partNames), r});
        }
    };
    if (emitCoatingShells) {
        emitConductors(false, "");           // bare copper conductor
        emitConductors(true, " coating");    // outer insulated footprint
    } else {
        emitConductors(paintCoating, "");
    }
    return out;
}

std::vector<NamedShape> MagneticBuilder::buildRealWindingTurnsNamed(
    const OpenMagnetics::Magnetic& magnetic,
    int wirePolygonSegments,
    int corePolygonSegments,
    bool paintCoating,
    bool femReady,
    bool diagnosticSkipCollisionCheck) const {
    // The core is built but NOT returned: the conductor builder needs the real core solids
    // to aim the terminal leads at the true window opening. Skipping them here would make a
    // viewer's leads land somewhere else than the exported assembly's for the same design —
    // the kind of divergence nobody notices until a STEP and a screenshot disagree.
    auto coreShapes = buildCoreNamed(magnetic.get_core(), corePolygonSegments);
    return buildRealWindingConductorsNamed(magnetic, coreShapes, wirePolygonSegments,
                                           paintCoating, /*emitCoatingShells=*/false, femReady,
                                           diagnosticSkipCollisionCheck);
}

std::vector<NamedShape> MagneticBuilder::buildTurnsNamedFromTurns(
    const std::vector<MAS::Turn>& turns,
    int wirePolygonSegments,
    bool paintCoating) const {
    std::vector<NamedShape> out;
    out.reserve(turns.size());
    TurnBuilder::clearCache();
    for (std::size_t i = 0; i < turns.size(); ++i) {
        TopoDS_Shape s = TurnBuilder::buildFromTurnAlone(turns[i], wirePolygonSegments, paintCoating);
        if (s.IsNull()) {
            throw std::runtime_error(
                "buildTurnsNamedFromTurns: TurnBuilder produced null shape for turn "
                + std::to_string(i));
        }
        const std::string& n = turns[i].get_name();
        out.emplace_back(s, n.empty() ? ("Turn_" + std::to_string(i)) : n, Role::Turn);
    }
    return out;
}

std::vector<ConductorBuilder::PathPolyline> MagneticBuilder::buildRealWindingPaths(
    const OpenMagnetics::Magnetic& magnetic) const {
    auto bobbinPd = getBobbinProcessed(magnetic.get_coil());
    patchBobbinDimensions(bobbinPd, magnetic.get_core());
    const bool toroidalCore = isCoreToroidal(magnetic.get_core());
    ConductorBuilder::Options copts;
    copts.femReady = true;
    copts.paintCoating = false;   // centrelines + copper radius; coating handled by the consumer
    auto cores = buildCoreNamed(magnetic.get_core(), DEFAULT_CORE_POLYGON_SEGMENTS);
    if (!toroidalCore) {
        for (const auto& ns : cores) copts.coreObstacles.push_back(ns.shape);
    }
    else {   // ABT #1248
        copts.toroidMounting = toroidMountingFrameOf(magnetic);
        copts.toroidCoreDepthAlongDown = coreDepthAlongDown(cores, copts.toroidMounting->down);
    }
    copts.woundColumnPerSection =
        resolveWoundColumnsPerSection(magnetic.get_coil(), magnetic.get_core(), bobbinPd);
    auto paths = ConductorBuilder::buildAllPaths(magnetic.get_coil(), bobbinPd, toroidalCore, copts);
    if (toroidalCore) {   // ABT #1248: into the exported frame
        const gp_Trsf& t = copts.toroidMounting->toExported;
        for (auto& pl : paths) {
            auto movePoint = [&](std::array<double, 3>& p) {
                gp_Pnt q(p[0], p[1], p[2]);
                q.Transform(t);
                p = {q.X(), q.Y(), q.Z()};
            };
            auto moveDir = [&](std::array<double, 3>& d) {
                gp_Vec v(d[0], d[1], d[2]);
                v.Transform(t);
                d = {v.X(), v.Y(), v.Z()};
            };
            for (auto& prim : pl.prims)
                for (auto& q : prim) movePoint(q);
            movePoint(pl.end0);
            movePoint(pl.end1);
            moveDir(pl.dir0);
            moveDir(pl.dir1);
        }
    }
    return paths;
}

ConductorBuilder::ToroidMounting MagneticBuilder::toroidMountingOf(
    const OpenMagnetics::Magnetic& magnetic) {
    if (!isCoreToroidal(magnetic.get_core()))
        throw std::runtime_error("toroidMountingOf: the core is not a toroid");
    return toroidMountingT<OpenMagnetics::Bobbin>(magnetic.get_coil().get_bobbin());
}

ConductorBuilder::ToroidMounting MagneticBuilder::toroidMountingOf(const MAS::Magnetic& magnetic) {
    if (!magnetic.get_core() || !magnetic.get_coil())
        throw std::runtime_error("toroidMountingOf: the magnetic needs a core and a coil");
    if (!isCoreToroidal(*magnetic.get_core()))
        throw std::runtime_error("toroidMountingOf: the core is not a toroid");
    return toroidMountingT<MAS::Bobbin>(magnetic.get_coil()->get_bobbin());
}

template<typename FunctionalDescriptionT>
static std::vector<std::string> windingOrderOf(const FunctionalDescriptionT& fd) {
    std::vector<std::string> order;
    for (const auto& w : fd) order.push_back(w.get_name());
    return order;
}

ConductorBuilder::ToroidMountingFrame MagneticBuilder::toroidMountingFrameOf(
    const OpenMagnetics::Magnetic& magnetic) {
    const auto mounting = toroidMountingOf(magnetic);
    const auto& coil = magnetic.get_coil();
    const auto turns = coil.get_turns_description();
    return ConductorBuilder::resolveToroidMountingFrame(
        mounting, turns ? *turns : std::vector<MAS::Turn>{},
        windingOrderOf(coil.get_functional_description()));
}

ConductorBuilder::ToroidMountingFrame MagneticBuilder::toroidMountingFrameOf(
    const MAS::Magnetic& magnetic) {
    const auto mounting = toroidMountingOf(magnetic);
    const MAS::Coil coil = magnetic.get_coil().value();
    const auto turns = coil.get_turns_description();
    return ConductorBuilder::resolveToroidMountingFrame(
        mounting, turns ? *turns : std::vector<MAS::Turn>{},
        windingOrderOf(coil.get_functional_description()));
}

} // namespace mvb
