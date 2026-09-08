#include "mvb/MagneticBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Symmetry.h"
#include "mvb/Utils.h"
#include "mvb/shapes/ShapeBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/ConductorBuilder.h"
#include "mvb/BobbinBuilder.h"
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
#include <stdexcept>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepOffset_Mode.hxx>
#include <GeomAbs_JoinType.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Precision.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <Geom_Surface.hxx>
#include <BRep_Tool.hxx>
#include <Standard_Failure.hxx>
#include <iostream>
#include <BRepBndLib.hxx>
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
                                      double scale) {
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
    exportSTEP(scaled, out.string());
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
    return drawMagneticCommon(named, outputPath, cfg.format, cfg.scale);
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
    return drawMagneticCommon(named, outputPath, cfg.format, cfg.scale);
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

        // Apply machining after rotation. Wrap each cut so that an
        // OCCT/std failure on one column doesn't abort the whole core —
        // the un-gapped shape is a strictly-better fallback than no shape.
        auto machiningOpt = piece.get_machining();
        if (machiningOpt) {
            for (const auto& mach : *machiningOpt) {
                shape = builder->applyMachining(shape, mach, dims);
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
    return BobbinBuilder::buildBobbin(bobbinPd, flangeThickness, !isCoreToroidal(core), polygonSegments);
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
        if (!s.IsNull()) out.push_back({s, "insulation_layer_" + std::to_string(i)});
    }
    return out;
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
                         single ? base : base + "_" + std::to_string(i));
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
        out.emplace_back(shapes[i], n);
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
        out.emplace_back(shapes[i], n);
    }
    return out;
}

NamedShape MagneticBuilder::buildBobbinNamed(const MAS::Coil& coil,
                                             const MAS::MagneticCore& core,
                                             int corePolygonSegments) const {
    NamedShape ns;
    ns.shape = buildBobbinShape_impl(coil, core, corePolygonSegments);
    ns.name = getBobbinNameT<MAS::Bobbin>(coil.get_bobbin(), "Bobbin");
    return ns;
}

NamedShape MagneticBuilder::buildBobbinNamed(const OpenMagnetics::Coil& coil,
                                             const MAS::MagneticCore& core,
                                             int corePolygonSegments) const {
    NamedShape ns;
    ns.shape = buildBobbinShape_impl(coil, core, corePolygonSegments);
    ns.name = getBobbinNameT<OpenMagnetics::Bobbin>(coil.get_bobbin(), "Bobbin");
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
        out.emplace_back(shapes[i], n);
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
                                                         bool femReady) const {
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
                if (!shell.IsNull()) coatings.push_back({shell, ns.name + " coating"});
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
                bobbin.shape = cut_bobbin(bobbin.shape, cutters);
                if (!bobbin.shape.IsNull()) all.push_back(bobbin);
            }
        }

        for (std::size_t i = 0; i < turnShapes.size(); ++i) {
            const std::string n = (i < turnNames.size() && !turnNames[i].empty())
                                    ? turnNames[i]
                                    : "Turn_" + std::to_string(i);
            all.emplace_back(turnShapes[i], n);
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
                all.emplace_back(fr4, "FR4Board");
            }
        }

        return apply_symmetry(std::move(all), symmetryPlanes);
    }

    OpenMagnetics::Magnetic enriched = magnetic_autocomplete_safe(magnetic, useRealWindingGeometry);
    return buildAllNamed(enriched, includeBobbin, symmetryPlanes,
                         wirePolygonSegments, corePolygonSegments, paintCoating, emitCoatingShells,
                         includeInsulation, coreCoatingThickness, useRealWindingGeometry, femReady);
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
                                                         bool femReady) const {
    auto all = buildCoreNamed(magnetic.get_core(), corePolygonSegments);

    if (coreCoatingThickness > 0.0) {   // conformal core-coating shells (offset core - core)
        std::vector<NamedShape> coatings;
        for (auto& ns : all) {
            auto shell = buildCoreCoatingShell(ns.shape, coreCoatingThickness);
            if (!shell.IsNull()) coatings.push_back({shell, ns.name + " coating"});
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
    if (useRealWindingGeometry) {
        for (auto& ns : buildRealWindingConductorsNamed(magnetic, all, wirePolygonSegments,
                                                        paintCoating, emitCoatingShells, femReady)) {
            turnShapes.push_back(ns.shape);
            turnNames.push_back(ns.name);
            turnPartNames.push_back(std::move(ns.partNames));
        }
    } else {
        turnShapes = buildTurnsImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(
            magnetic.get_coil(), magnetic.get_core(), &turnNames, wirePolygonSegments, paintCoating, emitCoatingShells);
    }

    if (includeBobbin) {
        auto bobbin = buildBobbinNamed(magnetic.get_coil(), magnetic.get_core(), corePolygonSegments);
        if (!bobbin.shape.IsNull()) {
            std::vector<TopoDS_Shape> cutters;
            for (const auto& ns : all) cutters.push_back(ns.shape);
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
                         /*femReady=*/false)) {
                    cutters.push_back(ns.shape);
                }
            } else {
                cutters.insert(cutters.end(), turnShapes.begin(), turnShapes.end());
            }
            bobbin.shape = cut_bobbin(bobbin.shape, cutters);
            if (!bobbin.shape.IsNull()) all.push_back(bobbin);
        }
    }

    // Append the already-built turns (no second build).
    for (std::size_t i = 0; i < turnShapes.size(); ++i) {
        const std::string n = (i < turnNames.size() && !turnNames[i].empty())
                                ? turnNames[i]
                                : "Turn_" + std::to_string(i);
        if (i < turnPartNames.size() && !turnPartNames[i].empty()) {
            all.emplace_back(turnShapes[i], n, turnPartNames[i]);   // ABT #685
        }
        else {
            all.emplace_back(turnShapes[i], n);
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
            all.emplace_back(fr4, "FR4Board");
        }
    }

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
    if (isToroidal)
        for (auto& ns : all) ns.shape = rotate_shape(ns.shape, -std::numbers::pi / 2.0, 0.0, 0.0);

    // FEM product only: re-express every surface as a B-spline (ABT #490 class, 2026-08-25).
    //
    // gmsh decides whether to send a face to its markedly more fragile periodic mesher by
    // reading the UNDERLYING SURFACE, not the face -- OCCFace.cpp:142,
    //   _periodic[0] = surface.IsUPeriodic();
    // so a 60-degree panel of a cylinder, or a 20-degree slice of a torus, is treated exactly
    // like a full closed one and routed to meshGeneratorPeriodic. That path produced 16 of the
    // 23 corpus mesh failures measured on 2026-08-24 ("Impossible to mesh periodic surface" x6
    // and "The 1D mesh seems not to be forming a closed loop" x10, both raised inside it).
    // Round wire is nothing but cylinders, tori and swept pipes, so nearly every conductor face
    // took it: 3938 cylinder + 1342 torus faces across the 30-design corpus.
    //
    // A conic IS exactly representable as a rational B-spline, so this conversion MOVES NO
    // POINT -- it is not a tolerance-based approximation and not a geometry change. It simply
    // drops the periodicity flag on faces that do not actually wrap, which is every conductor
    // face once the section is faceted. Measured: 03_buck 20 -> 0 periodic surfaces (226 faces
    // unchanged), 05_pfc toroid 2668 -> 0 (7288 faces unchanged), and both stay BRepCheck-valid.
    //
    // Deliberately NOT applied to the drawing product: exact analytic surfaces are smaller and
    // render better, and a viewer never meets gmsh. MVB_NO_NURBS=1 disables it for comparison.
    if (femReady && !std::getenv("MVB_NO_NURBS")) {
        for (auto& ns : all) {
            if (ns.shape.IsNull()) continue;
            // ONLY SOLIDS THAT ACTUALLY CARRY A PERIODIC FACE -- per SOLID, not per part.
            // NurbsConvert rewrites EVERY surface, so run unconditionally it re-expressed 9579
            // PLANAR faces of the faceted 14_dab as B-splines -- planes are never periodic, so
            // that bought nothing and made every downstream boolean heavier and hungrier (a
            // faceted design went 29 s -> 608 s). Per-PART gating was still too coarse: one
            // periodic lead pipe condemned every all-planar prism in the conductor to bspline
            // re-expression, and two abutting prisms' independently converted caps are what the
            // fragment then imprints as nm slivers and folded facets (10_emi's trimmed riser
            // PLC-failed ALONE; 12_boost's 32 fragment micro-edges). A planar solid stays
            // analytic: its junction caps remain the same exact planes as its neighbour's.
            auto hasPeriodicFace = [](const TopoDS_Shape& s) {
                for (TopExp_Explorer fx(s, TopAbs_FACE); fx.More(); fx.Next()) {
                    Handle(Geom_Surface) srf = BRep_Tool::Surface(TopoDS::Face(fx.Current()));
                    if (!srf.IsNull() && (srf->IsUPeriodic() || srf->IsVPeriodic())) return true;
                }
                return false;
            };
            if (!hasPeriodicFace(ns.shape)) continue;
            const bool wasValid = BRepCheck_Analyzer(ns.shape).IsValid();
            try {
                // Rebuild the part solid-by-solid, converting only the periodic-carrying ones.
                TopoDS_Shape out;
                if (ns.shape.ShapeType() == TopAbs_COMPOUND) {
                    TopoDS_Compound rebuilt;
                    BRep_Builder rb;
                    rb.MakeCompound(rebuilt);
                    for (TopoDS_Iterator it(ns.shape); it.More(); it.Next()) {
                        const TopoDS_Shape& child = it.Value();
                        if (hasPeriodicFace(child)) {
                            BRepBuilderAPI_NurbsConvert cc(child, /*Copy=*/Standard_True);
                            rb.Add(rebuilt, cc.Shape().IsNull() ? child : cc.Shape());
                        } else {
                            rb.Add(rebuilt, child);
                        }
                    }
                    out = rebuilt;
                } else {
                    BRepBuilderAPI_NurbsConvert conv(ns.shape, /*Copy=*/Standard_True);
                    out = conv.Shape();
                }
                if (out.IsNull()) continue;
                // Conversion is exact, so both volume and validity must survive it. Anything
                // else is OCC damage, not a representation change -- keep the original and say
                // so rather than shipping a silently degraded solid.
                // ADAPTIVE integration on BOTH sides. The default VolumeProperties integrates an
                // analytic face exactly and a B-spline face by fixed-order quadrature, so
                // comparing the two measures the INTEGRATION SCHEME, not the geometry: every part
                // came back ~0.17% "changed" (Core_0 84412.6 -> 84555.4 mm3) and was reverted,
                // which silently disabled this whole pass. The Eps overload refines until the
                // requested relative accuracy is met, so both numbers describe the same solid.
                // Judge the difference against the INTEGRATOR'S OWN reported accuracy, not a
                // chosen threshold. VolumeProperties(.., Eps) returns the relative error it
                // actually achieved; a difference inside the sum of the two reported errors says
                // nothing about the geometry, only about the quadrature. This is what makes the
                // test honest on both sides: analytic faces integrate almost exactly while
                // B-spline faces do not, so a fixed tolerance either reverts everything (1e-6 did)
                // or waves real changes through.
                GProp_GProps g0, g1;
                const double eps = 1e-7;
                const double err0 = std::abs(BRepGProp::VolumeProperties(ns.shape, g0, eps));
                const double err1 = std::abs(BRepGProp::VolumeProperties(out, g1, eps));
                const double v0 = std::abs(g0.Mass());
                // Two uncertainties, both physical, whichever is larger. (1) the integrators'
                // own reported error; (2) what the B-Rep can even express: every face is located
                // only to within the shape tolerance, so the volume is defined only to about
                // tolerance x surface area. A difference below that is not a geometry change, it
                // is the model's own resolution.
                GProp_GProps sp;
                BRepGProp::SurfaceProperties(ns.shape, sp);
                const double brepTol = Precision::Confusion() * std::abs(sp.Mass());
                const double noise = std::max((err0 + err1) * v0, brepTol);
                if (v0 > 0 && std::abs(g1.Mass() - v0) > noise) {
                    std::cerr << "[nurbs] '" << ns.name << "' volume changed "
                              << g0.Mass() * 1e9 << " -> " << g1.Mass() * 1e9
                              << " mm3 (beyond the model's own " << noise * 1e9
                              << " mm3); keeping the analytic surfaces\n";
                    continue;
                }
                if (wasValid && !BRepCheck_Analyzer(out).IsValid()) {
                    std::cerr << "[nurbs] '" << ns.name
                              << "' conversion produced an invalid solid; keeping the analytic"
                                 " surfaces\n";
                    continue;
                }
                ns.shape = out;
            } catch (const Standard_Failure& e) {
                std::cerr << "[nurbs] '" << ns.name << "' conversion threw ("
                          << e.GetMessageString() << "); keeping the analytic surfaces\n";
            }
        }
    }

    return apply_symmetry(std::move(all), symmetryPlanes);
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
    return NamedShape{geom, name};
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
        return NamedShape{TopoDS_Shape(), name};
    }
    return NamedShape{s, name};
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
    bool diagnosticSkipCollisionCheck) const {
    auto bobbinPd = getBobbinProcessed(magnetic.get_coil());
    patchBobbinDimensions(bobbinPd, magnetic.get_core());
    const bool toroidalCore = isCoreToroidal(magnetic.get_core());
    ConductorBuilder::Options copts;
    copts.wirePolygonSegments = wirePolygonSegments;
    copts.femReady = femReady;   // OM drawing -> fast compound; FEM export -> one-piece/conformal
    copts.diagnosticSkipCollisionCheck = diagnosticSkipCollisionCheck;
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

    std::vector<NamedShape> out;
    auto emitConductors = [&](bool coat, const std::string& suffix) {
        copts.paintCoating = coat;
        for (auto& ns : ConductorBuilder::buildAll(magnetic.get_coil(), bobbinPd,
                                                   toroidalCore, copts)) {
            // Carry the per-solid names through (ABT #685) — rebuilding the NamedShape from
            // {shape, name} alone silently dropped them, and the STEP went back to one unnamed
            // multi-solid product.
            out.push_back({ns.shape, ns.name + suffix, std::move(ns.partNames)});
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
        out.emplace_back(s, n.empty() ? ("Turn_" + std::to_string(i)) : n);
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
    if (!toroidalCore) {
        auto cores = buildCoreNamed(magnetic.get_core(), DEFAULT_CORE_POLYGON_SEGMENTS);
        for (const auto& ns : cores) copts.coreObstacles.push_back(ns.shape);
    }
    copts.woundColumnPerSection =
        resolveWoundColumnsPerSection(magnetic.get_coil(), magnetic.get_core(), bobbinPd);
    return ConductorBuilder::buildAllPaths(magnetic.get_coil(), bobbinPd, toroidalCore, copts);
}

} // namespace mvb
