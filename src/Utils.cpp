#include "mvb/Utils.h"
#include <algorithm>
#include <vector>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BOPAlgo_GlueEnum.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <gp_Pnt.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <variant>
#include <stdexcept>
#include <sstream>
#include <numbers>
#include <map>

// MKF headers for safe autocomplete
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"

namespace mvb {

double flatten_dimension(const MAS::Dimension& dim) {
    if (std::holds_alternative<MAS::DimensionWithTolerance>(dim)) {
        const auto& dwt = std::get<MAS::DimensionWithTolerance>(dim);
        if (dwt.get_nominal()) {
            return *dwt.get_nominal();
        }
        if (dwt.get_maximum().has_value() && dwt.get_minimum().has_value()) {
            return (*dwt.get_maximum() + *dwt.get_minimum()) / 2.0;
        }
        if (dwt.get_maximum().has_value()) {
            return *dwt.get_maximum();
        }
        if (dwt.get_minimum().has_value()) {
            return *dwt.get_minimum();
        }
        throw std::runtime_error("DimensionWithTolerance missing nominal value");
    }
    return std::get<double>(dim);
}

std::map<std::string, double> flatten_dimensions(const std::map<std::string, MAS::Dimension>& dims) {
    std::map<std::string, double> result;
    for (const auto& [key, dim] : dims) {
        result[key] = flatten_dimension(dim);
    }
    return result;
}

TopoDS_Wire build_polygon_circle(double radius, int segments) {
    if (segments <= 0) {
        gp_Circ circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), radius);
        BRepBuilderAPI_MakeEdge edge(circ);
        BRepBuilderAPI_MakeWire wire(edge);
        return wire.Wire();
    }
    BRepBuilderAPI_MakePolygon poly;
    for (int i = 0; i < segments; ++i) {
        double angle = 2.0 * std::numbers::pi * i / segments;
        poly.Add(gp_Pnt(radius * std::cos(angle), radius * std::sin(angle), 0.0));
    }
    poly.Close();
    return BRepBuilderAPI_MakeWire(poly.Wire()).Wire();
}

TopoDS_Shape build_polygon_cylinder(double height, double radius, int segments) {
    if (segments <= 0) {
        return BRepPrimAPI_MakeCylinder(radius, height).Shape();
    }
    TopoDS_Wire wire = build_polygon_circle(radius, segments);
    BRepBuilderAPI_MakeFace face(wire);
    gp_Vec vec(0, 0, height);
    return BRepPrimAPI_MakePrism(face.Face(), vec).Shape();
}

TopoDS_Shape build_polygon_ring(double turn_radius, double wire_radius,
                                 double y, int cross_segments,
                                 int revolution_segments) {
    // Fallback: analytic torus via MakeRevol of a polygonal cross-section.
    // Produces CYLINDRICAL/TOROIDAL_SURFACE in STEP; useful if you need the
    // old behaviour.
    if (revolution_segments <= 0) {
        gp_Ax2 profile_plane(gp_Pnt(turn_radius, y, 0.0),
                             gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        TopoDS_Wire prof_wire = build_polygon_circle(wire_radius, cross_segments);
        BRepBuilderAPI_Transform xform(
            prof_wire,
            gp_Trsf()  // identity — then rotate + translate
        );
        // Place the polygon face in `profile_plane`:
        gp_Trsf t;
        gp_Ax3 from(gp_Pnt(0,0,0), gp_Dir(0,0,1), gp_Dir(1,0,0));
        gp_Ax3 to(profile_plane);
        t.SetDisplacement(from, to);
        BRepBuilderAPI_Transform placer(prof_wire, t);
        TopoDS_Wire placed_wire = TopoDS::Wire(placer.Shape());
        BRepBuilderAPI_MakeFace face(placed_wire);
        gp_Ax1 rev_axis(gp_Pnt(0.0, y, 0.0), gp_Dir(0, 1, 0));
        return BRepPrimAPI_MakeRevol(face.Face(), rev_axis,
                                      2.0 * std::numbers::pi - 1e-6).Shape();
    }

    // ThruSections: place `revolution_segments` copies of the polygonal
    // cross-section around the Y axis at angles 0, 2π/N, ..., 2π — the
    // loft between consecutive sections uses RULED (planar) faces.
    BRepOffsetAPI_ThruSections gen(/*IsSolid*/ Standard_True, /*IsRuled*/ Standard_True);
    gen.SetSmoothing(Standard_False);

    TopoDS_Wire base = build_polygon_circle(wire_radius, cross_segments);
    // Place base at (turn_radius, y, 0), with the face's local Z along the
    // azimuthal (φ̂) direction at that angle — so rotating around Y sweeps
    // the polygon's plane through the ring.
    for (int i = 0; i <= revolution_segments; ++i) {
        double angle = 2.0 * std::numbers::pi * static_cast<double>(i)
                       / static_cast<double>(revolution_segments);
        double cx = turn_radius * std::cos(angle);
        double cz = turn_radius * std::sin(angle);
        // Plane at this station: normal = tangent (azimuthal), X axis = radial
        gp_Dir tangent(-std::sin(angle), 0.0, std::cos(angle));
        gp_Dir radial(std::cos(angle), 0.0, std::sin(angle));
        gp_Ax3 from(gp_Pnt(0,0,0), gp_Dir(0,0,1), gp_Dir(1,0,0));
        gp_Ax3 to(gp_Pnt(cx, y, cz), tangent, radial);
        gp_Trsf t;
        t.SetDisplacement(from, to);
        BRepBuilderAPI_Transform placer(base, t);
        TopoDS_Wire placed = TopoDS::Wire(placer.Shape());
        gen.AddWire(placed);
    }
    gen.Build();
    if (!gen.IsDone()) return TopoDS_Shape();
    return gen.Shape();
}

TopoDS_Shape rotate_shape(const TopoDS_Shape& shape, double rx, double ry, double rz) {
    TopoDS_Shape result = shape;
    if (rx != 0.0) {
        gp_Trsf trsf;
        trsf.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0)), rx);
        result = BRepBuilderAPI_Transform(result, trsf).Shape();
    }
    if (ry != 0.0) {
        gp_Trsf trsf;
        trsf.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,1,0)), ry);
        result = BRepBuilderAPI_Transform(result, trsf).Shape();
    }
    if (rz != 0.0) {
        gp_Trsf trsf;
        trsf.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,0,1)), rz);
        result = BRepBuilderAPI_Transform(result, trsf).Shape();
    }
    return result;
}

gp_Trsf translation_trsf(double x, double y, double z) {
    gp_Trsf trsf;
    trsf.SetTranslation(gp_Vec(x, y, z));
    return trsf;
}

TopoDS_Shape translate_shape(const TopoDS_Shape& shape, double x, double y, double z) {
    return BRepBuilderAPI_Transform(shape, translation_trsf(x, y, z)).Shape();
}

void patch_dimension_nominals(nlohmann::json& j) {
    if (j.is_object()) {
        // If this looks like a dimension object with min/max but no nominal, add one
        bool hasMin = j.contains("minimum") && !j["minimum"].is_null();
        bool hasMax = j.contains("maximum") && !j["maximum"].is_null();
        if (hasMin || hasMax) {
            if (!j.contains("nominal") || j["nominal"].is_null()) {
                if (hasMin && hasMax) {
                    j["nominal"] = (j["minimum"].get<double>() + j["maximum"].get<double>()) / 2.0;
                } else if (hasMax) {
                    j["nominal"] = j["maximum"].get<double>();
                } else {
                    j["nominal"] = j["minimum"].get<double>();
                }
            }
        }
        for (auto& [key, val] : j.items()) {
            patch_dimension_nominals(val);
        }
    } else if (j.is_array()) {
        for (auto& elem : j) {
            patch_dimension_nominals(elem);
        }
    }
}

std::string core_shape_family_to_string(MAS::CoreShapeFamily family) {
    switch (family) {
        case MAS::CoreShapeFamily::C: return "c";
        case MAS::CoreShapeFamily::DRUM: return "drum";
        case MAS::CoreShapeFamily::E: return "e";
        case MAS::CoreShapeFamily::EC: return "ec";
        case MAS::CoreShapeFamily::EFD: return "efd";
        case MAS::CoreShapeFamily::EI: return "ei";
        case MAS::CoreShapeFamily::EL: return "el";
        case MAS::CoreShapeFamily::ELP: return "elp";
        case MAS::CoreShapeFamily::EP: return "ep";
        case MAS::CoreShapeFamily::EPX: return "epx";
        case MAS::CoreShapeFamily::EQ: return "eq";
        case MAS::CoreShapeFamily::ER: return "er";
        case MAS::CoreShapeFamily::ETD: return "etd";
        case MAS::CoreShapeFamily::H: return "h";
        case MAS::CoreShapeFamily::LP: return "lp";
        case MAS::CoreShapeFamily::P: return "p";
        case MAS::CoreShapeFamily::PM: return "pm";
        case MAS::CoreShapeFamily::PQ: return "pq";
        case MAS::CoreShapeFamily::RM: return "rm";
        case MAS::CoreShapeFamily::ROD: return "rod";
        case MAS::CoreShapeFamily::T: return "t";
        case MAS::CoreShapeFamily::U: return "u";
        case MAS::CoreShapeFamily::UT: return "ut";
    }
    return "unknown";
}

OpenMagnetics::Magnetic magnetic_autocomplete_safe(const nlohmann::json& magneticJson) {
    using json = nlohmann::json;

    json coreJson = magneticJson.contains("core") ? magneticJson.at("core") : json::object();

    // Construct OpenMagnetics::Core and Coil directly from JSON.
    // The key fix: pass false to Coil constructor to skip wind(),
    // which crashes for raw MAS files with "Basic" bobbins.
    // Also: if the "coil" key is missing or empty, use a default Coil —
    // Coil(empty_json, false) throws "key 'bobbin' not found".
    OpenMagnetics::Core core;
    try {
        core = OpenMagnetics::Core(coreJson);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("magnetic_autocomplete_safe: Core ctor failed: ") + e.what());
    }
    OpenMagnetics::Coil coil;
    if (magneticJson.contains("coil") && !magneticJson.at("coil").is_null()
        && !magneticJson.at("coil").empty()) {
        try {
            coil = OpenMagnetics::Coil(magneticJson.at("coil"), false);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("magnetic_autocomplete_safe: Coil ctor failed: ") + e.what());
        }
    }

    OpenMagnetics::Magnetic om;
    om.set_core(core);
    om.set_coil(coil);

    if (magneticJson.contains("distributorsInfo")) {
        om.set_distributors_info(magneticJson.at("distributorsInfo").get<std::optional<std::vector<MAS::DistributorInfo>>>());
    }
    if (magneticJson.contains("manufacturerInfo")) {
        om.set_manufacturer_info(magneticJson.at("manufacturerInfo").get<std::optional<MAS::MagneticManufacturerInfo>>());
    }
    if (magneticJson.contains("rotation")) {
        om.set_rotation(magneticJson.at("rotation").get<std::optional<std::vector<double>>>());
    }

    // If the coil has no windings (core-only visualizer), calling the full
    // OpenMagnetics::magnetic_autocomplete crashes because it calls get_wire(0)
    // on an empty functional_description vector. Enrich the core directly instead.
    if (om.get_coil().get_functional_description().empty()) {
        auto& mutCore = om.get_mutable_core();

        // Guarantee required fields that MKF unconditionally dereferences.
        if (!mutCore.get_functional_description().get_number_stacks()) {
            mutCore.get_mutable_functional_description().set_number_stacks(int64_t(1));
        }

        mutCore.resolve_shape();

        auto shapeFamily = mutCore.get_shape_family();
        if (shapeFamily == MAS::CoreShapeFamily::T) {
            mutCore.get_mutable_functional_description().set_type(MAS::CoreType::TOROIDAL);
        } else {
            mutCore.get_mutable_functional_description().set_type(MAS::CoreType::TWO_PIECE_SET);
        }

        mutCore.resolve_material();

        if (!mutCore.get_processed_description()) {
            mutCore.process_data();
            mutCore.process_gap();
        }

        auto geoDesc = mutCore.create_geometrical_description();
        mutCore.set_geometrical_description(geoDesc);

        return om;
    }

    auto enriched = OpenMagnetics::magnetic_autocomplete(om, json{});

    // No-fallbacks rule: if the design has windings but MKF produced zero
    // turnsDescription (typically because the winding does not fit the
    // bobbin window), surface that as a hard failure instead of returning
    // a magnetic with an empty turns list that silently propagates as a
    // zero-turn build.
    {
        const auto& coilE = enriched.get_coil();
        const auto& fd    = coilE.get_functional_description();
        if (!fd.empty()) {
            int64_t expected = 0;
            for (const auto& w : fd) {
                expected += w.get_number_turns() * w.get_number_parallels();
            }
            const auto& td = coilE.get_turns_description();
            const size_t got = td ? td->size() : 0;
            if (expected > 0 && got == 0) {
                std::ostringstream s;
                s << "magnetic_autocomplete_safe: MKF produced 0 turns for a coil "
                     "that declares " << expected
                  << " turns total — the winding likely does not fit the bobbin "
                     "window. Refusing to return a magnetic with an empty "
                     "turnsDescription.";
                throw std::runtime_error(s.str());
            }
        }
    }

    return enriched;
}

OpenMagnetics::Magnetic magnetic_autocomplete_safe(const MAS::Magnetic& magnetic) {
    json j;
    to_json(j, magnetic);
    return magnetic_autocomplete_safe(j);
}

bool is_shape_usable(const TopoDS_Shape& shape) {
    if (shape.IsNull()) return false;
    Bnd_Box bb;
    BRepBndLib::Add(shape, bb);
    return !bb.IsVoid();
}

TopoDS_Shape cut_bobbin(const TopoDS_Shape& bobbin, const std::vector<TopoDS_Shape>& cutters) {
    if (std::getenv("MVB_NO_BOBBIN_CUT")) return bobbin;
    if (bobbin.IsNull()) return bobbin;

    // Scale into mm for the cut so coordinate magnitudes are well above
    // OCCT's 1e-7 numerical tolerance floor, then scale back to meters.
    const double S = 1000.0;
    gp_Trsf up;   up.SetScale(gp_Pnt(0, 0, 0), S);
    gp_Trsf down; down.SetScale(gp_Pnt(0, 0, 0), 1.0 / S);

    // Pre-scale tools once; reuse for every per-solid cut.
    struct ScaledTool { TopoDS_Shape shape; Bnd_Box box; };
    std::vector<ScaledTool> scaledTools;
    scaledTools.reserve(cutters.size());
    for (const auto& tool : cutters) {
        if (tool.IsNull()) continue;
        TopoDS_Shape s = BRepBuilderAPI_Transform(tool, up).Shape();
        Bnd_Box tb;
        BRepBndLib::Add(s, tb);
        if (tb.IsVoid()) continue;
        scaledTools.push_back({s, tb});
    }
    if (scaledTools.empty()) return bobbin;

    // Iterate solids in the bobbin. BobbinBuilder returns a compound of
    // body + topFlange + bottomFlange (3 solids) when flanges are present.
    // Per-solid AABB prefiltering means flange cuts only see the few turns
    // whose AABB actually pokes into the flange volume — typically 0 for a
    // tall winding. Drops `Cut` tool count from O(N) total to O(few) per
    // solid and is the difference between a multi-second BOP on a NURBS
    // bobbin against 50 polygonal turns and a near-instant compound result.
    std::vector<TopoDS_Shape> solids;
    for (TopExp_Explorer exp(bobbin, TopAbs_SOLID); exp.More(); exp.Next()) {
        solids.push_back(exp.Current());
    }
    if (solids.empty()) {
        // Bobbin isn't decomposable into solids — fall back to original
        // single-shape cut.
        solids.push_back(bobbin);
    }

    TopoDS_Compound resultComp;
    BRep_Builder cbld;
    cbld.MakeCompound(resultComp);

    // Chunk size for the BOP cut. Feeding more than ~10 tools to a single
    // BRepAlgoAPI_Cut blows the WASM heap (OCCT internally builds a full
    // bbox-pair table and intermediate shapes). 80-turn bobbin cuts then
    // crash with `Aborted()` mid-Build. Iterate `current <- current ∖ batch`
    // in fixed-size batches so peak BOP memory stays bounded regardless of
    // total tool count. Native (non-WASM) builds also benefit (linear vs
    // ~quadratic-on-tools growth in OCCT's internal data structures).
    constexpr std::size_t BOP_TOOL_CHUNK = 8;

    for (const auto& solid : solids) {
        TopoDS_Shape solidScaled = BRepBuilderAPI_Transform(solid, up).Shape();
        Bnd_Box solidBox;
        BRepBndLib::Add(solidScaled, solidBox);

        std::vector<TopoDS_Shape> relevant;
        relevant.reserve(scaledTools.size());
        for (const auto& st : scaledTools) {
            if (solidBox.IsOut(st.box)) continue;
            relevant.push_back(st.shape);
        }
        if (relevant.empty()) {
            cbld.Add(resultComp, solid);
            continue;
        }

        TopoDS_Shape current = solidScaled;
        bool ok = true;
        for (std::size_t i = 0; i < relevant.size(); i += BOP_TOOL_CHUNK) {
            TopTools_ListOfShape args, batchTools;
            args.Append(current);
            const std::size_t end = std::min(i + BOP_TOOL_CHUNK, relevant.size());
            for (std::size_t j = i; j < end; ++j) {
                batchTools.Append(relevant[j]);
            }

            BRepAlgoAPI_Cut cutter;
            cutter.SetArguments(args);
            cutter.SetTools(batchTools);
            cutter.SetUseOBB(true);
            cutter.SetRunParallel(true);
            cutter.SetGlue(BOPAlgo_GlueShift);
            cutter.SetCheckInverted(false);
            cutter.SetNonDestructive(false);
            cutter.Build();
            if (!cutter.IsDone()) { ok = false; break; }
            TopoDS_Shape step = cutter.Shape();
            if (step.IsNull())   { ok = false; break; }
            current = step;
        }
        if (!ok) {
            cbld.Add(resultComp, solid);
            continue;
        }

        bool inverted = false;
        for (TopExp_Explorer exp(current, TopAbs_SOLID); exp.More(); exp.Next()) {
            GProp_GProps props;
            BRepGProp::VolumeProperties(exp.Current(), props);
            if (props.Mass() < 0) { inverted = true; break; }
        }
        if (inverted) {
            cbld.Add(resultComp, solid);
            continue;
        }

        TopoDS_Shape downScaled = BRepBuilderAPI_Transform(current, down).Shape();
        cbld.Add(resultComp, downScaled);
    }

    return resultComp;
}

} // namespace mvb
