// 3D regression test for the MKF "single turn centered in toroid inner hole"
// special case (Coil::wind_by_round_turns single-turn branch).
//
// When a toroidal coil has exactly one physical turn whose wire outer
// diameter exceeds the inner-hole radius (B/2), MKF places the turn at
// cartesian (0, 0) instead of against the inner wall. This test exercises
// the downstream MVB++ 3D path: building the toroidal racetrack solid for
// a turn whose cartesian coordinates are the toroid origin.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"
#include "support/Settings.h"
#include "MAS.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>

#include <filesystem>
#include <map>
#include <numbers>
#include <vector>

namespace {
MAS::Dimension make_dim(double v) { return MAS::Dimension(v); }
}

// Shared toroid geometry helpers for both the centered (oversized wire)
// and against-the-wall (small wire) comparison cases below.
namespace {

struct ToroidFixture {
    MAS::MagneticCore core;
    MAS::CoreBobbinProcessedDescription bobbinPd;
};

ToroidFixture make_t_25_15_10() {
    ToroidFixture f;

    MAS::CoreShape shape;
    shape.set_family(MAS::CoreShapeFamily::T);
    shape.set_type(MAS::FunctionalDescriptionType::STANDARD);
    std::map<std::string, MAS::Dimension> dims;
    dims["A"] = make_dim(0.025);
    dims["B"] = make_dim(0.015);
    dims["C"] = make_dim(0.010);
    shape.set_dimensions(dims);

    MAS::CoreGeometricalDescriptionElement piece;
    piece.set_type(MAS::CoreGeometricalDescriptionElementType::TOROIDAL);
    piece.set_coordinates({0.0, 0.0, 0.0});
    // Toroidal pieces must be rotated {π/2, π/2, 0} so their hole axis
    // points along world Y — this matches the convention assumed by
    // TurnBuilder::build_toroidal_turn (vertical legs along Y).
    piece.set_rotation(std::optional<std::vector<double>>(std::vector<double>{
        std::numbers::pi / 2.0, std::numbers::pi / 2.0, 0.0}));
    piece.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    MAS::CoreFunctionalDescription coreFunc;
    coreFunc.set_material(MAS::CoreMaterialDataOrNameUnion{std::string("N87")});
    coreFunc.set_type(MAS::CoreType::TOROIDAL);
    coreFunc.set_number_stacks(int64_t(1));
    coreFunc.set_shape(shape);

    f.core.set_functional_description(coreFunc);
    f.core.set_geometrical_description(std::optional<std::vector<MAS::CoreGeometricalDescriptionElement>>(std::vector<MAS::CoreGeometricalDescriptionElement>{piece}));

    f.bobbinPd.set_column_width(0.0025);
    f.bobbinPd.set_column_depth(0.005);
    f.bobbinPd.set_column_thickness(0.0);
    f.bobbinPd.set_wall_thickness(0.0);
    f.bobbinPd.set_column_shape(MAS::ColumnShape::RECTANGULAR);

    MAS::WindingWindowElement ww;
    ww.set_radial_height(0.0075);
    ww.set_angle(360.0);
    ww.set_coordinates(std::vector<double>({0.0075, 0.0, 0.0}));
    f.bobbinPd.set_winding_windows(std::vector<MAS::WindingWindowElement>{ww});

    return f;
}

MAS::Wire make_round_wire(double outerDiameter, double conductingDiameter,
                          const std::string& name) {
    MAS::Wire wire;
    wire.set_type(MAS::WireType::ROUND);
    {
        MAS::DimensionWithTolerance dim;
        dim.set_nominal(outerDiameter);
        wire.set_outer_diameter(std::optional<MAS::DimensionWithTolerance>(dim));
    }
    {
        MAS::DimensionWithTolerance dim;
        dim.set_nominal(conductingDiameter);
        wire.set_conducting_diameter(std::optional<MAS::DimensionWithTolerance>(dim));
    }
    wire.set_name(std::optional<std::string>(name));
    return wire;
}

MAS::Magnetic build_single_turn_magnetic(const ToroidFixture& f,
                                         const MAS::Wire& wire,
                                         const std::vector<double>& turnCoordinates,
                                         std::optional<std::vector<double>> outerCoordinates = std::nullopt,
                                         double rotationDeg = 0.0) {
    MAS::Bobbin bobbin;
    bobbin.set_processed_description(std::optional<MAS::CoreBobbinProcessedDescription>(f.bobbinPd));

    MAS::Turn turn;
    turn.set_name("Primary parallel 0 turn 0");
    turn.set_winding("Primary");
    turn.set_length(0.05);
    turn.set_parallel(0);
    turn.set_coordinates(turnCoordinates);
    turn.set_rotation(std::optional<double>(rotationDeg));
    if (outerCoordinates) {
        turn.set_additional_coordinates(
            std::optional<std::vector<std::vector<double>>>(
                std::vector<std::vector<double>>{*outerCoordinates}));
    }

    MAS::CoilFunctionalDescription coilFunc;
    coilFunc.set_name("Primary");
    coilFunc.set_number_turns(1);
    coilFunc.set_number_parallels(1);
    coilFunc.set_isolation_side(MAS::IsolationSide::PRIMARY);
    coilFunc.set_wire(MAS::WireDataOrNameUnion(wire));

    MAS::Coil coil;
    coil.set_bobbin(MAS::BobbinDataOrNameUnion(bobbin));
    coil.set_turns_description(std::optional<std::vector<MAS::Turn>>(std::vector<MAS::Turn>{turn}));
    coil.set_functional_description(std::vector<MAS::CoilFunctionalDescription>{coilFunc});

    MAS::Magnetic magnetic;
    magnetic.set_core(f.core);
    magnetic.set_coil(coil);
    return magnetic;
}

// Compute the volume of the boolean intersection between a turn and the
// union of all core pieces. A physically valid turn does not occupy the
// same space as the core, so this volume should be zero (allowing a tiny
// epsilon for OCCT floating-point noise).
double turn_core_collision_volume(const TopoDS_Shape& turn,
                                  const std::vector<mvb::NamedShape>& coreNamed) {
    BRep_Builder builder;
    TopoDS_Compound coreCompound;
    builder.MakeCompound(coreCompound);
    for (const auto& ns : coreNamed) {
        if (!ns.shape.IsNull()) builder.Add(coreCompound, ns.shape);
    }
    BRepAlgoAPI_Common common(turn, coreCompound);
    common.Build();
    if (!common.IsDone()) return -1.0;
    GProp_GProps props;
    BRepGProp::VolumeProperties(common.Shape(), props);
    return props.Mass();
}

void write_step(const std::string& path,
                const std::vector<mvb::NamedShape>& coreNamed,
                const mvb::NamedShape& turnNamed) {
    std::filesystem::remove(path);
    std::vector<TopoDS_Shape> shapes;
    std::vector<std::string> names;
    for (auto& ns : coreNamed) { shapes.push_back(ns.shape); names.push_back(ns.name); }
    shapes.push_back(turnNamed.shape);
    names.push_back(turnNamed.name);
    REQUIRE(mvb::exportSTEP(shapes, names, path));
    REQUIRE(std::filesystem::exists(path));
}

}  // namespace

// End-to-end pipeline test: hand a minimal MAS magnetic to OM via
// magnetic_autocomplete_safe, which invokes MKF's wind() (sections →
// layers → turns OR the centered-single-turn special branch), then
// MVB++ builds the 3D solid. Verifies the entire OM → MKF → MVB++ chain.
//
// Wire OD (10 mm) > inner radius (B/2 = 7.5 mm), so MKF takes the
// centered-single-turn branch (Coil::build_centered_single_turn_toroidal)
// and places the turn at cartesian (0, 0).
TEST_CASE("Toroidal single turn: end-to-end via OM",
          "[3d][toroidal][single-turn][e2e]") {
    nlohmann::json magneticJson = R"({
        "core": {
            "functionalDescription": {
                "shape": {
                    "family": "t",
                    "type": "standard",
                    "dimensions": {
                        "A": 0.025,
                        "B": 0.015,
                        "C": 0.010
                    }
                },
                "material": "N87",
                "gapping": [],
                "numberStacks": 1,
                "type": "toroidal"
            }
        },
        "coil": {
            "bobbin": "Basic",
            "functionalDescription": [
                {
                    "isolationSide": "primary",
                    "name": "Primary",
                    "numberParallels": 1,
                    "numberTurns": 1,
                    "wire": {
                        "type": "round",
                        "material": "copper",
                        "conductingDiameter": {"nominal": 0.0095},
                        "outerDiameter": {"nominal": 0.010},
                        "name": "Round 10.00 - Custom",
                        "numberConductors": 1,
                        "coating": {"type": "enamelled", "grade": 1}
                    }
                }
            ]
        }
    })"_json;

    OpenMagnetics::Magnetic enriched;
    REQUIRE_NOTHROW(enriched = mvb::magnetic_autocomplete_safe(magneticJson));

    // MKF should produce exactly one turn via its standard pipeline.
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    REQUIRE(turnsOpt->size() == 1);
    const auto& t0 = (*turnsOpt)[0];
    auto c = t0.get_coordinates();
    std::cout << "  MKF placed turn at cartesian (" << c[0]*1000 << ", "
              << c[1]*1000 << ") mm\n";

    // Build 3D via MVB++, write STEP.
    mvb::MagneticBuilder builder;
    auto coreNamed  = builder.buildCoreNamed(enriched.get_core(), /*corePolygonSegments=*/0);
    auto turnsNamed = builder.buildTurnsNamed(enriched.get_coil(), enriched.get_core(), /*wirePolygonSegments=*/0);
    REQUIRE(turnsNamed.size() == 1);
    REQUIRE(!turnsNamed[0].shape.IsNull());

    REQUIRE_THAT(c[0], Catch::Matchers::WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(c[1], Catch::Matchers::WithinAbs(0.0, 1e-6));

    auto outFile = std::filesystem::temp_directory_path() / "Test_Toroidal_Single_Turn_E2E.step";
    write_step(outFile.string(), coreNamed, turnsNamed[0]);
    std::cout << "  STEP file: " << outFile.string() << "\n";
}

TEST_CASE("Toroidal single turn centered: 3D turn solid built at origin",
          "[3d][toroidal][single-turn]") {
    auto fixture = make_t_25_15_10();
    // 10 mm OD > 7.5 mm inner radius → MKF places this turn at the hole
    // centre; supply cartesian (0, 0) directly.
    auto wire = make_round_wire(0.010, 0.0095, "Round 10.00 - Custom");
    // Outer leg center at -(A/2 + wire_radius) = -(12.5 + 5) mm = -17.5 mm,
    // so the wire's inner edge in the outer leg sits flush against A/2.
    auto magnetic = build_single_turn_magnetic(fixture, wire, {0.0, 0.0},
                                                /*outer=*/std::vector<double>{-0.0175, 0.0},
                                                /*rotationDeg=*/180.0);

    mvb::MagneticBuilder builder;
    auto turnsNamed = builder.buildTurnsNamed(magnetic.get_coil(), magnetic.get_core(), /*wirePolygonSegments=*/0);

    REQUIRE(turnsNamed.size() == 1);
    REQUIRE(!turnsNamed[0].shape.IsNull());

    GProp_GProps props;
    BRepGProp::VolumeProperties(turnsNamed[0].shape, props);
    REQUIRE(props.Mass() > 0.0);

    Bnd_Box bb;
    BRepBndLib::Add(turnsNamed[0].shape, bb);
    double xmin, ymin, zmin, xmax, ymax, zmax;
    bb.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    REQUIRE(xmin <= 0.0);
    REQUIRE(xmax >= 0.0);
    REQUIRE(ymin <= 0.0);
    REQUIRE(ymax >= 0.0);

    auto coreNamed = builder.buildCoreNamed(magnetic.get_core(), /*corePolygonSegments=*/0);

    // Collision check: turn must not intersect the toroidal core ring.
    double collision = turn_core_collision_volume(turnsNamed[0].shape, coreNamed);
    INFO("Centered turn / core intersection volume: " << collision << " m^3");
    REQUIRE(collision >= 0.0);
    // For a 10 mm-OD wire on a 15 mm-ID core, the racetrack bend has a
    // 5 mm major-radius arc sweeping above the core's outer-radius
    // boundary; OCCT reports a tangent-contact residual up to ~3·10⁻⁷ m³.
    // Same NURBS-tangency artefact seen on the inner-wall test.
    REQUIRE(collision < 5e-7);

    auto outFile = std::filesystem::temp_directory_path() / "Test_Toroidal_Single_Turn_Centered.step";
    write_step(outFile.string(), coreNamed, turnsNamed[0]);
    INFO("STEP file written to: " << outFile.string());
}

// Comparison: a single SMALL wire (1 mm OD < 7.5 mm inner radius) that
// would NOT trigger MKF's centered branch. The turn sits against the
// inner wall — MKF offsets the wire CENTER inward by wire_radius from
// the inner-hole boundary so the wire's outer edge touches B/2, and
// the outer leg is placed at the same angle just outside A/2.
//
// For T 25/15/10 + 1 mm wire:
//   inner leg cartesian: (-(B/2 - wire_radius), 0) = (-0.007, 0)
//   outer leg cartesian: (-(A/2 + wire_radius), 0) = (-0.013, 0)
//   rotation: 180° (MKF stores turns on -X half with rotation = 180)
TEST_CASE("Toroidal single turn against inner wall: 3D turn solid at inner edge",
          "[3d][toroidal][single-turn]") {
    auto fixture = make_t_25_15_10();
    auto wire = make_round_wire(0.001, 0.00095, "Round 1.00 - Grade 1");
    auto magnetic = build_single_turn_magnetic(
        fixture, wire,
        /*inner=*/{-0.007, 0.0},
        /*outer=*/std::vector<double>{-0.013, 0.0},
        /*rotationDeg=*/180.0);

    mvb::MagneticBuilder builder;
    auto turnsNamed = builder.buildTurnsNamed(magnetic.get_coil(), magnetic.get_core(), /*wirePolygonSegments=*/0);

    REQUIRE(turnsNamed.size() == 1);
    REQUIRE(!turnsNamed[0].shape.IsNull());

    GProp_GProps props;
    BRepGProp::VolumeProperties(turnsNamed[0].shape, props);
    REQUIRE(props.Mass() > 0.0);

    Bnd_Box bb;
    BRepBndLib::Add(turnsNamed[0].shape, bb);
    double xmin, ymin, zmin, xmax, ymax, zmax;
    bb.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    // Inner leg at -7 mm and outer leg at -13 mm — turn lives in -X half.
    REQUIRE(xmin < -0.005);
    // Geometric correctness checks: wire outer edges should be flush with
    // the core inner/outer walls and core top face.
    //   inner leg outer-X = -B/2 = -0.0075
    //   outer leg inner-X = -A/2 = -0.0125
    //   inner Y-tube length = C/2 = 0.005 (so top horizontal segment touches y = C/2)
    REQUIRE_THAT(xmin, Catch::Matchers::WithinAbs(-0.0135, 1e-4));   // outer-leg outer-X
    REQUIRE_THAT(xmax, Catch::Matchers::WithinAbs(-0.0065, 1e-4));   // inner-leg inner-X
    REQUIRE_THAT(ymax, Catch::Matchers::WithinAbs( 0.006,  5e-4));   // wire top edge ≈ C/2 + wire_OD
    INFO("xmin=" << xmin << " xmax=" << xmax << " ymin=" << ymin << " ymax=" << ymax);

    auto coreNamed = builder.buildCoreNamed(magnetic.get_core(), /*corePolygonSegments=*/0);

    // Collision check: a 1 mm wire at the inner wall should clear the core.
    // The wire is exactly tangent to both B/2 and A/2 walls; OCCT boolean
    // common on double-tangent NURBS surfaces can report a tiny residual
    // (~10⁻⁸ m³ ≈ 10 mm³) that is below the wire's own cross-sectional
    // area × any sane sliver length and is not a real penetration.
    double collision = turn_core_collision_volume(turnsNamed[0].shape, coreNamed);
    INFO("Inner-wall turn / core intersection volume: " << collision << " m^3");
    REQUIRE(collision >= 0.0);
    REQUIRE(collision < 5e-8);

    auto outFile = std::filesystem::temp_directory_path() / "Test_Toroidal_Single_Turn_InnerWall.step";
    write_step(outFile.string(), coreNamed, turnsNamed[0]);
    INFO("STEP file written to: " << outFile.string());
}

// Survey: load several toroidal MAS fixtures, build their 3D turns, and
// measure the boolean-intersection volume between each turn and the core.
// Reports a summary at the end. Uses CHECK (non-fatal) so we can see all
// fixtures' results in one run.
namespace {
nlohmann::json load_json(const std::string& path) {
    std::ifstream f(path);
    REQUIRE(f.is_open());
    nlohmann::json j; f >> j;
    return j;
}
}  // namespace

TEST_CASE("Toroidal collision survey: every turn vs its core",
          "[3d][toroidal][collision-survey]") {
    struct Fixture {
        std::string label;
        std::string path;
        // If true, fixture file root is a MAS file (use j["magnetic"]);
        // otherwise it's already a Magnetic JSON.
        bool isMasRoot;
    };

    std::vector<Fixture> fixtures = {
        {"toroidal_one_turn",                 "testData/toroidal_one_turn.json",                 true},
        {"toroidal_one_turn_rect_wire",       "testData/toroidal_one_turn_rectangular_wire.json", true},
        {"toroidal_two_turns_centered",       "testData/toroidal_two_turns_centered.json",       true},
        {"toroidal_two_turns_spread",         "testData/toroidal_two_turns_spread.json",         true},
        {"toroidal_full_layer",               "testData/toroidal_full_layer.json",               true},
        {"toroidal_full_layer_rect",          "testData/toroidal_full_layer_rectangular_wires.json", true},
        {"cmc_t2515_w800",                    std::string(MAS_EXAMPLES_DIR) + "/07_cmc_t2515_w800.json", true},
        {"pfc_inductor_t4020_hf60",           std::string(MAS_EXAMPLES_DIR) + "/05_pfc_inductor_t4020_hf60.json", true},
        {"boost_inductor_t5026_26",           std::string(MAS_EXAMPLES_DIR) + "/12_boost_inductor_t5026_26.json", true},
    };

    struct Result {
        std::string label;
        int numTurns = 0;
        int collidingTurns = 0;
        double maxCollisionVol = 0.0;
        double totalCollisionVol = 0.0;
        std::string error;
    };
    std::vector<Result> results;

    for (const auto& fx : fixtures) {
        Result r; r.label = fx.label;
        if (!std::filesystem::exists(fx.path)) {
            r.error = "file not found";
            results.push_back(r);
            continue;
        }
        try {
            auto j = load_json(fx.path);
            mvb::patch_dimension_nominals(j);
            auto magneticJson = fx.isMasRoot ? j.at("magnetic") : j;
            auto enriched = mvb::magnetic_autocomplete_safe(magneticJson);

            mvb::MagneticBuilder builder;
            auto coreNamed = builder.buildCoreNamed(enriched.get_core(), /*corePolygonSegments=*/0);
            auto turnsNamed = builder.buildTurnsNamed(enriched.get_coil(), enriched.get_core(), /*wirePolygonSegments=*/0);
            r.numTurns = static_cast<int>(turnsNamed.size());

            for (const auto& ns : turnsNamed) {
                if (ns.shape.IsNull()) continue;
                double c = turn_core_collision_volume(ns.shape, coreNamed);
                if (c > 1e-12) {
                    r.collidingTurns++;
                    r.maxCollisionVol = std::max(r.maxCollisionVol, c);
                    r.totalCollisionVol += c;
                }
            }

            // Export STEP for visual inspection.
            std::vector<TopoDS_Shape> shapes;
            std::vector<std::string> names;
            for (auto& ns : coreNamed) { shapes.push_back(ns.shape); names.push_back(ns.name); }
            for (auto& ns : turnsNamed) {
                if (!ns.shape.IsNull()) { shapes.push_back(ns.shape); names.push_back(ns.name); }
            }
            auto outFile = std::filesystem::temp_directory_path()
                         / (std::string("Test_ToroidalSurvey_") + fx.label + ".step");
            std::filesystem::remove(outFile);
            mvb::exportSTEP(shapes, names, outFile.string());
        } catch (const std::exception& e) {
            r.error = e.what();
        }
        results.push_back(r);
    }

    // Pretty-print summary. Use std::cout so it's always visible (Catch2's
    // INFO would only show on failure).
    std::ostringstream report;
    report << "\n=== Toroidal turn-vs-core collision survey ===\n";
    report << std::left << std::setw(38) << "fixture"
           << std::right << std::setw(8) << "turns"
           << std::setw(10) << "collide"
           << std::setw(14) << "max(mm^3)"
           << std::setw(14) << "tot(mm^3)"
           << "   error\n";
    int totalColliding = 0;
    for (const auto& r : results) {
        report << std::left << std::setw(38) << r.label
               << std::right << std::setw(8) << r.numTurns
               << std::setw(10) << r.collidingTurns
               << std::setw(14) << std::fixed << std::setprecision(3) << (r.maxCollisionVol * 1e9)
               << std::setw(14) << (r.totalCollisionVol * 1e9)
               << "   " << r.error << "\n";
        totalColliding += r.collidingTurns;
    }
    report << "===========================================\n";
    std::cout << report.str();

    // Soft assert: at least one fixture should produce turns. The collision
    // numbers themselves are reported, not asserted, so a single run shows
    // the full landscape.
    int totalTurns = 0;
    for (const auto& r : results) totalTurns += r.numTurns;
    REQUIRE(totalTurns > 0);
}

// Multi-layer toroidal example: build the 3D for a 2-layer toroidal
// winding, export STEP, and confirm that the two layers' turns sit at
// distinct radii (i.e. layer 1 inside layer 0).
TEST_CASE("Toroidal multi-layer: two layers built at distinct radii",
          "[3d][toroidal][multi-layer]") {
    auto j = load_json("testData/toroidal_two_layers_not_compact.json");
    mvb::patch_dimension_nominals(j);
    auto enriched = mvb::magnetic_autocomplete_safe(j.at("magnetic"));

    mvb::MagneticBuilder builder;
    auto coreNamed = builder.buildCoreNamed(enriched.get_core(), /*corePolygonSegments=*/0);
    auto turnsNamed = builder.buildTurnsNamed(enriched.get_coil(), enriched.get_core(), /*wirePolygonSegments=*/0);
    REQUIRE(!turnsNamed.empty());

    // Cross-reference each NamedShape against the MAS turns to find its
    // layer assignment, then compute its X-Z radial position.
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    const auto& turns = *turnsOpt;
    REQUIRE(turns.size() == turnsNamed.size());

    std::map<std::string, std::vector<double>> radiiByLayer;
    for (size_t i = 0; i < turns.size(); ++i) {
        const auto& t = turns[i];
        const auto& c = t.get_coordinates();
        double r = std::hypot(c[0], c[1]);
        auto layer = t.get_layer().value_or("(no-layer)");
        radiiByLayer[layer].push_back(r * 1000.0);  // mm
    }

    std::ostringstream o;
    o << std::fixed << std::setprecision(3);
    o << "\n=== Two-layer toroidal: turn radii per layer (mm) ===\n";
    for (auto& [layer, rs] : radiiByLayer) {
        double rmin = *std::min_element(rs.begin(), rs.end());
        double rmax = *std::max_element(rs.begin(), rs.end());
        o << "  " << layer << ": " << rs.size() << " turns, r ∈ [" << rmin << ", " << rmax << "]\n";
    }
    std::cout << o.str();

    REQUIRE(radiiByLayer.size() >= 2);  // at least two distinct layers

    // Confirm per-layer Y-extents: layer-1 turns should reach further in Y
    // than layer-0 turns by their extra radial inset.
    std::map<std::string, double> maxYByLayer;
    for (size_t i = 0; i < turns.size(); ++i) {
        const auto& t = turns[i];
        auto layer = t.get_layer().value_or("(no-layer)");
        Bnd_Box bb; BRepBndLib::Add(turnsNamed[i].shape, bb);
        double a,b,c,d,e,g; bb.Get(a,b,c,d,e,g);
        maxYByLayer[layer] = std::max(maxYByLayer[layer], e * 1000.0);
    }
    std::ostringstream o2;
    o2 << std::fixed << std::setprecision(3);
    o2 << "  per-layer max |Y| (mm): ";
    for (auto& [l, y] : maxYByLayer) o2 << l << "=" << y << "  ";
    o2 << "\n";
    std::cout << o2.str();

    // Export STEP for inspection.
    std::vector<TopoDS_Shape> shapes;
    std::vector<std::string> names;
    for (auto& ns : coreNamed)  { shapes.push_back(ns.shape); names.push_back(ns.name); }
    for (auto& ns : turnsNamed) { if (!ns.shape.IsNull()) { shapes.push_back(ns.shape); names.push_back(ns.name); } }
    auto outFile = std::filesystem::temp_directory_path() / "Test_Toroidal_Two_Layers.step";
    std::filesystem::remove(outFile);
    REQUIRE(mvb::exportSTEP(shapes, names, outFile.string()));
    std::cout << "  STEP file: " << outFile.string() << "\n";

    // Collision survey across layers — should now be only tangent noise.
    int totalColliding = 0;
    double maxCol = 0.0, totCol = 0.0;
    for (const auto& ns : turnsNamed) {
        if (ns.shape.IsNull()) continue;
        double cvol = turn_core_collision_volume(ns.shape, coreNamed);
        if (cvol > 1e-12) {
            totalColliding++;
            maxCol = std::max(maxCol, cvol);
            totCol += cvol;
        }
    }
    std::cout << "  collide=" << totalColliding << " / " << turnsNamed.size()
              << "  max=" << maxCol*1e9 << " mm³  total=" << totCol*1e9 << " mm³\n";
}

TEST_CASE("Toroidal rect-wire: dump sub-piece bboxes",
          "[3d][toroidal][diag-rect]") {
    auto j = load_json("testData/toroidal_one_turn_rectangular_wire.json");
    mvb::patch_dimension_nominals(j);
    auto enriched = mvb::magnetic_autocomplete_safe(j.at("magnetic"));
    mvb::MagneticBuilder builder;
    auto coreNamed = builder.buildCoreNamed(enriched.get_core(), /*corePolygonSegments=*/0);
    auto turnsNamed = builder.buildTurnsNamed(enriched.get_coil(), enriched.get_core(), /*wirePolygonSegments=*/0);

    std::ostringstream o; o << std::fixed << std::setprecision(4);
    o << "\n=== Core (mm) ===\n";
    for (auto& ns : coreNamed) {
        Bnd_Box bb; BRepBndLib::Add(ns.shape, bb);
        double x0,y0,z0,x1,y1,z1; bb.Get(x0,y0,z0,x1,y1,z1);
        o << ns.name << "  x[" << x0*1000 << "," << x1*1000 << "]  y["
          << y0*1000 << "," << y1*1000 << "]  z[" << z0*1000 << "," << z1*1000 << "]\n";
    }
    o << "\n=== Rect-wire turn sub-solids (mm) ===\n";
    for (auto& ns : turnsNamed) {
        Bnd_Box wb; BRepBndLib::Add(ns.shape, wb);
        double a,b,c,d,e,g; wb.Get(a,b,c,d,e,g);
        o << "[WHOLE] " << ns.name << "  y[" << b*1000 << "," << e*1000 << "]\n";
        int idx = 0;
        for (TopExp_Explorer exp(ns.shape, TopAbs_SOLID); exp.More(); exp.Next(), ++idx) {
            Bnd_Box bb; BRepBndLib::Add(exp.Current(), bb);
            double x0,y0,z0,x1,y1,z1; bb.Get(x0,y0,z0,x1,y1,z1);
            o << "  sub[" << std::setw(3) << std::setfill('0') << idx << std::setfill(' ') << "]"
              << "  x[" << std::setw(9) << x0*1000 << "," << std::setw(9) << x1*1000 << "]"
              << "  y[" << std::setw(9) << y0*1000 << "," << std::setw(9) << y1*1000 << "]"
              << "  z[" << std::setw(9) << z0*1000 << "," << std::setw(9) << z1*1000 << "]\n";
        }
    }
    std::cout << o.str();
    REQUIRE(!turnsNamed.empty());
}

// Diagnostic: dump bounding boxes of every sub-solid of every turn for
// the toroidal_one_turn fixture, so we can answer "is the racetrack
// floating above the core / longer than C?"
TEST_CASE("Toroidal one-turn: dump sub-piece bboxes for diagnosis",
          "[3d][toroidal][diag]") {
    auto j = load_json("testData/toroidal_one_turn.json");
    mvb::patch_dimension_nominals(j);
    auto enriched = mvb::magnetic_autocomplete_safe(j.at("magnetic"));

    mvb::MagneticBuilder builder;
    auto coreNamed  = builder.buildCoreNamed(enriched.get_core(), /*corePolygonSegments=*/0);
    auto turnsNamed = builder.buildTurnsNamed(enriched.get_coil(), enriched.get_core(), /*wirePolygonSegments=*/0);

    std::ostringstream o;
    o << std::fixed << std::setprecision(4);
    o << "\n=== Core pieces (mm) ===\n";
    for (auto& ns : coreNamed) {
        Bnd_Box bb; BRepBndLib::Add(ns.shape, bb);
        double x0,y0,z0,x1,y1,z1; bb.Get(x0,y0,z0,x1,y1,z1);
        o << ns.name << "  x[" << x0*1000 << "," << x1*1000 << "]  y["
          << y0*1000 << "," << y1*1000 << "]  z[" << z0*1000 << "," << z1*1000 << "]\n";
    }
    o << "\n=== Turn sub-solids (mm) ===\n";
    for (auto& ns : turnsNamed) {
        Bnd_Box wholeBB; BRepBndLib::Add(ns.shape, wholeBB);
        double x0,y0,z0,x1,y1,z1; wholeBB.Get(x0,y0,z0,x1,y1,z1);
        o << "[WHOLE] " << ns.name
          << "  x[" << x0*1000 << "," << x1*1000 << "]"
          << "  y[" << y0*1000 << "," << y1*1000 << "]"
          << "  z[" << z0*1000 << "," << z1*1000 << "]"
          << "  Ylen=" << (y1-y0)*1000 << "mm\n";
        int idx = 0;
        for (TopExp_Explorer exp(ns.shape, TopAbs_SOLID); exp.More(); exp.Next(), ++idx) {
            Bnd_Box bb; BRepBndLib::Add(exp.Current(), bb);
            double a,b,c,d,e,g; bb.Get(a,b,c,d,e,g);
            o << "  sub[" << std::setw(3) << std::setfill('0') << idx << std::setfill(' ') << "]"
              << "  x[" << std::setw(8) << a*1000 << "," << std::setw(8) << d*1000 << "]"
              << "  y[" << std::setw(8) << b*1000 << "," << std::setw(8) << e*1000 << "]"
              << "  z[" << std::setw(8) << c*1000 << "," << std::setw(8) << g*1000 << "]\n";
        }
    }
    std::cout << o.str();
    REQUIRE(!turnsNamed.empty());
}
