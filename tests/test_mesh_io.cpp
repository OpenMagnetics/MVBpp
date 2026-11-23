// Round-trip gate for the mesh writers (ABT #1588): every format is read back and compared with
// the mesh it was written from. The mesh here is a hand-made cell, not a magnetic: it pins the
// format rules (ids, order, region names, units, refusals); the real meshes are round-tripped by
// tools/mvbpp_mesh_roundtrip on OMFEM's reference outputs.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "mvb/mesh/Mesh2d.h"
#include "mvb/mesh/MeshSupport.h"
#include "mvb/mesh/Mesher.h"
#include "mvb/mesh/SizeField.h"
#include "mvb/mesh/SurfaceCorridors.h"

#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

using namespace mvb::mesh;
using Catch::Matchers::ContainsSubstring;

namespace {

// Two tets sharing a face (one per volume region), a pyramid on the core, and the free faces of
// the air tet as a boundary region. Coordinates are not round in any decimal unit.
Mesh cell() {
    Mesh m;
    m.nodeIds = {1, 2, 3, 4, 5, 7, 8, 9};
    m.xyz = {0.0, 0.0, 0.0,
             0.0123456789012345678, 0.0, 0.0,
             0.0, 0.0098765432109876543, 0.0,
             0.0, 0.0, 0.0071234567890123457,
             0.0123456789012345678, 0.0098765432109876543, 0.0071234567890123457,
             -1.2345678901234567e-7, -0.0098765432109876543, 0.0,
             0.0123456789012345678, -0.0098765432109876543, 0.0,
             0.006, -0.005, -0.004};
    Region core{1, 3, "core", {MaterialKind::Core, nlohmann::json{
        {"name", "Test ferrite"}, {"density", 4800.0},
        {"heatConductivity", {{"nominal", 4.2}}}, {"heatCapacity", {{"minimum", 700.0}, {"maximum", 800.0}}},
        {"permeability", {{"initial", {{"value", 2200}}}}}, {"resistivity", {{{"value", 10.0}}}}}}};
    Region air{2, 3, "air", {MaterialKind::None, {}}};
    Region outer{4, 2, "outer", {MaterialKind::None, {}}};
    m.regions = {core, air, outer};
    m.blocks = {
        {ElementType::Tet4, 1, 11, {10}, {1, 2, 3, 4}},
        {ElementType::Pyramid5, 1, 11, {12}, {1, 2, 8, 7, 9}},
        {ElementType::Tet4, 2, 12, {20}, {2, 3, 4, 5}},
        {ElementType::Tri3, 4, 31, {40, 41, 42}, {2, 5, 3, 3, 5, 4, 4, 5, 2}},
    };
    return m;
}

// A BDF has no elementary entity: drop it and join the blocks that differed only by it, as the
// reader does, so the rest of the mesh compares exactly.
Mesh withoutEntities(const Mesh& m) {
    Mesh out = m;
    out.blocks.clear();
    for (auto b : m.blocks) {
        b.entity = 0;
        if (!out.blocks.empty() && out.blocks.back().type == b.type && out.blocks.back().region == b.region) {
            auto& last = out.blocks.back();
            last.ids.insert(last.ids.end(), b.ids.begin(), b.ids.end());
            last.nodes.insert(last.nodes.end(), b.nodes.begin(), b.nodes.end());
        } else {
            out.blocks.push_back(b);
        }
    }
    return out;
}

double tetVolume(const Mesh& m, const std::int64_t* n) {
    auto at = [&](std::int64_t id) {
        for (std::size_t i = 0; i < m.nodeIds.size(); ++i)
            if (m.nodeIds[i] == id) return &m.xyz[3 * i];
        FAIL("node " << id << " missing");
        return static_cast<const double*>(nullptr);
    };
    const double *a = at(n[0]), *b = at(n[1]), *c = at(n[2]), *d = at(n[3]);
    double u[3], v[3], w[3];
    for (int k = 0; k < 3; ++k) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; w[k] = d[k] - a[k]; }
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) +
            u[2] * (v[0] * w[1] - v[1] * w[0])) / 6.0;
}

// What a BDF cannot carry is the gmsh elementary entity; compare everything else, coordinates to
// the 16-character field's precision: 11 to 13 significant digits ("-1.234567890123-2" fills
// it), so 1e-11 relative, about 1e-13 m on a centimetre.
void requireBdfRoundTrip(const Mesh& original, const Mesh& back) {
    REQUIRE(back.nodeIds == original.nodeIds);
    for (std::size_t i = 0; i < original.xyz.size(); ++i) {
        const double x = original.xyz[i], y = back.xyz[i];
        INFO("coordinate " << i << ": " << x << " vs " << y);
        REQUIRE(std::abs(x - y) <= 1e-11 * std::max(std::abs(x), 1e-9));
    }
    Mesh a = withoutEntities(original), b = back;
    b.xyz = a.xyz;
    REQUIRE(diffMeshes(a, b) == "");
    for (std::size_t i = 0; i < a.blocks.size(); ++i)
        if (a.blocks[i].type == ElementType::Tet4)
            for (std::size_t e = 0; e < a.blocks[i].ids.size(); ++e) {
                const double v0 = tetVolume(original, &a.blocks[i].nodes[4 * e]);
                const double v1 = tetVolume(back, &back.blocks[i].nodes[4 * e]);
                REQUIRE((v0 > 0) == (v1 > 0));
                REQUIRE(std::abs(v1) >= 1e-18);
            }
}

}  // namespace

TEST_CASE("msh2 round-trips a mesh exactly", "[mesh-io]") {
    const Mesh m = cell();
    const std::string text = exportMesh(m, "msh2", "m");
    const Mesh back = importMesh(text, "msh2");
    REQUIRE(diffMeshes(m, back) == "");
    // And the comparator is not blind: one coordinate moved by one ulp is a difference.
    Mesh moved = back;
    moved.xyz[4] = std::nextafter(moved.xyz[4], 1.0);
    REQUIRE_THAT(diffMeshes(m, moved), ContainsSubstring("node 2 coordinate 1"));
}

TEST_CASE("msh2 refuses a unit it cannot state", "[mesh-io]") {
    REQUIRE_THROWS_WITH(exportMesh(cell(), "msh2", "mm"), ContainsSubstring("metres only"));
}

TEST_CASE("Nastran BDF round-trips in metres and millimetres", "[mesh-io]") {
    const Mesh m = cell();
    for (const char* unit : {"m", "mm"}) {
        INFO("unit " << unit);
        const std::string text = exportMesh(m, "bdf", unit);
        REQUIRE_THAT(text, ContainsSubstring(std::string("$ LENGTH UNIT: ") + unit));
        requireBdfRoundTrip(m, importMesh(text, "bdf"));
    }
}

TEST_CASE("Nastran BDF carries the MAS material, in the unit system of the length", "[mesh-io]") {
    const std::string m = exportMesh(cell(), "bdf", "m");
    // K, CP (resolved from the min/max range by MKF's resolver: 750) and RHO, SI.
    REQUIRE_THAT(m, ContainsSubstring("MAT4*   1               4.2+0           7.5+2           4.8+3"));
    REQUIRE_THAT(m, ContainsSubstring("$ MATERIAL core Test ferrite"));
    REQUIRE_THAT(m, ContainsSubstring("$ EM {\"permeability\""));
    // mm-t-s-K: K unchanged, CP x 1e6, RHO x 1e-12.
    const std::string mm = exportMesh(cell(), "bdf", "mm");
    REQUIRE_THAT(mm, ContainsSubstring("MAT4*   1               4.2+0           7.5+8           4.8-9"));
    REQUIRE_THAT(mm, ContainsSubstring("$ UNIT SYSTEM: mm-t-s-K"));
    // Air keeps its property and has no MAT card.
    REQUIRE_THAT(m, ContainsSubstring("$ MATERIAL none -\nPSOLID* 2               2\n$ REGION"));
}

TEST_CASE("Nastran BDF refuses a material MAS cannot complete", "[mesh-io]") {
    Mesh m = cell();
    m.regions[0].material.record.erase("heatConductivity");
    REQUIRE_THROWS_WITH(exportMesh(m, "bdf", "m"),
                        ContainsSubstring("MAS lacks heatConductivity"));
    // A conductor's conductivity is a temperature table: no reference temperature, no MAT4.
    m = cell();
    m.regions[0].material = {MaterialKind::Wire, nlohmann::json{
        {"name", "copper"}, {"permeability", 0.999994}, {"density", 8890.0},
        {"resistivity", {{"referenceValue", 1.678e-08}, {"referenceTemperature", 20}, {"temperatureCoefficient", 0.004041}}},
        {"thermalConductivity", {{{"temperature", 0}, {"value", 401}}, {{"temperature", 127}, {"value", 392}}}}}};
    REQUIRE_THROWS_WITH(exportMesh(m, "bdf", "m"), ContainsSubstring("options.temperature"));
    // With one, it still needs the specific heat.
    REQUIRE_THROWS_WITH(exportMesh(m, "bdf", "m", {{"temperature", 63.5}}), ContainsSubstring("MAS lacks specificHeat"));
}

TEST_CASE("Nastran BDF writes a conductor at the reference temperature, with its MAS table", "[mesh-io]") {
    Mesh m = cell();
    m.regions[0].material = {MaterialKind::Wire, nlohmann::json{
        {"name", "copper"}, {"permeability", 0.999994}, {"density", 8890.0}, {"specificHeat", 385.2},
        {"resistivity", {{"referenceValue", 1.678e-08}, {"referenceTemperature", 20}, {"temperatureCoefficient", 0.004041}}},
        {"thermalConductivity", {{{"temperature", 127}, {"value", 392}}, {{"temperature", 0}, {"value", 401}}}}}};
    // The ambient stands in when no temperature is asked for; an asked one wins.
    const std::string a = exportMesh(m, "bdf", "m", {{"ambientTemperature", 63.5}});
    REQUIRE_THAT(a, ContainsSubstring("MAT4*   1               3.965+2         3.852+2         8.89+3"));
    REQUIRE_THAT(a, ContainsSubstring("MATT4*  1               1\n"));
    REQUIRE_THAT(a, ContainsSubstring("TABLEM1*1\n*       0.              4.01+2          1.27+2          3.92+2\n*       ENDT"));
    REQUIRE_THAT(a, ContainsSubstring("reference temperature 63.5 C"));
    const std::string t = exportMesh(m, "bdf", "m", {{"ambientTemperature", 63.5}, {"temperature", 0.0}});
    REQUIRE_THAT(t, ContainsSubstring("MAT4*   1               4.01+2"));
    requireBdfRoundTrip(m, importMesh(a, "bdf"));
    REQUIRE_THROWS_WITH(exportMesh(m, "bdf", "m", {{"temperatur", 20.0}}), ContainsSubstring("unknown option"));
}

TEST_CASE("The ambient temperature comes from the MAS operating points, and only when they agree", "[mesh-io]") {
    auto mas = [](std::vector<double> ts) {
        nlohmann::json ops = nlohmann::json::array();
        for (double t : ts) ops.push_back({{"conditions", {{"ambientTemperature", t}}}});
        return nlohmann::json{{"inputs", {{"operatingPoints", ops}}}};
    };
    REQUIRE(ambientTemperature(mas({40.0, 40.0})) == 40.0);
    REQUIRE_THROWS_WITH(ambientTemperature(mas({25.0, 40.0})), ContainsSubstring("different ambients"));
    REQUIRE_THROWS_WITH(ambientTemperature(nlohmann::json::object()), ContainsSubstring("no inputs.operatingPoints"));
}

TEST_CASE("Exporters refuse a volume region with no material stated", "[mesh-io]") {
    Mesh m = cell();
    m.regions[1].material = {};
    REQUIRE_THROWS_WITH(exportMesh(m, "bdf", "m"), ContainsSubstring("'air' has no material set"));
    REQUIRE_THROWS_WITH(exportMaterials(m), ContainsSubstring("'air' has no material set"));
}

TEST_CASE("The materials companion lists every region with its MAS data", "[mesh-io]") {
    const auto j = exportMaterials(cell());
    REQUIRE(j.size() == 3);
    REQUIRE(j[0]["material"] == "Test ferrite");
    REQUIRE(j[0]["electromagnetic"].contains("permeability"));
    REQUIRE(j[0]["thermal"]["density"] == 4800.0);
    REQUIRE(j[1]["kind"] == "none");
    REQUIRE_FALSE(j[1].contains("material"));
}

TEST_CASE("Writers refuse a structurally broken mesh", "[mesh-io]") {
    Mesh m = cell();
    m.blocks[0].nodes[0] = 99;
    REQUIRE_THROWS_WITH(exportMesh(m, "msh2", "m"), ContainsSubstring("node 99, which does not exist"));
    m = cell();
    m.blocks[2].ids[0] = 10;
    REQUIRE_THROWS_WITH(exportMesh(m, "msh2", "m"), ContainsSubstring("element id 10 used twice"));
    m = cell();
    m.regions[2].tag = 1;
    REQUIRE_THROWS_WITH(exportMesh(m, "msh2", "m"), ContainsSubstring("tag 1 used twice"));
}

// Every winding is meshed at its own automatic size (Alf, 2026-11-21). two_switch_forward's 0.2 mm foil used to put its
// 0.12 mm on a 0.79 mm litz primary; a thick few-turn winding beside a thin many-turn one is the common case.
TEST_CASE("Each winding's automatic conductor target comes from its own wire", "[mesh-api][sizing]") {
    // a 0.79 mm litz bundle and a 0.2 mm foil (two_switch_forward's two wires) and a 0.355 mm round wire
    const std::map<std::string, double> target{{"Primary", wire_scale_target(0.3943e-3)},
                                               {"Secondary", foil_scale_target(0.2e-3)},
                                               {"Reset Winding", wire_scale_target(0.1775e-3)}};
    CHECK(wire_scale_target(0.1e-3) == 1.2 * 0.1e-3);
    CHECK(foil_scale_target(0.2e-3) == 2.0 * 0.2e-3);
    CHECK(foil_scale_target(0.2e-3) > 3.3 * wire_scale_target(0.1e-3));   // a foil is not sized as a wire of half its thickness
    REQUIRE_THROWS_AS(wire_scale_target(0.0), std::invalid_argument);
    REQUIRE_THROWS_AS(foil_scale_target(0.0), std::invalid_argument);
    REQUIRE_THROWS_AS(foil_scale_target(-1e-4), std::invalid_argument);

    const auto t = automatic_conductor_targets({"winding_Primary", "winding_Secondary", "winding_Reset Winding"}, target);
    REQUIRE(t.size() == 3);
    CHECK(t.at("winding_Primary") == target.at("Primary"));
    CHECK(t.at("winding_Secondary") == target.at("Secondary"));
    CHECK(t.at("winding_Reset Winding") == target.at("Reset Winding"));
    CHECK(t.at("winding_Primary") > 2.2 * t.at("winding_Reset Winding"));   // the thick winding is not sized by the thin one
    // per-turn regions are sized by the winding they belong to; the lumped region belongs to none
    const auto turns = automatic_conductor_targets({"turn_Primary_0", "turn_Primary_17", "turn_Reset Winding_3", "winding"}, target);
    REQUIRE(turns.size() == 2);
    CHECK(turns.at("winding_Primary") == t.at("winding_Primary"));
    CHECK(turns.at("winding_Reset Winding") == t.at("winding_Reset Winding"));
    CHECK(winding_of_conductor_region("turn_Reset Winding_3") == std::optional<std::string>("Reset Winding"));
    CHECK_FALSE(winding_of_conductor_region("winding").has_value());
    CHECK_FALSE(winding_of_conductor_region("core").has_value());
    // no wire, no size: neither a default nor the design's thinnest
    REQUIRE_THROWS_WITH(automatic_conductor_targets({"winding_Auxiliary"}, target), ContainsSubstring("no wire to be sized from"));
    // a corridor binds the coarser of its two surfaces
    CHECK(corridor_binds(0.3e-3, t.at("winding_Primary"), t.at("winding_Reset Winding")));               // under the primary's 0.47 mm
    CHECK_FALSE(corridor_binds(0.3e-3, t.at("winding_Reset Winding"), t.at("winding_Reset Winding")));   // over the reset wire's 0.21 mm
    CHECK_FALSE(corridor_binds(0.5e-3, t.at("winding_Primary"), t.at("winding_Secondary")));
}

// The chord bound, as rules (ABT #2460). One bound for every conductor: a facet of size h on a surface of radius rho sags
// h^2 / 8 rho, held to 0.32 of the gap.
TEST_CASE("The chord bound, the widest gap it can bind and what counts as facing", "[mesh-api][sizing]") {
    // a 0.5 mm wire 1.733 um from its neighbour (03_buck's tightest), and the sag that size leaves: 0.32 of the gap
    const double h = chord_bound(0.25e-3, 1.733e-6);
    CHECK(h == Catch::Approx(0.8 * std::sqrt(4.0 * 0.25e-3 * 1.733e-6)).epsilon(1e-12));
    CHECK(h * h / (8.0 * 0.25e-3) == Catch::Approx(0.32 * 1.733e-6).epsilon(1e-12));
    CHECK(chord_bound(1e-3, 0.0) == 0.0);
    REQUIRE_THROWS_AS(chord_bound(0.0, 1e-6), std::invalid_argument);
    REQUIRE_THROWS_AS(chord_bound(1e-3, -1e-9), std::invalid_argument);
    // the widest binding gap is where the bound reaches the target ...
    const double wide = widest_binding_gap(10.2e-3, 0.4e-3);   // a foil wrap of 10.2 mm, meshed at 0.4 mm
    CHECK(wide == Catch::Approx(6.127e-6).epsilon(1e-3));
    CHECK(chord_bound(10.2e-3, wide) == Catch::Approx(0.4e-3).epsilon(1e-12));
    // ... and never more than the surface can sag at all: a 50 um fillet under a 0.4 mm target stops at rho / 0.32
    CHECK(widest_binding_gap(50e-6, 0.4e-3) == Catch::Approx(50e-6 / 0.32).epsilon(1e-12));
    CHECK(widest_binding_gap(50e-6, 0.4e-3) < (0.4e-3 / 0.8) * (0.4e-3 / 0.8) / (4.0 * 50e-6));
    // a round wire at its own 1.2 r target never reaches that cap: 0.5625 r
    CHECK(widest_binding_gap(0.25e-3, wire_scale_target(0.25e-3)) == Catch::Approx(0.5625 * 0.25e-3).epsilon(1e-12));
    REQUIRE_THROWS_AS(widest_binding_gap(0.0, 1e-3), std::invalid_argument);
    REQUIRE_THROWS_AS(widest_binding_gap(1e-3, 0.0), std::invalid_argument);
    // facing: parallel walls always; a wall inclined by more than target / (2 radius) never; nothing past a right angle
    CHECK(surfaces_face(-1.0, 10.2e-3, 0.4e-3));
    CHECK(surfaces_face(-std::cos(0.019), 10.2e-3, 0.4e-3));         // limit 0.4 / 20.4 = 0.0196 rad
    CHECK_FALSE(surfaces_face(-std::cos(0.020), 10.2e-3, 0.4e-3));
    CHECK_FALSE(surfaces_face(-0.5, 1e-3, 0.4e-3));                  // a 60 degree re-entrant corner, limit 0.2 rad
    CHECK(surfaces_face(-std::cos(1.5), 50e-6, 0.4e-3));             // a small fillet: anything short of a right angle
    CHECK_FALSE(surfaces_face(0.0, 50e-6, 0.4e-3));
    CHECK_FALSE(surfaces_face(1.0, 50e-6, 0.4e-3));
    REQUIRE_THROWS_AS(surfaces_face(-1.0, 0.0, 1e-3), std::invalid_argument);
}

namespace {
// a foil layer about the Z axis: a rectangle [rIn, rOut] x [0, height] revolved through `angle`
TopoDS_Shape foilLayer(double rIn, double rOut, double height, double angle) {
    BRepBuilderAPI_MakePolygon section(gp_Pnt(rIn, 0, 0), gp_Pnt(rOut, 0, 0), gp_Pnt(rOut, 0, height), gp_Pnt(rIn, 0, height), true);
    return BRepPrimAPI_MakeRevol(BRepBuilderAPI_MakeFace(section.Wire(), true).Face(),
                                 gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), angle).Shape();
}
TopoDS_Shape plate(double x0, double y0, double z0, double x1, double y1, double z1) {
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape();
}
}  // namespace

// A rectangular wire or a foil has no radius of its own: its corridors are measured on its built surfaces, with the
// radius of curvature of the face that sees the gap (ABT #2460). The round-tube model it replaces put its refinement
// points inside the copper (two_switch_forward: 191 um inside a foil; 03_buck: 20 of 36 under the skin).
TEST_CASE("Rectangular conductors: corridors come from the curvature and the gaps of their built surfaces",
          "[mesh-api][corridors]") {
    const double mm = 1e-3, um = 1e-6, target = 0.4 * mm;

    SECTION("two foil layers about a column: the bound is the wrap radius's, the point is in the gap") {
        const double rOut = 10.2 * mm, gap = 2.0 * um;   // under the 6.13 um a 10.2 mm wrap can bind at 0.4 mm
        const std::vector<ConductorSolid> foils{{foilLayer(10.0 * mm, rOut, 5 * mm, 4.7), "Secondary parallel 0", "Secondary", true, target},
                                                {foilLayer(rOut + gap, rOut + gap + 0.2 * mm, 5 * mm, 4.7), "Secondary parallel 1", "Secondary", true, target}};
        const SurfaceCorridors sc = surface_corridor_points(foils, {});
        CHECK(sc.rectangularSolids == 2);
        CHECK(sc.curvedFaces == 4);
        REQUIRE(sc.points.size() > 1000);
        CHECK(sc.ownHits == 0);
        CHECK(sc.inclined == 0);
        CHECK(sc.roundHits == 0);
        for (const auto& c : sc.points) {
            const double r = std::hypot(c.p[0], c.p[1]);
            REQUIRE(r == Catch::Approx(rOut + 0.5 * gap).margin(1e-9));          // mid-gap: in air, never in copper
            REQUIRE(c.gap == Catch::Approx(gap).margin(1e-9));
            // seen from the inner layer's outer face (10.2 mm) or the outer layer's inner face (10.202 mm)
            REQUIRE(c.h >= chord_bound(rOut, c.gap) * (1.0 - 1e-9));
            REQUIRE(c.h <= chord_bound(rOut + gap, c.gap) * (1.0 + 1e-9));
            REQUIRE(c.h < target);
        }
        // not the round-tube bound: half the 0.2 mm thickness would ask for 23 um here, the wrap radius asks for 229 um
        CHECK(sc.points.front().h > 9.0 * chord_bound(0.1 * mm, gap));
    }
    SECTION("the same layers 10 um apart bind nothing: the wrap is too gently curved") {
        const double rOut = 10.2 * mm, gap = 10.0 * um;
        const std::vector<ConductorSolid> foils{{foilLayer(10.0 * mm, rOut, 5 * mm, 4.7), "Secondary parallel 0", "Secondary", true, target},
                                                {foilLayer(rOut + gap, rOut + gap + 0.2 * mm, 5 * mm, 4.7), "Secondary parallel 1", "Secondary", true, target}};
        const SurfaceCorridors sc = surface_corridor_points(foils, {});
        CHECK(sc.samples > 1000);   // the faces were measured
        CHECK(sc.points.empty());
    }
    SECTION("two turns of ONE foil body are seen as a gap between turns") {
        const double rOut = 10.2 * mm, gap = 2.0 * um;
        BRepAlgoAPI_Fuse bridge(foilLayer(10.0 * mm, rOut, 5 * mm, 4.7), plate(10.1 * mm, -0.3 * mm, 0, 10.3 * mm, 0, 5 * mm));
        REQUIRE(bridge.IsDone());
        BRepAlgoAPI_Fuse body(bridge.Shape(), foilLayer(rOut + gap, rOut + gap + 0.2 * mm, 5 * mm, 4.7));
        REQUIRE(body.IsDone());
        size_t solids = 0;
        for (TopExp_Explorer ex(body.Shape(), TopAbs_SOLID); ex.More(); ex.Next()) ++solids;
        REQUIRE(solids == 1);
        const SurfaceCorridors sc = surface_corridor_points({{body.Shape(), "Secondary parallel 0", "Secondary", true, target}}, {});
        REQUIRE(sc.points.size() > 1000);
        CHECK(sc.ownHits == sc.points.size());
        for (const auto& c : sc.points) REQUIRE(std::hypot(c.p[0], c.p[1]) == Catch::Approx(rOut + 0.5 * gap).margin(1e-9));
    }
    SECTION("flat faces 2 um apart have no sag and ask for nothing") {
        const std::vector<ConductorSolid> plates{{plate(0, 0, 0, 5 * mm, 5 * mm, 0.2 * mm), "A parallel 0", "A", true, target},
                                                 {plate(0, 0, 0.2 * mm + 2 * um, 5 * mm, 5 * mm, 0.4 * mm + 2 * um), "B parallel 0", "B", true, target}};
        const SurfaceCorridors sc = surface_corridor_points(plates, {});
        CHECK(sc.rectangularSolids == 2);
        CHECK(sc.curvedFaces == 0);
        CHECK(sc.samples == 0);
        CHECK(sc.points.empty());
    }
    SECTION("a re-entrant corner of one body is not a gap") {
        // a 2 mm lead sunk half its radius into a sheet, the body meshed at 2 mm: where they meet, the lead's surface
        // leaves the sheet at 60 degrees. The normals of the lead's first row of samples reach the sheet 0.85 mm away --
        // within the 1.56 mm a 1 mm radius can bind at that target, and the bound there (1.47 mm) is under the target --
        // but a wall inclined by 74 degrees is past the 57 degrees a 2 mm facet turns through: a corner, not a corridor.
        const double coarse = 2.0 * mm;
        BRepAlgoAPI_Fuse body(plate(-3 * mm, -3 * mm, -0.2 * mm, 3 * mm, 3 * mm, 0),
                              BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-2 * mm, 0, 0.5 * mm), gp_Dir(1, 0, 0)), 1.0 * mm, 4 * mm).Shape());
        REQUIRE(body.IsDone());
        const SurfaceCorridors sc = surface_corridor_points({{body.Shape(), "Secondary parallel 0", "Secondary", true, coarse}}, {});
        CHECK(sc.curvedFaces >= 1);
        REQUIRE(sc.hits > 0);               // normals near the junction do reach the sheet
        CHECK(sc.inclined == sc.hits);      // and every one of them meets it inclined
        CHECK(sc.points.empty());
        CHECK(chord_bound(1.0 * mm, 0.848 * mm) < coarse);   // the bound alone would have asked for refinement there
    }
    SECTION("a round wire over a foil: the wire's radius bounds the gap, measured to the foil's real face") {
        const double r = 0.25 * mm, gap = 5.0 * um, wireTarget = wire_scale_target(r);
        const std::vector<ConductorSolid> foil{{plate(0, 0, -0.2 * mm, 10 * mm, 5 * mm, 0), "Secondary parallel 0", "Secondary", true, target}};
        // along x over the foil's broad face, 3 mm in from its edges; and the same wire 1 mm above
        const std::vector<RoundCentreline> near{{"Primary parallel 0", r, wireTarget, {{{3 * mm, 2.5 * mm, r + gap}, {7 * mm, 2.5 * mm, r + gap}}}}};
        const std::vector<RoundCentreline> far{{"Primary parallel 0", r, wireTarget, {{{3 * mm, 2.5 * mm, r + 1 * mm}, {7 * mm, 2.5 * mm, r + 1 * mm}}}}};
        const SurfaceCorridors sc = surface_corridor_points(foil, near);
        CHECK(sc.curvedFaces == 0);
        REQUIRE(sc.roundHits == 32);        // 4 mm of wire sampled every half radius
        REQUIRE(sc.points.size() == 32);
        for (const auto& c : sc.points) {
            REQUIRE(c.p[2] == Catch::Approx(0.5 * gap).margin(1e-10));   // mid-gap, under the wire
            REQUIRE(c.p[1] == Catch::Approx(2.5 * mm).margin(1e-10));
            REQUIRE(c.gap == Catch::Approx(gap).margin(1e-10));
            REQUIRE(c.h == Catch::Approx(chord_bound(r, gap)).epsilon(1e-6));
        }
        CHECK(surface_corridor_points(foil, far).points.empty());
        // the foil's OWN lead is part of its body: its centreline lies a radius from its own surface, and that is not a gap
        const std::vector<RoundCentreline> ownLead{{"Secondary parallel 0 exit lead", r, target, {{{3 * mm, 2.5 * mm, r + gap}, {7 * mm, 2.5 * mm, r + gap}}}}};
        const SurfaceCorridors own = surface_corridor_points(foil, ownLead);
        CHECK(own.roundSamples == 0);
        CHECK(own.points.empty());
        // ... a lead of another parallel of the same winding is not
        const std::vector<RoundCentreline> otherLead{{"Secondary parallel 10 exit lead", r, target, {{{3 * mm, 2.5 * mm, r + gap}, {7 * mm, 2.5 * mm, r + gap}}}}};
        CHECK(surface_corridor_points(foil, otherLead).points.size() == 32);
    }
    SECTION("no rectangular conductor, no work; no target, no size") {
        const TopoDS_Shape wire = BRepPrimAPI_MakeCylinder(0.25 * mm, 5 * mm).Shape();
        const SurfaceCorridors none = surface_corridor_points({{wire, "Primary parallel 0", "Primary", false, target}}, {});
        CHECK(none.rectangularSolids == 0);
        CHECK(none.samples == 0);
        CHECK(none.points.empty());
        REQUIRE_THROWS_WITH(surface_corridor_points({{plate(0, 0, 0, 1 * mm, 1 * mm, 1 * mm), "Secondary parallel 0", "Secondary", true, 0.0}}, {}),
                            ContainsSubstring("no positive target"));
        REQUIRE_THROWS_WITH(surface_corridor_points({{plate(0, 0, 0, 1 * mm, 1 * mm, 1 * mm), "Secondary parallel 0", "Secondary", true, target}},
                                                    {{"Primary parallel 0", 0.0, target, {}}}),
                            ContainsSubstring("no positive radius"));
    }
}

// meshMagnetic takes every knob from its recipe: an ambient OMFEM_*/MVB_* is refused before any
// work, and a recipe the table rejects leaves no knob behind in the environment.
TEST_CASE("The outer box faces' surface algorithm meshes a fine zone that lies strictly inside a coarse face",
          "[mesh-api][sizing]") {
    using namespace mvb::mesh;
    // 128 mm2 at 0.67 mm is about 330 nodes; the algorithm the mesher gives its outer faces has to find them.
    const std::size_t found = nodes_in_enclosed_fine_zone(outer_face_algorithm_2d());
    INFO("nodes in the zone with the outer faces' algorithm: " << found);
    CHECK(found > 200);
    // What the choice rests on: gmsh's default (6, Frontal-Delaunay) never sees the zone -- 0 nodes of its own in
    // it (measured 2026-11-22, fork om-4.15.2). If a later gmsh fixes that, this fails and the choice can be
    // revisited; until then an outer face left on the default is isolated_buck's 180 tets below gamma 0.01.
    const std::size_t missed = nodes_in_enclosed_fine_zone(6);
    INFO("nodes in the zone with Frontal-Delaunay: " << missed);
    CHECK(missed < 20);
}

TEST_CASE("A faceted wire's surface algorithm leaves no flat wedge under the ridges of a bend",
          "[mesh-api][sizing]") {
    using namespace mvb::mesh;
    const std::size_t flat = flat_tets_in_faceted_bend(faceted_wire_face_algorithm_2d());
    INFO("tets below gamma 0.1 with the faceted-wire algorithm: " << flat);
    CHECK(flat == 0);
    // What the choice rests on: gmsh's default (6, Frontal-Delaunay) builds them on the same bend (5, measured
    // 2026-11-22, fork om-4.15.2, one thread). If a later gmsh stops, this fails and the choice can be revisited.
    const std::size_t withDefault = flat_tets_in_faceted_bend(6);
    INFO("tets below gamma 0.1 with Frontal-Delaunay: " << withDefault);
    CHECK(withDefault >= 3);
}

TEST_CASE("meshMagnetic refuses knobs from the environment and leaves none behind", "[mesh-api]") {
    const nlohmann::json magnetic = {{"core", nlohmann::json::object()}, {"coil", nlohmann::json::object()}};
    for (const char* var : {"OMFEM_NO_BOBBIN", "MVB_WELD_ALL"}) {
        setenv(var, "1", 1);
        REQUIRE_THROWS_WITH(meshMagnetic(magnetic, nlohmann::json::object()), ContainsSubstring(var));
        unsetenv(var);
    }
    // A known knob next to an unknown one: the table sets the first, then throws on the second.
    const nlohmann::json bad = {{"air", {{"target_m", 0.006}}}, {"air_typo", {{"x", 1}}}};
    REQUIRE_THROWS_WITH(meshMagnetic(magnetic, bad), ContainsSubstring("unknown key"));
    const nlohmann::json emptyTypo = {{"air_typo", nlohmann::json::object()}};   // an unknown section with nothing in it
    REQUIRE_THROWS_WITH(meshMagnetic(magnetic, emptyTypo), ContainsSubstring("unknown key(s): air_typo"));
    REQUIRE(std::getenv("OMFEM_AIR_TARGET") == nullptr);
    REQUIRE_THROWS_WITH(meshMagnetic(nlohmann::json::array(), nlohmann::json::object()),
                        ContainsSubstring("expects a MAS magnetic"));
}

// The mesher opens and closes its own gmsh session; nothing it leaves behind may call gmsh after
// the close (gmsh reports "Gmsh has not been initialized" on stderr, ABT #1659).
TEST_CASE("meshMagnetic refuses a mesh holding a tet below the recipe's quality floor", "[mesh-api]") {
    std::ifstream f(MAS_EXAMPLES_DIR "/000_debug.json");
    REQUIRE(f);
    const auto mas = nlohmann::json::parse(f);
    // No real mesh is regular throughout: a floor of 0.9 is below every mesh's worst tet. One attempt, so the
    // refusal is the mesh's own and not the last of a ladder of retries.
    const nlohmann::json strict = {{"mesher", {{"quality_floor", 0.9}, {"retries", 0}}}};
    REQUIRE_THROWS_WITH(meshMagnetic(mas.at("magnetic"), strict),
                        ContainsSubstring("below the quality floor 0.9") && ContainsSubstring("the worst has gamma"));
    const nlohmann::json negative = {{"mesher", {{"quality_floor", -1.0}, {"retries", 0}}}};
    REQUIRE_THROWS_WITH(meshMagnetic(mas.at("magnetic"), negative), ContainsSubstring("mesher.quality_floor must be a number >= 0"));
}

TEST_CASE("meshMagnetic meshes again in the same process after a mesh and after a refusal", "[mesh-api]") {
    std::ifstream f(MAS_EXAMPLES_DIR "/000_debug.json");
    REQUIRE(f);
    const auto mas = nlohmann::json::parse(f);
    const auto tetsOf = [](const MeshResult& r) {
        std::size_t n = 0;
        for (const auto& b : r.mesh.blocks)
            if (b.type == ElementType::Tet4) n += b.ids.size();
        return n;
    };
    // The mesher pinned the lead tips to the port plane through an environment variable it never took back, so the
    // second call found "its caller's" environment set and refused to start. Same recipe, same mesh, both times;
    // and a refusal in between (a mesh below an impossible quality floor) leaves nothing behind either.
    const std::size_t first = tetsOf(meshMagnetic(mas.at("magnetic"), nlohmann::json::object()));
    REQUIRE(first > 0);
    const nlohmann::json strict = {{"mesher", {{"quality_floor", 0.9}, {"retries", 0}}}};
    REQUIRE_THROWS_WITH(meshMagnetic(mas.at("magnetic"), strict), ContainsSubstring("below the quality floor"));
    const std::size_t second = tetsOf(meshMagnetic(mas.at("magnetic"), nlohmann::json::object()));
    CHECK(second == first);
}

TEST_CASE("meshMagnetic calls no gmsh after closing the session it opened", "[mesh-api]") {
    std::ifstream f(MAS_EXAMPLES_DIR "/000_debug.json");
    REQUIRE(f);
    const auto mas = nlohmann::json::parse(f);
    char path[] = "/tmp/mvbpp_mesh_stderr_XXXXXX";
    const int fd = mkstemp(path);
    REQUIRE(fd >= 0);
    std::fflush(stderr);
    const int saved = dup(2);
    dup2(fd, 2);
    std::size_t tets = 0;
    try {
        const auto r = meshMagnetic(mas.at("magnetic"), nlohmann::json::object());
        for (const auto& b : r.mesh.blocks)
            if (b.type == ElementType::Tet4) tets += b.ids.size();
    } catch (...) {
        std::fflush(stderr);
        dup2(saved, 2);
        throw;
    }
    std::fflush(stderr);
    dup2(saved, 2);
    close(saved);
    close(fd);
    std::ifstream e(path);
    const std::string log((std::istreambuf_iterator<char>(e)), std::istreambuf_iterator<char>());
    std::remove(path);
    REQUIRE(tets > 0);
    REQUIRE(log.find("[mesh3d]") != std::string::npos);   // the capture saw the mesher's own log
    REQUIRE(log.find("has not been initialized") == std::string::npos);
}

TEST_CASE("Abaqus INP and VTK round-trip a mesh exactly, ids, regions and entities included", "[mesh-io]") {
    const Mesh m = cell();
    for (const char* format : {"inp", "vtk"}) {
        INFO("format " << format);
        REQUIRE(diffMeshes(m, importMesh(exportMesh(m, format, "m"), format)) == "");
        // In millimetres the coordinates are scaled and back: equal to the last ulp or so.
        const Mesh mm = importMesh(exportMesh(m, format, "mm"), format);
        for (std::size_t i = 0; i < m.xyz.size(); ++i)
            REQUIRE(std::abs(m.xyz[i] - mm.xyz[i]) <= 4e-16 * std::max(std::abs(m.xyz[i]), 1e-9));
        Mesh same = mm;
        same.xyz = m.xyz;
        REQUIRE(diffMeshes(m, same) == "");
    }
    // The prism is the one cell whose VTK node order differs from gmsh's.
    Mesh p;
    p.nodeIds = {1, 2, 3, 4, 5, 6};
    p.xyz = {0, 0, 0, 1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3, 1e-3, 0, 1e-3, 0, 1e-3, 1e-3};
    p.regions = {Region{1, 3, "prism", {MaterialKind::None, {}}}};
    p.blocks = {{ElementType::Prism6, 1, 1, {7}, {1, 2, 3, 4, 5, 6}}};
    const std::string vtk = exportMesh(p, "vtk", "m");
    REQUIRE_THAT(vtk, ContainsSubstring("CELLS 1 7\n6 0 2 1 3 5 4\n"));
    REQUIRE(diffMeshes(p, importMesh(vtk, "vtk")) == "");
}

TEST_CASE("partsOnly leaves out the air and its box, keeps the parts and their terminals, and adds each part's skin", "[mesh-io]") {
    Mesh m = cell();
    // A terminal on the core: the tet's face {1, 2, 3}.
    m.regions.push_back(Region{5, 2, "term0", {MaterialKind::None, {}}});
    m.blocks.push_back({ElementType::Tri3, 5, 32, {43}, {1, 3, 2}});
    const Mesh p = partsOnly(m);

    std::vector<std::string> names;
    for (const auto& r : p.regions) names.push_back(r.name);
    CHECK(names == std::vector<std::string>{"core", "term0", "core skin"});
    CHECK(p.nodeIds == std::vector<std::int64_t>{1, 2, 3, 4, 7, 8, 9});   // node 5 was the air's alone
    REQUIRE(p.xyz.size() == 3 * p.nodeIds.size());
    CHECK(p.xyz[3 * 4] == m.xyz[3 * 5]);                                  // node 7, coordinates kept exactly
    // The core's elements and the terminal are kept as they were.
    REQUIRE(p.blocks.size() == 5);
    CHECK(p.blocks[0].ids == m.blocks[0].ids);
    CHECK(p.blocks[1].nodes == m.blocks[1].nodes);
    CHECK(p.blocks[2].ids == std::vector<std::int64_t>{43});
    // The skin: the tet's 4 faces and the pyramid's 4 triangles and base (they share none), new ids
    // after the mesh's last, every face's normal pointing out of the element it bounds.
    const auto& tris = p.blocks[3];
    const auto& quads = p.blocks[4];
    CHECK(tris.type == ElementType::Tri3);
    CHECK(quads.type == ElementType::Quad4);
    CHECK(tris.region == p.regions[2].tag);
    CHECK(tris.ids.size() == 8);
    CHECK(quads.ids.size() == 1);
    CHECK(tris.ids.front() == 44);
    auto at = [&](std::int64_t id) {
        for (std::size_t i = 0; i < m.nodeIds.size(); ++i)
            if (m.nodeIds[i] == id) return std::array<double, 3>{m.xyz[3 * i], m.xyz[3 * i + 1], m.xyz[3 * i + 2]};
        FAIL("no node " << id);
        return std::array<double, 3>{};
    };
    auto centre = [&](const std::vector<std::int64_t>& ids) {
        std::array<double, 3> c{0, 0, 0};
        for (auto id : ids)
            for (int k = 0; k < 3; ++k) c[k] += at(id)[k] / double(ids.size());
        return c;
    };
    const auto tet = centre({1, 2, 3, 4}), pyramid = centre({1, 2, 8, 7, 9});
    for (std::size_t e = 0; e < tris.ids.size(); ++e) {
        const std::vector<std::int64_t> f(tris.nodes.begin() + 3 * e, tris.nodes.begin() + 3 * e + 3);
        const bool ofTet = std::all_of(f.begin(), f.end(), [](std::int64_t id) { return id <= 4; });
        const auto a = at(f[0]), b = at(f[1]), c = at(f[2]);
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
        const auto mid = centre(f), owner = ofTet ? tet : pyramid;
        INFO("skin face " << f[0] << " " << f[1] << " " << f[2]);
        CHECK(n[0] * (mid[0] - owner[0]) + n[1] * (mid[1] - owner[1]) + n[2] * (mid[2] - owner[2]) > 0);
    }
    // It writes and reads back like any mesh; the export option is the same thing.
    CHECK(diffMeshes(p, importMesh(exportMesh(p, "msh2", "m"), "msh2")).empty());
    CHECK(exportMesh(m, "msh2", "m", {{"partsOnly", true}}) == exportMesh(p, "msh2", "m"));
    CHECK(exportMesh(m, "bdf", "mm", {{"partsOnly", true}}) == exportMesh(p, "bdf", "mm"));
    CHECK_THROWS_WITH(exportMesh(m, "msh2", "m", {{"partsOnly", 1}}), ContainsSubstring("true or false"));
}

TEST_CASE("partsOnly refuses a mesh with no air to leave out, and a skin name already taken", "[mesh-io]") {
    Mesh noAir = cell();
    noAir.regions[1].name = "gap air";
    CHECK_THROWS_WITH(partsOnly(noAir), ContainsSubstring("no volume region named 'air'"));
    Mesh taken = cell();
    taken.regions[2].name = "core skin";
    CHECK_THROWS_WITH(partsOnly(taken), ContainsSubstring("already has a region named 'core skin'"));
}

namespace {

// The 2-D regions a msh2 file names: its $PhysicalNames entries of dimension 2.
std::set<std::string> regions2d(const std::string& path) {
    std::ifstream f(path);
    REQUIRE(f);
    std::set<std::string> out;
    std::string line;
    bool in = false;
    while (std::getline(f, line)) {
        if (line == "$PhysicalNames") { in = true; std::getline(f, line); continue; }
        if (line == "$EndPhysicalNames") break;
        if (!in) continue;
        std::istringstream s(line);
        int dim = 0, tag = 0;
        std::string name;
        s >> dim >> tag >> name;
        if (dim == 2) out.insert(name.substr(1, name.size() - 2));
    }
    return out;
}

}  // namespace

// ABT #1314 (a): a turn wound with P parallels is P conductors in the 2-D section, one region per
// parallel ("turn_<w>_<t>_par<p>_<leg>"), so the solver can let the turn current split between
// them. The fixture is MKF's one-turn, 16-parallel planar winding (OMFEM's own acceptance case).
TEST_CASE("mesh2d_from_mas gives each parallel of a turn its own region", "[mesh-2d]") {
    std::ifstream f(MAS_COMPLETE_DIR "/planar_one_turn_16_parallels_100k.json");
    REQUIRE(f);
    const auto mas = nlohmann::json::parse(f);
    std::set<std::string> expected;
    for (const auto& t : mas.at("magnetic").at("coil").at("turnsDescription")) {
        REQUIRE(t.at("coordinates")[0].get<double>() > 0.0);   // all on the +x half the section keeps
        expected.insert("turn_Primary_0_par" + std::to_string(t.at("parallel").get<int>()) + "_plus");
    }
    REQUIRE(expected.size() == 16);

    const auto turns = [](const std::set<std::string>& regions) {
        std::set<std::string> out;
        for (const auto& r : regions)
            if (r.rfind("turn_", 0) == 0) out.insert(r);
        return out;
    };

    MeshOptions magnetic;
    magnetic.out_msh = "/tmp/mvbpp_mesh2d_parallels16.msh";
    magnetic.conductor_conducting_diameter = true;
    REQUIRE(mesh2d_from_mas(mas.at("magnetic"), magnetic) == magnetic.out_msh);
    CHECK(turns(regions2d(magnetic.out_msh)) == expected);

    std::remove(magnetic.out_msh.c_str());
}

// Mapped copper cuts a rectangular wire at its section rings; gmsh's addCurveLoop joins the
// four section curves in the order it is given. The ring found on 18_stacked came out sorted by
// tag as width, width, thickness, thickness (a fused junction had renumbered one edge), and gmsh
// refused it: "Curve loop is not closed".
TEST_CASE("walk_ring hands a section ring over end to end, whatever order it was found in", "[mesh-3d]") {
    // 154 and 161: the two 5 mm sides; 166 and 2051: the two 1 mm sides. Corners 1..4.
    const std::map<int, std::vector<int>> ends{{154, {1, 2}}, {161, {3, 4}}, {166, {2, 3}}, {2051, {4, 1}}};
    const auto walk = walk_ring({154, 161, 166, 2051}, ends);
    REQUIRE(walk.size() == 4);
    REQUIRE(std::set<int>(walk.begin(), walk.end()) == std::set<int>{154, 161, 166, 2051});
    const auto shareAnEnd = [&](int a, int b) {
        for (int p : ends.at(a)) for (int q : ends.at(b)) if (p == q) return true;
        return false;
    };
    for (std::size_t k = 0; k < walk.size(); ++k) CHECK(shareAnEnd(walk[k], walk[(k + 1) % walk.size()]));
    // Four curves that touch but do not close are refused, not handed to gmsh.
    const std::map<int, std::vector<int>> open{{1, {1, 2}}, {2, {2, 3}}, {3, {3, 4}}, {4, {4, 5}}};
    REQUIRE_THROWS_WITH(walk_ring({1, 2, 3, 4}, open), ContainsSubstring("does not return"));
}
