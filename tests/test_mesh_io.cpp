// Round-trip gate for the mesh writers (ABT #1588): every format is read back and compared with
// the mesh it was written from. The mesh here is a hand-made cell, not a magnetic: it pins the
// format rules (ids, order, region names, units, refusals); the real meshes are round-tripped by
// tools/mvbpp_mesh_roundtrip on OMFEM's reference outputs.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "mvb/mesh/Mesher.h"

#include <cmath>

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
    // A conductor: MAS has no specific heat field for wire materials at all.
    m = cell();
    m.regions[0].material = {MaterialKind::Wire, nlohmann::json{
        {"name", "copper"}, {"density", 8890.0}, {"thermalConductivity", {{{"temperature", 0}, {"value", 401}}}}}};
    REQUIRE_THROWS_WITH(exportMesh(m, "bdf", "m"), ContainsSubstring("specific heat"));
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

// meshMagnetic takes every knob from its recipe: an ambient OMFEM_*/MVB_* is refused before any
// work, and a recipe the table rejects leaves no knob behind in the environment.
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
    REQUIRE(std::getenv("OMFEM_AIR_TARGET") == nullptr);
    REQUIRE_THROWS_WITH(meshMagnetic(nlohmann::json::array(), nlohmann::json::object()),
                        ContainsSubstring("expects a MAS magnetic"));
}
