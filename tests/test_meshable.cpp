#include <catch2/catch_test_macros.hpp>

#include "mvb/MeshableModel.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"   // complete OpenMagnetics::Magnetic

#include <fstream>
#include <string>

#ifndef MAS_EXAMPLES_DIR
#define MAS_EXAMPLES_DIR "."
#endif

using namespace mvb;

// ── classify_region: the authoritative decoder of the build-time naming convention ──

TEST_CASE("classify_region: turn name -> winding/parallel/turn", "[meshable]") {
    RegionTag t = classify_region("Primary parallel 0 turn 3");
    CHECK(t.role == Role::Turn);
    CHECK(t.winding == "Primary");
    CHECK(t.parallel == 0);
    CHECK(t.turn == 3);
}

TEST_CASE("classify_region: turn with a multi-word winding name", "[meshable]") {
    RegionTag t = classify_region("Secondary winding parallel 2 turn 11");
    CHECK(t.role == Role::Turn);
    CHECK(t.winding == "Secondary winding");
    CHECK(t.parallel == 2);
    CHECK(t.turn == 11);
}

TEST_CASE("classify_region: turn split by a symmetry/section cut (_<int> suffix)", "[meshable]") {
    RegionTag t = classify_region("Primary parallel 0 turn 3_1");
    CHECK(t.role == Role::Turn);
    CHECK(t.winding == "Primary");
    CHECK(t.turn == 3);
}

TEST_CASE("classify_region: STEP label prefixes/suffixes are stripped", "[meshable]") {
    RegionTag a = classify_region("Shapes/Primary parallel 0 turn 1");
    CHECK(a.role == Role::Turn);
    CHECK(a.winding == "Primary");
    RegionTag b = classify_region("Primary parallel 0 turn 1/2/COMPOUND");
    CHECK(b.role == Role::Turn);
    CHECK(b.turn == 1);
}

TEST_CASE("classify_region: bobbin and FR4 board are Bobbin", "[meshable]") {
    CHECK(classify_region("Bobbin").role == Role::Bobbin);
    CHECK(classify_region("FR4Board").role == Role::Bobbin);
}

TEST_CASE("classify_region: core pieces (single and indexed)", "[meshable]") {
    RegionTag a = classify_region("ETD 34/17/11");   // a core shape name with slashes
    // first '/'-segment is taken -> "ETD 34" -> not a turn/bobbin -> Core
    CHECK(a.role == Role::Core);
    RegionTag b = classify_region("Core_1");
    CHECK(b.role == Role::Core);
    CHECK(b.corePiece == 1);
    RegionTag c = classify_region("Core");
    CHECK(c.role == Role::Core);
    CHECK(c.corePiece == 0);
}

// ── buildMeshable smoke test (enrich -> build -> tag) ───────────────────────────────

static MAS::Magnetic load_mas_example(const std::string& filename) {
    std::string path = std::string(MAS_EXAMPLES_DIR) + "/" + filename;
    std::ifstream f(path);
    REQUIRE(f.good());
    nlohmann::json j; f >> j;
    mvb::patch_dimension_nominals(j);                 // examples need nominal patching
    return j.at("magnetic").get<MAS::Magnetic>();     // files are full MAS docs
}

TEST_CASE("buildMeshable: tags a real inductor's regions", "[meshable][slow]") {
    OpenMagnetics::Magnetic enriched =
        magnetic_autocomplete_safe(load_mas_example("01_simple_inductor_etd34_n87.json"));

    MeshableOptions opt;            // full model, outer footprint
    MeshableModel m = buildMeshable(enriched, opt);

    REQUIRE(!m.regions.empty());
    CHECK(m.multiplicity == 1);
    CHECK(m.symmetry.empty());

    int cores = 0, turns = 0;
    for (const auto& r : m.regions) {
        if (r.tag.role == Role::Core) ++cores;
        if (r.tag.role == Role::Turn) {
            ++turns;
            CHECK_FALSE(r.tag.winding.empty());   // every turn must name its winding
            CHECK(r.tag.turn >= 0);
        }
    }
    CHECK(cores >= 1);
    CHECK(turns >= 1);
}

TEST_CASE("buildMeshable: half symmetry reports a plane + multiplicity 2", "[meshable][slow]") {
    OpenMagnetics::Magnetic enriched =
        magnetic_autocomplete_safe(load_mas_example("01_simple_inductor_etd34_n87.json"));

    MeshableOptions opt;
    opt.symmetryPlanes = 1;
    MeshableModel m = buildMeshable(enriched, opt);

    REQUIRE(!m.regions.empty());
    REQUIRE(m.symmetry.size() == 1);
    CHECK(m.multiplicity == 2);
    CHECK(m.symmetry[0].bc == WallBC::FluxTangential);
    CHECK(m.symmetry[0].plane != SymmetryPlane::None);
}
