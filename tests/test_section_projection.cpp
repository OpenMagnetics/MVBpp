// 2D-projection (cut2DFaces) element-count tests.
//
// These guard the section topology that downstream 2D FE consumers depend on: a solid with a
// through-hole must project to ONE face with an empty hole (not overlapping disks), and a
// toroid's turns must each cross the ring's mid-plane exactly twice (inner + outer). Both were
// silently wrong: the toroid annulus came out as two disks (outer + a plug over the bore), and
// the unfused toroidal turns sectioned into open wires (zero closed faces) or, once sectioned
// per sub-solid, into doubled coincident faces at the sub-piece joints.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/SectionBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include <nlohmann/json.hpp>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <numbers>
#include <string>

using json = nlohmann::json;
using Catch::Matchers::WithinRel;

#ifndef MAS_EXAMPLES_DIR
#define MAS_EXAMPLES_DIR "."
#endif

namespace {

MAS::Magnetic load_magnetic(const std::string& file) {
    std::ifstream in(std::string(MAS_EXAMPLES_DIR) + "/" + file);
    REQUIRE(in.is_open());
    json j; in >> j;
    mvb::patch_dimension_nominals(j);
    return j.at("magnetic").get<MAS::Magnetic>();
}

bool name_has(const std::string& n, const char* sub) {
    std::string a = n, b = sub;
    std::transform(a.begin(), a.end(), a.begin(), [](unsigned char c){ return std::tolower(c); });
    std::transform(b.begin(), b.end(), b.begin(), [](unsigned char c){ return std::tolower(c); });
    return a.find(b) != std::string::npos;
}

double face_area(const TopoDS_Shape& s) {
    GProp_GProps g; BRepGProp::SurfaceProperties(s, g); return std::abs(g.Mass());
}

}  // namespace

// A toroid lies flat in XY; its XY mid-plane section is the ring's full magnetic-path cross
// section. The annulus must be ONE face (empty bore), and each turn threads the ring so it
// crosses the plane exactly twice.
TEST_CASE("Toroid XY section: one annular core face, two crossings per turn",
          "[section][toroid]") {
    auto magnetic = load_magnetic("12_boost_inductor_t5026_26.json");   // T 50/30/19, 52 turns
    mvb::MagneticBuilder builder;
    // True-circle core (0 core segments) so the annulus area is exact.
    auto solids = builder.buildAllNamed(magnetic, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                        /*wireSegs=*/0, /*coreSegs=*/0, /*paintCoating=*/false);
    REQUIRE_FALSE(solids.empty());

    int num_turn_solids = 0;
    for (const auto& ns : solids) if (name_has(ns.name, "turn")) ++num_turn_solids;
    REQUIRE(num_turn_solids == 52);

    auto faces = mvb::SectionBuilder::cut2DFaces(solids, mvb::parseSectionPlane("XY"), 0.0);

    int n_core = 0, n_turn = 0;
    double core_area = 0.0;
    for (const auto& ns : faces) {
        if (name_has(ns.name, "turn")) { ++n_turn; }
        else if (name_has(ns.name, "core")) { ++n_core; core_area += face_area(ns.shape); }
    }

    // The donut projects as ONE annular face — not an outer disk plus an inner disk plugging
    // the bore.
    REQUIRE(n_core == 1);

    // Every turn threads the ring -> exactly two crossings (inner hole-side + outer).
    REQUIRE(n_turn == 2 * num_turn_solids);

    // Area = pi/4 (OD^2 - ID^2), T 50/30/19: OD=50mm, ID=30mm -> ~1257 mm^2; crucially NOT the
    // outer disk (pi/4 OD^2 ~ 1963 mm^2), which is what the two-disk bug produced.
    const double expected = std::numbers::pi / 4.0 * (0.050 * 0.050 - 0.030 * 0.030);
    REQUIRE_THAT(core_area, WithinRel(expected, 0.05));
}

// Regression guard for the per-sub-solid + dedup change on a NON-toroid core: a gapped ETD
// inductor must still section to a sane, connected core (no spurious extra or missing faces).
TEST_CASE("Non-toroid XY section keeps a connected core and its turns", "[section][etd]") {
    auto magnetic = load_magnetic("01_simple_inductor_etd34_n87.json");
    mvb::MagneticBuilder builder;
    auto solids = builder.buildAllNamed(magnetic, /*includeBobbin=*/false, /*symmetryPlanes=*/0);
    REQUIRE_FALSE(solids.empty());

    auto faces = mvb::SectionBuilder::cut2DFaces(solids, mvb::parseSectionPlane("XY"), 0.0);

    int n_core = 0, n_turn = 0;
    for (const auto& ns : faces) {
        if (name_has(ns.name, "turn")) ++n_turn;
        else if (name_has(ns.name, "core")) ++n_core;
    }
    // A 1 mm-gapped ETD centre leg splits into a small number of core pieces (not zero, not
    // dozens); turns are present.
    REQUIRE(n_core >= 1);
    REQUIRE(n_core <= 4);
    REQUIRE(n_turn > 0);
}
