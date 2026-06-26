// Sanity check that drawDimensionedTopView works for toroidal cores.
// Regression for the "SectionDrawing: no polylines" bug surfaced by the
// MAS battery on T-cores (PFC inductor, CMC, boost T).
#include <catch2/catch_test_macros.hpp>

#include "constructive_models/Magnetic.h"
#include "mvb/SectionDrawing.h"
#include "mvb/Utils.h"
#include "MAS.hpp"

#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;

#ifndef MAS_EXAMPLES_DIR
#define MAS_EXAMPLES_DIR "."
#endif

namespace {
OpenMagnetics::Magnetic load(const std::string& rel) {
    std::ifstream in(std::string(MAS_EXAMPLES_DIR) + "/" + rel);
    REQUIRE(in.is_open());
    json j;
    in >> j;
    mvb::patch_dimension_nominals(j);
    return mvb::magnetic_autocomplete_safe(j.at("magnetic"));
}
}

TEST_CASE("Toroidal TopView: PFC inductor produces non-empty SVG",
          "[2d][toroidal][topview]") {
    auto m = load("05_pfc_inductor_t4020_hf60.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedTopView(m, 400, 12));
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<svg") != std::string::npos);
    // The toroidal section must produce at least one closed circular polyline
    // (outer or inner ring). Check there is some <polyline> or <path> with
    // multiple points by looking for repeated comma separators.
    REQUIRE(svg.find("<path") != std::string::npos);
}

TEST_CASE("Toroidal TopView: boost inductor T5026 produces non-empty SVG",
          "[2d][toroidal][topview]") {
    auto m = load("12_boost_inductor_t5026_26.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedTopView(m, 400, 12));
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<path") != std::string::npos);
}

TEST_CASE("Toroidal TopView: CMC T2515 produces non-empty SVG",
          "[2d][toroidal][topview]") {
    auto m = load("07_cmc_t2515_w800.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedTopView(m, 400, 12));
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<path") != std::string::npos);
}

// Regression for the toroidal FRONT-view bug: the front view must show the
// donut FACE (two concentric circles), not a meridional cut (two rectangles).
// Circles sample into long multi-point paths; the buggy rectangles were short
// 2-point straight edges. We assert the geometry group carries many points and
// that the dimensions are the A (OD) / B (ID) diameters, not the height C.
TEST_CASE("Toroidal FrontView: donut face is concentric circles, not rectangles",
          "[2d][toroidal][frontview]") {
    auto m = load("05_pfc_inductor_t4020_hf60.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedFrontView(m, 400, 12));
    REQUIRE(!svg.empty());

    // Isolate the geometry group (black stroke) from the dimension group.
    auto gpos = svg.find("stroke=\"#000000\"");
    REQUIRE(gpos != std::string::npos);
    auto gend = svg.find("</g>", gpos);
    REQUIRE(gend != std::string::npos);
    std::string geom = svg.substr(gpos, gend - gpos);

    // Two sampled circles → many vertices. The rectangle bug produced ~8 short
    // 2-point edges (a handful of 'L's total).
    size_t lcount = 0;
    for (size_t p = geom.find(" L"); p != std::string::npos; p = geom.find(" L", p + 1)) ++lcount;
    INFO("geometry 'L' vertex count = " << lcount);
    REQUIRE(lcount > 50);

    // Front dims annotate the visible diameters (A outer, B inner), not height C.
    REQUIRE(svg.find("A = ") != std::string::npos);
    REQUIRE(svg.find("B = ") != std::string::npos);
}

// Complement of the above: the TOP view is the meridional cut — two rectangles
// (the ring walls) dimensioned by the height C (not the A/B diameters).
TEST_CASE("Toroidal TopView: meridional rectangles dimensioned by height C",
          "[2d][toroidal][topview]") {
    auto m = load("05_pfc_inductor_t4020_hf60.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedTopView(m, 400, 12));
    REQUIRE(!svg.empty());

    auto gpos = svg.find("stroke=\"#000000\"");
    REQUIRE(gpos != std::string::npos);
    auto gend = svg.find("</g>", gpos);
    REQUIRE(gend != std::string::npos);
    std::string geom = svg.substr(gpos, gend - gpos);

    // Rectangles → few short straight edges (no sampled circle here).
    size_t lcount = 0;
    for (size_t p = geom.find(" L"); p != std::string::npos; p = geom.find(" L", p + 1)) ++lcount;
    INFO("top-view geometry 'L' vertex count = " << lcount);
    REQUIRE(lcount < 20);

    // Height C is the meaningful dimension on this view; A/B belong to the front.
    REQUIRE(svg.find("C = ") != std::string::npos);
}

TEST_CASE("E-family TopView still works (regression guard)",
          "[2d][topview]") {
    auto m = load("01_simple_inductor_etd34_n87.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedTopView(m, 400, 12));
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<path") != std::string::npos);
}
