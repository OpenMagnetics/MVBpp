// Mirrors MVB Python's gapping tests in test_builder.py:
//   - test_all_subtractive_gapped_cores_generated
//   - test_all_subtractive_distributed_gapped_cores_generated
//   - test_all_additive_gapped_cores_generated
//
// For each shape in MAS/data/core_shapes.ndjson (excluding ui/ut/pqi/t —
// same exclusions as MVB Python), build a minimal Magnetic with the
// requested gapping spec, let MKF compute geometricalDescription/machining,
// then build the core via MagneticBuilder and verify a positive volume.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Core.h"
#include <nlohmann/json.hpp>
#include <Standard_Failure.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

#ifndef MAS_DATA_DIR
#define MAS_DATA_DIR "."
#endif

// Families whose gaps are STRUCTURAL, not user-specified, so a user gapping template is invalid
// for them by definition -- not a build failure. "drumRing" joined this list when ShapeDrum
// registered the family: before that it was unknown to the factory and skipped as unrecognised, so
// the two DR+SRI shapes were never attempted. Now they are, and MVB++ correctly rejects them with
// "drumRing cores cannot carry user gapping: their two annular clearance gaps are structural,
// derived from A/K/D/F". Classify it alongside the others rather than letting a declared property
// read as a regression.
static const std::set<std::string> EXCLUDED = {"ui", "ut", "pqi", "t", "drumRing"};

namespace {

double total_volume(const std::vector<mvb::NamedShape>& pieces) {
    double total = 0.0;
    for (const auto& ns : pieces) {
        if (ns.shape.IsNull()) continue;
        GProp_GProps p;
        BRepGProp::VolumeProperties(ns.shape, p);
        total += p.Mass();
    }
    return total;
}

struct GappingResult {
    bool ok = false;
    int n_pieces = 0;
    double volume = 0.0;
    std::string error;
};

// Build just the core with the given gapping spec — bypasses the coil to
// match MVB Python's PyMKF.calculate_core_data flow (which is core-only).
GappingResult build_with_gapping(const json& shape, const json& gapping,
                                  int numberStacks = 1) {
    GappingResult r;
    json coreJson;
    coreJson["functionalDescription"] = {
        {"name", "dummy"},
        {"type", "two-piece set"},
        {"material", "N97"},
        {"shape", shape},
        {"gapping", gapping},
        {"numberStacks", numberStacks},
    };

    try {
        // Constructing OpenMagnetics::Core triggers gap processing and
        // populates geometricalDescription with the per-column machining.
        OpenMagnetics::Core core(coreJson);
        // If geometricalDescription is empty (e.g., open-circuit families
        // where MKF couldn't compute a gap layout), build the core without
        // gapping — this matches MVB Python's behaviour of just rendering
        // the un-gapped solid in such cases.
        if (!core.get_geometrical_description() ||
            core.get_geometrical_description()->empty()) {
            json noGap = coreJson;
            noGap["functionalDescription"]["gapping"] = json::array();
            core = OpenMagnetics::Core(noGap);
        }
        mvb::MagneticBuilder builder;
        auto pieces = builder.buildCoreNamed(core);
        r.n_pieces = static_cast<int>(pieces.size());
        r.volume = total_volume(pieces);
        r.ok = (r.n_pieces > 0 && r.volume > 0.0);
        if (!r.ok) r.error = "n_pieces=" + std::to_string(r.n_pieces) +
                              " vol=" + std::to_string(r.volume);
    } catch (const Standard_Failure& e) {
        r.error = std::string("OCCT: ") + e.GetMessageString();
    } catch (const std::exception& e) {
        r.error = e.what();
    } catch (...) {
        r.error = "unknown exception";
    }
    return r;
}

void run_gap_test(const std::string& label, const json& gappingTemplate,
                   int numberStacks = 1) {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());

    int total = 0, skipped = 0;
    std::vector<std::string> failures;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json shape = json::parse(line, nullptr, false);
        if (shape.is_discarded()) continue;
        std::string family = shape.value("family", "");
        std::string name = shape.value("name", "?");
        if (EXCLUDED.count(family)) { ++skipped; continue; }

        mvb::patch_dimension_nominals(shape);

        // Apply the gapping template — MKF computes per-column coordinates
        // during enrichment, so we just hand it the lengths/types.
        auto result = build_with_gapping(shape, gappingTemplate, numberStacks);
        if (!result.ok) failures.push_back(name + ": " + result.error);
        ++total;
    }

    std::cerr << "[" << label << "] total=" << total << " skipped=" << skipped
              << " failed=" << failures.size() << "\n";
    for (const auto& f : failures) std::cerr << "  FAIL: " << f << "\n";

    REQUIRE(total > 400);  // ~461 shapes after excluding T/UI/UT/PQI
    REQUIRE(failures.empty());
}

} // namespace

TEST_CASE("All shapes build with subtractive gapping",
          "[shapes][gapping][subtractive]") {
    json gapping = json::array({
        {{"length", 0.001}, {"type", "subtractive"}},
        {{"length", 0.002}, {"type", "subtractive"}},
        {{"length", 0.0},   {"type", "subtractive"}},
    });
    run_gap_test("subtractive_gapping", gapping, 3);
}

TEST_CASE("All shapes build with subtractive distributed gapping",
          "[shapes][gapping][distributed]") {
    json gapping = json::array({
        {{"length", 0.001},   {"type", "subtractive"}},
        {{"length", 0.0005},  {"type", "subtractive"}},
        {{"length", 0.002},   {"type", "subtractive"}},
        {{"length", 0.00005}, {"type", "residual"}},
        {{"length", 0.00005}, {"type", "residual"}},
    });
    run_gap_test("subtractive_distributed_gapping", gapping, 1);
}

TEST_CASE("All shapes build with additive gapping",
          "[shapes][gapping][additive]") {
    json gapping = json::array({
        {{"length", 0.0001}, {"type", "additive"}},
        {{"length", 0.0001}, {"type", "additive"}},
        {{"length", 0.0001}, {"type", "additive"}},
    });
    run_gap_test("additive_gapping", gapping, 1);
}

// EVERY SUBTRACTIVE GAP MUST CUT. Found 2026-09-05 on the OMFEM 3D pipeline with 00_debug
// (PQ 20/16, one turn, three 5 um residual gaps retyped to subtractive so the FEM sees them):
// ShapeP::applyMachining -- which PQ, RM and PM inherit -- returned the piece untouched for any
// side-column machining ("no outer column to cut"), and MagneticBuilder accepted a cut that
// removed nothing. The core reached the mesher gapped on the centre post only and the 3D
// inductance read +82% against MKF; the fault was chased through the mesher and the solver
// before the geometry was checked. This test builds the same core and asks the only question
// that matters: did each machining remove material?
//   - centre gap alone must remove ~ (column area) x gap
//   - all three must remove strictly more than the centre alone
// and MagneticBuilder now throws on a cut that removes nothing, so a regression fails loudly
// either way.
TEST_CASE("PQ 20/16: a subtractive gap on each leg is machined, not just the centre post",
          "[shapes][gapping][sidegap][abt-omfem-00debug]") {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());
    json pq;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json shape = json::parse(line, nullptr, false);
        if (!shape.is_discarded() && shape.value("name", "") == "PQ 20/16") { pq = shape; break; }
    }
    REQUIRE(!pq.is_null());
    mvb::patch_dimension_nominals(pq);

    const double gap = 5e-6;
    // The gap layout MKF itself emits for an all-legs ground core: one per column, with
    // coordinates (00_debug's gapping, retyped).
    json all = json::array({
        {{"length", gap}, {"type", "subtractive"}, {"coordinates", {0.0, 0.0, 0.0}}, {"shape", "round"},
         {"area", 6.0821e-05}, {"sectionDimensions", {0.0088, 0.0088}},
         {"distanceClosestNormalSurface", 0.0051475}, {"distanceClosestParallelSurface", 0.0046}},
        {{"length", gap}, {"type", "subtractive"}, {"coordinates", {0.010122107, 0.0, 0.0}}, {"shape", "irregular"},
         {"area", 3.1419e-05}, {"sectionDimensions", {0.002244214, 0.014}},
         {"distanceClosestNormalSurface", 0.0051475}, {"distanceClosestParallelSurface", 0.0046}},
        {{"length", gap}, {"type", "subtractive"}, {"coordinates", {-0.010122107, 0.0, 0.0}}, {"shape", "irregular"},
         {"area", 3.1419e-05}, {"sectionDimensions", {0.002244214, 0.014}},
         {"distanceClosestNormalSurface", 0.0051475}, {"distanceClosestParallelSurface", 0.0046}},
    });
    json centreOnly = json::array({all[0]});
    json none = json::array();

    auto vNone = build_with_gapping(pq, none);
    auto vCentre = build_with_gapping(pq, centreOnly);
    auto vAll = build_with_gapping(pq, all);
    INFO("ungapped: " << vNone.error << " centre: " << vCentre.error << " all: " << vAll.error);
    REQUIRE(vNone.ok); REQUIRE(vCentre.ok); REQUIRE(vAll.ok);

    const double removedCentre = vNone.volume - vCentre.volume;   // both halves
    const double removedAll = vNone.volume - vAll.volume;
    INFO("removed by the centre gap: " << removedCentre * 1e9 << " mm3, by all three: " << removedAll * 1e9 << " mm3");
    // centre column: F = 8.8 mm round -> 60.8 mm2 x 5 um (2.5 um per half, two halves) ~ 0.30 mm3;
    // the 12-segment polygon column is a little smaller than the circle.
    CHECK(removedCentre > 0.25e-9);
    CHECK(removedCentre < 0.35e-9);
    // the two legs add ~ 2 x 31.4 mm2 x 5 um ~ 0.31 mm3 -- at least half of that must be there
    CHECK(removedAll - removedCentre > 0.15e-9);
}

// ---------------------------------------------------------------------------------------------
// ABT #1184. "This cut removed nothing" answers two different questions, and only one of them is
// a defect. Both are pinned here, so neither half of the distinction can be lost again.
namespace {

// A core of the named shape with one subtractive gap, taken through MKF so the machining that
// comes back is the real thing rather than something hand-written.
OpenMagnetics::Core one_gap_core(const std::string& shapeName, double gapLength) {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());
    std::string line;
    json shape;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json candidate = json::parse(line, nullptr, false);
        if (candidate.is_discarded()) continue;
        if (candidate.value("name", "") == shapeName) { shape = candidate; break; }
    }
    REQUIRE(!shape.is_null());
    mvb::patch_dimension_nominals(shape);
    json coreJson;
    coreJson["functionalDescription"] = {
        {"name", "dummy"}, {"type", "two-piece set"}, {"material", "N97"}, {"shape", shape},
        {"gapping", json::array({ json{{"length", gapLength}, {"type", "subtractive"}} })},
        {"numberStacks", 1},
    };
    return OpenMagnetics::Core(coreJson);
}

} // namespace

TEST_CASE("A gap tool that misses the column is an error", "[shapes][gapping][machining]") {
    auto core = one_gap_core("E 42/21/20", 0.001);
    auto geometry = core.get_geometrical_description().value();

    mvb::MagneticBuilder builder;
    REQUIRE(builder.buildCoreNamed(core).size() == 2);   // as MKF places it, the gap is cut

    // The same gap with its tool moved 50 mm along the column axis, far outside the piece.
    // Nothing is removed because nothing is there, and that silently-ungapped core is what this
    // guard exists to catch: a FEM run on it reports the gapless inductance with no sign
    // anything is wrong.
    int toolsMoved = 0;
    for (auto& piece : geometry) {
        if (!piece.get_machining()) continue;
        auto machining = *piece.get_machining();
        for (auto& operation : machining) {
            auto coordinates = operation.get_coordinates();
            coordinates[1] = (coordinates[1] >= 0 ? 0.05 : -0.05);
            operation.set_coordinates(coordinates);
            ++toolsMoved;
        }
        piece.set_machining(machining);
    }
    REQUIRE(toolsMoved > 0);

    auto missed = core;
    missed.set_geometrical_description(geometry);
    REQUIRE_THROWS_WITH(builder.buildCoreNamed(missed),
                        Catch::Matchers::ContainsSubstring("core machining removed nothing"));
}

TEST_CASE("A gap already cut by an overlapping neighbour is not a missing gap",
          "[shapes][gapping][machining]") {
    auto core = one_gap_core("E 42/21/20", 0.001);
    mvb::MagneticBuilder builder;
    const double singleCutVolume = total_volume(builder.buildCoreNamed(core));
    REQUIRE(singleCutVolume > 0.0);

    // Two cuts in the same place. The second finds the material gone -- yet the gap it asks for
    // is there, ground by the first. MKF produces exactly this whenever a distributed gapping's
    // gaps are longer than the columnHeight / (numberGaps + 1) spacing it lays them out on,
    // which is why 102 catalogue shapes used to fail the distributed case above.
    auto geometry = core.get_geometrical_description().value();
    int piecesDuplicated = 0;
    for (auto& piece : geometry) {
        if (!piece.get_machining()) continue;
        auto machining = *piece.get_machining();
        machining.push_back(machining.front());
        piece.set_machining(machining);
        ++piecesDuplicated;
    }
    REQUIRE(piecesDuplicated > 0);

    auto overlapping = core;
    overlapping.set_geometrical_description(geometry);
    auto pieces = builder.buildCoreNamed(overlapping);
    REQUIRE(pieces.size() == 2);
    // cutting the same slab twice leaves exactly the geometry cutting it once does
    CHECK(total_volume(pieces) == Catch::Approx(singleCutVolume));
}
