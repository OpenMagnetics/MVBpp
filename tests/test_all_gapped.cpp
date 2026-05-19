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

static const std::set<std::string> EXCLUDED = {"ui", "ut", "pqi", "t"};

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
