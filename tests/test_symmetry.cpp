#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Symmetry.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include <nlohmann/json.hpp>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <fstream>
#include <string>
#include <unordered_map>

using Catch::Matchers::WithinRel;
using json = nlohmann::json;

#ifndef MAS_EXAMPLES_DIR
#define MAS_EXAMPLES_DIR "."
#endif

// ── helpers ──────────────────────────────────────────────────────────────────

static double total_volume(const std::vector<mvb::NamedShape>& shapes) {
    double v = 0.0;
    for (const auto& ns : shapes) {
        if (ns.shape.IsNull()) continue;
        GProp_GProps props;
        BRepGProp::VolumeProperties(ns.shape, props);
        v += props.Mass();
    }
    return v;
}

static MAS::Magnetic load_mas_example(const std::string& filename) {
    std::string path = std::string(MAS_EXAMPLES_DIR) + "/" + filename;
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("Cannot open " + path);
    json j;
    f >> j;
    mvb::patch_dimension_nominals(j);
    return j.at("magnetic").get<MAS::Magnetic>();
}

// Per-fixture cache of the full (sym=0, includeBobbin=true) assembly.
// Boolean ops on ETD34 take ~7s per buildAllNamed call; the 8 tests in this
// file together do 9 ETD34 builds and 2 PQ3230 builds. Caching the full
// build collapses that to one build per fixture; cut variants are then
// produced cheaply via apply_symmetry(cached_full, N).
//
// Test execution is single-threaded under Catch2's default runner, so a
// plain unordered_map needs no synchronisation. The shapes are immutable
// after construction — callers must not mutate the returned vector.
static const std::vector<mvb::NamedShape>&
full_assembly(const std::string& filename) {
    static std::unordered_map<std::string, std::vector<mvb::NamedShape>> cache;
    auto it = cache.find(filename);
    if (it != cache.end()) return it->second;
    mvb::MagneticBuilder builder;
    auto shapes = builder.buildAllNamed(load_mas_example(filename),
                                        /*includeBobbin=*/true,
                                        /*symmetryPlanes=*/0);
    auto [ins, _] = cache.emplace(filename, std::move(shapes));
    return ins->second;
}

// analyze_symmetry is a pure function of the assembly; for the ETD34
// fixture it costs ~12s in boolean Common ops because every winding turn
// straddles all three symmetry planes (an inductor wraps the central
// column). Four of the eight tests in this file need the analyze result
// for the same fixture, so we memoize per filename.
static const mvb::SymmetryResult& cached_symmetry(const std::string& filename) {
    static std::unordered_map<std::string, mvb::SymmetryResult> cache;
    auto it = cache.find(filename);
    if (it != cache.end()) return it->second;
    auto res = mvb::analyze_symmetry(full_assembly(filename));
    auto [ins, _] = cache.emplace(filename, std::move(res));
    return ins->second;
}

// Cut result per (fixture, nPlanes). apply_symmetry runs analyze + per-plane
// cuts; both are deterministic for a given assembly.
static const std::vector<mvb::NamedShape>&
cached_cut(const std::string& filename, int nPlanes) {
    static std::unordered_map<std::string, std::vector<mvb::NamedShape>> cache;
    std::string key = filename + "#" + std::to_string(nPlanes);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    auto shapes = mvb::apply_symmetry(full_assembly(filename), nPlanes);
    auto [ins, _] = cache.emplace(std::move(key), std::move(shapes));
    return ins->second;
}

// ── symmetry detection ────────────────────────────────────────────────────────

TEST_CASE("ETD34 inductor has at least two symmetry planes", "[symmetry]") {
    const auto& shapes = full_assembly("01_simple_inductor_etd34_n87.json");
    REQUIRE_FALSE(shapes.empty());

    const auto& result = cached_symmetry("01_simple_inductor_etd34_n87.json");
    INFO("Valid planes found: " << result.valid_planes.size());
    CHECK(result.valid_planes.size() >= 2);
}

// ── volume ratios via symmetryPlanes parameter ────────────────────────────────

TEST_CASE("ETD34: symmetryPlanes=1 gives half total volume", "[symmetry]") {
    const auto& full = full_assembly("01_simple_inductor_etd34_n87.json");
    double V_full = total_volume(full);
    REQUIRE(V_full > 0.0);

    double V_half = total_volume(cached_cut("01_simple_inductor_etd34_n87.json", 1));
    REQUIRE(V_half > 0.0);

    CHECK_THAT(V_half / V_full, WithinRel(0.5, 0.01));
}

TEST_CASE("ETD34: symmetryPlanes=2 gives quarter total volume", "[symmetry]") {
    const auto& full = full_assembly("01_simple_inductor_etd34_n87.json");
    double V_full = total_volume(full);
    REQUIRE(V_full > 0.0);

    double V_quarter = total_volume(cached_cut("01_simple_inductor_etd34_n87.json", 2));
    REQUIRE(V_quarter > 0.0);

    CHECK_THAT(V_quarter / V_full, WithinRel(0.25, 0.01));
}

TEST_CASE("PQ3230 inductor: symmetryPlanes=1 gives half volume", "[symmetry]") {
    const auto& full = full_assembly("03_buck_inductor_pq3230_n95.json");
    double V_full = total_volume(full);
    REQUIRE(V_full > 0.0);

    double V_half = total_volume(cached_cut("03_buck_inductor_pq3230_n95.json", 1));
    REQUIRE(V_half > 0.0);

    CHECK_THAT(V_half / V_full, WithinRel(0.5, 0.01));
}

TEST_CASE("Cut shapes are all non-null and non-zero volume", "[symmetry]") {
    const auto& shapes = cached_cut("01_simple_inductor_etd34_n87.json", 1);
    REQUIRE_FALSE(shapes.empty());
    for (const auto& ns : shapes) {
        CHECK_FALSE(ns.shape.IsNull());
        GProp_GProps props;
        BRepGProp::VolumeProperties(ns.shape, props);
        CHECK(props.Mass() > 0.0);
    }
}

TEST_CASE("symmetryPlanes=0 returns full assembly", "[symmetry]") {
    const auto& shapes = full_assembly("01_simple_inductor_etd34_n87.json");
    // At minimum: 2 core halves + 1 bobbin
    CHECK(shapes.size() >= 3);
}

// ── direct API: cut_to_region ─────────────────────────────────────────────────

TEST_CASE("Explicit cut_to_region gives half volume", "[symmetry]") {
    const auto& full = full_assembly("01_simple_inductor_etd34_n87.json");
    double V_full = total_volume(full);
    REQUIRE(V_full > 0.0);

    const auto& sym = cached_symmetry("01_simple_inductor_etd34_n87.json");
    REQUIRE_FALSE(sym.valid_planes.empty());

    auto bbox = mvb::aggregate_bbox(full);
    std::vector<std::pair<mvb::SymmetryPlane, mvb::SymmetryHalf>> cuts = {
        {sym.valid_planes[0], mvb::SymmetryHalf::Positive}
    };
    auto cut = mvb::cut_to_region(full, cuts, bbox);
    double V_cut = total_volume(cut);

    CHECK_THAT(V_cut / V_full, WithinRel(0.5, 0.01));
}

TEST_CASE("Explicit cut_to_region two planes gives quarter volume", "[symmetry]") {
    const auto& full = full_assembly("01_simple_inductor_etd34_n87.json");
    double V_full = total_volume(full);
    REQUIRE(V_full > 0.0);

    const auto& sym = cached_symmetry("01_simple_inductor_etd34_n87.json");
    REQUIRE(sym.valid_planes.size() >= 2);

    auto bbox = mvb::aggregate_bbox(full);
    std::vector<std::pair<mvb::SymmetryPlane, mvb::SymmetryHalf>> cuts = {
        {sym.valid_planes[0], mvb::SymmetryHalf::Positive},
        {sym.valid_planes[1], mvb::SymmetryHalf::Positive},
    };
    auto cut = mvb::cut_to_region(full, cuts, bbox);
    double V_cut = total_volume(cut);

    CHECK_THAT(V_cut / V_full, WithinRel(0.25, 0.01));
}
