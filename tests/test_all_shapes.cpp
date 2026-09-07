// Mirrors MVB Python's test_builder.py::test_all_shapes_generated.
// Iterates over MAS/data/core_shapes.ndjson and verifies every shape
// (excluding ui/pqi/ut, same exclusions as MVB Python) builds a non-empty
// solid via the C++ shape factory.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/catch_approx.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "mvb/shapes/ShapeBuilder.h"
#include "MAS.hpp"
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

static const std::set<std::string> EXCLUDED_FAMILIES = {"ui", "pqi", "ut"};

#ifndef MAS_DATA_DIR
#define MAS_DATA_DIR "."
#endif

TEST_CASE("All MAS core shapes build a non-empty solid",
          "[shapes][all_shapes]") {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());

    int total = 0, skipped = 0, failed = 0;
    std::vector<std::string> failures;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json j = json::parse(line, nullptr, false);
        if (j.is_discarded()) continue;

        std::string family = j.value("family", "");
        std::string name = j.value("name", "?");
        if (EXCLUDED_FAMILIES.count(family)) { ++skipped; continue; }

        // Patch null dim subfields (same as the step generator does)
        mvb::patch_dimension_nominals(j);

        MAS::CoreShape shape;
        try {
            shape = j.get<MAS::CoreShape>();
        } catch (const std::exception& e) {
            failures.push_back(name + " (parse: " + e.what() + ")");
            ++failed; ++total;
            continue;
        }

        auto subtype = shape.get_family_subtype().value_or("");
        auto builder = mvb::shapes::createShapeBuilder(shape.get_family(), subtype);
        if (!builder) {
            failures.push_back(name + " (no builder)");
            ++failed; ++total;
            continue;
        }

        TopoDS_Shape piece;
        try {
            piece = builder->buildPiece(shape);
        } catch (const std::exception& e) {
            failures.push_back(name + " (build: " + e.what() + ")");
            ++failed; ++total;
            continue;
        }

        if (piece.IsNull()) {
            failures.push_back(name + " (null)");
            ++failed; ++total;
            continue;
        }

        GProp_GProps props;
        BRepGProp::VolumeProperties(piece, props);
        if (props.Mass() <= 0.0) {
            failures.push_back(name + " (zero volume)");
            ++failed;
        }
        ++total;
    }

    std::cerr << "[all_shapes] total=" << total << " skipped=" << skipped
              << " failed=" << failed << "\n";
    for (const auto& f : failures) std::cerr << "  FAIL: " << f << "\n";

    REQUIRE(total > 800);          // sanity: should hit ~882-952 shapes
    REQUIRE(failures.empty());
}

// ABT #1126. "Non-empty solid" above is a weak test: RM 7LP passed it for years while
// rendering as a mangled polyhedron, because its MAS record had no C and ShapeRM defaulted
// c to 0, extruding the profile with zero half-depth. Two things pin the fix.
//
// First, the data. TDK's RM 7 LP (B65819P) and Ferroxcube's RM7/ILP publish the same core
// to every digit -- l/A 0.52 mm^-1, le 23.5 mm, Ae 45.3 mm^2, Amin 39.6 mm^2, Ve 1060 mm^3 --
// and differ only in the height tolerance. So the two pieces must come out the same size.
// With C absent, RM 7LP's volume was a small fraction of its twin's.
//
// Second, the code. A shape whose profile needs a dimension it has not got must say so
// rather than build something from zero.
TEST_CASE("RM 7LP is the same core as its twin, and a missing dimension is refused",
          "[shapes][abt1126]") {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());

    auto volumeOf = [](json j) {
        mvb::patch_dimension_nominals(j);
        auto shape = j.get<MAS::CoreShape>();
        auto builder = mvb::shapes::createShapeBuilder(shape.get_family(),
                                                      shape.get_family_subtype().value_or(""));
        REQUIRE(builder != nullptr);
        auto piece = builder->buildPiece(shape);
        REQUIRE_FALSE(piece.IsNull());
        GProp_GProps props;
        BRepGProp::VolumeProperties(piece, props);
        return props.Mass();
    };

    json lowProfile, twin;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json j = json::parse(line, nullptr, false);
        if (j.is_discarded()) continue;
        if (j.value("name", "") == "RM 7LP") lowProfile = j;
        if (j.value("name", "") == "RM 7/ILP") twin = j;
    }
    REQUIRE_FALSE(lowProfile.is_null());
    REQUIRE_FALSE(twin.is_null());

    double lowProfileVolume = volumeOf(lowProfile);
    double twinVolume = volumeOf(twin);
    INFO("RM 7LP = " << lowProfileVolume * 1e9 << " mm3, RM 7/ILP = " << twinVolume * 1e9 << " mm3");
    // 5%: the two records differ only in D's upper tolerance, so the pieces are the same
    // size to well within this. The degenerate build was out by far more.
    CHECK(lowProfileVolume == Catch::Approx(twinVolume).epsilon(0.05));

    // Take C away again and the builder must refuse, naming what is missing.
    json crippled = lowProfile;
    crippled["dimensions"].erase("C");
    mvb::patch_dimension_nominals(crippled);
    auto shape = crippled.get<MAS::CoreShape>();
    auto builder = mvb::shapes::createShapeBuilder(shape.get_family(),
                                                  shape.get_family_subtype().value_or(""));
    REQUIRE(builder != nullptr);
    REQUIRE_THROWS_WITH(builder->buildPiece(shape),
                        Catch::Matchers::ContainsSubstring("RM 7LP")
                            && Catch::Matchers::ContainsSubstring("dimension C"));
}

// Mirrors MVB Python's test_get_families: every family in core_shapes.ndjson
// (except ui/pqi) must be routable via the factory.
TEST_CASE("All MAS families are supported by the factory",
          "[shapes][get_families]") {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());

    std::set<std::string> families;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json j = json::parse(line, nullptr, false);
        if (j.is_discarded()) continue;
        families.insert(j.value("family", ""));
    }

    std::vector<std::string> unsupported;
    for (const auto& fam : families) {
        if (fam == "ui" || fam == "pqi") continue;  // MVB Python excludes these
        // Look up via a minimal CoreShape just to test factory routing
        json stub = {{"family", fam}, {"name", "stub"},
                     {"type", "standard"}, {"dimensions", json::object()}};
        MAS::CoreShape s;
        try { s = stub.get<MAS::CoreShape>(); }
        catch (...) { unsupported.push_back(fam + " (parse failed)"); continue; }
        auto b = mvb::shapes::createShapeBuilder(s.get_family(), "");
        if (!b) unsupported.push_back(fam);
    }

    std::cerr << "[get_families] families=" << families.size()
              << " unsupported=" << unsupported.size() << "\n";
    for (const auto& u : unsupported) std::cerr << "  UNSUPPORTED: " << u << "\n";

    REQUIRE(unsupported.empty());
}
