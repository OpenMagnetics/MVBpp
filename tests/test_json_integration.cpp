#include <catch2/catch_test_macros.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "mvb/SectionDrawing.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs.hxx>

using json = nlohmann::json;

TEST_CASE("Build complete magnetic from concentric_rectangular_column_one_turn.json", "[json][integration][.]" ) {
    std::ifstream f("testData/concentric_rectangular_column_one_turn.json");
    REQUIRE(f.is_open());
    json j;
    try {
        f >> j;
    } catch (const std::exception& e) {
        FAIL("JSON read failed: " << e.what());
    }

    try {
        mvb::patch_dimension_nominals(j);
    } catch (const std::exception& e) {
        FAIL("JSON patch failed: " << e.what());
    }

    MAS::Magnetic magnetic;
    try {
        magnetic = j.at("magnetic").get<MAS::Magnetic>();
    } catch (const std::exception& e) {
        FAIL("MAS parse failed: " << e.what());
    } catch (...) {
        FAIL("MAS parse failed with unknown exception type");
    }

    mvb::MagneticBuilder builder;
    std::vector<mvb::NamedShape> corePieces;
    mvb::NamedShape bobbinNamed;
    std::vector<mvb::NamedShape> turnPieces;
    try {
        corePieces = builder.buildCoreNamed(magnetic.get_core());
        bobbinNamed = builder.buildBobbinNamed(magnetic.get_coil(), magnetic.get_core());
        turnPieces = builder.buildTurnsNamed(magnetic.get_coil(), magnetic.get_core());
    } catch (const std::exception& e) {
        FAIL("Build failed: " << e.what());
    }

    REQUIRE(!corePieces.empty());
    // bobbinNamed.shape may legitimately be null: this fixture has
    // column_thickness=0 and wall_thickness=0, which means MAS declares no
    // bobbin material. Per the no-fallback policy we don't synthesise one.
    REQUIRE(!turnPieces.empty());

    int solids = 0;
    double totalVolume = 0.0;
    for (const auto& ns : corePieces) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(ns.shape, props);
        totalVolume += props.Mass();
        for (TopExp_Explorer exp(ns.shape, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    }
    if (!bobbinNamed.shape.IsNull()) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(bobbinNamed.shape, props);
        totalVolume += props.Mass();
        for (TopExp_Explorer exp(bobbinNamed.shape, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    }
    for (const auto& ns : turnPieces) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(ns.shape, props);
        totalVolume += props.Mass();
        for (TopExp_Explorer exp(ns.shape, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    }

    REQUIRE(totalVolume > 0.0);
    REQUIRE(solids >= 3); // at minimum: 2 core halves + 1 turn
}

// Simulate what Core2DVisualizer.vue passes to drawDimensionedFrontView:
// a magnetic with core-only (no coil), no processedDescription, no geometricalDescription.
// This exercises the magnetic_autocomplete_safe empty-coil path.
TEST_CASE("Core2DVisualizer: drawDimensionedFrontView with core-only E shape", "[2d][empty-coil]") {
    // Minimal JSON that Core2DVisualizer passes: shape by name, material, gapping, no coil
    json magneticJson = R"({
        "core": {
            "functionalDescription": {
                "shape": "E 32/6/20",
                "material": "N87",
                "gapping": [
                    {"type": "subtractive", "length": 0},
                    {"type": "subtractive", "length": 0},
                    {"type": "residual", "length": 0.0001}
                ],
                "numberStacks": 1,
                "type": "two-piece set"
            },
            "geometricalDescription": null,
            "processedDescription": null
        }
    })"_json;

    // Step 1: magnetic_autocomplete_safe should not throw and should set geometricalDescription
    OpenMagnetics::Magnetic enriched;
    REQUIRE_NOTHROW(enriched = mvb::magnetic_autocomplete_safe(magneticJson));
    INFO("geometricalDescription set: " << (enriched.get_core().get_geometrical_description() ? "yes" : "no"));
    REQUIRE(enriched.get_core().get_geometrical_description().has_value());
    REQUIRE(!enriched.get_core().get_geometrical_description()->empty());

    // Step 2: drawDimensionedFrontView should not throw and should return non-empty SVG
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedFrontView(enriched, 400, 12));
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<svg") != std::string::npos);
}

TEST_CASE("Core2DVisualizer: drawDimensionedFrontView with core-only E shape (gapped)", "[2d][empty-coil][gap]") {
    json magneticJson = R"({
        "core": {
            "functionalDescription": {
                "shape": "E 32/6/20",
                "material": "N87",
                "gapping": [
                    {"type": "subtractive", "length": 0.0002},
                    {"type": "subtractive", "length": 0.0002},
                    {"type": "residual", "length": 0.0001}
                ],
                "numberStacks": 1,
                "type": "two-piece set"
            },
            "geometricalDescription": null,
            "processedDescription": null
        }
    })"_json;

    OpenMagnetics::Magnetic enriched;
    REQUIRE_NOTHROW(enriched = mvb::magnetic_autocomplete_safe(magneticJson));
    REQUIRE(enriched.get_core().get_geometrical_description().has_value());

    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedFrontView(enriched, 400, 12));
    REQUIRE(!svg.empty());
}

TEST_CASE("Core3DVisualizer: buildCoreSTL with core-only E shape", "[3d][empty-coil]") {
    json magneticJson = R"({
        "core": {
            "functionalDescription": {
                "shape": "E 32/6/20",
                "material": "N87",
                "gapping": [
                    {"type": "subtractive", "length": 0},
                    {"type": "subtractive", "length": 0},
                    {"type": "residual", "length": 0.0001}
                ],
                "numberStacks": 1,
                "type": "two-piece set"
            },
            "geometricalDescription": null,
            "processedDescription": null
        }
    })"_json;

    OpenMagnetics::Magnetic enriched;
    REQUIRE_NOTHROW(enriched = mvb::magnetic_autocomplete_safe(magneticJson));
    REQUIRE(enriched.get_core().get_geometrical_description().has_value());

    mvb::MagneticBuilder builder;
    std::vector<mvb::NamedShape> pieces;
    REQUIRE_NOTHROW(pieces = builder.buildCoreNamed(enriched.get_core()));
    REQUIRE(!pieces.empty());
}
