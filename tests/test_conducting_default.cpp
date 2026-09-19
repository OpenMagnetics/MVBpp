// ABT #1261: conductors are drawn at their CONDUCTING diameter by default; the coated envelope
// is an opt-in (paintCoating=true, used by the web viewer).
//
// The default is asserted by BEHAVIOUR, not by reading a default argument: a build with no
// coating argument must put exactly as much copper in the model as an explicit
// paintCoating=false build, and strictly less than a paintCoating=true one. Flip any of the
// defaults back to true and the first comparison fails.
//
//   ./mvb_tests "[abt1261]"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "constructive_models/Magnetic.h"
#include "mvb/ConductorBuilder.h"
#include "mvb/MagneticBuilder.h"
#include "mvb/NamedShape.h"
#include "mvb/Utils.h"

#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#ifndef MAS_COMPLETE_DIR
#define MAS_COMPLETE_DIR MAS_EXAMPLES_DIR "/complete"
#endif

namespace {

using json = nlohmann::json;
using Catch::Matchers::WithinRel;

json load_magnetic(const std::string& fixture) {
    const std::filesystem::path path = std::filesystem::path(MAS_COMPLETE_DIR) / fixture;
    REQUIRE(std::filesystem::exists(path));
    std::ifstream in(path);
    json mas; in >> mas;
    return mas.contains("magnetic") ? mas.at("magnetic") : mas;
}

// Total volume of the conductor solids (Role::Turn) of an assembly [m^3].
double copper_volume(const std::vector<mvb::NamedShape>& all) {
    double v = 0.0;
    int n = 0;
    for (const auto& ns : all) {
        if (ns.role != mvb::Role::Turn) continue;
        GProp_GProps props;
        BRepGProp::VolumeProperties(ns.shape, props);
        v += props.Mass();
        ++n;
    }
    REQUIRE(n > 0);
    return v;
}

}  // namespace

TEST_CASE("Option structs default to the conducting diameter", "[abt1261]") {
    CHECK_FALSE(mvb::DrawConfig{}.paintCoating);
    CHECK_FALSE(mvb::ConductorBuilder::Options{}.paintCoating);
}

TEST_CASE("Per-turn build: the default draws the conducting diameter", "[abt1261]") {
    const auto enriched = mvb::magnetic_autocomplete_safe(load_magnetic("flyback_transformer_complete.json"));
    mvb::MagneticBuilder builder;

    const double byDefault = copper_volume(builder.buildAllNamed(enriched, /*includeBobbin=*/false));
    const double conducting = copper_volume(builder.buildAllNamed(
        enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
        mvb::DEFAULT_CORE_POLYGON_SEGMENTS, /*paintCoating=*/false));
    const double coated = copper_volume(builder.buildAllNamed(
        enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
        mvb::DEFAULT_CORE_POLYGON_SEGMENTS, /*paintCoating=*/true));

    INFO("default " << byDefault << " m^3, conducting " << conducting << " m^3, coated " << coated << " m^3");
    CHECK_THAT(byDefault, WithinRel(conducting, 1e-9));
    CHECK(byDefault < coated);
}

TEST_CASE("Real-winding build: the default draws the conducting diameter", "[abt1261][realwinding]") {
    const auto enriched = mvb::magnetic_autocomplete_safe(load_magnetic("buck_inductor_complete.json"),
                                                          /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;

    const double byDefault = copper_volume(builder.buildRealWindingTurnsNamed(enriched));
    const double conducting = copper_volume(builder.buildRealWindingTurnsNamed(
        enriched, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS, mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
        /*paintCoating=*/false));
    const double coated = copper_volume(builder.buildRealWindingTurnsNamed(
        enriched, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS, mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
        /*paintCoating=*/true));

    INFO("default " << byDefault << " m^3, conducting " << conducting << " m^3, coated " << coated << " m^3");
    CHECK_THAT(byDefault, WithinRel(conducting, 1e-9));
    CHECK(byDefault < coated);
}
