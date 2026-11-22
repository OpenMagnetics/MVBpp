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
#include "mvb/SectionBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"

#include <BRepAlgoAPI_Common.hxx>
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

// ABT #2511. MKF flattens a litz bundle that would not fit otherwise: a stadium of the round bundle's area, its short
// side k x the outer diameter (here k = 0.99: Litz 270x0.16 served, 3.7555 mm round, wound 3.7179 x 3.7772 mm on a
// T 50/30/15 whose faces the first layer touches). The copper footprint of such a turn was the UNDEFORMED bare bundle,
// 3.7228 mm: wider than the wire it is the inside of, and 2.4 um into the ferrite at every first-layer turn -- the 2D
// mesher refused the section ("0.000293 mm^2 ... claimed by core_0, turn_Boost inductor_35"). The copper of a deformed
// bundle is its own outline taken in by the serving.
TEST_CASE("A flattened litz turn's copper lies inside the flattened wire and clear of the core it touches", "[abt2511]") {
    const auto enriched = mvb::magnetic_autocomplete_safe(load_magnetic("abt2511_t50_flattened_litz.json"));
    OpenMagnetics::Coil coil = enriched.get_coil();   // mutable copy: resolve_wire is non-const
    const auto wire = coil.resolve_wire(0);
    const auto turns = coil.get_turns_description();   // held: the getter returns a copy
    REQUIRE(turns.has_value());
    size_t flattened = 0;
    for (const auto& turn : *turns) {
        // The premise: MKF wound this bundle flattened.
        const auto shape = turn.get_cross_sectional_shape();
        REQUIRE(shape.has_value());
        REQUIRE(*shape == MAS::TurnCrossSectionalShape::OVAL);
        const auto dimensions = turn.get_dimensions();
        REQUIRE(dimensions.has_value());
        REQUIRE(dimensions->size() == 2);
        const auto& wound = *dimensions;
        const auto copper = mvb::TurnBuilder::wireDimensions(wire, turn, /*paintCoating=*/false);
        INFO(turn.get_name() << ": copper " << copper.first * 1e3 << " x " << copper.second * 1e3 << " mm in a wire wound "
                             << wound[0] * 1e3 << " x " << wound[1] * 1e3 << " mm");
        CHECK(copper.first < wound[0]);
        CHECK(copper.second < wound[1]);
        ++flattened;
    }
    REQUIRE(flattened == 49);

    // What the 2D mesher sections: true circles for the core, the copper footprint for the turns.
    mvb::MagneticBuilder builder;
    const auto solids = builder.buildAllNamed(enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                              mvb::DEFAULT_WIRE_POLYGON_SEGMENTS, /*coreSegments=*/0, /*paintCoating=*/false);
    const auto faces = mvb::SectionBuilder::cut2DFaces(solids, mvb::parseSectionPlane("XY"), 0.0);
    const mvb::NamedShape* core = nullptr;
    size_t turnFaces = 0, inTheCore = 0;
    for (const auto& ns : faces)
        if (ns.role == mvb::Role::Core) core = &ns;
    REQUIRE(core != nullptr);
    for (const auto& ns : faces) {
        if (ns.role != mvb::Role::Turn) continue;
        ++turnFaces;
        BRepAlgoAPI_Common common(ns.shape, core->shape);
        REQUIRE(common.IsDone());
        GProp_GProps props;
        BRepGProp::SurfaceProperties(common.Shape(), props);
        if (props.Mass() > 0.0) {
            ++inTheCore;
            UNSCOPED_INFO(ns.name << " shares " << props.Mass() * 1e6 << " mm2 with the core");
        }
    }
    REQUIRE(turnFaces == 98);   // every turn crosses the mid-plane twice: in the bore and outside the ring
    CHECK(inTheCore == 0);
}
