// ABT #1249: the bobbin's plastic pin rails (the bars that hold the pins) are drawn, fused with
// the bottom flange, exactly where MKF's Bobbin::get_pin_rails() puts them.
//
// Fixture boost_inductor_complete.json: "Bobbin PQ 26/25" (Miles-Platts PQ0040), rails from the
// record's labels a, b, a1, b1, H1, H3 -- four blocks, two per pin row, from the bottom flange's
// outer face (y = -H1/2 = -7.735 mm) down to the pin standoff (y = -(H1/2 + H3) = -14.165 mm),
// below the PQ 26/25 core (half height 12.375 mm) and outside its depth (|z| >= 9.905 mm > 9.5).
//
// Pins are separate solids that START on the rail underside: they touch the rail, they do not
// pass through it, so the pin cutters of the bobbin cut remove nothing from a rail and no volume
// is shared in the FEM assembly. The rails are part of the bottom-flange solid of the bobbin
// compound, under the bobbin's own name (Role::Bobbin), which OMFEM classifies as bobbin.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "mvb/BobbinBuilder.h"
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <nlohmann/json.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

json load_fixture_magnetic(const std::string& file) {
    auto path = std::filesystem::path(MAS_COMPLETE_DIR) / file;
    REQUIRE(std::filesystem::exists(path));
    std::ifstream in(path);
    json j;
    in >> j;
    return j.at("magnetic");
}

double volume_of(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}

double common_volume(const TopoDS_Shape& a, const TopoDS_Shape& b) {
    Bnd_Box ba, bb;
    BRepBndLib::Add(a, ba);
    BRepBndLib::Add(b, bb);
    if (ba.IsOut(bb)) return 0.0;
    BRepAlgoAPI_Common common(a, b);
    REQUIRE(common.IsDone());
    return volume_of(common.Shape());
}

OpenMagnetics::Bobbin enriched_bobbin(const json& magneticJson) {
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson);
    const auto* bobbin = std::get_if<OpenMagnetics::Bobbin>(&enriched.get_coil().get_bobbin());
    REQUIRE(bobbin != nullptr);
    return *bobbin;
}

// Every drawn pin starts on the underside of exactly one rail block and lies inside its footprint.
void check_pins_on_rails(const std::vector<mvb::NamedShape>& all, const std::vector<OpenMagnetics::Bobbin::PinRailBlock>& rails) {
    std::size_t pins = 0;
    for (const auto& ns : all) {
        if (ns.role != mvb::Role::Pin) continue;
        ++pins;
        Bnd_Box box;
        BRepBndLib::AddOptimal(ns.shape, box, false, false);
        double xlo, ylo, zlo, xhi, yhi, zhi;
        box.Get(xlo, ylo, zlo, xhi, yhi, zhi);
        std::size_t holders = 0;
        for (const auto& rail : rails) {
            if (xlo >= rail.centre[0] - rail.halfExtents[0] - 1e-9 && xhi <= rail.centre[0] + rail.halfExtents[0] + 1e-9 &&
                zlo >= rail.centre[2] - rail.halfExtents[2] - 1e-9 && zhi <= rail.centre[2] + rail.halfExtents[2] + 1e-9 &&
                std::abs(yhi - (rail.centre[1] - rail.halfExtents[1])) < 1e-7)
                ++holders;
        }
        UNSCOPED_INFO(ns.name << " y [" << ylo << ", " << yhi << "]");
        CHECK(holders == 1u);
    }
    CHECK(pins > 0u);
}

}  // namespace

TEST_CASE("boost_inductor_complete draws the PQ 26/25 pin rails under the core, holding every pin", "[pinrail][pins][abt1249]") {
    const json magneticJson = load_fixture_magnetic("boost_inductor_complete.json");
    const auto bobbin = enriched_bobbin(magneticJson);
    const auto rails = bobbin.get_pin_rails();
    REQUIRE(rails.size() == 4u);
    const auto core = OpenMagnetics::Core(mvb::magnetic_autocomplete_safe(magneticJson).get_core());
    const double coreHalfHeight = core.get_height() / 2;

    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(magneticJson.get<MAS::Magnetic>(), /*includeBobbin=*/true);
    const mvb::NamedShape* drawnBobbin = nullptr;
    std::size_t bobbins = 0;
    for (const auto& ns : all)
        if (ns.role == mvb::Role::Bobbin) { drawnBobbin = &ns; ++bobbins; }
    REQUIRE(bobbins == 1u);
    CHECK(drawnBobbin->name == bobbin.get_name().value());
    std::size_t solids = 0;
    for (TopExp_Explorer exp(drawnBobbin->shape, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    CHECK(solids == 3u);   // tube, top flange, bottom flange + rails

    // The drawn bobbin reaches down to the pin standoff, below the core.
    Bnd_Box box;
    BRepBndLib::AddOptimal(drawnBobbin->shape, box, false, false);
    double xlo, ylo, zlo, xhi, yhi, zhi;
    box.Get(xlo, ylo, zlo, xhi, yhi, zhi);
    CHECK(std::abs(ylo - (-(0.01547 / 2 + 0.00643))) < 1e-7);
    CHECK(ylo < -coreHalfHeight);
    CHECK(std::abs(zhi - 0.014605) < 1e-7);   // b/2: the rails are the widest part across the rows

    // Every rail block is inside the drawn (cut) bobbin whole: nothing of it was cut away.
    const auto railShapes = mvb::BobbinBuilder::buildPinRailsNamed(rails, drawnBobbin->name);
    for (const auto& rail : railShapes) {
        INFO(rail.name);
        CHECK(std::abs(common_volume(rail.shape, drawnBobbin->shape) - volume_of(rail.shape)) < 1e-13);
    }
    // No overlap between a rail and a core piece, a conductor, or a pin.
    for (const auto& rail : railShapes) {
        for (const auto& ns : all) {
            if (ns.role == mvb::Role::Bobbin) continue;
            UNSCOPED_INFO(rail.name << " vs " << ns.name);
            CHECK(common_volume(rail.shape, ns.shape) < 1e-15);
        }
    }
    check_pins_on_rails(all, rails);
}

TEST_CASE("Pin rails add exactly their own volume, fused into the bottom flange; no pins, no change", "[pinrail][abt1249]") {
    const auto bobbin = enriched_bobbin(load_fixture_magnetic("boost_inductor_complete.json"));
    const auto pd = bobbin.get_processed_description().value();
    const auto rails = bobbin.get_pin_rails();
    const auto bare = mvb::BobbinBuilder::buildBobbin(pd, pd.get_wall_thickness(), true, mvb::DEFAULT_CORE_POLYGON_SEGMENTS);
    const auto railed = mvb::BobbinBuilder::buildBobbin(pd, pd.get_wall_thickness(), true, mvb::DEFAULT_CORE_POLYGON_SEGMENTS, rails);
    double railVolume = 0.0;
    for (const auto& rail : rails) railVolume += 8 * rail.halfExtents[0] * rail.halfExtents[1] * rail.halfExtents[2];
    CHECK(std::abs(volume_of(railed) - volume_of(bare) - railVolume) < 1e-12);

    // A bobbin without pins is drawn exactly as before (flyback_transformer_complete: Bobbin ETD 34).
    const auto noPins = enriched_bobbin(load_fixture_magnetic("flyback_transformer_complete.json"));
    REQUIRE(noPins.get_pin_rails().empty());
    const auto noPinsPd = noPins.get_processed_description().value();
    const auto before = mvb::BobbinBuilder::buildBobbin(noPinsPd, noPinsPd.get_wall_thickness(), true, mvb::DEFAULT_CORE_POLYGON_SEGMENTS);
    const auto after = mvb::BobbinBuilder::buildBobbin(noPinsPd, noPinsPd.get_wall_thickness(), true, mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                                        noPins.get_pin_rails());
    CHECK(volume_of(after) == volume_of(before));

    // A rail that does not hang from the flange face is refused, never drawn floating.
    auto floating = rails;
    floating[0].centre[1] -= 0.001;
    CHECK_THROWS_WITH(mvb::BobbinBuilder::buildBobbin(pd, pd.get_wall_thickness(), true, mvb::DEFAULT_CORE_POLYGON_SEGMENTS, floating),
                      Catch::Matchers::ContainsSubstring("not on the drawn bottom flange face"));
}

TEST_CASE("A quick E 42/21/15 bobbin with synthesised pins draws its table rails", "[pinrail][abt1249]") {
    OpenMagnetics::settings.reset();
    OpenMagnetics::settings.set_coil_quick_bobbin_generate_pins(true);
    json coilJson;
    coilJson["bobbin"] = "Basic";
    coilJson["functionalDescription"] = json::array();
    json winding;
    winding["name"] = "Primary";
    winding["numberTurns"] = 20;
    winding["numberParallels"] = 1;
    winding["isolationSide"] = "primary";
    winding["wire"] = "Round 0.5 - Grade 1";
    coilJson["functionalDescription"].push_back(winding);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagnetics::Core(json::parse(R"({"functionalDescription": {"type": "two-piece set", "material": "N87",
        "shape": "E 42/21/15", "gapping": [], "numberStacks": 1}})")));
    magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
    magnetic = OpenMagnetics::magnetic_autocomplete(magnetic);
    OpenMagnetics::settings.reset();
    const auto* bobbin = std::get_if<OpenMagnetics::Bobbin>(&magnetic.get_coil().get_bobbin());
    REQUIRE(bobbin != nullptr);
    REQUIRE(bobbin->get_processed_description()->get_pins());
    const auto rails = bobbin->get_pin_rails();
    REQUIRE(rails.size() == 2u);

    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(magnetic, /*includeBobbin=*/true);
    const mvb::NamedShape* drawnBobbin = nullptr;
    for (const auto& ns : all)
        if (ns.role == mvb::Role::Bobbin) drawnBobbin = &ns;
    REQUIRE(drawnBobbin != nullptr);
    for (const auto& rail : mvb::BobbinBuilder::buildPinRailsNamed(rails, drawnBobbin->name)) {
        INFO(rail.name);
        CHECK(std::abs(common_volume(rail.shape, drawnBobbin->shape) - volume_of(rail.shape)) < 1e-13);
        for (const auto& ns : all) {
            if (ns.role == mvb::Role::Bobbin) continue;
            UNSCOPED_INFO(rail.name << " vs " << ns.name);
            CHECK(common_volume(rail.shape, ns.shape) < 1e-15);
        }
    }
    check_pins_on_rails(all, rails);
}
