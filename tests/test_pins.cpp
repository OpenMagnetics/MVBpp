// ABT #1171 (WP2 "Bobbin pins"): a catalogue bobbin with a pin footprint draws its pins where
// MKF placed them.
//
// MKF's Bobbin::expand_pinout turns a record's pinout (count, rows, pitches, pinDescription,
// orientation) into processedDescription.pins[]; MVB++'s PinBuilder only turns each MAS::Pin into
// a solid. The fixture is boost_inductor_complete.json, whose "Bobbin PQ 26/25" is one of the
// catalogue records with a complete footprint (12 round THT pins, 6+6, row distance 25.4 mm since MAS 901af03, ABT #1209;
// the rails holding them are drawn since ABT #1249, see test_pinrail.cpp).
// flyback_transformer_complete.json's "Bobbin ETD 34" states no footprint, so MKF places no pins
// there and none may be drawn -- MVB++ must not invent any.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "mvb/MagneticBuilder.h"
#include "mvb/PinBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <nlohmann/json.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <numbers>
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

// The pins MKF placed, keyed by pin name, from the same autocomplete buildAllNamed runs.
std::map<std::string, MAS::Pin> mkf_pins(const json& magneticJson, std::string& bobbinName) {
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson);
    const auto& bobbinVariant = enriched.get_coil().get_bobbin();
    const auto* bobbin = std::get_if<OpenMagnetics::Bobbin>(&bobbinVariant);
    REQUIRE(bobbin != nullptr);
    REQUIRE(bobbin->get_name().has_value());
    bobbinName = bobbin->get_name().value();
    REQUIRE(bobbin->get_processed_description().has_value());
    std::map<std::string, MAS::Pin> out;
    const auto pins = bobbin->get_processed_description()->get_pins();
    if (!pins) return out;
    for (const auto& pin : *pins) {
        REQUIRE(pin.get_name().has_value());
        out.emplace(pin.get_name().value(), pin);
    }
    return out;
}

MAS::Pin round_pin(std::vector<double> dimensions, std::optional<std::vector<double>> coordinates) {
    MAS::Pin pin;
    pin.set_name("7");
    pin.set_shape(MAS::PinShape::ROUND);
    pin.set_type(MAS::PinDescriptionType::THT);
    pin.set_dimensions(std::move(dimensions));
    pin.set_coordinates(std::move(coordinates));
    return pin;
}

}  // namespace

TEST_CASE("A catalogue bobbin with a footprint draws every pin where MKF placed it", "[pins][abt1171]") {
    const json magneticJson = load_fixture_magnetic("boost_inductor_complete.json");
    std::string bobbinName;
    const auto expected = mkf_pins(magneticJson, bobbinName);
    INFO("bobbin '" << bobbinName << "', MKF pins: " << expected.size());
    REQUIRE(expected.size() == 12u);   // Bobbin PQ 26/25: numberPins 12

    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(magneticJson.get<MAS::Magnetic>(), /*includeBobbin=*/true);

    std::vector<mvb::NamedShape> pins;
    for (const auto& ns : all)
        if (ns.role == mvb::Role::Pin) pins.push_back(ns);
    INFO("solids: " << all.size() << ", pin solids: " << pins.size());
    REQUIRE(pins.size() == expected.size());

    const std::string prefix = bobbinName + " pin ";
    for (const auto& pin : pins) {
        REQUIRE(pin.name.rfind(prefix, 0) == 0);
        const std::string pinName = pin.name.substr(prefix.size());
        auto it = expected.find(pinName);
        REQUIRE(it != expected.end());
        const auto& mkf = it->second;
        const auto coordinates = mkf.get_coordinates().value();

        // Centred on MKF's coordinate (bobbin.json: coordinates are the pin's CENTRE).
        GProp_GProps props;
        BRepGProp::VolumeProperties(pin.shape, props);
        const gp_Pnt centre = props.CentreOfMass();
        UNSCOPED_INFO(pin.name << " centre (" << centre.X() << ", " << centre.Y() << ", " << centre.Z()
                      << ") vs MKF (" << coordinates[0] << ", " << coordinates[1] << ", " << coordinates[2] << ")");
        CHECK(std::abs(centre.X() - coordinates[0]) <= 1e-6);
        CHECK(std::abs(centre.Y() - coordinates[1]) <= 1e-6);
        CHECK(std::abs(centre.Z() - coordinates[2]) <= 1e-6);

        // A round pin is a cylinder of the stated diameter and length along Y.
        REQUIRE(mkf.get_shape() == MAS::PinShape::ROUND);
        const auto& d = mkf.get_dimensions();
        CHECK(std::abs(volume_of(pin.shape) - std::numbers::pi * d[0] * d[0] / 4 * d[2]) <= 1e-12);
        Bnd_Box box;
        BRepBndLib::AddOptimal(pin.shape, box, /*useTriangulation=*/false, /*useShapeTolerance=*/false);
        double xlo, ylo, zlo, xhi, yhi, zhi;
        box.Get(xlo, ylo, zlo, xhi, yhi, zhi);
        CHECK(std::abs((xhi - xlo) - d[0]) <= 1e-9);
        CHECK(std::abs((zhi - zlo) - d[0]) <= 1e-9);
        CHECK(std::abs((yhi - ylo) - d[2]) <= 1e-9);
    }

    // Gate: no pin shares volume with a conductor, a core piece, or the (cut) bobbin.
    std::size_t turns = 0, cores = 0, bobbins = 0;
    for (const auto& ns : all) {
        if (ns.role == mvb::Role::Turn) ++turns;
        if (ns.role == mvb::Role::Core) ++cores;
        if (ns.role == mvb::Role::Bobbin) ++bobbins;
    }
    REQUIRE(turns > 0);
    REQUIRE(cores > 0);
    REQUIRE(bobbins == 1u);
    for (const auto& pin : pins) {
        for (const auto& ns : all) {
            if (ns.role != mvb::Role::Turn && ns.role != mvb::Role::TurnCoating &&
                ns.role != mvb::Role::Core && ns.role != mvb::Role::Bobbin) continue;
            BRepAlgoAPI_Common common(pin.shape, ns.shape);
            REQUIRE(common.IsDone());
            UNSCOPED_INFO(pin.name << " vs " << ns.name);
            CHECK(volume_of(common.Shape()) < 1e-15);
        }
    }
}

TEST_CASE("A bobbin whose record states no footprint draws no pins", "[pins][abt1171]") {
    const json magneticJson = load_fixture_magnetic("flyback_transformer_complete.json");
    std::string bobbinName;
    REQUIRE(mkf_pins(magneticJson, bobbinName).empty());

    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(magneticJson.get<MAS::Magnetic>(), /*includeBobbin=*/true);
    std::size_t bobbins = 0;
    for (const auto& ns : all) {
        CHECK(ns.role != mvb::Role::Pin);
        if (ns.role == mvb::Role::Bobbin) ++bobbins;
    }
    REQUIRE(bobbins == 1u);   // the bobbin itself IS drawn, so "no pins" is not "no bobbin"
}

TEST_CASE("PinBuilder refuses a pin it cannot draw instead of dropping it", "[pins][abt1171]") {
    using Catch::Matchers::ContainsSubstring;
    CHECK_THROWS_WITH(mvb::PinBuilder::buildPin(round_pin({0.0009, 0.0009}, std::vector<double>{0, 0, 0})),
                      ContainsSubstring("'7'") && ContainsSubstring("2 dimensions"));
    CHECK_THROWS_WITH(mvb::PinBuilder::buildPin(round_pin({0.0009, 0.0009, 0.007}, std::nullopt)),
                      ContainsSubstring("'7'") && ContainsSubstring("no coordinates"));

    MAS::CoreBobbinProcessedDescription pd;
    pd.set_pins(std::vector<MAS::Pin>{round_pin({0.0009, 0.0009, 0.007}, std::vector<double>{0.001, -0.01, 0.0}),
                                      round_pin({0.0009, 0.0009, 0.007}, std::nullopt)});
    CHECK_THROWS_WITH(mvb::PinBuilder::buildPinsNamed(pd, "Bobbin X"), ContainsSubstring("no coordinates"));

    MAS::CoreBobbinProcessedDescription none;
    CHECK(mvb::PinBuilder::buildPinsNamed(none, "Bobbin X").empty());
}
