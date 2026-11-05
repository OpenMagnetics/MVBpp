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
#include <optional>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "mvb/MagneticBuilder.h"
#include "mvb/PinBuilder.h"
#include "mvb/BobbinBuilder.h"
#include "mvb/ConductorBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <BRepAlgoAPI_Common.hxx>
#include <TopExp_Explorer.hxx>
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

// ---- ABT #1172 (WP3): every terminal lead ends in a wrap around the pin MKF assigned it ----------

namespace {

struct AssignedPin {
    std::string name;
    double x = 0, z = 0, radius = 0, baseY = 0, tipY = 0;
};

// connections[] of the enriched coil: "<winding> parallel <p>" + entrance(0)/exit(1) -> the pin.
std::map<std::pair<std::string, int>, AssignedPin> assigned_pins(const OpenMagnetics::Magnetic& enriched) {
    const auto& coil = enriched.get_coil();
    const auto* bobbin = std::get_if<OpenMagnetics::Bobbin>(&coil.get_bobbin());
    REQUIRE(bobbin != nullptr);
    std::map<std::string, MAS::Pin> pins;
    const auto pinList = bobbin->get_processed_description()->get_pins();
    REQUIRE(pinList.has_value());
    for (const auto& pin : *pinList) pins[pin.get_name().value()] = pin;
    std::map<std::pair<std::string, int>, AssignedPin> out;
    const auto& windings = coil.get_functional_description();
    for (const auto& w : windings) {
        const auto connections = w.get_connections();
        if (!connections) continue;
        for (const auto& c : *connections) {
            if (!c.get_pin_name() || !c.get_end() || c.get_end().value() == MAS::End::TAP) continue;
            const int terminal = c.get_end().value() == MAS::End::START ? 0 : 1;
            const auto& pin = pins.at(c.get_pin_name().value());
            AssignedPin a;
            a.name = c.get_pin_name().value();
            const auto xyz = pin.get_coordinates().value();
            a.x = xyz[0];
            a.z = xyz[2];
            a.radius = pin.get_dimensions()[0] / 2;
            a.baseY = xyz[1] + pin.get_dimensions()[2] / 2;
            a.tipY = xyz[1] - pin.get_dimensions()[2] / 2;
            for (int64_t p = 0; p < w.get_number_parallels(); ++p) {
                if (c.get_parallel() && c.get_parallel().value() != p) continue;
                out[{w.get_name() + " parallel " + std::to_string(p), terminal}] = a;
            }
        }
    }
    return out;
}

// Settings coil_connect_leads_to_pins (MKF, off by default): set for the scope of a test, and
// always put back, so no other test inherits it.
struct ConnectLeadsToPins {
    explicit ConnectLeadsToPins(bool on) { OpenMagnetics::settings.set_coil_connect_leads_to_pins(on); }
    ~ConnectLeadsToPins() { OpenMagnetics::settings.set_coil_connect_leads_to_pins(false); }
};

OpenMagnetics::Magnetic enriched_boost() {
    return mvb::magnetic_autocomplete_safe(load_fixture_magnetic("boost_inductor_complete.json"), true);
}

MAS::CoreBobbinProcessedDescription bobbin_pd(const OpenMagnetics::Magnetic& magnetic) {
    const auto* bobbin = std::get_if<OpenMagnetics::Bobbin>(&magnetic.get_coil().get_bobbin());
    REQUIRE(bobbin != nullptr);
    REQUIRE(bobbin->get_processed_description().has_value());
    return bobbin->get_processed_description().value();
}

// The conductor builder's options as MagneticBuilder sets them for a single-window core. The
// collision gate and the pin gate stay ON.
mvb::ConductorBuilder::Options builder_options(const OpenMagnetics::Magnetic& magnetic) {
    mvb::MagneticBuilder builder;
    mvb::ConductorBuilder::Options opts;
    opts.femReady = true;
    opts.paintCoating = false;
    for (const auto& ns : builder.buildCoreNamed(magnetic.get_core())) opts.coreObstacles.push_back(ns.shape);
    return opts;
}

OpenMagnetics::Magnetic without_pins(OpenMagnetics::Magnetic magnetic) {
    auto bobbin = std::get<OpenMagnetics::Bobbin>(magnetic.get_coil().get_bobbin());
    auto pd = bobbin.get_processed_description().value();
    pd.set_pins(std::nullopt);
    bobbin.set_processed_description(pd);
    auto coil = magnetic.get_coil();
    auto fd = coil.get_functional_description();
    for (auto& w : fd) w.set_connections(std::nullopt);
    coil.set_functional_description(fd);
    coil.set_bobbin(bobbin);
    magnetic.set_coil(coil);
    return magnetic;
}

}  // namespace

TEST_CASE("coil_connect_leads_to_pins off: the terminals are the free-air leads, pins or no pins", "[pins][pinroute][abt1172]") {
    // Pins assigned (the setting on for autocomplete), then drawn with the setting OFF (its default):
    // every primitive, terminal runs included, equals the drawing of the same design with no pins
    // at all, i.e. exactly what MVB++ drew before WP3.
    OpenMagnetics::Magnetic enriched;
    {
        ConnectLeadsToPins on(true);
        enriched = enriched_boost();
    }
    REQUIRE(assigned_pins(enriched).size() == 4u);
    REQUIRE(!OpenMagnetics::settings.get_coil_connect_leads_to_pins());
    const auto noPins = without_pins(enriched);
    const auto off = mvb::ConductorBuilder::buildAllPaths(enriched.get_coil(), bobbin_pd(enriched), false,
                                                          builder_options(enriched));
    const auto bare = mvb::ConductorBuilder::buildAllPaths(noPins.get_coil(), bobbin_pd(noPins), false,
                                                           builder_options(noPins));
    REQUIRE(off.size() == bare.size());
    for (std::size_t c = 0; c < off.size(); ++c) {
        INFO(off[c].name);
        CHECK(off[c].prims == bare[c].prims);
        CHECK(off[c].end0 == bare[c].end0);
        CHECK(off[c].end1 == bare[c].end1);
        // Both terminals on the common tip plane, not on a pin.
        CHECK(std::abs(off[c].end0[2] - off[c].end1[2]) <= 1e-12);
    }
    // And the default autocomplete assigns nothing (MKF's side of the flag).
    CHECK(assigned_pins(enriched_boost()).empty());
}
