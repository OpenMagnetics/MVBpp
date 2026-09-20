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

// Where MKF's route of each terminal ends (on the pin axis, under the rail): the wrap's top.
std::map<std::pair<std::string, int>, double> route_end_heights(const OpenMagnetics::Magnetic& magnetic) {
    mvb::ConductorBuilder::LeadBendPolicy bendPolicy;   // see route_exits (ABT #1172)
    auto coil = magnetic.get_coil();
    std::map<std::pair<std::string, int>, double> out;
    for (const auto& r : coil.get_connection_layout().routes) {
        if (r.pinName.empty()) continue;
        const int terminal = r.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE ? 0 : 1;
        const auto& end = terminal == 0 ? r.pinWaypoints.front() : r.pinWaypoints.back();
        out[{r.winding + " parallel " + std::to_string(r.parallel), terminal}] = end[1];
    }
    return out;
}

// MKF's window exit of each pin lead (the window end of route.pinWaypoints), keyed like
// assigned_pins: "<winding> parallel <p>" + entrance(0)/exit(1).
std::map<std::pair<std::string, int>, std::array<double, 3>> route_exits(const OpenMagnetics::Magnetic& magnetic) {
    // ABT #1172: MKF plans the pin runs when asked for the layout, for the bend the consumer
    // declares. Read them under the SAME policy the builder draws with, or the expectations here
    // are sharp-planned routes and disagree with the drawing by the leg offset (6.9 um on the boost).
    mvb::ConductorBuilder::LeadBendPolicy bendPolicy;
    auto coil = magnetic.get_coil();
    std::map<std::pair<std::string, int>, std::array<double, 3>> out;
    for (const auto& r : coil.get_connection_layout().routes) {
        if (r.pinName.empty()) continue;
        const int terminal = r.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE ? 0 : 1;
        const auto& q = terminal == 0 ? r.pinWaypoints.back() : r.pinWaypoints.front();
        out[{r.winding + " parallel " + std::to_string(r.parallel), terminal}] = {q[0], q[1], q[2]};
    }
    return out;
}

// The drawn straight in-window run of a lead: a piece at constant (x, y) = (its slot x, MKF's row)
// running along Z up to MKF's exit depth. Its slot x is MKF's exit x, or that x moved OUTWARD by at
// most one wire radius where the 3-D lead is not proven clear at MKF's own lane (ABT #1237: MKF's
// lanes are one coated OD apart with no margin). Returns the drawn x of such a run.
std::optional<double> run_x_near(const mvb::ConductorBuilder::PathPolyline& path, const std::array<double, 3>& exit) {
    for (const auto& prim : path.prims) {
        if (prim.size() < 2) continue;
        const double x0 = prim.front()[0];
        if (std::abs(x0 - exit[0]) > path.wireRadius + 1e-12) continue;
        bool onLine = true;
        double zLo = 1e9, zHi = -1e9;
        for (const auto& q : prim) {
            if (std::abs(q[0] - x0) > 1e-9 || std::abs(q[1] - exit[1]) > 1e-9) { onLine = false; break; }
            zLo = std::min(zLo, q[2]);
            zHi = std::max(zHi, q[2]);
        }
        // its outer end stops at the exit depth, or short of it by the corner fillet into the pin run
        if (onLine && zHi - zLo > 1e-4 && zLo >= exit[2] - 1e-9 && zLo - exit[2] < 2e-3) return x0;
    }
    return std::nullopt;
}

void require_runs_at_mkf_slots(const OpenMagnetics::Magnetic& magnetic,
                               const std::vector<mvb::ConductorBuilder::PathPolyline>& paths,
                               const std::map<std::pair<std::string, int>, double>& expectedX) {
    const auto exits = route_exits(magnetic);
    REQUIRE(exits.size() == expectedX.size());
    for (const auto& [key, x] : expectedX) {
        INFO(key.first << " terminal " << key.second << ": expected exit x " << x * 1e3 << " mm");
        REQUIRE(exits.count(key));
        const auto& exit = exits.at(key);
        CHECK(std::abs(exit[0] - x) <= 1e-7);   // MKF's lane arithmetic, to the nm of the coated OD
        std::optional<double> drawnX;
        for (const auto& path : paths)
            if (path.name == key.first) drawnX = run_x_near(path, exit);
        REQUIRE(drawnX.has_value());
        INFO("drawn at x " << *drawnX * 1e3 << " mm, " << (*drawnX - exit[0]) * 1e3 << " mm off MKF's exit");
        // never moved INWARD, towards the terminal plane (a lead on the plane may go either side)
        if (std::abs(exit[0]) > 1e-12) CHECK(std::abs(*drawnX) >= std::abs(exit[0]) - 1e-12);
    }
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

TEST_CASE("Every terminal lead of the boost ends in a wrap around the pin MKF assigned it", "[pins][pinroute][abt1172]") {
    ConnectLeadsToPins on(true);
    auto enriched = enriched_boost();
    const auto assigned = assigned_pins(enriched);
    for (const auto& [key, pin] : assigned)
        UNSCOPED_INFO(key.first << " terminal " << key.second << " -> pin " << pin.name);
    REQUIRE(assigned.size() == 4u);   // Primary, 2 parallels, start + finish
    const auto ends = route_end_heights(enriched);
    REQUIRE(ends.size() == 4u);

    const auto paths = mvb::ConductorBuilder::buildAllPaths(enriched.get_coil(), bobbin_pd(enriched), false,
                                                            builder_options(enriched));
    REQUIRE(paths.size() == 2u);
    // Every run is drawn at MKF's exit slot (Coil::terminal_exit_slots, coated OD 0.943 mm): the
    // entrance bundle side by side at the plane, the exit bundle one lane out of its dragbacks.
    require_runs_at_mkf_slots(enriched, paths,
                              {{{"Primary parallel 0", 0}, 0.0}, {{"Primary parallel 1", 0}, 0.943e-3},
                               {{"Primary parallel 0", 1}, 0.943e-3}, {{"Primary parallel 1", 1}, 1.886e-3}});
    // dfm_rules.json R4.wrapTurns (2), drawn as half-turn pieces.
    const std::size_t wrapPieces = 4;
    for (const auto& path : paths) {
        REQUIRE(path.prims.size() > 2 * wrapPieces);
        for (int terminal = 0; terminal < 2; ++terminal) {
            const auto& pin = assigned.at({path.name, terminal});
            const double wrapTop = ends.at({path.name, terminal});
            const double wrapRadius = pin.radius + path.wireRadius;
            INFO(path.name << " terminal " << terminal << " on pin " << pin.name << ", wrap radius " << wrapRadius);
            double yMin = 1e9, yMax = -1e9;
            for (std::size_t k = 0; k < wrapPieces; ++k) {
                const auto& prim = terminal == 0 ? path.prims[k] : path.prims[path.prims.size() - 1 - k];
                for (const auto& q : prim) {
                    CHECK(std::abs(std::hypot(q[0] - pin.x, q[2] - pin.z) - wrapRadius) <= 1e-9);
                    yMin = std::min(yMin, q[1]);
                    yMax = std::max(yMax, q[1]);
                }
            }
            // Starts where MKF's route ends (at least a wire radius under the pin base) and advances
            // along the pin, one wire OD per turn, staying on the pin.
            CHECK(std::abs(yMax - wrapTop) <= 1e-9);
            CHECK(yMax <= pin.baseY - path.wireRadius + 1e-9);
            CHECK(std::abs((yMax - yMin) - 2 * 2 * path.wireRadius) <= 1e-9);
            CHECK(yMin > pin.tipY);
            // The free end (the FEM port) is the wrap's end, on the wrap circle.
            const auto& end = terminal == 0 ? path.end0 : path.end1;
            CHECK(std::abs(std::hypot(end[0] - pin.x, end[2] - pin.z) - wrapRadius) <= 1e-9);
            CHECK(std::abs(end[1] - yMin) <= 1e-9);
        }
    }
}

TEST_CASE("A pin changes only the terminal leads: the winding copper is byte-identical without one", "[pins][pinroute][abt1172]") {
    ConnectLeadsToPins on(true);
    auto enriched = enriched_boost();
    const auto noPins = without_pins(enriched);
    const auto withPins = mvb::ConductorBuilder::buildAllPaths(enriched.get_coil(), bobbin_pd(enriched), false,
                                                               builder_options(enriched));
    const auto without = mvb::ConductorBuilder::buildAllPaths(noPins.get_coil(), bobbin_pd(noPins), false,
                                                              builder_options(noPins));
    REQUIRE(withPins.size() == without.size());
    // A lead on a pin leaves the window at MKF's lane (Coil::terminal_exit_slots), not where the
    // free-air fan puts it, so the LEAD pieces differ by design. The TURN copper -- every piece that
    // is not terminal-lead copper by the lead-length rule (PathPolyline::primIsLead) -- must not
    // move because a pin exists: same pieces, same order, point for point. "Same point" is
    // floating-point residue, not a geometric allowance: the terminal stubs the fillet cuts at each
    // lead end carry at most 8.6e-13 m of summation noise (boost, measured 2026-09-19); a turn that
    // really moved -- by the 1.75 um to 182 um the pin slots shift the leads -- is 6+ orders above.
    constexpr double kSameCopper = 1e-12;   // metres
    auto turnPieces = [](const mvb::ConductorBuilder::PathPolyline& path) {
        REQUIRE(path.primIsLead.size() == path.prims.size());
        std::vector<std::vector<std::array<double, 3>>> out;
        for (std::size_t k = 0; k < path.prims.size(); ++k)
            if (!path.primIsLead[k]) out.push_back(path.prims[k]);
        return out;
    };
    for (std::size_t c = 0; c < withPins.size(); ++c) {
        REQUIRE(withPins[c].name == without[c].name);
        const auto a = turnPieces(withPins[c]);
        const auto b = turnPieces(without[c]);
        INFO(withPins[c].name << ": " << a.size() << " turn pieces with pins, " << b.size() << " without");
        REQUIRE(!b.empty());
        REQUIRE(a.size() == b.size());
        for (std::size_t k = 0; k < b.size(); ++k) {
            INFO(withPins[c].name << " turn piece " << k);
            REQUIRE(a[k].size() == b[k].size());
            double worst = 0.0;
            for (std::size_t q = 0; q < b[k].size(); ++q)
                worst = std::max(worst, std::hypot(a[k][q][0] - b[k][q][0], a[k][q][1] - b[k][q][1],
                                                   a[k][q][2] - b[k][q][2]));
            INFO("largest point deviation " << worst << " m");
            CHECK(worst <= kSameCopper);
        }
    }
}

TEST_CASE("The PQ 32/30 flyback draws each lead at MKF's slot, primary and secondary on opposite rows", "[pins][pinroute][abt1172]") {
    ConnectLeadsToPins on(true);
    const json magneticJson = json::parse(R"({
        "core": {"functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/30",
                                           "gapping": [], "numberStacks": 1}},
        "coil": {"bobbin": "Bobbin PQ 32/30", "functionalDescription": [
            {"name": "Primary", "numberTurns": 40, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 0.5 - Grade 1"},
            {"name": "Secondary", "numberTurns": 6, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 0.5 - Grade 1"}]}
    })");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);
    const auto assigned = assigned_pins(enriched);
    REQUIRE(assigned.size() == 4u);
    CHECK(assigned.at({"Primary parallel 0", 0}).z * assigned.at({"Secondary parallel 0", 0}).z < 0.0);   // opposite rows
    const auto paths = mvb::ConductorBuilder::buildAllPaths(enriched.get_coil(), bobbin_pd(enriched), false,
                                                            builder_options(enriched));
    require_runs_at_mkf_slots(enriched, paths,
                              {{{"Primary parallel 0", 0}, 0.0}, {{"Primary parallel 0", 1}, 0.534e-3},
                               {{"Secondary parallel 0", 0}, 0.0}, {{"Secondary parallel 0", 1}, 0.0}});
    for (const auto& path : paths)
        for (int terminal = 0; terminal < 2; ++terminal) {
            const auto& pin = assigned.at({path.name, terminal});
            const auto& end = terminal == 0 ? path.end0 : path.end1;
            INFO(path.name << " terminal " << terminal << " on pin " << pin.name);
            CHECK(std::abs(std::hypot(end[0] - pin.x, end[2] - pin.z) - (pin.radius + path.wireRadius)) <= 1e-9);
        }
}

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

TEST_CASE("A MAS round-tripped coil draws the same leads to the same pins", "[pins][pinroute][abt1172]") {
    ConnectLeadsToPins on(true);
    auto enriched = enriched_boost();
    json coilJson;
    to_json(coilJson, enriched.get_coil());
    const MAS::Coil masCoil = coilJson.get<MAS::Coil>();

    auto opts = builder_options(enriched);
    opts.femReady = false;   // the fast compound: the part names are what is compared
    const auto bobbinPd = bobbin_pd(enriched);

    const auto fromMkf = mvb::ConductorBuilder::buildAll(enriched.get_coil(), bobbinPd, false, opts);
    const auto fromMas = mvb::ConductorBuilder::buildAll(masCoil, bobbinPd, false, opts);
    REQUIRE(fromMkf.size() == fromMas.size());
    std::size_t wraps = 0;
    for (std::size_t c = 0; c < fromMkf.size(); ++c) {
        INFO(fromMkf[c].name);
        CHECK(fromMkf[c].partNames == fromMas[c].partNames);
        for (const auto& part : fromMas[c].partNames)
            if (part.find("lead wrap") != std::string::npos) ++wraps;
    }
    CHECK(wraps == 2u * 2u * 4u);   // 2 conductors x 2 terminals x 2 turns in half-turn pieces
}

TEST_CASE("The boost's FEM assembly: wraps and pin runs keep out of the pins, the core and the bobbin", "[pins][pinroute][abt1172]") {
    ConnectLeadsToPins on(true);
    auto enriched = enriched_boost();
    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                     mvb::DEFAULT_WIRE_POLYGON_SEGMENTS, mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                     /*paintCoating=*/false, /*emitCoatingShells=*/false, /*includeInsulation=*/false,
                                     /*coreCoatingThickness=*/0.0, /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto uncutBobbin = builder.buildBobbinNamed(enriched.get_coil(), enriched.get_core());
    REQUIRE(!uncutBobbin.shape.IsNull());

    std::vector<const mvb::NamedShape*> obstacles;
    std::size_t terminals = 0, conductors = 0, wraps = 0, runs = 0;
    for (const auto& ns : all) {
        if (ns.role == mvb::Role::Pin || ns.role == mvb::Role::Core) obstacles.push_back(&ns);
        // The FEM port caps: planar faces named "<conductor> terminal <k>" (they arrive as
        // Role::Terminal since buildAllNamed keeps the conductor builder's roles; counted by name
        // and shape type here, which is what this test is about).
        if (ns.name.find(" terminal ") != std::string::npos && ns.shape.ShapeType() == TopAbs_FACE) ++terminals;
    }
    REQUIRE(obstacles.size() > 12u);
    for (const auto& ns : all) {
        if (ns.role != mvb::Role::Turn || ns.partNames.empty()) continue;
        ++conductors;
        std::size_t k = 0;
        for (TopExp_Explorer e(ns.shape, TopAbs_SOLID); e.More(); e.Next(), ++k) {
            REQUIRE(k < ns.partNames.size());
            const auto& part = ns.partNames[k];
            const bool wrap = part.find("lead wrap") != std::string::npos;
            const bool run = part.find("lead pin run") != std::string::npos;
            if (!wrap && !run) continue;
            (wrap ? wraps : runs)++;
            std::vector<const mvb::NamedShape*> against = obstacles;
            against.push_back(&uncutBobbin);
            for (const auto* other : against) {
                BRepAlgoAPI_Common common(e.Current(), other->shape);
                REQUIRE(common.IsDone());
                UNSCOPED_INFO(part << " vs " << other->name);
                CHECK(volume_of(common.Shape()) < 1e-15);
            }
        }
    }
    CHECK(conductors == 2u);
    CHECK(terminals == 2u * conductors);   // the FEM ports: one planar cap per free end, unchanged
    CHECK(wraps == 2u * 4u * conductors);
    CHECK(runs > 0u);
}

TEST_CASE("A design whose bobbin has no pins draws no wraps and no pin runs", "[pins][pinroute][abt1172]") {
    ConnectLeadsToPins on(true);
    auto enriched = mvb::magnetic_autocomplete_safe(load_fixture_magnetic("flyback_transformer_complete.json"), true);
    mvb::MagneticBuilder builder;
    const auto conductors = builder.buildRealWindingTurnsNamed(enriched);
    REQUIRE(!conductors.empty());
    for (const auto& ns : conductors)
        for (const auto& part : ns.partNames) {
            CHECK(part.find("lead wrap") == std::string::npos);
            CHECK(part.find("pin run") == std::string::npos);
        }
}


// ABT #1172/#1237: MKF plans a pin run's corners for the bend radius the CONSUMER declares
// (Settings::coil_lead_bend_radius_factor), offsetting each leg off the obstacle faces so a bend of
// exactly that radius clears the edge. MVB++ declares its policy for as long as it holds the coil --
// the routes are planned when the BUILDER asks for them, not during enrichment; declaring it only
// around the enrichment left MKF planning sharp corners and put drawn copper inside the pin rail.
//
// Drawn at the COATED footprint, which is the radius MKF plans on: there the drawn bend equals the
// planned one, so any missing leg offset is copper in the rail. (At the conducting default the drawn
// bend is smaller than planned and clears by luck, which is why the gate stayed green while the
// geometry was wrong -- and why this measures the VOLUME rather than trusting the gate's 1e-15 m^3
// threshold, which a 7 um sliver slips under.)
TEST_CASE("The boost's pin runs leave room for the bend MVB++ draws: no copper in the pin rail", "[pins][pinroute][abt1172]") {
    ConnectLeadsToPins on(true);
    auto enriched = enriched_boost();
    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                     mvb::DEFAULT_WIRE_POLYGON_SEGMENTS, mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                     /*paintCoating=*/true, /*emitCoatingShells=*/false, /*includeInsulation=*/false,
                                     /*coreCoatingThickness=*/0.0, /*useRealWindingGeometry=*/true, /*femReady=*/false);

    const auto& bobbinVariant = enriched.get_coil().get_bobbin();
    const auto* bobbin = std::get_if<OpenMagnetics::Bobbin>(&bobbinVariant);
    REQUIRE(bobbin != nullptr);
    const auto rails = mvb::BobbinBuilder::buildPinRailsNamed(bobbin->get_pin_rails(),
                                                              bobbin->get_name().value_or("Bobbin"));
    REQUIRE(!rails.empty());

    std::size_t conductors = 0;
    for (const auto& ns : all) {
        // A conductor PAINTED at its outer diameter is still the conductor and stays Role::Turn;
        // only the separate shell pass (emitCoatingShells, " coating" suffix) is a TurnCoating
        // (ABT #1245). This case paints, so the copper to measure is the Turns.
        if (ns.role != mvb::Role::Turn) continue;
        ++conductors;
        for (const auto& rail : rails) {
            BRepAlgoAPI_Common common(rail.shape, ns.shape);
            REQUIRE(common.IsDone());
            GProp_GProps props;
            BRepGProp::VolumeProperties(common.Shape(), props);
            INFO(rail.name << " vs " << ns.name << ": " << props.Mass() * 1e9 << " mm^3");
            CHECK(props.Mass() == 0.0);   // exactly none: the legs are offset for the drawn bend
        }
    }
    CHECK(conductors == 2u);
}
