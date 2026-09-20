// ABT #1169 (WP0): NamedShape::role — every solid buildAllNamed produces says what it IS,
// not only what it is called.
//
// Names are the channel that survives a STEP round-trip, so OMFEM still classifies imported
// geometry by name. The role is the in-process channel for consumers that hold the
// NamedShape vector. If a producer forgets to set it, the shape silently claims to be a core
// piece (the enum's default) — exactly the failure mode WP0 exists to remove — so this test
// asserts the role of EVERY solid of two complete MAS fixtures against the role its name
// implies.
//
//   ./mvb_tests "[roles]"

#include <catch2/catch_test_macros.hpp>

#include "constructive_models/Magnetic.h"
#include "mvb/MagneticBuilder.h"
#include "mvb/NamedShape.h"
#include "mvb/Utils.h"
#include "MAS.hpp"

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

bool ends_with(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// The role a solid's NAME implies, derived independently of the builder so the test is a
// real cross-check rather than a restatement of the production code. Mirrors D2's name
// families plus the ones MVB++ already emits.
mvb::Role role_from_name(const std::string& name, const std::string& bobbinName) {
    if (name == "FR4Board")                              return mvb::Role::FR4;
    if (ends_with(name, " solder joint"))                return mvb::Role::Solder;
    if (name.find(" terminal ") != std::string::npos)    return mvb::Role::Terminal;
    if (name.rfind("insulation_layer", 0) == 0)          return mvb::Role::Insulation;
    if (name.rfind("Spacer_", 0) == 0)                   return mvb::Role::Spacer;
    if (name.rfind("Shunt_", 0) == 0)                    return mvb::Role::Shunt;
    if (ends_with(name, " sleeve"))                      return mvb::Role::Sleeve;
    if (name.find(" divider ") != std::string::npos)     return mvb::Role::Divider;
    if (name.find(" pin ") != std::string::npos)         return mvb::Role::Pin;
    if (!bobbinName.empty() && name.rfind(bobbinName, 0) == 0) return mvb::Role::Bobbin;
    if (name.rfind("Bobbin", 0) == 0)                    return mvb::Role::Bobbin;
    // A conductor: "<winding> parallel <p>[ turn <i>]", possibly " coating".
    if (name.find(" parallel ") != std::string::npos || name.rfind("Turn_", 0) == 0)
        return ends_with(name, " coating") ? mvb::Role::TurnCoating : mvb::Role::Turn;
    return ends_with(name, " coating") ? mvb::Role::CoreCoating : mvb::Role::Core;
}

std::string bobbin_name_of(const OpenMagnetics::Magnetic& m) {
    const auto& bob = m.get_coil().get_bobbin();
    if (const auto* s = std::get_if<std::string>(&bob)) return *s;
    if (const auto* b = std::get_if<OpenMagnetics::Bobbin>(&bob)) {
        const auto nm = b->get_name();
        if (nm && !nm->empty()) return *nm;
    }
    return "Bobbin";
}

// expectedBobbins: a CONCENTRIC part gets one former solid; a toroid gets none, because MKF's
// BobbinTDataProcessor is a zero-thickness virtual bobbin and buildBobbinNamed returns a null
// shape for it (the gap WP4/D1 fills with a real toroid base). Stated per fixture rather than
// assumed, so the count stays a real assertion instead of a guess.
// realWinding: the femReady real-winding build (ABT #1245), the only path that emits terminal caps
// (and, for a foil, solder joints) beside the copper. Drawn without a bobbin, as OMFEM builds it.
void check_roles(const std::string& fixture, int expectedBobbins, bool accessories = false,
                 bool realWinding = false, bool expectSolder = false) {
    const std::filesystem::path path = std::filesystem::path(MAS_COMPLETE_DIR) / fixture;
    REQUIRE(std::filesystem::exists(path));
    std::ifstream in(path);
    json mas; in >> mas;
    const json magneticJson = mas.contains("magnetic") ? mas.at("magnetic") : mas;

    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/realWinding);
    const std::string bobbinName = bobbin_name_of(enriched);

    mvb::MagneticBuilder builder;
    // `accessories` turns on the three producers the default call never reaches: the conformal
    // core-coating shell, the per-turn coating shells and the insulation layers. Without it the
    // test would only ever see Core/Turn/Bobbin and could not catch a producer that forgot to
    // set the role on a coating or an insulation layer.
    auto all = realWinding
        ? builder.buildAllNamed(enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                /*useRealWindingGeometry=*/true, /*femReady=*/true)
        : accessories
        ? builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                /*paintCoating=*/false, /*emitCoatingShells=*/true,
                                /*includeInsulation=*/true,
                                /*coreCoatingThickness=*/50e-6)
        : builder.buildAllNamed(enriched, /*includeBobbin=*/true);
    REQUIRE(all.size() > 1);

    int cores = 0, turns = 0, bobbins = 0, coreCoatings = 0, turnCoatings = 0;
    int terminals = 0, solders = 0;
    for (const auto& ns : all) {
        const mvb::Role expected = role_from_name(ns.name, bobbinName);
        INFO(fixture << ": solid '" << ns.name << "' role=" << mvb::role_name(ns.role)
                     << " expected=" << mvb::role_name(expected));
        CHECK(ns.role == expected);
        if (ns.role == mvb::Role::Core)   ++cores;
        if (ns.role == mvb::Role::Turn)   ++turns;
        if (ns.role == mvb::Role::Bobbin) ++bobbins;
        if (ns.role == mvb::Role::CoreCoating) ++coreCoatings;
        if (ns.role == mvb::Role::TurnCoating) ++turnCoatings;
        if (ns.role == mvb::Role::Terminal) ++terminals;
        if (ns.role == mvb::Role::Solder)   ++solders;
    }
    if (realWinding) {
        // Two caps per (winding, parallel) conductor. Demanding them makes the name cross-check
        // above bite: with none present it could not see a cap mislabelled as a Turn.
        CHECK(terminals > 0);
        CHECK(terminals == 2 * turns);
        if (expectSolder) CHECK(solders > 0);
    }
    if (accessories) {
        CHECK(coreCoatings == cores);    // one conformal shell per core piece
        CHECK(turnCoatings == turns);    // one outer footprint per conductor
    }
    // A fixture whose roles are ALL the enum default would pass a pure name-agreement check
    // if the name map were also degenerate; demand the three roles we know must be present.
    CHECK(cores   > 0);
    CHECK(turns   > 0);
    CHECK(bobbins == expectedBobbins);
}

}  // namespace

TEST_CASE("buck_inductor_complete: every solid carries its role", "[roles][abt1169]") {
    check_roles("buck_inductor_complete.json", /*expectedBobbins=*/0);   // toroid: virtual bobbin
}

TEST_CASE("flyback_transformer_complete: every solid carries its role", "[roles][abt1169]") {
    check_roles("flyback_transformer_complete.json", /*expectedBobbins=*/1);
}

TEST_CASE("buck_inductor_complete: coating shells and insulation layers carry their roles",
          "[roles][abt1169]") {
    check_roles("buck_inductor_complete.json", /*expectedBobbins=*/0, /*accessories=*/true);
}

TEST_CASE("buck_inductor_complete real winding: terminal caps are Terminal, not Turn",
          "[roles][abt1245]") {
    check_roles("buck_inductor_complete.json", /*expectedBobbins=*/0, /*accessories=*/false,
                /*realWinding=*/true);
}

TEST_CASE("two_switch_forward real winding: foil solder joints are Solder, caps are Terminal",
          "[roles][abt1245]") {
    check_roles("two_switch_forward_transformer_complete.json", /*expectedBobbins=*/0,
                /*accessories=*/false, /*realWinding=*/true, /*expectSolder=*/true);
}

TEST_CASE("role_name covers every Role", "[roles][abt1169]") {
    const mvb::Role every[] = {
        mvb::Role::Core, mvb::Role::CoreCoating, mvb::Role::Turn, mvb::Role::TurnCoating,
        mvb::Role::Insulation, mvb::Role::Bobbin, mvb::Role::Pin, mvb::Role::Divider,
        mvb::Role::Base, mvb::Role::Spacer, mvb::Role::Shunt, mvb::Role::Sleeve,
        mvb::Role::Terminal, mvb::Role::Solder, mvb::Role::FR4};
    for (mvb::Role r : every) {
        INFO("role " << static_cast<int>(r));
        CHECK(std::string(mvb::role_name(r)) != "unknown");
    }
}

// The real-winding assembly: the conductor builder roles every solid it emits (the conductor, its
// FEM terminal caps, a foil's solder bodies, a lead's sleeve), and buildAllNamed must hand those
// roles through. It used to re-derive them from the name, so every terminal cap reached the
// assembly as a Turn (found on ABT #1172) -- and a sleeve would have been meshed as copper.
TEST_CASE("flyback_transformer_complete real winding: terminal caps keep Role::Terminal", "[roles][abt1169][abt1174]") {
    const std::filesystem::path path = std::filesystem::path(MAS_COMPLETE_DIR) / "flyback_transformer_complete.json";
    REQUIRE(std::filesystem::exists(path));
    std::ifstream in(path);
    json mas; in >> mas;
    const json magneticJson = mas.contains("magnetic") ? mas.at("magnetic") : mas;
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    const std::string bobbinName = bobbin_name_of(enriched);

    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                     mvb::DEFAULT_WIRE_POLYGON_SEGMENTS, mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                     /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                     /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                     /*useRealWindingGeometry=*/true, /*femReady=*/true);
    int terminals = 0, turns = 0;
    for (const auto& ns : all) {
        const mvb::Role expected = role_from_name(ns.name, bobbinName);
        INFO("solid '" << ns.name << "' role=" << mvb::role_name(ns.role) << " expected=" << mvb::role_name(expected));
        CHECK(ns.role == expected);
        if (ns.role == mvb::Role::Terminal) ++terminals;
        if (ns.role == mvb::Role::Turn) ++turns;
    }
    CHECK(turns > 0);
    CHECK(terminals == 2 * turns);   // two FEM port caps per conductor
}
