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
void check_roles(const std::string& fixture, int expectedBobbins, bool accessories = false) {
    const std::filesystem::path path = std::filesystem::path(MAS_COMPLETE_DIR) / fixture;
    REQUIRE(std::filesystem::exists(path));
    std::ifstream in(path);
    json mas; in >> mas;
    const json magneticJson = mas.contains("magnetic") ? mas.at("magnetic") : mas;

    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson);
    const std::string bobbinName = bobbin_name_of(enriched);

    mvb::MagneticBuilder builder;
    // `accessories` turns on the three producers the default call never reaches: the conformal
    // core-coating shell, the per-turn coating shells and the insulation layers. Without it the
    // test would only ever see Core/Turn/Bobbin and could not catch a producer that forgot to
    // set the role on a coating or an insulation layer.
    auto all = accessories
        ? builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                /*paintCoating=*/false, /*emitCoatingShells=*/true,
                                /*includeInsulation=*/true,
                                /*coreCoatingThickness=*/50e-6)
        : builder.buildAllNamed(enriched, /*includeBobbin=*/true);
    REQUIRE(all.size() > 1);

    int cores = 0, turns = 0, bobbins = 0, coreCoatings = 0, turnCoatings = 0;
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
