// The enamel gate's verdict is handed back to the caller of buildAllNamed, so a consumer can
// rely on the proof (OMFEM skips welding a winding's parallels only when the gate CERTIFIED that
// no two conductors share volume). A wrong Certified would silently skip a weld that is needed,
// so every non-proof outcome is pinned here:
//   - a certified build reports Certified;
//   - the same build with the checks disabled reports Skipped;
//   - the same out-variable reused for a certified then a skipped build reads Skipped (no stale
//     proof survives into the next build);
//   - a build that throws, or never reaches the gate, leaves NotRun.
//
//   ./mvb_tests "[gateverdict]"

#include <catch2/catch_test_macros.hpp>

#include "constructive_models/Magnetic.h"
#include "mvb/ConductorBuilder.h"
#include "mvb/MagneticBuilder.h"
#include "mvb/NamedShape.h"
#include "mvb/Utils.h"

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
using mvb::EnamelGateVerdict;

json load_magnetic(const std::string& fixture) {
    const std::filesystem::path path = std::filesystem::path(MAS_COMPLETE_DIR) / fixture;
    REQUIRE(std::filesystem::exists(path));
    std::ifstream in(path);
    json mas; in >> mas;
    return mas.contains("magnetic") ? mas.at("magnetic") : mas;
}

std::vector<mvb::NamedShape> build(const mvb::MagneticBuilder& builder,
                                   const OpenMagnetics::Magnetic& magnetic,
                                   bool realWinding, bool skipChecks, EnamelGateVerdict* out) {
    return builder.buildAllNamed(magnetic, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                 mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                 mvb::DEFAULT_CORE_POLYGON_SEGMENTS, /*paintCoating=*/false,
                                 /*emitCoatingShells=*/false, /*includeInsulation=*/false,
                                 /*coreCoatingThickness=*/0.0, realWinding, /*femReady=*/false,
                                 skipChecks, out);
}

const OpenMagnetics::Magnetic& buck() {
    static const OpenMagnetics::Magnetic m = mvb::magnetic_autocomplete_safe(
        load_magnetic("buck_inductor_complete.json"), /*useRealWindingGeometry=*/true);
    return m;
}

}  // namespace

TEST_CASE("Enamel gate verdict: a certified real-winding build reports Certified",
          "[gateverdict][realwinding]") {
    mvb::MagneticBuilder builder;
    EnamelGateVerdict v = EnamelGateVerdict::NotRun;
    const auto all = build(builder, buck(), /*realWinding=*/true, /*skipChecks=*/false, &v);
    REQUIRE_FALSE(all.empty());
    CHECK(v == EnamelGateVerdict::Certified);
}

TEST_CASE("Enamel gate verdict: the same build with checks disabled reports Skipped",
          "[gateverdict][realwinding]") {
    mvb::MagneticBuilder builder;
    EnamelGateVerdict v = EnamelGateVerdict::NotRun;
    build(builder, buck(), /*realWinding=*/true, /*skipChecks=*/true, &v);
    CHECK(v == EnamelGateVerdict::Skipped);
}

TEST_CASE("Enamel gate verdict: no stale Certified survives into the next build",
          "[gateverdict][realwinding]") {
    mvb::MagneticBuilder builder;
    EnamelGateVerdict v = EnamelGateVerdict::NotRun;
    build(builder, buck(), /*realWinding=*/true, /*skipChecks=*/false, &v);
    REQUIRE(v == EnamelGateVerdict::Certified);
    build(builder, buck(), /*realWinding=*/true, /*skipChecks=*/true, &v);
    CHECK(v == EnamelGateVerdict::Skipped);
}

TEST_CASE("Enamel gate verdict: a per-turn build runs no gate and reports NotRun",
          "[gateverdict]") {
    mvb::MagneticBuilder builder;
    const auto perTurn = mvb::magnetic_autocomplete_safe(load_magnetic("buck_inductor_complete.json"));
    EnamelGateVerdict v = EnamelGateVerdict::Certified;   // as a previous build left it
    build(builder, perTurn, /*realWinding=*/false, /*skipChecks=*/false, &v);
    CHECK(v == EnamelGateVerdict::NotRun);
}

TEST_CASE("Enamel gate verdict: a build that throws leaves NotRun, not a previous Certified",
          "[gateverdict]") {
    mvb::MagneticBuilder builder;
    EnamelGateVerdict v = EnamelGateVerdict::Certified;   // as a previous build left it
    MAS::Magnetic noCoil;                                 // no core, no coil: refused on entry
    CHECK_THROWS(builder.buildAllNamed(noCoil, false, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, false, false, false,
                                       0.0, /*useRealWindingGeometry=*/true, false, false, &v));
    CHECK(v == EnamelGateVerdict::NotRun);
}
