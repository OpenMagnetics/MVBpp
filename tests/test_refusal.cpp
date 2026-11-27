// A refusal carries its reason (mvb::RefusalReason), so nobody has to read its message to know what it is.
#include <catch2/catch_test_macros.hpp>
#include "constructive_models/Magnetic.h"
#include "mvb/Refusal.h"
#include "mvb/Utils.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>

using json = nlohmann::json;

TEST_CASE("Every refusal reason has one fixed name, and the name gives the reason back", "[refusal]") {
    size_t seen = 0;
    for (const auto& [reason, name] : mvb::detail::kRefusalReasonNames) {
        INFO(name);
        CHECK(mvb::to_string(reason) == name);
        const auto back = mvb::refusalReasonFromString(name);
        REQUIRE(back.has_value());
        CHECK(*back == reason);
        ++seen;
    }
    CHECK(seen == 10);
    CHECK_FALSE(mvb::refusalReasonFromString("").has_value());
    CHECK_FALSE(mvb::refusalReasonFromString("planar real winding").has_value());
}

TEST_CASE("A planar construction asked for real winding is refused as PlanarRealWindingUnsupported", "[refusal]") {
    std::ifstream f(std::string(MAS_EXAMPLES_DIR) + "/09_planar_xfmr_er2510_3c94.json");
    REQUIRE(f.good());
    json j = json::parse(f);
    mvb::patch_dimension_nominals(j);
    const json magnetic = j.contains("magnetic") ? j.at("magnetic") : j;
    bool refused = false;
    try {
        mvb::magnetic_autocomplete_safe(magnetic, /*useRealWindingGeometry=*/true);
    } catch (const mvb::Refusal& e) {
        refused = true;
        CHECK(e.reason() == mvb::RefusalReason::PlanarRealWindingUnsupported);
        CHECK_FALSE(std::string(e.what()).empty());
    }
    CHECK(refused);
}
