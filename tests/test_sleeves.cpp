// WP5 lead sleeving (ABT #1174): MVB++ draws the sleeve MKF decided and reserved; it decides
// nothing itself. Plan: /home/alf/OpenMagnetics/docs/2026-09-13_wp5_mvbpp_plan.md.
#include <catch2/catch_test_macros.hpp>

#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "support/Utils.h"

#include <filesystem>
#include <fstream>

#ifndef MAS_EXAMPLES_DIR
#define MAS_EXAMPLES_DIR "."
#endif

using json = nlohmann::json;

namespace {

// The margin-wound flyback of the MAS examples (PQ 32/30, reinforced, IEC 60664-1 + IEC 62368-1):
// MKF's mas_autocomplete runs the insulation coordinator and records each lead's sleeve on the
// winding's connections[]; the magnetic it returns carries them.
OpenMagnetics::Magnetic sleeved_flyback() {
    std::ifstream file(std::filesystem::path(MAS_EXAMPLES_DIR) / "24_margin_interleaved_flyback_pq3230_3c94.json");
    REQUIRE(file.good());
    const json masJson = json::parse(file);
    // Processed inputs (the excitations' waveforms from the design requirements), as MKF's own
    // lead-sleeve tests load this file; the coordinator reads the insulation requirement from them.
    OpenMagnetics::Mas mas;
    mas.set_inputs(OpenMagnetics::Inputs(masJson.at("inputs"), true));
    mas.set_magnetic(OpenMagnetics::Magnetic(masJson.at("magnetic")));
    const auto completed = OpenMagnetics::mas_autocomplete(mas, false);
    json magneticJson;
    OpenMagnetics::to_json(magneticJson, completed.get_magnetic());
    return mvb::magnetic_autocomplete_safe(magneticJson, true);
}

}  // namespace

TEST_CASE("The margin flyback carries MKF's sleeved terminal routes through the MAS round trip", "[sleeve][abt1174]") {
    const auto magnetic = sleeved_flyback();
    auto coil = magnetic.get_coil();
    size_t sleeved = 0, terminals = 0;
    for (const auto& r : coil.get_connection_layout().routes) {
        const bool terminal = r.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE ||
                              r.kind == OpenMagnetics::ConnectionKind::TERMINAL_EXIT;
        if (!terminal) continue;
        ++terminals;
        if (r.sleeveOuterDiameter) {
            ++sleeved;
            WARN(r.winding << " parallel " << r.parallel << " "
                           << (r.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE ? "entrance" : "exit")
                           << ": sleeve OD " << *r.sleeveOuterDiameter * 1e3 << " mm");
        }
    }
    WARN(sleeved << " of " << terminals << " terminal routes sleeved");
    CHECK(sleeved > 0);
}
