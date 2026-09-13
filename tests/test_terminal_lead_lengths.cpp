// ABT #1215 ([leadlength]): the terminal-lead copper MVB++ draws, per winding, measured on the
// finished real-winding centreline. OMFEM's 3D R_dc includes these leads; MKF's reference must be
// validated against the same number, so every piece is checked independently of primLength and
// every end is proven to run contiguously from the free terminal tip into the conductor.
#include <catch2/catch_test_macros.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"
#include "json.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <string>

using json = nlohmann::json;

namespace {

json loadLeadFixture(const std::string& name) {
    std::ifstream f("testData/" + name);
    if (!f.good()) f = std::ifstream("tests/realwinding_fixtures/" + name);
    if (!f.good()) f = std::ifstream("tests/mas_complete_fixtures/" + name);
    if (!f.good()) {
        const auto root = std::filesystem::path{__FILE__}.parent_path();
        f = std::ifstream(root / "mas_complete_fixtures" / name);
        if (!f.good()) f = std::ifstream(root / "realwinding_fixtures" / name);
    }
    REQUIRE(f.good());
    json j = json::parse(f);
    return j.contains("magnetic") ? j.at("magnetic") : j;
}

double dist(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                     (a[2] - b[2]) * (a[2] - b[2]));
}

struct Measured {
    std::map<std::string, mvb::ConductorBuilder::TerminalLeadLength> leads;
    std::vector<mvb::ConductorBuilder::PathPolyline> paths;
};

// The FEM product OMFEM meshes: bare copper, femReady -- the settings buildRealWindingPaths uses,
// so the path free ends can be compared with the lead tips.
Measured measure(const std::string& fixture) {
    auto enriched = mvb::magnetic_autocomplete_safe(loadLeadFixture(fixture),
                                                   /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    Measured m;
    m.leads = builder.measureTerminalLeadLengths(enriched, /*paintCoating=*/false,
                                                 /*femReady=*/true);
    m.paths = builder.buildRealWindingPaths(enriched);
    return m;
}

// Independent checks on every winding: each piece's length recomputed from its printed geometry
// (not through primLength), ends summed by hand, every end contiguous from the terminal tip, and
// the turn copper the lead hands over to NOT among the pieces.
void requireConsistent(const Measured& m, const std::string& fixture) {
    REQUIRE(!m.leads.empty());
    for (const auto& [winding, t] : m.leads) {
        INFO(fixture << " winding '" << winding << "'");
        REQUIRE(t.parallels >= 1);
        REQUIRE(t.per_end.size() == 2 * t.parallels);
        double handTotal = 0.0;
        for (size_t k = 0; k < t.per_end.size(); ++k) {
            const auto& e = t.per_end[k];
            INFO("parallel " << e.parallel << " " << e.end);
            CHECK(e.parallel == k / 2);
            CHECK(e.end == (k % 2 == 0 ? "entrance" : "exit"));
            REQUIRE(!e.pieces.empty());
            double handEnd = 0.0;
            for (size_t i = 0; i < e.pieces.size(); ++i) {
                const auto& pc = e.pieces[i];
                INFO("piece " << i << " '" << pc.label << "' " << pc.kind);
                double independent = 0.0;
                if (pc.kind == "SEG") {
                    independent = dist(pc.start, pc.end);
                }
                else if (pc.kind == "ARC3") {
                    // From the chord, not from radius*sweep: 2 R asin(c / 2R), valid below pi.
                    REQUIRE(pc.sweep_rad < std::numbers::pi);
                    REQUIRE(pc.radius_m > 0.0);
                    independent = 2.0 * pc.radius_m *
                                  std::asin(std::min(1.0, dist(pc.start, pc.end) / (2.0 * pc.radius_m)));
                }
                else {
                    FAIL("lead piece kind " << pc.kind << " has no independent check here");
                }
                CHECK(std::fabs(independent - pc.length_m) < 1e-9);
                // A lead piece is never a wrap / turn transition. The one exception by name is a
                // terminal fillet arc: it is cut from the wrap's terminal stub and keeps the stub's
                // label with " fillet arc N" appended.
                if (pc.label.find("(terminal stub) fillet arc") == std::string::npos) {
                    CHECK(pc.label.find("->") == std::string::npos);
                    CHECK(pc.label.rfind("wrap", 0) != 0);
                }
                if (i + 1 < e.pieces.size())
                    CHECK(dist(pc.end, e.pieces[i + 1].start) < 1e-9);   // contiguous
                handEnd += independent;
            }
            CHECK(std::fabs(handEnd - e.length_m) < 1e-9);
            handTotal += handEnd;

            // The end starts (entrance) / finishes (exit) at the conductor's free terminal tip.
            const std::string name = winding + " parallel " + std::to_string(e.parallel);
            const mvb::ConductorBuilder::PathPolyline* path = nullptr;
            for (const auto& p : m.paths)
                if (p.name == name) path = &p;
            REQUIRE(path != nullptr);
            if (e.end == "entrance")
                CHECK(dist(e.pieces.front().start, path->end0) < 1e-9);
            else
                CHECK(dist(e.pieces.back().end, path->end1) < 1e-9);

            // THE TURN IS EXCLUDED. The primitive the lead hands over to is turn copper: find it
            // in the full sampled centreline (it starts where the entrance lead ends / ends where
            // the exit lead starts) and prove it is not one of the pieces.
            const auto& junction = e.end == "entrance" ? e.pieces.back().end : e.pieces.front().start;
            bool found = false;
            for (size_t q = 0; q < path->prims.size(); ++q) {
                const auto& prim = path->prims[q];
                const auto& near = e.end == "entrance" ? prim.front() : prim.back();
                const auto& far = e.end == "entrance" ? prim.back() : prim.front();
                if (dist(near, junction) > 1e-9 || dist(far, junction) < 1e-9) continue;
                bool isPiece = false;
                for (const auto& pc : e.pieces)
                    if (dist(pc.start, prim.front()) < 1e-9 && dist(pc.end, prim.back()) < 1e-9)
                        isPiece = true;
                if (isPiece) continue;
                found = true;
                double turnLen = 0.0;
                for (size_t s = 1; s < prim.size(); ++s)
                    turnLen += dist(prim[s - 1], prim[s]);
                INFO("adjacent turn primitive " << q << ", sampled length " << turnLen * 1e3 << " mm");
                CHECK(turnLen > 0.0);
                break;
            }
            CHECK(found);
        }
        CHECK(std::fabs(handTotal - t.total_m) < 1e-9);
    }
}

void report(const Measured& m, const std::string& fixture) {
    std::cout << std::setprecision(12);
    for (const auto& [winding, t] : m.leads) {
        std::cout << "[leadlength] " << fixture << " '" << winding << "': total " << t.total_m * 1e3
                  << " mm, parallels " << t.parallels << "\n";
        for (const auto& e : t.per_end) {
            std::cout << "[leadlength]   parallel " << e.parallel << " " << e.end << ": "
                      << e.length_m * 1e3 << " mm\n";
            for (const auto& pc : e.pieces)
                std::cout << "[leadlength]     " << pc.kind << " " << pc.length_m * 1e3 << " mm  '"
                          << pc.label << "'  (" << pc.start[0] * 1e3 << "," << pc.start[1] * 1e3
                          << "," << pc.start[2] * 1e3 << ") -> (" << pc.end[0] * 1e3 << ","
                          << pc.end[1] * 1e3 << "," << pc.end[2] * 1e3 << ") mm"
                          << (pc.kind == "ARC3" ? "  R " + std::to_string(pc.radius_m * 1e3) +
                                                      " mm sweep " + std::to_string(pc.sweep_rad) +
                                                      " rad"
                                                : std::string())
                          << "\n";
        }
    }
}

} // namespace

TEST_CASE("Terminal leads: the 3-parallel buck toroid reports its six drawn leads (ABT #1215)",
          "[leadlength][realwinding]") {
    const auto m = measure("buck_inductor_complete.json");
    report(m, "buck_inductor_complete");
    REQUIRE(m.leads.size() == 1);
    const auto& t = m.leads.begin()->second;
    CHECK(t.parallels == 3);
    CHECK(t.per_end.size() == 6);
    requireConsistent(m, "buck_inductor_complete");
    // Every toroid terminal finishes in its -Y drop: the entrance's first piece and the exit's
    // last piece are straight legs running along Y, and all six tips share one plane.
    double planeY = std::numeric_limits<double>::quiet_NaN();
    for (const auto& e : t.per_end) {
        const auto& drop = e.end == "entrance" ? e.pieces.front() : e.pieces.back();
        const auto& tip = e.end == "entrance" ? drop.start : drop.end;
        const auto& root = e.end == "entrance" ? drop.end : drop.start;
        CHECK(drop.kind == "SEG");
        CHECK(std::fabs(tip[0] - root[0]) < 1e-12);
        CHECK(std::fabs(tip[2] - root[2]) < 1e-12);
        CHECK(tip[1] < root[1]);
        if (std::isnan(planeY)) planeY = tip[1];
        CHECK(std::fabs(tip[1] - planeY) < 1e-9);
        // The chain the router builds in the default VERTICAL mounting (ABT #1248; pieces are in the
        // exported frame): drop, corner, axial. Before #1248 (the ring flat, a radial run over the
        // face) it was drop, corner, radial, corner, axial, and this asserted the same >= 3.
        CHECK(e.pieces.size() >= 3);
    }
}

TEST_CASE("Terminal leads: the common-mode choke reports both windings (ABT #1215)",
          "[leadlength][realwinding]") {
    const auto m = measure("common_mode_choke_complete.json");
    report(m, "common_mode_choke_complete");
    CHECK(m.leads.size() == 2);
    requireConsistent(m, "common_mode_choke_complete");
}

TEST_CASE("Terminal leads: a concentric flyback's exits are measured (ABT #1215)",
          "[leadlength][realwinding]") {
    const auto m = measure("flyback_transformer_complete.json");
    report(m, "flyback_transformer_complete");
    CHECK(m.leads.size() >= 2);
    requireConsistent(m, "flyback_transformer_complete");
    // Every concentric terminal here bends off its wrap through TerminalFillet's arcs, and those
    // arcs are lead copper: each end must carry at least one, next to a straight lead leg.
    for (const auto& [winding, t] : m.leads) {
        for (const auto& e : t.per_end) {
            INFO(winding << " parallel " << e.parallel << " " << e.end);
            size_t fillets = 0, legs = 0;
            for (const auto& pc : e.pieces) {
                if (pc.label.find("(terminal stub) fillet arc") != std::string::npos) ++fillets;
                else if (pc.kind == "SEG") ++legs;
            }
            CHECK(fillets >= 1);
            CHECK(legs >= 1);
        }
    }
}

TEST_CASE("Terminal leads: the sidecar next to the STEP round-trips the measurement (ABT #1215)",
          "[leadlength][realwinding]") {
    const auto m = measure("buck_inductor_complete.json");
    const auto dir = std::filesystem::path{__FILE__}.parent_path().parent_path() / "output";
    std::filesystem::create_directories(dir);
    const std::string step = (dir / "leadlength_buck.step").string();
    CHECK(mvb::MagneticBuilder::terminalLeadSidecarPath(step) ==
          (dir / "leadlength_buck.leads.json").string());
    std::filesystem::remove(mvb::MagneticBuilder::terminalLeadSidecarPath(step));
    const std::string written = mvb::MagneticBuilder::writeTerminalLeadSidecar(m.leads, step);
    REQUIRE(std::filesystem::exists(written));
    std::ifstream f(written);
    const json j = json::parse(f);
    REQUIRE(j.size() == m.leads.size());
    for (const auto& [winding, t] : m.leads) {
        const std::string key = "winding_" + winding;
        REQUIRE(j.contains(key));
        CHECK(j.at(key).at("terminal_lead_length_m").get<double>() == t.total_m);   // bit-exact
        CHECK(j.at(key).at("parallels").get<size_t>() == t.parallels);
        REQUIRE(j.at(key).at("ends").size() == t.per_end.size());
        double sum = 0.0;
        for (size_t k = 0; k < t.per_end.size(); ++k) {
            const auto& je = j.at(key).at("ends").at(k);
            CHECK(je.at("parallel").get<size_t>() == t.per_end[k].parallel);
            CHECK(je.at("end").get<std::string>() == t.per_end[k].end);
            CHECK(je.at("length_m").get<double>() == t.per_end[k].length_m);
            CHECK(je.at("pieces").size() == t.per_end[k].pieces.size());
            sum += je.at("length_m").get<double>();
        }
        CHECK(std::fabs(sum - t.total_m) < 1e-12);
    }
}
