// mvbpp_pathprobe -- ABT #871 verification. Builds the real-winding conductor CENTRELINES and
// checks the contract that matters for multi-column placement: every MKF turn keeps its exact
// drawn position, i.e. each turn's crossing (x, y) is a point the finished conductor actually
// passes through. Prints each path's extent and the worst per-turn miss.
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "Magnetic.h"
#include "Utils.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

using json = nlohmann::json;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: mvbpp_pathprobe <input.json>\n";
        return 2;
    }
    std::ifstream f(argv[1]);
    if (!f.is_open()) {
        std::cerr << "cannot open " << argv[1] << "\n";
        return 2;
    }
    json j;
    f >> j;
    mvb::patch_dimension_nominals(j);
    json magneticJson = j.contains("magnetic") ? j.at("magnetic") : j;

    OpenMagnetics::Magnetic enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);
    mvb::MagneticBuilder builder;
    const auto paths = builder.buildRealWindingPaths(enriched);

    auto coil = enriched.get_mutable_coil();
    const auto turns = coil.get_turns_description().value();

    std::cout << std::fixed << std::setprecision(4);
    for (const auto& path : paths) {
        double lowX = 1e9, highX = -1e9, lowY = 1e9, highY = -1e9, lowZ = 1e9, highZ = -1e9;
        std::vector<std::array<double, 3>> points;
        for (const auto& primitive : path.prims) {
            for (const auto& p : primitive) {
                points.push_back(p);
                lowX = std::min(lowX, p[0]); highX = std::max(highX, p[0]);
                lowY = std::min(lowY, p[1]); highY = std::max(highY, p[1]);
                lowZ = std::min(lowZ, p[2]); highZ = std::max(highZ, p[2]);
            }
        }
        std::cout << path.name << ": end0 (" << path.end0[0]*1e3 << "," << path.end0[1]*1e3 << "," << path.end0[2]*1e3
                  << ") dir0 (" << path.dir0[0] << "," << path.dir0[1] << "," << path.dir0[2] << ")  end1 ("
                  << path.end1[0]*1e3 << "," << path.end1[1]*1e3 << "," << path.end1[2]*1e3 << ") dir1 ("
                  << path.dir1[0] << "," << path.dir1[1] << "," << path.dir1[2] << ")  [mm; lead ends = terminal caps]\n";
        std::cout << path.name << ": x [" << lowX * 1e3 << "," << highX * 1e3 << "] y ["
                  << lowY * 1e3 << "," << highY * 1e3 << "] z [" << lowZ * 1e3 << ","
                  << highZ * 1e3 << "] mm, " << points.size() << " sampled points\n";

        // The contract MVB++ owes MKF: every drawn turn keeps its exact position. For a
        // conductor wrapping a column whose axis is `axisX`, that means the racetrack's own
        // x-extent at the turn's height reaches exactly |turnX - axisX| from the axis. The axis
        // is the midpoint of the wrap's x-extent (a closed loop around the leg is symmetric
        // about it), which for a main-column conductor is 0 and for a leg is the leg's own x.
        const double axisX = 0.5 * (lowX + highX);
        double worst = 0.0;
        std::string worstTurn;
        for (const auto& turn : turns) {
            const std::string owner = turn.get_winding() + " parallel " +
                                      std::to_string(turn.get_parallel());
            if (owner != path.name) continue;
            const double turnY = turn.get_coordinates()[1];
            const double drawnRadial = std::abs(turn.get_coordinates()[0] - axisX);
            // Closest approach of the sampled centreline to the drawn crossing, measured in the
            // (radial, axial) half-plane the crossing is drawn in.
            double closest = std::numeric_limits<double>::max();
            for (const auto& p : points) {
                closest = std::min(closest, std::hypot(std::abs(p[0] - axisX) - drawnRadial,
                                                       p[1] - turnY));
            }
            if (closest > worst) { worst = closest; worstTurn = turn.get_name(); }
        }
        std::cout << "    axis x = " << axisX * 1e3 << " mm; worst turn-radial miss "
                  << worst * 1e6 << " um"
                  << (worstTurn.empty() ? "" : " (" + worstTurn + ")") << "\n";
    }

    // ABT #1215: the terminal-lead pieces of the same centreline (bare copper, femReady -- the
    // settings buildRealWindingPaths uses), one line each, so the per-winding total can be summed
    // by hand from the printed geometry.
    const auto leads = builder.measureTerminalLeadLengths(enriched, /*paintCoating=*/false,
                                                          /*femReady=*/true);
    std::cout << std::setprecision(17);   // round-trips: a hand sum to 1e-9 m needs every digit
    for (const auto& [winding, t] : leads) {
        std::cout << "LEADS '" << winding << "': total " << t.total_m << " m, parallels "
                  << t.parallels << "\n";
        for (const auto& e : t.per_end) {
            std::cout << "  parallel " << e.parallel << " " << e.end << ": " << e.length_m << " m\n";
            for (const auto& pc : e.pieces) {
                std::cout << "    " << pc.kind << " len " << pc.length_m << " m  a=(" << pc.start[0]
                          << "," << pc.start[1] << "," << pc.start[2] << ") b=(" << pc.end[0] << ","
                          << pc.end[1] << "," << pc.end[2] << ")";
                if (pc.kind != "SEG")
                    std::cout << " R=" << pc.radius_m << " sweep=" << pc.sweep_rad;
                std::cout << "  '" << pc.label << "'\n";
            }
        }
    }
    return 0;
}
