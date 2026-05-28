#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "Magnetic.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <cmath>
#include <numbers>

using json = nlohmann::json;

int main(int argc, char* argv[]) {
    if (argc < 2) { std::cerr << "usage: " << argv[0] << " input.json [filter]\n"; return 1; }
    std::ifstream f(argv[1]);
    if (!f.is_open()) { std::cerr << "cannot open " << argv[1] << "\n"; return 1; }
    json j; f >> j;
    mvb::patch_dimension_nominals(j);
    json magneticJson = j.contains("magnetic") ? j.at("magnetic") : j;
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson);

    std::string filter = (argc > 2) ? argv[2] : "";
    auto turnsOpt = enriched.get_coil().get_turns_description();
    if (!turnsOpt || turnsOpt->empty()) {
        std::cerr << "No turns_description\n"; return 1;
    }
    for (size_t i = 0; i < turnsOpt->size(); ++i) {
        const auto& t = (*turnsOpt)[i];
        std::string name = t.get_name();
        if (!filter.empty() && name.find(filter) == std::string::npos) continue;
        const auto& coords = t.get_coordinates();
        double cx = coords.size() > 0 ? coords[0] : 0.0;
        double cy = coords.size() > 1 ? coords[1] : 0.0;
        double r = std::sqrt(cx*cx + cy*cy);
        double angle = std::atan2(cy, cx) * 180.0 / std::numbers::pi;
        double rotation = t.get_rotation().value_or(0.0);
        std::cout << "[" << i << "] " << name
                  << "  coords=(" << cx << ", " << cy << ")"
                  << "  r=" << r << "  angle=" << angle << "°"
                  << "  rotation=" << rotation << "°";
        auto addCoords = t.get_additional_coordinates();
        if (addCoords && !addCoords->empty() && (*addCoords)[0].size() >= 2) {
            double ax = (*addCoords)[0][0], ay = (*addCoords)[0][1];
            double or_ = std::sqrt(ax*ax + ay*ay);
            double oa = std::atan2(ay, ax) * 180.0 / std::numbers::pi;
            std::cout << "  outer=(" << ax << ", " << ay << ") r=" << or_ << " angle=" << oa << "°";
        }
        std::cout << "\n";
    }
    return 0;
}
