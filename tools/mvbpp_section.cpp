// mvbpp_section: builds the core from a MAS JSON, then writes an SVG
// for the XY-plane and XZ-plane sections of the fused core.
#include "mvb/SectionBuilder.h"
#include "mvb/SectionDrawing.h"
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Core.h"

#include <BRep_Builder.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <TopoDS_Compound.hxx>

#include <fstream>
#include <iostream>
#include <filesystem>
#include <nlohmann/json.hpp>

namespace fs = std::filesystem;
using json = nlohmann::json;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <input.json> <output_dir>\n"
                  << "Writes <output_dir>/<basename>_section_xy.svg and _xz.svg\n";
        return 1;
    }
    fs::path inPath  = argv[1];
    fs::path outDir  = argv[2];
    fs::create_directories(outDir);

    std::ifstream f(inPath);
    if (!f.is_open()) { std::cerr << "Cannot open " << inPath << "\n"; return 1; }
    json j; f >> j;
    mvb::patch_dimension_nominals(j);

    json magneticJson = j.contains("magnetic") ? j.at("magnetic") : j;

    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson);
    mvb::MagneticBuilder builder;
    auto namedPieces = builder.buildCoreNamed(enriched.get_core());
    if (namedPieces.empty()) {
        std::cerr << "No core pieces built\n";
        return 1;
    }

    // Fuse core halves so the section captures the assembled cross-section.
    TopoDS_Shape core = namedPieces[0].shape;
    for (size_t i = 1; i < namedPieces.size(); ++i) {
        BRepAlgoAPI_Fuse fuser(core, namedPieces[i].shape);
        if (fuser.IsDone()) core = fuser.Shape();
    }

    std::string base = inPath.stem().string();
    auto outXY = (outDir / (base + "_section_xy.svg")).string();
    auto outXZ = (outDir / (base + "_section_xz.svg")).string();

    try {
        mvb::SectionBuilder::writeSectionSvg(core, mvb::SectionPlane::XY, outXY);
        std::cout << "Wrote " << outXY << "\n";
    } catch (const std::exception& e) {
        std::cerr << "XY section failed: " << e.what() << "\n";
    }
    try {
        mvb::SectionBuilder::writeSectionSvg(core, mvb::SectionPlane::XZ, outXZ);
        std::cout << "Wrote " << outXZ << "\n";
    } catch (const std::exception& e) {
        std::cerr << "XZ section failed: " << e.what() << "\n";
    }

    // Dimensioned FrontView + TopView, both reading from the enriched
    // Magnetic (needs gap coordinates and processedDescription populated
    // by magnetic_autocomplete).
    auto outFront = (outDir / (base + "_dimensioned_front.svg")).string();
    try {
        mvb::SectionDrawing::writeDimensionedFrontView(enriched, outFront);
        std::cout << "Wrote " << outFront << "\n";
    } catch (const std::exception& e) {
        std::cerr << "FrontView drawing failed: " << e.what() << "\n";
    }
    auto outTop = (outDir / (base + "_dimensioned_top.svg")).string();
    try {
        mvb::SectionDrawing::writeDimensionedTopView(enriched, outTop);
        std::cout << "Wrote " << outTop << "\n";
    } catch (const std::exception& e) {
        std::cerr << "TopView drawing failed: " << e.what() << "\n";
    }
    return 0;
}
