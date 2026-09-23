#include "mvb/MagneticBuilder.h"
#include "mvb/StepExporter.h"
#include <map>
#include <optional>
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "Magnetic.h"
#include "Mas.h"
#include "Utils.h"
#include "support/Painter.h"
#include "support/Settings.h"
#include <BRepBuilderAPI_Transform.hxx>
#include <OSD_ThreadPool.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numbers>
#include <filesystem>

namespace fs = std::filesystem;
using json = nlohmann::json;

static void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog << " [options] <input.json>\n"
              << "Options:\n"
              << "  -o, --output <path>   Output STEP file path\n"
              << "  -d, --output-dir <dir> Output directory (batch mode)\n"
              << "  --no-mkf              Skip MKF enrichment\n"
              << "  --real                Real winding: continuous conductor per (winding, parallel)\n"
              << "  --skip-collision-check  DIAGNOSTIC: build the conductors even when two of\n"
                 "                        them overlap. Interpenetrating copper: a picture, never\n"
                 "                        a part. Also via MVB_SKIP_COLLISION_CHECK.\n"
              << "  --fem                 FEM geometry: one-piece / conformal conductors (slow); "
                 "default is the fast drawing compound.\n"
                 "                        With --real, also writes <output>.leads.json: the\n"
                 "                        terminal-lead copper length per winding (ABT #1215)\n"
              << "  --coated              Draw conductors at the OUTER (insulation) envelope.\n"
                 "                        DEFAULT IS THE CONDUCTING (copper) footprint (ABT #1261):\n"
                 "                        true copper cross-section, and the enamel MKF's one-OD\n"
                 "                        turn pitch reserves becomes real air between neighbours.\n"
                 "                        --copper is accepted and ignored (it is the default now).\n"
              << "  --segments <N>        Wire AND core polygon segments (0 = exact analytic\n"
              << "                        curves). The two facet in LOCKSTEP: a faceted wire\n"
              << "                        against an exact core wall touches at every polygon\n"
              << "                        vertex (measured, 03_buck).\n"
              << "  --core-segments <N>   Core polygon segments explicitly (overrides the above)\n"
              << "  --min-bend-radius <m> Minimum radius for ANY drawn corner/fillet, in metres\n"
              << "                        (also via MVB_MIN_BEND_RADIUS). Floors every corner's\n"
              << "                        policy radius; never tightens one.\n"
              << "  --toroid-mounting <vertical|horizontal>  How a toroid is mounted (ABT #1248;\n"
              << "                        sets MKF Settings toroid_mounting, default vertical; a\n"
              << "                        toroid base record's base.mounting still overrides it)\n"
              << "  -h, --help            Show this help\n";
}


// ABT #685 (Alf, 2026-08-17): "always produce step and 2D SVG together." The 2D view is the
// faster diagnosis for most 3D faults -- a wrap collision is usually a layout problem visible
// in the cross-section -- so the STEP never ships without it. Painted from the SAME enriched
// magnetic the 3D was built from, so the two cannot describe different coils.
//
// Failures here are reported and swallowed: a projection the Painter does not implement for
// this core (the whole-magnetic views are two-piece-set + rectangular-window only) must not
// cost the caller the STEP it asked for.
static void paintProjections(const MAS::Magnetic& rawMagnetic, bool useRealWinding,
                             const fs::path& stepPath) {
    using OpenMagnetics::Painter;
    using OpenMagnetics::PainterProjection;
    OpenMagnetics::Magnetic magnetic;
    try {
        auto enriched = useRealWinding
            ? mvb::magnetic_autocomplete_safe(nlohmann::json(rawMagnetic), true)
            : mvb::magnetic_autocomplete_safe(nlohmann::json(rawMagnetic));
        magnetic = enriched;
    } catch (const std::exception& e) {
        std::cerr << "2D painter: could not enrich for painting: " << e.what() << "\n";
        return;
    }
    const std::vector<std::pair<PainterProjection, std::string>> views{
        {PainterProjection::XY, ""},          // the winding-window cross-section
        {PainterProjection::XZ, ".top"},      // looking down the column: lanes and bumps
    };
    for (const auto& [projection, suffix] : views) {
        fs::path svg = stepPath;
        svg.replace_extension("");
        svg += suffix + std::string(".svg");
        try {
            fs::remove(svg);
            Painter painter(svg);
            painter.paint_magnetic(magnetic, projection);
            painter.export_svg();
            if (fs::exists(svg)) {
                std::cout << "Generated: " << svg << "\n";
            }
        } catch (const std::exception& e) {
            std::cerr << "2D painter (" << (suffix.empty() ? "XY" : "XZ") << "): " << e.what()
                      << "\n";
        }
    }
}

static bool processFile(const fs::path& inputPath, const fs::path& outputPath, bool useMkf,
                        bool useRealWinding, int segments, int coreSegments, bool coatedFootprint = false,
                        bool femReady = false) {
    try {
        // Read JSON
        std::ifstream f(inputPath);
        if (!f.is_open()) {
            std::cerr << "Error: Cannot open " << inputPath << "\n";
            return false;
        }
        
        json j;
        f >> j;
        
        // Patch dimensions
        mvb::patch_dimension_nominals(j);
        
        // Get magnetic data
        MAS::Magnetic magnetic;
        if (j.contains("magnetic")) {
            magnetic = j.at("magnetic").get<MAS::Magnetic>();
        } else {
            magnetic = j.get<MAS::Magnetic>();
        }
        
        // Generate STEP
        mvb::MagneticBuilder builder;
        const std::string format       = "step";
        const bool includeBobbin       = true;
        const double scale             = 1.0;     // exporters emit mm natively (ABT #317); 1000 double-scaled 1e6x
        const int symmetryPlanes       = 0;

        std::string result;
        // ABT #1215: the FEM product's terminal-lead lengths, measured with the settings the
        // conductors were drawn with, written next to the STEP once it has its final name.
        std::optional<std::map<std::string, mvb::ConductorBuilder::TerminalLeadLength>> leads;
        if (useRealWinding) {
            // Real winding requires the MKF wind (turn blocking on); no fallback.
            auto enriched = mvb::magnetic_autocomplete_safe(magnetic, true);
            // MVB_DUMP_TURNS=<path>: dump the ENRICHED turnsDescription (exactly what MKF's wind/
            // blocking placed and what ConductorBuilder consumes) BEFORE drawing, so turn placement
            // can be analysed even when the 3D build throws on a collision. This is the "MKF painter"
            // view: turn centres + per-winding wire, in window cross-section coordinates.
            if (const char* dump = std::getenv("MVB_DUMP_TURNS")) {
                nlohmann::json enr; OpenMagnetics::to_json(enr, enriched);
                nlohmann::json out;
                out["turnsDescription"] = enr["coil"].contains("turnsDescription")
                                              ? enr["coil"]["turnsDescription"] : nlohmann::json::array();
                out["windings"] = nlohmann::json::array();
                OpenMagnetics::Coil coil = enriched.get_coil();
                size_t wi = 0;
                for (const auto& w : enr["coil"]["functionalDescription"]) {
                    const auto wire = coil.resolve_wire(wi++);
                    nlohmann::json ww; ww["name"] = w["name"];
                    if (wire.get_conducting_diameter())
                        ww["conductingDiameter"] = OpenMagnetics::resolve_dimensional_values(wire.get_conducting_diameter().value());
                    if (wire.get_outer_diameter())
                        ww["outerDiameter"] = OpenMagnetics::resolve_dimensional_values(wire.get_outer_diameter().value());
                    if (wire.get_conducting_width())
                        ww["conductingWidth"] = OpenMagnetics::resolve_dimensional_values(wire.get_conducting_width().value());
                    if (wire.get_conducting_height())
                        ww["conductingHeight"] = OpenMagnetics::resolve_dimensional_values(wire.get_conducting_height().value());
                    out["windings"].push_back(ww);
                }
                std::ofstream df(dump); df << out.dump(1);
                std::cerr << "[dump-turns] wrote " << dump << "\n";
            }
            mvb::DrawConfig cfg{format, includeBobbin, scale, symmetryPlanes};
            cfg.useRealWindingGeometry = true;
            cfg.paintCoating = coatedFootprint;   // CD by default (ABT #1261); --coated => OUTER envelope
            cfg.femReady = femReady;              // --fem => slow one-piece/conformal meshable geometry
            if (segments >= 0) cfg.wirePolygonSegments = segments;
            // --core-segments defaults to EXACT for the FEM product. Faceting the core is not
            // free the way faceting the wire is: a round centre column discretised to an
            // inscribed n-gon loses 2.55% of its area at n=16, and the core's cross-section is
            // exactly what sets reluctance -- so the mesh knob would quietly bias inductance.
            // The winding needs faceting to stay off gmsh's periodic mesher; the core does not,
            // because the NURBS pass in buildAllNamed already strips periodicity from its
            // cylinders WITHOUT moving a point. Pass --core-segments explicitly to override.
            // CORE FACETING FOLLOWS WIRE FACETING, always (Alf, 2026-08-26). The two are a
            // phase-locked pair: a faceted wire's inner-corner vertices reach EXACTLY the radius
            // where they meet a matching faceted wall, keeping the conducting-vs-outer gap at
            // every azimuth. Forcing the core EXACT under --fem (the 2026-08-25 "A_e is physics"
            // change) broke that contract: measured on 03_buck seg=16, the wire polygon's
            // vertices landed at r=6.7255 -- the exact cylinder wall -- 50 um deeper than the
            // intended inner face, giving zero-distance copper-core tangency that OCC's fragment
            // then corrupts. The A_e bias of a faceted core is the FEM post-processor's business;
            // geometry consistency is this tool's.
            if (coreSegments >= 0)      cfg.corePolygonSegments = coreSegments;
            else if (segments >= 0)     cfg.corePolygonSegments = segments;
            result = builder.drawMagnetic(enriched, outputPath.parent_path().string(), cfg);
            if (femReady && !result.empty()) {
                leads = builder.measureTerminalLeadLengths(
                    enriched, cfg.paintCoating, cfg.femReady, cfg.wirePolygonSegments,
                    cfg.corePolygonSegments,
                    mvb::MagneticBuilder::declaredCoreCoatingThickness(enriched.get_core()));
            }
        } else if (useMkf) {
            try {
                // Use MVB++'s safe MKF wrapper to avoid Coil::wind() crashes on raw MAS files
                auto enriched = mvb::magnetic_autocomplete_safe(magnetic);
                result = builder.drawMagnetic(enriched, outputPath.parent_path().string(),
                                               format, includeBobbin, scale, symmetryPlanes);
            } catch (const std::exception& e) {
                std::cerr << "MKF enrichment failed: " << e.what() << "\n";
                std::cerr << "Falling back to raw magnetic...\n";
                result = builder.drawMagnetic(magnetic, outputPath.parent_path().string(),
                                               format, includeBobbin, scale, symmetryPlanes);
            }
        } else {
            result = builder.drawMagnetic(magnetic, outputPath.parent_path().string(),
                                           format, includeBobbin, scale, symmetryPlanes);
        }
        
        if (!result.empty()) {
            // drawMagnetic writes to <parent>/magnetic.step; rename to requested output path
            fs::path generated = outputPath.parent_path() / "magnetic.step";
            if (fs::exists(generated) && generated != outputPath) {
                fs::rename(generated, outputPath);
            }
            std::cout << "Generated: " << outputPath << "\n";
            if (leads) {
                const std::string sidecar =
                    mvb::MagneticBuilder::writeTerminalLeadSidecar(*leads, outputPath.string());
                std::cout << "Terminal leads: " << sidecar << "\n";
            }
            paintProjections(magnetic, useRealWinding, outputPath);
            return true;
        } else {
            std::cerr << "Failed to generate STEP for " << inputPath << "\n";
            return false;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error processing " << inputPath << ": " << e.what() << "\n";
        return false;
    }
}

// OCC runs its boolean edge/face intersections on a thread pool and each worker allocates its own
// Extrema sampling grid, so PEAK MEMORY SCALES WITH THREAD COUNT, not just with geometry. On the
// faceted product that is the difference between fitting a memory budget and not: 14_dab at
// --segments 16 peaked at 15.1 GB and died against a 16 GB cap. This box is shared with other
// agents' builds, and a job that takes the whole machine OOM-kills somebody else's work.
// MVB_OCC_THREADS caps the pool; unset or <= 0 leaves OCC's default (all cores).
// FACETED GEOMETRY RUNS ITS BOOLEANS SERIALLY (2026-09-06, ABT #1111). At --segments 12 the
// corpus crashed and hung inside OCCT's parallel boolean workers: single_switch 28 min in
// BOPAlgo_PaveFiller::PerformEF -> Extrema_GenExtPS::BuildGrid, then (after a Common was
// removed) SIGSEGV in NCollection_Array2::Resize from BOPAlgo_VertexFace under cut_bobbin;
// 14_dab the same Resize crash from PerformEF. With the pool at ONE thread the same
// single_switch run completes in 470 s. Every worker builds its own Extrema grid on the
// faceted B-spline strips, and that is what does not fit (this file already records 14_dab
// at 15.1 GB against a 16 GB cap at 16 segments). Exact round geometry (segments 0) keeps
// OCC's default pool: it has run the whole corpus clean on it. MVB_OCC_THREADS overrides both.
static void applyOccThreadCap(int segments) {
    const char* v = std::getenv("MVB_OCC_THREADS");
    int n = 0;
    if (v) n = std::atoi(v);
    else if (segments > 0) n = 1;
    if (n <= 0) return;
    OSD_ThreadPool::DefaultPool()->Init(n);
    std::cerr << "[occ] boolean thread pool capped at " << n << " thread(s)"
              << (v ? " (MVB_OCC_THREADS)" : " (faceted geometry: serial booleans)") << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }
    
    fs::path inputPath;
    fs::path outputPath;
    fs::path outputDir;
    bool useMkf = true;
    bool useRealWinding = false;
    // ABT #1261: the conductor is drawn at its CONDUCTING footprint by default. --coated asks
    // for the OUTER (insulation) envelope, which is what the web viewer draws.
    bool coatedFootprint = false;
    bool femReady = false;
    int segments = -1;  // -1 = builder default; 0 = exact analytic curves
    int coreSegments = -1;  // -1 = follow --segments (drawing) / stay EXACT (--fem)
    double minBendRadius = std::getenv("MVB_MIN_BEND_RADIUS")
                               ? std::atof(std::getenv("MVB_MIN_BEND_RADIUS")) : 0.0;

    // Parse arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-o" || arg == "--output") {
            if (++i < argc) outputPath = argv[i];
        } else if (arg == "-d" || arg == "--output-dir") {
            if (++i < argc) outputDir = argv[i];
        } else if (arg == "--no-mkf") {
            useMkf = false;
        } else if (arg == "--real") {
            useRealWinding = true;
        } else if (arg == "--copper") {
            // Kept so existing scripts keep working: the copper footprint is now the default.
        } else if (arg == "--coated") {
            coatedFootprint = true;
        } else if (arg == "--fem") {
            femReady = true;
        } else if (arg == "--skip-collision-check") {
            // Diagnosis only: builds the conductors even when two of them overlap, so the
            // overlap can be LOOKED AT instead of only read about in an exception. The STEP
            // that comes out is interpenetrating copper -- a picture, never a part.
            setenv("MVB_SKIP_COLLISION_CHECK", "1", 1);
        } else if (arg == "--core-segments") {
            if (i + 1 < argc) coreSegments = std::stoi(argv[++i]);
        } else if (arg == "--min-bend-radius") {
            if (i + 1 < argc) {
                minBendRadius = std::atof(argv[++i]);
                // ConductorBuilder::Options reads MVB_MIN_BEND_RADIUS as its default, so the
                // flag reaches every construction site without threading a parameter through
                // four signatures. setenv BEFORE any Options is constructed.
                setenv("MVB_MIN_BEND_RADIUS", argv[i], 1);
            }
        } else if (arg == "--toroid-mounting") {
            if (++i >= argc) { printUsage(argv[0]); return 1; }
            const std::string m = argv[i];
            if (m == "vertical")
                OpenMagnetics::Settings::GetInstance().set_toroid_mounting(MAS::OrientationEnum::VERTICAL);
            else if (m == "horizontal")
                OpenMagnetics::Settings::GetInstance().set_toroid_mounting(MAS::OrientationEnum::HORIZONTAL);
            else {
                std::cerr << "Error: --toroid-mounting takes vertical or horizontal, got '" << m << "'\n";
                return 1;
            }
        } else if (arg == "--segments") {
            if (++i < argc) segments = std::stoi(argv[i]);
        } else if (arg[0] != '-') {
            inputPath = arg;
        }
    }
    applyOccThreadCap(segments);   // after --segments: faceted geometry gets a serial pool

    if (inputPath.empty()) {
        std::cerr << "Error: No input file specified\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Determine output path
    if (outputPath.empty()) {
        fs::path dir = outputDir.empty() ? fs::current_path() : outputDir;
        fs::create_directories(dir);
        outputPath = dir / (inputPath.stem().string() + "_mvbpp.step");
    }
    
    bool success = processFile(inputPath, outputPath, useMkf, useRealWinding, segments, coreSegments,
                               coatedFootprint, femReady);
    return success ? 0 : 1;
}
