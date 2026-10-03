// mvbpp_cadbench <mas.json> [out.step] [segments] — CAD-only build benchmark + volume manifest.
// Moved verbatim from OMFEM tools/omfem_cadbench.cpp at faaf17a (ABT #1588, step 8); renamed; OMFEM's binary-freshness check dropped (it guards OMFEM's link to MVB++; this tool is part of MVB++). No logic changed.
//
// PHASE 2 of the 3D plan (Alf, 2026-08-10): make geometry generation faster while PROVING the
// result did not change. This tool exists because the meshing path must stay out of the way:
// it links MVB++ only (no gmsh, no mesher), does exactly what the CAD sweep does — enrich with
// MKF, build the named solids, optionally verify and export — and TIMES EACH PHASE, so an
// optimisation can be attributed instead of guessed at.
//
// Equivalence: --volumes writes a manifest of every named solid's volume (and the assembly
// total) to stdout. Byte-identical STEP is the strongest proof and the sweep's sha256 compare
// already does it; volumes are the fallback when a variant legitimately changes the B-Rep
// (different faceting for the fast web-frontend build) but must still describe the same metal.
//
// THREE INDEPENDENT AXES (Alf, 2026-08-10). Keep them independent — each has a different
// owner, and collapsing any two of them into one flag is how a viewer ends up drawing
// geometry nobody asked for:
//
//   PRODUCT   --viz  visualisation grade: per-piece bodies, no femReady welding, EXACT
//                    surfaces. Exact is not a compromise here: it is both faster and truer
//                    than faceting (a true cylinder is ONE face where a 12-gon is twelve,
//                    and the 12-gon under-fills the wire by 4.51%).
//             --fem  FEM grade: one continuous welded body per (winding, parallel).
//
//   SURFACES  exact (0) or faceted (N segments). FIXED exact for --viz; a FREE CHOICE for
//             --fem because it changes what the mesher has to chew on — that decision
//             belongs to the meshing phase.
//
//   WINDING   --winding real|ideal. A USER SETTING, not a product decision: the WebFrontend
//             exposes it for both 2D and 3D so a user can look at the real winding (leads,
//             pitch, dragbacks) or the idealised closed rings. Applies to BOTH products.
//
// Env: OMFEM_BENCH_REPEAT=n  repeat the build n times and report the best (warm-cache) time.
#include "mvb/MagneticBuilder.h"
#include "constructive_models/Magnetic.h"
#include "mvb/StepExporter.h"
#include "mvb/Utils.h"
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <nlohmann/json.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double since(const Clock::time_point& t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

double shapeVolume(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
                     "usage: %s <mas.json> [out.step] [segments=0] [--volumes] [--verify]\n"
                     "  segments: wire/core polygon segments; 0 = exact surfaces (FEM/CAD grade)\n",
                     argv[0]);
        return 2;
    }
    // Drives MVB++/MKF directly, bypassing the guarded bridge functions (ABT #683).
    const std::string masPath = argv[1];
    std::string outStep;
    int segments = 0;
    bool wantVolumes = false;
    bool wantVerify = false;
    bool fastVariant = false;   // --viz: visualisation grade (no femReady welding)
    bool realWinding = true;    // --winding real|ideal (the WebFrontend's user setting)
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--volumes") wantVolumes = true;
        else if (a == "--verify") wantVerify = true;
        else if (a == "--fast" || a == "--viz") { fastVariant = true; segments = 0; }
        else if (a == "--fem") fastVariant = false;
        else if (a == "--winding" && i + 1 < argc) {
            const std::string w = argv[++i];
            if (w == "real") realWinding = true;
            else if (w == "ideal") realWinding = false;
            else throw std::runtime_error("--winding takes 'real' or 'ideal', got '" + w + "'");
        }
        else if (!a.empty() && (a[0] == '-' || std::isdigit(static_cast<unsigned char>(a[0]))))
            segments = std::atoi(a.c_str());
        else outStep = a;
    }
    if (fastVariant) {
        segments = 0;   // the visualisation variant is EXACT by definition
    }
    const int repeat = std::getenv("OMFEM_BENCH_REPEAT")
                           ? std::max(1, std::atoi(std::getenv("OMFEM_BENCH_REPEAT"))) : 1;

    try {
        const auto tParse = Clock::now();
        std::ifstream in(masPath);
        if (!in) { std::fprintf(stderr, "cannot open %s\n", masPath.c_str()); return 2; }
        std::stringstream ss; ss << in.rdbuf();
        const nlohmann::json mas = nlohmann::json::parse(ss.str());
        const nlohmann::json magneticJson = mas.contains("magnetic") ? mas.at("magnetic") : mas;
        const double parseSecs = since(tParse);

        double bestEnrich = 1e30, bestBuild = 1e30, bestVerify = 0.0, bestExport = 0.0;
        size_t solidCount = 0;
        std::vector<std::pair<std::string, double>> volumes;
        double totalVolume = 0.0;
        int badShapes = 0;

        for (int pass = 0; pass < repeat; ++pass) {
            const auto tEnrich = Clock::now();
            // The enrich MUST use the same winding model as the build: with
            // useRealWindingGeometry the coil is wound with connection-lead turn blocking, and
            // ConductorBuilder then honours those positions verbatim. Enriching ideal and
            // building real produces a DIFFERENT layout — measured on 24_margin: an ideal
            // enrich collides at the gate where the real one is clean.
            // REAL WINDING IN BOTH VARIANTS. The visualisation product must show the same
            // physical winding as the FEM one — leads, pitch, dragbacks — otherwise the
            // WebFrontend draws idealised rings that do not exist. What --viz drops is the
            // femReady WELDING of the per-run pieces into one body, which only the mesher
            // needs; the pieces still touch, so the picture is identical.
            OpenMagnetics::Magnetic enriched =
                mvb::magnetic_autocomplete_safe(magneticJson,
                                                /*useRealWindingGeometry=*/true);
            const double enrichSecs = since(tEnrich);

            const auto tBuild = Clock::now();
            mvb::MagneticBuilder builder;
            // Match the sweep EXACTLY (mesh3d_from_mas): real winding + femReady, which fuses
            // each conductor's per-run pieces into ONE continuous body. That fusing is the
            // expensive part and the whole point of the FEM-grade variant; --fast drops it to
            // measure the visualisation-grade cost (per-piece compound, no welding).
            const std::vector<mvb::NamedShape> named =
                builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                      segments, segments, /*paintCoating=*/false,
                                      /*emitCoatingShells=*/false, /*includeInsulation=*/false,
                                      /*coreCoatingThickness=*/0.0,
                                      /*useRealWindingGeometry=*/true,
                                      /*femReady=*/!fastVariant);
            const double buildSecs = since(tBuild);

            double verifySecs = 0.0;
            if (wantVerify) {
                const auto tVerify = Clock::now();
                badShapes = 0;
                for (const auto& ns : named) {
                    if (ns.name.find(" terminal ") != std::string::npos) continue;
                    int openShells = 0, solids = 0;
                    for (TopExp_Explorer ex(ns.shape, TopAbs_SOLID); ex.More(); ex.Next()) {
                        ++solids;
                        for (TopExp_Explorer sh(ex.Current(), TopAbs_SHELL); sh.More(); sh.Next())
                            if (!TopoDS::Shell(sh.Current()).Closed() &&
                                !BRep_Tool::IsClosed(sh.Current()))
                                ++openShells;
                    }
                    const bool ok = BRepCheck_Analyzer(ns.shape).IsValid() && openShells == 0 &&
                                    solids >= 1 && shapeVolume(ns.shape) > 1e-15;
                    if (!ok) ++badShapes;
                }
                verifySecs = since(tVerify);
            }

            double exportSecs = 0.0;
            if (!outStep.empty()) {
                std::vector<mvb::NamedShape> solids;
                solids.reserve(named.size());
                for (const auto& ns : named)
                    if (ns.name.find(" terminal ") == std::string::npos) solids.push_back(ns);
                const auto tExport = Clock::now();
                std::ostringstream sink;
                std::streambuf* old = std::cout.rdbuf(sink.rdbuf());
                const bool ok = mvb::exportSTEP(solids, outStep);
                std::cout.rdbuf(old);
                if (!ok) throw std::runtime_error("exportSTEP failed");
                exportSecs = since(tExport);
            }

            if (pass + 1 == repeat) {   // manifest from the last pass
                solidCount = named.size();
                volumes.clear();
                totalVolume = 0.0;
                for (const auto& ns : named) {
                    const double v = shapeVolume(ns.shape);
                    volumes.emplace_back(ns.name, v);
                    totalVolume += v;
                }
            }
            bestEnrich = std::min(bestEnrich, enrichSecs);
            bestBuild = std::min(bestBuild, buildSecs);
            bestVerify = std::max(bestVerify, verifySecs);
            bestExport = std::max(bestExport, exportSecs);
        }

        std::fprintf(stderr,
                     "[cadbench] %s segments=%d%s: parse %.3fs | enrich(MKF) %.3fs | "
                     "build(MVB++) %.3fs | verify %.3fs | export %.3fs | TOTAL %.3fs | "
                     "%zu solids\n",
                     masPath.c_str(), segments, (std::string(fastVariant ? " VIZ" : " FEM") + (realWinding ? "/real" : "/ideal")).c_str(), parseSecs, bestEnrich, bestBuild, bestVerify,
                     bestExport, parseSecs + bestEnrich + bestBuild + bestVerify + bestExport,
                     solidCount);
        if (wantVerify) {
            std::fprintf(stderr, "[cadbench] verify: %d bad shape(s)\n", badShapes);
        }
        if (wantVolumes) {
            // Manifest: one line per named solid, then the assembly total. Stable ordering is
            // the builder's own emission order, so two runs diff line by line.
            for (const auto& [name, volume] : volumes) {
                std::printf("%.9g\t%s\n", volume * 1e9, name.c_str());   // mm^3
            }
            std::printf("%.9g\tTOTAL\n", totalVolume * 1e9);
        }
        return badShapes == 0 ? 0 : 1;
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "ERROR: %s\n", e.what());
        return 1;
    }
}
