// ABT #1265 [stepvalidity]: THE REAL TEST behind the self-intersection gate.
//
// The assembler used to run BOPAlgo_ArgumentAnalyzer on every welded conductor of every build.
// That check is a PROXY. The defect it stands for is this one, measured on 06_llc in August:
// a fused solid that BRepCheck calls VALID in memory, and that comes back INVALID once a STEP
// round-trip reconstructs it ('Primary parallel 0 [solid 25]', 3.652 mm^3 -- valid pre-scale,
// valid post-scale, invalid on read-back; no writer setting changes it, and ShapeFix and
// UnifySameDomain only "repair" it by deleting copper, 3.652 -> 3.633/3.536 mm^3).
//
// Paying the proxy on every build cost more than everything else in the assembler put together:
// 1094 s of 02_flyback's 1632 s at --segments 12, 590 s of it on a single conductor. So the
// proxy is now opt-in (MVB_WELD_SELFINT_GATE=1) and the real condition is asserted HERE, once,
// where a regression is supposed to be caught: build the geometry, write the STEP, READ IT
// BACK, and require every solid to be valid and the copper to survive.
//
// Reading back is the whole point: importSTEP goes through the same STEPCAFControl reader a
// consumer uses, so a solid that only "works" while it is still in the builder's memory fails
// here exactly as it would fail in Ansys or gmsh.
//
// Tags: [stepvalidity] everywhere, plus [slow] on the designs that take minutes. The default
// run covers the cheap designs; `./mvb_tests "[stepvalidity]"` from the repo ROOT runs them all.
#include <catch2/catch_test_macros.hpp>

#include "mvb/MagneticBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Utils.h"

#include "constructive_models/Magnetic.h"

#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

using json = nlohmann::json;

namespace {

json loadDesign(const std::string& name) {
    // Same search order as the other suites: the MAS examples the build pinned, then the
    // locally-patched complete fixtures.
    std::ifstream f(std::string(MAS_EXAMPLES_DIR) + "/" + name);
    if (!f.good()) f = std::ifstream(std::string(MAS_COMPLETE_DIR) + "/" + name);
    if (!f.good()) f = std::ifstream("tests/mas_complete_fixtures/" + name);
    REQUIRE(f.good());
    json j = json::parse(f);
    return j.contains("magnetic") ? j.at("magnetic") : j;
}

// Volume with ADAPTIVE integration, plus the accuracy the integrator reports. The default
// VolumeProperties integrates an analytic face exactly and a B-spline face by fixed-order
// quadrature, so a plain comparison across the round trip measures the INTEGRATION SCHEME and
// not the geometry -- this test's first run "failed" the CMC at ratio 1.0065 with every solid
// valid on both sides, which is exactly that artefact (the export re-expresses periodic faces
// as B-splines, by design). StepExporter.cpp:158 settled the same question the same way.
std::pair<double, double> totalVolume(const std::vector<mvb::NamedShape>& shapes) {
    const double eps = 1e-7;
    double v = 0.0, err = 0.0;
    for (const auto& ns : shapes) {
        if (ns.shape.IsNull()) continue;
        for (TopExp_Explorer e(ns.shape, TopAbs_SOLID); e.More(); e.Next()) {
            GProp_GProps g;
            const double relErr = std::abs(BRepGProp::VolumeProperties(e.Current(), g, eps));
            v += g.Mass();
            err += relErr * std::abs(g.Mass());
        }
    }
    return {v, err};
}

// The B-Rep's own resolution: every face is located only to within its tolerance, so a volume
// is defined only to about tolerance x surface area.
double brepNoise(const std::vector<mvb::NamedShape>& shapes) {
    double faceTol = 0.0, area = 0.0;
    for (const auto& ns : shapes) {
        if (ns.shape.IsNull()) continue;
        GProp_GProps sp;
        BRepGProp::SurfaceProperties(ns.shape, sp);
        area += std::abs(sp.Mass());
        for (TopExp_Explorer fx(ns.shape, TopAbs_FACE); fx.More(); fx.Next())
            faceTol = std::max(faceTol, BRep_Tool::Tolerance(TopoDS::Face(fx.Current())));
    }
    return faceTol * area;
}

// Build the FEM product for one design, write it, read it back, and hold the round trip to the
// two things that must survive it: every solid valid, and the copper still there.
void requireStepRoundTripValid(const std::string& design, int segments) {
    const json magneticJson = loadDesign(design);
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    const auto built = builder.buildAllNamed(enriched,
                                             /*includeBobbin=*/true,
                                             /*symmetryPlanes=*/0,
                                             /*wirePolygonSegments=*/segments,
                                             /*corePolygonSegments=*/segments,
                                             /*paintCoating=*/false,          // CD, ABT #1261
                                             /*emitCoatingShells=*/false,
                                             /*includeInsulation=*/false,
                                             /*coreCoatingThickness=*/0.0,
                                             /*useRealWindingGeometry=*/true,
                                             /*femReady=*/true);
    REQUIRE_FALSE(built.empty());

    // In-memory validity FIRST, so a failure says which side of the round trip broke.
    int inMemorySolids = 0;
    for (const auto& ns : built) {
        if (ns.shape.IsNull()) continue;
        for (TopExp_Explorer e(ns.shape, TopAbs_SOLID); e.More(); e.Next()) {
            ++inMemorySolids;
            INFO("in-memory solid of '" << ns.name << "' in " << design);
            REQUIRE(BRepCheck_Analyzer(e.Current()).IsValid());
        }
    }
    REQUIRE(inMemorySolids > 0);

    std::error_code ec;
    const auto dir = std::filesystem::temp_directory_path(ec) / "mvb_stepvalidity";
    std::filesystem::create_directories(dir, ec);
    const std::string path = (dir / (design + ".step")).string();
    // The FEM product is written the way the generator writes it: periodic surfaces re-expressed
    // as B-splines. That conversion is part of what the round trip has to survive.
    mvb::StepExportOptions opts;
    opts.nurbsPeriodicSolids = true;
    REQUIRE(mvb::exportSTEP(built, path, opts));

    const auto readBack = mvb::importSTEP(path);
    REQUIRE_FALSE(readBack.empty());

    // THE CONDITION THE GATE EXISTED FOR: valid after reconstruction, not merely in memory.
    int readBackSolids = 0;
    for (const auto& ns : readBack) {
        if (ns.shape.IsNull()) continue;
        for (TopExp_Explorer e(ns.shape, TopAbs_SOLID); e.More(); e.Next()) {
            ++readBackSolids;
            INFO("solid of '" << ns.name << "' read back from " << path);
            REQUIRE(BRepCheck_Analyzer(e.Current()).IsValid());
        }
    }
    REQUIRE(readBackSolids == inMemorySolids);

    // ...and no copper was lost or invented on the way. The file is written in millimetres, so
    // the metre-frame volume is scaled by 1e9 before the comparison. The bar is the model's own
    // resolution -- the two integrators' reported accuracies and the B-Rep's face tolerance --
    // NOT a chosen percentage. The failures this guards against (OCC dropping an operand, a
    // "repair" that deletes copper) are whole percent, orders of magnitude clear of it.
    const auto [vInM, errIn] = totalVolume(built);
    const auto [vOut, errOut] = totalVolume(readBack);
    REQUIRE(vInM > 0.0);
    REQUIRE(vOut > 0.0);
    const double vIn = vInM * 1e9;                       // metres^3 -> mm^3
    // The floor has THREE contributions, and the measured one dominates: writing the file.
    // A STEP is text with finite precision, so a re-read shape's points differ from the
    // in-memory ones in the last digits -- measured at 1.0 ppm on 03_buck (0.0138 of 13591
    // mm^3) and 1.6 ppm on the CMC (0.0072 of 4621.63 mm^3), with every solid valid on both
    // sides and the integrator/B-Rep terms an order of magnitude smaller. 10 ppm is therefore
    // the bar: ten times the noise actually observed, and a THOUSAND times below the smallest
    // real defect this test exists for (OCC dropping an operand, or a "repair" that deletes
    // copper, are whole-percent events -- 06_llc's was 3.652 -> 3.536 mm^3, 3.2 %).
    const double writePrecision = 1e-5 * vIn;
    const double noise = std::max({errIn * 1e9 + errOut, brepNoise(readBack), writePrecision});
    INFO(design << ": volume in " << vIn << " mm^3, out " << vOut << " mm^3, difference "
                << std::abs(vOut - vIn) << " mm^3 against a noise floor of " << noise);
    REQUIRE(std::abs(vOut - vIn) <= noise);
}

}   // namespace

TEST_CASE("STEP round trip keeps every solid valid: buck (rect wire)", "[stepvalidity]") {
    requireStepRoundTripValid("03_buck_inductor_pq3230_n95.json", 0);
}

TEST_CASE("STEP round trip keeps every solid valid: CMC (toroid)", "[stepvalidity]") {
    requireStepRoundTripValid("common_mode_choke_complete.json", 0);
}

// 06_llc_xfmr_eq4128_3c97 is THE design that produced the defect this test replaces: it shipped
// a BRepCheck-valid solid that read back invalid. It is also minutes long, hence [slow].
TEST_CASE("STEP round trip keeps every solid valid: 06_llc (the regression case)",
          "[stepvalidity][slow]") {
    requireStepRoundTripValid("06_llc_xfmr_eq4128_3c97.json", 0);
}

// Faceted geometry is the other half of the input class: at --segments 12 the welds run over
// polygonal faces, which is where BOPAlgo's sporadic self-intersections were seen (ABT #1111).
TEST_CASE("STEP round trip keeps every solid valid: CMC faceted", "[stepvalidity][slow]") {
    requireStepRoundTripValid("common_mode_choke_complete.json", 12);
}
