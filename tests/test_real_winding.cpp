// Real-winding geometry ([realwinding]): ONE continuous conductor per (winding, parallel)
// replacing the per-turn closed loops, with every MKF turn position honoured exactly.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <filesystem>
#include <set>
#include <functional>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "mvb/StepExporter.h"
#include "constructive_models/Magnetic.h"
#include "json.hpp"
#include "support/Settings.h"
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <TopoDS_Vertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Bnd_Box.hxx>
#include <BRepTools.hxx>
#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <GeomAbs_CurveType.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Torus.hxx>
#include "mvb/WireAssembler.h"
#include <gp_Dir.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BOPAlgo_CheckerSI.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt.hxx>
#include <cmath>
#include <array>
#include <cstdlib>
#include <limits>
#include <map>
#include <string>
#include <memory>
#include <fstream>
#include <sstream>
#include <numbers>
#include <variant>

using json = nlohmann::json;

namespace {

// Every artifact a test writes goes to <repo>/output, never the working directory: run from the
// build tree or the repo root, the STEPs land in one place and the project root stays clean
// (Alf, 2026-08-14: "can you clean the root of the project and just put all generated files in
// output?"). The path is derived from THIS source file, so it does not depend on the CWD.
std::string outputPath(const std::string& name) {
    const auto dir = std::filesystem::path{__FILE__}.parent_path().parent_path() / "output";
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    return (dir / name).string();
}

json loadFixture(const std::string& name) {
    std::ifstream f("testData/" + name);
    if (!f.good()) f = std::ifstream("tests/realwinding_fixtures/" + name);
    if (!f.good()) f = std::ifstream("tests/mas_complete_fixtures/" + name);
    REQUIRE(f.good());
    json j = json::parse(f);
    return j.contains("magnetic") ? j.at("magnetic") : j;
}

// A FEM conductor must be ONE connected solid. Counting solids is the check that was
// missing: a compound of DISCONNECTED per-turn bodies has the right name and a perfectly
// good volume, so name/volume assertions passed while the geometry was unusable — the
// winding could not be meshed at all ("1D mesh seems not to be forming a closed loop",
// ABT #332, root cause: the junction fuse ran at a 1e-7 fuzzy value and left invalid
// seams, so the builder silently fell back to a compound).
int solidCount(const TopoDS_Shape& s) {
    int n = 0;
    for (TopExp_Explorer e(s, TopAbs_SOLID); e.More(); e.Next()) ++n;
    return n;
}

// CONNECTED components among a shape's solids: two solids connect when they touch or
// overlap (BRepExtrema distance ~ 0). This -- not solidCount == 1 -- is the actual FEM
// requirement: the meshing fragment welds touching/overlapping same-region solids into one
// conformal region, so a compound whose pieces all touch meshes as ONE winding. Demanding a
// single BREP solid was retired 2026-08: the whole-spine single-body sweep that produced it
// folds onto itself at junctions (measured: 3.7 um between adjacent patches at ETD34's
// wrap->lead fillet; 4.7 um at every e138 racetrack corner), which NO element size can
// discretise, while welded/touching exact pieces mesh cleanly (ETD34: 3.66M tets, READY).
// A DISCONNECTED conductor is still rejected -- that was the original point of the check.
int connectedSolidComponents(const TopoDS_Shape& s) {
    std::vector<TopoDS_Shape> solids;
    for (TopExp_Explorer e(s, TopAbs_SOLID); e.More(); e.Next()) solids.push_back(e.Current());
    if (solids.size() <= 1) return (int)solids.size();
    std::vector<int> parent(solids.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = (int)i;
    std::function<int(int)> find = [&](int a) { while (parent[a] != a) a = parent[a] = parent[parent[a]]; return a; };
    // Bounding-box prefilter: all-pairs BRepExtrema on a ~150-solid conformal compound is ~11k
    // exact distance queries (~1 h); boxes farther apart than the touch tolerance can never touch.
    std::vector<Bnd_Box> boxes(solids.size());
    for (size_t i = 0; i < solids.size(); ++i) BRepBndLib::Add(solids[i], boxes[i]);
    for (size_t i = 0; i < solids.size(); ++i)
        for (size_t j = i + 1; j < solids.size(); ++j) {
            if (find((int)i) == find((int)j)) continue;
            if (boxes[i].Distance(boxes[j]) > 1e-6) continue;
            BRepExtrema_DistShapeShape d(solids[i], solids[j]);
            if (d.IsDone() && d.Value() < 1e-6) parent[find((int)i)] = find((int)j);
        }
    std::set<int> roots;
    for (size_t i = 0; i < solids.size(); ++i) roots.insert(find((int)i));
    return (int)roots.size();
}

bool pointStrictlyInsideShape(const TopoDS_Shape& shape, const gp_Pnt& p);

// The femReady round-wire conductor contract: a CONFORMAL mitre compound. One connected chain,
// every solid individually valid (BRepCheck), and consecutive solids abut WITHOUT volumetric
// overlap (their contact is coincident faces, which the meshing fragment welds conformally).
// Interpenetration is checked by CLASSIFYING A GRID over the junction region (the bbox
// intersection), NOT with BRepAlgoAPI_Common: OCC booleans on abutting BSpline pipe solids are
// the documented pathological input class this architecture exists to avoid, and Common ground
// for >10 min PER PAIR in 2d-extrema root-finding on the toroid's hole-threading pipes
// (measured via gdb backtrace). A point strictly inside BOTH neighbours = real overlap; the
// historical per-run compound overlaps (~0.5-1 mm^3 per joint) light up dozens of grid points.
void requireConformalConductor(const TopoDS_Shape& shape) {
    std::vector<TopoDS_Shape> solids;
    for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next()) solids.push_back(e.Current());
    REQUIRE(!solids.empty());
    REQUIRE(connectedSolidComponents(shape) == 1);
    for (const auto& s : solids) REQUIRE(BRepCheck_Analyzer(s).IsValid());
    int stride = std::max<int>(1, static_cast<int>(solids.size()) / 15);
    for (size_t i = 0; i + 1 < solids.size(); i += static_cast<size_t>(stride)) {
        Bnd_Box ba, bb;
        BRepBndLib::Add(solids[i], ba);
        BRepBndLib::Add(solids[i + 1], bb);
        double ax0, ay0, az0, ax1, ay1, az1, bx0, by0, bz0, bx1, by1, bz1;
        ba.Get(ax0, ay0, az0, ax1, ay1, az1);
        bb.Get(bx0, by0, bz0, bx1, by1, bz1);
        const double x0 = std::max(ax0, bx0), x1 = std::min(ax1, bx1);
        const double y0 = std::max(ay0, by0), y1 = std::min(ay1, by1);
        const double z0 = std::max(az0, bz0), z1 = std::min(az1, bz1);
        if (x0 >= x1 || y0 >= y1 || z0 >= z1) continue;   // disjoint boxes: nothing to probe
        constexpr int N = 6;
        int inBoth = 0;
        gp_Pnt firstHit;
        for (int gx = 0; gx < N; ++gx)
            for (int gy = 0; gy < N; ++gy)
                for (int gz = 0; gz < N; ++gz) {
                    const gp_Pnt p(x0 + (x1 - x0) * (gx + 0.5) / N,
                                   y0 + (y1 - y0) * (gy + 0.5) / N,
                                   z0 + (z1 - z0) * (gz + 0.5) / N);
                    if (pointStrictlyInsideShape(solids[i], p) && pointStrictlyInsideShape(solids[i + 1], p)) {
                        if (inBoth == 0) firstHit = p;
                        ++inBoth;
                    }
                }
        GProp_GProps gpa, gpb;
        BRepGProp::VolumeProperties(solids[i], gpa);
        BRepGProp::VolumeProperties(solids[i + 1], gpb);
        if (inBoth > 0) {   // offline forensics (CWD): the exact solids the probe flagged
            BRepTools::Write(solids[i], ("mitre_overlap_A_" + std::to_string(i) + ".brep").c_str());
            BRepTools::Write(solids[i + 1], ("mitre_overlap_B_" + std::to_string(i + 1) + ".brep").c_str());
        }
        INFO("mitre neighbours [" << i << "," << (i + 1) << "]: " << inBoth
             << " junction-grid points inside BOTH solids (first at ("
             << firstHit.X() << "," << firstHit.Y() << "," << firstHit.Z() << ")); A centroid=("
             << gpa.CentreOfMass().X() << "," << gpa.CentreOfMass().Y() << ","
             << gpa.CentreOfMass().Z() << ") vol=" << gpa.Mass() << "; B centroid=("
             << gpb.CentreOfMass().X() << "," << gpb.CentreOfMass().Y() << ","
             << gpb.CentreOfMass().Z() << ") vol=" << gpb.Mass());
        REQUIRE(inBoth == 0);
    }
}

double shapeVolume(const TopoDS_Shape& s) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props);
    return props.Mass();
}

double commonVolume(const TopoDS_Shape& a, const TopoDS_Shape& b) {
    BRepAlgoAPI_Common common(a, b);
    // A boolean that cannot run is NOT evidence of zero overlap: returning -1 here made
    // every 'overlap <= tol' assertion pass vacuously whenever OCC choked on the operands
    // (exact-quadric cylinder pairs do exactly that). Return a loud sentinel that fails
    // any sane tolerance instead.
    if (!common.IsDone()) return 1e9;
    return shapeVolume(common.Shape());
}


bool pointInsideShape(const TopoDS_Shape& shape, const gp_Pnt& p, double tol = 1e-9) {
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        BRepClass3d_SolidClassifier cls(TopoDS::Solid(exp.Current()), p, tol);
        if (cls.State() == TopAbs_IN || cls.State() == TopAbs_ON) return true;
    }
    return false;
}

// STRICT interior (TopAbs_IN only): the conformal junction probe must NOT count TopAbs_ON --
// points on the coincident abutment faces of a tangent mitre junction classify ON for BOTH
// neighbours (that contact IS the conformal design), and counting them read as interpenetration
// (16/216 false hits on the 12-turn toroid; a 24^3 strict-IN census of the same pair found 0).
bool pointStrictlyInsideShape(const TopoDS_Shape& shape, const gp_Pnt& p) {
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        BRepClass3d_SolidClassifier cls(TopoDS::Solid(exp.Current()), p, 1e-9);
        if (cls.State() == TopAbs_IN) return true;
    }
    return false;
}

// Exact OCCT self-intersection check on the emitted body. Consecutive-wrap CONTACT is by
// design (a spring); only genuine face-face interference counts. When the conductor is a
// fused/swept single solid, any residual self-intersection is a modelling defect.
bool hasSelfIntersections(const TopoDS_Shape& s) {
    BOPAlgo_CheckerSI checker;
    TopTools_ListOfShape args;
    args.Append(s);
    checker.SetArguments(args);
    checker.Perform();
    return checker.HasErrors();
}

// Facet-wedge bound: the cores are polygon-faceted (n-gon) approximations whose flats dip
// inside the true round window bore by sag = r*(1-cos(pi/n)); a lead legitimately ending
// at the true window border can therefore interpenetrate a facet by up to a wedge of
// volume ~ pi*wireRadius^2 * sag, at each of the two lead ends.
double coreFacetWedgeBound(double wireRadius, double borderRadius, int coreSegments) {
    double sag = borderRadius * (1.0 - std::cos(std::numbers::pi / coreSegments));
    return 2.0 * std::numbers::pi * wireRadius * wireRadius * sag;
}

// All-pairs boolean interference among named bodies (skipping the bobbin, which is
// deliberately cut to yield). Tolerance covers polygon-facet slivers of the cores.
void requireNoPairwiseOverlap(const std::vector<mvb::NamedShape>& named, double tol) {
    // Decompose every body into its solids once (with bboxes) and run booleans only on solid
    // pairs whose boxes actually come near: Common on whole CONFORMAL COMPOUNDS (~150 solids
    // since the femReady mitre default) took tens of minutes per pair and timed out the suite;
    // near-pair pruning keeps the check exact while touching only real contacts.
    struct Body {
        std::vector<TopoDS_Shape> solids;
        std::vector<Bnd_Box> boxes;
    };
    std::vector<Body> bodies(named.size());
    for (size_t i = 0; i < named.size(); ++i)
        for (TopExp_Explorer e(named[i].shape, TopAbs_SOLID); e.More(); e.Next()) {
            Bnd_Box bb;
            BRepBndLib::Add(e.Current(), bb);
            bodies[i].solids.push_back(e.Current());
            bodies[i].boxes.push_back(bb);
        }
    for (size_t i = 0; i < named.size(); ++i) {
        if (named[i].name.find("Bobbin") != std::string::npos) continue;
        for (size_t j = i + 1; j < named.size(); ++j) {
            if (named[j].name.find("Bobbin") != std::string::npos) continue;
            double v = 0.0;
            std::string verify;
            for (size_t a = 0; a < bodies[i].solids.size(); ++a)
                for (size_t b = 0; b < bodies[j].solids.size(); ++b) {
                    if (bodies[i].boxes[a].Distance(bodies[j].boxes[b]) > 1e-9) continue;
                    double vv = commonVolume(bodies[i].solids[a], bodies[j].solids[b]);
                    // Second opinion WITHOUT booleans: classify the common region's centroid in
                    // both bodies. OCC booleans on quadric pairs have been caught returning
                    // "empty" for genuinely overlapping solids (and fabricating the reverse), so
                    // an interference verdict must not rest on one algorithm.
                    if (vv > tol && vv < 1e8) {
                        BRepAlgoAPI_Common common(bodies[i].solids[a], bodies[j].solids[b]);
                        if (common.IsDone()) {
                            GProp_GProps gp_;
                            BRepGProp::VolumeProperties(common.Shape(), gp_);
                            const gp_Pnt c = gp_.CentreOfMass();
                            const bool inA = pointInsideShape(bodies[i].solids[a], c);
                            const bool inB = pointInsideShape(bodies[j].solids[b], c);
                            verify += std::string(" [") + std::to_string(a) + "," +
                                      std::to_string(b) + " centroid-in-A=" + (inA ? "yes" : "no") +
                                      " centroid-in-B=" + (inB ? "yes" : "no") + "]";
                            if (!(inA && inB)) vv = 0.0;   // boolean fabricated the overlap
                        }
                    }
                    v += vv;
                }
            INFO("pairwise overlap '" << named[i].name << "' vs '" << named[j].name
                                      << "' = " << v << verify);
            REQUIRE(v <= tol);
        }
    }
}

} // namespace

// REGRESSION (ABT #332): E 13/7/4 with 0.4 mm round wire — the EXACT geometry that failed.
// Verified to be a real guard: it FAILS at the old MVB_FUSE_FUZZY=1e-7 and passes at the 1e-6
// default. (A coarser E 32/16/9 fixture fuses fine even at 1e-7, so it would not have caught this.)
// A RECT-COLUMN core with round wire takes the analytic
// per-primitive path, whose junction fuse must weld into ONE solid. At the old 1e-7 fuzzy
// value the union left invalid seams, BRepCheck rejected it, and the builder SILENTLY
// returned a compound of disconnected per-turn solids — geometry that names and volumes
// alone cannot distinguish from the real thing, and that gmsh cannot mesh at all.
// ABT #615 / Alf (custom_magnetic 37, 2026-08-21): AN INTER-SECTION CONNECTION MUST HAVE A WAY
// OUT OF ITS OWN SECTION.
//
// E 16/6/5, Primary 22t x 2p interleaved. Primary section 0 is the FULL window height (6.5 mm)
// and its single layer is completely full: parallel 0 holds y = -3.094..-0.163 mm and parallel 1
// holds +0.163..+3.094, 20 turns at a 0.326 mm pitch. Parallel 0's last turn therefore sits in
// the MIDDLE of the stack, at the boundary with its sibling.
//
// MKF then routes that turn to the next section on an inter-section BAND at a window edge, choosing
// the edge with `exitTurn.y >= windowCentre`. At -0.163 that reads "below centre" and sends the band
// to the BOTTOM, so the wire must descend 1.58 mm back through its own turns 4..8 -- and by ABT
// #615's own rule a continuation reserves nothing in its own section, so nothing is blocked out of
// its way. The 3D gate refuses it as a hard bare-copper collision, correctly.
//
// There is no good edge: down crosses parallel 0's turns, up crosses parallel 1's. The section is
// full, so both are copper. The defect is that the layout does not pay for the escape at all --
// MKF reports the wind as fitting, then hands over a route that cannot be drawn. Fixing the edge
// choice alone is NOT enough (measured: it moves the failure to parallel 1); either the band must
// sit at the exit turn's own row, or a full section owing an inter-section exit must be one turn
// over capacity and refused.
//
// This test pins the CONTRACT: a design MKF accepts must be buildable. It was written failing
// ([!shouldfail]) and flipped on 2026-08-22 when the design first built CERTIFIED CLEAR (0 nm,
// every pair proven). What it took, in order (ABT #849): MKF's N-filar law (placement f5cb139f +
// distribution 1f897b5c) so every layer holds all parallels side by side and the exit turn is at
// the layer's end, not mid-stack; whole-winding parallel-order continuity (c318b004, kept through
// the L-shape revert 90c253bf); and on the MVB++ side the band chain riding one OD off the
// destination face (b16b2ff), the chain ending where the next wrap begins (no 180-degree fold at
// the crossing), and EXIT LANES -- parallels leave side by side on the +X face, the last wrap
// ending at its lane, the stub one OD off the face (Alf: 'as parallels they should be going out
// side by side'). Any regression in that chain shows here as the gate's own refusal.
TEST_CASE("Real winding: an interleaved N-filar design builds CERTIFIED CLEAR (cm37)",
          "[realwinding][connectivity][abt615][abt849]") {
    auto magneticJson = loadFixture("realwinding_interleaved_full_section_e16.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    // EXACTLY the CLI's path and config (mvbpp_step_generator --real), because that is where the
    // defect shows: buildAllNamed with polygonSegments=0 builds this design without complaint,
    // so a test written against it would have reported the bug fixed while the tool still refused.
    mvb::DrawConfig cfg{"step", /*includeBobbin=*/true, /*scale=*/1.0, /*symmetryPlanes=*/0};
    cfg.useRealWindingGeometry = true;
    cfg.paintCoating = true;
    // The gate is the assertion: it throws on any bare-copper overlap, so a clean build IS the
    // contract. Nothing is relaxed here -- see ABT #839, the gate is never weakened to pass.
    REQUIRE_NOTHROW(builder.drawMagnetic(enriched, outputPath(""), cfg));
}

TEST_CASE("Real winding: rect-column conductor fuses into ONE connected solid",
          "[realwinding][connectivity]") {
    auto magneticJson = loadFixture("realwinding_e138_rectcolumn.json");  // E 13/7/4, 6 turns, 0.4 mm round
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                       /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);

    int conductors = 0;
    for (const auto& ns : named) {
        // conductors are "<winding> parallel <k>"; terminal leads are
        // "<winding> parallel <k> terminal <n>" and are not solids — skip them.
        if (ns.name.find("parallel") == std::string::npos) continue;
        if (ns.name.find("terminal") != std::string::npos) continue;
        ++conductors;
        INFO("conductor: " << ns.name << "  solids=" << solidCount(ns.shape));
        REQUIRE(!ns.shape.IsNull());
        REQUIRE(shapeVolume(ns.shape) > 0.0);
        // ONE CONNECTED conductor (see connectedSolidComponents): welded/touching exact
        // pieces are FEM-equivalent to one solid; only DISCONNECTION is a failure.
        REQUIRE(connectedSolidComponents(ns.shape) == 1);
    }
    REQUIRE(conductors > 0);
}

// The web case Alf hit on 2026-08-13 (E 16/8/8, 19 turns of litz 31x0.1, 3 stacks, U winding
// order and ~2.92 mm of margin tape at the top). U alone built, margin alone built, and the two
// together died in ConductorBuilder with the entrance lead running through a turn's wrap — the
// 3D view showed a core and a bobbin with no copper at all (MKF ABT #682: the overflowing layer
// was centred on a span it did not fit and landed on the lead's own row; it now gives back the
// deeper reservation instead). Kept as a fixture because it is the combination, not either half,
// that regressed: a design carrying BOTH is the only thing that would have caught it.
// TEMP sweep: every MAS example through the real-winding builder, checked for what a mesh needs
// — every solid valid (BRepCheck), ONE connected component per conductor, and no volumetric
// interpenetration at the junctions.
TEST_CASE("Tmp mesh sweep", "[meshsweep]") {
    namespace fs = std::filesystem;
    std::vector<fs::path> files;
    for (const auto& e : fs::directory_iterator(MAS_EXAMPLES_DIR)) {
        if (e.path().extension() == ".json" && std::isdigit(e.path().filename().string()[0])) {
            files.push_back(e.path());
        }
    }
    std::sort(files.begin(), files.end());
    for (const auto& f : files) {
        std::ifstream in(f.string());
        json j = json::parse(in);
        json magneticJson = j.contains("magnetic") ? j.at("magnetic") : j;
        std::string name = f.filename().string();
        try {
            auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);
            mvb::MagneticBuilder builder;
            auto named = builder.buildAllNamed(enriched, false, 0,
                                               mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                               mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false,
                                               false, 0.0, true, true);
            int conductors = 0, bad = 0;
            std::string why;
            for (const auto& ns : named) {
                if (ns.name.find(" parallel ") == std::string::npos) continue;
                std::vector<TopoDS_Shape> solids;
                if (!ns.shape.IsNull())
                    for (TopExp_Explorer e(ns.shape, TopAbs_SOLID); e.More(); e.Next())
                        solids.push_back(e.Current());
                if (solids.empty()) continue;   // coating shells / empty helpers, not the copper
                ++conductors;
                if (connectedSolidComponents(ns.shape) != 1) { ++bad; why += " disconnected"; }
                for (const auto& sol : solids)
                    if (!BRepCheck_Analyzer(sol).IsValid()) { ++bad; why += " invalid"; break; }
            }
            std::cout << "[mesh] " << name << " conductors=" << conductors
                      << (bad == 0 ? "  OK (valid, connected)" : ("  PROBLEM:" + why)) << std::endl;
        }
        catch (const std::exception& e) {
            std::cout << "[mesh] " << name << "  THREW: " << std::string(e.what()).substr(0, 90) << std::endl;
        }
    }
}

TEST_CASE("Real winding: U order with margin tape builds a continuous conductor",
          "[realwinding][abt682]") {
    auto magneticJson = loadFixture("realwinding_u_order_margin_e16.json");
    // The fixture is the design as the web hands it over: already wound, U on the bobbin's
    // winding window, margin on the section.
    REQUIRE(magneticJson.at("coil").at("bobbin").at("processedDescription")
                        .at("windingWindows").at(0).at("windingOrder") == "U");

    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       /*paintCoating=*/true, /*emitCoatingShells=*/false,
                                       /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);

    // The layers must connect HORIZONTALLY (MKF ABT #683): U exists so the next layer starts
    // level with the turn it connects to. Read it off the geometry MVB++ actually builds from.
    {
        auto& enrichedCoil = enriched.get_mutable_coil();
        auto enrichedTurnsOptional = enrichedCoil.get_turns_description();
        REQUIRE(enrichedTurnsOptional);
        const auto& enrichedTurns = enrichedTurnsOptional.value();
        std::vector<std::string> layerOrder;
        std::map<std::string, std::vector<double>> heightsPerLayer;
        for (const auto& turn : enrichedTurns) {
            if (!turn.get_layer()) continue;
            const std::string layerName = turn.get_layer().value();   // BY VALUE: the optional comes back by value, so .value() would dangle
            if (!heightsPerLayer.count(layerName)) layerOrder.push_back(layerName);
            heightsPerLayer[layerName].push_back(turn.get_coordinates()[1]);
        }
        REQUIRE(layerOrder.size() >= 2);
        for (size_t i = 0; i + 1 < layerOrder.size(); ++i) {
            double arrival = heightsPerLayer[layerOrder[i]].back();
            double landing = heightsPerLayer[layerOrder[i + 1]].front();
            INFO(layerOrder[i] << " leaves at " << arrival * 1000 << " mm, "
                 << layerOrder[i + 1] << " starts at " << landing * 1000 << " mm");
            CHECK(std::abs(arrival - landing) < 0.0002);
        }
    }

    const mvb::NamedShape* conductor = nullptr;
    for (const auto& ns : named) {
        if (ns.name == "Primary parallel 0") {
            REQUIRE(conductor == nullptr);
            conductor = &ns;
        }
    }
    REQUIRE(conductor != nullptr);
    REQUIRE(!conductor->shape.IsNull());
    REQUIRE(shapeVolume(conductor->shape) > 0.0);
    REQUIRE(connectedSolidComponents(conductor->shape) == 1);

    // The margin is geometry, not decoration (MKF ABT #676): no copper above the tape's inner
    // face. Read the face from the design rather than hard-coding it, so the fixture stays the
    // source of truth.
    const auto& window = magneticJson.at("coil").at("bobbin").at("processedDescription")
                                     .at("windingWindows").at(0);
    const double windowTop = window.at("coordinates").at(1).get<double>()
                           + window.at("height").get<double>() / 2;
    double margin = 0.0;
    for (const auto& section : magneticJson.at("coil").at("sectionsDescription")) {
        if (section.at("type") != "conduction" || !section.contains("margin")) continue;
        const auto& m = section.at("margin");
        margin = std::max(margin, m.is_array() ? m.at(0).get<double>()
                                               : m.at("topOrLeftWidth").get<double>());
    }
    REQUIRE(margin > 0.0);
    Bnd_Box box;
    BRepBndLib::Add(conductor->shape, box);
    double xMin, yMin, zMin, xMax, yMax, zMax;
    box.Get(xMin, yMin, zMin, xMax, yMax, zMax);
    // Bounding box is in mm; a faceting tolerance of one polygon sagitta is expected on round wire.
    INFO("copper top " << yMax << " mm against the margin's inner face "
                       << (windowTop - margin) * 1000 << " mm");
    REQUIRE(yMax <= (windowTop - margin) * 1000 + 0.05);
}

TEST_CASE("Real winding: single-parallel PQ33 becomes one continuous conductor",
          "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_round_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       /*paintCoating=*/true, /*emitCoatingShells=*/false,
                                       /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);

    // Exactly one conductor for the single (winding, parallel); no per-turn ring solids.
    const mvb::NamedShape* conductor = nullptr;
    for (const auto& ns : named) {
        INFO(ns.name);
        REQUIRE_THAT(ns.name, !Catch::Matchers::ContainsSubstring(" turn "));
        if (ns.name == "Primary parallel 0") {
            REQUIRE(conductor == nullptr);
            conductor = &ns;
        }
    }
    REQUIRE(conductor != nullptr);
    REQUIRE(!conductor->shape.IsNull());
    REQUIRE(shapeVolume(conductor->shape) > 0.0);
    // CONNECTIVITY: one connected component — a continuous conductor, not a bag of loose
    // turns. Touching/overlapping pieces count as connected (the meshing fragment welds
    // them into one conformal region); see connectedSolidComponents.
    REQUIRE(connectedSolidComponents(conductor->shape) == 1);

    // The "nothing moves" regression guard: every MKF turn station lies INSIDE the
    // conductor's copper. Helical wraps pass the exact station at their start phase, so
    // probe the full station RING (a thin torus at (r, y) around the column axis) — it
    // must intersect the copper.
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    REQUIRE(!turnsOpt->empty());
    // The station RING (circle at the turn's exact (r, y)) must enter the copper: its
    // minimum distance to the conductor is 0 when the wire centreline passes through the
    // station somewhere (extrema is far more robust than boolean common on compounds).
    for (const auto& turn : *turnsOpt) {
        const auto& c = turn.get_coordinates();
        REQUIRE(c.size() >= 2);
        gp_Circ stationRing(gp_Ax2(gp_Pnt(0.0, c[1], 0.0), gp_Dir(0, 1, 0)), c[0]);
        TopoDS_Edge ringEdge = BRepBuilderAPI_MakeEdge(stationRing).Edge();
        BRepExtrema_DistShapeShape dist(ringEdge, conductor->shape);
        REQUIRE(dist.IsDone());
        INFO("turn " << turn.get_name() << " station (r=" << c[0] << ", y=" << c[1]
                     << ") ring-distance=" << dist.Value());
        REQUIRE(dist.Value() <= 1e-9);
    }

    // "Absolutely no body collides with another or itself":
    // (a) all-pairs boolean interference across every emitted body (tolerance 1e-10 m^3
    //     for the cores' polygon-facet slivers — same class [battery] tolerates at 1e-7);
    // (b) exact OCCT self-intersection check when the conductor fused to a single solid.
    //     When OCCT's tangent-contact fuse defect forces the per-run compound fallback,
    //     the consecutive pieces legitimately overlap at their junctions (the wire's own
    //     crossovers) — for that case the capsule gate + station probes + (a) are the
    //     collision guarantee.
    // Tolerance: the facet-wedge bound for this fixture (wire radius 0.4795 mm, window
    // border ~13.75 mm, 16-gon cores).
    requireNoPairwiseOverlap(named, coreFacetWedgeBound(0.0004795, 0.01375, 16));
    int conductorSolids = 0;
    for (TopExp_Explorer exp(conductor->shape, TopAbs_SOLID); exp.More(); exp.Next())
        ++conductorSolids;
    if (conductorSolids == 1) {
        REQUIRE(!hasSelfIntersections(conductor->shape));
    }
}

TEST_CASE("Real winding: flag off keeps the per-turn loops unchanged", "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_round_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/false);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched);

    size_t turnSolids = 0;
    bool sawConductor = false;
    for (const auto& ns : named) {
        if (ns.name.find(" turn ") != std::string::npos) ++turnSolids;
        if (ns.name == "Primary parallel 0") sawConductor = true;
    }
    REQUIRE(turnSolids == 16);   // 16 turns x 1 parallel, one closed loop each (2 layers)
    REQUIRE(!sawConductor);
}

TEST_CASE("Real winding: two parallels become two independent conductors", "[realwinding]") {
    auto magneticJson = loadFixture("round_2p_1layer.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);

    const mvb::NamedShape* p0 = nullptr;
    const mvb::NamedShape* p1 = nullptr;
    for (const auto& ns : named) {
        if (ns.name == "Primary parallel 0") p0 = &ns;
        if (ns.name == "Primary parallel 1") p1 = &ns;
    }
    REQUIRE(p0 != nullptr);
    REQUIRE(p1 != nullptr);
    REQUIRE(shapeVolume(p0->shape) > 0.0);
    REQUIRE(shapeVolume(p1->shape) > 0.0);

    // Parallel conductors are fully independent copper bodies: zero overlap (contact only),
    // no pairwise interference anywhere, and no self-intersections when fused single.
    double v = commonVolume(p0->shape, p1->shape);
    INFO("parallel-parallel common volume = " << v);
    REQUIRE(v <= 1e-12);
    requireNoPairwiseOverlap(named, 1e-10);
    for (const mvb::NamedShape* c : {p0, p1}) {
        int solids = 0;
        for (TopExp_Explorer exp(c->shape, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
        if (solids == 1) REQUIRE(!hasSelfIntersections(c->shape));
    }
}

TEST_CASE("Real winding: multi-layer multi-parallel builds collision-free",
          "[realwinding]") {
    // HISTORY — this fixture characterized two generations of a collision, both fixed. First MKF
    // drew every parallel's terminal lead at the SAME edge row (coincident copper) until ABT
    // #229/#240 gave each parallel its own row. The build then still collided: the U layer links
    // landed LEVEL with the previous layer's last turn, so the parallels' landing revolutions
    // overlapped. The ABT #608 final form fixed that too — MKF places each non-first U layer's
    // first station below the tangential arrival (as far as the window allows) and the landing
    // wrap descends, chunk included — so the 8t x 2p multi-layer U fixture now builds valid,
    // meshable copper (verified ALL WATERTIGHT in the full battery). The gate throwing here again
    // means a REGRESSION in one of those two fixes.
    auto magneticJson = loadFixture("realwinding_round_2p.json");   // 8t x 2p -> multi-layer U
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true);

    const mvb::NamedShape* p0 = nullptr;
    const mvb::NamedShape* p1 = nullptr;
    for (const auto& ns : named) {
        if (ns.name == "Primary parallel 0") p0 = &ns;
        if (ns.name == "Primary parallel 1") p1 = &ns;
    }
    REQUIRE(p0 != nullptr);
    REQUIRE(p1 != nullptr);
    REQUIRE(shapeVolume(p0->shape) > 0.0);
    REQUIRE(shapeVolume(p1->shape) > 0.0);
    // Same guarantee the single-layer 2p case asserts: the parallels are independent copper.
    double v = commonVolume(p0->shape, p1->shape);
    INFO("parallel-parallel common volume = " << v);
    REQUIRE(v <= 1e-12);
}

namespace {

// Test-side mirror of MagneticBuilder's bobbin resolution (getBobbinProcessed +
// patchBobbinDimensions) so probes use the same column dimensions as the builder.
MAS::CoreBobbinProcessedDescription bobbinPdOf(const OpenMagnetics::Magnetic& m) {
    auto bobbinVar = m.get_coil().get_bobbin();
    auto* b = std::get_if<OpenMagnetics::Bobbin>(&bobbinVar);
    REQUIRE(b != nullptr);
    auto pd = b->get_processed_description();
    REQUIRE(pd.has_value());
    MAS::CoreBobbinProcessedDescription bobbinPd = *pd;
    if (bobbinPd.get_column_width().value_or(0.0) <= 0.0) {
        auto corePd = m.get_core().get_processed_description();
        REQUIRE(corePd.has_value());
        REQUIRE(!corePd->get_columns().empty());
        const auto& col = corePd->get_columns()[0];
        double wall = bobbinPd.get_wall_thickness();
        if (std::isnan(wall) || wall < 0.0) wall = 0.0;
        bobbinPd.set_column_width(col.get_width() / 2.0 + wall);
        bobbinPd.set_column_depth(col.get_depth() / 2.0 + wall);
        bobbinPd.set_column_shape(col.get_shape());
    }
    return bobbinPd;
}

// A crossing must lie ON THE CENTERLINE of the copper: the crossing itself and four
// probes at ±0.99·wireRadius along the two directions perpendicular to the wire's travel
// must all be inside the solid. If the centerline missed the crossing by more than
// 0.01·wireRadius, the probe opposite the offset would fall outside the pipe. (A
// boundary-distance check is unsound here: at the conductor's free ends and at run
// junctions internal cap faces pass exactly through the crossing.)
void requireCrossingOnCenterline(const TopoDS_Shape& conductor, const gp_Pnt& crossing,
                                 double wireRadius, const gp_Dir& perpA, const gp_Dir& perpB,
                                 const std::string& what) {
    INFO(what << " at (" << crossing.X() << "," << crossing.Y() << "," << crossing.Z()
              << "), wire radius " << wireRadius);
    REQUIRE(pointInsideShape(conductor, crossing, 1e-9));
    // The emitted section is an INSCRIBED n-gon, so its flats lie at the apothem
    // r*cos(pi/n), not at r (0.981*r at the default 16 segments). Probing at 0.99*r would
    // land outside the copper along a face normal even for a perfectly centred wire, so the
    // offset is measured against the apothem -- still a tight centring bound (the probe sits
    // within ~1% of the real material boundary), just an honest one for a faceted section.
    // ABT #860: segments <= 0 is the EXACT CIRCLE -- no flats, so the apothem IS the radius.
    // Dividing pi by a zero segment count gave cos(inf) = NaN, and a NaN probe is inside
    // nothing: the moment the default became analytic this helper failed four designs while
    // the geometry was correct.
    const double sectionApothem =
        mvb::DEFAULT_WIRE_POLYGON_SEGMENTS > 0
            ? wireRadius * std::cos(std::numbers::pi / mvb::DEFAULT_WIRE_POLYGON_SEGMENTS)
            : wireRadius;
    for (const gp_Dir* d : {&perpA, &perpB}) {
        for (double sgn : {1.0, -1.0}) {
            gp_Pnt probe(crossing.XYZ() + d->XYZ() * (sgn * 0.99 * sectionApothem));
            INFO("perpendicular probe at (" << probe.X() << "," << probe.Y() << ","
                                            << probe.Z() << ")");
            REQUIRE(pointInsideShape(conductor, probe, 1e-9));
        }
    }
}

// Outer-footprint wire radius of an MKF-enriched turn (turn.dimensions = outer w/h).
double turnWireRadius(const MAS::Turn& turn) {
    auto dims = turn.get_dimensions();
    REQUIRE(dims.has_value());
    REQUIRE(dims->size() >= 2);
    return std::min((*dims)[0], (*dims)[1]) / 2.0;
}

const mvb::NamedShape* findConductor(const std::vector<mvb::NamedShape>& named,
                                     const std::string& name) {
    const mvb::NamedShape* found = nullptr;
    for (const auto& ns : named) {
        REQUIRE_THAT(ns.name, !Catch::Matchers::ContainsSubstring(" turn "));
        if (ns.name == name) {
            REQUIRE(found == nullptr);
            found = &ns;
        }
    }
    REQUIRE(found != nullptr);
    REQUIRE(!found->shape.IsNull());
    REQUIRE(shapeVolume(found->shape) > 0.0);
    return found;
}

} // namespace

TEST_CASE("Real winding: rectangular-column E core zigzag racetrack conductor",
          "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_rect_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");

    // Every MKF crossing must lie ON the copper. For a rectangular column the crossing
    // sits at the -Z transition-face centre: (0, y, -(x + columnDepth - columnWidth)).
    auto bobbinPd = bobbinPdOf(enriched);
    REQUIRE(bobbinPd.get_column_shape() == MAS::ColumnShape::RECTANGULAR);
    double zoff = bobbinPd.get_column_depth() - bobbinPd.get_column_width().value();
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    for (const auto& turn : *turnsOpt) {
        const auto& c = turn.get_coordinates();
        REQUIRE(c.size() >= 2);
        requireCrossingOnCenterline(conductor->shape, gp_Pnt(0.0, c[1], -(c[0] + zoff)),
                                    turnWireRadius(turn), gp_Dir(0, 1, 0), gp_Dir(0, 0, 1),
                                    "crossing " + turn.get_name());
    }

    // E-core window walls are planar (no facet sag): tangent contact only.
    requireNoPairwiseOverlap(named, 1e-10);
    // CONFORMAL CONTRACT (femReady round wire): the conductor is a mitre-jointed compound --
    // one CONNECTED chain of individually valid solids whose neighbours abut without volumetric
    // overlap. (The old fused-ONE-solid contract is retired: OCC booleans on winding chains are
    // the documented self-interference failure class -- ABT #490; conformal-by-construction is
    // the architecture. BOPAlgo_CheckerSI is NOT run: coincident abutting faces are the intended
    // conformal contact that gmsh's fragment+glue welds.)
    requireConformalConductor(conductor->shape);
}

TEST_CASE("Real winding: oblong-column EP core stadium conductor", "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_oblong_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");

    // Oblong crossing sits at the -Z cap apex: same z = -(x + D - W) mapping.
    auto bobbinPd = bobbinPdOf(enriched);
    REQUIRE(bobbinPd.get_column_shape() == MAS::ColumnShape::OBLONG);
    double zoff = std::max(0.0, bobbinPd.get_column_depth() -
                                    bobbinPd.get_column_width().value());
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    for (const auto& turn : *turnsOpt) {
        const auto& c = turn.get_coordinates();
        REQUIRE(c.size() >= 2);
        requireCrossingOnCenterline(conductor->shape, gp_Pnt(0.0, c[1], -(c[0] + zoff)),
                                    turnWireRadius(turn), gp_Dir(0, 1, 0), gp_Dir(0, 0, 1),
                                    "crossing " + turn.get_name());
    }

    requireNoPairwiseOverlap(named, coreFacetWedgeBound(0.00028, 0.0089, 16));
    int solids = 0;
    for (TopExp_Explorer exp(conductor->shape, TopAbs_SOLID); exp.More(); exp.Next())
        ++solids;
    if (solids == 1) REQUIRE(!hasSelfIntersections(conductor->shape));
}

TEST_CASE("Real winding: toroidal conductor threads the exact inner and outer crossings",
          "[realwinding]") {
    for (auto mounting : {MAS::OrientationEnum::VERTICAL, MAS::OrientationEnum::HORIZONTAL}) {
    OpenMagnetics::SettingsGuard<MAS::OrientationEnum> mountingGuard(
        OpenMagnetics::Settings::GetInstance(), &OpenMagnetics::Settings::get_toroid_mounting,
        &OpenMagnetics::Settings::set_toroid_mounting, mounting);
    auto magneticJson = loadFixture("realwinding_toroid.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    const auto mountFrame = mvb::MagneticBuilder::toroidMountingFrameOf(enriched);
    auto placed = [&](double cx, double cy) { return gp_Pnt(cx, 0.0, cy).Transformed(mountFrame.toExported); };
    const gp_Dir inPlaneA = gp_Dir(1, 0, 0).Transformed(mountFrame.toExported);
    const gp_Dir inPlaneB = gp_Dir(0, 0, 1).Transformed(mountFrame.toExported);

    mvb::MagneticBuilder builder;
    // Exact core (segments=0): the wall-adjacent ring's tubes touch the bore tangentially
    // along their whole length, so any polygon-faceted bore would interpenetrate them by
    // its facet sag; the exact annulus makes true tangency testable at boolean tolerance.
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");

    // ABT #1248: the assembly is placed for its mounting. A hole-plane crossing (cx, cy) is the build
    // point (cx, 0, cy) moved by the mounting's rigid motion (VERTICAL default: the MAS frame, ring in
    // XY, turned about Z so the terminals sit at the bottom; before #1248 it was always the MAS frame
    // with no turn, (cx, cy, 0)). Checked in both mountings: the inner crossing of every turn, the
    // outer (additionalCoordinates) of every wrapped turn (the last entry's outer is not wrapped).
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    const auto& turns = *turnsOpt;
    for (size_t i = 0; i < turns.size(); ++i) {
        const auto& c = turns[i].get_coordinates();
        REQUIRE(c.size() >= 2);
        double wr = turnWireRadius(turns[i]);
        requireCrossingOnCenterline(conductor->shape, placed(c[0], c[1]), wr,
                                    inPlaneA, inPlaneB,
                                    "turn " + turns[i].get_name() + " inner crossing");
        if (i + 1 < turns.size()) {
            auto add = turns[i].get_additional_coordinates();
            REQUIRE(add.has_value());
            REQUIRE(!add->empty());
            requireCrossingOnCenterline(conductor->shape,
                                        placed((*add)[0][0], (*add)[0][1]), wr,
                                        inPlaneA, inPlaneB,
                                        "turn " + turns[i].get_name() + " outer crossing");
        }
    }

    // Exact bore: wall contact is true tangency, zero interference within OCCT booleans.
    requireNoPairwiseOverlap(named, 1e-12);
    // CONFORMAL CONTRACT (femReady round wire): a connected mitre compound of valid solids,
    // no volumetric overlap between neighbours (see requireConformalConductor). The historical
    // single-MakePipe body is retired for femReady -- conformal-by-construction, ABT #490.
    requireConformalConductor(conductor->shape);

    // FEM terminal faces: the two free ends of the conductor are flat PLANAR discs (the
    // swept lead cylinder's end cap, no sphere), so downstream FEM can assign a current
    // BC on a planar surface. The wire is round and the whole conductor is otherwise
    // cylinders/tori/revolves/spheres, so the only planar faces are the two terminals.
    int planarFaces = 0;
    for (TopExp_Explorer exp(conductor->shape, TopAbs_FACE); exp.More(); exp.Next()) {
        BRepAdaptor_Surface sa(TopoDS::Face(exp.Current()));
        if (sa.GetType() == GeomAbs_Plane) ++planarFaces;
    }
    INFO("planar (terminal) faces on the conductor = " << planarFaces);
    REQUIRE(planarFaces >= 2);
    }
}

TEST_CASE("Real winding: toroidal RECTANGULAR wire threads the crossings", "[realwinding]") {
    // Rectangular wire on a toroid (single layer): each turn is built from per-primitive rect
    // solids (prisms + revolved poloidal elbows) oriented on the local AZIMUTHAL axis, then fused.
    // Every MKF inner/outer crossing must still lie on the copper (the section's inscribed circle
    // is min(w,h)/2 = turnWireRadius, so the 0.99*r probes stay inside whatever the orientation).
    for (auto mounting : {MAS::OrientationEnum::VERTICAL, MAS::OrientationEnum::HORIZONTAL}) {
    OpenMagnetics::SettingsGuard<MAS::OrientationEnum> mountingGuard(
        OpenMagnetics::Settings::GetInstance(), &OpenMagnetics::Settings::get_toroid_mounting,
        &OpenMagnetics::Settings::set_toroid_mounting, mounting);
    auto magneticJson = loadFixture("realwinding_toroid_rect.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    const auto mountFrame = mvb::MagneticBuilder::toroidMountingFrameOf(enriched);
    auto placed = [&](double cx, double cy) { return gp_Pnt(cx, 0.0, cy).Transformed(mountFrame.toExported); };
    const gp_Dir inPlaneA = gp_Dir(1, 0, 0).Transformed(mountFrame.toExported);
    const gp_Dir inPlaneB = gp_Dir(0, 0, 1).Transformed(mountFrame.toExported);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(shapeVolume(conductor->shape) > 0.0);

    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    const auto& turns = *turnsOpt;
    for (size_t i = 0; i < turns.size(); ++i) {
        const auto& c = turns[i].get_coordinates();
        REQUIRE(c.size() >= 2);
        double wr = turnWireRadius(turns[i]);
        requireCrossingOnCenterline(conductor->shape, placed(c[0], c[1]), wr,
                                    inPlaneA, inPlaneB,
                                    "turn " + turns[i].get_name() + " inner crossing");
        if (i + 1 < turns.size()) {
            auto add = turns[i].get_additional_coordinates();
            REQUIRE(add.has_value());
            REQUIRE(!add->empty());
            requireCrossingOnCenterline(conductor->shape,
                                        placed((*add)[0][0], (*add)[0][1]), wr,
                                        inPlaneA, inPlaneB,
                                        "turn " + turns[i].get_name() + " outer crossing");
        }
    }
    // Flat-wire-on-round-bore sagitta: a rectangular wire's FLAT inner face can't sit flush against
    // the round bore the way a round wire's tangent does -- placed tangent at its centre, its
    // corners dip into the core by ~(height/2)^2 / (2*boreRadius). That is a real, tiny (< 1e-3
    // mm^3) geometric artifact of flat wire on a curved bore, not an interference to fix, so the
    // core<->conductor tolerance here is looser than the round-wire toroid's exact-tangency 1e-12.
    requireNoPairwiseOverlap(named, 1e-9);
    }
}

// Count the solids in a named conductor body.
static int conductorSolidCount(const TopoDS_Shape& shape) {
    int n = 0;
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) ++n;
    return n;
}

TEST_CASE("Real winding: LITZ wire builds ONE continuous body", "[realwinding]") {
    // Litz flows through the round path as a bare bundle. Round column -> ONE single solid.
    auto magneticJson = loadFixture("realwinding_litz_round.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(shapeVolume(conductor->shape) > 0.0);
    REQUIRE(connectedSolidComponents(conductor->shape) == 1);    // FEM-ready CONNECTED conductor
    REQUIRE(BRepCheck_Analyzer(conductor->shape).IsValid());
}

TEST_CASE("Real winding: round-column RECTANGULAR wire is ONE body", "[realwinding]") {
    // 03_buck (PQ32, 3x0.5 mm rectangular wire) -> fixed-binormal single solid.
    auto magneticJson = loadFixture("realwinding_rect_wire_round.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(shapeVolume(conductor->shape) > 0.0);
    // ONE CONNECTED conductor: the per-run rect compound welds/touches into one component
    // (the whole-spine single body creased at the wrap->lead fillet; see 03_buck).
    REQUIRE(connectedSolidComponents(conductor->shape) == 1);
    REQUIRE(BRepCheck_Analyzer(conductor->shape).IsValid());
}

TEST_CASE("Real winding: rectangular-column RECTANGULAR wire builds", "[realwinding]") {
    // 18_stacked (E70, 5x1 mm): the per-turn racetrack solids. Its copper turns TOUCH, so the fuse
    // would short them into a brick -- the per-turn compound is kept (correct), so it is multi-solid
    // BY DESIGN. Assert only that it builds valid positive-volume copper for every turn.
    auto magneticJson = loadFixture("realwinding_rect_wire_rect.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(shapeVolume(conductor->shape) > 0.0);
    REQUIRE(conductorSolidCount(conductor->shape) >= 1);
}

// ABT #1271. THE LEAD CORNER MUST CLOSE ONTO THE WRAP STRAIGHT IT MEETS.
// The entrance/exit lead corner arc is built FLAT (constant axial coordinate) while the wrap
// straight it joins carries the helical advance (dy/ds = 6.394e-3 on this design, 0.37 deg).
// Each piece's cap is square to its OWN axis, so unless the coplanar-cap shear reaches this
// junction the two cap planes cross at mid-thickness and leave a re-entrant wedge: the junction's
// single 1.000 mm THICKNESS edge comes out as TWO 0.500 mm halves (twelve of them over the two
// leads), the section never closes into a 4-curve W-T-W-T ring, and OMFEM's mapped-hex path
// correctly refuses the design (8 incomplete rings, 2 blocks of 12 faces instead of 6).
// The wedge is ~0.0043 mm^3 against a solid of tens of thousands, so a volume or
// builds-without-throwing check CANNOT see it and passes with the shear reverted -- which is why
// this asserts the junction CLOSES: no straight edge of the conductor measures half the wire
// thickness. Fixture: realwinding_rect_wire_rect.json is 18_stacked itself (E 70/33/32 x2,
// 30 turns of 5 x 1 mm rectangular wire), the design the defect was measured on.
TEST_CASE("Real winding: the lead corner closes onto the wrap straight (no half-thickness edges)",
          "[realwinding][abt1271]") {
    auto magneticJson = loadFixture("realwinding_rect_wire_rect.json");
    const double wireThickness = magneticJson.at("coil")
                                     .at("functionalDescription")[0]
                                     .at("wire")
                                     .at("conductingHeight")
                                     .at("nominal")
                                     .get<double>();
    REQUIRE(wireThickness == Catch::Approx(0.001));   // else the 0.5 mm signature below is wrong
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    // paintCoating=false: the CONDUCTING footprint, which is what the FEM product is drawn at
    // (ABT #1261) and what the half-thickness signature below is stated in. Drawn at the coated
    // envelope the same defect measures half of the OUTER height instead, and a 0.5 mm test would
    // look green over it.
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, /*paintCoating=*/false,
                                       false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(shapeVolume(conductor->shape) > 0.0);

    TopTools_IndexedMapOfShape edges;
    TopExp::MapShapes(conductor->shape, TopAbs_EDGE, edges);
    REQUIRE(edges.Extent() > 0);
    std::ostringstream halves;
    int halfCount = 0, fullCount = 0;
    for (int i = 1; i <= edges.Extent(); ++i) {
        const TopoDS_Edge& e = TopoDS::Edge(edges(i));
        if (BRep_Tool::Degenerated(e)) continue;
        BRepAdaptor_Curve curve(e);
        if (curve.GetType() != GeomAbs_Line) continue;   // only the straight section edges
        GProp_GProps lp;
        BRepGProp::LinearProperties(e, lp);
        // Whole thickness edges must dominate -- if the conductor were drawn at some other
        // thickness (the coated envelope, say) the half-thickness signature below would be
        // measured against the wrong number and the test would pass over the defect.
        if (std::abs(lp.Mass() - wireThickness) < 2e-6) ++fullCount;
        // Half the THICKNESS, to 2 um -- the split-thickness edge the crossing caps produce.
        if (std::abs(lp.Mass() - wireThickness / 2.0) > 2e-6) continue;
        ++halfCount;
        const gp_Pnt c = lp.CentreOfMass();
        halves << "\n  " << lp.Mass() * 1e3 << " mm edge at (" << c.X() * 1e3 << ","
               << c.Y() * 1e3 << "," << c.Z() * 1e3 << ") mm";
    }
    INFO("conductor drawn with " << fullCount << " full-thickness ("
                                 << wireThickness * 1e3 << " mm) straight edges");
    REQUIRE(fullCount > 100);   // 30 turns x 4 corners: the section thickness edges are there
    INFO("half-thickness (" << wireThickness / 2.0 * 1e3 << " mm) straight edges: " << halfCount
                            << halves.str());
    REQUIRE(halfCount == 0);
}

TEST_CASE("Real winding: MULTI-LAYER spread 3-winding toroidal CMC builds clean",
          "[realwinding]") {
    // A spread 3-winding toroidal CMC with TWO rings per 120-degree section
    // (realwinding_cmc_3w_2layer: 18 turns of 1.4 mm wire per winding). Multilayer
    // toroidal builds once MKF's non-physical outer crossings are corrected on the
    // builder side (MKF ABT #231): each outer ring's outer crossing is re-placed at the
    // physical radial stack (ring 0 outer + ring*OD) AND at the inner crossing's azimuth
    // — MKF staggers outer-ring outer angles out of sequence, which would cross
    // consecutive turns' top chords (and the gate exempts consecutive wraps, so it slips
    // through). Ring returns are depth-staggered under the core. The three windings must
    // be fully independent bodies with two layers each.
    auto magneticJson = loadFixture("realwinding_cmc_3w_2layer.json");
    // Fixture sanity via the CLASSIC path (the real-winding enrichment refuses this layout,
    // asserted below): 2 rings per section, else the test is vacuous.
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, false);
    {
        auto layers = enriched.get_coil().get_layers_description();
        REQUIRE(layers.has_value());
        int primaryConductionLayers = 0;
        for (const auto& L : *layers)
            if (L.get_type() == MAS::ElectricalType::CONDUCTION &&
                !L.get_partial_windings().empty() &&
                L.get_partial_windings()[0].get_winding() == "Primary")
                ++primaryConductionLayers;
        REQUIRE(primaryConductionLayers == 2);
    }

    // TWO-LAYER sections bury the inner layer's entrance: ring 2 nests one wire OD inside
    // ring 1 (crossings 1.47 mm apart at a 1.4 mm envelope -- the TURNS are legal), but every
    // straight-dressed lead path from ring 1's start crossing is then blocked by ring 2's
    // returns. A real winder lays the lead FIRST and winds layer 2 around it -- rigid MAS turn
    // geometry cannot express that (lead-space reservation is the layout's job, MKF ABT #187),
    // so the routed-lead builder must REFUSE loudly rather than emit overlapping copper.
    // FULL 2-layer build, end-to-end validated. The chain that makes it possible (all in MKF's
    // layout, per Alf's ruling that the final 3D must come out non-crossing from the winder):
    // (1) the INPUT-CONNECTION ANGULAR CORRIDOR -- rings after the first surrender the
    // connection parallels' angular slots at the section-start edge (span shrunk + shifted),
    // so no later ring places a station behind the connection; (2) per-ring capacity measured
    // at the ring's own radius (kills the ABT #563 overhang); (3) the outer-crossing sweep
    // discards candidates whose implied 3D runs would cross a connection vertical. Measured
    // layout here: ring 0 holds 10 stations from the entrance, ring 1 starts one corridor
    // later -- and the leads route with the classic 90-degree drop.
    enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    // CMC spread windings stay on the fast drawing compound (femReady=false): the test only asserts
    // the three windings don't overlap EACH OTHER, which the per-run compound already satisfies.
    auto named = builder.buildAllNamed(enriched, false, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/false);
    const char* names[] = {"Primary parallel 0", "Secondary parallel 0",
                           "Tertiary parallel 0"};
    std::vector<const mvb::NamedShape*> conductors;
    for (const char* n : names) conductors.push_back(findConductor(named, n));
    for (size_t i = 0; i < conductors.size(); ++i)
        for (size_t j = i + 1; j < conductors.size(); ++j) {
            double v = commonVolume(conductors[i]->shape, conductors[j]->shape);
            INFO("winding-winding overlap '" << names[i] << "' vs '" << names[j]
                                             << "' = " << v);
            REQUIRE(v <= 1e-12);
        }
    requireNoPairwiseOverlap(named, 1e-12);
}

TEST_CASE("Real winding: single-layer spread 3-winding toroidal CMC builds clean",
          "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_cmc_3w_1layer.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);

    mvb::MagneticBuilder builder;
    // CMC spread windings stay on the fast drawing compound (femReady=false): the test only asserts
    // the three windings don't overlap EACH OTHER, which the per-run compound already satisfies.
    auto named = builder.buildAllNamed(enriched, false, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/false);

    // One continuous conductor per winding (three windings, each spread over its own
    // ~120-degree arc), all fully independent copper bodies.
    const char* names[] = {"Primary parallel 0", "Secondary parallel 0",
                           "Tertiary parallel 0"};
    std::vector<const mvb::NamedShape*> conductors;
    for (const char* n : names) conductors.push_back(findConductor(named, n));
    for (size_t i = 0; i < conductors.size(); ++i)
        for (size_t j = i + 1; j < conductors.size(); ++j) {
            double v = commonVolume(conductors[i]->shape, conductors[j]->shape);
            INFO("winding-winding overlap '" << names[i] << "' vs '" << names[j]
                                             << "' = " << v);
            REQUIRE(v <= 1e-12);
        }
    requireNoPairwiseOverlap(named, 1e-12);
}

TEST_CASE("Real winding: pre-enriched input with the flag on throws", "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_round_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, false);

    // Round-trip to a MAS::Magnetic that carries the geometricalDescription.
    json enrichedJson;
    to_json(enrichedJson, enriched);
    MAS::Magnetic masMagnetic = enrichedJson.get<MAS::Magnetic>();
    REQUIRE(masMagnetic.get_core().value().get_geometrical_description().has_value());

    mvb::MagneticBuilder builder;
    REQUIRE_THROWS_WITH(
        builder.buildAllNamed(masMagnetic, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                              mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                              /*useRealWindingGeometry=*/true),
        Catch::Matchers::ContainsSubstring("geometricalDescription"));
}

// ABT #646 — a design that could not be routed until the layout was re-wound.
//
// A real user design: E 16/8/8, N87-class 95, one winding, 12 turns of Synthesized Litz
// 31x0.100mm, 2 conduction layers. Real winding used to refuse it, because the Primary
// entrance lead and an inter-layer dragback were COINCIDENT — centreline distance 8.3e-19 m,
// i.e. exactly zero, against 0.000802 m required.
//
// The cause was upstream and not in the router: MKF's magnetic_autocomplete only wound a coil
// that had no turnsDescription, so a design arriving already wound kept a layout laid out
// WITHOUT the connection corridors. Layer 1 spanned the full 10.2 mm window with a turn sitting
// exactly on the input connection's reserved rectangle, and the router was handed turn positions
// that had never reserved a slot for the lead it then had to route. MKF 90594876 re-winds when
// real winding is asked for; layer 1 comes back one wire slot shorter at the bottom and the
// conductor routes.
//
// This asserts the ROUTE, not the refusal: a design that regresses to a collision here is a
// layout regression, and the exception text will say exactly which two runs met.
TEST_CASE("Real winding: the ABT #646 litz design routes once the layout is re-wound",
          "[realwinding][abt646]") {
    auto magneticJson = loadFixture("realwinding_e16_litz_2layer_leadcollision.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;

    // Idealised loops were never affected: the defect was in the routed conductor.
    REQUIRE_NOTHROW(builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                          mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true));

    std::vector<mvb::NamedShape> named;
    REQUIRE_NOTHROW(named = builder.buildAllNamed(enriched, true, 0, /*wireSeg=*/0,
                                                  mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false,
                                                  false, 0.0, /*useRealWindingGeometry=*/true));
    REQUIRE_FALSE(named.empty());

    // The re-wound layer must actually be shorter than the window — that is the fix, and a
    // full-height layer would mean the corridors were not reserved even if nothing collided.
    auto enrichedCoil = enriched.get_coil();
    auto layers = enrichedCoil.get_layers_description().value();
    double windowHeight = enrichedCoil.resolve_bobbin()
                              .get_processed_description().value()
                              .get_winding_windows()[0].get_height().value();
    bool anyShortened = false;
    for (const auto& layer : layers) {
        if (layer.get_type() != MAS::ElectricalType::CONDUCTION) continue;
        if (layer.get_dimensions()[1] < windowHeight - 1e-9) anyShortened = true;
    }
    INFO("window height " << windowHeight);
    CHECK(anyShortened);
}

// The routed conductors, exported so the winding can be inspected in CAD. No collision gate is
// skipped here any more: this design routes, so the file is real geometry.
TEST_CASE("Real winding: export the ABT #646 design",
          "[realwinding][abt646][diagnostic]") {
    auto magneticJson = loadFixture("realwinding_e16_litz_2layer_leadcollision.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(enriched, /*includeBobbin=*/true, 0, /*wireSeg=*/0,
                                     mvb::DEFAULT_CORE_POLYGON_SEGMENTS, /*paintCoating=*/true,
                                     false, false, 0.0, /*useRealWindingGeometry=*/true);
    REQUIRE_FALSE(all.empty());
    const std::string path = outputPath("abt646_e16_litz_realwinding.step");
    REQUIRE(mvb::exportSTEP(all, path));
    WARN("ABT #646 real-winding geometry written to " << path);
}

TEST_CASE("Real winding: FEM dense toroid is a conformal (non-overlapping) mitre compound",
          "[realwinding]") {
    // A 60-turn toroid is too dense for the single-solid MakePipe sweep to close on its packed hole
    // spine, so femReady=true builds the CONFORMAL mitre-jointed compound instead: each primitive is
    // its own round solid, and neighbours are sliced on their shared angle-bisector plane so they
    // ABUT on a coincident elliptical face rather than interpenetrating. This is the FEM-meshable
    // (no double material) form of a winding that cannot be a single solid.
    auto magneticJson = loadFixture("realwinding_toroid_3in.json");
    // At the fixture's full 60 turns the hole packs TWO layers. Historically this was a hard
    // refusal (no lead corridor anywhere); with MKF's input-connection angular corridor the
    // second ring now starts one corridor past the entrance, the layout re-flows, and the FULL
    // 60-turn build goes through validated -- assert exactly that.
    {
        auto denseEnriched =
            mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
        mvb::MagneticBuilder denseBuilder;
        auto denseNamed =
            denseBuilder.buildAllNamed(denseEnriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       0, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
        const auto* denseConductor = findConductor(denseNamed, "Primary parallel 0");
        REQUIRE(denseConductor != nullptr);
        requireConformalConductor(denseConductor->shape);
    }
    // The conformal mitre-compound structure is exercised at the densest ROUTABLE packing of
    // the same core/wire: 30 turns = a full single layer (rim gap 0.3 mm).
    magneticJson["coil"]["functionalDescription"][0]["numberTurns"] = 30;
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(conductor != nullptr);

    int nSolids = solidCount(conductor->shape);
    // WAS `nSolids > 1`: MakePipe cannot close this dense toroid, so the assembler left it a
    // multi-solid conformal compound and the test pinned that. ABT #1265 then closed the
    // assembly with a glued fuse ("an object per wire to simulate in FEM"), and this conductor
    // now comes out as ONE body -- 132 pieces into 1 on the CMC, 69 on the buck. The old
    // assertion was pinning the emission strategy; what the test is FOR is that the conductor is
    // not a heap of disconnected per-turn loops (see solidCount's comment above), and one body
    // satisfies that completely. So: one body when the close succeeds, and when it does not, the
    // compound must still be conformal -- which is what the next line checks either way.
    REQUIRE(nSolids >= 1);
    // Connected + per-solid valid + neighbours abut WITHOUT interpenetration. The overlap probe
    // is the junction-grid classifier, NOT BRepAlgoAPI_Common: Common on abutting BSpline pipe
    // pairs of this toroid ground >10 min/pair in 2d-extrema root-finding or returned !IsDone
    // (both measured here) -- the exact OCC-boolean pathology the conformal build avoids.
    if (nSolids > 1) requireConformalConductor(conductor->shape);
}

TEST_CASE("Real winding: export the 8t x 2p multi-parallel design",
          "[realwinding][abt685][diagnostic]") {
    // ABT #685 before/after comparison (Alf: "can you show me one STEP from before and one from
    // after"). MVB_STEP_OUT names the file so the same binary can write both sides of a
    // rebuild; MVB_LEAD_NO_VALIDATE skips the collision gate when the point IS to look at a
    // collision.
    auto magneticJson = loadFixture("realwinding_round_2p.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(enriched, /*includeBobbin=*/true, 0, /*wireSeg=*/0,
                                     mvb::DEFAULT_CORE_POLYGON_SEGMENTS, /*paintCoating=*/true,
                                     false, false, 0.0, /*useRealWindingGeometry=*/true);
    REQUIRE_FALSE(all.empty());
    const char* out = std::getenv("MVB_STEP_OUT");
    const std::string path = out ? out : outputPath("abt685_multiparallel.step");
    REQUIRE(mvb::exportSTEP(all, path));
    WARN("8t x 2p real-winding geometry written to " << path);
}

TEST_CASE("TMP failing-case export and overlap location", "[tmpfail]") {
    auto magneticJson = loadFixture("realwinding_round_2p.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);
    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                                       true);
    const mvb::NamedShape* p0 = nullptr;
    const mvb::NamedShape* p1 = nullptr;
    for (const auto& ns : named) {
        if (ns.name == "Primary parallel 0") p0 = &ns;
        if (ns.name == "Primary parallel 1") p1 = &ns;
    }
    REQUIRE(p0);
    REQUIRE(p1);
    // Solids in EXPORT ORDER — the order a STEP viewer numbers them ("Primary parallel 0",
    // then 001, 002, ...). Prim indices are not it: degenerate slivers are pruned first.
    auto solidsOf = [](const TopoDS_Shape& sh) {
        std::vector<TopoDS_Shape> out;
        for (TopExp_Explorer e(sh, TopAbs_SOLID); e.More(); e.Next()) out.push_back(e.Current());
        return out;
    };
    const auto s0 = solidsOf(p0->shape);
    const auto s1 = solidsOf(p1->shape);
    auto nameOf = [](const char* base, size_t k) {
        char buf[64];
        if (k == 0) std::snprintf(buf, sizeof buf, "%s", base);
        else std::snprintf(buf, sizeof buf, "%s%02zu", base, k);
        return std::string(buf);
    };
    std::fprintf(stderr, "[solids] Primary parallel 0: %zu solids, Primary parallel 1: %zu\n",
                 s0.size(), s1.size());
    for (size_t a = 0; a < s0.size(); ++a) {
        for (size_t b = 0; b < s1.size(); ++b) {
            Bnd_Box ba, bb;
            BRepBndLib::Add(s0[a], ba);
            BRepBndLib::Add(s1[b], bb);
            if (ba.IsOut(bb)) continue;
            BRepAlgoAPI_Common common(s0[a], s1[b]);
            if (!common.IsDone()) continue;
            const double v = shapeVolume(common.Shape());
            if (v <= 1e-15) continue;
            std::fprintf(stderr, "[hit] %-24s  vs  %-24s   common = %.4e mm3\n",
                         nameOf("Primary parallel 0", a).c_str(),
                         nameOf("Primary parallel 1", b).c_str(), v * 1e9);
        }
    }
}

TEST_CASE("TMP fixture probe", "[tmpfix]") {
    const char* name = std::getenv("MVB_FIXTURE");
    REQUIRE(name != nullptr);
    std::ifstream f(std::string(MAS_EXAMPLES_DIR) + "/" + name);
    REQUIRE(f.good());
    json j = json::parse(f);
    auto magneticJson = j.contains("magnetic") ? j.at("magnetic") : j;
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);
    mvb::MagneticBuilder builder;
    try {
        auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                           mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false,
                                           0.0, true);
        std::fprintf(stderr, "[fix] %s BUILT ok\n", name);
        const char* out = std::getenv("MVB_STEP_OUT");
        if (out) REQUIRE(mvb::exportSTEP(named, out));
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[fix] %s THROW: %s\n", name, e.what());
    }
}

// ABT #353 -- NO ROUTED COPPER IN THE FLANGE. The WE-TI 7447720470 drum is the design the ticket
// was filed on: its winding groove runs y = [-3.2488, +3.2488] mm and its flanges fill the bands
// beyond that (the bottom one out to r = 3.898 mm, the top one to r = 2.498 mm), so anything a
// router places "just below the winding" is inside ferrite. The Z end-run planner did exactly
// that -- its duck plane came from the copper ring inventory alone and sat 0.23 mm below the
// groove floor, with the wire 0.44 mm into the flange -- and it now classifies every candidate
// lane against the core solids instead (planZEndRuns). This is the design-level pin: what MVB++
// routes on this drum stays in the groove.
//
// SCOPE (Alf, 2026-08-23): copper touching the core is a filling-factor consequence of what the
// user asked for, not something the builder vetoes -- there is deliberately no general
// copper-vs-core refusal in the build. What is pinned here is MVB++'s own PLACEMENT on a design
// whose turns fit their window with room to spare; a failure means a router put copper in the
// ferrite, and is not to be relaxed away.
//
// Deliberately on the PATH (centreline) API: it exercises the same routing as the solid build
// without depending on wire assembly, so a defect in the assembler cannot mask -- or fake -- a
// routing regression.
TEST_CASE("Real winding: the drum's connections stay out of the core flange (ABT #353)",
          "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_drum_we_ti.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    // Hands the core solids down to ConductorBuilder as routing obstacles.
    auto paths = builder.buildRealWindingPaths(enriched);
    REQUIRE(!paths.empty());
    auto coreShapes = builder.buildCoreNamed(enriched.get_core());
    REQUIRE(!coreShapes.empty());

    // The obstacle, measured from the CORE SOLID itself rather than from any dimension in the
    // fixture: the flange radius is the core's own bounding radius, and the groove floor/ceiling
    // are found by walking the axis of the winding band until the classifier reports material.
    // Nothing here can be defused by editing a number in the JSON.
    double flangeR = 0.0;
    for (const auto& ns : coreShapes) {
        Bnd_Box bx;
        BRepBndLib::Add(ns.shape, bx);
        double x0, y0, z0, x1, y1, z1;
        bx.Get(x0, y0, z0, x1, y1, z1);
        flangeR = std::max({flangeR, std::abs(x0), std::abs(x1), std::abs(z0), std::abs(z1)});
    }
    REQUIRE(flangeR > 0.0);

    std::vector<std::unique_ptr<BRepClass3d_SolidClassifier>> cls;
    for (const auto& ns : coreShapes)
        for (TopExp_Explorer se(ns.shape, TopAbs_SOLID); se.More(); se.Next())
            cls.push_back(std::make_unique<BRepClass3d_SolidClassifier>(se.Current()));
    REQUIRE(!cls.empty());
    auto material = [&](double x, double y, double z) {
        for (auto& c : cls) {
            c->Perform(gp_Pnt(x, y, z), 1e-9);
            if (c->State() == TopAbs_IN) return true;
        }
        return false;
    };

    // The groove and the two flanges, MEASURED from the solid (this drum's flanges differ:
    // the bottom one is the full 7.8 mm OD, the top one only 5.0 mm). Walk the post's own
    // radius up and down until material starts -- that is the groove face -- then walk radially
    // inside each flange band to find how far that flange reaches.
    const double step = flangeR / 2000.0;
    double postR = 0.0;
    for (double r = step; r < flangeR; r += step)
        if (material(r, 0.0, 0.0)) postR = r;
    REQUIRE(postR > 0.0);
    const double rGroove = postR + 4.0 * step;
    REQUIRE(!material(rGroove, 0.0, 0.0));
    double yFloor = 0.0, yCeil = 0.0;
    for (double y = 0.0; y > -flangeR * 4.0; y -= step)
        if (material(rGroove, y, 0.0)) { yFloor = y + step; break; }
    for (double y = 0.0; y < flangeR * 4.0; y += step)
        if (material(rGroove, y, 0.0)) { yCeil = y - step; break; }
    REQUIRE(yFloor < 0.0);
    REQUIRE(yCeil > 0.0);
    double rBottomFlange = 0.0, rTopFlange = 0.0;
    for (double r = rGroove; r < flangeR + step; r += step) {
        if (material(r, yFloor - 4.0 * step, 0.0)) rBottomFlange = r;
        if (material(r, yCeil + 4.0 * step, 0.0)) rTopFlange = r;
    }
    REQUIRE(rBottomFlange > postR);
    REQUIRE(rTopFlange > postR);
    std::fprintf(stderr,
                 "[ABT353] drum measured from the solid: post r=%.4f mm, groove y=[%.4f, %.4f] "
                 "mm, bottom flange to r=%.4f mm, top flange to r=%.4f mm\n",
                 postR * 1e3, yFloor * 1e3, yCeil * 1e3, rBottomFlange * 1e3, rTopFlange * 1e3);

    int inFlange = 0, inCore = 0;
    double deepest = 0.0;   // worst excursion past a groove face, over the flange (metres)
    for (const auto& p : paths) {
        for (const auto& prim : p.prims) {
            for (const auto& q : prim) {
                const double r = std::hypot(q[0], q[2]);
                const double ux = r > 1e-12 ? q[0] / r : 1.0;
                const double uz = r > 1e-12 ? q[2] / r : 0.0;
                // (a) the ticket's own measurement: the wire TUBE reaching past the groove floor
                //     or ceiling while still under that side's flange is copper inside ferrite.
                {
                    const double over =
                        std::max(r <= rBottomFlange ? yFloor - (q[1] - p.wireRadius) : -1.0,
                                 r <= rTopFlange ? (q[1] + p.wireRadius) - yCeil : -1.0);
                    if (over > 1e-9) {
                        INFO("copper at r=" << r * 1e3 << " mm y=" << q[1] * 1e3
                                            << " mm reaches " << over * 1e3
                                            << " mm past the groove face (floor " << yFloor * 1e3
                                            << ", ceiling " << yCeil * 1e3 << " mm, flanges r="
                                            << rBottomFlange * 1e3 << "/" << rTopFlange * 1e3
                                            << " mm) in '" << p.name << "'");
                        ++inFlange;
                        deepest = std::max(deepest, over);
                    }
                }
                // (b) the general statement: no part of the wire tube is inside a core solid.
                const std::array<std::array<double, 3>, 5> probes{{
                    {q[0], q[1], q[2]},
                    {q[0], q[1] - p.wireRadius, q[2]},
                    {q[0], q[1] + p.wireRadius, q[2]},
                    {q[0] - ux * p.wireRadius, q[1], q[2] - uz * p.wireRadius},
                    {q[0] + ux * p.wireRadius, q[1], q[2] + uz * p.wireRadius}}};
                for (const auto& pr : probes) {
                    if (!material(pr[0], pr[1], pr[2])) continue;
                    INFO("copper at (" << q[0] * 1e3 << "," << q[1] * 1e3 << "," << q[2] * 1e3
                                       << ") mm has tube material inside the core, in '" << p.name
                                       << "'");
                    ++inCore;
                    break;
                }
            }
        }
    }
    std::fprintf(stderr, "[ABT353] deepest excursion into a flange band: %.4f mm\n",
                 deepest * 1e3);
    CHECK(inFlange == 0);
    CHECK(inCore == 0);
}

// A SOLDERED WIRE TERMINAL NEVER TOUCHES THE FOIL (the mesher's veto on an exact tangency).
//
// The terminal was drawn tangent to its sheet -- bare copper on bare copper, no standoff -- so
// that the bump the next turn rides is exactly one OD and the wire cannot reach the turn above
// (Alf, 2026-09-04). The CAD gate liked it: WATERTIGHT, NO OVERLAPS, and the enamel rule
// CERTIFIED at 0 nm. It is nevertheless unmeshable, and no element size fixes it: a cylinder
// laid on a plane meets it at a ZERO dihedral angle, and every way of decomposing that corner --
// two solids in contact, one fused body, fillets that stop short of the tangent line or run all
// the way to it -- still contains a feature of zero thickness. tetgen said exactly that on
// two_switch_forward:
//     PLC Error: A segment and a facet intersect at point (-0.267, 5.183, 15.9845) mm
// which is the bottom generatrix of 'Secondary parallel 0 exit lead' on the face of 'Secondary
// parallel 0' -- the tangent line itself.
//
// The physics agrees with the mesher: a soldered wire floats on the joint, it never touches bare
// copper. So the wire is lifted by a solder film and the solder fills its whole footprint. Alf's
// reason for asking for tangency is kept exactly -- the LANE the stack rides over the terminals
// grew by the film too, so the wire still cannot reach the turn above.
//
// This test pins the invariant that the tangency must not come back: a strictly positive gap
// between the terminal and its sheet, of the film's size, BRIDGED by the solder. Revert the film
// to zero and the first REQUIRE fails (the gap goes to 0).
TEST_CASE("Real winding: a foil terminal floats on its solder film, never tangent to the sheet",
          "[realwinding][foil]") {
    auto magneticJson = loadFixture("two_switch_forward_transformer_complete.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    // --copper: the terminal is drawn on its CONDUCTING diameter, which is the surface that is
    // soldered. The outer diameter is what the bump RESERVES, not what is drawn.
    auto named = builder.buildAllNamed(enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                       /*wirePolygonSegments=*/0,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                       /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);

    auto find = [&](const std::string& n) {
        for (const auto& ns : named)
            if (ns.name == n) return ns.shape;
        FAIL("no solid named '" << n << "' in the build");
        return TopoDS_Shape{};
    };
    auto gap = [](const TopoDS_Shape& a, const TopoDS_Shape& b) {
        BRepExtrema_DistShapeShape d(a, b);
        REQUIRE(d.IsDone());
        return d.Value();
    };
    // A failure here is a CONTACT, and the only useful thing to say about a contact is where it
    // is: the point names the feature that is touching.
    auto whereClosest = [](const TopoDS_Shape& a, const TopoDS_Shape& b) {
        BRepExtrema_DistShapeShape d(a, b);
        std::ostringstream o;
        o.precision(6);
        if (d.IsDone() && d.NbSolution() > 0) {
            const gp_Pnt p1 = d.PointOnShape1(1), p2 = d.PointOnShape2(1);
            o << "(" << p1.X() * 1e3 << ", " << p1.Y() * 1e3 << ", " << p1.Z() * 1e3 << ") mm"
              << " <-> (" << p2.X() * 1e3 << ", " << p2.Y() * 1e3 << ", " << p2.Z() * 1e3 << ") mm";
        }
        return o.str();
    };

    int joints = 0;
    for (const std::string& which : {std::string("entrance"), std::string("exit")}) {
        const std::string lead = "Secondary parallel 0 " + which + " lead";
        const TopoDS_Shape sheet = find("Secondary parallel 0");
        const TopoDS_Shape joint = find(lead + " solder joint");
        // The lead is an L: the SOLDERED RUN down the sheet's height, and a radial piece that
        // elbows out through the margin to the tip. Only the soldered run is measured here. The
        // radial piece grazes the sheet's top edge on its way out and that is not a defect --
        // they are the same conductor, and the meshing fragment welds a same-region contact into
        // one volume. What killed the mesher was the SOLDERED RUN lying on the face.
        TopoDS_Shape wire;
        {
            double best = -1.0;
            for (TopExp_Explorer e(find(lead), TopAbs_SOLID); e.More(); e.Next()) {
                Bnd_Box b;
                BRepBndLib::Add(e.Current(), b);
                double x0, y0, z0, x1, y1, z1;
                b.Get(x0, y0, z0, x1, y1, z1);
                if (y1 - y0 > best) { best = y1 - y0; wire = e.Current(); }
            }
            REQUIRE(best > 0.0);
        }

        // The wire's drawn diameter, measured, so the check is independent of the unit the
        // builder works in: the soldered run is a cylinder, its X extent is one diameter.
        Bnd_Box wb;
        BRepBndLib::Add(wire, wb);
        double xlo, ylo, zlo, xhi, yhi, zhi;
        wb.Get(xlo, ylo, zlo, xhi, yhi, zhi);
        const double diameter = xhi - xlo;
        REQUIRE(diameter > 0.0);

        const double film = gap(sheet, wire);
        INFO(lead << ": film " << film << ", diameter " << diameter
                  << ", ratio " << film / diameter
                  << ", closest " << whereClosest(sheet, wire));
        // NOT TANGENT. This is the assertion the PLC error asked for; at tangency it is 0.
        REQUIRE(film > 0.1 * diameter * 0.5);
        // ...and it is the film, not an accidental drift: 0.05 mm under a 0.501 mm wire.
        REQUIRE(film / diameter == Catch::Approx(0.05 / 0.501).epsilon(0.05));
        // BRIDGED. A gap that nothing crosses would be an open circuit, not a joint.
        REQUIRE(gap(joint, sheet) < 1e-9 * diameter + 1e-12);
        REQUIRE(gap(joint, wire) < 1e-9 * diameter + 1e-12);
        ++joints;
    }
    REQUIRE(joints == 2);
}

// ABT #1111 (2026-09-06). The FEM product's periodic-face -> B-spline pass used to run in
// buildAllNamed, in METRES; exportSTEP then scaled the poles x1000 and left the knots alone, so
// every length-parametrised B-spline direction (a facet strip's generatrix, a planar cap) spanned
// 1000 mm of geometry per parameter unit. BOPAlgo converts 3D tolerances to parameter space
// through exactly that ratio and reported sporadic self-intersections on faceted revolves: 11 of
// 38 corpus designs at --segments 12, each one clean again once rescaled to metres. The pass now
// runs inside the exporter on the millimetre shape. This reads the written file back and checks
// the invariant directly: no B-spline face may span more than a few millimetres of geometry per
// parameter unit (a well-parametrised conversion spans ~1 mm/unit along a length, ~r mm/rad
// along an angle), and every solid must pass BOPAlgo's self-intersection analysis.
TEST_CASE("Real winding: the FEM STEP's B-splines are parametrised in the millimetres they are "
          "written in",
          "[realwinding][step][nurbs]") {
    auto magneticJson = loadFixture("realwinding_e138_rectcolumn.json");  // E 13/7/4, 6 turns, 0.4 mm round
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                       /*wirePolygonSegments=*/12,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                       /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    REQUIRE_FALSE(named.empty());

    const std::string path = outputPath("abt1111_e138_seg12_fem.step");
    mvb::StepExportOptions opts;
    opts.nurbsPeriodicSolids = true;
    REQUIRE(mvb::exportSTEP(named, path, opts));

    const auto back = mvb::importSTEP(path);
    REQUIRE_FALSE(back.empty());
    int bsplineFaces = 0, periodicFaces = 0, solids = 0, faulty = 0;
    double worstMmPerUnit = 0.0;
    std::string worstFace;
    for (const auto& ns : back) {
        for (TopExp_Explorer se(ns.shape, TopAbs_SOLID); se.More(); se.Next(), ++solids) {
            int fidx = 0;
            for (TopExp_Explorer fe(se.Current(), TopAbs_FACE); fe.More(); fe.Next(), ++fidx) {
                const TopoDS_Face f = TopoDS::Face(fe.Current());
                Handle(Geom_Surface) srf = BRep_Tool::Surface(f);
                if (!srf.IsNull() && (srf->IsUPeriodic() || srf->IsVPeriodic())) ++periodicFaces;
                BRepAdaptor_Surface ad(f, false);
                if (ad.GetType() != GeomAbs_BSplineSurface) continue;
                ++bsplineFaces;
                // Resolution(1 mm) is the parameter step that moves the surface by 1 mm, i.e.
                // 1 / (mm per parameter unit). Metre-scale knots under mm poles give 1e-3.
                for (double res : {ad.UResolution(1.0), ad.VResolution(1.0)}) {
                    if (!(res > 0.0)) continue;
                    const double mmPerUnit = 1.0 / res;
                    if (mmPerUnit > worstMmPerUnit) {
                        worstMmPerUnit = mmPerUnit;
                        worstFace = ns.name + " solid " + std::to_string(solids) + " face " +
                                    std::to_string(fidx);
                    }
                }
            }
            BOPAlgo_ArgumentAnalyzer an;
            an.SetShape1(se.Current());
            an.SelfInterMode() = Standard_True;
            an.ArgumentTypeMode() = Standard_False;
            an.Perform();
            if (an.HasFaulty()) ++faulty;
        }
    }
    INFO("solids " << solids << ", B-spline faces " << bsplineFaces << ", periodic faces left "
                   << periodicFaces << ", worst parametrisation " << worstMmPerUnit
                   << " mm per parameter unit at " << worstFace);
    // The design carries faceted revolves at 12 segments, so the pass had something to convert.
    REQUIRE(bsplineFaces > 0);
    CHECK(periodicFaces == 0);
    // Length-parametrised directions span ~1 mm/unit, angles ~r mm/rad (r < 1 mm here). The
    // metre-frame conversion this guards against sits at 1000.
    CHECK(worstMmPerUnit < 20.0);
    CHECK(faulty == 0);
    std::filesystem::remove(path);
}

// ABT #1111 (2026-09-06): in faceted mode an arc bent tighter than kFacetTightArcExactRatio wire
// radii is revolved on the exact round profile (a torus), because the 12-gon's innermost facet
// strip on such a piece is micrometres long and BOPAlgo rejects it. The in-memory FEM product
// keeps its analytic surfaces (the B-spline re-expression happens in the exporter), so the rule
// is visible directly: the winding must contain tori, every one of them tighter than the ratio,
// and each with the wire's own minor radius.
TEST_CASE("Real winding: faceted mode revolves tight arcs exactly, and only tight arcs",
          "[realwinding][tightarc]") {
    // A toroid: its window corners are ARC3 pieces bent at kRoundCornerBendFactor (1.05 r), the
    // tightest arcs the assembler emits. (A rect-column design's corners are pitched SPIRALs and
    // would exercise nothing here.)
    auto magneticJson = loadFixture("realwinding_toroid.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                       /*wirePolygonSegments=*/12,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                       /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    REQUIRE_FALSE(named.empty());

    int tori = 0, loose = 0;
    double minorMin = 1e300, minorMax = 0.0, ratioMax = 0.0;
    for (const auto& ns : named) {
        if (ns.name.find(" parallel ") == std::string::npos || ns.name.find(" terminal ") != std::string::npos)
            continue;
        for (TopExp_Explorer fe(ns.shape, TopAbs_FACE); fe.More(); fe.Next()) {
            BRepAdaptor_Surface ad(TopoDS::Face(fe.Current()), false);
            if (ad.GetType() != GeomAbs_Torus) continue;
            ++tori;
            const gp_Torus t = ad.Torus();
            const double ratio = t.MajorRadius() / t.MinorRadius();
            ratioMax = std::max(ratioMax, ratio);
            minorMin = std::min(minorMin, t.MinorRadius());
            minorMax = std::max(minorMax, t.MinorRadius());
            if (ratio >= mvb::kFacetTightArcExactRatio) ++loose;
        }
    }
    INFO("torus faces " << tori << ", tightest-to-loosest bend ratio max " << ratioMax
                        << ", minor radius " << minorMin * 1e6 << ".." << minorMax * 1e6 << " um");
    // The terminal-corner fillets sit at 1.05 r, so the rule must have fired somewhere...
    CHECK(tori > 0);
    // ...and nowhere else: a faceted design has no other reason to carry a torus.
    CHECK(loose == 0);
    // The revolve is of the wire itself: one radius for every torus, and a wire-sized one.
    CHECK(minorMax - minorMin < 1e-9);
    CHECK(minorMin > 0.05e-3);
    CHECK(minorMax < 1.0e-3);
}

// ---------------------------------------------------------------------------------------
// TOROID TERMINALS DROP TO ONE PLANE (Alf, 2026-09-12; ABT #1159, superseding the #1155 x/z tip
// pinning). ABT #1248 moved these checks to the EXPORTED frame (buildRealWindingPaths now returns
// it, the frame of the STEP) and runs them in BOTH mountings: before #1248 the paths came back in
// the build frame, where the drops ran -Y while the STEP showed them along +Z. The -Y direction
// and the one-plane assertions are unchanged; the "not gratuitously lower" bound now also accepts
// a plane set by the CORE (a standing ring's bare rim can reach below the copper). Every toroidal terminal lead ends with a -Y drop past the rim, and every tip of the
// component lies on ONE XZ plane below the lowest copper surface by max(2 OD of the thinnest
// wire, 1 OD of the thickest) -- like a real toroid dressed for the board. Under the old rule the tips pointed
// radially from their own crossing azimuths, sat on a circle, and common_mode_choke_complete was
// refused outright ("runs nearly parallel to the common terminal plane").
namespace {
struct ToroidPlaneReport { size_t tips; double planeY; double lowestCopperY; double clearance; };

ToroidPlaneReport requireToroidTerminalsOnOnePlane(const std::string& fixture,
                                                   MAS::OrientationEnum mounting) {
    OpenMagnetics::SettingsGuard<MAS::OrientationEnum> mountingGuard(
        OpenMagnetics::Settings::GetInstance(), &OpenMagnetics::Settings::get_toroid_mounting,
        &OpenMagnetics::Settings::set_toroid_mounting, mounting);
    auto magneticJson = loadFixture(fixture);
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    std::vector<mvb::ConductorBuilder::PathPolyline> paths;
    // Includes the collision gate: the now-parallel drops are checked against each other.
    REQUIRE_NOTHROW(paths = builder.buildRealWindingPaths(enriched));
    REQUIRE(!paths.empty());

    struct Tip { std::string owner; std::array<double, 3> p, u; };
    std::vector<Tip> tips;
    double odMax = 0.0, odMin = std::numeric_limits<double>::max();
    for (const auto& path : paths) {
        tips.push_back({path.name + " end0", path.end0, path.dir0});
        tips.push_back({path.name + " end1", path.end1, path.dir1});
        odMax = std::max(odMax, 2.0 * path.wireRadius);
        odMin = std::min(odMin, 2.0 * path.wireRadius);
    }
    const double clearance = std::max(2.0 * odMin, odMax);
    double planeLo = std::numeric_limits<double>::max(), planeHi = std::numeric_limits<double>::lowest();
    for (const auto& t : tips) {
        INFO(fixture << ": " << t.owner << " at (" << t.p[0] * 1e3 << ", " << t.p[1] * 1e3 << ", "
                     << t.p[2] * 1e3 << ") mm dir (" << t.u[0] << ", " << t.u[1] << ", " << t.u[2] << ")");
        // Every terminal points straight down: the port normal is -Y for all of them.
        CHECK(std::fabs(t.u[0]) < 1e-9);
        CHECK(std::fabs(t.u[2]) < 1e-9);
        CHECK(t.u[1] < -(1.0 - 1e-9));
        planeLo = std::min(planeLo, t.p[1]);
        planeHi = std::max(planeHi, t.p[1]);
    }
    INFO(fixture << ": tip plane y spans [" << planeLo * 1e3 << ", " << planeHi * 1e3 << "] mm");
    CHECK(planeHi - planeLo < 1e-9);

    // The plane sits `clearance` below the lowest copper that is NOT a drop. A sampled point
    // whose envelope reaches below (plane + clearance) must belong to a drop: same (x, z) as a tip.
    double lowest = std::numeric_limits<double>::max();
    for (const auto& path : paths) {
        for (const auto& prim : path.prims) {
            for (const auto& q : prim) {
                bool onDrop = false;
                for (const auto& t : tips)
                    if (std::hypot(q[0] - t.p[0], q[2] - t.p[2]) < 1e-6) { onDrop = true; break; }
                if (onDrop) continue;
                lowest = std::min(lowest, q[1] - path.wireRadius);
            }
        }
    }
    // The core in the exported frame: the plane clears it too (ABT #1248).
    Bnd_Box coreBox;
    const auto frame = mvb::MagneticBuilder::toroidMountingFrameOf(enriched);
    for (const auto& ns : builder.buildCoreNamed(enriched.get_core()))
        BRepBndLib::AddOptimal(BRepBuilderAPI_Transform(ns.shape, frame.toExported, true).Shape(),
                               coreBox, false, false);
    REQUIRE_FALSE(coreBox.IsVoid());
    const double coreBottom = coreBox.CornerMin().Y();
    INFO(fixture << (mounting == MAS::OrientationEnum::VERTICAL ? " VERTICAL" : " HORIZONTAL")
                 << ": lowest non-drop copper surface " << lowest * 1e3 << " mm, core bottom "
                 << coreBottom * 1e3 << " mm, plane "
                 << planeLo * 1e3 << " mm, clearance max(2 OD_thin, 1 OD_thick) = "
                 << clearance * 1e3 << " mm");
    CHECK(planeLo <= lowest - clearance + 1e-9);
    CHECK(planeLo <= coreBottom - clearance + 1e-9);
    lowest = std::min(lowest, coreBottom);
    // ... and not gratuitously lower. The lowest non-drop copper is the rim fillet's underside,
    // and the consumer polyline samples that arc coarsely enough to miss its extremum by ~10 % of
    // an OD (measured 0.12 / 0.06 / 0.40 mm on CMC / buck / CT), so the bound carries a
    // quarter-OD (of the thickest wire) of slack; the builder measures the same arcs finer.
    CHECK(planeLo >= lowest - clearance - 0.25 * odMax);
    return {tips.size(), planeLo, lowest, clearance};
}
} // namespace

TEST_CASE("Real winding: CMC toroid terminals all drop to one plane below the part (ABT #1159)",
          "[realwinding][toroidtips]") {
    // Two windings on opposite sides of the ring: under the radial-tip rule no common plane
    // could hold their four tips. Now all four drop in -Y onto one plane.
    for (auto m : {MAS::OrientationEnum::VERTICAL, MAS::OrientationEnum::HORIZONTAL}) {
        const auto r = requireToroidTerminalsOnOnePlane("common_mode_choke_complete.json", m);
        CHECK(r.tips == 4);
    }
}

TEST_CASE("Real winding: the 3-parallel buck toroid's six terminals share one plane (ABT #1155)",
          "[realwinding][toroidtips]") {
    // Six leads at six azimuths, six drops, one plane, and the gate proves the parallel drops
    // clear each other.
    for (auto m : {MAS::OrientationEnum::VERTICAL, MAS::OrientationEnum::HORIZONTAL}) {
        const auto r = requireToroidTerminalsOnOnePlane("buck_inductor_complete.json", m);
        CHECK(r.tips == 6);
    }
}

TEST_CASE("Real winding: single-turn toroid primary drops both terminals without collision",
          "[realwinding][toroidtips]") {
    // A one-turn bore-through conductor: entrance and exit share the crossing azimuth, so the
    // exit's drop from the top face would run straight through the entrance lead below unless
    // the router moves it -- and the gate must SEE that pair (same turn ordinal, different
    // terminals).
    for (auto m : {MAS::OrientationEnum::VERTICAL, MAS::OrientationEnum::HORIZONTAL}) {
        const auto r = requireToroidTerminalsOnOnePlane("current_transformer_complete.json", m);
        CHECK(r.tips >= 4);
    }
}
