// mvbpp_stepcheck <file.step> [options] — is this exported assembly WATERTIGHT and MESHABLE?
// Moved verbatim from OMFEM tools/omfem_stepcheck.cpp at 8ed8a92 (ABT #1588, step 8); renamed; SizeField.hpp -> mvb; StepAudit.hpp -> mvb; omfem::cad:: -> mvb; omfem::mesh:: -> mvb. No logic changed.
//
// The CAD sweep proves MVB++ can *write* a STEP; omfem_step_intersect proves the solids do not
// interpenetrate. Neither answers the question a mesher asks: can OCC hand every solid to gmsh
// as a closed, valid, conformal volume, and does a tet mesh actually come out the other end?
// This tool answers exactly that, on the artifact ON DISK — not on the in-memory shapes the
// builder happened to have — because the file is what a customer (and OMFEM's own
// mesh3d_from_mas STEP round-trip) actually consumes.
//
// STAGE A — B-Rep validity, per named solid (OCCT only, no booleans):
//   * BRepCheck_Analyzer (with BOP-consistency) — the topology/geometry checker.
//   * volume > 0 and finite (a zero/negative-volume solid meshes to nothing or inverts).
//   * WATERTIGHT: every edge of the solid is shared by exactly two faces. One face = a free
//     edge = an open shell = the mesher's "cannot classify volume". More than two = a
//     non-manifold junction, which tets cannot represent either. Degenerated edges (cone/
//     sphere apex) are exempt: they bound no material.
//   * closed shells (BRep_Tool::IsClosed) and shell count (>1 = internal void; legal, reported).
//
// STAGE B — meshability, the real gmsh path (--mesh):
//   import (OCCTargetUnit M, as MasMesher does) -> optional air box -> fragment the whole
//   assembly (BOPAlgo GLUE + non-destructive: the conformal winding abuts on exactly coincident
//   faces) -> 3D tets. Reports volumes in/out, whether fragment PRESERVED total volume, element
//   count, and min SICN quality. A negative SICN is an inverted element: meshed, but unusable.
//
// A file is FEM-ready when stage A is clean, fragment conserves volume, and the mesh has no
// inverted elements. Exit 0 = ready, 1 = not, 2 = usage/read error.
#include <sys/stat.h>   // --contacts-cache keys on the STEP's size+mtime

#include <thread>
#include <atomic>

#include "mvb/mesh/StepAudit.h"   // --audit: the overlap audit, on the SAME import

#include <STEPCAFControl_Reader.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDataStd_Name.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepAdaptor_Surface.hxx>        // surface type of a faulty chunk's faces
#include <GeomAbs_SurfaceType.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>   // self-intersection: BRepCheck cannot see it
#include <BOPAlgo_CheckResult.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>
#include <gmsh.h>
#include "mvb/mesh/SizeField.h"

#include <algorithm>
#include <chrono>
#include <csignal>
#include <unistd.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

struct SolidInfo {
    std::string name;
    double volume_mm3 = 0;      // STEP declares mm; OCCT reads the file's own unit
    bool   analyzerOk = true;
    // BRepCheck does NOT test self-intersection -- only BOPAlgo_ArgumentAnalyzer does. A fuse
    // can be BRepCheck-VALID, watertight, overlap-clean and STILL be self-intersecting, and
    // that solid is a trap: valid in memory, often INVALID after the STEP round trip, and dear
    // or impossible to boolean. Measured 2026-09-01: 02_flyback -- a design this checker calls
    // WATERTIGHT and meshes 588/588 -- ships at least 6 self-intersecting copper solids
    // (Primary parallel 0, solids 7/14/25/34/61/68). Nothing in the sweep looked for it, so it
    // shipped green. It is a defect of the geometry, so it belongs in the CAD gate.
    bool   selfIntersects = false;
    int    selfIntCount = 0;
    // Small (sub-tolerance) edges, the OTHER finding the same analyzer makes. Nothing gated
    // these either, and they are a classic source of boolean instability. REPORTED, NOT GATED:
    // unlike a self-intersection, a small edge is not automatically a defect (a genuinely tiny
    // feature can be legitimate), and turning it into a stage-A failure would change the
    // corpus verdict on Alf's behalf. Surfaced so the decision can be made on evidence --
    // measured 2026-09-01: 08_pq carries 2, which is why stepheal called it FAULTY while its
    // self-intersection count is zero.
    int    smallEdgeCount = 0;
    int    freeEdges = 0;       // shared by exactly one face
    int    nonManifoldEdges = 0;// shared by three or more
    int    shells = 0;
    int    openShells = 0;
    int    faces = 0;
    double cx = 0, cy = 0, cz = 0;   // centre of mass [mm] -- names an unmeshed gmsh volume
    // Bbox extents [mm], longest-to-shortest. For a self-intersecting chunk this is the datum
    // that says WHICH defect it is: a piece whose longest extent is no bigger than the wire
    // diameter is a stub too short to sweep, whereas a long thin one folded on itself means the
    // PATH curved tighter than the wire radius. Both are emitted defects, but they have
    // different owners (chunk construction vs the drawn route).
    double ext0 = 0, ext1 = 0, ext2 = 0;
    // Surface-type histogram. Says whether a faulty chunk is a STRAIGHT prism (planar sides,
    // which cannot fold on its own) or a CURVED sweep (BSpline sides, which folds when the path
    // turns tighter than the section). Different causes, different fixes.
    std::string surfMix;
};

double solidVolume(const TopoDS_Shape& s, double* cx = nullptr, double* cy = nullptr,
                   double* cz = nullptr) {
    GProp_GProps p;
    BRepGProp::VolumeProperties(s, p);
    if (cx) { const gp_Pnt c = p.CentreOfMass(); *cx = c.X(); *cy = c.Y(); *cz = c.Z(); }
    return p.Mass();
}

// Watertightness of ONE solid, boolean-free: an edge of a closed manifold solid is used by
// exactly two face-sides. TopExp counts each USE, so a cylinder's seam (one face, two uses)
// correctly reports 2 and is not flagged.
void edgeManifoldCheck(const TopoDS_Shape& solid, int& freeEdges, int& nonManifold) {
    freeEdges = 0; nonManifold = 0;
    TopTools_IndexedDataMapOfShapeListOfShape map;
    TopExp::MapShapesAndUniqueAncestors(solid, TopAbs_EDGE, TopAbs_FACE, map, /*useOrientation*/false);
    // MapShapesAndUniqueAncestors collapses the seam's two uses into one face, so count uses
    // ourselves: walk the faces and tally every edge occurrence.
    std::map<int, int> uses;   // edge index in map -> number of face uses
    for (TopExp_Explorer fx(solid, TopAbs_FACE); fx.More(); fx.Next()) {
        for (TopExp_Explorer ex(fx.Current(), TopAbs_EDGE); ex.More(); ex.Next()) {
            const int idx = map.FindIndex(ex.Current());
            if (idx > 0) uses[idx]++;
        }
    }
    for (const auto& kv : uses) {
        const TopoDS_Edge& e = TopoDS::Edge(map.FindKey(kv.first));
        if (BRep_Tool::Degenerated(e)) continue;   // apex edge: bounds no material
        if (kv.second == 1) ++freeEdges;
        else if (kv.second > 2) ++nonManifold;
    }
}

bool hasArg(int argc, char** argv, const char* flag) {
    for (int i = 2; i < argc; ++i) if (!std::strcmp(argv[i], flag)) return true;
    return false;
}
double argVal(int argc, char** argv, const char* flag, double dflt) {
    for (int i = 2; i + 1 < argc; ++i) if (!std::strcmp(argv[i], flag)) return std::atof(argv[i + 1]);
    return dflt;
}
std::string argStr(int argc, char** argv, const char* flag, const std::string& dflt) {
    for (int i = 2; i + 1 < argc; ++i) if (!std::strcmp(argv[i], flag)) return argv[i + 1];
    return dflt;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: %s <file.step> [--mesh] [--airbox] [--max <m>] [--min <m>] [--curv <n>]\n"
            "                      [--quiet] [--msh <out.msh>]\n"
            "  --mesh    also run the gmsh import/fragment/tet stage (stage B)\n"
            "  --airbox  surround the assembly with the FEM air box before fragmenting\n"
            "  --max/--min  mesh size bounds [m] (default: derived from the bbox/feature size)\n"
            "  --curv    elements per 2*pi of curvature (default 8; 0 disables)\n", argv[0]);
        return 2;
    }
    const std::string path = argv[1];
    const bool doMesh  = hasArg(argc, argv, "--mesh");
    const bool airbox  = hasArg(argc, argv, "--airbox");
    const bool quiet   = hasArg(argc, argv, "--quiet");
    // --audit runs the pairwise overlap audit on the solids this tool has ALREADY read, instead
    // of paying for a second import of the same 10-20 MB STEP in omfem_step_intersect. The
    // import dominated both checks (Alf, 2026-09-03: "is there no way to make them faster,
    // without losing accuracy") and the audit itself is byte-for-byte the same code.
    const bool doAudit = hasArg(argc, argv, "--audit");
    // --list: one line per solid with its name, volume and bounding box. The way to VERIFY a
    // STEP's geometry from the STEP itself (Alf: "have you checked out the step?") -- the
    // generator's SVGs are MKF's 2D painter output, not a rendering of the solids.
    const bool doList  = hasArg(argc, argv, "--list");
    const double auditTol = argVal(argc, argv, "--audit-tol", 1e-6);
    const double curv  = argVal(argc, argv, "--curv", 8.0);
    std::string mshOut;
    for (int i = 2; i + 1 < argc; ++i) if (!std::strcmp(argv[i], "--msh")) mshOut = argv[i + 1];

    // ---------------- STAGE A: read the file back and check every solid ----------------
    auto t0 = std::chrono::steady_clock::now();
    Handle(TDocStd_Document) doc;
    XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
    STEPCAFControl_Reader reader;
    reader.SetNameMode(true);
    if (reader.ReadFile(path.c_str()) != IFSelect_RetDone || !reader.Transfer(doc)) {
        std::printf("STEPCHECK %s READ_FAIL\n", path.c_str());
        return 2;
    }
    Handle(XCAFDoc_ShapeTool) tool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
    TDF_LabelSequence labels;
    tool->GetFreeShapes(labels);

    std::vector<SolidInfo> solids;
    std::vector<mvb::cad::AuditPart> auditParts;   // only filled under --audit
    Bnd_Box assemblyBox;
    double minFeature_mm = std::numeric_limits<double>::max();
    // Every solid's own smallest bbox extent. minFeature is the MINIMUM of these, which any
    // single trim sliver or short riser can drive to a fraction of the wire; the percentiles
    // below describe the population the mesh actually has to resolve.
    std::vector<double> solidFeature_mm;

    // ------------------------------------------------------------------------------------
    // Stage A runs the per-solid checks IN PARALLEL. Each solid's verdict depends only on that
    // solid -- BRepCheck_Analyzer, the self-intersection analyzer, the edge-manifold count and
    // the surface histogram all read one shape and share nothing -- so the pass is embarrassingly
    // parallel, and it was the single most expensive thing either checker did (measured on the
    // 2026-09-03 corpus sweep: 538 s of stage A on complete_current_transformer's 961 solids,
    // 481 s on 12_boost's 466). The answers are IDENTICAL to the serial pass; only the wall clock
    // changes. OMFEM_CAD_JOBS sets the width (default 4, the same cap the OCC thread budget uses
    // on this shared box); OMFEM_CAD_JOBS=1 restores the serial order for debugging.
    // ------------------------------------------------------------------------------------
    struct SolidItem { TopoDS_Shape shape; std::string name, body; };
    std::vector<SolidItem> items;
    for (Standard_Integer i = 1; i <= labels.Length(); ++i) {
        Handle(TDataStd_Name) nm;
        std::string name = "shape" + std::to_string(i);
        if (labels.Value(i).FindAttribute(TDataStd_Name::GetID(), nm))
            name = std::string(TCollection_AsciiString(nm->Get()).ToCString());
        TopoDS_Shape sh = tool->GetShape(labels.Value(i));
        if (sh.IsNull()) continue;
        int k = 0;
        for (TopExp_Explorer ex(sh, TopAbs_SOLID); ex.More(); ex.Next(), ++k)
            items.push_back({ex.Current(),
                             name + (k ? " [solid " + std::to_string(k) + "]" : ""), name});
    }
    solids.assign(items.size(), SolidInfo{});
    {
        int jobs = std::getenv("OMFEM_CAD_JOBS") ? std::atoi(std::getenv("OMFEM_CAD_JOBS")) : 4;
        jobs = std::max(1, std::min<int>(jobs, static_cast<int>(items.size())));
        std::atomic<size_t> next{0};
        auto worker = [&]() {
            for (;;) {
                const size_t idx = next.fetch_add(1);
                if (idx >= items.size()) return;
                const SolidItem& item = items[idx];
                const TopoDS_Shape& sol = item.shape;
                SolidInfo& si = solids[idx];
        si.name = item.name;
        si.volume_mm3 = solidVolume(sol, &si.cx, &si.cy, &si.cz);
        si.analyzerOk = BRepCheck_Analyzer(sol, /*bopConsistency*/true).IsValid();
        // Self-intersection, the check BRepCheck above cannot make. On by default because a
        // self-intersecting solid is a defect, not a preference; OMFEM_NO_SELFINT=1 skips it
        // for the rare case where the analyzer itself is the thing being debugged (13's
        // --selfint has been measured taking >2400 s, so the cost is real and is reported).
        if (!std::getenv("OMFEM_NO_SELFINT")) {
            try {
                BOPAlgo_ArgumentAnalyzer an;
                an.SetShape1(sol);
                an.ArgumentTypeMode() = Standard_True;
                an.SelfInterMode()    = Standard_True;
                an.SmallEdgeMode()    = Standard_True;
                an.Perform();
                if (an.HasFaulty()) {
                    for (const auto& r : an.GetCheckResult()) {
                        if (r.GetCheckStatus() == BOPAlgo_SelfIntersect) ++si.selfIntCount;
                        else if (r.GetCheckStatus() == BOPAlgo_TooSmallEdge)
                            ++si.smallEdgeCount;
                    }
                    si.selfIntersects = si.selfIntCount > 0;
                }
            } catch (const Standard_Failure&) {
            }
        }
        edgeManifoldCheck(sol, si.freeEdges, si.nonManifoldEdges);
        for (TopExp_Explorer sx(sol, TopAbs_SHELL); sx.More(); sx.Next()) {
            ++si.shells;
            if (!BRep_Tool::IsClosed(sx.Current())) ++si.openShells;
        }
        for (TopExp_Explorer fx(sol, TopAbs_FACE); fx.More(); fx.Next()) ++si.faces;
        if (si.selfIntersects) {
            std::map<std::string, int> tm;
            for (TopExp_Explorer fx(sol, TopAbs_FACE); fx.More(); fx.Next()) {
                BRepAdaptor_Surface bs(TopoDS::Face(fx.Current()), Standard_False);
                const char* n = "other";
                switch (bs.GetType()) {
                    case GeomAbs_Plane:            n = "plane";    break;
                    case GeomAbs_Cylinder:         n = "cyl";      break;
                    case GeomAbs_Cone:             n = "cone";     break;
                    case GeomAbs_Sphere:           n = "sph";      break;
                    case GeomAbs_Torus:            n = "torus";    break;
                    case GeomAbs_BezierSurface:    n = "bezier";   break;
                    case GeomAbs_BSplineSurface:   n = "bspline";  break;
                    case GeomAbs_SurfaceOfRevolution: n = "revol"; break;
                    case GeomAbs_SurfaceOfExtrusion:  n = "extr";  break;
                    default: break;
                }
                tm[n]++;
            }
            for (const auto& kv : tm)
                si.surfMix += (si.surfMix.empty() ? "" : ",") + kv.first + ":" +
                              std::to_string(kv.second);
        }
        {
            Bnd_Box sb;
            BRepBndLib::Add(sol, sb);
            if (!sb.IsVoid()) {
                double a0, b0, c0, a1, b1, c1;
                sb.Get(a0, b0, c0, a1, b1, c1);
                double e[3] = {a1 - a0, b1 - b0, c1 - c0};
                std::sort(e, e + 3, std::greater<double>());
                si.ext0 = e[0]; si.ext1 = e[1]; si.ext2 = e[2];
            }
        }
            }
        };
        if (jobs == 1) worker();
        else {
            std::vector<std::thread> pool;
            for (int j = 0; j < jobs; ++j) pool.emplace_back(worker);
            for (auto& t : pool) t.join();
        }
    }
    for (size_t idx = 0; idx < items.size(); ++idx) {
        const TopoDS_Shape& sol = items[idx].shape;
        if (doAudit) {
            mvb::cad::AuditPart ap;
            ap.shape = sol;
            ap.body = items[idx].body;
            ap.name = items[idx].name;
            BRepBndLib::Add(sol, ap.box);
            ap.volume = solids[idx].volume_mm3;
            auditParts.push_back(std::move(ap));
        }
        Bnd_Box b;
        BRepBndLib::Add(sol, b);
        assemblyBox.Add(b);
        double x0, y0, z0, x1, y1, z1;
        b.Get(x0, y0, z0, x1, y1, z1);
        if (doList)
            std::printf("SOLID %-70s vol=%10.4f mm3  x[%8.3f,%8.3f] y[%8.3f,%8.3f] z[%8.3f,%8.3f]\n",
                        items[idx].name.c_str(), solids[idx].volume_mm3, x0, x1, y0, y1, z0, z1);
        // smallest extent of the smallest solid ~ the wire/feature the mesher must resolve
        const double feat = std::max(1e-6, std::min({x1 - x0, y1 - y0, z1 - z0}));
        solidFeature_mm.push_back(feat);
        minFeature_mm = std::min(minFeature_mm, feat);
    }
    const double readSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    if (solids.empty()) {
        std::printf("STEPCHECK %s NO_SOLIDS\n", path.c_str());
        return 1;
    }

    int badAnalyzer = 0, openSolids = 0, nonManifoldSolids = 0, zeroVol = 0, selfIntSolids = 0;
    int smallEdgeSolids = 0;
    double totalVol = 0;
    for (const auto& s : solids) {
        totalVol += s.volume_mm3;
        if (!s.analyzerOk) ++badAnalyzer;
        if (s.selfIntersects) ++selfIntSolids;
        if (s.smallEdgeCount) ++smallEdgeSolids;
        if (s.freeEdges || s.openShells) ++openSolids;
        if (s.nonManifoldEdges) ++nonManifoldSolids;
        if (!(s.volume_mm3 > 1e-12) || !std::isfinite(s.volume_mm3)) ++zeroVol;
    }
    if (!quiet) {
        for (const auto& s : solids) {
            const bool bad = !s.analyzerOk || s.freeEdges || s.openShells ||
                             s.nonManifoldEdges || !(s.volume_mm3 > 1e-12) || s.selfIntersects;
            if (bad)
                std::printf("  DEFECT %-44s vol=%.6f mm3 faces=%d ext=%.4fx%.4fx%.4f mm %s%s"
                            "%s%s%s%s%s\n",
                            s.name.c_str(), s.volume_mm3, s.faces, s.ext0, s.ext1, s.ext2,
                            s.surfMix.empty() ? "" : ("[" + s.surfMix + "] ").c_str(),
                            s.selfIntersects
                                ? ("[" + std::to_string(s.selfIntCount) +
                                   " SELF-INTERSECTION(s), BOPAlgo] ").c_str() : "",
                            s.analyzerOk ? "" : "[invalid B-Rep] ",
                            s.freeEdges ? ("[" + std::to_string(s.freeEdges) + " free edge(s)] ").c_str() : "",
                            s.openShells ? ("[" + std::to_string(s.openShells) + " open shell(s)] ").c_str() : "",
                            s.nonManifoldEdges ? ("[" + std::to_string(s.nonManifoldEdges) + " non-manifold edge(s)] ").c_str() : "",
                            (s.volume_mm3 > 1e-12) ? "" : "[zero/negative volume] ");
        }
    }

    double bx0, by0, bz0, bx1, by1, bz1;
    assemblyBox.Get(bx0, by0, bz0, bx1, by1, bz1);
    const double diag_mm = std::sqrt((bx1 - bx0) * (bx1 - bx0) + (by1 - by0) * (by1 - by0) +
                                     (bz1 - bz0) * (bz1 - bz0));

    // selfint counts as a stage-A defect: a self-intersecting solid is not sound geometry, and
    // calling such a file WATERTIGHT is exactly how 02_flyback shipped green for days.
    const bool stageA = (badAnalyzer == 0 && openSolids == 0 && nonManifoldSolids == 0 &&
                         zeroVol == 0 && selfIntSolids == 0);
    std::printf("STEPCHECK %s CAD %s solids=%zu vol=%.3f mm3 bboxdiag=%.2f mm minfeat=%.4f mm "
                "invalid=%d open=%d nonmanifold=%d zerovol=%d selfint=%d smalledge=%d "
                "read=%.1fs\n",
                path.c_str(), stageA ? "WATERTIGHT" : "DEFECTIVE", solids.size(), totalVol,
                diag_mm, minFeature_mm, badAnalyzer, openSolids, nonManifoldSolids, zeroVol,
                selfIntSolids, smallEdgeSolids, readSec);

    // ROBUST FEATURE SCALE. minFeature is an extremum: one sliver sets it for the whole
    // assembly. The percentile of the same population is what the bulk of the geometry
    // actually is -- for a wound part, the conductor cross-section.
    auto featurePercentile = [&](double q) {
        if (solidFeature_mm.empty()) return minFeature_mm;
        std::vector<double> v = solidFeature_mm;
        std::sort(v.begin(), v.end());
        const size_t idx = std::min(v.size() - 1,
                                    static_cast<size_t>(q * static_cast<double>(v.size() - 1)));
        return v[idx];
    };
    const double featP10_mm = featurePercentile(0.10);
    const double featP50_mm = featurePercentile(0.50);
    if (std::getenv("OMFEM_SIZE_DIAG"))
        std::printf("  SIZEDIAG solids=%zu minfeat=%.4f p10=%.4f p50=%.4f diag=%.3f "
                    "bboxRule=%.4f cap2xMin=%.4f\n",
                    solidFeature_mm.size(), minFeature_mm, featP10_mm, featP50_mm, diag_mm,
                    diag_mm / 25.0, 2.0 * minFeature_mm);

    // The overlap audit runs here: after the stage-A verdict (so a defective read still reports
    // both findings) and BEFORE the --mesh early return, because the sweep asks for the CAD
    // check and the audit together on a file it must not import twice.
    bool auditClean = true;
    if (doAudit) {
        std::printf("[intersect] %s:\n", path.c_str());
        std::fflush(stdout);
        auditClean = mvb::cad::audit_overlaps(auditParts, auditTol) == 0;
        std::fflush(stdout);
    }

    if (!doMesh) return (stageA && auditClean) ? 0 : 1;

    // ---------------- STAGE B: the gmsh path a FEM run actually takes ----------------
    // Sizes are metres here: the file declares mm and we ask OCC to convert on import, exactly
    // as MasMesher does (ABT #317) -- otherwise the geometry comes back 1000x too big.
    const double minFeature_m = minFeature_mm * 1e-3;
    const double diag_m = diag_mm * 1e-3;
    // Default szMax: bbox/25, but never coarser than 2x the smallest feature. Elements must
    // RESOLVE the wire: 07_cmc (0.8 mm wire, minfeat 0.49 mm) got szMax 1.635 mm from the
    // bbox rule and tetgen's recovery folded facets on the faceted tubes -- PLC errors, 1/235
    // meshed; capped at 0.8 mm the same file meshes 235/235 (2026-08-28). The failures scale
    // with domain size because --isolate shrank the bbox and "fixed" it, which is how a
    // geometry-looking failure turned out to be a sizing rule.
    // SIZING IS NOT WHAT 13/14 DIE OF -- do not "fix" them here (measured 2026-08-31):
    //   13_current_sense  --nofrag: MESHABLE 152/152 in  275 s   | with fragment: >7200 s
    //   14_dab            --nofrag: MESHABLE 132/132 in  162 s   | with fragment: >3000 s
    // Both mesh CORRECTLY at the sizes this rule already picks; the cost is entirely in OCC's
    // fragment, which runs BEFORE any meshing and is independent of szMax/szMin. Proof: with
    // a 7200 s budget 13 never printed "[stage] fragment+sync done" -- that line goes to
    // stderr (unbuffered), so its absence means the boolean never returned. Their failure
    // belongs with the fragment/geometry, not with mesh sizing.
    //
    // THE 2x CAP STAYS, because a design still needs it (also measured 2026-08-31):
    //   07_cmc  (the design the cap was added for) now meshes 303/303 at BOTH the capped
    //           0.9808 mm AND the uncapped bbox rule 1.635 mm (205 s) -- its old PLC failure
    //           was the coplanar-facet osculation class, fixed at the source by the junction
    //           stagger work, so 07 alone would no longer justify the cap.
    //   10_emi  DOES justify it: at 0.8455 mm (p50 rule) gmsh reports "Could not recover
    //           boundary mesh: error 2" and meshes 5/45; at the capped 0.3923 mm it is 45/45.
    //           Its thin solids are 0.196 mm -- 0.845 mm elements cannot resolve them. That is
    //           physics, not tuning, so the extremum-based cap is the safe rule and remains
    //           the DEFAULT.
    // OMFEM_SIZE_RULE=p50 selects the bulk-scale variant below for measurement only; it is
    // NOT the default precisely because it regresses 10_emi.
    //
    // BULK feature scale, not the extremum (OMFEM_SIZE_RULE=p50). Measured across the corpus
    // 2026-08-31 (per-solid smallest bbox extent, mm):
    //     design  minfeat   p10     p50     -> the wire is p50; minfeat is whatever sliver is
    //     02      0.3923  0.3950  0.4009      smallest anywhere.
    //     04      0.3925  0.3985  0.4997
    //     05      0.7846  0.7884  0.7968
    //     07      0.4904  0.4928  0.4977
    //     08      0.6179  1.4956  1.4956
    //     10      0.1962  0.1996  0.4227
    //     12      0.9808  0.9850  0.9957
    //     13      0.0981  0.0981  0.2264   <- 2.3x spread: sliver class vs real conductor
    // Where the winding is uniform (02/05/07/12) minfeat == p50 to ~2% and the rule is
    // unchanged. Where a sliver class exists (10/13) capping on minfeat forced the WHOLE
    // domain to the sliver scale: 13 got szMax 0.196 mm over an 18.59 mm part and blew a
    // 3000 s mesh budget. szMin stays tied to minFeature so small features are still resolved
    // LOCALLY (gmsh refines from boundary/curvature); only the global ceiling changes.
    // The 2x multiplier came from 07_cmc's PLC failure at 1.635 mm (pre-stagger). That failure
    // no longer reproduces -- 07 now meshes 303/303 at 1.635 mm -- so the multiplier is held
    // by 10_emi instead: its thinnest solids are 0.196 mm and it needs szMax 0.3923 mm (2x) to
    // recover their boundary, failing 5/45 at 0.8455 mm.
    const double bulkFeature_m =
        (std::getenv("OMFEM_SIZE_RULE") && std::string(std::getenv("OMFEM_SIZE_RULE")) == "p50")
            ? featP50_mm * 1e-3
            : minFeature_m;
    const double szMax = argVal(argc, argv, "--max",
                                std::min(std::max(diag_m / 25.0, bulkFeature_m),
                                         2.0 * bulkFeature_m));
    const double szMin = argVal(argc, argv, "--min", minFeature_m / 4.0);

    int rc = (stageA && auditClean) ? 0 : 1;
    auto tm = std::chrono::steady_clock::now();
    // FAIL FAST (Alf, 2026-08-27: "the moment something breaks or collides, stop all the
    // checking and investigate; don't run for hours when you know that there is something
    // wrong already"). A healthy assembly fragments and meshes in seconds to a few minutes;
    // the 03_buck lockstep check instead ground for 2 h inside OCC intersections before its
    // timeout, telling nobody anything a 10-minute cutoff would not have. --budget <sec>
    // (default 600) arms a hard watchdog over stage B: on expiry it prints an explicit
    // verdict naming the stage and exits 3 -- distinguishable from NOT_MESHABLE (1), so the
    // caller knows this was "too slow = investigate", not "meshed and failed".
    {
        const double budget = argVal(argc, argv, "--budget", 600.0);
        if (budget > 0) {
            static std::string budgetPath = path;
            static double budgetSec = budget;
            std::signal(SIGALRM, [](int) {
                std::printf("STEPCHECK %s MESH BUDGET_EXCEEDED after %.0fs -- fail fast: the "
                            "geometry is the suspect, not the clock\n",
                            budgetPath.c_str(), budgetSec);
                std::fflush(stdout);
                std::_Exit(3);
            });
            alarm(static_cast<unsigned>(budget));
        }
    }
    try {
        gmsh::initialize();
        // First error ends the run: grinding past a failure only buries the signal.
        // OMFEM_GMSH_RAW is the diagnostic exception: tetgen prints a PLC failure as an Error
        // line FOLLOWED by the offending coordinates at Info level, so aborting on the error
        // throws the one datum we need away.
        gmsh::option::setNumber("General.AbortOnError",
                                std::getenv("OMFEM_GMSH_RAW") ? 0 : 2);
        // Capture gmsh's own log. The per-volume "unmeshed" tally OVERSTATES the defect: gmsh
        // stops the 3D stage after a 2D failure, so ONE surface with overlapping facets can
        // leave 400 volumes empty. The surface-level errors are the atomic defects; count them.
        gmsh::logger::start();
        gmsh::option::setNumber("General.Terminal", std::getenv("OMFEM_GMSH_VERBOSE") ? 1 : 0);
        gmsh::option::setNumber("General.Verbosity", std::getenv("OMFEM_GMSH_VERBOSE") ? 99 : 1);
        gmsh::option::setString("Geometry.OCCTargetUnit", "M");
        gmsh::model::add("stepcheck");

        gmsh::vectorpair imported;
        gmsh::model::occ::importShapes(path, imported);
        gmsh::model::occ::synchronize();

        // --isolate N[,M,...]: keep only those volumes. Separates "this solid cannot be meshed
        // at all" from "it cannot be meshed in company" (a fragment/neighbour interaction).
        std::set<int> keep;
        for (int i = 2; i + 1 < argc; ++i) {
            if (std::strcmp(argv[i], "--isolate")) continue;
            std::string spec = argv[i + 1];
            size_t pos = 0;
            while (pos < spec.size()) {
                const size_t c = spec.find(',', pos);
                const int t = std::atoi(spec.substr(pos, c - pos).c_str());
                if (t > 0) keep.insert(t);
                if (c == std::string::npos) break;
                pos = c + 1;
            }
        }
        if (!keep.empty()) {
            gmsh::vectorpair all;
            gmsh::model::occ::getEntities(all, 3);
            gmsh::vectorpair drop;
            for (auto& dt : all) if (!keep.count(dt.second)) drop.push_back(dt);
            gmsh::model::occ::remove(drop, true);
            gmsh::model::occ::synchronize();
        }

        gmsh::vectorpair vols;
        gmsh::model::occ::getEntities(vols, 3);
        double volIn = 0;
        for (auto& dt : vols) { double m = 0; gmsh::model::occ::getMass(3, dt.second, m); volIn += m; }
        const size_t nIn = vols.size();

        if (airbox) {
            double x0, y0, z0, x1, y1, z1;
            gmsh::model::getBoundingBox(-1, -1, x0, y0, z0, x1, y1, z1);
            const double mx = 0.5 * (x0 + x1), my = 0.5 * (y0 + y1), mz = 0.5 * (z0 + z1);
            const double hx = 1.5 * (x1 - x0), hy = 1.5 * (y1 - y0), hz = 1.5 * (z1 - z0);
            const int box = gmsh::model::occ::addBox(mx - hx, my - hy, mz - hz, 2 * hx, 2 * hy, 2 * hz);
            vols.push_back({3, box});
            gmsh::model::occ::synchronize();
        }

        // Fragment the assembly into ONE conformal domain. GLUE + non-destructive is the mode
        // OCCT specifies for inputs that abut on exactly coincident faces with no true
        // intersections -- which is what conformal-by-construction windings are.
        const bool noglue = hasArg(argc, argv, "--noglue");
        const bool nofrag = hasArg(argc, argv, "--nofrag");
        // --touchfrag: fragment only the parts that actually TOUCH, never disjoint ones. The
        // conformality argument for the fragment only exists across a CONTACT: parts separated
        // by real gaps (the bare-copper product's enamel clearances) share no boundary, so
        // feeding them through the n-ary boolean buys nothing and costs everything -- measured
        // on 03_buck: 4 clean solids (audit: NO OVERLAPS, each meshes alone) came out of the
        // full fragment as 12 volumes with unrecoverable edges, 0 meshed. Touch-groups are
        // found by pairwise BRepExtrema distance <= tol on the gmsh shapes' OCC handles;
        // each group of >= 2 is fragmented on its own, singletons pass through untouched.
        const bool touchfrag = hasArg(argc, argv, "--touchfrag");
        if (touchfrag && !nofrag && vols.size() > 1) {
            std::vector<int> uf(vols.size());
            for (size_t i = 0; i < uf.size(); ++i) uf[i] = (int)i;
            std::function<int(int)> find = [&](int a) {
                while (uf[a] != a) a = uf[a] = uf[uf[a]];
                return a;
            };
            const double touchTol = 1e-7;   // metres; contact means contact, not proximity
            for (size_t i = 0; i < vols.size(); ++i)
                for (size_t j = i + 1; j < vols.size(); ++j) {
                    double d = 1e30, xx1, yy1, zz1, xx2, yy2, zz2;
                    try {
                        gmsh::model::occ::getDistance(3, vols[i].second, 3, vols[j].second,
                                                      d, xx1, yy1, zz1, xx2, yy2, zz2);
                    } catch (...) { continue; }
                    if (std::getenv("OMFEM_TOUCH_DEBUG")) {
                        // Name both volumes by nearest stage-A centroid, as the UNMESHED
                        // reporter does -- reading gmsh tags back to pieces by hand across the
                        // toroid counter-rotation cost a debugging hour and two wrong frames.
                        auto nameOf = [&](int tag) -> std::string {
                            double cx = 0, cy = 0, cz = 0;
                            try { gmsh::model::occ::getCenterOfMass(3, tag, cx, cy, cz); }
                            catch (...) { return "?"; }
                            cx *= 1e3; cy *= 1e3; cz *= 1e3;
                            const SolidInfo* best = nullptr;
                            double bd = 1e30;
                            for (const auto& si : solids) {
                                const double dd = std::fabs(si.cx - cx) + std::fabs(si.cy - cy) +
                                                  std::fabs(si.cz - cz);
                                if (dd < bd) { bd = dd; best = &si; }
                            }
                            return best ? best->name : "?";
                        };
                        std::printf("  touchdist vol%d('%s')-vol%d('%s'): %.6f mm  "
                                    "A(%.4f, %.4f, %.4f) mm\n",
                                    vols[i].second, nameOf(vols[i].second).c_str(),
                                    vols[j].second, nameOf(vols[j].second).c_str(), d * 1e3,
                                    xx1 * 1e3, yy1 * 1e3, zz1 * 1e3);
                    }
                    if (d <= touchTol) uf[find((int)i)] = find((int)j);
                }
            std::map<int, gmsh::vectorpair> groups;
            for (size_t i = 0; i < vols.size(); ++i) groups[find((int)i)].push_back(vols[i]);
            gmsh::option::setNumber("Geometry.OCCBooleanGlue", noglue ? 0 : 1);
            gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", noglue ? 0 : 1);
            size_t fragged_groups = 0;
            for (auto& kv : groups) {
                if (kv.second.size() < 2) continue;
                gmsh::vectorpair od;
                std::vector<gmsh::vectorpair> omap;
                gmsh::vectorpair obj{kv.second.front()},
                    tls(kv.second.begin() + 1, kv.second.end());
                try {
                    gmsh::model::occ::fragment(obj, tls, od, omap);
                    ++fragged_groups;
                } catch (const std::exception& e) {
                    std::printf("  touchfrag: group fragment failed (%s)\n", e.what());
                }
            }
            gmsh::model::occ::synchronize();
            gmsh::option::setNumber("Geometry.OCCBooleanGlue", 0);
            gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", 0);
            std::printf("  touchfrag: %zu group(s) of touching parts fragmented, %zu part(s) disjoint\n",
                        fragged_groups, groups.size() - fragged_groups);
        }
        // --chainfrag: fragment ALONG THE TOUCH GRAPH, one contact at a time, instead of one
        // n-ary fragment over everything.
        //
        // WHY (measured on 21_interleaved_flyback, 2026-08-31). The single n-ary fragment tests
        // every pair whose bounding boxes overlap. A helical turn's bbox is mostly AIR: in a
        // SECTIONED winding a turn only box-overlaps the chain neighbours it actually touches,
        // but INTERLEAVING threads the other winding's turns through that empty space, so boxes
        // interpenetrate while the SURFACES NEVER MEET. OCC then runs a full face-vs-face search
        // (~55x55 faces) to PROVE each such pair does not intersect -- pure waste, quadratic in
        // how many interleaved solids share the call. Full 21: 400 bbox pairs, only 72 TOUCHING,
        // 328 APART -- 82% of the pair work futile. Marginal cost per apart pair measured at a
        // near-constant 15.9 / 19.6 / 15.8 s; the late region grows n^2.6 (47/158/334/586 s at
        // n=10/15/20/26) while the early region, which has ZERO apart pairs, is flat n^0.74.
        //
        // Conformality is only required ACROSS A CONTACT, and contacts form a chain of O(n) (72
        // for 21's 76 solids), not the O(n^2) candidate set. Fragmenting each contact on its own
        // means two non-touching solids are NEVER arguments to the same boolean, so the futile
        // work cannot arise. Partitioning proof: 51..76 as two 13-solid groups = 227 s vs 586 s
        // together (2.6x).
        //
        // NOT the same as --touchfrag, which groups by touch-CONNECTIVITY: a winding is ONE
        // connected chain, so that never partitions and the futile pairs all stay in one call.
        //
        // *** OPT-IN, AND IT MUST STAY OPT-IN UNTIL THE PROBE IS CHEAP. *** The fragmenting is a
        // clear win (21: 3934 s -> 407 s, 9.7x, 76/76 meshed, volDev +0.0001%), but FINDING the
        // contacts costs more than it saves on this design (probe 1890 s) and is O(n^2) exact
        // BRepExtrema, so it gets catastrophically worse as solids multiply:
        //     21: 76 solids ->    400 bbox pairs -> 1890 s probe (4.7 s/pair)
        //     02: 588 solids -> ~25000 bbox pairs -> MEASURED >7780 s and still not done;
        //         extrapolates past 30 h. Two runs were killed at 7780 s / 2852 s.
        // 02 is a GREEN design through the default n-ary path, so switching it to --chainfrag
        // would turn minutes into hours. Make this the default only once the contact set is
        // obtained cheaply -- via --contacts-cache (a sweep re-checks the same STEP repeatedly),
        // a parallel probe (the pairs are independent; omfem_step_intersect already forks a pool
        // for exactly this), or a screen that actually discriminates. NOTE the face-box screen
        // below does NOT: it rejected 6 of 400 on 21, because a helical turn's faces are long
        // strips whose boxes are nearly the solid's own.
        const bool chainfrag = hasArg(argc, argv, "--chainfrag");
        std::string fragErrFirst;
        size_t chainContacts = 0, chainFutile = 0, chainDone = 0, chainFailed = 0;
        if (chainfrag && !nofrag && !touchfrag && vols.size() > 1) {
            // Bbox pre-filter first: only pairs whose boxes overlap can touch, and the exact
            // distance call is far dearer than the box test (2850 pairs -> ~400 on 21).
            std::vector<double> bb(6 * vols.size(), 0.0);
            for (size_t i = 0; i < vols.size(); ++i) {
                try {
                    gmsh::model::getBoundingBox(3, vols[i].second, bb[6 * i + 0], bb[6 * i + 1],
                                                bb[6 * i + 2], bb[6 * i + 3], bb[6 * i + 4],
                                                bb[6 * i + 5]);
                } catch (...) {}
            }
            auto boxOverlap = [&](size_t a, size_t b) {
                for (int k = 0; k < 3; ++k) {
                    if (bb[6 * a + k] > bb[6 * b + k + 3] + 1e-9) return false;
                    if (bb[6 * b + k] > bb[6 * a + k + 3] + 1e-9) return false;
                }
                return true;
            };
            const double touchTol = 1e-7;   // metres; contact means contact, not proximity

            // FACE-BOX SCREEN, between the solid-box test and the exact distance call.
            // MEASURED (21, 2026-08-31): the exact probe cost 1932 s of a 2337 s run -- 4.7x the
            // 407 s of fragmenting it exists to guard -- because on a pair that does NOT touch
            // BRepExtrema must search exhaustively to prove the minimum. Replacing each cheap
            // futile intersection with a dear exact test gave back most of the win.
            // A SOLID box is far too coarse here: a helical turn's box is mostly air, which is
            // the whole reason interleaving generates futile pairs. But two solids can only
            // touch if some FACE box of one meets some FACE box of the other -- pure arithmetic,
            // and it rejects the futile pairs before they ever reach OCC.
            std::vector<std::vector<double>> fbb(vols.size());
            for (size_t i = 0; i < vols.size(); ++i) {
                gmsh::vectorpair bnd;
                try { gmsh::model::getBoundary({vols[i]}, bnd, false, false, false); }
                catch (...) { continue; }
                fbb[i].reserve(6 * bnd.size());
                for (auto& f : bnd) {
                    double a0, a1, a2, a3, a4, a5;
                    try {
                        gmsh::model::getBoundingBox(2, std::abs(f.second), a0, a1, a2, a3, a4, a5);
                    } catch (...) { continue; }
                    fbb[i].insert(fbb[i].end(), {a0, a1, a2, a3, a4, a5});
                }
            }
            auto faceBoxesMeet = [&](size_t a, size_t b) {
                const std::vector<double>&A = fbb[a], &B = fbb[b];
                if (A.empty() || B.empty()) return true;   // no face data: defer to the exact test
                for (size_t p = 0; p + 5 < A.size(); p += 6)
                    for (size_t q = 0; q + 5 < B.size(); q += 6) {
                        bool sep = false;
                        for (int k = 0; k < 3 && !sep; ++k) {
                            if (A[p + k] > B[q + k + 3] + touchTol) sep = true;
                            if (B[q + k] > A[p + k + 3] + touchTol) sep = true;
                        }
                        if (!sep) return true;
                    }
                return false;
            };

            std::vector<std::pair<size_t, size_t>> contacts;
            size_t screened = 0;
            const auto tProbe0 = std::chrono::steady_clock::now();

            // --contacts-cache <file>: the contact set is a pure function of the STEP, but
            // proving it costs more than the fragmenting it enables (measured on 21: probe
            // 1890 s vs fragment 407 s -- BRepExtrema has to search exhaustively to prove a
            // minimum for a pair that does NOT touch, and a face-box screen cannot reject those
            // because a helical turn's faces are long strips with near-solid-sized boxes, 6 of
            // 400 rejected). A sweep re-checks the same STEP many times, so the probe is worth
            // remembering. The cache is keyed on the file's size+mtime+solid count and is
            // IGNORED (recomputed, never trusted blind) on any mismatch.
            const std::string cachePath = argStr(argc, argv, "--contacts-cache", "");
            bool cacheHit = false;
            struct stat stbuf;
            const bool haveStat = (::stat(path.c_str(), &stbuf) == 0);
            if (!cachePath.empty() && haveStat) {
                if (std::FILE* cf = std::fopen(cachePath.c_str(), "r")) {
                    long long csize = -1, cmtime = -1;
                    long cn = -1;
                    char line[512];
                    std::vector<std::pair<size_t, size_t>> got;
                    while (std::fgets(line, sizeof line, cf)) {
                        long long a, b;
                        if (std::sscanf(line, "size %lld", &a) == 1) csize = a;
                        else if (std::sscanf(line, "mtime %lld", &a) == 1) cmtime = a;
                        else if (std::sscanf(line, "solids %lld", &a) == 1) cn = (long)a;
                        else if (std::sscanf(line, "c %lld %lld", &a, &b) == 2)
                            got.push_back({(size_t)a, (size_t)b});
                    }
                    std::fclose(cf);
                    if (csize == (long long)stbuf.st_size && cmtime == (long long)stbuf.st_mtime &&
                        cn == (long)vols.size()) {
                        contacts = got;
                        cacheHit = true;
                        std::fprintf(stderr, "  chainfrag: %zu contact(s) loaded from cache %s\n",
                                     contacts.size(), cachePath.c_str());
                    } else {
                        std::fprintf(stderr, "  chainfrag: cache %s STALE (size/mtime/solids "
                                     "mismatch) -- reprobing\n", cachePath.c_str());
                    }
                }
            }

            if (!cacheHit)
            for (size_t i = 0; i < vols.size(); ++i)
                for (size_t j = i + 1; j < vols.size(); ++j) {
                    if (!boxOverlap(i, j)) continue;
                    if (!faceBoxesMeet(i, j)) { ++chainFutile; ++screened; continue; }
                    double d = 1e30, x1, y1, z1, x2, y2, z2;
                    try {
                        gmsh::model::occ::getDistance(3, vols[i].second, 3, vols[j].second,
                                                      d, x1, y1, z1, x2, y2, z2);
                    } catch (...) { continue; }
                    if (d <= touchTol) contacts.push_back({i, j});
                    else ++chainFutile;
                }
            chainContacts = contacts.size();
            const double probeSec = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - tProbe0).count();
            if (!cachePath.empty() && !cacheHit && haveStat) {
                if (std::FILE* cf = std::fopen(cachePath.c_str(), "w")) {
                    std::fprintf(cf, "# omfem chainfrag contacts v1\n");
                    std::fprintf(cf, "step %s\n", path.c_str());
                    std::fprintf(cf, "size %lld\n", (long long)stbuf.st_size);
                    std::fprintf(cf, "mtime %lld\n", (long long)stbuf.st_mtime);
                    std::fprintf(cf, "solids %zu\n", vols.size());
                    for (auto& c : contacts)
                        std::fprintf(cf, "c %zu %zu\n", c.first, c.second);
                    std::fclose(cf);
                    std::fprintf(stderr, "  chainfrag: wrote %zu contact(s) to cache %s\n",
                                 contacts.size(), cachePath.c_str());
                }
            }

            // Each original solid may be SPLIT by a contact, so track the live pieces it became
            // and feed all of them to its next contact -- otherwise a later contact would be
            // fragmented against a stale tag that the boolean already consumed.
            std::vector<gmsh::vectorpair> cur(vols.size());
            for (size_t i = 0; i < vols.size(); ++i) cur[i] = {vols[i]};

            gmsh::option::setNumber("Geometry.OCCBooleanGlue", noglue ? 0 : 1);
            gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", noglue ? 0 : 1);
            for (auto& c : contacts) {
                if (cur[c.first].empty() || cur[c.second].empty()) continue;
                gmsh::vectorpair obj = cur[c.first], tls = cur[c.second], od;
                std::vector<gmsh::vectorpair> omap;
                try {
                    gmsh::model::occ::fragment(obj, tls, od, omap);
                } catch (const std::exception& e) {
                    ++chainFailed;
                    if (fragErrFirst.empty()) fragErrFirst = e.what();
                    continue;
                }
                ++chainDone;
                gmsh::vectorpair na, nb;
                for (size_t k = 0; k < omap.size(); ++k)
                    for (auto& t : omap[k]) (k < obj.size() ? na : nb).push_back(t);
                auto dedupe = [](gmsh::vectorpair& v) {
                    std::sort(v.begin(), v.end());
                    v.erase(std::unique(v.begin(), v.end()), v.end());
                };
                dedupe(na); dedupe(nb);
                if (!na.empty()) cur[c.first] = na;
                if (!nb.empty()) cur[c.second] = nb;
            }
            gmsh::model::occ::synchronize();
            gmsh::option::setNumber("Geometry.OCCBooleanGlue", 0);
            gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", 0);
            const double fragSec = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - tProbe0).count() - probeSec;
            std::fprintf(stderr,
                        "  chainfrag: %zu contact(s) fragmented in %.1fs (%zu failed); "
                        "%zu futile pair(s) SKIPPED that one n-ary fragment would have tested; "
                        "contact probe %.1fs over %zu bbox pair(s) (%zu rejected by face-box "
                        "screen, %zu needed the exact test)\n",
                        chainDone, fragSec, chainFailed, chainFutile,
                        probeSec, chainContacts + chainFutile, screened,
                        chainContacts + chainFutile - screened);
        }

        std::string fragErr = fragErrFirst;
        bool fragged = touchfrag || chainfrag;
        if (vols.size() > 1 && !nofrag && !touchfrag && !chainfrag) {
            gmsh::option::setNumber("Geometry.OCCBooleanGlue", noglue ? 0 : 1);
            gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", noglue ? 0 : 1);
            // --fuzzy F pins the boolean tolerance [m] instead of climbing the ladder. MVB++'s
            // ConductorBuilder welds its chain primitives at fuzzy=1e-7 m because consecutive
            // sections are "co-located but not bit-coincident"; the consumer's fragment needs
            // the same headroom to see those junctions as one.
            const double pinned = argVal(argc, argv, "--fuzzy", -1.0);
            const double ladder[] = {0.0, 1e-9, 1e-7};
            std::vector<double> tolLadder;
            if (pinned >= 0) tolLadder.push_back(pinned);
            else tolLadder.assign(std::begin(ladder), std::end(ladder));
            for (double tol : tolLadder) {
                try {
                    gmsh::option::setNumber("Geometry.ToleranceBoolean", tol);
                    gmsh::vectorpair od;
                    std::vector<gmsh::vectorpair> omap;
                    gmsh::vectorpair obj{vols[0]}, tools(vols.begin() + 1, vols.end());
                    gmsh::model::occ::fragment(obj, tools, od, omap);
                    gmsh::model::occ::synchronize();
                    fragged = true;
                    break;
                } catch (const std::exception& e) {
                    fragErr = e.what();
                }
            }
            gmsh::option::setNumber("Geometry.ToleranceBoolean", 0.0);
            gmsh::option::setNumber("Geometry.OCCBooleanGlue", 0);
            gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", 0);
        } else {
            fragged = true;
        }
        if (!fragged) {
            std::printf("STEPCHECK %s MESH FRAGMENT_FAIL (%s)\n", path.c_str(), fragErr.c_str());
            gmsh::finalize();
            return 1;
        }

        gmsh::vectorpair vout;
        gmsh::model::occ::getEntities(vout, 3);
        double volOut = 0;
        for (auto& dt : vout) { double m = 0; gmsh::model::occ::getMass(3, dt.second, m); volOut += m; }
        // Fragment must conserve material. With the air box the total legitimately grows to the
        // box volume, so compare only the assembly's own volume when there is no box.
        const double volDev = (!airbox && volIn > 0) ? (volOut - volIn) / volIn : 0.0;

        // --heal-frag: repair the FRAGMENT'S OUTPUT before meshing. Measured on 03_buck: the
        // input STEP is clean, but the boolean hands gmsh two faces whose wires are misoriented
        // (surface area integrates NEGATIVE, -1.4e-2 mm2 over the whole winding band) and edge
        // recovery then fails on exactly those faces. healShapes runs OCCT's ShapeFix over the
        // model -- fixing wire orientation is precisely its job. Applied only on request so the
        // default path still measures the raw boolean.
        if (hasArg(argc, argv, "--heal-frag")) {
            try {
                gmsh::vectorpair healed;
                gmsh::model::occ::healShapes(healed);
                gmsh::model::occ::synchronize();
                std::printf("  healed fragment: %zu top-level shape(s)\n", healed.size());
            } catch (const std::exception& e) {
                std::printf("  heal-frag failed: %s\n", e.what());
            }
            gmsh::vectorpair vh;
            gmsh::model::occ::getEntities(vh, 3);
            vout = vh;
            volOut = 0;
            for (auto& dt : vout) { double m = 0; gmsh::model::occ::getMass(3, dt.second, m); volOut += m; }
        }

        gmsh::option::setNumber("Mesh.MeshSizeMax", szMax);
        gmsh::option::setNumber("Mesh.MeshSizeMin", szMin);
        gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", curv);
        gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 1);
        // --wire-refine [n]: ADAPTIVE SIZING AROUND THE CONDUCTORS (Alf, 2026-09-03: "can we
        // have finer mesh around the wires?"). The uniform ceiling this gate derives from the
        // part's minimum feature is what folds the surface mesh on round wire: 01_simple gets
        // szMax 1.6 mm for a 1 mm tube, and gmsh lays triangles that overlap ("invalid boundary
        // mesh"). OMFEM's own FEM mesher never does this -- mvb::mesh::SizeFieldBuilder
        // composes skin-depth, gap and component sizes into one Min field -- so the gate now
        // borrows the same builder for the one feature it can identify without part names:
        // a solid whose smallest bounding-box extent is far under the part's scale is a
        // conductor piece, and its own faces get n elements across that thickness (default 3),
        // graded back to the ceiling. Everything else keeps the ceiling.
        int wireRefineN = 0;
        for (int i = 2; i < argc; ++i) {
            if (std::strcmp(argv[i], "--wire-refine")) continue;
            wireRefineN = 3;
            if (i + 1 < argc && argv[i + 1][0] != '-') wireRefineN = std::atoi(argv[i + 1]);
            break;
        }
        if (wireRefineN > 0) {
            // One implementation, shared with the FEM mesher: omfem::mesh composes every size
            // source into a single Min background field (SizeField.hpp). Here the only feature
            // identifiable without part names is the conductor's own thinness.
            const double thinCut = 3.0 * minFeature_m;
            try {
                mvb::mesh::SizeFieldBuilder sb(3, szMax);
                sb.add_ceiling(szMax);
                const auto r = mvb::mesh::add_thin_solid_refinement(sb, thinCut, wireRefineN);
                if (r.faces == 0) {
                    std::printf("  wire-refine: no solid thinner than %.4f mm; nothing to refine\n",
                                thinCut * 1e3);
                }
                else {
                    sb.finalize(r.size_min);
                    std::printf("  wire-refine: %zu thin solid(s), thinnest %.4f mm -> %d elements "
                                "across (size %.4f mm), %zu face(s) refined, ceiling %.4f mm\n",
                                r.solids, r.thinnest * 1e3, wireRefineN, r.size_min * 1e3, r.faces,
                                szMax * 1e3);
                }
            } catch (const std::exception& e) {
                std::printf("  wire-refine FAILED: %s\n", e.what());
            }
        }
        // --dump-frag: write the POST-FRAGMENT geometry so it can be diagnosed with the same
        // OCCT battery as the input. This separates "the STEP shipped a face gmsh cannot mesh"
        // from "the boolean CREATED one" -- measured on 03_buck, whose input has zero seam
        // faces yet whose mesh fails on a periodic surface.
        for (int i = 2; i + 1 < argc; ++i) {
            if (std::strcmp(argv[i], "--dump-frag")) continue;
            try { gmsh::write(argv[i + 1]); std::printf("  wrote post-fragment geometry %s\n", argv[i + 1]); }
            catch (const std::exception& e) { std::printf("  dump-frag failed: %s\n", e.what()); }
        }
        // OMFEM_ALGO2D overrides the 2D algorithm (6 = frontal-Delaunay default; 5 = Delaunay,
        // 1 = MeshAdapt) -- frontal-Delaunay can fold facets on trimmed bspline patches the
        // fragment carves at exact-abutment junctions (12_boost surface 4560).
        {
            const char* a2 = std::getenv("OMFEM_ALGO2D");
            gmsh::option::setNumber("Mesh.Algorithm", a2 ? std::atof(a2) : 6);
        }
        // OMFEM_ALGO3D overrides the 3D algorithm (1 = Delaunay default; 10 = HXT, parallel
        // and much faster in boundary recovery -- never MMG3D by default).
        {
            const char* a3 = std::getenv("OMFEM_ALGO3D");
            gmsh::option::setNumber("Mesh.Algorithm3D", a3 ? std::atof(a3) : 1);
        }
        // A volume the mesher cannot fill makes gmsh throw, which would hide every OTHER
        // volume's verdict. The interesting answer is "which ones and how many", so swallow
        // the throw and let the per-volume audit below report the whole picture.
        std::string genErr;
        // OMFEM_CURVE_PROBE="c1,c2,...": after the fragment, print each curve's bbox and its
        // owning surfaces/volumes (named by stage-A centroid) and exit -- for chasing gmsh's
        // "There are N intersections in the 1D mesh (curves a b)" retry loop back to parts.
        if (const char* cp = std::getenv("OMFEM_CURVE_PROBE")) {
            std::string list(cp);
            size_t pos = 0;
            while (pos < list.size()) {
                const int cTag = std::atoi(list.c_str() + pos);
                double bx0, by0, bz0, bx1, by1, bz1;
                try {
                    gmsh::model::getBoundingBox(1, cTag, bx0, by0, bz0, bx1, by1, bz1);
                    std::printf("  CURVE %d bbox x[%.4f,%.4f] y[%.4f,%.4f] z[%.4f,%.4f] mm\n",
                                cTag, bx0 * 1e3, bx1 * 1e3, by0 * 1e3, by1 * 1e3, bz0 * 1e3, bz1 * 1e3);
                    std::vector<int> up, down;
                    gmsh::model::getAdjacencies(1, cTag, up, down);
                    for (int s : up) {
                        std::vector<int> vup, vdown;
                        gmsh::model::getAdjacencies(2, s, vup, vdown);
                        for (int v : vup) {
                            double ccx = 0, ccy = 0, ccz = 0;
                            try { gmsh::model::occ::getCenterOfMass(3, v, ccx, ccy, ccz); } catch (...) {}
                            ccx *= 1e3; ccy *= 1e3; ccz *= 1e3;
                            const SolidInfo* best = nullptr;
                            double bestd = std::numeric_limits<double>::max();
                            for (const auto& si : solids) {
                                const double d = std::fabs(si.cx - ccx) + std::fabs(si.cy - ccy) +
                                                 std::fabs(si.cz - ccz);
                                if (d < bestd) { bestd = d; best = &si; }
                            }
                            std::printf("    on surface %d of volume %d -> '%s'\n", s, v,
                                        best ? best->name.c_str() : "?");
                        }
                    }
                } catch (const std::exception& e) {
                    std::printf("  CURVE %d probe failed: %s\n", cTag, e.what());
                }
                pos = list.find(',', pos);
                if (pos == std::string::npos) break;
                ++pos;
            }
            gmsh::finalize();
            return 0;
        }
        // OMFEM_MESH_TRACE: stream gmsh's own Info lines live to the terminal so a hung run's
        // log tail names the exact surface/volume being meshed when the watchdog kills it.
        if (std::getenv("OMFEM_MESH_TRACE")) {
            gmsh::option::setNumber("General.Terminal", 1);
            gmsh::option::setNumber("General.Verbosity", 4);
        }
        // Stage timestamps: a BUDGET_EXCEEDED that names no stage cost an hour of guessing
        // whether 04's hang was the fragment or the mesh. stderr, flushed, so the watchdog
        // kill cannot lose them.
        const auto tFrag = std::chrono::steady_clock::now();
        std::fprintf(stderr, "[stage] fragment+sync done, starting 2D mesh\n");
        // CURVATURE RETRY LADDER (2026-08-31, measured on 17_cllc).
        // The default 8 elements per 2*pi of curvature leaves a chord deviation on swept
        // lateral faces big enough that adjacent pieces' SURFACE meshes interpenetrate near
        // tangent junctions -- tetgen then dies in boundary recovery with "PLC Error: a
        // segment and a facet intersect", even though the CAD is watertight and the overlap
        // audit is clean. 17_cllc measured: curv 8/10/12 fail, 16/20/24 give 378/378, with the
        // element size unchanged -- so curvature, not size, is the operative variable.
        // A blanket increase is the wrong fix: it costs ~7x on designs that already pass
        // (10_emi 202 s -> ~25 min at curv 16) and would worsen the ones that are cost-bound.
        // So retry ONLY on a boundary-recovery failure: designs that pass at the default never
        // pay, and the gate answers the honest question -- is this geometry meshable at a
        // reasonable setting? -- instead of failing it on a discretisation default.
        auto recoveryFailure = [](const std::string& e) {
            return e.find("PLC Error") != std::string::npos ||
                   e.find("recover") != std::string::npos ||
                   e.find("Invalid boundary mesh") != std::string::npos;
        };
        const double curvLadder[] = {0.0, 16.0, 24.0};   // 0 = keep the caller's setting
        for (size_t attempt = 0; attempt < sizeof(curvLadder) / sizeof(curvLadder[0]); ++attempt) {
            if (attempt > 0) {
                if (curv <= 0 || curvLadder[attempt] <= curv) break;   // caller already finer
                std::fprintf(stderr, "[stage] retrying at curv=%.0f\n", curvLadder[attempt]);
                // REBUILD THE MODEL, do not just clear the mesh. mesh::clear() leaves state
                // that makes the retry land one volume short of a clean run (measured on
                // 17_cllc: in-process retry 377/378 vs 378/378 for a fresh --curv 16 run), so
                // the retry re-imports and re-fragments exactly as the first attempt did.
                // Cheap where it matters: the ladder only fires on boundary-recovery failures,
                // and a design whose FRAGMENT is the bottleneck (13/14) never reaches here.
                try {
                    gmsh::clear();
                    gmsh::vectorpair reimported;
                    gmsh::model::occ::importShapes(path, reimported);
                    gmsh::model::occ::synchronize();
                    gmsh::vectorpair vols2;
                    gmsh::model::occ::getEntities(vols2, 3);
                    if (vols2.size() > 1 && !nofrag) {
                        gmsh::option::setNumber("Geometry.OCCBooleanGlue", noglue ? 0 : 1);
                        gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", noglue ? 0 : 1);
                        gmsh::vectorpair od2;
                        std::vector<gmsh::vectorpair> omap2;
                        gmsh::vectorpair obj2{vols2[0]}, tools2(vols2.begin() + 1, vols2.end());
                        gmsh::model::occ::fragment(obj2, tools2, od2, omap2);
                        gmsh::model::occ::synchronize();
                        gmsh::option::setNumber("Geometry.OCCBooleanGlue", 0);
                        gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", 0);
                    }
                    vout.clear();
                    gmsh::model::occ::getEntities(vout, 3);
                    gmsh::option::setNumber("Mesh.MeshSizeMax", szMax);
                    gmsh::option::setNumber("Mesh.MeshSizeMin", szMin);
                    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 1);
                    gmsh::option::setNumber("Mesh.Algorithm", 6);
                    gmsh::option::setNumber("Mesh.Algorithm3D", 1);
                } catch (const std::exception& e) {
                    std::fprintf(stderr, "[stage] retry rebuild failed: %s\n", e.what());
                    break;
                }
                gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", curvLadder[attempt]);
                genErr.clear();
            }
            try {
                gmsh::model::mesh::generate(2);
                if (attempt == 0)
                    std::fprintf(stderr, "[stage] 2D mesh done in %.1fs, starting 3D\n",
                                 std::chrono::duration<double>(
                                     std::chrono::steady_clock::now() - tFrag).count());
                gmsh::model::mesh::generate(3);
            } catch (const std::exception& e) {
                genErr = e.what();
            }
            // Retry only when the failure is the recoverable class AND volumes are missing.
            std::vector<int> et0;
            std::vector<std::vector<std::size_t>> tg0, nd0;
            gmsh::model::mesh::getElements(et0, tg0, nd0, 3);
            std::size_t tets0 = 0;
            for (size_t i = 0; i < et0.size(); ++i)
                if (et0[i] == 4) tets0 += tg0[i].size();
            std::vector<std::string> glog0;
            gmsh::logger::get(glog0);
            bool recoverable = recoveryFailure(genErr);
            for (const auto& ln : glog0)
                if (recoveryFailure(ln)) { recoverable = true; break; }
            if (tets0 == 0) break;
            // DEGENERATE ELEMENTS ARE THE SAME QUESTION, ASKED LATER. A run can finish with no
            // gmsh error and every volume filled and still be unusable: two flat tets inside a
            // 0.25 mm-radius terminal wire (two_switch_forward, curv 8, SICN -4e-15) are the
            // same under-resolved curvature that makes tetgen die outright elsewhere -- the
            // wire is the smallest cylinder in the model and 8 elements per 2*pi leaves a
            // 0.196 mm chord under a 0.094 mm element size. curv 16 clears it (0 inverted,
            // minSICN +0.0004). So the ladder answers for quality too, on the same terms: only
            // a run that is actually defective pays for the finer rung.
            std::size_t inverted0 = 0;
            {
                std::vector<std::size_t> tets0Tags;
                for (size_t i = 0; i < et0.size(); ++i)
                    if (et0[i] == 4) tets0Tags.insert(tets0Tags.end(), tg0[i].begin(), tg0[i].end());
                std::vector<double> q0;
                gmsh::model::mesh::getElementQualities(tets0Tags, q0, "minSICN");
                for (double v : q0)
                    if (v <= 0) ++inverted0;
            }
            if (!recoverable && inverted0 == 0) break;
            std::set<int> meshed0;
            for (auto& dt : vout) {
                std::vector<int> et1;
                std::vector<std::vector<std::size_t>> tg1, nd1;
                gmsh::model::mesh::getElements(et1, tg1, nd1, 3, dt.second);
                for (size_t i = 0; i < et1.size(); ++i)
                    if (et1[i] == 4 && !tg1[i].empty()) { meshed0.insert(dt.second); break; }
            }
            // Complete AND sound: nothing to retry for. Incomplete or degenerate: climb.
            if (meshed0.size() == vout.size() && inverted0 == 0) break;
            if (inverted0 > 0)
                std::fprintf(stderr, "[stage] %zu degenerate tet(s) at this curvature\n", inverted0);
        }

        std::vector<int> etypes;
        std::vector<std::vector<std::size_t>> etags, enodes;
        gmsh::model::mesh::getElements(etypes, etags, enodes, 3);
        std::size_t tets = 0;
        std::vector<std::size_t> allTets;
        for (size_t i = 0; i < etypes.size(); ++i)
            if (etypes[i] == 4) { tets += etags[i].size(); allTets.insert(allTets.end(), etags[i].begin(), etags[i].end()); }

        double minQ = 1.0;
        std::size_t inverted = 0;
        if (!allTets.empty()) {
            std::vector<double> q;
            gmsh::model::mesh::getElementQualities(allTets, q, "minSICN");
            for (size_t qi = 0; qi < q.size(); ++qi) {
                minQ = std::min(minQ, q[qi]);
                if (q[qi] <= 0) {
                    ++inverted;
                    // WHERE is it? A degenerate tet at a coincident unmerged interface is a
                    // geometry verdict, and its position names the junction (04_forward: one
                    // zero-volume tet was the entire NOT_MESHABLE).
                    if (inverted <= 6) {
                        try {
                            int etype = 0, edim = 0, etag = 0;
                            std::vector<std::size_t> enodes2;
                            gmsh::model::mesh::getElement(allTets[qi], etype, enodes2, edim, etag);
                            double cx = 0, cy = 0, cz = 0;
                            for (std::size_t nid : enodes2) {
                                std::vector<double> coord, pcoord;
                                int ndim = 0, ntag = 0;
                                gmsh::model::mesh::getNode(nid, coord, pcoord, ndim, ntag);
                                if (coord.size() >= 3) { cx += coord[0]; cy += coord[1]; cz += coord[2]; }
                            }
                            const double n4 = std::max<size_t>(1, enodes2.size());
                            std::printf("  INVERTED tet SICN=%.3e at (%.4f, %.4f, %.4f) mm in volume %d\n",
                                        q[qi], cx / n4 * 1e3, cy / n4 * 1e3, cz / n4 * 1e3, etag);
                        } catch (...) {
                        }
                    }
                }
            }
        }
        // How many of the fragmented volumes actually received tets? A volume the mesher
        // silently skipped is a hole in the FEM domain, not a success. Name each empty one by
        // matching its centre of mass back to the STEP solid it came from.
        std::set<int> meshedVols;
        std::vector<int> emptyVols;
        for (auto& dt : vout) {
            std::vector<int> et;
            std::vector<std::vector<std::size_t>> tg, nd;
            gmsh::model::mesh::getElements(et, tg, nd, 3, dt.second);
            bool got = false;
            for (size_t i = 0; i < et.size(); ++i)
                if (et[i] == 4 && !tg[i].empty()) { got = true; break; }
            if (got) meshedVols.insert(dt.second);
            else emptyVols.push_back(dt.second);
        }
        for (int v : emptyVols) {
            double m = 0;
            gmsh::model::occ::getMass(3, v, m);
            double ccx = 0, ccy = 0, ccz = 0;
            gmsh::model::occ::getCenterOfMass(3, v, ccx, ccy, ccz);
            ccx *= 1e3; ccy *= 1e3; ccz *= 1e3;    // metres (import unit) -> mm (stage A unit)
            const SolidInfo* best = nullptr;
            double bestd = std::numeric_limits<double>::max();
            for (const auto& s : solids) {
                const double d = std::fabs(s.cx - ccx) + std::fabs(s.cy - ccy) + std::fabs(s.cz - ccz);
                if (d < bestd) { bestd = d; best = &s; }
            }
            std::printf("  UNMESHED volume %d  %.6f mm3  at (%.3f, %.3f, %.3f) mm -> '%s' (match %.4f mm)\n",
                        v, m * 1e9, ccx, ccy, ccz, best ? best->name.c_str() : "?", bestd);
        }
        if (!genErr.empty()) std::printf("  gmsh: %s\n", genErr.c_str());

        // Distinct gmsh errors, deduplicated: "3 x Invalid boundary mesh (overlapping facets)"
        // is the honest defect count, not the 400 volumes it left empty downstream.
        std::vector<std::string> glog;
        gmsh::logger::get(glog);
        std::map<std::string, int> errKinds;
        int nErr = 0;
        for (const auto& line : glog) {
            // gmsh logs as "Error: msg" through the logger and "Error   : msg" on the terminal
            const size_t colon = line.find(':');
            if (colon == std::string::npos) continue;
            std::string lvl = line.substr(0, colon);
            while (!lvl.empty() && lvl.back() == ' ') lvl.pop_back();
            if (lvl != "Error") continue;
            ++nErr;
            // OMFEM_GMSH_RAW: print every gmsh Error verbatim, before the collapsing below drops
            // the coordinates. The aggregate is right for a sweep, but a PLC intersection is only
            // actionable with the POINT it names -- that is what locates the offending solids.
            if (std::getenv("OMFEM_GMSH_RAW")) std::printf("  GMSH_RAW %s\n", line.c_str());
            std::string msg = line.substr(colon + 1);
            while (!msg.empty() && msg.front() == ' ') msg.erase(msg.begin());
            // strip trailing entity tags so the same failure class collapses to one key
            const size_t on = msg.find(" on surface ");
            if (on != std::string::npos) {
                // Before collapsing, name the CULPRIT: gmsh's per-volume failures list the
                // first VICTIMS, but the failing surface pins the actual defect. Resolve the
                // surface tag to its adjacent volumes and name them by stage-A centroid --
                // 12_boost bisected clean pairwise precisely because the reported empty
                // volumes were downstream of the one poisoned surface.
                // The message can name TWO surfaces ("... on surface A surface B"): facets of
                // A overlap facets of B. Resolve both.
                std::vector<int> sTags;
                {
                    const char* p = msg.c_str() + on;
                    while ((p = std::strstr(p, "surface ")) != nullptr) {
                        const int t = std::atoi(p + 8);
                        if (t > 0) sTags.push_back(t);
                        p += 8;
                    }
                }
                for (const int sTag : sTags) {
                    std::vector<int> up, down;
                    try { gmsh::model::getAdjacencies(2, sTag, up, down); } catch (...) {}
                    for (int v : up) {
                        double ccx = 0, ccy = 0, ccz = 0;
                        try { gmsh::model::occ::getCenterOfMass(3, v, ccx, ccy, ccz); }
                        catch (...) { continue; }
                        ccx *= 1e3; ccy *= 1e3; ccz *= 1e3;
                        const SolidInfo* best = nullptr;
                        double bestd = std::numeric_limits<double>::max();
                        for (const auto& s : solids) {
                            const double d = std::fabs(s.cx - ccx) + std::fabs(s.cy - ccy) +
                                             std::fabs(s.cz - ccz);
                            if (d < bestd) { bestd = d; best = &s; }
                        }
                        double scx = 0, scy = 0, scz = 0;
                        try { gmsh::model::occ::getCenterOfMass(2, sTag, scx, scy, scz); }
                        catch (...) {}
                        double sArea = 0;
                        try { gmsh::model::occ::getMass(2, sTag, sArea); } catch (...) {}
                        double bx0 = 0, by0 = 0, bz0 = 0, bx1 = 0, by1 = 0, bz1 = 0;
                        try { gmsh::model::getBoundingBox(2, sTag, bx0, by0, bz0, bx1, by1, bz1); }
                        catch (...) {}
                        std::printf("  CULPRIT surface %d at (%.3f, %.3f, %.3f) mm area=%.6f mm2 "
                                    "bbox x[%.3f,%.3f] y[%.3f,%.3f] z[%.3f,%.3f] mm "
                                    "bounds volume %d -> '%s' (match %.4f mm)\n",
                                    sTag, scx * 1e3, scy * 1e3, scz * 1e3, sArea * 1e6,
                                    bx0 * 1e3, bx1 * 1e3, by0 * 1e3, by1 * 1e3, bz0 * 1e3, bz1 * 1e3,
                                    v, best ? best->name.c_str() : "?", bestd);
                    }
                }
                msg = msg.substr(0, on) + " on surface <n>";
            }
            else if (msg.compare(0, 20, "No elements in volume") == 0) msg = "No elements in volume <n>";
            errKinds[msg]++;
        }
        for (const auto& kv : errKinds)
            std::printf("  GMSH_ERR %4d x %s\n", kv.second, kv.first.c_str());
        gmsh::logger::stop();

        if (!mshOut.empty()) {
            gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
            gmsh::write(mshOut);
        }
        const double meshSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - tm).count();

        const bool ok = tets > 0 && inverted == 0 && meshedVols.size() == vout.size() &&
                        std::fabs(volDev) < 1e-3;
        double unmeshedVol = 0;
        for (int v : emptyVols) { double m = 0; gmsh::model::occ::getMass(3, v, m); unmeshedVol += m; }
        std::printf("STEPCHECK %s MESH %s volsIn=%zu volsOut=%zu meshedVols=%zu unmeshedVol=%.4fmm3 "
                    "gmshErrs=%d tets=%zu minSICN=%.4f inverted=%zu volDev=%+.4f%% "
                    "szMin=%.4gmm szMax=%.4gmm t=%.1fs\n",
                    path.c_str(), ok ? "MESHABLE" : "NOT_MESHABLE", nIn, vout.size(),
                    meshedVols.size(), unmeshedVol * 1e9, nErr, tets, minQ, inverted,
                    volDev * 100.0, szMin * 1e3, szMax * 1e3, meshSec);
        if (!ok) rc = 1;
        gmsh::finalize();
    } catch (const std::exception& e) {
        std::printf("STEPCHECK %s MESH EXCEPTION %s\n", path.c_str(), e.what());
        try { gmsh::finalize(); } catch (...) {}
        return 1;
    }
    return rc;
}
