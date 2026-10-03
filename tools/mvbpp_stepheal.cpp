// mvbpp_stepheal <in.step> [--diag] [--out out.step] [heal ops...]
// Moved verbatim from OMFEM tools/omfem_stepheal.cpp at d92be53 (ABT #1588, step 8); renamed. No logic changed.
//
// WHY. The corpus audit (2026-08-24) found every produced STEP watertight yet only 7/30
// conformally meshable. The mesher's complaints -- "Impossible to mesh periodic surface",
// "The 1D mesh seems not to be forming a closed loop", "Invalid boundary mesh (overlapping
// facets)" -- are SYMPTOMS reported against a gmsh entity tag, which says nothing about the
// CAD feature that caused them. This tool measures the CAD features those symptoms come from,
// and can apply the standard OCCT repairs so a proposed fix can be MEASURED (heal -> re-run
// omfem_stepcheck --mesh) instead of argued.
//
// --diag reports, per file and per named solid:
//   * SEAM faces: a face whose surface closes on itself (full-revolution cylinder/torus/
//     BSpline pipe) carries a seam edge. gmsh cannot mesh a periodic face -- this is the
//     direct cause of the "Impossible to mesh periodic surface" class.
//   * MICRO edges / SLIVER faces: sub-tolerance features. A 1D loop that "seems not to be
//     forming a closed loop" is normally a wire whose tiny edge fell below the mesh size and
//     collapsed, breaking the loop.
//   * surface-type histogram, so a class can be attributed to a primitive kind.
//
// Heal ops (applied in this order, each optional, names preserved through XCAF):
//   --divide-closed   ShapeUpgrade_ShapeDivideClosed: split every closed/periodic face into
//                     non-periodic pieces. The canonical fix for the periodic class.
//   --fix <tol>       ShapeFix_Shape: the general repair (wires, orientation, small edges).
//   --sew <tol>       BRepBuilderAPI_Sewing: re-sew faces into shells at a wider tolerance.
//   --unify           ShapeUpgrade_UnifySameDomain: merge co-surface faces/edges, removing
//                     the imprint clutter a partial weld leaves behind.
//   --nurbs           BRepBuilderAPI_NurbsConvert: re-express every analytic surface as a
//                     B-spline. A conic IS exactly representable as a rational B-spline, so
//                     this moves NO point -- but a cylinder/cone/torus face trimmed to less
//                     than a full turn loses the IsUPeriodic flag that gmsh dispatches on
//                     (OCCFace.cpp:142), taking it off the fragile periodic mesher.
// Exit 0 on success, 2 on usage/read/write error.
#include <STEPCAFControl_Reader.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDataStd_Name.hxx>
#include <TCollection_ExtendedString.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Status.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <ShapeUpgrade_ShapeDivideClosed.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <ShapeFix_Shape.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace {

const char* surfName(GeomAbs_SurfaceType t) {
    switch (t) {
        case GeomAbs_Plane:              return "plane";
        case GeomAbs_Cylinder:           return "cylinder";
        case GeomAbs_Cone:               return "cone";
        case GeomAbs_Sphere:             return "sphere";
        case GeomAbs_Torus:              return "torus";
        case GeomAbs_BezierSurface:      return "bezier";
        case GeomAbs_BSplineSurface:     return "bspline";
        case GeomAbs_SurfaceOfRevolution:return "revolution";
        case GeomAbs_SurfaceOfExtrusion: return "extrusion";
        case GeomAbs_OffsetSurface:      return "offset";
        default:                         return "other";
    }
}

struct Stats {
    int faces = 0, edges = 0;
    int seamFaces = 0;          // face carrying a seam edge -> periodic, gmsh cannot mesh it
    int periodicSurf = 0;       // surface declares U or V periodicity (seam or not)
    int wrapFaces = 0;          // face whose parametric range spans the FULL period: a 360-degree
                                // wrap with no seam split. THIS is what gmsh refuses with
                                // "Impossible to mesh periodic surface" -- a cylinder face can be
                                // periodic-by-surface yet only span 90 degrees, which meshes fine.
    int microEdges = 0;         // shorter than edgeTol
    int sliverFaces = 0;        // area below faceTol
    // SHAPE TOLERANCE. OCC stores a per-vertex/edge/face tolerance and every boolean works to
    // the WIDEST tolerance it meets. A shape written at 1e-7 m but carrying 1e-3 mm edges makes
    // the intersector search a fat band instead of a curve -- slow, and prone to finding
    // "intersections" that are not there. Worth reporting because nothing else here does:
    // BRepCheck only flags InvalidToleranceValue (a violated invariant), not a merely INFLATED
    // one, so a shape can be VALID, WATERTIGHT and still be dear to intersect.
    double maxVertTol = 0.0, maxEdgeTol = 0.0, maxFaceTol = 0.0;
    double sumEdgeTol = 0.0;
    int fatEdges = 0;           // edge tolerance above 1e-6 mm
    int fatFaces = 0;
    std::map<std::string, int> byType;        // surface type histogram
    std::map<std::string, int> seamByType;    // seam faces by surface type
    std::map<std::string, int> wrapByType;    // full-period wrapped faces by surface type
};

double edgeLength(const TopoDS_Edge& e) {
    if (BRep_Tool::Degenerated(e)) return 0.0;
    BRepAdaptor_Curve c(e);
    return GCPnts_AbscissaPoint::Length(c);
}

double faceArea(const TopoDS_Face& f) {
    GProp_GProps p;
    BRepGProp::SurfaceProperties(f, p);
    return p.Mass();
}

void collect(const TopoDS_Shape& sh, Stats& s, double edgeTol, double faceTol) {
    TopTools_IndexedMapOfShape fm, em;
    TopExp::MapShapes(sh, TopAbs_FACE, fm);
    TopExp::MapShapes(sh, TopAbs_EDGE, em);
    s.faces += fm.Extent();
    s.edges += em.Extent();
    // Tolerance sweep (mm; the STEP declares millimetres and OCCT reads the file's own unit).
    {
        TopTools_IndexedMapOfShape vm;
        TopExp::MapShapes(sh, TopAbs_VERTEX, vm);
        for (int i = 1; i <= vm.Extent(); ++i)
            s.maxVertTol = std::max(s.maxVertTol,
                                    (double)BRep_Tool::Tolerance(TopoDS::Vertex(vm(i))));
        for (int i = 1; i <= em.Extent(); ++i) {
            const double t = BRep_Tool::Tolerance(TopoDS::Edge(em(i)));
            s.maxEdgeTol = std::max(s.maxEdgeTol, t);
            s.sumEdgeTol += t;
            if (t > 1e-6) ++s.fatEdges;
        }
        for (int i = 1; i <= fm.Extent(); ++i) {
            const double t = BRep_Tool::Tolerance(TopoDS::Face(fm(i)));
            s.maxFaceTol = std::max(s.maxFaceTol, t);
            if (t > 1e-6) ++s.fatFaces;
        }
    }
    for (int i = 1; i <= fm.Extent(); ++i) {
        const TopoDS_Face& f = TopoDS::Face(fm(i));
        BRepAdaptor_Surface bs(f, Standard_False);
        const std::string tn = surfName(bs.GetType());
        s.byType[tn]++;
        Handle(Geom_Surface) gs = BRep_Tool::Surface(f);
        const bool per = !gs.IsNull() && (gs->IsUPeriodic() || gs->IsVPeriodic());
        if (per) ++s.periodicSurf;
        bool seam = false;
        for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
            if (BRepTools::IsReallyClosed(TopoDS::Edge(ex.Current()), f)) { seam = true; break; }
        }
        if (seam) { ++s.seamFaces; s.seamByType[tn]++; }
        if (!gs.IsNull()) {
            double u0, u1, v0, v1;
            BRepTools::UVBounds(f, u0, u1, v0, v1);
            const bool uWrap = gs->IsUPeriodic() && (u1 - u0) >= gs->UPeriod() * (1.0 - 1e-6);
            const bool vWrap = gs->IsVPeriodic() && (v1 - v0) >= gs->VPeriod() * (1.0 - 1e-6);
            if (uWrap || vWrap) { ++s.wrapFaces; s.wrapByType[tn]++; }
        }
        if (faceArea(f) < faceTol) {
            ++s.sliverFaces;
            // Name the sliver: without a location nobody can find which CONTACT produced it,
            // and slivers are exactly what kills gmsh edge recovery ("Unable to recover the
            // edge"): a face thinner than the mesh size whose two long edges the 2D mesher
            // cannot keep apart.
            Bnd_Box fb;
            BRepBndLib::Add(f, fb);
            if (!fb.IsVoid()) {
                double x0,y0,z0,x1,y1,z1; fb.Get(x0,y0,z0,x1,y1,z1);
                std::printf("  SLIVER area=%.3e mm2 %s bbox x[%.4f,%.4f] y[%.4f,%.4f] z[%.4f,%.4f]\n",
                            faceArea(f), tn.c_str(), x0,x1,y0,y1,z0,z1);
            }
        }
    }
    for (int i = 1; i <= em.Extent(); ++i) {
        const TopoDS_Edge& e = TopoDS::Edge(em(i));
        if (BRep_Tool::Degenerated(e)) continue;
        if (edgeLength(e) < edgeTol) {
            ++s.microEdges;
            // Locate the first few: 34 anonymous micro-edges in 12_boost's fragment were
            // undiagnosable without their positions.
            if (s.microEdges <= 6) {
                Bnd_Box bb;
                BRepBndLib::Add(e, bb);
                if (!bb.IsVoid()) {
                    double x0, y0, z0, x1, y1, z1;
                    bb.Get(x0, y0, z0, x1, y1, z1);
                    std::printf("  MICROEDGE len=%.3e mm at x[%.4f,%.4f] y[%.4f,%.4f] z[%.4f,%.4f]\n",
                                edgeLength(e), x0, x1, y0, y1, z0, z1);
                }
            }
        }
    }
}

bool hasArg(int argc, char** argv, const char* f) {
    for (int i = 2; i < argc; ++i) if (!std::strcmp(argv[i], f)) return true;
    return false;
}
double argVal(int argc, char** argv, const char* f, double d) {
    for (int i = 2; i + 1 < argc; ++i) if (!std::strcmp(argv[i], f)) return std::atof(argv[i + 1]);
    return d;
}
const char* argStr(int argc, char** argv, const char* f, const char* d) {
    for (int i = 2; i + 1 < argc; ++i) if (!std::strcmp(argv[i], f)) return argv[i + 1];
    return d;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: %s <in.step> [--diag] [--per-solid] [--out out.step]\n"
            "         [--divide-closed] [--fix <tol>] [--sew <tol>] [--unify]\n"
            "         [--edge-tol <mm>] [--face-tol <mm2>]\n", argv[0]);
        return 2;
    }
    const std::string in = argv[1];
    const bool diag     = hasArg(argc, argv, "--diag");
    const bool perSolid = hasArg(argc, argv, "--per-solid");
    const std::string out = argStr(argc, argv, "--out", "");
    const bool divClosed = hasArg(argc, argv, "--divide-closed");
    const bool unify     = hasArg(argc, argv, "--unify");
    const bool nurbs     = hasArg(argc, argv, "--nurbs");
    const double fixTol  = argVal(argc, argv, "--fix", -1.0);
    const double sewTol  = argVal(argc, argv, "--sew", -1.0);
    // STEP declares millimetres, and OCCT reads the file's own unit, so these are mm / mm2.
    const double edgeTol = argVal(argc, argv, "--edge-tol", 1e-3);
    const double faceTol = argVal(argc, argv, "--face-tol", 1e-6);

    Handle(TDocStd_Document) doc;
    XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
    STEPCAFControl_Reader reader;
    reader.SetNameMode(true);
    if (reader.ReadFile(in.c_str()) != IFSelect_RetDone || !reader.Transfer(doc)) {
        std::printf("STEPHEAL %s READ_FAIL\n", in.c_str());
        return 2;
    }
    Handle(XCAFDoc_ShapeTool) tool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
    TDF_LabelSequence labels;
    tool->GetFreeShapes(labels);

    struct Part { std::string name; TopoDS_Shape shape; };
    std::vector<Part> parts;
    for (Standard_Integer i = 1; i <= labels.Length(); ++i) {
        Handle(TDataStd_Name) nm;
        std::string name = "shape" + std::to_string(i);
        if (labels.Value(i).FindAttribute(TDataStd_Name::GetID(), nm))
            name = std::string(TCollection_AsciiString(nm->Get()).ToCString());
        TopoDS_Shape sh = tool->GetShape(labels.Value(i));
        if (!sh.IsNull()) parts.push_back({name, sh});
    }
    if (parts.empty()) { std::printf("STEPHEAL %s NO_SHAPES\n", in.c_str()); return 2; }

    // --why: for every BRepCheck-INVALID solid, name the offending sub-shapes and their exact
    // BRepCheck status. "invalid=1" alone cannot be acted on; the status says whether it is a
    // self-intersection, a bad orientation, a tolerance violation, a non-closed shell, ...
    if (hasArg(argc, argv, "--why")) {
        auto statusName = [](BRepCheck_Status s) -> const char* {
            switch (s) {
                case BRepCheck_NoError: return "NoError";
                case BRepCheck_InvalidPointOnCurve: return "InvalidPointOnCurve";
                case BRepCheck_InvalidPointOnCurveOnSurface: return "InvalidPointOnCurveOnSurface";
                case BRepCheck_InvalidPointOnSurface: return "InvalidPointOnSurface";
                case BRepCheck_No3DCurve: return "No3DCurve";
                case BRepCheck_Multiple3DCurve: return "Multiple3DCurve";
                case BRepCheck_Invalid3DCurve: return "Invalid3DCurve";
                case BRepCheck_NoCurveOnSurface: return "NoCurveOnSurface";
                case BRepCheck_InvalidCurveOnSurface: return "InvalidCurveOnSurface";
                case BRepCheck_InvalidCurveOnClosedSurface: return "InvalidCurveOnClosedSurface";
                case BRepCheck_InvalidSameRangeFlag: return "InvalidSameRangeFlag";
                case BRepCheck_InvalidSameParameterFlag: return "InvalidSameParameterFlag";
                case BRepCheck_InvalidDegeneratedFlag: return "InvalidDegeneratedFlag";
                case BRepCheck_FreeEdge: return "FreeEdge";
                case BRepCheck_InvalidMultiConnexity: return "InvalidMultiConnexity";
                case BRepCheck_InvalidRange: return "InvalidRange";
                case BRepCheck_EmptyWire: return "EmptyWire";
                case BRepCheck_RedundantEdge: return "RedundantEdge";
                case BRepCheck_SelfIntersectingWire: return "SelfIntersectingWire";
                case BRepCheck_NoSurface: return "NoSurface";
                case BRepCheck_InvalidWire: return "InvalidWire";
                case BRepCheck_RedundantWire: return "RedundantWire";
                case BRepCheck_IntersectingWires: return "IntersectingWires";
                case BRepCheck_InvalidImbricationOfWires: return "InvalidImbricationOfWires";
                case BRepCheck_EmptyShell: return "EmptyShell";
                case BRepCheck_RedundantFace: return "RedundantFace";
                case BRepCheck_UnorientableShape: return "UnorientableShape";
                case BRepCheck_NotClosed: return "NotClosed";
                case BRepCheck_NotConnected: return "NotConnected";
                case BRepCheck_SubshapeNotInShape: return "SubshapeNotInShape";
                case BRepCheck_BadOrientation: return "BadOrientation";
                case BRepCheck_BadOrientationOfSubshape: return "BadOrientationOfSubshape";
                case BRepCheck_InvalidPolygonOnTriangulation: return "InvalidPolygonOnTriangulation";
                case BRepCheck_InvalidToleranceValue: return "InvalidToleranceValue";
                case BRepCheck_CheckFail: return "CheckFail";
                default: return "Unknown";
            }
        };
        auto kindName = [](TopAbs_ShapeEnum t) -> const char* {
            switch (t) {
                case TopAbs_VERTEX: return "vertex";
                case TopAbs_EDGE: return "edge";
                case TopAbs_WIRE: return "wire";
                case TopAbs_FACE: return "face";
                case TopAbs_SHELL: return "shell";
                case TopAbs_SOLID: return "solid";
                default: return "shape";
            }
        };
        for (const auto& pt : parts) {
            int k = 0;
            for (TopExp_Explorer ex(pt.shape, TopAbs_SOLID); ex.More(); ex.Next(), ++k) {
                BRepCheck_Analyzer ana(ex.Current());
                if (ana.IsValid()) continue;
                std::printf("  WHY %s [solid %d] is INVALID:\n", pt.name.c_str(), k);
                int shown = 0;
                const TopAbs_ShapeEnum kinds[] = {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_WIRE,
                                                  TopAbs_FACE, TopAbs_SHELL};
                for (TopAbs_ShapeEnum kind : kinds) {
                    for (TopExp_Explorer sx(ex.Current(), kind); sx.More() && shown < 20; sx.Next()) {
                        Handle(BRepCheck_Result) res;
                        try { res = ana.Result(sx.Current()); } catch (...) { continue; }
                        if (res.IsNull()) continue;
                        for (BRepCheck_ListIteratorOfListOfStatus it(res->Status()); it.More();
                             it.Next()) {
                            if (it.Value() == BRepCheck_NoError) continue;
                            Bnd_Box bb;
                            BRepBndLib::Add(sx.Current(), bb);
                            double x0 = 0, y0 = 0, z0 = 0, x1 = 0, y1 = 0, z1 = 0;
                            if (!bb.IsVoid()) bb.Get(x0, y0, z0, x1, y1, z1);
                            std::printf("    %-6s %-32s at (%.4f, %.4f, %.4f)..(%.4f, %.4f, %.4f)\n",
                                        kindName(kind), statusName(it.Value()),
                                        x0, y0, z0, x1, y1, z1);
                            ++shown;
                            break;
                        }
                    }
                }
                if (shown == 0) std::printf("    (analyzer reports no sub-shape status)\n");
            }
        }
    }

    // --selfint: OCC's own self-interference verdict per named solid (BOPAlgo_ArgumentAnalyzer,
    // the check BOP runs on its arguments). The 04_forward welds are BRepCheck-VALID and
    // watertight yet unmeshable -- this is the test that can tell a sound weld from one that
    // carries internal face-face intersections, so it is the candidate ACCEPTANCE test for
    // "always weld when possible" (Alf, 2026-08-26).
    if (hasArg(argc, argv, "--selfint")) {
        int bad = 0;
        for (const auto& pt : parts) {
            int k = 0;
            for (TopExp_Explorer ex(pt.shape, TopAbs_SOLID); ex.More(); ex.Next(), ++k) {
                BOPAlgo_ArgumentAnalyzer an;
                an.SetShape1(ex.Current());
                an.ArgumentTypeMode() = Standard_True;
                an.SelfInterMode() = Standard_True;
                an.SmallEdgeMode() = Standard_True;
                an.Perform();
                if (an.HasFaulty()) {
                    ++bad;
                    int selfInt = 0, smallEdge = 0, other = 0;
                    int shown = 0;
                    for (const auto& r : an.GetCheckResult()) {
                        if (r.GetCheckStatus() == BOPAlgo_SelfIntersect) ++selfInt;
                        else if (r.GetCheckStatus() == BOPAlgo_TooSmallEdge) ++smallEdge;
                        else ++other;
                        // Name the first few faulty regions: 1228 anonymous self-intersections
                        // on 04_forward's primary were unlocatable without this.
                        if (r.GetCheckStatus() == BOPAlgo_SelfIntersect && shown < 8) {
                            for (const auto& fs : r.GetFaultyShapes1()) {
                                Bnd_Box bb;
                                BRepBndLib::Add(fs, bb);
                                if (bb.IsVoid()) continue;
                                double x0, y0, z0, x1, y1, z1;
                                bb.Get(x0, y0, z0, x1, y1, z1);
                                std::printf("    faulty %s bbox x[%.4f,%.4f] y[%.4f,%.4f] z[%.4f,%.4f] mm\n",
                                            fs.ShapeType() == TopAbs_FACE ? "face" :
                                            fs.ShapeType() == TopAbs_EDGE ? "edge" : "shape",
                                            x0 * 1e3, x1 * 1e3, y0 * 1e3, y1 * 1e3, z0 * 1e3, z1 * 1e3);
                                if (++shown >= 8) break;
                            }
                        }
                    }
                    std::printf("  SELFINT %s%s: selfIntersections=%d smallEdges=%d other=%d\n",
                                pt.name.c_str(),
                                k ? (" [solid " + std::to_string(k) + "]").c_str() : "",
                                selfInt, smallEdge, other);
                }
            }
        }
        std::printf("STEPHEAL %s SELFINT %s (%d faulty solid(s))\n", in.c_str(),
                    bad ? "FAULTY" : "CLEAN", bad);
    }

    if (diag) {
        Stats all;
        for (const auto& p : parts) {
            Stats s;
            collect(p.shape, s, edgeTol, faceTol);
            if (perSolid && (s.seamFaces || s.wrapFaces || s.microEdges || s.sliverFaces))
                std::printf("  PART %-44s faces=%d seam=%d wrap=%d micro=%d sliver=%d\n",
                            p.name.c_str(), s.faces, s.seamFaces, s.wrapFaces, s.microEdges, s.sliverFaces);
            all.faces += s.faces; all.edges += s.edges;
            all.seamFaces += s.seamFaces; all.periodicSurf += s.periodicSurf;
            all.wrapFaces += s.wrapFaces;
            all.microEdges += s.microEdges; all.sliverFaces += s.sliverFaces;
            all.maxVertTol = std::max(all.maxVertTol, s.maxVertTol);
            all.maxEdgeTol = std::max(all.maxEdgeTol, s.maxEdgeTol);
            all.maxFaceTol = std::max(all.maxFaceTol, s.maxFaceTol);
            all.sumEdgeTol += s.sumEdgeTol;
            all.fatEdges += s.fatEdges; all.fatFaces += s.fatFaces;
            for (const auto& kv : s.byType) all.byType[kv.first] += kv.second;
            for (const auto& kv : s.seamByType) all.seamByType[kv.first] += kv.second;
            for (const auto& kv : s.wrapByType) all.wrapByType[kv.first] += kv.second;
        }
        std::string hist, seamHist;
        for (const auto& kv : all.byType)
            hist += (hist.empty() ? "" : ",") + kv.first + ":" + std::to_string(kv.second);
        for (const auto& kv : all.wrapByType)
            seamHist += (seamHist.empty() ? "" : ",") + kv.first + ":" + std::to_string(kv.second);
        std::printf("STEPHEAL %s DIAG parts=%zu faces=%d edges=%d seamFaces=%d periodicSurf=%d "
                    "wrapFaces=%d microEdges(<%.4gmm)=%d sliverFaces(<%.4gmm2)=%d types=[%s] wrapTypes=[%s]\n",
                    in.c_str(), parts.size(), all.faces, all.edges, all.seamFaces,
                    all.periodicSurf, all.wrapFaces, edgeTol, all.microEdges, faceTol, all.sliverFaces,
                    hist.c_str(), seamHist.empty() ? "-" : seamHist.c_str());
        // Tolerance is reported separately because it is NOT a validity question: a shape can be
        // VALID and WATERTIGHT and still be dear to intersect if its edges carry a fat tolerance.
        std::printf("STEPHEAL %s TOL maxVert=%.3gmm maxEdge=%.3gmm maxFace=%.3gmm "
                    "meanEdge=%.3gmm fatEdges(>1e-6mm)=%d/%d fatFaces=%d/%d\n",
                    in.c_str(), all.maxVertTol, all.maxEdgeTol, all.maxFaceTol,
                    all.edges ? all.sumEdgeTol / all.edges : 0.0,
                    all.fatEdges, all.edges, all.fatFaces, all.faces);
    }

    if (out.empty()) return 0;

    // ---- healing, name-preserving -------------------------------------------------------
    Handle(TDocStd_Document) odoc;
    XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", odoc);
    Handle(XCAFDoc_ShapeTool) otool = XCAFDoc_DocumentTool::ShapeTool(odoc->Main());
    int healed = 0, brokeValidity = 0;
    for (const auto& p : parts) {
        TopoDS_Shape s = p.shape;
        const bool wasValid = BRepCheck_Analyzer(s).IsValid();
        try {
            if (divClosed) {
                ShapeUpgrade_ShapeDivideClosed d(s);
                d.SetNbSplitPoints(1);
                if (d.Perform() && !d.Result().IsNull()) s = d.Result();
            }
            if (fixTol > 0) {
                Handle(ShapeFix_Shape) sf = new ShapeFix_Shape(s);
                sf->SetPrecision(fixTol);
                sf->SetMaxTolerance(fixTol * 10);
                if (sf->Perform() && !sf->Shape().IsNull()) s = sf->Shape();
            }
            if (sewTol > 0) {
                BRepBuilderAPI_Sewing sew(sewTol);
                sew.Add(s);
                sew.Perform();
                if (!sew.SewedShape().IsNull()) s = sew.SewedShape();
            }
            if (nurbs) {
                BRepBuilderAPI_NurbsConvert nc(s, /*Copy*/ Standard_True);
                if (!nc.Shape().IsNull()) s = nc.Shape();
            }
            if (unify) {
                ShapeUpgrade_UnifySameDomain u(s, true, true, false);
                u.Build();
                if (!u.Shape().IsNull()) s = u.Shape();
            }
        } catch (const std::exception& e) {
            std::printf("  HEAL_FAIL %-40s %s\n", p.name.c_str(), e.what());
            s = p.shape;
        }
        if (wasValid && !BRepCheck_Analyzer(s).IsValid()) {
            // Never ship a repair that breaks a solid that was already valid.
            std::printf("  HEAL_REVERTED %-40s (repair made a valid solid invalid)\n", p.name.c_str());
            s = p.shape;
            ++brokeValidity;
        } else if (!s.IsEqual(p.shape)) {
            ++healed;
        }
        TDF_Label lab = otool->AddShape(s, false);
        TDataStd_Name::Set(lab, TCollection_ExtendedString(p.name.c_str()));
    }
    STEPCAFControl_Writer writer;
    writer.SetNameMode(true);
    if (!writer.Transfer(odoc, STEPControl_AsIs) || writer.Write(out.c_str()) != IFSelect_RetDone) {
        std::printf("STEPHEAL %s WRITE_FAIL %s\n", in.c_str(), out.c_str());
        return 2;
    }
    std::printf("STEPHEAL %s HEALED parts=%zu changed=%d reverted=%d -> %s\n",
                in.c_str(), parts.size(), healed, brokeValidity, out.c_str());
    return 0;
}
