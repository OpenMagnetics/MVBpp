// mvbpp_faceknots <file.step> [topN]
// Moved verbatim from OMFEM tools/omfem_faceknots.cpp at c12864b (ABT #1588, step 8); renamed. No logic changed.
// Which faces will make OCCT's booleans explode? BOPAlgo_PaveFiller::PerformEF projects edge points
// onto faces with GeomAPI_ProjectPointOnSurf -> Extrema_GenExtPS::BuildGrid, whose sample grid scales
// with the B-spline's knot counts. A faceted sweep (--segments N) turns every lateral strip into a
// B-spline; a strip with thousands of knots is a projection that never returns (single_switch at 12
// segments: 28 min in BuildGrid) or an Array2 that overflows (14_dab: SIGSEGV in Resize). This lists
// the worst faces by knots (U x V), with their solid, so the sweep that made them can be found.
#include <STEPCAFControl_Reader.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDataStd_Name.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Trsf.hxx>
#include <BOPAlgo_ListOfCheckResult.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopExp.hxx>
#include <BRepTools.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <cmath>
#include <TCollection_AsciiString.hxx>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

struct Rec { std::string solid; int faceIdx; std::string type; int uk, vk, up, vp, ud, vd; double bb[6]; long long score; double u0, u1, v0, v1, area; bool wild; };

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s <file.step> [topN]\n", argv[0]); return 2; }
    const int topN = argc > 2 ? std::atoi(argv[2]) : 15;
    // FACEKNOTS_SOLID=<name-substring>:<A>-<B> lists EVERY face of solids A..B of the parts whose
    // name contains the substring (file order, no sorting) and runs BOPAlgo's self-intersection
    // analyzer on each, naming the faulty faces -- the per-face view of a stepcheck DEFECT line.
    std::string filtName; int filtA = -1, filtB = -1;
    if (const char* fs = std::getenv("FACEKNOTS_SOLID")) {
        std::string v(fs); auto c = v.rfind(':'); std::string rng = v;
        if (c != std::string::npos) { filtName = v.substr(0, c); rng = v.substr(c + 1); }
        if (std::sscanf(rng.c_str(), "%d-%d", &filtA, &filtB) == 1) filtB = filtA;
    }
    Handle(TDocStd_Document) doc; XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
    STEPCAFControl_Reader rd; rd.SetNameMode(true);
    if (rd.ReadFile(argv[1]) != IFSelect_RetDone || !rd.Transfer(doc)) { std::fprintf(stderr, "read failed\n"); return 2; }
    Handle(XCAFDoc_ShapeTool) st = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
    TDF_LabelSequence free; st->GetFreeShapes(free);
    std::vector<Rec> recs; std::map<std::string, int> typeCount; long long totalFaces = 0;
    for (int i = 1; i <= free.Length(); ++i) {
        TopoDS_Shape sh = st->GetShape(free.Value(i));
        // FACEKNOTS_SCALE=<f> rescales the geometry before analysis (poles move, knots stay), to
        // test whether a fault is a parametrisation-vs-geometry scale artefact of the m->mm export.
        if (const char* sc = std::getenv("FACEKNOTS_SCALE")) {
            gp_Trsf t; t.SetScale(gp_Pnt(0, 0, 0), std::atof(sc));
            sh = BRepBuilderAPI_Transform(sh, t, Standard_True).Shape();
        }
        std::string nm = "?"; Handle(TDataStd_Name) n;
        if (free.Value(i).FindAttribute(TDataStd_Name::GetID(), n)) nm = TCollection_AsciiString(n->Get()).ToCString();
        int sidx = 0;
        for (TopExp_Explorer se(sh, TopAbs_SOLID); se.More(); se.Next(), ++sidx) {
            const bool inFilter = filtA >= 0 && sidx >= filtA && sidx <= filtB &&
                                  (filtName.empty() || nm.find(filtName) != std::string::npos);
            if (filtA >= 0 && !inFilter) continue;
            if (inFilter) {
                TopTools_IndexedMapOfShape fmap; TopExp::MapShapes(se.Current(), TopAbs_FACE, fmap);
                BOPAlgo_ArgumentAnalyzer an; an.SetShape1(se.Current());
                an.SelfInterMode() = Standard_True; an.SmallEdgeMode() = Standard_True;
                an.ArgumentTypeMode() = Standard_False; an.StopOnFirstFaulty() = Standard_False;
                an.Perform();
                std::printf("== %s [solid %d]: %d faces, BOPAlgo %s\n", nm.c_str(), sidx, fmap.Extent(), an.HasFaulty() ? "FAULTY" : "clean");
                for (BOPAlgo_ListIteratorOfListOfCheckResult it(an.GetCheckResult()); it.More(); it.Next()) {
                    const BOPAlgo_CheckResult& cr = it.Value();
                    std::printf("   status %d faces:", (int)cr.GetCheckStatus());
                    for (TopTools_ListIteratorOfListOfShape fi(cr.GetFaultyShapes1()); fi.More(); fi.Next()) {
                        const TopoDS_Shape& fx = fi.Value();
                        if (fx.ShapeType() == TopAbs_FACE) std::printf(" F%d", fmap.FindIndex(fx) - 1);
                        else std::printf(" (%d)", (int)fx.ShapeType());
                    }
                    std::printf("\n");
                }
            }
            int fidx = 0;
            for (TopExp_Explorer fe(se.Current(), TopAbs_FACE); fe.More(); fe.Next(), ++fidx) {
                ++totalFaces;
                const TopoDS_Face f = TopoDS::Face(fe.Current());
                BRepAdaptor_Surface ad(f, false);
                Rec r; r.solid = nm + " [solid " + std::to_string(sidx) + "]"; r.faceIdx = fidx;
                r.uk = r.vk = r.up = r.vp = r.ud = r.vd = 0; r.score = 1;
                switch (ad.GetType()) {
                    case GeomAbs_Plane: r.type = "plane"; break;
                    case GeomAbs_Cylinder: r.type = "cylinder"; break;
                    case GeomAbs_Torus: r.type = "torus"; break;
                    case GeomAbs_Cone: r.type = "cone"; break;
                    case GeomAbs_Sphere: r.type = "sphere"; break;
                    case GeomAbs_BSplineSurface: {
                        r.type = "bspline";
                        Handle(Geom_BSplineSurface) bs = ad.BSpline();
                        r.uk = bs->NbUKnots(); r.vk = bs->NbVKnots(); r.up = bs->NbUPoles(); r.vp = bs->NbVPoles();
                        r.ud = bs->UDegree(); r.vd = bs->VDegree();
                        // Extrema's grid grows like (UDegree+1)*NbUKnots x (VDegree+1)*NbVKnots
                        r.score = (long long)((r.ud + 1) * r.uk) * (long long)((r.vd + 1) * r.vk);
                        break; }
                    case GeomAbs_BezierSurface: r.type = "bezier"; break;
                    case GeomAbs_SurfaceOfExtrusion: r.type = "extrusion"; break;
                    case GeomAbs_SurfaceOfRevolution: r.type = "revolution"; break;
                    default: r.type = "other"; break;
                }
                typeCount[r.type]++;
                // PARAMETRISATION SANITY: Extrema_GenExtPS sizes its grid from the face's UV range; a
                // face with a wild or non-finite range is the crash in Array2::Resize (2026-09-06,
                // single_switch / 14_dab at --segments 12, inside three different booleans).
                BRepTools::UVBounds(f, r.u0, r.u1, r.v0, r.v1);
                GProp_GProps sp; BRepGProp::SurfaceProperties(f, sp); r.area = sp.Mass();
                const double span = std::max(std::fabs(r.u1 - r.u0), std::fabs(r.v1 - r.v0));
                r.wild = !std::isfinite(r.u0) || !std::isfinite(r.u1) || !std::isfinite(r.v0) || !std::isfinite(r.v1) ||
                         span > 1e4 || !(r.area > 0.0) || !std::isfinite(r.area);
                if (inFilter) { Bnd_Box b; BRepBndLib::Add(f, b); b.Get(r.bb[0], r.bb[1], r.bb[2], r.bb[3], r.bb[4], r.bb[5]); recs.push_back(r); continue; }
                if (r.type == "bspline" || r.wild) {
                    Bnd_Box b; BRepBndLib::Add(f, b); b.Get(r.bb[0], r.bb[1], r.bb[2], r.bb[3], r.bb[4], r.bb[5]);
                    if (r.wild) r.score = (long long)1e15;      // sort to the top
                    recs.push_back(r);
                }
            }
        }
    }
    std::printf("faces=%lld  by type:", totalFaces);
    for (auto& kv : typeCount) std::printf(" %s=%d", kv.first.c_str(), kv.second);
    std::printf("\n");
    if (filtA < 0) std::sort(recs.begin(), recs.end(), [](const Rec& a, const Rec& b) { return a.score > b.score; });
    long long over1e4 = 0, wild = 0; double minArea = 1e300; for (auto& r : recs) { if (r.score > 10000 && !r.wild) ++over1e4; if (r.wild) ++wild; minArea = std::min(minArea, r.area); }
    std::printf("listed faces=%zu  grid>1e4: %lld  WILD parametrisation/area: %lld  smallest face area %.3g mm2\n", recs.size(), over1e4, wild, minArea * 1e6);
    for (int i = 0; i < (filtA >= 0 ? (int)recs.size() : std::min<int>(topN, (int)recs.size())); ++i) {
        const Rec& r = recs[i];
        std::printf("  %s%-8s U[%.3g,%.3g] V[%.3g,%.3g] area=%.3g mm2 knots U=%d V=%d deg U=%d V=%d  %s face %d  x[%.3f,%.3f] y[%.3f,%.3f] z[%.3f,%.3f] mm\n",
                    r.wild ? "WILD " : "", r.type.c_str(), r.u0, r.u1, r.v0, r.v1, r.area * 1e6, r.uk, r.vk, r.ud, r.vd, r.solid.c_str(), r.faceIdx,
                    r.bb[0] * 1e3, r.bb[3] * 1e3, r.bb[1] * 1e3, r.bb[4] * 1e3, r.bb[2] * 1e3, r.bb[5] * 1e3);
    }
    return 0;
}
