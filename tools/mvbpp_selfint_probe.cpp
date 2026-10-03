// mvbpp_selfint_probe: WHICH faces of one solid cross each other.
// Moved verbatim from OMFEM tools/omfem_selfint_probe.cpp at 00740b3 (ABT #1588, step 8); renamed. No logic changed.
//
// omfem_stepcheck counts self-intersections per solid (BOPAlgo_ArgumentAnalyzer) and names the
// solid the way it labels it -- "<free shape name> [solid k]" -- but says nothing about where
// the defect is. This tool takes exactly that label and prints the solid's faces (surface type,
// area, extents) and the analyzer's faulty sub-shapes for every self-intersection it reports, so
// a defect class can be tied to a construction (a mitre cap, a swept lateral, a sphere elbow)
// instead of guessed at. Diagnostic only: it changes nothing and writes nothing.
//
//   mvbpp_selfint_probe <file.step> "<free shape name>" <k>
//
// k is the 0-based solid index inside that free shape, as in stepcheck's "[solid k]"; omit the
// suffix for k = 0.
#include <STEPCAFControl_Reader.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDataStd_Name.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopExp.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <GeomAbs_CurveType.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BOPAlgo_CheckResult.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Trsf.hxx>
#include <GeomLProp_SLProps.hxx>
#include <BRepTools.hxx>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace {

const char* surfaceName(GeomAbs_SurfaceType t) {
    switch (t) {
        case GeomAbs_Plane: return "plane";
        case GeomAbs_Cylinder: return "cylinder";
        case GeomAbs_Cone: return "cone";
        case GeomAbs_Sphere: return "sphere";
        case GeomAbs_Torus: return "torus";
        case GeomAbs_BezierSurface: return "bezier";
        case GeomAbs_BSplineSurface: return "bspline";
        case GeomAbs_SurfaceOfRevolution: return "revolution";
        case GeomAbs_SurfaceOfExtrusion: return "extrusion";
        case GeomAbs_OffsetSurface: return "offset";
        default: return "other";
    }
}

const char* curveName(GeomAbs_CurveType t) {
    switch (t) {
        case GeomAbs_Line: return "line";
        case GeomAbs_Circle: return "circle";
        case GeomAbs_Ellipse: return "ellipse";
        case GeomAbs_BSplineCurve: return "bspline";
        case GeomAbs_BezierCurve: return "bezier";
        default: return "other";
    }
}

const char* typeName(TopAbs_ShapeEnum t) {
    switch (t) {
        case TopAbs_FACE: return "FACE";
        case TopAbs_EDGE: return "EDGE";
        case TopAbs_VERTEX: return "VERTEX";
        case TopAbs_SOLID: return "SOLID";
        case TopAbs_SHELL: return "SHELL";
        case TopAbs_WIRE: return "WIRE";
        default: return "shape";
    }
}

void describe(const TopoDS_Shape& s, const TopTools_IndexedMapOfShape& faces,
              const TopTools_IndexedMapOfShape& edges, const char* indent) {
    Bnd_Box bb;
    BRepBndLib::Add(s, bb);
    double x0 = 0, y0 = 0, z0 = 0, x1 = 0, y1 = 0, z1 = 0;
    if (!bb.IsVoid()) bb.Get(x0, y0, z0, x1, y1, z1);
    if (s.ShapeType() == TopAbs_FACE) {
        BRepAdaptor_Surface surf(TopoDS::Face(s));
        GProp_GProps g;
        BRepGProp::SurfaceProperties(s, g);
        std::printf("%sFACE #%d %-10s area=%.6f mm2  ext=%.4f x %.4f x %.4f mm\n", indent,
                    faces.FindIndex(s), surfaceName(surf.GetType()), g.Mass(), x1 - x0, y1 - y0,
                    z1 - z0);
    }
    else if (s.ShapeType() == TopAbs_EDGE) {
        BRepAdaptor_Curve c(TopoDS::Edge(s));
        GProp_GProps g;
        BRepGProp::LinearProperties(s, g);
        std::printf("%sEDGE #%d %-10s len=%.6f mm  ext=%.4f x %.4f x %.4f mm\n", indent,
                    edges.FindIndex(s), curveName(c.GetType()), g.Mass(), x1 - x0, y1 - y0,
                    z1 - z0);
    }
    else {
        std::printf("%s%s  ext=%.4f x %.4f x %.4f mm\n", indent, typeName(s.ShapeType()),
                    x1 - x0, y1 - y0, z1 - z0);
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <file.step> \"<free shape name>\" [k]\n", argv[0]);
        return 2;
    }
    const std::string path = argv[1];
    const std::string want = argv[2];
    const int wantIndex = argc > 3 ? std::atoi(argv[3]) : 0;
    // PAIR MODE: a fourth argument names a second solid of the same free shape; the tool then
    // reports WHERE the two overlap (the boolean COMMON's volume and extents, in the mm frame
    // the booleans need), which tells an axial slab from a lateral wedge.
    // Same-shape pair: <step> "<A>" k j (argc 5). Cross-shape pair, for core-vs-copper:
    // <step> "<A>" k "<B>" j (argc 6).
    const bool crossPair = argc > 5;
    const int pairIndex = (!crossPair && argc > 4) ? std::atoi(argv[4]) : -1;
    const std::string wantB = crossPair ? argv[4] : want;
    const int pairIndexB = crossPair ? std::atoi(argv[5]) : pairIndex;

    Handle(TDocStd_Document) doc;
    XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
    STEPCAFControl_Reader reader;
    reader.SetNameMode(true);
    if (reader.ReadFile(path.c_str()) != IFSelect_RetDone || !reader.Transfer(doc)) {
        std::printf("READ_FAIL %s\n", path.c_str());
        return 2;
    }
    Handle(XCAFDoc_ShapeTool) tool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
    TDF_LabelSequence labels;
    tool->GetFreeShapes(labels);

    TopoDS_Shape target, other;
    for (Standard_Integer i = 1; i <= labels.Length(); ++i) {
        Handle(TDataStd_Name) nm;
        std::string name = "shape" + std::to_string(i);
        if (labels.Value(i).FindAttribute(TDataStd_Name::GetID(), nm))
            name = std::string(TCollection_AsciiString(nm->Get()).ToCString());
        if (name != want && name != wantB) continue;
        TopoDS_Shape sh = tool->GetShape(labels.Value(i));
        int k = 0;
        for (TopExp_Explorer ex(sh, TopAbs_SOLID); ex.More(); ex.Next(), ++k) {
            if (name == want && k == wantIndex) target = ex.Current();
            if (name == wantB && pairIndexB >= 0 && k == pairIndexB) other = ex.Current();
        }
    }
    if (pairIndexB >= 0 && (argc > 4)) {
        if (target.IsNull() || other.IsNull()) {
            std::printf("NOT FOUND: '%s' solid %d / '%s' solid %d\n", want.c_str(), wantIndex,
                        wantB.c_str(), pairIndexB);
            return 1;
        }
        auto report = [&](const char* tag, const TopoDS_Shape& sh) {
            Bnd_Box bb; BRepBndLib::Add(sh, bb);
            double x0, y0, z0, x1, y1, z1; bb.Get(x0, y0, z0, x1, y1, z1);
            GProp_GProps g; BRepGProp::VolumeProperties(sh, g);
            std::printf("%s vol=%.6f mm3  box x[%.4f,%.4f] y[%.4f,%.4f] z[%.4f,%.4f]\n", tag,
                        g.Mass(), x0, x1, y0, y1, z0, z1);
        };
        report("A", target);
        report("B", other);
        // Every PLANAR face of each solid: centroid and normal. Two pieces that should abut
        // share a cap plane; the axial distance between the matching caps is the overlap.
        auto caps = [&](const char* tag, const TopoDS_Shape& sh) {
            for (TopExp_Explorer fe(sh, TopAbs_FACE); fe.More(); fe.Next()) {
                const TopoDS_Face f = TopoDS::Face(fe.Current());
                BRepAdaptor_Surface surf(f);
                // Every face under 1.5 mm2 is a cap-sized face: planes AND the BSpline caps a
                // NURBS export pass or a pipe sweep leaves behind.
                GProp_GProps ga; BRepGProp::SurfaceProperties(f, ga);
                if (surf.GetType() != GeomAbs_Plane && ga.Mass() > 1.5) continue;
                GProp_GProps g; BRepGProp::SurfaceProperties(f, g);
                const gp_Pnt c = g.CentreOfMass();
                double u0, u1, v0, v1; BRepTools::UVBounds(f, u0, u1, v0, v1);
                GeomLProp_SLProps props(BRep_Tool::Surface(f), 0.5 * (u0 + u1), 0.5 * (v0 + v1), 1, 1e-9);
                gp_Dir n = props.IsNormalDefined() ? props.Normal() : gp_Dir(0, 0, 1);
                if (f.Orientation() == TopAbs_REVERSED) n.Reverse();
                std::printf("%s cap %-8s area=%.6f centre=(%.6f,%.6f,%.6f) normal=(%.6f,%.6f,%.6f)\n",
                            tag, surfaceName(surf.GetType()), g.Mass(), c.X(), c.Y(), c.Z(), n.X(), n.Y(), n.Z());
            }
        };
        caps("A", target);
        caps("B", other);
        // COMMON in the millimetre frame (the STEP is already in mm here).
        BRepAlgoAPI_Common common(target, other);
        if (!common.IsDone() || common.Shape().IsNull()) {
            std::printf("COMMON: boolean failed\n");
            return 1;
        }
        report("COMMON", common.Shape());
        int nf = 0;
        for (TopExp_Explorer fe(common.Shape(), TopAbs_FACE); fe.More(); fe.Next()) ++nf;
        std::printf("COMMON faces=%d\n", nf);
        TopTools_IndexedMapOfShape cf, ce;
        TopExp::MapShapes(common.Shape(), TopAbs_FACE, cf);
        TopExp::MapShapes(common.Shape(), TopAbs_EDGE, ce);
        for (int f = 1; f <= cf.Extent(); ++f) describe(cf(f), cf, ce, "  ");
        return 0;
    }
    if (target.IsNull()) {
        std::printf("NOT FOUND: free shape '%s' solid %d\n", want.c_str(), wantIndex);
        return 1;
    }

    TopTools_IndexedMapOfShape faces, edges;
    TopExp::MapShapes(target, TopAbs_FACE, faces);
    TopExp::MapShapes(target, TopAbs_EDGE, edges);
    GProp_GProps vg;
    BRepGProp::VolumeProperties(target, vg);
    std::printf("SOLID '%s' [solid %d]: vol=%.6f mm3 faces=%d edges=%d\n", want.c_str(), wantIndex,
                vg.Mass(), faces.Extent(), edges.Extent());
    for (int f = 1; f <= faces.Extent(); ++f) {
        describe(faces(f), faces, edges, "  ");
        // The face's own boundary, so a cap can be told from a lateral by its edges.
        for (TopExp_Explorer ee(faces(f), TopAbs_EDGE); ee.More(); ee.Next())
            describe(ee.Current(), faces, edges, "      ");
    }

    BOPAlgo_ArgumentAnalyzer an;
    an.SetShape1(target);
    an.ArgumentTypeMode() = Standard_True;
    an.SelfInterMode()    = Standard_True;
    an.SmallEdgeMode()    = Standard_True;
    an.Perform();
    if (!an.HasFaulty()) {
        std::printf("ANALYZER: no faults\n");
        return 0;
    }
    int n = 0;
    for (const auto& r : an.GetCheckResult()) {
        ++n;
        const char* what = r.GetCheckStatus() == BOPAlgo_SelfIntersect ? "SELF-INTERSECTION"
                         : r.GetCheckStatus() == BOPAlgo_TooSmallEdge ? "SMALL-EDGE"
                         : "OTHER";
        std::printf("FAULT %d: %s\n", n, what);
        std::printf("  shapes1:\n");
        for (TopTools_ListIteratorOfListOfShape it(r.GetFaultyShapes1()); it.More(); it.Next())
            describe(it.Value(), faces, edges, "    ");
        std::printf("  shapes2:\n");
        for (TopTools_ListIteratorOfListOfShape it(r.GetFaultyShapes2()); it.More(); it.Next())
            describe(it.Value(), faces, edges, "    ");
    }
    return 0;
}
