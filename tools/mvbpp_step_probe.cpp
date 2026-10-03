// step_probe <file.step> <nameSubstring> [nameSubstringB] — inspect named solids in a STEP:
// Moved verbatim from OMFEM tools/omfem_step_probe.cpp at 14c2bde (ABT #1588, step 8); renamed. No logic changed.
// volume + bounding box; with a second filter, also the common shape's volume and bbox, so a
// reported "overlap" can be located and judged (a tangency sliver is thin and sits ON a face;
// a real interpenetration has depth).
#include <STEPCAFControl_Reader.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDataStd_Name.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <gp_Pnt.hxx>
#include <limits>
#include <BRepClass3d_SolidClassifier.hxx>
#include <cstdlib>
#include <BRepExtrema_DistShapeShape.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <TDataStd_Name.hxx>
#include <TCollection_ExtendedString.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <iostream>
#include <vector>
#include <string>

namespace {
struct Part { TopoDS_Shape shape; std::string name; double volume; Bnd_Box box; };

double vol(const TopoDS_Shape& s) {
    GProp_GProps p; BRepGProp::VolumeProperties(s, p); return p.Mass();
}
void printBox(const Bnd_Box& b) {
    double x0,y0,z0,x1,y1,z1; b.Get(x0,y0,z0,x1,y1,z1);
    printf("bbox x[%.3f,%.3f] y[%.3f,%.3f] z[%.3f,%.3f] mm", x0,x1,y0,y1,z0,z1);
}
}

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: step_probe <file.step> <nameA> [nameB]\n"; return 2; }
    Handle(TDocStd_Document) doc;
    XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
    STEPCAFControl_Reader reader; reader.SetNameMode(true);
    if (reader.ReadFile(argv[1]) != IFSelect_RetDone || !reader.Transfer(doc)) {
        std::cerr << "ERROR reading\n"; return 2;
    }
    Handle(XCAFDoc_ShapeTool) tool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
    TDF_LabelSequence labels; tool->GetFreeShapes(labels);
    std::vector<Part> parts;
    for (Standard_Integer i = 1; i <= labels.Length(); ++i) {
        Handle(TDataStd_Name) nm; std::string name = "shape" + std::to_string(i);
        if (labels.Value(i).FindAttribute(TDataStd_Name::GetID(), nm)) {
            name = std::string(TCollection_AsciiString(nm->Get()).ToCString());
        }
        TopoDS_Shape sh = tool->GetShape(labels.Value(i));
        if (sh.IsNull()) continue;
        int k = 0;
        for (TopExp_Explorer ex(sh, TopAbs_SOLID); ex.More(); ex.Next(), ++k) {
            Part p; p.shape = ex.Current();
            p.name = name + (k ? " [solid " + std::to_string(k) + "]" : "");
            p.volume = vol(p.shape); BRepBndLib::Add(p.shape, p.box);
            parts.push_back(p);
        }
    }
    // Point mode: step_probe <file> point x y z — classify one point against EVERY solid.
    // Boolean-free ground truth for "is this copper inside that part".
    if (std::string(argv[2]) == "point" && argc >= 6) {
        gp_Pnt probe(std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]));
        printf("point (%.4f, %.4f, %.4f):\n", probe.X(), probe.Y(), probe.Z());
        for (const auto& p : parts) {
            if (p.box.IsOut(probe)) continue;
            BRepClass3d_SolidClassifier cls(p.shape);
            cls.Perform(probe, 1e-7);
            const char* state = cls.State() == TopAbs_IN ? "IN" :
                                cls.State() == TopAbs_ON ? "ON" :
                                cls.State() == TopAbs_OUT ? "out" : "unknown";
            if (cls.State() == TopAbs_IN || cls.State() == TopAbs_ON) {
                printf("   %-42s %s\n", p.name.c_str(), state);
            }
        }
        return 0;
    }
    // Common mode: step_probe <file> common <nameA> <nameB> <out.step> — export the actual
    // intersection solids (plus the two parents, for context) as their own named STEP, so the
    // interpenetration can be opened and measured directly instead of argued about.
    if (std::string(argv[2]) == "common" && argc >= 6) {
        const std::string fa = argv[3], fb = argv[4];
        const std::string outPath = argv[5];
        Handle(TDocStd_Document) outDoc;
        XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", outDoc);
        Handle(XCAFDoc_ShapeTool) outTool = XCAFDoc_DocumentTool::ShapeTool(outDoc->Main());
        auto addNamed = [&](const TopoDS_Shape& sh, const std::string& nm) {
            TDF_Label lab = outTool->AddShape(sh, false);
            TDataStd_Name::Set(lab, TCollection_ExtendedString(nm.c_str()));
        };
        size_t found = 0;
        double total = 0;
        for (const auto& p : parts) {
            if (p.name.find(fa) == std::string::npos) continue;
            for (const auto& q : parts) {
                if (q.name.find(fb) == std::string::npos || q.name == p.name) continue;
                if (p.box.IsOut(q.box)) continue;
                BRepAlgoAPI_Common common(p.shape, q.shape);
                if (!common.IsDone() || common.Shape().IsNull()) continue;
                double v = 0, heaviest = 0;
                gp_Pnt centroid(0, 0, 0);
                for (TopExp_Explorer ex(common.Shape(), TopAbs_SOLID); ex.More(); ex.Next()) {
                    GProp_GProps gp; BRepGProp::VolumeProperties(ex.Current(), gp);
                    v += gp.Mass();
                    if (gp.Mass() > heaviest) { heaviest = gp.Mass(); centroid = gp.CentreOfMass(); }
                }
                if (v <= 1e-9) continue;
                // Same verification the audit uses: a common whose own centre of mass is not
                // inside both parents is an OCCT artefact, not shared material — exporting it
                // would put phantom geometry in front of a reviewer.
                BRepClass3d_SolidClassifier ca(p.shape), cb(q.shape);
                ca.Perform(centroid, 1e-7); cb.Perform(centroid, 1e-7);
                const bool inA = ca.State() == TopAbs_IN || ca.State() == TopAbs_ON;
                const bool inB = cb.State() == TopAbs_IN || cb.State() == TopAbs_ON;
                if (!inA || !inB) {
                    printf("rejected ARTEFACT %.6f mm3: '%s' x '%s' (common centroid outside %s)\n",
                           v, p.name.c_str(), q.name.c_str(), !inA ? p.name.c_str() : q.name.c_str());
                    continue;
                }
                ++found; total += v;
                printf("exporting COMMON %.6f mm3: '%s' x '%s'\n", v, p.name.c_str(), q.name.c_str());
                addNamed(common.Shape(), "COMMON " + p.name + " x " + q.name);
                addNamed(p.shape, "PARENT " + p.name);
                addNamed(q.shape, "PARENT " + q.name);
            }
        }
        if (!found) { printf("no common found for '%s' x '%s'\n", fa.c_str(), fb.c_str()); return 1; }
        STEPCAFControl_Writer writer;
        writer.SetNameMode(true);
        if (!writer.Transfer(outDoc, STEPControl_AsIs) ||
            writer.Write(outPath.c_str()) != IFSelect_RetDone) {
            printf("ERROR writing %s\n", outPath.c_str());
            return 2;
        }
        printf("wrote %s (%zu common solid group(s), total %.6f mm3)\n", outPath.c_str(), found, total);
        return 0;
    }
    // Extract mode: step_probe <file> extract <nameA> <nameB> <out.step> — write just the
    // named solids to their own STEP, so a contact/collision can be inspected on its own
    // instead of hunting for two pieces among thousands.
    if (std::string(argv[2]) == "extract" && argc >= 6) {
        Handle(TDocStd_Document) outDoc;
        XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", outDoc);
        Handle(XCAFDoc_ShapeTool) outTool = XCAFDoc_DocumentTool::ShapeTool(outDoc->Main());
        size_t written = 0;
        auto matches = [](const std::string& name, const std::string& filter) {
            // A leading '=' means EXACT: "Secondary parallel 0" must not drag in every
            // "Secondary parallel 0 [solid N]" when it is the bare first solid we want.
            if (!filter.empty() && filter[0] == '=') return name == filter.substr(1);
            return name.find(filter) != std::string::npos;
        };
        for (const auto& p : parts) {
            if (!matches(p.name, argv[3]) && !matches(p.name, argv[4])) {
                continue;
            }
            TDF_Label lab = outTool->AddShape(p.shape, false);
            TDataStd_Name::Set(lab, TCollection_ExtendedString(p.name.c_str()));
            printf("extracted %-42s vol=%.6f mm3\n", p.name.c_str(), p.volume);
            ++written;
        }
        if (!written) { printf("nothing matched\n"); return 1; }
        STEPCAFControl_Writer writer;
        writer.SetNameMode(true);
        if (!writer.Transfer(outDoc, STEPControl_AsIs) ||
            writer.Write(argv[5]) != IFSelect_RetDone) { printf("ERROR writing\n"); return 2; }
        printf("wrote %s (%zu solid(s))\n", argv[5], written);
        return 0;
    }
    // Distance mode: step_probe <file> dist <nameA> <nameB> — true minimum distance between
    // solids. Positive means they do NOT interpenetrate, whatever a boolean claims.
    if (std::string(argv[2]) == "dist" && argc >= 5) {
        const std::string fa = argv[3], fb = argv[4];
        for (const auto& p : parts) {
            if (p.name.find(fa) == std::string::npos) continue;
            for (const auto& q : parts) {
                if (q.name.find(fb) == std::string::npos || q.name == p.name) continue;
                BRepExtrema_DistShapeShape measure(p.shape, q.shape);
                if (!measure.IsDone()) continue;
                const double d = measure.Value();
                if (d > 0.5) continue;   // only report near neighbours
                printf("%-40s <-> %-40s  distance = %.6f mm (%.2f um)\n",
                       p.name.c_str(), q.name.c_str(), d, d * 1000);
            }
        }
        return 0;
    }
    const std::string a = argv[2];
    const std::string b = argc > 3 ? argv[3] : "";
    for (const auto& p : parts) {
        if (p.name.find(a) == std::string::npos) continue;
        // Vertex extents: the exact geometry, with none of Bnd_Box's tolerance inflation.
        double vy0 = std::numeric_limits<double>::max(), vy1 = -vy0;
        double vx0 = vy0, vx1 = -vy0, vz0 = vy0, vz1 = -vy0;
        for (TopExp_Explorer vx(p.shape, TopAbs_VERTEX); vx.More(); vx.Next()) {
            gp_Pnt pt = BRep_Tool::Pnt(TopoDS::Vertex(vx.Current()));
            vx0 = std::min(vx0, pt.X()); vx1 = std::max(vx1, pt.X());
            vy0 = std::min(vy0, pt.Y()); vy1 = std::max(vy1, pt.Y());
            vz0 = std::min(vz0, pt.Z()); vz1 = std::max(vz1, pt.Z());
        }
        BRepCheck_Analyzer analyzer(p.shape);
        printf("%-42s vol=%.5f mm3 %s\n    vertices x[%.4f,%.4f] y[%.4f,%.4f] z[%.4f,%.4f]  ",
               p.name.c_str(), p.volume, analyzer.IsValid() ? "VALID" : "*** INVALID ***",
               vx0, vx1, vy0, vy1, vz0, vz1);
        printBox(p.box); printf("\n");
        if (b.empty()) continue;
        for (const auto& q : parts) {
            if (q.name.find(b) == std::string::npos || q.name == p.name) continue;
            if (p.box.IsOut(q.box)) continue;
            BRepAlgoAPI_Common common(p.shape, q.shape);
            if (!common.IsDone() || common.Shape().IsNull()) continue;
            double v = 0; Bnd_Box cb;
            for (TopExp_Explorer ex(common.Shape(), TopAbs_SOLID); ex.More(); ex.Next()) {
                v += vol(ex.Current()); BRepBndLib::Add(ex.Current(), cb);
            }
            if (v <= 0) continue;
            double x0,y0,z0,x1,y1,z1; cb.Get(x0,y0,z0,x1,y1,z1);
            printf("    COMMON with %-34s vol=%.6f mm3 extent %.4f x %.4f x %.4f mm  ",
                   q.name.c_str(), v, x1-x0, y1-y0, z1-z0);
            printBox(cb); printf("\n");
        }
    }
    return 0;
}
