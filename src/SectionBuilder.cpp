#include "mvb/SectionBuilder.h"

#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRep_Builder.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <GCPnts_TangentialDeflection.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_HSequenceOfShape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Pln.hxx>
#include <BRep_Tool.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <ElSLib.hxx>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace mvb {

namespace {

gp_Pln planeFor(SectionPlane plane) {
    switch (plane) {
        case SectionPlane::XY: return gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
        case SectionPlane::XZ: return gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0));
        case SectionPlane::YZ: return gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0));
    }
    return gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
}

gp_Pln planeFor(SectionPlane plane, double offset) {
    switch (plane) {
        case SectionPlane::XY: return gp_Pln(gp_Pnt(0, 0, offset), gp_Dir(0, 0, 1));
        case SectionPlane::XZ: return gp_Pln(gp_Pnt(0, offset, 0), gp_Dir(0, 1, 0));
        case SectionPlane::YZ: return gp_Pln(gp_Pnt(offset, 0, 0), gp_Dir(1, 0, 0));
    }
    return gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
}

// Project a 3D point onto the chosen plane's 2D coordinates.
// Convention: SVG +x = right, +y = down, so we flip the vertical axis.
std::pair<double, double> project(const gp_Pnt& p, SectionPlane plane) {
    switch (plane) {
        case SectionPlane::XY: return {p.X(), p.Y()};
        case SectionPlane::XZ: return {p.X(), p.Z()};
        case SectionPlane::YZ: return {p.Y(), p.Z()};
    }
    return {p.X(), p.Y()};
}

// Sample one edge into a polyline using OCCT's tangential deflection
// adaptive sampler — handles lines, circles, b-splines, etc.
std::vector<gp_Pnt> sampleEdge(const TopoDS_Edge& edge,
                                double angularDeflection = 0.05,
                                double curvatureDeflection = 1e-4) {
    std::vector<gp_Pnt> pts;
    BRepAdaptor_Curve curve(edge);
    GCPnts_TangentialDeflection sampler(curve, angularDeflection,
                                          curvatureDeflection);
    pts.reserve(sampler.NbPoints());
    for (int i = 1; i <= sampler.NbPoints(); ++i) {
        pts.push_back(sampler.Value(i));
    }
    return pts;
}

} // namespace

TopoDS_Shape SectionBuilder::sectionCore(const TopoDS_Shape& solid,
                                          SectionPlane plane) {
    if (solid.IsNull()) return TopoDS_Shape();
    BRepAlgoAPI_Section section(solid, planeFor(plane), Standard_False);
    section.ComputePCurveOn1(Standard_True);
    section.Approximation(Standard_True);
    section.Build();
    if (!section.IsDone()) return TopoDS_Shape();
    return section.Shape();
}

std::string SectionBuilder::edgesToSvg(const TopoDS_Shape& edges,
                                        SectionPlane plane,
                                        double width_px,
                                        double margin_px,
                                        double stroke_width,
                                        const std::string& stroke_color) {
    if (edges.IsNull()) {
        throw std::runtime_error("edgesToSvg: empty section");
    }

    // Sample every edge to a polyline and collect bbox in the projected plane.
    std::vector<std::vector<std::pair<double, double>>> polylines;
    double xmin = +1e30, xmax = -1e30, ymin = +1e30, ymax = -1e30;
    bool any = false;

    for (TopExp_Explorer exp(edges, TopAbs_EDGE); exp.More(); exp.Next()) {
        TopoDS_Edge edge = TopoDS::Edge(exp.Current());
        std::vector<gp_Pnt> pts;
        try { pts = sampleEdge(edge); }
        catch (...) { continue; }
        if (pts.size() < 2) continue;

        std::vector<std::pair<double, double>> poly;
        poly.reserve(pts.size());
        for (const auto& p : pts) {
            auto [u, v] = project(p, plane);
            poly.emplace_back(u, v);
            if (u < xmin) xmin = u; if (u > xmax) xmax = u;
            if (v < ymin) ymin = v; if (v > ymax) ymax = v;
            any = true;
        }
        polylines.push_back(std::move(poly));
    }

    if (!any) {
        throw std::runtime_error("edgesToSvg: section has no sampled edges");
    }

    // Scale-to-fit, with the SVG y-axis flipped (model +y is up, SVG +y is down).
    double model_w = xmax - xmin;
    double model_h = ymax - ymin;
    if (model_w <= 0 || model_h <= 0) {
        throw std::runtime_error("edgesToSvg: degenerate section bbox");
    }
    double inner_w = width_px - 2.0 * margin_px;
    double scale = inner_w / model_w;
    double height_px = model_h * scale + 2.0 * margin_px;

    auto toSvg = [&](double u, double v) -> std::pair<double, double> {
        double x = (u - xmin) * scale + margin_px;
        double y = height_px - ((v - ymin) * scale + margin_px);
        return {x, y};
    };

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(3);
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "width=\"" << width_px << "\" height=\"" << height_px << "\" "
        << "viewBox=\"0 0 " << width_px << " " << height_px << "\">\n";
    svg << "  <g fill=\"none\" stroke=\"" << stroke_color
        << "\" stroke-width=\"" << stroke_width
        << "\" stroke-linecap=\"round\" stroke-linejoin=\"round\">\n";

    for (const auto& poly : polylines) {
        if (poly.size() < 2) continue;
        svg << "    <path d=\"";
        for (size_t i = 0; i < poly.size(); ++i) {
            auto [x, y] = toSvg(poly[i].first, poly[i].second);
            svg << (i == 0 ? "M" : "L") << x << "," << y << " ";
        }
        svg << "\"/>\n";
    }

    svg << "  </g>\n";
    svg << "</svg>\n";
    return svg.str();
}

void SectionBuilder::writeSectionSvg(const TopoDS_Shape& solid,
                                      SectionPlane plane,
                                      const std::string& outputPath) {
    auto edges = sectionCore(solid, plane);
    auto svg = edgesToSvg(edges, plane);
    std::ofstream f(outputPath);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open " + outputPath + " for writing");
    }
    f << svg;
}

SectionPlane parseSectionPlane(const std::string& s) {
    std::string up = s;
    std::transform(up.begin(), up.end(), up.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    if (up == "XY") return SectionPlane::XY;
    if (up == "XZ") return SectionPlane::XZ;
    if (up == "YZ") return SectionPlane::YZ;
    throw std::runtime_error("Invalid section plane '" + s + "' (expected XY, XZ or YZ)");
}

std::vector<NamedShape> SectionBuilder::cut2DFaces(const std::vector<NamedShape>& shapes,
                                                    SectionPlane plane,
                                                    double offset) {
    std::vector<NamedShape> out;
    out.reserve(shapes.size());

    gp_Pln pln = planeFor(plane, offset);

    for (const auto& ns : shapes) {
        if (ns.shape.IsNull()) continue;

        BRepAlgoAPI_Section section(ns.shape, pln, Standard_False);
        section.ComputePCurveOn1(Standard_True);
        section.Approximation(Standard_True);
        section.Build();
        if (!section.IsDone()) continue;
        TopoDS_Shape edgesShape = section.Shape();
        if (edgesShape.IsNull()) continue;

        // Collect edges
        Handle(TopTools_HSequenceOfShape) edgeSeq = new TopTools_HSequenceOfShape();
        for (TopExp_Explorer exp(edgesShape, TopAbs_EDGE); exp.More(); exp.Next()) {
            edgeSeq->Append(exp.Current());
        }
        if (edgeSeq->IsEmpty()) continue;

        // Connect edges into wires
        Handle(TopTools_HSequenceOfShape) wireSeq;
        ShapeAnalysis_FreeBounds::ConnectEdgesToWires(edgeSeq, 1e-6, Standard_False, wireSeq);

        BRep_Builder builder;
        TopoDS_Compound compound;
        builder.MakeCompound(compound);
        std::vector<TopoDS_Shape> children;

        if (!wireSeq.IsNull()) {
            // Separate closed wires (candidate faces) from open ones.
            std::vector<TopoDS_Wire> cw;
            for (Standard_Integer i = 1; i <= wireSeq->Length(); ++i) {
                TopoDS_Wire wire = TopoDS::Wire(wireSeq->Value(i));
                if (wire.IsNull()) continue;
                if (wire.Closed()) cw.push_back(wire);
                else children.push_back(wire);   // open wire: keep as-is
            }
            // Build a face + interior test point (a vertex projected to the section plane)
            // for each closed wire, then nest CONTAINED wires as HOLES. A single solid still
            // sections into an outer loop PLUS inner loop(s) when it has a through-hole (a
            // toroid annulus, a closed winding window); making one face per wire would yield
            // overlapping disks where the inner disk fills the hole. Even containment depth =
            // solid outer boundary, odd = hole (punched out of its enclosing outer face).
            std::vector<TopoDS_Face> face(cw.size());
            std::vector<gp_Pnt2d>    pt(cw.size());
            std::vector<bool>        ok(cw.size(), false);
            std::vector<int>         depth(cw.size(), 0);
            for (std::size_t i = 0; i < cw.size(); ++i) {
                BRepBuilderAPI_MakeFace mk(pln, cw[i], Standard_True);
                TopExp_Explorer ve(cw[i], TopAbs_VERTEX);
                if (!mk.IsDone() || !ve.More()) { children.push_back(cw[i]); continue; }
                face[i] = mk.Face();
                gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(ve.Current()));
                Standard_Real u, v; ElSLib::Parameters(pln, p, u, v);
                pt[i] = gp_Pnt2d(u, v);
                ok[i] = true;
            }
            for (std::size_t i = 0; i < cw.size(); ++i) {
                if (!ok[i]) continue;
                for (std::size_t j = 0; j < cw.size(); ++j) {
                    if (j == i || !ok[j]) continue;
                    BRepClass_FaceClassifier cls(face[j], pt[i], 1e-7);
                    if (cls.State() == TopAbs_IN) depth[i]++;
                }
            }
            for (std::size_t i = 0; i < cw.size(); ++i) {
                if (!ok[i] || (depth[i] % 2) != 0) continue;   // skip holes (odd depth)
                // Subtract the directly-nested hole faces with a boolean cut (orientation-
                // robust: adding a wire as a hole depends on its winding, which differs
                // between polygon and true-circle section loops and silently inflated the
                // area). Cut the inner disks out of the outer disk -> the annulus.
                TopoDS_Shape result = face[i];
                for (std::size_t j = 0; j < cw.size(); ++j) {
                    if (j == i || !ok[j] || depth[j] != depth[i] + 1) continue;
                    BRepClass_FaceClassifier cls(face[i], pt[j], 1e-7);
                    if (cls.State() != TopAbs_IN) continue;
                    BRepAlgoAPI_Cut cut(result, face[j]);
                    if (cut.IsDone() && !cut.Shape().IsNull()) result = cut.Shape();
                }
                children.push_back(result);
            }
        }
        if (children.empty()) {
            // Fallback: keep raw edges
            for (TopExp_Explorer exp(edgesShape, TopAbs_EDGE); exp.More(); exp.Next()) {
                children.push_back(exp.Current());
            }
        }
        if (children.empty()) continue;

        // Emit one NamedShape per child so STEPCAFControl_Writer attaches
        // the source name to every product, not just the compound parent
        // (which OCCT auto-decomposes on write).  Multi-face sections get a
        // suffix so each face is uniquely identifiable.
        if (children.size() == 1) {
            out.push_back(NamedShape{children[0], ns.name});
        } else {
            for (std::size_t k = 0; k < children.size(); ++k) {
                out.push_back(NamedShape{children[k],
                                          ns.name + "_" + std::to_string(k)});
            }
        }
    }

    return out;
}

} // namespace mvb
