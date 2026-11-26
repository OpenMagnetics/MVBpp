#include "mvb/shapes/ShapeEpx.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

// EPX (EP variant with elongated/stadium central column).
// Mirrors MVB.js EPXShape.
TopoDS_Face ShapeEpx::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A") / 2.0;
    const double c = required(dims, "C");
    const double k = required(dims, "K");   // the oblong post's length: what makes an EPX
    const double f = required(dims, "F");

    double columnLength = k + f / 2.0;
    double topC = c - columnLength / 2.0;
    double bottomC = columnLength / 2.0;

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a, topC, 0.0));
    poly.Add(gp_Pnt( a, topC, 0.0));
    poly.Add(gp_Pnt( a, -bottomC, 0.0));
    poly.Add(gp_Pnt(-a, -bottomC, 0.0));
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

TopoDS_Shape ShapeEpx::buildWindingWindow(const std::map<std::string, double>& dims) const {
    const double b = required(dims, "B");
    const double c = required(dims, "C");
    const double d = required(dims, "D");
    const double e = required(dims, "E");
    const double f = required(dims, "F");
    const double k = required(dims, "K");
    // G is the slot above the post; a record without it has none (MKF reads it the same way).
    const double g = featureDimension(dims, "G").value_or(0.0);
    if (!(k > f / 2.0)) {
        throw std::runtime_error("shape '" + m_shapeName + "' has K = " + std::to_string(k) + " and F = " + std::to_string(f) +
                                 ": an EPX post is oblong, K must exceed F / 2");
    }

    double zBase = b - d;
    double rectangularPartWidth = k - f / 2.0;
    double columnWidth = k + f / 2.0;

    // Stadium-shaped central column: rectangle F × rectangularPartWidth + cylinders at ±rpw/2 in Y
    TopoDS_Shape centralCol;
    {
        gp_Pnt cn(-f / 2.0, -rectangularPartWidth / 2.0, zBase);
        TopoDS_Shape rect = BRepPrimAPI_MakeBox(cn, f, rectangularPartWidth, d).Shape();
        TopoDS_Shape topCyl = BRepPrimAPI_MakeCylinder(f / 2.0, d).Shape();
        topCyl = translate_shape(topCyl, 0.0, rectangularPartWidth / 2.0, zBase);
        TopoDS_Shape botCyl = BRepPrimAPI_MakeCylinder(f / 2.0, d).Shape();
        botCyl = translate_shape(botCyl, 0.0, -rectangularPartWidth / 2.0, zBase);
        centralCol = done(BRepAlgoAPI_Fuse(rect, topCyl), "joining the oblong post's first end");
        centralCol = done(BRepAlgoAPI_Fuse(centralCol, botCyl), "joining the oblong post's second end");
    }

    // Outer winding-window cylinder, offset by +rpw/2 in Y
    TopoDS_Shape ww = BRepPrimAPI_MakeCylinder(e / 2.0, d).Shape();
    ww = translate_shape(ww, 0.0, rectangularPartWidth / 2.0, zBase);
    ww = done(BRepAlgoAPI_Cut(ww, centralCol), "leaving the post in the winding window");

    // Top G-slot cutout if applicable
    if (g > 0.0) {
        gp_Pnt cn(-g / 2.0, c / 2.0 + columnWidth / 2.0 - c / 2.0, zBase);
        // MVB.js: translate [0, width/2 + columnWidth/2, ...] where width=C
        // so Y center = C/2 + columnWidth/2; box centered Y, so Y range = (C/2+columnWidth/2 - C/2, C/2+columnWidth/2 + C/2) = (columnWidth/2, C+columnWidth/2)
        // Actually box spans Y = (yc - C/2, yc + C/2) where yc = C/2 + columnWidth/2.
        gp_Pnt cn2(-g / 2.0, columnWidth / 2.0, zBase);
        TopoDS_Shape topBox = BRepPrimAPI_MakeBox(cn2, g, c, d).Shape();
        ww = done(BRepAlgoAPI_Fuse(ww, topBox), "adding the G slot to the winding window");
    }

    // Bottom rectangle E × C × D, offset by -C/2 + rpw/2 in Y, minus central column
    {
        // MVB.js: translate [0, -width/2 + rpw/2, ...] where width=C
        // Y center = -C/2 + rpw/2; box spans Y = (-C + rpw/2, rpw/2)
        gp_Pnt cn(-e / 2.0, -c + rectangularPartWidth / 2.0, zBase);
        TopoDS_Shape bottomBox = BRepPrimAPI_MakeBox(cn, e, c, d).Shape();
        bottomBox = done(BRepAlgoAPI_Cut(bottomBox, centralCol), "leaving the post in the open side of the window");
        ww = done(BRepAlgoAPI_Fuse(ww, bottomBox), "adding the open side to the winding window");
    }

    return ww;
}

} // namespace shapes
} // namespace mvb
