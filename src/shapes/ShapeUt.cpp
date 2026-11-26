#include "mvb/shapes/ShapeUt.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace mvb {
namespace shapes {

// Mirrors MVB Python CadQueryBuilder.Ut: rectangle A×C profile, box WW,
// then two asymmetric columns added on top.
TopoDS_Shape ShapeUt::buildPiece(const MAS::CoreShape& shapeData) const {
    const auto dims = dimensionsOf(shapeData);

    const double a = required(dims, "A");
    const double b = required(dims, "B");
    const double c = required(dims, "C");
    const double d = required(dims, "D");
    const double e = required(dims, "E");
    const double fdim = required(dims, "F");

    double ha = a / 2.0, hc = c / 2.0;
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-ha,  hc, 0.0));
    poly.Add(gp_Pnt( ha,  hc, 0.0));
    poly.Add(gp_Pnt( ha, -hc, 0.0));
    poly.Add(gp_Pnt(-ha, -hc, 0.0));
    poly.Close();
    TopoDS_Face face = BRepBuilderAPI_MakeFace(poly.Wire()).Face();
    TopoDS_Shape piece = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, b)).Shape();

    // Cut WW: box(2A, 2C, D) translated to z = B/2 (centered at z=B/2)
    {
        gp_Pnt cn(-a, -c, b / 2.0 - d / 2.0);
        TopoDS_Shape ww = BRepPrimAPI_MakeBox(cn, 2.0 * a, 2.0 * c, d).Shape();
        piece = done(BRepAlgoAPI_Cut(piece, ww), "cutting the winding window");
    }

    // Top column: F × C × D at x = -A/2 + F/2, z = B/2
    {
        double xc = -a / 2.0 + fdim / 2.0;
        gp_Pnt cn(xc - fdim / 2.0, -hc, b / 2.0 - d / 2.0);
        TopoDS_Shape topCol = BRepPrimAPI_MakeBox(cn, fdim, c, d).Shape();
        piece = done(BRepAlgoAPI_Fuse(piece, topCol), "adding the first column");
    }

    // Bottom column: (A-E-F) × C × D at x = A/2 - bcw/2, z = B/2
    const double bcw = a - e - fdim;
    if (!(bcw > 0.0)) {
        throw std::runtime_error("shape '" + m_shapeName + "' leaves no second column: A - E - F = " + std::to_string(bcw));
    }
    {
        double xc = a / 2.0 - bcw / 2.0;
        gp_Pnt cn(xc - bcw / 2.0, -hc, b / 2.0 - d / 2.0);
        TopoDS_Shape botCol = BRepPrimAPI_MakeBox(cn, bcw, c, d).Shape();
        piece = done(BRepAlgoAPI_Fuse(piece, botCol), "adding the second column");
    }

    piece = translate_shape(piece, 0.0, 0.0, -b);
    piece = rotate_shape(piece, -std::numbers::pi / 2.0, 0.0, 0.0);
    return piece;
}

} // namespace shapes
} // namespace mvb
