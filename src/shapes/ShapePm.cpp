#include "mvb/shapes/ShapePm.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Trsf.hxx>
#include <gp_Ax1.hxx>
#include <cmath>
#include <numbers>

namespace mvb {
namespace shapes {

// Mirrors MVB.js PMShape.getShapeExtras(): cuts two wedge openings on
// ±Y, fuses the central F-radius column, optionally cuts the central H
// hole, then translates so the bottom is at Z = -B (matching ShapeP).
TopoDS_Shape ShapePm::applyExtras(const std::map<std::string, double>& dims,
                                  const TopoDS_Shape& piece) const {
    const double a = required(dims, "A") / 2.0;
    const double b = required(dims, "B");
    const double c = required(dims, "C") / 2.0;
    const double f = required(dims, "F") / 2.0;
    // H is the centre hole; a record without it has a solid post.
    const double h = featureDimension(dims, "H").value_or(0.0) / 2.0;
    // The opening angle of the two wedges, in degrees: every PM record carries it.
    const double alpha = required(dims, "alpha");
    if (!(alpha < 180.0)) {
        throw std::runtime_error("shape '" + m_shapeName + "' has a wedge angle alpha of " + std::to_string(alpha) + " degrees; it must be below 180");
    }

    double halfAlphaRad = (alpha / 2.0) * std::numbers::pi / 180.0;
    double wedgeLength = a * 1.5;
    double wedgeHalfWidth = wedgeLength * std::tan(halfAlphaRad);

    // Build wedge profile in XY, extrude by B in +Z.
    auto makeWedge = [&]() {
        BRepBuilderAPI_MakePolygon poly;
        if (familySubtype_ == "1") {
            poly.Add(gp_Pnt(0.0, 0.0, 0.0));
            poly.Add(gp_Pnt(wedgeHalfWidth, wedgeLength, 0.0));
            poly.Add(gp_Pnt(-wedgeHalfWidth, wedgeLength, 0.0));
        } else {
            poly.Add(gp_Pnt(-f, 0.0, 0.0));
            poly.Add(gp_Pnt(f, 0.0, 0.0));
            poly.Add(gp_Pnt(wedgeHalfWidth, wedgeLength, 0.0));
            poly.Add(gp_Pnt(-wedgeHalfWidth, wedgeLength, 0.0));
        }
        poly.Close();
        BRepBuilderAPI_MakeFace face(poly.Wire());
        return BRepPrimAPI_MakePrism(face.Face(), gp_Vec(0, 0, b)).Shape();
    };

    TopoDS_Shape topWedge = translate_shape(makeWedge(), 0.0, c, 0.0);
    TopoDS_Shape result = piece;
    result = done(BRepAlgoAPI_Cut(result, topWedge), "cutting the first wedge opening");

    // Bottom wedge: rotate 180° around Z then translate to -c
    TopoDS_Shape bottomWedge = makeWedge();
    {
        gp_Trsf rot;
        rot.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,0,1)), std::numbers::pi);
        bottomWedge = BRepBuilderAPI_Transform(bottomWedge, rot).Shape();
    }
    bottomWedge = translate_shape(bottomWedge, 0.0, -c, 0.0);
    result = done(BRepAlgoAPI_Cut(result, bottomWedge), "cutting the second wedge opening");

    // Add central column at full height B
    {
        TopoDS_Shape column = BRepPrimAPI_MakeCylinder(f, b).Shape();
        result = done(BRepAlgoAPI_Fuse(result, column), "adding the centre post");
    }

    // Optional central hole H
    if (h > 0.0) {
        TopoDS_Shape hole = BRepPrimAPI_MakeCylinder(h, b).Shape();
        result = done(BRepAlgoAPI_Cut(result, hole), "drilling the centre hole");
    }

    return translate_shape(result, 0.0, 0.0, -b);
}

} // namespace shapes
} // namespace mvb
