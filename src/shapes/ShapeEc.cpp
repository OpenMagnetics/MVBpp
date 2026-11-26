#include "mvb/shapes/ShapeEc.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

// EC: E-shape with two side dents (between ±(t+s) and ±a, height ±s).
// Mirrors MVB.js ECShape.getShapeBase.
TopoDS_Face ShapeEc::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A") / 2.0;
    const double c = required(dims, "C") / 2.0;
    const double t = required(dims, "T") / 2.0;
    const double sv = required(dims, "s") / 2.0;

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a,  c, 0.0));
    poly.Add(gp_Pnt( a,  c, 0.0));
    poly.Add(gp_Pnt( a,  sv, 0.0));
    poly.Add(gp_Pnt( t + sv,  sv, 0.0));
    poly.Add(gp_Pnt( t + sv, -sv, 0.0));
    poly.Add(gp_Pnt( a, -sv, 0.0));
    poly.Add(gp_Pnt( a, -c, 0.0));
    poly.Add(gp_Pnt(-a, -c, 0.0));
    poly.Add(gp_Pnt(-a, -sv, 0.0));
    poly.Add(gp_Pnt(-(t + sv), -sv, 0.0));
    poly.Add(gp_Pnt(-(t + sv),  sv, 0.0));
    poly.Add(gp_Pnt(-a,  sv, 0.0));
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

} // namespace shapes
} // namespace mvb
