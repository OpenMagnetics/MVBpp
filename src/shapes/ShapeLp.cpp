#include "mvb/shapes/ShapeLp.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

// LP shares EP's winding window topology (donut + top G-slot + bottom E-rect),
// but uses a centered A × C rectangular profile (ER-like).
TopoDS_Face ShapeLp::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A") / 2.0;
    const double c = required(dims, "C") / 2.0;

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a,  c, 0.0));
    poly.Add(gp_Pnt( a,  c, 0.0));
    poly.Add(gp_Pnt( a, -c, 0.0));
    poly.Add(gp_Pnt(-a, -c, 0.0));
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

} // namespace shapes
} // namespace mvb
