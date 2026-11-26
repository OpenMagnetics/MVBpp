#include "mvb/shapes/ShapeT.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>

namespace mvb {
namespace shapes {

TopoDS_Face ShapeT::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A") / 2.0; // outer diameter
    const double b = required(dims, "B") / 2.0; // inner diameter
    if (!(b < a)) {
        throw std::runtime_error("shape '" + m_shapeName + "' has an inner diameter B that is not below its outer diameter A");
    }

    TopoDS_Wire outer = build_polygon_circle(a, m_corePolygonSegments);
    // The bore holds the winding: circumscribe so faceting never bites into it.
    TopoDS_Wire inner = build_polygon_circle(b, m_corePolygonSegments, /*circumscribed=*/true);
    inner.Reverse(); // inner hole must wind opposite to outer for correct face orientation

    BRepBuilderAPI_MakeFace face(outer);
    face.Add(inner);
    return face.Face();
}

TopoDS_Shape ShapeT::applyExtras(const std::map<std::string, double>& dims,
                                 const TopoDS_Shape& piece) const {
    const double c = required(dims, "C");
    return translate_shape(piece, 0.0, 0.0, -c / 2.0);
}

} // namespace shapes
} // namespace mvb
