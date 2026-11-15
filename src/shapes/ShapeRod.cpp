#include "mvb/shapes/ShapeRod.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <stdexcept>

namespace mvb {
namespace shapes {

// Build frame (pre intrinsic -90degX rotation): the disc of diameter A extruded along Z from 0 to B.

TopoDS_Face ShapeRod::buildProfile(const std::map<std::string, double>& dims) const {
    if (!dims.count("A") || !(dims.at("A") > 0.0))
        throw std::runtime_error("ShapeRod: dimension A (rod diameter) is missing or non-positive");
    if (!dims.count("B") || !(dims.at("B") > 0.0))
        throw std::runtime_error("ShapeRod: dimension B (rod length) is missing or non-positive");
    TopoDS_Wire outer = build_polygon_circle(dims.at("A") / 2.0, m_corePolygonSegments);
    BRepBuilderAPI_MakeFace face(outer);
    return face.Face();
}

TopoDS_Shape ShapeRod::applyExtras(const std::map<std::string, double>& dims,
                                   const TopoDS_Shape& piece) const {
    const double b = dims.at("B");
    const double h = dims.count("H") ? dims.at("H") : 0.0;
    TopoDS_Shape result = piece;
    if (h > 0.0) {
        if (!(h < dims.at("A")))
            throw std::runtime_error("ShapeRod: bore H is not smaller than the rod diameter A");
        TopoDS_Shape bore = build_polygon_cylinder(b, h / 2.0, m_corePolygonSegments, /*circumscribed=*/true);
        BRepAlgoAPI_Cut cutter(result, bore);
        if (!cutter.IsDone()) throw std::runtime_error("ShapeRod: the bore cut failed (OCCT boolean)");
        result = cutter.Shape();
    }
    // Centred on the origin: after the intrinsic rotation the rod axis lands on Y.
    return translate_shape(result, 0.0, 0.0, -b / 2.0);
}

} // namespace shapes
} // namespace mvb
