#include "mvb/shapes/ShapeEpc.h"
#include "mvb/Utils.h"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>

namespace mvb {
namespace shapes {

namespace {

// The pole's section in the plane z = 0: `length` along X, `thickness` along Y, each end a half
// disc of diameter `thickness`. Exact arcs when segments <= 0; otherwise each half disc is an
// inscribed polygon (the pole is core material and must not bite into the winding's room).
TopoDS_Face stadiumFace(double length, double thickness, int segments) {
    const double r = thickness / 2.0;
    const double s = (length - thickness) / 2.0;   // half the straight part
    if (segments <= 0) {
        BRepBuilderAPI_MakeWire wire;
        if (s > 0.0) wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(-s, -r, 0.0), gp_Pnt(s, -r, 0.0)).Edge());
        wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(gp_Pnt(s, -r, 0.0), gp_Pnt(s + r, 0.0, 0.0), gp_Pnt(s, r, 0.0)).Value()).Edge());
        if (s > 0.0) wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(s, r, 0.0), gp_Pnt(-s, r, 0.0)).Edge());
        wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(gp_Pnt(-s, r, 0.0), gp_Pnt(-s - r, 0.0, 0.0), gp_Pnt(-s, -r, 0.0)).Value()).Edge());
        return BRepBuilderAPI_MakeFace(wire.Wire()).Face();
    }
    const int perEnd = std::max(2, segments / 2);
    BRepBuilderAPI_MakePolygon poly;
    for (int end = 0; end < 2; ++end) {
        const double centre = end == 0 ? s : -s;
        const double start = end == 0 ? -std::numbers::pi / 2.0 : std::numbers::pi / 2.0;
        // A round pole (no straight part) shares its end points between the two halves.
        const int last = s > 0.0 ? perEnd : perEnd - 1;
        for (int i = 0; i <= last; ++i) {
            const double angle = start + std::numbers::pi * i / perEnd;
            poly.Add(gp_Pnt(centre + r * std::cos(angle), r * std::sin(angle), 0.0));
        }
    }
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

}  // namespace

// The footprint on the mating face: A wide and C deep, centred on the pole's axis.
TopoDS_Face ShapeEpc::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A") / 2.0;
    const double c = required(dims, "C") / 2.0;
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a, c, 0.0));
    poly.Add(gp_Pnt(a, c, 0.0));
    poly.Add(gp_Pnt(a, -c, 0.0));
    poly.Add(gp_Pnt(-a, -c, 0.0));
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

TopoDS_Shape ShapeEpc::buildWindingWindow(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A");
    const double b = required(dims, "B");
    const double c = required(dims, "C");
    const double d = required(dims, "D");
    const double e = required(dims, "E");
    const double f = required(dims, "F");
    const double f2 = required(dims, "F2");
    // The same reading of the record as MKF's CorePieceEpc: a pole that is a stadium long along A,
    // inside a closed E-like section.
    if (f2 > f) {
        throw std::runtime_error("shape '" + m_shapeName + "': the pole thickness F2 = " + std::to_string(f2) +
                                 " is larger than its length F = " + std::to_string(f) + "; an EPC pole is a stadium long along A");
    }
    if (f2 > c || f >= e || e >= a || d >= b) {
        throw std::runtime_error("shape '" + m_shapeName + "' is not a closed E-like section: need F2 <= C, F < E < A and D < B");
    }

    // The window between the outer legs, E wide and D high, through the whole depth (2C: past
    // both faces of the A x C profile), minus the pole.
    const double zBase = b - d;
    TopoDS_Shape window = BRepPrimAPI_MakeBox(gp_Pnt(-e / 2.0, -c, zBase), e, 2.0 * c, d).Shape();
    TopoDS_Shape pole = BRepPrimAPI_MakePrism(stadiumFace(f, f2, m_corePolygonSegments), gp_Vec(0.0, 0.0, d)).Shape();
    pole = translate_shape(pole, 0.0, 0.0, zBase);
    return done(BRepAlgoAPI_Cut(window, pole), "leaving the stadium pole in the winding window");
}

TopoDS_Shape ShapeEpc::applyExtras(const std::map<std::string, double>& dims,
                                   const TopoDS_Shape& piece) const {
    return translate_shape(piece, 0.0, 0.0, -required(dims, "B"));
}

} // namespace shapes
} // namespace mvb
