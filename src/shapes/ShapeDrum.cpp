#include "mvb/shapes/ShapeDrum.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <algorithm>

namespace mvb {
namespace shapes {

// Build frame (pre intrinsic -90degX rotation): extruded along Z from 0 to B, bottom
// flange at the base. The groove is cut as an annular band, exactly the ShapeP idiom.

TopoDS_Face ShapeDrum::buildProfile(const std::map<std::string, double>& dims) const {
    // A2 is the top flange's diameter where it differs from the bottom one's (A).
    const double flange = std::max(required(dims, "A"), featureDimension(dims, "A2").value_or(0.0));
    TopoDS_Wire outer = build_polygon_circle(flange / 2.0, m_corePolygonSegments);
    BRepBuilderAPI_MakeFace face(outer);
    return face.Face();
}

TopoDS_Shape ShapeDrum::buildWindingWindow(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A");
    const double c = required(dims, "C");
    const double e = required(dims, "E");
    const double f = required(dims, "F");

    // Annular groove: from the post surface out past the flange edge, spanning z in [F, F+E].
    TopoDS_Shape outerCyl = build_polygon_cylinder(e, a, m_corePolygonSegments);  // radius a: safely past any flange
    outerCyl = translate_shape(outerCyl, 0.0, 0.0, f);
    TopoDS_Shape post = build_polygon_cylinder(e, c / 2.0, m_corePolygonSegments);
    post = translate_shape(post, 0.0, 0.0, f);

    return done(BRepAlgoAPI_Cut(outerCyl, post), "leaving the post in the winding groove");
}

TopoDS_Shape ShapeDrum::applyExtras(const std::map<std::string, double>& dims,
                                    const TopoDS_Shape& piece) const {
    const double a = required(dims, "A");
    const double a2 = featureDimension(dims, "A2").value_or(a);
    const double b = required(dims, "B");
    const double d = required(dims, "D");
    const double f = required(dims, "F");
    // H is the bore through the drum; a record without it is solid.
    const double h = featureDimension(dims, "H").value_or(0.0);

    TopoDS_Shape result = piece;

    // Asymmetric flanges: the profile used the LARGER disc; trim the other flange's band
    // down to its own diameter. Convention: A = bottom flange, A2 = top flange.
    double larger = std::max(a, a2);
    if (a2 < larger || a < larger) {
        double smaller = std::min(a, a2);
        double bandBase = (a2 < a) ? (b - d) : 0.0;   // trim the top band if A2 is smaller, else the bottom
        double bandHeight = (a2 < a) ? d : f;
        TopoDS_Shape trimOuter = build_polygon_cylinder(bandHeight, larger, m_corePolygonSegments);
        trimOuter = translate_shape(trimOuter, 0.0, 0.0, bandBase);
        TopoDS_Shape trimInner = build_polygon_cylinder(bandHeight, smaller / 2.0, m_corePolygonSegments);
        trimInner = translate_shape(trimInner, 0.0, 0.0, bandBase);
        TopoDS_Shape ring = done(BRepAlgoAPI_Cut(trimOuter, trimInner), "making the ring that trims the smaller flange");
        result = done(BRepAlgoAPI_Cut(result, ring), "trimming the smaller flange");
    }

    // Bore through the whole piece.
    if (h > 0.0) {
        TopoDS_Shape bore =
            build_polygon_cylinder(b, h / 2.0, m_corePolygonSegments, /*circumscribed=*/true);
        result = done(BRepAlgoAPI_Cut(result, bore), "drilling the bore");
    }

    // Centre the single piece on the origin (open-shape cores are emitted as one CLOSED
    // element at {0,0,z}); after the intrinsic rotation the drum axis lands on Y.
    return translate_shape(result, 0.0, 0.0, -b / 2.0);
}

} // namespace shapes
} // namespace mvb
