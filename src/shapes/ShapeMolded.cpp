#include "mvb/shapes/ShapeMolded.h"
#include "mvb/Utils.h"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>

namespace mvb {
namespace shapes {

// Build frame (pre intrinsic -90degX rotation): the block is extruded along Z from 0 to B by
// ShapeBuilder::buildPiece, so Z is the coil axis here and the footprint lies in XY.
TopoDS_Face ShapeMolded::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A") / 2.0;
    const double c = required(dims, "C") / 2.0;

    // A rectangular body footprint: A wide, C deep, centred on the axis. Moulded inductors
    // are plain blocks — the shaping is all internal.
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a, -c, 0.0));
    poly.Add(gp_Pnt( a, -c, 0.0));
    poly.Add(gp_Pnt( a,  c, 0.0));
    poly.Add(gp_Pnt(-a,  c, 0.0));
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

TopoDS_Shape ShapeMolded::buildWindingWindow(const std::map<std::string, double>& dims) const {
    const double b = required(dims, "B");
    const double d = required(dims, "D");   // cavity height
    const double e = required(dims, "E");   // cavity outer diameter
    // The cavity is the annulus the winding occupies: outer diameter E, inner diameter F,
    // height D, centred on the body so the composite closes over and under the coil. F is the
    // post: a record without it (or with F = 0) has none, and the cavity is a plain cylinder.
    const double f = featureDimension(dims, "F").value_or(0.0);

    // A cavity that does not fit its body is refused: the solid block drawn (and meshed) in its
    // place would be a core with no room for the coil inside it.
    if (e <= f || d >= b) {
        throw std::runtime_error("shape '" + m_shapeName + "' has a cavity that does not fit: E = " + std::to_string(e) +
                                 ", F = " + std::to_string(f) + ", D = " + std::to_string(d) + ", B = " + std::to_string(b) +
                                 " (need E > F and D < B)");
    }

    double base = (b - d) / 2.0;
    // Circumscribe the cavity: faceting must never leave a sliver of composite inside the
    // volume the winding occupies, or the turns intersect the core in the render.
    TopoDS_Shape cavity = build_polygon_cylinder(d, e / 2.0, m_corePolygonSegments,
                                                 /*circumscribed=*/true);
    cavity = translate_shape(cavity, 0.0, 0.0, base);
    if (f <= 0.0) {
        return cavity;
    }

    // The post is inscribed for the mirrored reason: it is composite that must not eat into
    // the winding volume.
    TopoDS_Shape post = build_polygon_cylinder(d, f / 2.0, m_corePolygonSegments);
    post = translate_shape(post, 0.0, 0.0, base);

    return done(BRepAlgoAPI_Cut(cavity, post), "leaving the post in the cavity");
}

TopoDS_Shape ShapeMolded::applyExtras(const std::map<std::string, double>& dims,
                                      const TopoDS_Shape& piece) const {
    // Centre the body on the origin, the same closing move ShapeDrum makes.
    //
    // ShapeBuilder::buildPiece extrudes the profile from z = 0 to B, so a piece that does not
    // translate is left sitting entirely on one side of the plane. MKF places the winding
    // window at y = 0 (CorePieceMolded sets the window coordinates with no vertical offset),
    // so an uncentred body renders half a body-height away from its own coil — the block
    // floating above the winding rather than moulded around it.
    //
    // Every other family here ends with this translation; ShapeMolded originally had no
    // applyExtras at all, which is how it was missed.
    const double b = required(dims, "B");
    return translate_shape(piece, 0.0, 0.0, -b / 2.0);
}

} // namespace shapes
} // namespace mvb
