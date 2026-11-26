#include "mvb/shapes/ShapePQ.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <gp_Pnt.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <cmath>
#include <numbers>

namespace mvb {
namespace shapes {

static TopoDS_Shape buildLateralPiece(double a, double c, double e, double sin_g, double cos_g,
                                       double j, double l, double height) {
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(j / 4.0, l / 4.0, 0.0));
    poly.Add(gp_Pnt(j / 4.0, -l / 4.0, 0.0));
    poly.Add(gp_Pnt(j / 2.0, -l / 2.0, 0.0));
    // When e*cos_g > c the arc extremum lies outside the flat-edge boundary,
    // producing a self-intersecting polygon.  Only include it when it's inside.
    if (e * cos_g < c) {
        poly.Add(gp_Pnt(e * sin_g, -e * cos_g, 0.0));
    }
    poly.Add(gp_Pnt(e * sin_g, -c, 0.0));
    poly.Add(gp_Pnt(a, -c, 0.0));
    poly.Add(gp_Pnt(a, c, 0.0));
    poly.Add(gp_Pnt(e * sin_g, c, 0.0));
    if (e * cos_g < c) {
        poly.Add(gp_Pnt(e * sin_g, e * cos_g, 0.0));
    }
    poly.Add(gp_Pnt(j / 2.0, l / 2.0, 0.0));
    poly.Close();

    BRepBuilderAPI_MakeWire wire(poly.Wire());
    BRepBuilderAPI_MakeFace face(wire.Wire());
    gp_Vec vec(0, 0, height);
    return BRepPrimAPI_MakePrism(face.Face(), vec).Shape();
}

TopoDS_Shape ShapePQ::buildPiece(const MAS::CoreShape& shapeData) const {
    const auto dims = dimensionsOf(shapeData);

    // Every one of these shapes the piece: an absent G used to build the lateral legs at an
    // angle made up from E and F instead of the opening the record should have given
    // (ABT #1143: PQ 32/12 had no G).
    const double a = required(dims, "A") / 2.0;
    const double b = required(dims, "B");
    const double c = required(dims, "C") / 2.0;
    const double d = required(dims, "D");
    const double e = required(dims, "E") / 2.0;
    const double f = required(dims, "F") / 2.0;
    const double g = required(dims, "G");
    if (!(g <= e * 2.0)) {
        throw std::runtime_error("shape '" + m_shapeName + "' has an opening G = " + std::to_string(g) +
                                 " wider than its window E = " + std::to_string(e * 2.0));
    }

    // J and L describe the lateral legs' inner flank. Most PQ records lack both, and MKF then
    // takes J = F / 2 and L = F + (C - F) / 3 (CorePiecePq, its documented approximation, the
    // one its Ae and le are computed with). The core is drawn with MKF's flank, never another.
    const auto recordJ = featureDimension(dims, "J");
    const auto recordL = featureDimension(dims, "L");
    if (recordJ.has_value() != recordL.has_value()) {
        throw std::runtime_error("shape '" + m_shapeName + "' carries only one of J and L; the lateral flank needs both or neither");
    }
    const double j = recordJ.value_or(f);
    const double l = recordL.value_or((f * 2.0) + (c * 2.0 - f * 2.0) / 3.0);

    const double g_angle = std::asin(g / (e * 2.0));
    double sin_g = std::sin(g_angle);
    double cos_g = std::cos(g_angle);

    // Central structure: two cylinders stacked.
    // 1) Connector disk (radius E/2, height B−D): solid yoke base that joins the
    //    column to the lateral wings — this is the "solid cylinder joining the two
    //    parts" of the PQ (column side + wing side).
    // 2) Real column (radius F/2, height D): the cylinder that protrudes into the
    //    winding space.  The winding-window cut will remove the E/2→F/2 annulus at
    //    this height, leaving only this column.
    double yoke_h = b - d;
    TopoDS_Shape central;
    {
        TopoDS_Shape connector = build_polygon_cylinder(
            (yoke_h > 1e-9 ? yoke_h : b), f, m_corePolygonSegments);
        TopoDS_Shape column    = build_polygon_cylinder(d, f, m_corePolygonSegments);
        // Translate column to sit on top of the connector disk.
        if (yoke_h > 1e-9) {
            gp_Trsf shift; shift.SetTranslation(gp_Vec(0, 0, yoke_h));
            column = BRepBuilderAPI_Transform(column, shift).Shape();
        }
        central = yoke_h > 1e-9 ? done(BRepAlgoAPI_Fuse(connector, column), "joining the column to the yoke disk") : connector;
    }

    // Right lateral piece
    TopoDS_Shape rightLateral = buildLateralPiece(a, c, e, sin_g, cos_g, j, l, b);

    // Left lateral piece (mirror across YZ plane)
    gp_Trsf mirror;
    mirror.SetMirror(gp_Ax2(gp_Pnt(0, 0, b / 2.0), gp_Dir(1, 0, 0)));
    TopoDS_Shape leftLateral = BRepBuilderAPI_Transform(rightLateral, mirror).Shape();

    // Fuse all three
    TopoDS_Shape fused = done(BRepAlgoAPI_Fuse(central, rightLateral), "joining the first lateral leg");
    fused = done(BRepAlgoAPI_Fuse(fused, leftLateral), "joining the second lateral leg");

    // Apply winding window cut
    TopoDS_Shape window = buildWindingWindow(dims);
    fused = done(BRepAlgoAPI_Cut(fused, window), "cutting the winding window");

    // Apply extras (translate by -B)
    fused = applyExtras(dims, fused);

    // Intrinsic rotation
    fused = rotate_shape(fused, -std::numbers::pi / 2.0, 0.0, 0.0);
    return fused;
}

} // namespace shapes
} // namespace mvb
