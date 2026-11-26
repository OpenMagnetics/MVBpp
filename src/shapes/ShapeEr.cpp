#include <stdexcept>
#include "mvb/shapes/ShapeEr.h"
#include "mvb/Utils.h"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Trsf.hxx>
#include <gp_Pnt.hxx>
#include <cmath>
#include <numbers>

namespace mvb {
namespace shapes {

TopoDS_Shape ShapeEr::buildWindingWindow(const std::map<std::string, double>& dims) const {
    const double b = required(dims, "B");
    const double c = required(dims, "C");
    const double d = required(dims, "D");
    const double e = required(dims, "E");
    const double f = required(dims, "F");
    // G is the slot between the outer legs of an ER / EQ: an ETD or EER record has none, and MKF
    // reads an absent or zero G the same way (no aperture).
    const double g = featureDimension(dims, "G").value_or(0.0);

    double zCenter = b - d / 2.0;

    // Window-carving tool: its outer boundary becomes the window wall the winding
    // faces — circumscribe so the carved window contains the nominal window.
    TopoDS_Shape outerCyl =
        build_polygon_cylinder(d, e / 2.0, m_corePolygonSegments, /*circumscribed=*/true);
    outerCyl = translate_shape(outerCyl, 0.0, 0.0, zCenter - d / 2.0);

    TopoDS_Shape innerCyl = build_polygon_cylinder(d, f / 2.0, m_corePolygonSegments);
    innerCyl = translate_shape(innerCyl, 0.0, 0.0, zCenter - d / 2.0);

    TopoDS_Shape ww = done(BRepAlgoAPI_Cut(outerCyl, innerCyl), "leaving the round column in the winding window");

    // ER/EQ variant: when G > F and C > F, add a rectangular slot of
    // (G × C × D) centred on the column, minus the central column cylinder,
    // to the winding window cut. Matches MVB Python Er.get_negative_winding_window.
    if (g > f && c > f) {
        gp_Pnt boxCorner(-g / 2.0, -c / 2.0, zCenter - d / 2.0);
        TopoDS_Shape cube = BRepPrimAPI_MakeBox(boxCorner, g, c, d).Shape();
        TopoDS_Shape slot = done(BRepAlgoAPI_Cut(cube, innerCyl), "leaving the round column in the G slot");
        ww = done(BRepAlgoAPI_Fuse(ww, slot), "adding the G slot to the winding window");
    }
    return ww;
}

TopoDS_Shape ShapeEr::applyMachining(const TopoDS_Shape& piece,
                                     const MAS::Machining& machining,
                                     const std::map<std::string, double>& dims) const {
    const std::vector<double>& coords = machining.get_coordinates();
    if (coords.size() < 2) return piece;

    double gapLength = machining.get_length();
    // Skip zero-length residual gaps to avoid degenerate cylinder/box tools.
    if (std::abs(gapLength) < 1e-12) return piece;
    double xCoord = coords[0];
    double yCoord = coords[1];

    if (std::abs(xCoord) < 1e-12) {
        const double f = required(dims, "F");

        // Center column: polygon-faceted cylinder along Y axis, CENTERED on
        // yCoord. build_polygon_cylinder creates a base-at-z=0 prism along Z;
        // rotate -90° about X to point along +Y, then translate.
        // Gap-cutting tool: circumscribe so no faceted ridge is left in the gap.
        TopoDS_Shape toolZ = build_polygon_cylinder(gapLength, f / 2.0,
                                                    m_corePolygonSegments,
                                                    /*circumscribed=*/true);
        if (toolZ.IsNull()) return piece;
        TopoDS_Shape toolY = rotate_shape(toolZ, -std::numbers::pi / 2.0, 0.0, 0.0);
        TopoDS_Shape tool  = translate_shape(toolY, 0.0,
                                             yCoord - gapLength / 2.0, 0.0);
        if (tool.IsNull()) return piece;

        BRepAlgoAPI_Cut cutter(piece, tool);
        if (!cutter.IsDone()) throw std::runtime_error("applyMachining: the gap cut did not complete (OCCT boolean failed) at x=" + std::to_string(coords[0]) + " y=" + std::to_string(coords[1]) + " length=" + std::to_string(gapLength));
    return cutter.Shape();
    }

    // Side columns: use generic rectangular tool from base class
    return ShapeE::applyMachining(piece, machining, dims);
}

} // namespace shapes
} // namespace mvb
