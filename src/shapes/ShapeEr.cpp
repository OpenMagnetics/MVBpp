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
    double b = 0.0, c = 0.0, d = 0.0, e = 0.0, f = 0.0, g = 0.0;
    auto it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;
    it = dims.find("G"); if (it != dims.end()) g = it->second;

    if (b == 0.0 || d == 0.0 || e == 0.0) return TopoDS_Shape();

    double zCenter = b - d / 2.0;

    TopoDS_Shape outerCyl = build_polygon_cylinder(d, e / 2.0, m_corePolygonSegments);
    outerCyl = translate_shape(outerCyl, 0.0, 0.0, zCenter - d / 2.0);

    TopoDS_Shape innerCyl = build_polygon_cylinder(d, f / 2.0, m_corePolygonSegments);
    innerCyl = translate_shape(innerCyl, 0.0, 0.0, zCenter - d / 2.0);

    BRepAlgoAPI_Cut ringCut(outerCyl, innerCyl);
    if (!ringCut.IsDone()) return TopoDS_Shape();
    TopoDS_Shape ww = ringCut.Shape();

    // ER/EQ variant: when G > F and C > F, add a rectangular slot of
    // (G × C × D) centred on the column, minus the central column cylinder,
    // to the winding window cut. Matches MVB Python Er.get_negative_winding_window.
    if (g > f && c > f) {
        gp_Pnt boxCorner(-g / 2.0, -c / 2.0, zCenter - d / 2.0);
        TopoDS_Shape cube = BRepPrimAPI_MakeBox(boxCorner, g, c, d).Shape();
        BRepAlgoAPI_Cut cubeMinus(cube, innerCyl);
        if (cubeMinus.IsDone()) {
            BRepAlgoAPI_Fuse fuser(ww, cubeMinus.Shape());
            if (fuser.IsDone()) ww = fuser.Shape();
        }
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
        double f = 0.0;
        auto it = dims.find("F");
        if (it != dims.end()) f = it->second;
        if (f <= 0.0) return piece;

        // Center column: polygon-faceted cylinder along Y axis, CENTERED on
        // yCoord. build_polygon_cylinder creates a base-at-z=0 prism along Z;
        // rotate -90° about X to point along +Y, then translate.
        TopoDS_Shape toolZ = build_polygon_cylinder(gapLength, f / 2.0,
                                                    m_corePolygonSegments);
        if (toolZ.IsNull()) return piece;
        TopoDS_Shape toolY = rotate_shape(toolZ, -std::numbers::pi / 2.0, 0.0, 0.0);
        TopoDS_Shape tool  = translate_shape(toolY, 0.0,
                                             yCoord - gapLength / 2.0, 0.0);
        if (tool.IsNull()) return piece;

        BRepAlgoAPI_Cut cutter(piece, tool);
        return cutter.IsDone() ? cutter.Shape() : piece;
    }

    // Side columns: use generic rectangular tool from base class
    return ShapeE::applyMachining(piece, machining, dims);
}

} // namespace shapes
} // namespace mvb
