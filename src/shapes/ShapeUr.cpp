#include <stdexcept>
#include "mvb/shapes/ShapeUr.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace mvb {
namespace shapes {

// Mirrors MVB.js URShape: 4 subtypes with rectangular profile + cylindrical
// winding/lateral columns.
TopoDS_Shape ShapeUr::buildPiece(const MAS::CoreShape& shapeData) const {
    const auto dims = dimensionsOf(shapeData);

    // The subtype says which letters make the piece: 1 and 2 have a winding column of diameter C,
    // 3 and 4 one of diameter F; 1 and 3 have a rectangular return leg H wide, 2 and 4 a round one.
    if (!shapeData.get_family_subtype()) {
        throw std::runtime_error("shape '" + m_shapeName + "' has no familySubtype; UR geometry is defined for types 1 to 4");
    }
    const std::string subtype = *shapeData.get_family_subtype();
    if (subtype != "1" && subtype != "2" && subtype != "3" && subtype != "4") {
        throw std::runtime_error("shape '" + m_shapeName + "' has familySubtype '" + subtype + "'; UR geometry is defined for types 1 to 4");
    }
    const bool columnIsF = subtype == "3" || subtype == "4";
    const bool legIsBox = subtype == "1" || subtype == "3";

    const double aFull = required(dims, "A");
    const double b = required(dims, "B");
    const double cFull = required(dims, "C");
    const double d = required(dims, "D");
    const double fFull = columnIsF ? required(dims, "F") : 0.0;
    const double legWidth = legIsBox ? required(dims, "H") : 0.0;
    double c = cFull / 2.0;

    double windingColumnWidth = 0.0;
    bool isHalfRect = (subtype == "1" || subtype == "3");

    if (subtype == "1" || subtype == "2") windingColumnWidth = cFull;
    else /* "3" or "4" */                windingColumnWidth = fFull;

    double leftA;
    if (isHalfRect) leftA = aFull - windingColumnWidth / 2.0;
    else            leftA = aFull - windingColumnWidth;

    BRepBuilderAPI_MakePolygon poly;
    if (isHalfRect) {
        poly.Add(gp_Pnt(0.0, c, 0.0));
        poly.Add(gp_Pnt(-leftA, c, 0.0));
        poly.Add(gp_Pnt(-leftA, -c, 0.0));
        poly.Add(gp_Pnt(0.0, -c, 0.0));
        poly.Add(gp_Pnt(windingColumnWidth / 2.0, -c, 0.0));
        poly.Add(gp_Pnt(windingColumnWidth / 2.0, c, 0.0));
    } else {
        poly.Add(gp_Pnt(0.0, c, 0.0));
        poly.Add(gp_Pnt(-leftA, c, 0.0));
        poly.Add(gp_Pnt(-leftA, -c, 0.0));
        poly.Add(gp_Pnt(0.0, -c, 0.0));
    }
    poly.Close();

    TopoDS_Face face = BRepBuilderAPI_MakeFace(poly.Wire()).Face();
    TopoDS_Shape piece = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, b)).Shape();

    // Winding window cutout: box(2A, 2C, D) at top
    {
        gp_Pnt corner(-aFull, -cFull, b - d);
        TopoDS_Shape ww = BRepPrimAPI_MakeBox(corner, 2.0 * aFull, 2.0 * cFull, d).Shape();
        piece = done(BRepAlgoAPI_Cut(piece, ww), "cutting the winding window");
    }

    // Add winding & lateral columns
    auto fuse = [&](TopoDS_Shape& base, const TopoDS_Shape& add) {
        base = done(BRepAlgoAPI_Fuse(base, add), "adding a column");
    };

    if (subtype == "1") {
        TopoDS_Shape wcCol = BRepPrimAPI_MakeCylinder(cFull / 2.0, d).Shape();
        wcCol = translate_shape(wcCol, 0.0, 0.0, b - d);
        fuse(piece, wcCol);
        // lateral box H × C × D at -(A - wcw/2 - H/2)
        const double h = legWidth;
        {
            double xOff = -(aFull - cFull / 2.0 - h / 2.0);
            gp_Pnt cn(xOff - h / 2.0, -cFull / 2.0, b - d);
            TopoDS_Shape latBox = BRepPrimAPI_MakeBox(cn, h, cFull, d).Shape();
            fuse(piece, latBox);
        }
    } else if (subtype == "2") {
        TopoDS_Shape wcCol = BRepPrimAPI_MakeCylinder(cFull / 2.0, b).Shape();
        fuse(piece, wcCol);
        TopoDS_Shape latCol = BRepPrimAPI_MakeCylinder(cFull / 2.0, b).Shape();
        latCol = translate_shape(latCol, -(aFull - cFull), 0.0, 0.0);
        fuse(piece, latCol);
    } else if (subtype == "3") {
        TopoDS_Shape wcCol = BRepPrimAPI_MakeCylinder(fFull / 2.0, d).Shape();
        wcCol = translate_shape(wcCol, 0.0, 0.0, b - d);
        fuse(piece, wcCol);
        const double h = legWidth;
        {
            double xOff = -(aFull - fFull / 2.0 - h / 2.0);
            gp_Pnt cn(xOff - h / 2.0, -cFull / 2.0, b - d);
            TopoDS_Shape latBox = BRepPrimAPI_MakeBox(cn, h, cFull, d).Shape();
            fuse(piece, latBox);
        }
    } else if (subtype == "4") {
        TopoDS_Shape wcCol = BRepPrimAPI_MakeCylinder(fFull / 2.0, b).Shape();
        fuse(piece, wcCol);
        TopoDS_Shape latCol = BRepPrimAPI_MakeCylinder(fFull / 2.0, b).Shape();
        latCol = translate_shape(latCol, -(aFull - fFull), 0.0, 0.0);
        fuse(piece, latCol);
    }

    piece = translate_shape(piece, 0.0, 0.0, -b);
    piece = rotate_shape(piece, -std::numbers::pi / 2.0, 0.0, 0.0);
    return piece;
}

// UR winding column is cylindrical (radius F/2 for subtypes 3/4, C/2 for 1/2).
// MKF emits gap coordinates with xCoord=0 referring to the winding column
// centre. After the piece is rotated -90° about X, the column axis lies along
// world Y, so the gap tool is a Y-aligned cylinder cut from the piece.
TopoDS_Shape ShapeUr::applyMachining(const TopoDS_Shape& piece,
                                     const MAS::Machining& machining,
                                     const std::map<std::string, double>& dims) const {
    const std::vector<double>& coords = machining.get_coordinates();
    if (coords.size() < 2) return piece;

    double gapLength = machining.get_length();
    if (std::abs(gapLength) < 1e-12) return piece;

    double xCoord = coords[0];
    double yCoord = coords[1];

    // Winding-column gap (xCoord == 0): cylindrical tool aligned with Y.
    if (std::abs(xCoord) < 1e-12) {
        // The winding column's diameter: F where the record carries it (types 3 and 4), else C.
        const double radius = featureDimension(dims, "F").value_or(required(dims, "C")) / 2.0;

        // Gap-cutting tool: circumscribe so no faceted ridge is left in the gap.
        TopoDS_Shape toolZ = build_polygon_cylinder(gapLength, radius,
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

    // Lateral-column gap: rectangular box across the lateral column.
    const double h = required(dims, "H");
    const double c = required(dims, "C");

    TopoDS_Shape tool = makeBox(h, gapLength, c);
    tool = translate_shape(tool, xCoord, yCoord, 0.0);
    BRepAlgoAPI_Cut cutter(piece, tool);
    if (!cutter.IsDone()) throw std::runtime_error("applyMachining: the gap cut did not complete (OCCT boolean failed) at x=" + std::to_string(coords[0]) + " y=" + std::to_string(coords[1]) + " length=" + std::to_string(gapLength));
    return cutter.Shape();
}

} // namespace shapes
} // namespace mvb
