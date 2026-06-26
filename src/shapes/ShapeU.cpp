#include "mvb/shapes/ShapeU.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

TopoDS_Face ShapeU::buildProfile(const std::map<std::string, double>& dims) const {
    double a = 0.0, c = 0.0, e = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second / 2.0;
    it = dims.find("E"); if (it != dims.end()) e = it->second;

    // Coordinate origin sits INSIDE the winding column (one leg), not at the
    // core centre — this is the MAS/MKF convention the bobbin and turns are
    // expressed in, so the core must follow it. The winding column is centred at
    // X=0 (width windingColumnWidth); the window and the return leg extend to one
    // side. Both legs are equal width = (A-E)/2.
    double windingColumnWidth = (a - e) / 2.0;
    double leftA = a - windingColumnWidth / 2.0;  // outer edge of the far (return) leg
    double rightA = windingColumnWidth / 2.0;      // outer edge of the winding column

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(rightA, c, 0.0));
    poly.Add(gp_Pnt(-leftA, c, 0.0));
    poly.Add(gp_Pnt(-leftA, -c, 0.0));
    poly.Add(gp_Pnt(rightA, -c, 0.0));
    poly.Close();

    BRepBuilderAPI_MakeWire wire(poly.Wire());
    BRepBuilderAPI_MakeFace face(wire.Wire());
    return face.Face();
}

TopoDS_Shape ShapeU::buildWindingWindow(const std::map<std::string, double>& dims) const {
    double aFull = 0.0, b = 0.0, cFull = 0.0, d = 0.0, e = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) aFull = it->second;
    it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("C"); if (it != dims.end()) cFull = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;

    if (b == 0.0 || d == 0.0) return TopoDS_Shape();

    // The window occupies the gap between the winding column (centred at X=0,
    // inner edge at -windingColumnWidth/2) and the return leg, i.e.
    // X in [-(windingColumnWidth/2 + E), -windingColumnWidth/2]. The previous
    // -(wcw/2 + E/2) placed it half a window too far, cutting into the winding
    // column and leaving unequal legs.
    (void)aFull;
    double windingColumnWidth = (aFull - e) / 2.0;
    gp_Pnt corner(-(windingColumnWidth / 2.0 + e), -cFull, b - d);
    TopoDS_Shape ww = BRepPrimAPI_MakeBox(corner, e, 2.0 * cFull, d).Shape();
    return ww;
}

TopoDS_Shape ShapeU::applyExtras(const std::map<std::string, double>& dims,
                                 const TopoDS_Shape& piece) const {
    double b = 0.0;
    auto it = dims.find("B");
    if (it != dims.end()) b = it->second;
    return translate_shape(piece, 0.0, 0.0, -b);
}

TopoDS_Shape ShapeU::applyMachining(const TopoDS_Shape& piece,
                                    const MAS::Machining& machining,
                                    const std::map<std::string, double>& dims) const {
    const std::vector<double>& coords = machining.get_coordinates();
    if (coords.size() < 2) return piece;

    double gapLength = machining.get_length();
    // Zero-length residual gaps (emitted by MKF as part of multi-gap layouts)
    // would build a degenerate box and throw Standard_DomainError. Skip them.
    if (std::abs(gapLength) < 1e-12) return piece;
    double xCoord = coords[0];
    double yCoord = coords[1];

    double a = 0.0, c = 0.0, e = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;

    double windingColumnWidth = (a - e) / 2.0;
    TopoDS_Shape tool = makeBox(windingColumnWidth, gapLength, c);
    tool = translate_shape(tool, xCoord, yCoord, 0.0);

    BRepAlgoAPI_Cut cutter(piece, tool);
    return cutter.IsDone() ? cutter.Shape() : piece;
}

} // namespace shapes
} // namespace mvb
