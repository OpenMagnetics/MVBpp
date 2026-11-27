#include <cmath>
#include <stdexcept>
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

ShapeU::Legs ShapeU::legs(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A");
    const double e = required(dims, "E");
    if (!(e < a)) {
        throw std::runtime_error("shape '" + m_shapeName + "' has a window E that is not narrower than its width A");
    }
    // A record that carries H gives the wound leg's width outright, and that is the width MKF
    // sizes the bobbin and the turns around (CorePieceU::process_columns). The leg must not be
    // derived from A and E when H is there: a catalogue may store E as a minimum or as a
    // reference value (UI 93/76: E is a minimum of 34.6 mm, so (A - E) / 2 gives 29.2 mm for a
    // leg that is 28.0 mm on the datasheet). The return leg is what A, E and H leave: MKF winds
    // into a window E wide from the wound leg's face, so the return leg starts where that window
    // ends, and MKF's shape constants take the same A - E - H. Without H both legs are (A - E) / 2.
    const auto h = featureDimension(dims, "H");
    if (!h) {
        return {(a - e) / 2.0, (a - e) / 2.0};
    }
    const double returning = a - e - *h;
    if (!(returning > 0.0)) {
        throw std::runtime_error("shape '" + m_shapeName + "' leaves no return leg: A - E - H = " +
                                 std::to_string(returning));
    }
    return {*h, returning};
}

TopoDS_Face ShapeU::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A");
    const double c = required(dims, "C") / 2.0;
    const Legs leg = legs(dims);

    // Coordinate origin sits INSIDE the winding column (one leg), not at the
    // core centre — this is the MAS/MKF convention the bobbin and turns are
    // expressed in, so the core must follow it. The winding column is centred at
    // X=0; the window and the return leg extend to one side, and the outer width is A.
    double leftA = a - leg.winding / 2.0;  // outer edge of the far (return) leg
    double rightA = leg.winding / 2.0;     // outer edge of the winding column

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
    const double b = required(dims, "B");
    const double cFull = required(dims, "C");
    const double d = required(dims, "D");
    const double e = required(dims, "E");

    // The window occupies the gap between the winding column (centred at X=0,
    // inner edge at -windingColumnWidth/2) and the return leg, i.e.
    // X in [-(windingColumnWidth/2 + E), -windingColumnWidth/2].
    const double windingColumnWidth = legs(dims).winding;
    gp_Pnt corner(-(windingColumnWidth / 2.0 + e), -cFull, b - d);
    TopoDS_Shape ww = BRepPrimAPI_MakeBox(corner, e, 2.0 * cFull, d).Shape();
    return ww;
}

TopoDS_Shape ShapeU::applyExtras(const std::map<std::string, double>& dims,
                                 const TopoDS_Shape& piece) const {
    const double b = required(dims, "B");
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

    const double a = required(dims, "A");
    const double c = required(dims, "C");
    const double e = required(dims, "E");

    // The gap is as wide as the leg it sits in. MKF puts the wound leg at x = 0 and the return
    // leg at (A + E) / 2, so the nearer of the two centres names the leg.
    const Legs leg = legs(dims);
    const double returnLegCentre = (a + e) / 2.0;
    const bool inWoundLeg = std::abs(xCoord) < std::abs(std::abs(xCoord) - returnLegCentre);
    TopoDS_Shape tool = makeBox(inWoundLeg ? leg.winding : leg.returning, gapLength, c);
    tool = translate_shape(tool, xCoord, yCoord, 0.0);

    BRepAlgoAPI_Cut cutter(piece, tool);
    if (!cutter.IsDone()) throw std::runtime_error("applyMachining: the gap cut did not complete (OCCT boolean failed) at x=" + std::to_string(coords[0]) + " y=" + std::to_string(coords[1]) + " length=" + std::to_string(gapLength));
    return cutter.Shape();
}

} // namespace shapes
} // namespace mvb
