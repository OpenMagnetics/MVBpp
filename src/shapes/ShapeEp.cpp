#include <stdexcept>
#include "mvb/shapes/ShapeEp.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

// EP cores: asymmetric profile (height C, top portion C-K, bottom K),
// circular winding window on top half + rectangular cut-out on bottom +
// optional G-slot above the central column. Mirrors MVB.js EPShape.
TopoDS_Face ShapeEp::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A") / 2.0;
    const double c = required(dims, "C");
    // K places the post in the depth of the pot. A record without it has the post half its own
    // diameter from the open face: MKF's rule for the same records (CorePieceEp, K = F / 2), so
    // the core is drawn where MKF computed it.
    const double k = featureDimension(dims, "K").value_or(required(dims, "F") / 2.0);

    double topC = c - k;
    double bottomC = k;

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a, topC, 0.0));
    poly.Add(gp_Pnt( a, topC, 0.0));
    poly.Add(gp_Pnt( a, -bottomC, 0.0));
    poly.Add(gp_Pnt(-a, -bottomC, 0.0));
    poly.Close();
    BRepBuilderAPI_MakeFace face(poly.Wire());
    return face.Face();
}

TopoDS_Shape ShapeEp::buildWindingWindow(const std::map<std::string, double>& dims) const {
    const double b = required(dims, "B");
    const double c = required(dims, "C");
    const double d = required(dims, "D");
    const double e = required(dims, "E");
    const double f = required(dims, "F");
    // G is the slot above the post; most EP records have none, and MKF reads an absent or zero G
    // the same way (no aperture).
    const double g = featureDimension(dims, "G").value_or(0.0);

    double zBase = b - d;

    // Outer cylinder (winding window)
    TopoDS_Shape ww = BRepPrimAPI_MakeCylinder(e / 2.0, d).Shape();
    ww = translate_shape(ww, 0.0, 0.0, zBase);

    // Subtract central column
    TopoDS_Shape centralCol = BRepPrimAPI_MakeCylinder(f / 2.0, d).Shape();
    centralCol = translate_shape(centralCol, 0.0, 0.0, zBase);
    ww = done(BRepAlgoAPI_Cut(ww, centralCol), "leaving the post in the winding window");

    // Top G-cube cutout (above the central column)
    if (g > 0.0) {
        gp_Pnt corner(-g / 2.0, f / 2.0, zBase);
        TopoDS_Shape topCube = BRepPrimAPI_MakeBox(corner, g, c, d).Shape();
        ww = done(BRepAlgoAPI_Fuse(ww, topCube), "adding the G slot to the winding window");
    }

    // Bottom rectangle E × C × D, minus central column (stays in WW set)
    {
        gp_Pnt corner(-e / 2.0, -c, zBase);
        TopoDS_Shape bottomCube = BRepPrimAPI_MakeBox(corner, e, c, d).Shape();
        // Cut central column from this bottom box
        TopoDS_Shape centralCol2 = BRepPrimAPI_MakeCylinder(f / 2.0, d).Shape();
        centralCol2 = translate_shape(centralCol2, 0.0, 0.0, zBase);
        bottomCube = done(BRepAlgoAPI_Cut(bottomCube, centralCol2), "leaving the post in the open side of the window");
        ww = done(BRepAlgoAPI_Fuse(ww, bottomCube), "adding the open side to the winding window");
    }

    return ww;
}

// EP center column is round (diameter F). The generic gap cutter (ShapeBuilder::
// applyMachining) sizes the center-column tool depth to C, a rectangular-column
// assumption; for EP, C >> F, so the box reaches past the round post into the closed pot
// back wall and clips the column corner at the -Z face. Size the cutter to F instead so
// it matches the round column exactly (its corners then fall in the bore = air, harmless).
// Side-column gaps (not used by round EP) defer to the generic handler.
TopoDS_Shape ShapeEp::applyMachining(const TopoDS_Shape& piece,
                                     const MAS::Machining& machining,
                                     const std::map<std::string, double>& dims) const {
    const std::vector<double>& coords = machining.get_coordinates();
    if (coords.size() < 2) return piece;
    const double gapLength = machining.get_length();
    if (std::abs(gapLength) < 1e-12) return piece;       // zero-gap no-op
    const double xCoord = coords[0];
    if (std::abs(xCoord) >= 1e-12)
        return ShapeBuilder::applyMachining(piece, machining, dims);  // side column
    const double f = required(dims, "F");
    TopoDS_Shape tool = makeBox(f, gapLength, f);         // F (width) x gap x F (round depth)
    tool = translate_shape(tool, 0.0, coords[1], 0.0);
    BRepAlgoAPI_Cut cutter(piece, tool);
    if (!cutter.IsDone()) throw std::runtime_error("applyMachining: the gap cut did not complete (OCCT boolean failed) at x=" + std::to_string(coords[0]) + " y=" + std::to_string(coords[1]) + " length=" + std::to_string(gapLength));
    return cutter.Shape();
}

TopoDS_Shape ShapeEp::applyExtras(const std::map<std::string, double>& dims,
                                  const TopoDS_Shape& piece) const {
    const double b = required(dims, "B");
    return translate_shape(piece, 0.0, 0.0, -b);
}

} // namespace shapes
} // namespace mvb
