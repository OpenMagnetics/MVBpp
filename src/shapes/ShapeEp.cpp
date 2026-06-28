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
    double a = 0.0, c = 0.0, k = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second / 2.0;
    it = dims.find("C"); if (it != dims.end()) c = it->second;
    it = dims.find("K"); if (it != dims.end()) k = it->second;

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
    double b = 0.0, c = 0.0, d = 0.0, e = 0.0, f = 0.0, g = 0.0;
    auto it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;
    it = dims.find("G"); if (it != dims.end()) g = it->second;

    if (b == 0.0 || d == 0.0 || e == 0.0) return TopoDS_Shape();

    double zBase = b - d;

    // Outer cylinder (winding window)
    TopoDS_Shape ww = BRepPrimAPI_MakeCylinder(e / 2.0, d).Shape();
    ww = translate_shape(ww, 0.0, 0.0, zBase);

    // Subtract central column
    TopoDS_Shape centralCol = BRepPrimAPI_MakeCylinder(f / 2.0, d).Shape();
    centralCol = translate_shape(centralCol, 0.0, 0.0, zBase);
    {
        BRepAlgoAPI_Cut cut(ww, centralCol);
        if (cut.IsDone()) ww = cut.Shape();
    }

    // Top G-cube cutout (above the central column)
    if (g > 0.0) {
        gp_Pnt corner(-g / 2.0, f / 2.0, zBase);
        TopoDS_Shape topCube = BRepPrimAPI_MakeBox(corner, g, c, d).Shape();
        BRepAlgoAPI_Fuse fuser(ww, topCube);
        if (fuser.IsDone()) ww = fuser.Shape();
    }

    // Bottom rectangle E × C × D, minus central column (stays in WW set)
    {
        gp_Pnt corner(-e / 2.0, -c, zBase);
        TopoDS_Shape bottomCube = BRepPrimAPI_MakeBox(corner, e, c, d).Shape();
        // Cut central column from this bottom box
        TopoDS_Shape centralCol2 = BRepPrimAPI_MakeCylinder(f / 2.0, d).Shape();
        centralCol2 = translate_shape(centralCol2, 0.0, 0.0, zBase);
        BRepAlgoAPI_Cut cut(bottomCube, centralCol2);
        if (cut.IsDone()) bottomCube = cut.Shape();
        BRepAlgoAPI_Fuse fuser(ww, bottomCube);
        if (fuser.IsDone()) ww = fuser.Shape();
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
    const double f = dims.count("F") ? dims.at("F") : 0.0;
    TopoDS_Shape tool = makeBox(f, gapLength, f);         // F (width) x gap x F (round depth)
    tool = translate_shape(tool, 0.0, coords[1], 0.0);
    BRepAlgoAPI_Cut cutter(piece, tool);
    return cutter.IsDone() ? cutter.Shape() : piece;
}

TopoDS_Shape ShapeEp::applyExtras(const std::map<std::string, double>& dims,
                                  const TopoDS_Shape& piece) const {
    double b = 0.0;
    auto it = dims.find("B"); if (it != dims.end()) b = it->second;
    return translate_shape(piece, 0.0, 0.0, -b);
}

} // namespace shapes
} // namespace mvb
