#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include "mvb/shapes/ShapeBuilder.h"
#include "mvb/Utils.h"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <cmath>
#include <numbers>

namespace mvb {
namespace shapes {

TopoDS_Shape ShapeBuilder::buildPiece(const MAS::CoreShape& shapeData) const {
    auto dimsOpt = shapeData.get_dimensions();
    if (!dimsOpt) return TopoDS_Shape();
    auto dims = flatten_dimensions(*dimsOpt);

    // Build base profile and extrude
    TopoDS_Face profile = buildProfile(dims);
    double height = dims.count("B") ? dims.at("B") : 0.0;
    if (isToroidal()) {
        height = dims.count("C") ? dims.at("C") : height;
    }
    TopoDS_Shape piece = extrude(profile, height);

    // Cut winding window
    TopoDS_Shape window = buildWindingWindow(dims);
    if (!window.IsNull()) {
        piece = BRepAlgoAPI_Cut(piece, window).Shape();
    }

    // Apply extras
    piece = applyExtras(dims, piece);

    // Intrinsic rotation: concentric cores are built in XY-extruded-Z then rotated -90° around X
    if (!isToroidal()) {
        piece = rotate_shape(piece, -std::numbers::pi / 2.0, 0.0, 0.0);
    }

    return piece;
}

TopoDS_Shape ShapeBuilder::applyMachining(const TopoDS_Shape& piece,
                                          const MAS::Machining& machining,
                                          const std::map<std::string, double>& dims) const {
    // Generic rectangular machining tool oriented along Y (column axis)
    const std::vector<double>& coords = machining.get_coordinates();
    if (coords.size() < 2) throw std::runtime_error("applyMachining: machining has no coordinates");

    double gapLength = machining.get_length();
    if (std::abs(gapLength) < 1e-12) return piece; // zero-gap is a no-op
    double xCoord = coords[0];
    double yCoord = coords[1];

    TopoDS_Shape tool;
    double f = dims.count("F") ? dims.at("F") : 0.0;
    double c = dims.count("C") ? dims.at("C") : 0.0;
    double a = dims.count("A") ? dims.at("A") : 0.0;
    // For EFD cores the central column depth is F2 (not C). Use F2 when
    // available so the gap tool matches MVB Python's EFD.apply_machining.
    double centerColumnDepth = dims.count("F2") ? dims.at("F2") : c;

    if (std::abs(xCoord) < 1e-12) {
        // Center column: F (width) × gap × centerColumnDepth
        tool = makeBox(f, gapLength, centerColumnDepth);
    } else {
        // Side column: A/2 x gap x C (or A for EFD-family, matching MVB Python) -- OVERSIZED.
        // The leg's outer face sits exactly at x = A/2 and its depth is exactly C, so a tool of
        // that size has three faces coincident with the piece's own faces, and the parting plane
        // makes a fourth. That is the case OCCT's boolean drops silently: IsDone() true, piece
        // returned untouched. Measured 2026-09-05 on PQ 20/16: the centre cut (tool strictly
        // inside the section) removed 0.1452 mm3, each lateral cut removed 0.00000 mm3, and the
        // core reached the FEM with two of its three gaps missing. A knife must overhang what it
        // cuts: extend the tool beyond the leg in x and z, and away from the piece in y (nothing
        // is there to over-cut; the slab thickness on the piece's side is still exactly the gap).
        const double over = 1e-3;                             // 1 mm past every free face
        double sideDepth = dims.count("F2") ? a : c;
        const double sideW = a / 2.0 + over;                  // from x=0 out past A/2
        double sideX = (xCoord < 0) ? -sideW / 2.0 : sideW / 2.0;
        tool = makeBox(sideW, gapLength + over, sideDepth + 2.0 * over);
        // the extra `over` in y goes to the side AWAY from the piece: the piece lies on the side
        // of the parting plane that yCoord points to, so shift the box the other way by over/2
        tool = translate_shape(tool, sideX, (yCoord > 0 ? -over : over) / 2.0, 0.0);

        // Avoid cutting into center column (with matching depth)
        TopoDS_Shape centerTool = makeBox(f * 1.001, gapLength, centerColumnDepth * 1.001);
        BRepAlgoAPI_Cut cutter(tool, centerTool);
        if (!cutter.IsDone()) throw std::runtime_error("applyMachining: side-column tool minus centre column failed");
        tool = cutter.Shape();
    }

    tool = translate_shape(tool, 0.0, yCoord, 0.0);
    if (std::getenv("MVB_CORE_DIAG")) {
        Bnd_Box bp, bt; BRepBndLib::Add(piece, bp); BRepBndLib::Add(tool, bt);
        double p0[6], t0[6]; bp.Get(p0[0],p0[1],p0[2],p0[3],p0[4],p0[5]); bt.Get(t0[0],t0[1],t0[2],t0[3],t0[4],t0[5]);
        BRepAlgoAPI_Common com(piece, tool); GProp_GProps g; if (com.IsDone()) BRepGProp::VolumeProperties(com.Shape(), g);
        std::fprintf(stderr, "[machining-tool] x=%.4g y=%.4g: piece x[%.3f,%.3f] y[%.4f,%.4f] z[%.3f,%.3f] mm | tool x[%.3f,%.3f] y[%.4f,%.4f] z[%.3f,%.3f] mm | common %.5f mm3 (done=%d)\n",
            xCoord*1e3, yCoord*1e3, p0[0]*1e3,p0[3]*1e3,p0[1]*1e3,p0[4]*1e3,p0[2]*1e3,p0[5]*1e3,
            t0[0]*1e3,t0[3]*1e3,t0[1]*1e3,t0[4]*1e3,t0[2]*1e3,t0[5]*1e3, g.Mass()*1e9, (int)com.IsDone());
    }
    BRepAlgoAPI_Cut cutter(piece, tool);
    if (!cutter.IsDone()) throw std::runtime_error("applyMachining: the gap cut did not complete (OCCT boolean failed) at x=" + std::to_string(coords[0]) + " y=" + std::to_string(coords[1]) + " length=" + std::to_string(gapLength));
    return cutter.Shape();
}

TopoDS_Shape ShapeBuilder::buildWindingWindow(const std::map<std::string, double>&) const {
    return TopoDS_Shape();
}

TopoDS_Shape ShapeBuilder::applyExtras(const std::map<std::string, double>&,
                                        const TopoDS_Shape& piece) const {
    return piece;
}

TopoDS_Shape ShapeBuilder::extrude(const TopoDS_Face& face, double height) {
    gp_Vec vec(0, 0, height);
    return BRepPrimAPI_MakePrism(face, vec).Shape();
}

TopoDS_Shape ShapeBuilder::makeBox(double x, double y, double z) {
    gp_Pnt corner(-x / 2.0, -y / 2.0, -z / 2.0);
    return BRepPrimAPI_MakeBox(corner, x, y, z).Shape();
}

TopoDS_Shape ShapeBuilder::makeCylinder(double height, double radius, int segments) {
    return build_polygon_cylinder(height, radius, segments);
}

} // namespace shapes
} // namespace mvb
