#include <stdexcept>
#include <algorithm>
#include "mvb/shapes/ShapeP.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <gp_Pnt.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

namespace mvb {
namespace shapes {

TopoDS_Face ShapeP::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A") / 2.0;

    TopoDS_Wire outer = build_polygon_circle(a, m_corePolygonSegments);
    BRepBuilderAPI_MakeFace face(outer);
    return face.Face();
}

TopoDS_Shape ShapeP::buildWindingWindow(const std::map<std::string, double>& dims) const {
    const double b = required(dims, "B");
    const double d = required(dims, "D");
    const double e = required(dims, "E");
    const double f = required(dims, "F");

    double zCenter = b - d / 2.0;

    // build_polygon_cylinder creates a base-at-zero prism; shift down by d/2 to center it at zCenter
    // Window-carving tool: its outer boundary becomes the window wall the winding
    // faces — circumscribe so the carved window contains the nominal window.
    TopoDS_Shape outerCyl =
        build_polygon_cylinder(d, e / 2.0, m_corePolygonSegments, /*circumscribed=*/true);
    outerCyl = translate_shape(outerCyl, 0.0, 0.0, zCenter - d / 2.0);

    TopoDS_Shape innerCyl = build_polygon_cylinder(d, f / 2.0, m_corePolygonSegments);
    innerCyl = translate_shape(innerCyl, 0.0, 0.0, zCenter - d / 2.0);

    return done(BRepAlgoAPI_Cut(outerCyl, innerCyl), "leaving the centre post in the winding window");
}

TopoDS_Shape ShapeP::applyExtras(const std::map<std::string, double>& dims,
                                 const TopoDS_Shape& piece) const {
    const double b = required(dims, "B");
    return translate_shape(piece, 0.0, 0.0, -b);
}

// Round-column gap cut: cylinder of radius F/2, height = gap length,
// centered at machining.coordinates (matches MVB.js ERShape.applyMachining).
TopoDS_Shape ShapeP::applyMachining(const TopoDS_Shape& piece,
                                    const MAS::Machining& machining,
                                    const std::map<std::string, double>& dims) const {
    const auto& coords = machining.get_coordinates();
    if (coords.size() < 2) return piece;
    double gapLength = machining.get_length();
    if (gapLength <= 0.0) return piece;

    double xCoord = coords[0];
    double yCoord = coords[1];
    if (std::abs(xCoord) > 1e-12) {
        // SIDE-COLUMN GAP. This used to `return piece` with the note that a P-family core has
        // "no outer column to cut". PQ, RM and PM all derive from this class and all have outer
        // legs, and MKF emits one Machining per leg for an all-legs-ground core -- so every one
        // of those gaps was silently dropped and the FEM saw a core gapped on the centre post
        // only. Measured 2026-09-05 on PQ 20/16 (3 x 5 um): centre cut 0.1452 mm3, each lateral
        // 0.00000 mm3, and the core reached the mesher missing two thirds of its reluctance.
        // The knife: everything OUTSIDE the winding window (r > E/2) on this gap's side of the
        // core, one gap slab thick at the parting plane. That is the leg for PQ/RM/PM and the
        // skirt for a pot core -- whatever the outer structure is, the window cylinder is the
        // one boundary all of them share. Half-space per gap so each Machining cuts ITS leg and
        // the second cut is not a no-op. The knife overhangs every free face (1 mm) and extends
        // AWAY from the piece in y: nothing is there to over-cut, and coincident faces are what
        // makes OCCT's boolean drop a cut silently.
        // THE KNIFE IS SIZED FROM THE PIECE ITSELF (ABT #1720). It used to be max(A, C)/2: right
        // for PQ/RM/PM, whose footprint is A x C, but a POT core is round -- its footprint is A
        // and the MAS shape carries no C (P 59/36: A,B,D,E,F,G,H,M,r1) -- so every pot core with
        // side gaps threw "needs dimensions A, C and E". That is every pot core OMFEM meshes
        // without a functional gap: retype_residual_gaps turns MKF's 5 um residual gaps
        // (centre + one per side) into machined ones. The piece's bounding box is the outer
        // extent of whatever the family is, so no dimension is read for it at all.
        const double e = required(dims, "E");   // the winding window diameter
        Bnd_Box bb; BRepBndLib::Add(piece, bb);
        if (bb.IsVoid()) throw std::runtime_error("applyMachining: P-family side gap on an empty piece");
        double xmin, ymin, zmin, xmax, ymax, zmax; bb.Get(xmin, ymin, zmin, xmax, ymax, zmax);
        const double over = 1e-3;
        const double R = std::max({std::abs(xmin), std::abs(xmax), std::abs(zmin), std::abs(zmax)}) + over;   // past every outer face
        TopoDS_Shape knife = makeBox(R, gapLength + over, 2.0 * R);  // x in [-R/2, R/2] before shift
        knife = translate_shape(knife, (xCoord > 0 ? R : -R) / 2.0,  // x in [0, R] on the gap's side
                                yCoord + (yCoord > 0 ? -over : over) / 2.0, 0.0);
        gp_Ax2 winAxis(gp_Pnt(0.0, yCoord - gapLength - over, 0.0), gp_Dir(0, 1, 0));
        TopoDS_Shape window = BRepPrimAPI_MakeCylinder(winAxis, e / 2.0, 2.0 * (gapLength + over)).Shape();
        BRepAlgoAPI_Cut mk(knife, window);
        if (!mk.IsDone()) throw std::runtime_error("applyMachining: side knife minus window cylinder failed");
        BRepAlgoAPI_Cut cutter(piece, mk.Shape());
        if (!cutter.IsDone())
            throw std::runtime_error("applyMachining: the side gap cut did not complete (OCCT boolean failed) at x=" +
                                     std::to_string(xCoord) + " y=" + std::to_string(yCoord) + " length=" + std::to_string(gapLength));
        return cutter.Shape();
    }

    // Centre-column gap: cylinder along Y (the column axis after the
    // intrinsic -90°X rotation), radius F/2, height gapLength, centred at
    // Y = yCoord. Using default BRepPrimAPI_MakeCylinder orients along Z
    // and cuts a transverse slot through the column — same mistake Python
    // avoids by its coord system and ShapeEr handles via an explicit axis.
    const double f = required(dims, "F");
    gp_Ax2 cylAxis(gp_Pnt(0.0, yCoord - gapLength / 2.0, 0.0), gp_Dir(0, 1, 0));
    TopoDS_Shape tool = BRepPrimAPI_MakeCylinder(cylAxis, f / 2.0, gapLength).Shape();
    BRepAlgoAPI_Cut cutter(piece, tool);
    if (!cutter.IsDone()) throw std::runtime_error("applyMachining: the gap cut did not complete (OCCT boolean failed) at x=" + std::to_string(coords[0]) + " y=" + std::to_string(coords[1]) + " length=" + std::to_string(gapLength));
    return cutter.Shape();
}

} // namespace shapes
} // namespace mvb
