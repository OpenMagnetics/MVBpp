#include <stdexcept>
#include <algorithm>
#include "mvb/shapes/ShapeP.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

TopoDS_Face ShapeP::buildProfile(const std::map<std::string, double>& dims) const {
    double a = 0.0;
    auto it = dims.find("A");
    if (it != dims.end()) a = it->second / 2.0;

    TopoDS_Wire outer = build_polygon_circle(a, m_corePolygonSegments);
    BRepBuilderAPI_MakeFace face(outer);
    return face.Face();
}

TopoDS_Shape ShapeP::buildWindingWindow(const std::map<std::string, double>& dims) const {
    double b = 0.0, d = 0.0, e = 0.0, f = 0.0;
    auto it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;

    if (b == 0.0 || d == 0.0 || e == 0.0) return TopoDS_Shape();

    double zCenter = b - d / 2.0;

    // build_polygon_cylinder creates a base-at-zero prism; shift down by d/2 to center it at zCenter
    // Window-carving tool: its outer boundary becomes the window wall the winding
    // faces — circumscribe so the carved window contains the nominal window.
    TopoDS_Shape outerCyl =
        build_polygon_cylinder(d, e / 2.0, m_corePolygonSegments, /*circumscribed=*/true);
    outerCyl = translate_shape(outerCyl, 0.0, 0.0, zCenter - d / 2.0);

    TopoDS_Shape innerCyl = build_polygon_cylinder(d, f / 2.0, m_corePolygonSegments);
    innerCyl = translate_shape(innerCyl, 0.0, 0.0, zCenter - d / 2.0);

    BRepAlgoAPI_Cut cutter(outerCyl, innerCyl);
    return cutter.IsDone() ? cutter.Shape() : TopoDS_Shape();
}

TopoDS_Shape ShapeP::applyExtras(const std::map<std::string, double>& dims,
                                 const TopoDS_Shape& piece) const {
    double b = 0.0;
    auto it = dims.find("B");
    if (it != dims.end()) b = it->second;
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

    double f = 0.0;
    auto it = dims.find("F"); if (it != dims.end()) f = it->second;
    if (f <= 0.0) return piece;

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
        const double e = dims.count("E") ? dims.at("E") : 0.0;
        const double a = dims.count("A") ? dims.at("A") : 0.0;
        const double c = dims.count("C") ? dims.at("C") : 0.0;
        if (e <= 0.0 || a <= 0.0 || c <= 0.0)
            throw std::runtime_error("applyMachining: P-family side gap needs dimensions A, C and E");
        const double over = 1e-3;
        const double R = std::max(a, c) / 2.0 + over;              // past every outer face
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
    gp_Ax2 cylAxis(gp_Pnt(0.0, yCoord - gapLength / 2.0, 0.0), gp_Dir(0, 1, 0));
    TopoDS_Shape tool = BRepPrimAPI_MakeCylinder(cylAxis, f / 2.0, gapLength).Shape();
    BRepAlgoAPI_Cut cutter(piece, tool);
    if (!cutter.IsDone()) throw std::runtime_error("applyMachining: the gap cut did not complete (OCCT boolean failed) at x=" + std::to_string(coords[0]) + " y=" + std::to_string(coords[1]) + " length=" + std::to_string(gapLength));
    return cutter.Shape();
}

} // namespace shapes
} // namespace mvb
