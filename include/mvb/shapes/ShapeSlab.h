#pragma once

#include "mvb/shapes/ShapeP.h"

namespace mvb {
namespace shapes {

// Slab cores DS / HS (MKF CorePieceSlab, ABT #263): a pot core with two sides cut away by
// flats. A = outline diameter, C = across the flats, B = piece height, D = window height,
// E = window outer diameter, F = post diameter, H (optional) = bore through the post (HS).
// Built as the P profile, window and post, then the two flats and the bore are cut before
// the piece is dropped to the parting plane. The flats are normal to profile Y, which the
// intrinsic -90 deg X rotation turns into the core's depth: width A, depth C, as MKF's
// CorePieceSlab reports them (set_width(A), set_depth(C)).
class ShapeSlab : public ShapeP {
protected:
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;
public:
    // The bore H through the post, when the shape has one (HS, RS). Shared by both RS halves.
    static TopoDS_Shape drillBore(const std::map<std::string, double>& dims, const TopoDS_Shape& piece);
};

// RS cores (MKF CorePieceRs, ABT #263 / #1443): ONE slab half paired with ONE plain pot-core
// round ("A slab piece can be paired with a standard pot core round to make an RS combination",
// Magnetics ferrite catalogue, slab cores). The RS shape's dimension table is the slab's: C (across
// the flats) and G (the wire opening) exist only on that piece, so the RS piece on its own is the
// slab, built exactly as DS/HS. The round half of a set is the same pot (A B D E F, bore H) without
// flats and without wire slots, which is how MKF's CorePieceRs models it
// (symmetric_set_constants(slab, flats=false)).
//
// Which half is which comes from MKF: Core::create_geometrical_description emits the TOP half
// (rotation {pi, pi, 0}) and the BOTTOM half (rotation {0, 0, 0}), and MKF's thermal model
// (Temperature.cpp calculateConcentricCoreExteriorAreas) books the flat-cut plate on the top yoke
// and the full round plate on the bottom yoke. So: top = slab, bottom = round. buildCoreShapes
// draws the bottom half with buildRoundHalf.
class ShapeRs : public ShapeSlab {
public:
    // The set's plain pot round: the same pot without the flats, drilled with H.
    TopoDS_Shape buildRoundHalf(const MAS::CoreShape& shapeData) const;
};

} // namespace shapes
} // namespace mvb
