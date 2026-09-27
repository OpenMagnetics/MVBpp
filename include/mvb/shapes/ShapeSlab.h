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
};

} // namespace shapes
} // namespace mvb
