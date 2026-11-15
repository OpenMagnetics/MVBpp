#pragma once

#include "mvb/shapes/ShapeBuilder.h"

namespace mvb {
namespace shapes {

// ROD core: a bare ferrite cylinder, no flanges and no return limb; the coil is wound straight
// on it and the circuit closes through the air. Dimension convention, MKF's CorePieceRod:
// A rod diameter, B rod length, H optional axial bore. C is not read (a rod is its own post).
class ShapeRod : public ShapeBuilder {
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
    // Bore, then centred on the origin along its axis, as ShapeDrum does.
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;
};

} // namespace shapes
} // namespace mvb
