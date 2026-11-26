#pragma once

#include "mvb/shapes/ShapeE.h"

namespace mvb {
namespace shapes {

// EPC (TDK's low-profile series). An EPC is not an EP: seen on its mating face it is an E core, a
// back plate, two outer legs and a centre pole whose section is a STADIUM, F long along the
// width A and F2 thick across the depth C, with end discs of diameter F2. That is the piece MKF
// models and winds (CorePieceEpc: a rectangular column F x F2 with corner radius F2 / 2, centred
// on the axis), so it is the piece drawn. G, the opening between the ears above the window, is
// recorded and not modelled by MKF, and is not drawn either (ABT #2586).
//
// The gap cutters are ShapeBuilder's: their centre tool is F x gap x F2 when the record carries
// F2, which is the pole's bounding box.
class ShapeEpc : public ShapeE {
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;
};

} // namespace shapes
} // namespace mvb
