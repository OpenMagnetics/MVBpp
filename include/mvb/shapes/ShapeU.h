#pragma once

#include "mvb/shapes/ShapeBuilder.h"

namespace mvb {
namespace shapes {

class ShapeU : public ShapeBuilder {
protected:
    /// The two legs of a U: the one that is wound (centred on x = 0) and the one that returns the flux.
    struct Legs { double winding; double returning; };
    /// MKF's reading of A, E and H (CorePieceU, ABT #2587). Throws when A and E leave no leg.
    Legs legs(const std::map<std::string, double>& dims) const;

    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;
    TopoDS_Shape applyMachining(const TopoDS_Shape& piece,
                                const MAS::Machining& machining,
                                const std::map<std::string, double>& dims) const override;
};

} // namespace shapes
} // namespace mvb
