#pragma once

#include "mvb/shapes/ShapeE.h"

namespace mvb {
namespace shapes {

class ShapeEp : public ShapeE {
public:
    // EP center column is ROUND (diameter F). The generic gap cutter sizes its depth to
    // C (a rectangular-column assumption), which for EP is far larger than F and so the
    // box clips the closed pot back wall / column corner at the -Z face. Override to size
    // the cutter to the round column (depth F).
    TopoDS_Shape applyMachining(const TopoDS_Shape& piece,
                                const MAS::Machining& machining,
                                const std::map<std::string, double>& dims) const override;
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;
};

} // namespace shapes
} // namespace mvb
