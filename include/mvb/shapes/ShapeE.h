#pragma once

#include "mvb/shapes/ShapeBuilder.h"

namespace mvb {
namespace shapes {

class ShapeE : public ShapeBuilder {
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;

private:
    // The record describes an EFD pole: shallower than it is wide (F2 < F), or chamfered (q).
    bool hasEfdPole(const std::map<std::string, double>& dims) const;
};

} // namespace shapes
} // namespace mvb
