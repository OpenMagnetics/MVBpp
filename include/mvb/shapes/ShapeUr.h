#pragma once

#include "mvb/shapes/ShapeBuilder.h"

namespace mvb {
namespace shapes {

class ShapeUr : public ShapeBuilder {
public:
    TopoDS_Shape buildPiece(const MAS::CoreShape& shapeData) const override;
    TopoDS_Shape applyMachining(const TopoDS_Shape& piece,
                                const MAS::Machining& machining,
                                const std::map<std::string, double>& dims) const override;
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>&) const override {
        return TopoDS_Face();  // unused — buildPiece is fully overridden
    }
};

} // namespace shapes
} // namespace mvb
