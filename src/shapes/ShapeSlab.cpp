#include <stdexcept>
#include <string>
#include "mvb/shapes/ShapeSlab.h"
#include "mvb/Utils.h"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>

namespace mvb {
namespace shapes {

TopoDS_Shape ShapeSlab::applyExtras(const std::map<std::string, double>& dims,
                                    const TopoDS_Shape& piece) const {
    for (const char* letter : {"A", "B", "C", "F"}) {
        auto it = dims.find(letter);
        if (it == dims.end() || !(it->second > 0.0))
            throw std::runtime_error(std::string("ShapeSlab: dimension ") + letter +
                                     " is missing or non-positive");
    }
    const double a = dims.at("A"), b = dims.at("B"), c = dims.at("C"), f = dims.at("F");
    if (c <= f || c > a)
        throw std::runtime_error("ShapeSlab: the flats C must clear the post F and lie within the outline A");

    // Two flats. Each knife starts exactly at |y| = C/2 and overhangs every other face by a
    // margin, so no knife face is coincident with a piece face (OCCT drops such cuts quietly).
    const double over = 1e-3;
    const double knifeDepth = a / 2.0 + over;                      // reaches past the outline
    TopoDS_Shape shape = piece;
    for (double sign : {1.0, -1.0}) {
        TopoDS_Shape knife = makeBox(a + 2.0 * over, knifeDepth, b + 2.0 * over);
        knife = translate_shape(knife, 0.0, sign * (c / 2.0 + knifeDepth / 2.0), b / 2.0);
        BRepAlgoAPI_Cut cut(shape, knife);
        if (!cut.IsDone())
            throw std::runtime_error("ShapeSlab: cutting the flats failed (OCCT boolean)");
        shape = cut.Shape();
    }

    shape = drillBore(dims, shape);

    return ShapeP::applyExtras(dims, shape);
}

TopoDS_Shape ShapeSlab::drillBore(const std::map<std::string, double>& dims,
                                  const TopoDS_Shape& piece) {
    auto h = dims.find("H");
    if (h == dims.end() || !(h->second > 0.0)) return piece;   // DS: solid post, no H
    const double b = dims.at("B"), f = dims.at("F");
    if (h->second >= f)
        throw std::runtime_error("ShapeSlab: bore H must be smaller than the post F");
    const double over = 1e-3;
    gp_Ax2 axis(gp_Pnt(0.0, 0.0, -over), gp_Dir(0, 0, 1));
    TopoDS_Shape bore = BRepPrimAPI_MakeCylinder(axis, h->second / 2.0, b + 2.0 * over).Shape();
    BRepAlgoAPI_Cut cut(piece, bore);
    if (!cut.IsDone())
        throw std::runtime_error("ShapeSlab: drilling the bore failed (OCCT boolean)");
    return cut.Shape();
}

namespace {

// The round half of an RS set: ShapeP's pot (outline A, window E x D, post F) plus the bore H.
// No wire slots: MKF's CorePieceRs models the round as CorePieceP with n = 0.
class RsRoundHalf : public ShapeP {
protected:
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override {
        for (const char* letter : {"A", "B", "D", "E", "F"}) {
            auto it = dims.find(letter);
            if (it == dims.end() || !(it->second > 0.0))
                throw std::runtime_error(std::string("ShapeRs: round half dimension ") + letter +
                                         " is missing or non-positive");
        }
        if (!(dims.at("F") < dims.at("E") && dims.at("E") < dims.at("A")))
            throw std::runtime_error("ShapeRs: round half needs F < E < A (post inside the window inside the outline)");
        if (!(dims.at("D") < dims.at("B")))
            throw std::runtime_error("ShapeRs: round half window height D must be smaller than the half height B");
        return ShapeP::applyExtras(dims, ShapeSlab::drillBore(dims, piece));
    }
};

} // namespace

TopoDS_Shape ShapeRs::buildRoundHalf(const MAS::CoreShape& shapeData) const {
    RsRoundHalf round;
    round.setCorePolygonSegments(corePolygonSegments());
    TopoDS_Shape piece = round.buildPiece(shapeData);
    if (piece.IsNull())
        throw std::runtime_error("ShapeRs: the round half of '" + shapeData.get_name().value_or(std::string("?")) +
                                 "' built no solid (the shape has no dimensions)");
    return piece;
}

} // namespace shapes
} // namespace mvb
