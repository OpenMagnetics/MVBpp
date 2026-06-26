#include "mvb/BobbinBuilder.h"
#include "mvb/Utils.h"
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <algorithm>
#include <cmath>
#include <numbers>

namespace mvb {

// Build a rectangular prism (w x d, extruded `h` along +Z from `zBottom`) whose
// four vertical corners are rounded with radius `r` (true circular arcs, so STEP
// stays clean). `r <= 0` (or a degenerate radius) falls back to a plain box, so
// callers can pass r=0 for non-rounded shapes. `r` is clamped to < min(w,d)/2.
static TopoDS_Shape makeRoundedRectPrism(double w, double d, double h,
                                         double r, double zBottom) {
    const double hw = w / 2.0, hd = d / 2.0;
    const double rmax = std::min(hw, hd);
    if (r > rmax - 1e-9) r = rmax - 1e-9;
    if (r <= 1e-9 || hw <= 0.0 || hd <= 0.0) {
        gp_Pnt corner(-hw, -hd, zBottom);
        return BRepPrimAPI_MakeBox(corner, w, d, h).Shape();
    }
    const double pi = std::numbers::pi;
    auto line = [](gp_Pnt a, gp_Pnt b) { return BRepBuilderAPI_MakeEdge(a, b).Edge(); };
    auto arc = [&](double cx, double cy, double a1, double a2) {
        gp_Circ c(gp_Ax2(gp_Pnt(cx, cy, 0.0), gp_Dir(0, 0, 1)), r);
        return BRepBuilderAPI_MakeEdge(c, a1, a2).Edge();
    };
    BRepBuilderAPI_MakeWire wire;
    wire.Add(line(gp_Pnt(hw, -(hd - r), 0), gp_Pnt(hw, hd - r, 0)));    // right
    wire.Add(arc(hw - r, hd - r, 0.0, pi / 2));                          // top-right
    wire.Add(line(gp_Pnt(hw - r, hd, 0), gp_Pnt(-(hw - r), hd, 0)));     // top
    wire.Add(arc(-(hw - r), hd - r, pi / 2, pi));                        // top-left
    wire.Add(line(gp_Pnt(-hw, hd - r, 0), gp_Pnt(-hw, -(hd - r), 0)));   // left
    wire.Add(arc(-(hw - r), -(hd - r), pi, 3 * pi / 2));                 // bottom-left
    wire.Add(line(gp_Pnt(-(hw - r), -hd, 0), gp_Pnt(hw - r, -hd, 0)));   // bottom
    wire.Add(arc(hw - r, -(hd - r), 3 * pi / 2, 2 * pi));                // bottom-right
    TopoDS_Face face = BRepBuilderAPI_MakeFace(wire.Wire()).Face();
    TopoDS_Shape prism = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, h)).Shape();
    return translate_shape(prism, 0.0, 0.0, zBottom);
}

static void checkIsDone(const BRepAlgoAPI_Cut& op, const char* ctx) {
    if (!op.IsDone()) throw std::runtime_error(std::string("BobbinBuilder: Cut failed: ") + ctx);
}
static void checkIsDone(const BRepAlgoAPI_Fuse& op, const char* ctx) {
    if (!op.IsDone()) throw std::runtime_error(std::string("BobbinBuilder: Fuse failed: ") + ctx);
}

TopoDS_Shape BobbinBuilder::buildBobbin(const MAS::CoreBobbinProcessedDescription& bobbin, double flangeThickness, bool axisIsY, int polygonSegments) {
    if (flangeThickness < 0.0 || std::isnan(flangeThickness)) {
        throw std::invalid_argument("BobbinBuilder: flangeThickness is invalid");
    }
    bool hasFlanges = flangeThickness > 0.0;

    double colWidth = bobbin.get_column_width().value_or(0.0);
    double colDepth = bobbin.get_column_depth();
    double wallThickness = bobbin.get_wall_thickness();
    if (wallThickness < 0.0 || std::isnan(wallThickness)) {
        throw std::invalid_argument("BobbinBuilder: wallThickness is invalid");
    }

    if (colWidth <= 0.0 || colDepth <= 0.0) {
        throw std::invalid_argument("BobbinBuilder: column_width and column_depth must be > 0");
    }

    const auto& wwList = bobbin.get_winding_windows();
    double wwWidth = 0.0;
    double wwHeight = 0.0;
    double wwDepth = 0.0;
    if (!wwList.empty()) {
        wwWidth = wwList[0].get_width().value_or(0.0);
        wwHeight = wwList[0].get_height().value_or(0.0);
        wwDepth = wwList[0].get_radial_height().value_or(wwWidth);
    }
    if (wwWidth <= 0.0) wwWidth = colWidth;
    if (wwHeight <= 0.0) wwHeight = colDepth;
    if (wwDepth <= 0.0) wwDepth = wwWidth;

    double height = wwHeight;
    if (height <= 0.0) return TopoDS_Shape();

    TopoDS_Shape bobbinShape;

    if (bobbin.get_column_shape() == MAS::ColumnShape::ROUND) {
        double outerR = colWidth;
        double innerR = colWidth - bobbin.get_column_thickness();
        if (innerR <= 0.0) innerR = outerR - wallThickness;

        // Build the cylinder along +Z then translate down by height/2 so it
        // straddles z=0 (matches the legacy MakeCylinder placement).
        auto make_cyl = [&](double radius, double h, double zBottom) -> TopoDS_Shape {
            TopoDS_Shape c = build_polygon_cylinder(h, radius, polygonSegments);
            return translate_shape(c, 0.0, 0.0, zBottom);
        };

        TopoDS_Shape outerCyl = make_cyl(outerR, height, -height / 2.0);

        TopoDS_Shape body;
        if (innerR + 1e-9 >= outerR) {
            body = outerCyl;
        } else {
            double innerHeight = height * 1.1;
            TopoDS_Shape innerCyl = make_cyl(innerR, innerHeight, -innerHeight / 2.0);
            BRepAlgoAPI_Cut cutOp(outerCyl, innerCyl);
            checkIsDone(cutOp, "round body hollow");
            body = cutOp.Shape();
        }

        if (hasFlanges) {
            double flangeOuterR = outerR + wwWidth;
            double holeHeight = flangeThickness * 1.2;
            double holeMargin = (holeHeight - flangeThickness) / 2.0;

            TopoDS_Shape topFlangeDisc = make_cyl(flangeOuterR, flangeThickness, height / 2.0);
            TopoDS_Shape topHole = make_cyl(innerR, holeHeight, height / 2.0 - holeMargin);
            BRepAlgoAPI_Cut topCut(topFlangeDisc, topHole);
            checkIsDone(topCut, "round top flange");

            TopoDS_Shape bottomFlangeDisc = make_cyl(flangeOuterR, flangeThickness,
                                                     -(height / 2.0 + flangeThickness));
            TopoDS_Shape bottomHole = make_cyl(innerR, holeHeight,
                                               -(height / 2.0 + flangeThickness) - holeMargin);
            BRepAlgoAPI_Cut bottomCut(bottomFlangeDisc, bottomHole);
            checkIsDone(bottomCut, "round bottom flange");

            // Skip the Fuse: return body + flanges as a compound of 3 solids.
            // Downstream cut_bobbin can then AABB-prefilter tools per solid;
            // most turns sit in the winding window and don't intersect the
            // flange AABBs, so per-solid cuts run against far fewer tools.
            TopoDS_Compound comp;
            BRep_Builder bld;
            bld.MakeCompound(comp);
            bld.Add(comp, body);
            bld.Add(comp, topCut.Shape());
            bld.Add(comp, bottomCut.Shape());
            bobbinShape = comp;
        } else {
            bobbinShape = body;
        }
    } else {
        // Rectangular (or irregular = EFD) column.
        double columnThickness = bobbin.get_column_thickness();
        double outerWidth = 2.0 * colWidth;
        double outerDepth = 2.0 * colDepth;
        double holeWidth = 2.0 * std::max(0.0, colWidth - columnThickness);
        double holeDepth = 2.0 * std::max(0.0, colDepth - columnThickness);
        if (holeWidth <= 0.0) holeWidth = outerWidth * 0.8;
        if (holeDepth <= 0.0) holeDepth = outerDepth * 0.8;
        double eps = wallThickness > 0 ? wallThickness * 0.1 : 1e-6;

        // Rounded corners (radius = 1/4 of the smaller yoke dimension): real U/E
        // bobbins are moulded with a corner radius on the OUTER winding surface
        // rather than a sharp 90°. The inner bore (and the flange holes) stay
        // SQUARE so they still match the square core central column — rounding
        // the bore would push bobbin material into the column corners and
        // overlap it. Only for true RECTANGULAR columns; round/oblong/irregular
        // pass r=0 (box).
        const bool roundCorners = (bobbin.get_column_shape() == MAS::ColumnShape::RECTANGULAR);
        const double rOuter = roundCorners ? 0.25 * std::min(outerWidth, outerDepth) : 0.0;
        const double rInner = 0.0;  // bore matches the (square) core column

        TopoDS_Shape outer = makeRoundedRectPrism(outerWidth, outerDepth, height, rOuter, -height / 2.0);
        TopoDS_Shape central = makeRoundedRectPrism(holeWidth, holeDepth, height + 2.0 * eps, rInner, -height / 2.0 - eps);
        BRepAlgoAPI_Cut bodyCut(outer, central);
        checkIsDone(bodyCut, "rect body hollow");

        if (hasFlanges) {
            double flangeWidth = 2.0 * (colWidth + wwWidth);
            double flangeDepth = 2.0 * (colDepth + wwWidth);
            double holeEps = flangeThickness * 0.1;
            const double rFlange = roundCorners ? 0.25 * std::min(flangeWidth, flangeDepth) : 0.0;

            TopoDS_Shape topFlange = makeRoundedRectPrism(flangeWidth, flangeDepth, flangeThickness, rFlange, height / 2.0);
            TopoDS_Shape topHole = makeRoundedRectPrism(holeWidth, holeDepth, flangeThickness + 2.0 * holeEps, rInner, height / 2.0 - holeEps);
            BRepAlgoAPI_Cut topCut(topFlange, topHole);
            checkIsDone(topCut, "rect top flange");

            TopoDS_Shape bottomFlange = makeRoundedRectPrism(flangeWidth, flangeDepth, flangeThickness, rFlange, -(height / 2.0 + flangeThickness));
            TopoDS_Shape bottomHole = makeRoundedRectPrism(holeWidth, holeDepth, flangeThickness + 2.0 * holeEps, rInner, -(height / 2.0 + flangeThickness) - holeEps);
            BRepAlgoAPI_Cut bottomCut(bottomFlange, bottomHole);
            checkIsDone(bottomCut, "rect bottom flange");

            // Skip the Fuse: compound body + flanges as 3 solids (see round path).
            TopoDS_Compound comp;
            BRep_Builder bld;
            bld.MakeCompound(comp);
            bld.Add(comp, bodyCut.Shape());
            bld.Add(comp, topCut.Shape());
            bld.Add(comp, bottomCut.Shape());
            bobbinShape = comp;
        } else {
            bobbinShape = bodyCut.Shape();
        }
    }

    if (axisIsY && !bobbinShape.IsNull()) {
        gp_Trsf rot;
        rot.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0)), -std::numbers::pi / 2.0);
        bobbinShape = BRepBuilderAPI_Transform(bobbinShape, rot).Shape();
    }

    return bobbinShape;
}

} // namespace mvb
