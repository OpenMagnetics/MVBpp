#include "mvb/FR4Builder.h"

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>

namespace mvb {

namespace {
// An OCCT boolean that fails would otherwise leave the uncut solid: a board with no core cut-out,
// sitting in the ferrite. Refuse it by name instead (ABT #1729).
template <class Op>
TopoDS_Shape booleanResult(Op& op, const char* what)
{
    if (!op.IsDone()) {
        throw std::runtime_error(std::string("FR4Builder::buildFR4Board: OCCT failed to ") + what +
                                 "; the board would keep no core cut-out");
    }
    return op.Shape();
}
}  // namespace

double FR4Builder::coreCutoutClearance(const MAS::Group& group, double columnFace)
{
    const auto& coords = group.get_coordinates();
    const auto& dims   = group.get_dimensions();
    if (coords.empty() || dims.empty()) {
        throw std::runtime_error(
            "FR4Builder: the printed group has no coordinates/dimensions, so the core cut-out "
            "of its board cannot be placed");
    }
    // The copper region (the group) starts where MKF placed it: coreToLayerDistance (or the
    // PCB's coreToTrack) out from the column face. That inset IS the board's clearance to the
    // column when the group carries no PCB description; it is never invented here.
    const double groupInnerEdge = coords[0] - dims[0] / 2.0;
    const double groupInset     = groupInnerEdge - columnFace;
    if (groupInset < 0.0) {
        throw std::runtime_error(
            "FR4Builder: the printed group's inner edge (" + std::to_string(groupInnerEdge * 1e3) +
            " mm) lies inside the central column face (" + std::to_string(columnFace * 1e3) +
            " mm): its copper would sit in the ferrite");
    }
    const auto& pcb = group.get_pcb();
    if (!pcb) return groupInset;

    // MAS pcbOutline.coreCutoutClearance: "Gap between the core column faces and the cut-out
    // edges. Defaults to coreToTrack when absent (stated rule)."
    const auto stated = pcb->get_outline().get_core_cutout_clearance();
    const double clearance = stated ? *stated : pcb->get_design_rules().get_core_to_track();
    if (clearance > groupInset) {
        throw std::runtime_error(
            "FR4Builder: the PCB's core cut-out clearance (" + std::to_string(clearance * 1e3) +
            " mm) puts the board edge beyond the printed group's inner edge (" +
            std::to_string(groupInset * 1e3) + " mm from the column face): the copper would hang "
            "off the board");
    }
    return clearance;
}

TopoDS_Shape FR4Builder::buildFR4Board(
    const std::vector<MAS::Group>& groups,
    const MAS::CoreBobbinProcessedDescription& bobbinPd)
{
    if (groups.empty()) return TopoDS_Shape();

    const auto& group0 = groups.front();
    if (group0.get_type() != MAS::WiringTechnology::PRINTED) return TopoDS_Shape();

    const auto& groupCoords = group0.get_coordinates();
    const auto& groupDims   = group0.get_dimensions();
    if (groupCoords.size() < 2 || groupDims.size() < 2) return TopoDS_Shape();

    const double groupX      = groupCoords[0];                // radial
    const double groupZ      = (groupCoords.size() > 1) ? groupCoords[1] : 0.0; // height
    const double groupWidth  = groupDims[0];
    const double groupHeight = groupDims[1];

    // Widest layer across all groups — used as the Y (depth) extent of
    // rectangular/oblong boards. Round boards size purely from group 0's
    // own dimensions.
    double maxLayerWidth = 0.0;
    for (const auto& g : groups) {
        const auto& d = g.get_dimensions();
        if (!d.empty()) maxLayerWidth = std::max(maxLayerWidth, d[0]);
    }

    // Bobbin column metrics — required, must be populated by the caller.
    if (!bobbinPd.get_column_width()) {
        throw std::runtime_error(
            "FR4Builder::buildFR4Board: bobbinPd.column_width is unset. "
            "The caller must populate the bobbin processed description "
            "(see MagneticBuilder::patchBobbinDimensions) before drawing "
            "the FR4 board.");
    }
    const double columnWidth          = *bobbinPd.get_column_width();
    const double columnDepth          = bobbinPd.get_column_depth();
    const MAS::ColumnShape columnShape = bobbinPd.get_column_shape();

    // Board thickness: the group height minus a 4 µm clearance so the FR4 plane doesn't
    // Z-fight with the turn planes. No thickness floor: a floor thicker than the group pushed
    // the board out of the winding window into the yokes (ABT #1729).
    const double zFightClearance = 4e-6;
    const double fr4Thickness    = groupHeight - zFightClearance;
    if (fr4Thickness <= 0.0) {
        throw std::runtime_error(
            "FR4Builder::buildFR4Board: the printed group is " + std::to_string(groupHeight * 1e6) +
            " um high, too thin to hold a board");
    }

    // The board's outer edge is the group's own outer edge (the group is the board footprint,
    // already inset from the lateral legs). The core cut-out is the column face plus the
    // board-to-column clearance from the data, so the board never reaches into the ferrite
    // (ABT #1729: the hole used to be the group edge minus an invented 250 µm, which put the
    // board 0.15 mm inside a 5.1 mm column face and recessed the column in every planar mesh).
    const double outerRadius = groupX + groupWidth / 2.0;
    const double clearance   = coreCutoutClearance(group0, columnWidth);

    // Y-depth of the board: use bobbin column depth + widest layer width.
    const double groupDepth = columnDepth + maxLayerWidth;
    if (!(groupDepth > 0.0)) {
        throw std::runtime_error(
            "FR4Builder::buildFR4Board: column depth + widest layer width is not positive; the "
            "board depth is undefined");
    }

    TopoDS_Shape board;

    if (columnShape == MAS::ColumnShape::ROUND) {
        // Annulus: outer disc minus inner hole cylinder, oriented along Z,
        // centred at groupZ.
        gp_Ax2 axOuter(gp_Pnt(0.0, 0.0, groupZ - fr4Thickness / 2.0), gp_Dir(0, 0, 1));
        TopoDS_Shape outerDisc = BRepPrimAPI_MakeCylinder(axOuter, outerRadius, fr4Thickness).Shape();

        // Hole extends slightly beyond so the cut is clean.
        const double holeExtra = 2e-3;
        gp_Ax2 axHole(gp_Pnt(0.0, 0.0, groupZ - fr4Thickness / 2.0 - holeExtra / 2.0), gp_Dir(0, 0, 1));
        TopoDS_Shape hole = BRepPrimAPI_MakeCylinder(axHole, columnWidth + clearance,
                                                       fr4Thickness + holeExtra).Shape();
        BRepAlgoAPI_Cut cutter(outerDisc, hole);
        board = booleanResult(cutter, "cut the round column hole");
    }
    else if (columnShape == MAS::ColumnShape::OBLONG) {
        // Rectangular plate, stadium-shaped hole.
        const double boardW = outerRadius * 2.0;
        const double boardD = groupDepth * 2.0;
        gp_Pnt boardCorner(-boardW / 2.0, -boardD / 2.0, groupZ - fr4Thickness / 2.0);
        TopoDS_Shape plate = BRepPrimAPI_MakeBox(boardCorner, boardW, boardD, fr4Thickness).Shape();

        const double fullColDepth = columnDepth * 2.0;
        const double fullColWidth = columnWidth * 2.0;
        const double holeRadius   = columnWidth + clearance;
        const double straightLen  = std::max(fullColDepth - fullColWidth, 0.0);
        const double holeExtra    = 2e-3;
        const double holeThickness = fr4Thickness + holeExtra;

        // Central rectangle
        gp_Pnt rectCorner(-holeRadius, -straightLen / 2.0,
                          groupZ - holeThickness / 2.0);
        TopoDS_Shape rect = BRepPrimAPI_MakeBox(rectCorner, 2.0 * holeRadius,
                                                  straightLen, holeThickness).Shape();

        // Two semicircle caps
        gp_Ax2 axTop(gp_Pnt(0.0, straightLen / 2.0, groupZ - holeThickness / 2.0),
                     gp_Dir(0, 0, 1));
        TopoDS_Shape capTop = BRepPrimAPI_MakeCylinder(axTop, holeRadius, holeThickness).Shape();
        gp_Ax2 axBot(gp_Pnt(0.0, -straightLen / 2.0, groupZ - holeThickness / 2.0),
                     gp_Dir(0, 0, 1));
        TopoDS_Shape capBot = BRepPrimAPI_MakeCylinder(axBot, holeRadius, holeThickness).Shape();

        TopoDS_Shape hole = rect;
        { BRepAlgoAPI_Fuse f(hole, capTop); hole = booleanResult(f, "fuse the oblong hole's top cap"); }
        { BRepAlgoAPI_Fuse f(hole, capBot); hole = booleanResult(f, "fuse the oblong hole's bottom cap"); }

        BRepAlgoAPI_Cut cutter(plate, hole);
        board = booleanResult(cutter, "cut the oblong column hole");
    }
    else {
        // Rectangular column → rectangular board with rectangular hole.
        const double boardW = outerRadius * 2.0;
        const double boardD = groupDepth * 2.0;
        gp_Pnt boardCorner(-boardW / 2.0, -boardD / 2.0, groupZ - fr4Thickness / 2.0);
        TopoDS_Shape plate = BRepPrimAPI_MakeBox(boardCorner, boardW, boardD, fr4Thickness).Shape();

        const double holeW = (columnWidth + clearance) * 2.0;
        const double holeD = (columnDepth + clearance) * 2.0;
        const double holeExtra    = 2e-3;
        const double holeThickness = fr4Thickness + holeExtra;

        if (!(holeW > 0.0 && holeD > 0.0)) {
            throw std::runtime_error(
                "FR4Builder::buildFR4Board: the rectangular column hole is not positive (" +
                std::to_string(holeW * 1e3) + " x " + std::to_string(holeD * 1e3) +
                " mm); the column dimensions are invalid");
        }
        gp_Pnt holeCorner(-holeW / 2.0, -holeD / 2.0,
                          groupZ - holeThickness / 2.0);
        TopoDS_Shape hole = BRepPrimAPI_MakeBox(holeCorner, holeW, holeD, holeThickness).Shape();
        BRepAlgoAPI_Cut cutter(plate, hole);
        board = booleanResult(cutter, "cut the rectangular column hole");
    }

    // The board is built with the column axis along Z (flat face in XY).
    // The rest of the magnetic — bobbin (BobbinBuilder.cpp:172), turns,
    // spacers — uses Y as the column axis. Rotate -π/2 around X to bring
    // the board into the Y-axis convention so its flat face lands in XZ.
    if (!board.IsNull()) {
        gp_Trsf rot;
        rot.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)),
                        -std::numbers::pi / 2.0);
        board = BRepBuilderAPI_Transform(board, rot).Shape();
    }

    return board;
}

} // namespace mvb
