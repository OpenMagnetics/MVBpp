#pragma once

#include "MAS.hpp"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
#include <optional>
#include <utility>
#include <vector>

namespace mvb {

// The CONDUCTING cross-section of a MAS wire, {width, height} in metres: round -> {d, d};
// rectangular/planar -> {w, h}; LITZ -> the bare bundle treated as one solid conductor, its
// diameter from MKF's Wire::get_outer_diameter_bare_litz (capped by the catalogue outer
// diameter). The single source of truth for what MVB++ draws as copper and what OMFEM's
// level-set mesher carves; throws (no silent fallback) when the wire lacks the data.
std::pair<double, double> conducting_dimensions(const MAS::Wire& wire);

class TurnBuilder {
public:
    // Context for a turn wound around a NON-MAIN core column (MAS multi-column placement:
    // section/winding `windingWindow` -> winding window `column` edge -> columnData). The turn
    // loop is built with the same shape builders but RELATIVE to this column's axis, then
    // translated to the leg position. Half dimensions include the lateral bobbin wall
    // (columnThickness), mirroring how the main-column path uses the bobbin's columnWidth
    // (core half-width + wall). Absent spec = main column at the origin, byte-identical
    // legacy behaviour.
    struct WoundColumnSpec {
        double axisX = 0.0;        // MAS x of the wound column axis
        double halfWidth = 0.0;    // column width/2 + bobbin columnThickness
        double halfDepth = 0.0;    // column depth/2 + bobbin columnThickness
        MAS::ColumnShape shape = MAS::ColumnShape::RECTANGULAR;
    };

    // GLOBAL setting: fuse each turn's sub-solid compound (straight segments + corner pieces) into
    // ONE solid. A round/rect turn is built as a COMPOUND of overlapping cylinders/boxes + corner
    // tori/swept-rects -- fast, and perfectly fine for STL/STEP VISUALISATION (the web frontend), so
    // the default is OFF. But the overlapping sub-solids make a tetrahedral MESH fail ("overlapping
    // facets"), so the FEM meshing path must enable fusing to get a single clean conductor per turn.
    // Fusing is a boolean per turn (slower) -- hence the toggle. Also honoured via env MVB_FUSE_TURNS.
    static void setFuseTurnParts(bool on);
    static bool fuseTurnParts();

    // paintCoating: true  → turn solid drawn at the OUTER (insulation) footprint
    //                       (default; unchanged behaviour for visualisation).
    //               false → turn solid drawn at the CONDUCTING (copper) footprint
    //                       — what FEM winding-loss extraction must mesh. For
    //                       LITZ this is the bare bundle treated as a single
    //                       solid conductor (diameter from MKF).
    // woundColumn: present when the turn wraps a non-main column (see WoundColumnSpec).
    static TopoDS_Shape buildTurn(const MAS::Turn& turn,
                                  const MAS::Wire& wire,
                                  const MAS::CoreBobbinProcessedDescription& bobbin,
                                  bool isToroidal,
                                  int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                  int wireRevolutionSegments = DEFAULT_WIRE_REVOLUTION_SEGMENTS,
                                  bool paintCoating = true,
                                  const std::optional<WoundColumnSpec>& woundColumn = std::nullopt);

    // Build a turn using ONLY the data on the Turn itself: coordinates,
    // dimensions, cross_sectional_shape, additional_coordinates and
    // rotation. No Wire or Bobbin lookup. Used by drawTurns/drawWinding
    // when the caller does not provide (or want) that surrounding context.
    //
    // Topology heuristic:
    //   - additional_coordinates present → toroidal racetrack turn
    //   - otherwise                       → concentric round-column turn
    //                                       (turn radius = |coordinates|)
    //
    // Throws `std::runtime_error` if any required field is missing.
    // paintCoating must stay true here: a standalone Turn carries only its
    // outer footprint, so requesting the conducting cross-section throws.
    static TopoDS_Shape buildFromTurnAlone(const MAS::Turn& turn,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                            bool paintCoating = true);

    // Resolve the wire cross-section {width, height} (outer footprint or copper per
    // paintCoating). Used by ConductorBuilder, which sweeps whole conductors instead of
    // per-turn solids. Prefer the Turn overload: turn.dimensions carries the OUTER
    // footprint, whereas the context-less overload can only fall back to the wire's own
    // diameters (a bare round wire then resolves to its CONDUCTING diameter — too thin
    // for paintCoating=true). Throws on missing wire data (no fallbacks).
    static std::pair<double, double> wireDimensions(const MAS::Wire& wire, bool paintCoating);
    static std::pair<double, double> wireDimensions(const MAS::Wire& wire,
                                                    const MAS::Turn& turn,
                                                    bool paintCoating);

    // Fuse all SOLIDs of a compound into one — the same routine the per-turn fuse toggle
    // uses; exposed for ConductorBuilder's per-conductor compounds.
    static TopoDS_Shape fuseSolids(const TopoDS_Shape& compound);

    static void clearCache();
};

} // namespace mvb
