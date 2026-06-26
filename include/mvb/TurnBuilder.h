#pragma once

#include "MAS.hpp"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
#include <vector>

namespace mvb {

class TurnBuilder {
public:
    // paintCoating: true  → turn solid drawn at the OUTER (insulation) footprint
    //                       (default; unchanged behaviour for visualisation).
    //               false → turn solid drawn at the CONDUCTING (copper) footprint
    //                       — what FEM winding-loss extraction must mesh. For
    //                       LITZ this is the bare bundle treated as a single
    //                       solid conductor (diameter from MKF).
    static TopoDS_Shape buildTurn(const MAS::Turn& turn,
                                  const MAS::Wire& wire,
                                  const MAS::CoreBobbinProcessedDescription& bobbin,
                                  bool isToroidal,
                                  int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                  int wireRevolutionSegments = DEFAULT_WIRE_REVOLUTION_SEGMENTS,
                                  bool paintCoating = true);

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

    static void clearCache();
};

} // namespace mvb
