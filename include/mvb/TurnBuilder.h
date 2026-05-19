#pragma once

#include "MAS.hpp"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
#include <vector>

namespace mvb {

class TurnBuilder {
public:
    static TopoDS_Shape buildTurn(const MAS::Turn& turn,
                                  const MAS::Wire& wire,
                                  const MAS::CoreBobbinProcessedDescription& bobbin,
                                  bool isToroidal,
                                  int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                  int wireRevolutionSegments = DEFAULT_WIRE_REVOLUTION_SEGMENTS);

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
    static TopoDS_Shape buildFromTurnAlone(const MAS::Turn& turn,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS);

    static void clearCache();
};

} // namespace mvb
