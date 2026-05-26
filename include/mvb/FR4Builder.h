#pragma once

#include "MAS.hpp"
#include <TopoDS_Shape.hxx>
#include <string>

namespace mvb {

// Builds the FR4 PCB board for a planar transformer coil. Mirrors the logic
// in MVB.js ReplicadBuilder::getFR4Board.
//
// Inputs are in MAS metres. The output is in metres as well (no implicit
// SCALE applied — callers apply `scale` uniformly in the WASM wrapper).
//
// Returns a null shape when:
//   - the coil has no groupsDescription, or
//   - the first group is not a PCB (WiringTechnology::PRINTED).
//
// Parameters:
//   borderToWireDistance : extra margin on the outer edge of the PCB
//   coreToLayerDistance  : margin between the column cut-out and the PCB
//                          inner edge
class FR4Builder {
public:
    static constexpr double DEFAULT_BORDER_TO_WIRE_DISTANCE = 90e-6;   // 90 µm
    static constexpr double DEFAULT_CORE_TO_LAYER_DISTANCE  = 250e-6;  // 250 µm
    static constexpr double MIN_FR4_THICKNESS               = 0.5e-3;  // 0.5 mm

    // The caller is responsible for handing in:
    //   - groups: the coil's `groupsDescription` (non-empty; first group's
    //             type must be PRINTED — else this returns a null shape)
    //   - bobbinPd: the bobbin's processed description, ALREADY PATCHED
    //               (column_shape / column_width / column_depth populated
    //               — typically by reading the core's processed centre
    //               column when the variant came from MKF as a string).
    // No internal fallbacks. Throws if column_width is unset on bobbinPd
    // for a planar group — the caller must populate it upstream.
    static TopoDS_Shape buildFR4Board(
        const std::vector<MAS::Group>& groups,
        const MAS::CoreBobbinProcessedDescription& bobbinPd,
        double borderToWireDistance = DEFAULT_BORDER_TO_WIRE_DISTANCE,
        double coreToLayerDistance  = DEFAULT_CORE_TO_LAYER_DISTANCE);
};

} // namespace mvb
