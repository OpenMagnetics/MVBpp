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
// The core cut-out is the column face plus the board-to-column clearance from the data
// (coreCutoutClearance below); the board never reaches into the ferrite (ABT #1729).
class FR4Builder {
public:
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
        const MAS::CoreBobbinProcessedDescription& bobbinPd);

    // Gap between the central column face (`columnFace`, measured from the column axis) and the
    // board's cut-out edge. With a `pcb` on the group: pcb.outline.coreCutoutClearance, or
    // pcb.designRules.coreToTrack when absent (the MAS-stated rule). Without one: the group's own
    // inset from the column face -- the coreToLayerDistance MKF wound the group with. Throws when
    // the group's inner edge is inside the column, or when the PCB clearance would leave the
    // group's copper off the board.
    static double coreCutoutClearance(const MAS::Group& group, double columnFace);
};

} // namespace mvb
