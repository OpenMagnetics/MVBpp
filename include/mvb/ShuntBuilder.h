#pragma once

#include "MAS.hpp"
#include "mvb/NamedShape.h"
#include <TopoDS_Shape.hxx>
#include <string>
#include <vector>

namespace OpenMagnetics {
class Magnetic;
}

namespace mvb {

// Magnetic shunts (MAS-RFC 0015, ABT #1176 WP7): the permeable sheets listed in
// `magnetic.shunts[]`, drawn as solids of their own.
//
// CONTRACT WITH MKF. MKF owns the shunt model (`physical_models/MagneticShunt.*`) and its
// geometric validation; MVB++ draws only what MKF accepts, and re-derives nothing:
//   - `outsideWindow` throws (MagneticShuntModel::check_supported_placements);
//   - an `inWindow` / `betweenSections` sheet is validated by
//     MagneticShuntModel::extract_network_inputs: E-type core with three RECTANGULAR columns
//     (a round/oblong column throws NotImplementedException, "a sheet across a round or oblong
//     column is an annulus, not modelled" -- so no annulus is drawn either), the sheet inside a
//     column only through that column's gap, `gapToColumns` equal to the DRAWN clearance, and
//     `segments` adding up to the width inside one window;
//   - an `onColumn` sheet is validated by MagneticShuntModel::apply_shunts_to_gapping, the call
//     MKF's magnetizing inductance makes (it resolves the sheet material too).
// Every MKF exception propagates unchanged; nothing is drawn for a magnetic MKF refuses.
//
// GEOMETRY. MAS `coordinates` is the sheet CENTRE in the main-column frame and `dimensions` is
// [width along x, height along the column axis y, depth along z] -- the same frame MKF places
// the core halves, the gaps and the turns in, and the frame MVB++ draws them in, so no axis is
// swapped. A missing z coordinate is 0 (MKF reads it the same way).
// `gapToColumns` needs NO cut: MKF requires the box to be drawn with that clearance already
// (1 um / 1 % agreement) and throws otherwise, so the solid IS the `dimensions` box.
//
// SEGMENTS. `segments[k] = {length, gap}` is one piece followed by the gap after it. Pieces are
// laid along +x from the box's minimum-x face: piece k spans [x0 + sum_{j<k}(length_j + gap_j),
// that + length_k]. MKF has checked sum(length + gap) == width; a trailing gap is empty box.
//
// NAMES (decision, ABT #1176). `Shunt_<i>` with i the index in `magnetic.shunts[]`, and
// `Shunt_<i>_<k>` for piece k of a segmented sheet. The MAS `name` is NOT used: it is optional,
// free text and not unique, while the index is the one key a consumer can follow back to the
// shunt's record (placement, material). OMFEM's strict classify() maps the `Shunt_` prefix to
// the shunt region. Role::Shunt; `materialName` = the material string, or the `name` of the
// inline core-material record.
class ShuntBuilder {
public:
    // Validated, drawn shunts of the magnetic, in `shunts[]` order (segments in piece order).
    // Empty when the magnetic has no shunts.
    static std::vector<NamedShape> buildShuntsNamed(const OpenMagnetics::Magnetic& magnetic);
    static std::vector<NamedShape> buildShuntsNamed(const MAS::Magnetic& magnetic);
};

// Hook body for MagneticBuilder::appendAccessorySolids: appends the shunts to `all` and cuts
// every already-present Spacer solid the shunts overlap (a sheet held in an additive gap takes
// that part of the gap: MKF models the gap as g - t of air/shim plus t of sheet). A spacer cut
// that fails or does not remove exactly the overlap volume throws.
void appendShuntSolids(std::vector<NamedShape>& all, const OpenMagnetics::Magnetic& magnetic);
void appendShuntSolids(std::vector<NamedShape>& all, const MAS::Magnetic& magnetic);

// Take out of `host` the volume every Shunt solid in `all` shares with it, with the same strict
// check (a failed boolean, or a cut that removes anything but the shared volume, throws). Used for
// the spacers (above) and for the planar FR4 board, which MVB++ draws as one slab over the whole
// printed group -- copper and the interface insulation where a sheet sits included.
void cutByShunts(NamedShape& host, const std::vector<NamedShape>& all);

// The shunt collision gate, run on the finished assembly: a Shunt solid that shares volume with
// any other solid (turn, turn coating, core, core coating, insulation, FR4 board, pin, bobbin,
// spacer, another shunt) throws, naming both. The bobbin, the spacers and the FR4 board are checked
// too, AFTER their cut: a cut that silently did not happen is caught here. Terminal faces carry no volume.
void checkShuntCollisions(const std::vector<NamedShape>& all);

}  // namespace mvb
