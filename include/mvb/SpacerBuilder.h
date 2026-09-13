#pragma once

#include "MAS.hpp"
#include "mvb/NamedShape.h"
#include <TopoDS_Shape.hxx>
#include <string>
#include <vector>

namespace mvb {

// One plastic shim of an additive-gapped core, as MAS describes it.
//
// MAS `core/spacer.json`: `dimensions` is "the cube defining the spacer" and MKF fills it
// (Core.cpp, ABT #1170) as [X, Y, Z] in the SAME frame the half-sets are placed in — Y is
// the column axis, so dimensions[1] IS the shim thickness (the additive gap). `coordinates`
// is the CENTRE of that cube on all three axes. The old header comment here claimed
// "(dimensions[0], dimensions[2], dimensions[1])" and the .cpp comment claimed
// "[width(X), thickness(Z), depth(Y)]"; both contradicted the code, which was right.
struct Spacer {
    TopoDS_Shape shape;
    // `insulationMaterial` off the MAS element — a material NAME, or the name of the inline
    // insulation-material record. Carried through so OMFEM can give the region a dielectric
    // constant and a thermal conductivity instead of guessing.
    std::string  insulationMaterialName;
};

class SpacerBuilder {
public:
    // Build one box per SPACER element of the geometrical description, in MVB++'s
    // column-along-Y frame, centred on all three coordinates.
    //
    // Throws on a SPACER element that cannot be drawn (missing/short `dimensions`, a
    // non-positive extent, fewer than three coordinates). A shim that MAS declares and MVB++
    // silently drops is a part that leaves the factory with a gap nobody modelled — there is
    // no safe default here.
    //
    // Returns an empty vector when the core has no spacers (no additive gap): that is a
    // legitimate core, not a failure.
    static std::vector<Spacer> buildSpacers(
        const std::vector<MAS::CoreGeometricalDescriptionElement>& geometricalDescription);

    // Convenience: fuse all spacer shapes into a single compound (or the
    // single spacer itself if only one exists). Returns a null shape when
    // the list is empty.
    static TopoDS_Shape buildSpacersCompound(
        const std::vector<MAS::CoreGeometricalDescriptionElement>& geometricalDescription);
};

// Append the core's spacers to an assembly as `Spacer_<i>`, role `Spacer`.
//
// WP0 (ABT #1169) adds `MagneticBuilder::appendAccessorySolids(all, magnetic, opts)` to both
// `buildAllNamed` overloads; that hook calls this. Kept as a free function so WP1 owns no
// part of the hook itself and the two work packages do not collide in MagneticBuilder.cpp.
void appendSpacerSolids(std::vector<NamedShape>& all, const MAS::MagneticCore& core);

} // namespace mvb
