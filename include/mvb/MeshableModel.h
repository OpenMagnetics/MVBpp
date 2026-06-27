#pragma once

#include "mvb/NamedShape.h"
#include "mvb/Symmetry.h"

#include <string>
#include <vector>

// A structured, ready-to-mesh description of a magnetic for FEM consumers (e.g. OMFEM).
// It carries the SAME geometry the builders produce, but with structured region
// identity (so consumers classify by data, not by re-parsing name strings) and — when
// a symmetry reduction is applied — which planes were used, the per-plane FE boundary
// condition, and the energy/flux multiplicity. This keeps ALL component geometry +
// symmetry knowledge in MVB++; the consumer only meshes, applies the labelled BCs, and
// scales by the multiplicity.
namespace OpenMagnetics { class Magnetic; }

namespace mvb {

// Logical role of a built region.
enum class Role { Core, Bobbin, Turn, Other };

// Structured identity decoded from the build-time naming convention (which this
// library owns and produces). winding/parallel/turn are set for Role::Turn;
// corePiece for Role::Core.
struct RegionTag {
    Role        role = Role::Other;
    std::string winding;        // set when role == Turn
    int         parallel = -1;  // turn parallel index (-1 if n/a)
    int         turn     = -1;  // turn index within the winding (-1 if n/a)
    int         corePiece = -1; // core piece index (0-based; -1 if n/a)
    std::string name;           // canonicalized source name (debugging / compat)
};

// FE wall condition for a symmetry cut face in a MAGNETICS model. A magnetic driven
// symmetrically about a geometric mirror plane has its flux crossing that plane
// normally, so the in-plane vector potential vanishes there -> FluxTangential
// (essential A = 0). This is the standard magnetics mirror BC and the default here;
// FluxNormal (natural / do-nothing) is provided for completeness.
enum class WallBC { FluxTangential, FluxNormal };

struct AppliedSymmetry {
    SymmetryPlane plane = SymmetryPlane::None;
    SymmetryHalf  half  = SymmetryHalf::Positive;
    WallBC        bc    = WallBC::FluxTangential;
};

struct TaggedRegion {
    NamedShape shape;
    RegionTag  tag;
};

struct MeshableModel {
    std::vector<TaggedRegion>    regions;       // reduced solids if symmetry applied
    std::vector<AppliedSymmetry> symmetry;      // empty => full model (no reduction)
    int                          multiplicity = 1;  // 2^(#applied planes)
};

struct MeshableOptions {
    bool includeBobbin       = false;
    int  symmetryPlanes      = 0;     // 0..3 (full / half / quarter / eighth)
    int  wirePolygonSegments = 16;
    int  corePolygonSegments = 16;
    bool paintCoating        = true;  // false => conducting (copper) footprint (FEM loss)
};

// Decode the build-time naming convention into a structured tag. Authoritative: the
// naming convention is produced by the builders in this library, so the decoder lives
// here as the single source of truth. Handles STEP-label prefixes ("Shapes/"),
// "/<dup>/COMPOUND" suffixes, and "_<int>" children added by sectioning/symmetry cuts.
RegionTag classify_region(const std::string& name);

// One call: build the named 3D solids (honouring opt), optionally apply
// opt.symmetryPlanes of the detected symmetry (reporting which planes/halves + the
// multiplicity + the per-plane BC), and tag every region. 2D sectioning (cut2DFaces)
// and the FEM truncation air box are the consumer's job and stay there.
MeshableModel buildMeshable(const OpenMagnetics::Magnetic& magnetic,
                            const MeshableOptions& opt);

} // namespace mvb
