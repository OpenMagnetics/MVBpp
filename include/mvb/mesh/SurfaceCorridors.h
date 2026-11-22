#pragma once
// Corridors of rectangular conductors, measured on their built surfaces. See the block comment in the .cpp.
#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <TopoDS_Shape.hxx>

namespace mvb::mesh {

// One corridor: where the fine size is asked for [m], the size [m] and the gap it was measured in [m].
struct CorridorCandidate { double p[3]; double h, gap; };
// A conductor solid as it will be meshed: the built conductor it is ("Secondary parallel 0"), the winding it belongs
// to, whether its section is rectangular (rectangular wire or foil) and the size its own surface is meshed at [m].
struct ConductorSolid {
    TopoDS_Shape shape; std::string conductor, winding; bool rectangular = false; double target = 0.0;
};
// The built centreline of a ROUND conductor: its name, its radius [m], its own target [m], and its sampled pieces [m].
// A foil's round leads are centrelines of their own ("Secondary parallel 0 entrance lead") and are PART of the foil's
// body: a centreline is never measured against the solid whose conductor name it carries or extends.
struct RoundCentreline {
    std::string name;
    double radius = 0.0, target = 0.0;
    std::vector<std::vector<std::array<double, 3>>> pieces;
};
struct SurfaceCorridors {
    std::vector<CorridorCandidate> points;
    size_t rectangularSolids = 0, curvedFaces = 0, samples = 0, hits = 0, ownHits = 0, inclined = 0, unbound = 0;
    size_t roundSamples = 0, roundHits = 0;
    double seconds = 0.0;
};

// The corridors that involve a rectangular conductor. Empty when no solid is rectangular. Throws std::runtime_error on
// a solid without a positive target, a round centreline without a positive radius or target, or solids with no faces.
SurfaceCorridors surface_corridor_points(const std::vector<ConductorSolid>& solids,
                                         const std::vector<RoundCentreline>& roundCentrelines);

}  // namespace mvb::mesh
