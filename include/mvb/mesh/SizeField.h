#pragma once
#include <cstddef>

#include <stdexcept>
#include <string>
#include <vector>

// Shared, composable mesh-SIZE policy for the OMFEM magnetics mesher (2D and 3D).
//
// The whole point: a magnetic needs DIFFERENT element sizes in different places --
//   * skin depth delta in the conductors (eddy/proximity), delta = 1/sqrt(pi f mu sigma)
//   * several elements across each air GAP (fringing / reluctance -> inductance)
//   * a size scaled to each COMPONENT's feature dimension (thin core legs, bobbin, air)
//   * curvature on round wires / posts
// gmsh resolves these by taking the MINIMUM of every active size source at each point. The
// historical OMFEM bug was calling setAsBackgroundMesh TWICE (conductor then coating) so only the
// LAST field survived. This module fixes that structurally: every feature contributes a sub-field,
// they are composed into ONE `Min` field, registered exactly once, and all competing size sources
// (MeshSizeFromPoints / ExtendFromBoundary) are disabled. Works identically in 2D and 3D.
//
// No silent fallbacks: skin_depth() and the builder throw on missing/invalid inputs rather than
// substitute a default (per project rule). The .cpp wraps the gmsh field API; this header is
// gmsh-free so typed consumers stay light.

namespace mvb::mesh {

// Tunable sizing policy. Every field has a documented default; none is silently substituted.
struct SizingPolicy {
    double N_skin   = 2.5;   // target elements per skin depth near a conductor surface
    double N_gap    = 6.0;   // target elements across each air gap (validated value)
    int    N_curv   = 16;    // elements per 2*pi on round wires/posts (2D only; 3D uses hausd)
    double dist_band = 4.0;  // refinement extends to dist_band * size_min away from the feature
    double grade    = 1.3;   // max neighbour element-size ratio (informational / MMG -hgrad)
};

// Skin depth [m]: delta = 1 / sqrt(pi * f * mu0 * mur * sigma). Copper ~ 66.1/sqrt(f[Hz]) mm.
// Throws std::invalid_argument on non-positive f / sigma (no fallback).
double skin_depth(double f_hz, double sigma, double mur = 1.0);

// Accumulates feature-driven sub-fields, then composes ONE Min background field.
// dim = 2 or 3. air_target is the coarse far-field ceiling [m].
class SizeFieldBuilder {
public:
    SizeFieldBuilder(int dim, double air_target, const SizingPolicy& policy = {});

    // Graded refinement toward CAD entities: size ramps from size_min at the entity to air_target
    // at dist_max away (StopAtDistMax so the far field defers to the ceiling). For dim==2 pass
    // CURVE tags; for dim==3 pass SURFACE tags. dist_min/dist_max default to size_min and
    // dist_band*size_min when <=0. Skipped (no field) if entities is empty.
    void add_distance_refinement(const std::vector<double>& entities, double size_min,
                                 double dist_min = 0.0, double dist_max = 0.0, int sampling = 100);

    // Graded refinement toward POINTS (0-dim CAD entities): a small halo around a singular corner,
    // size_min at the point ramping to air_target at dist_max. Use this, not a band along a whole
    // curve, when only the corner is singular (a gap mouth). Skipped (no field) if points is empty.
    void add_point_refinement(const std::vector<double>& points, double size_min,
                              double dist_min, double dist_max);

    // A size of its OWN for a region (2D: surface tags, 3D: volume tags): size_in inside and on its
    // boundary, the air ceiling elsewhere. For a material whose loss integrand needs a resolution of its
    // own (the ferrite's B^beta), not a halo around a feature. Skipped (no field) if regions is empty.
    void add_region_size(const std::vector<double>& regions, double size_in);

    // Axis-aligned box refined to size_in inside, air_target outside, with a graded `thickness`
    // transition band. Use for planar/E-core gaps and winding-window regions.
    void add_box(double size_in, double x0, double y0, double z0,
                 double x1, double y1, double z1, double thickness);

    // Cylinder (round centre-post gap): size_in inside a cylinder of given radius about an axis.
    void add_cylinder(double size_in, double radius, double thickness,
                      double cx, double cy, double cz, double ax, double ay, double az);

    // LOCAL CHORD-BOUND REFINEMENT (2026-09-20). A magnetic's mesh is bounded by the narrow air
    // CORRIDORS between nearly-touching conductor surfaces -- and only THERE. A winding revolved
    // about an axis turns each corridor of the 2D layout cross-section into a CIRCLE about that
    // axis, so the region that needs the fine size is a thin torus, not the whole winding window.
    // One RingBand is one such torus: `size_in` out to `tube` from the circle, then graded over
    // `thickness` to the air ceiling. All bands are composed into ONE MathEval field (gmsh
    // serialises MathEval evaluation, so one field with a nested min() is far cheaper than N).
    // Bands whose size_in is at or above the air ceiling are skipped -- they refine nothing.
    struct RingBand {
        double size_in;    // element size demanded inside the band [m]
        double radius;     // the circle's distance from the revolution axis [m]
        double axis_pos;   // the circle's position ALONG the revolution axis [m]
        double tube;       // fine out to this distance from the circle [m]
        double thickness;  // graded transition from `tube` out to the air ceiling [m]
    };
    // axis: 'X', 'Y' or 'Z' -- the axis the winding is revolved about ('Y' for MVB++ E/ETD/PQ).
    // Throws on a non-positive size/tube/thickness/radius (no silent skip of a real constraint).
    void add_ring_bands(char axis, const std::vector<RingBand>& bands);

    // Constant ceiling size everywhere (always include one so unconstrained air is sane).
    void add_ceiling(double size);

    bool empty() const { return fields_.empty(); }

    // Register a single Min over all sub-fields as the background mesh, disable competing size
    // sources, and set Mesh.MeshSizeMin/Max. min_floor<=0 => derived from the finest sub-field.
    // Returns the Min field tag. Throws if no sub-fields were added.
    int finalize(double min_floor = 0.0);

private:
    int dim_;
    double air_target_;
    SizingPolicy policy_;
    std::vector<int> fields_;
    double finest_ = 1e30;   // smallest size_in seen, for the derived MeshSizeMin floor
};

// WHAT THE WIRE NEEDS, WITHOUT PART NAMES (2026-09-03). A STEP that has been through OCC's
// fragment carries no names, so a mesher working from geometry alone cannot ask "which solid is
// copper". It does not have to: a conductor piece is the only thing in a magnetic whose smallest
// extent is a fraction of a millimetre while the part spans tens. Refining around those solids is
// what turns an unmeshable exact-surface assembly into a meshed one -- measured on
// 01_simple_inductor_etd34_n87 at segments=0: a uniform 1.6 mm ceiling meshed 5 of 30 volumes and
// left 571 mm3 of copper out ("invalid boundary mesh, overlapping facets"), while refining the
// thin solids to 3 elements across meshed 38 of 38 with no unmeshed volume and no gmsh error.
// Reads the CURRENT gmsh model, so call it after the geometry is built and fragmented.
struct ThinSolidRefinement {
    std::size_t solids = 0;    // solids judged thin (conductor pieces)
    std::size_t faces = 0;     // faces handed to the distance field
    double thinnest = 0.0;     // smallest bounding-box extent found [m]
    double size_min = 0.0;     // element size applied at those faces [m]
};
// thin_cut [m]: a solid counts as thin when its smallest bounding-box extent is at or under this.
// A solid must ALSO be elongated -- longest extent at least `min_aspect` times the shortest -- or
// the test degenerates on parts that have no thin features at all: there thin_cut, derived from
// the smallest solid present, admits the core itself (measured on 03_buck, whose four chunky
// solids were all "thin" and got the core refined for nothing). Wire is long and thin; a core
// half or a bobbin is not.
// elements_across: how many elements to put across the thinnest solid (>= 1).
ThinSolidRefinement add_thin_solid_refinement(SizeFieldBuilder& builder, double thin_cut,
                                              int elements_across, double min_aspect = 3.0);

// 3D meshing-algorithm fallback chain (Mesh.Algorithm3D ids), most-robust first.
// 1 = Delaunay (robust boundary recovery, honours size fields + embedded entities),
// 10 = HXT (fast/parallel for large clean models). NEVER 7 (MMG3D) for first-pass generation --
// MMG3D ignores general size fields and is an adaptive backend, not a robust generator.
std::vector<int> algorithm_chain_3d();

// 2D chain (Mesh.Algorithm ids): 6 = Frontal-Delaunay (best quality), 5 = Delaunay (fallback).
// 7 = BAMG only when an anisotropic / boundary-layer metric is requested.
std::vector<int> algorithm_chain_2d(bool boundary_layer = false);

// Generate a `dim`-D mesh trying each algorithm in `chain` until one yields a non-empty mesh.
// Sets Mesh.Algorithm (2D) / Mesh.Algorithm3D (3D) per attempt, catches gmsh failures, escalates
// Mesh.ToleranceInitialDelaunay on retry. Returns the algorithm id that succeeded. Throws
// std::runtime_error with the collected failures if every algorithm fails (no silent bad mesh).
int generate_robust(int dim, const std::vector<int>& chain);

// VALIDATION GATE: throw std::runtime_error if any physical group of dimension `topdim` (the
// material regions: faces in 2D, volumes in 3D) has zero mesh elements. A region tagged in the CAD
// but swallowed by the boolean fragment / mis-classified would otherwise produce a SILENTLY wrong
// solve (e.g. a turn with no copper). Fails loud instead. Boundary groups (lower dim) are ignored.
void validate_regions_nonempty(int topdim);

// Run gmsh mesh optimisation (Mesh.Optimize, + Netgen for 3D) to repair slivers before writing.
// Best-effort: logs but does not throw (the region gate is the hard check).
void optimize_mesh(int topdim);

}  // namespace mvb::mesh
