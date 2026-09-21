#pragma once
// ABT #685 (Alf, 2026-08-16) — THE CENTRELINE VOCABULARY AND THE ONE ASSEMBLER.
//
// Alf's standing rule for real-winding Coil geometry: the six approved chunk primitives
// (Turn, Turn-with-bump, U connection, Z connection, inter-section connection, Terminal)
// are the ONLY producers of geometry, and they reach OCCT through exactly ONE assembler —
// this one. No other translation unit may build a solid. A seventh primitive, or a second
// way of turning a centreline into copper, needs Alf's explicit approval first.
//
// The split is: a chunk emits a CENTRELINE (a chain of exact analytic pieces — straight,
// arc, spiral), and assembleWire turns that chain into copper, joining consecutive pieces
// with one of the two approved corner constructions. Everything OCCT knows about wire lives
// below this line; everything the winder knows about lives above it.
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_XYZ.hxx>
#include <gp_XY.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

namespace mvb {

constexpr double kPi = std::numbers::pi;
constexpr double kTwoPi = 2.0 * std::numbers::pi;
// Curved primitives are collision-checked as sampled polylines; the sampling step is
// chosen so each chord sags inward by at most this fraction of the wire radius. This is
// a DISCRETIZATION-DENSITY parameter of the measurement (like polygonSegments), not a
// clearance policy.
constexpr double kMaxSagFraction = 0.02;
// ABT #685: a ROUND corner's bend radius, as a multiple of the wire radius. It must be
// strictly greater than 1: a centreline radius EQUAL to the wire radius sweeps a HORN TORUS,
// whose tube touches its own axis of revolution, so no valid solid exists. It is also the
// physics — a wire does not bend tighter than its own radius. Where a corner has no room for
// this, it takes the BisectionMitre construction instead, which has no radius to violate.
//
// 1.05 is the smallest margin that is buildable, and it is deliberately the SAME floor
// appendFilletedPolyline already applies to its fillets — one number for one constraint. It is
// not slack to be spent: a round corner intrudes (factor - 1) * wireRadius into its own inside,
// and on a toroid that space belongs to the terminal lead's drop. Measured on the common-mode
// choke, a 1.5 factor put the poloidal corner 0.5 * wireRadius = 219 um into the entrance
// lead's corridor and the router (correctly) refused the design.
constexpr double kRoundCornerBendFactor = 1.05;
// A DRAGBACK corner's fillet radius, as a multiple of the dragback's own bend radius (not of the
// wire radius). Strictly above 1 for the same horn-torus reason as above; 1.5 rather than the
// buildable floor because real wire cannot bend to its own radius without crushing the inner
// fibre, and a dragback is a formed bend in the finished part. Named here, not repeated as a bare
// literal at the site (Alf, 2026-09-19: no magic numbers); the value is unchanged from the site's.
constexpr double kDragbackFilletBendFactor = 1.5;
// FACETED MODE ONLY (segments > 0): an arc bent tighter than this many wire radii is revolved
// on the EXACT round profile, not the polygon (ABT #1111, 2026-09-06; Alf's call: "do option
// 1"). A 12-gon revolved about an axis 1.05 r away puts its innermost facet ~0.08 r from the
// axis; on a short remainder arc (the terminal fillet's 5-degree third piece) that strip is
// under 2 um long and BOPAlgo reports too-small edges and self-intersections -- real
// micro-geometry, faulty at any scale (single_switch Reset Winding solid 4), unlike the
// parametrisation defect the exporter now fixes. The torus of the same arc has no strip. 1.5
// covers every corner fillet and toroid window corner (1.05-1.15 r); the bump-riser biarcs sit
// at 9-140 r and stay faceted.
constexpr double kFacetTightArcExactRatio = 1.5;
// The measurement's own error bound: a sampled polyline lies within this sag of the true
// curve, so a measured polyline-polyline distance can under-read the true centreline
// distance by at most the two prims' sag bounds. Derived from the sampling rule actually
// applied in samplePrim/curveSampleCount -- never a policy constant.
inline double samplingSag(double wireRadius) {
    return std::max(1e-6, kMaxSagFraction * wireRadius);
}

// Azimuth convention (OCCT right-handed rotation about +Y, the column axis), generalised
// to a vertical axis through the horizontal point (cx, cz):
//   pos(r, y, az) = (cx + r cos az, y, cz - r sin az)
inline gp_Pnt azPointC(double cx, double cz, double r, double y, double az) {
    return gp_Pnt(cx + r * std::cos(az), y, cz - r * std::sin(az));
}

// Rodrigues rotation of v about the unit axis by ang.
inline gp_XYZ rotateXYZ(const gp_XYZ& v, const gp_XYZ& axisUnit, double ang) {
    return v * std::cos(ang) + axisUnit.Crossed(v) * std::sin(ang) +
           axisUnit * (axisUnit.Dot(v) * (1.0 - std::cos(ang)));
}

// ---------------------------------------------------------------------------------------
// Path model.
//   SEG    — straight segment (leads, racetrack faces, toroidal tubes and chords).
//   ARC3   — exact circular arc about an arbitrary axis: P(t) = c + rot(axis, t)·v0 for
//            t in [0, sweep] (racetrack corners, oblong caps, toroidal elbows).
//   SPIRAL — spiral about a VERTICAL axis through (cx, cz): azimuth linear in the
//            parameter, radius/height linear (round-column wraps) or cosine-blended
//            (blend=true, oblong cap ramps) so the end tangents are purely azimuthal.
//   BLEND  — S-curve between two points whose tangent is parallel to a common direction
//            at BOTH ends (the rectangular-column transition ramp): the lateral offset
//            follows a cosine blend, so the chain stays tangent-continuous through the
//            crossings while passing EXACTLY through them. (A straight ramp leaves
//            ~15-degree kinks at its junctions; OCCT's pipe-shell corner transitions
//            flare unboundedly on such kinks — observed: 9 mm spikes.)
struct Seg {
    gp_Pnt a, b;
};
struct Arc3 {
    gp_Pnt c;            // arc centre
    gp_XYZ axis{0, 1, 0};  // unit rotation axis (right-handed sweep)
    gp_XYZ v0{0, 0, 0};    // centre -> start vector, |v0| = bend radius
    double sweep = 0;
};
struct Spiral {
    double cx = 0, cz = 0;   // vertical axis position
    double r0 = 0, y0 = 0, az0 = 0;
    double r1 = 0, y1 = 0, az1 = 0;
    bool blend = false;      // cosine-blend r/y (end tangents purely azimuthal)
};
struct Blend {
    gp_Pnt a, b;
    gp_XYZ u{1, 0, 0};       // unit tangent direction at both ends
};
struct Primitive {
    enum Kind { SEG, ARC3, SPIRAL, BLEND } kind = SEG;
    Seg seg{};
    Arc3 arc{};
    Spiral spiral{};
    Blend blendc{};
    std::string label;
    // Electrical turn ordinal this primitive belongs to (entrance lead = first turn's,
    // exit lead = last turn's). A continuous wire legitimately contacts itself only
    // between CONSECUTIVE turns; the gate exempts same-conductor pairs with
    // |ordinal diff| <= 1 and checks everything farther apart.
    size_t turnOrdinal = 0;
    // Lead primitives are the replayed MKF connection routes. Emission splits runs at
    // lead<->wrap boundaries: those junctions are the only 90-degree corners of the path
    // (wrap chains are tangent-continuous), and OCCT's pipe-shell mitre flares unboundedly
    // across them, while lead chains sweep robustly as exact cylinders + sphere elbows.
    bool isLead = false;
    // An inter-layer link (MKF's blue ConnectionReservedSpace box): the straight radial
    // step that carries the wire from one concentric layer out to the next. Like a lead it
    // meets its neighbouring wraps at 90-degree corners, so it is isolated into its own run
    // and swept piecewise (cylinder + sphere elbows) rather than through the pipe-shell.
    bool isConnection = false;
    // Which TERMINAL a lead primitive belongs to: 0 = entrance, 1 = exit, -1 = not a terminal
    // lead (turn copper, connections, concentric leads that do not tag themselves). The
    // collision gate exempts same-conductor pairs of adjacent turn ordinals (they touch by
    // construction); the two terminals of a one- or two-turn conductor share those ordinals
    // yet are NOT each other's continuation, so a tagged pair of different terminals is gated.
    int terminal = -1;
    // ABT #1215 — MEASUREMENT TAG ONLY, read by nothing that builds geometry. True on the arcs
    // filletTerminalCorners inserts where a concentric wrap's terminal stub meets its lead. They
    // are copied from the stub (so they carry isLead = false and must keep it: isLead changes how
    // a run is split and emitted), yet they exist only because the lead leaves the wrap there, so
    // the terminal-lead length measurement counts them as lead copper.
    bool terminalFillet = false;
};

// THE TWO APPROVED CORNER CONSTRUCTIONS (Alf, ABT #685: "any corner can be made in two
// ways: round corner, or bisection mitre"). Nothing else may join two centreline pieces.
enum class CornerStyle {
    // Both sides are grown past the joint and sliced on the ANGLE-BISECTOR plane, so they
    // share one identical face — the default, and the only construction that is the same
    // for every parallel regardless of how long its legs are.
    BisectionMitre,
    // An exact tangent fillet arc replaces the joint. Needs a leg longer than the wire
    // radius (a fillet AT the wire radius sweeps a horn torus, which is not a valid solid),
    // so it is only available where the geometry has room for it.
    Round,
};

// ---- the centreline -> copper assembler --------------------------------------------------
// Builds one solid per centreline piece, abutting exactly: tangent junctions meet on flush
// perpendicular discs, corner junctions are grown and sliced per `corners`. No boolean fuse,
// no overlap: coincidence is exact by construction, which is the input class OCCT's boolean
// spec and gmsh's fragment actually support.
// `primIndexPerSolid`, when given, receives one entry per SOLID of the returned compound: the
// index into `centreline` of the piece that solid was built from. That is the only correct way to
// name the parts -- a solid count can exceed the primitive count (a mitre trim may fragment a
// piece), so neither an index nor a centroid match is reliable.
// `skipWeld` (ABT #1265): emit the pieces UNWELDED, as the weld-lens fallback already does.
// For a cutting tool -- the coated-envelope build that carves the bobbin's lead slots -- the
// weld is pure waste: overlapping pieces cut exactly the same volume, and that build was
// costing as much as the product one (measured on 01_etd34 at --segments 4: two identical
// conformal assemblies, 126 s total).
// The ruled prism between two cap polygons (vertex i of `start` joined to vertex i of `end`),
// the solid every faceted straight is. Lofted in the MILLIMETRE frame and returned in metres at
// the model's own confusion -- see the definition for why a metre-frame loft is not (ABT #1266).
// Returns a null shape, with the reason in *why, when OCC cannot build a valid solid.
TopoDS_Shape loftRuledPrism(const std::vector<gp_Pnt>& start, const std::vector<gp_Pnt>& end,
                            std::string* why = nullptr);

TopoDS_Shape assembleWire(const std::vector<const Primitive*>& centreline, double wireRadius,
                          int polygonSegments, CornerStyle corners = CornerStyle::BisectionMitre,
                          std::vector<size_t>* primIndexPerSolid = nullptr,
                          bool skipWeld = false);

// Sampling / geometry queries on a centreline piece, shared by the chunk builders and the
// collision gate. Pure measurement — they build no copper.
int curveSampleCount(double radius, double azSpan, double wireRadius);
// A SPIRAL is not an arc: its radius (and height) run while its azimuth turns, so the
// azimuthal rule above under-samples it. See the definition for the derivation.
int spiralSampleCount(const Spiral& sp, double wireRadius);
std::vector<gp_Pnt> samplePrim(const Primitive& p, double wireRadius);
std::pair<gp_Pnt, gp_Pnt> primEndpoints(const Primitive& p);
// EXACT centreline arc length of one piece, metres (ABT #1215). SEG: chord; ARC3: the radius
// perpendicular to the axis times |sweep|; linear SPIRAL: the closed form of
// integral sqrt(dr^2 + dy^2 + (r(t) daz)^2) dt. A BLENDED spiral and a BLEND have no closed form
// (their length is an elliptic integral) and THROW — no quadrature estimate is returned.
double primLength(const Primitive& p);
gp_Dir primFwdStart(const Primitive& p, double r);
gp_Dir primFwdEnd(const Primitive& p, double r);

// Profile wires/faces and the single-piece edge, used by the assembler and (until they are
// retired) by the legacy sweep paths in ConductorBuilder.
TopoDS_Wire wireProfileWire(const gp_Pnt& center, const gp_Dir& normal, double radius, int segments);
TopoDS_Face wireProfile(const gp_Pnt& center, const gp_Dir& normal, double radius, int segments);
TopoDS_Wire wireProfileWireSplit(const gp_Pnt& center, const gp_Dir& normal, double radius, int segments);
TopoDS_Wire rectProfileWire(const gp_Pnt& center, const gp_Dir& tangent, const gp_Dir& axialAxis,
                            double width, double height);
// Growth (overA/overB, in arc length) is folded into the pcurve for spirals that support it --
// see primEdgeGrowsAnalytically. Everything else ignores it and is grown by the caller.
TopoDS_Edge primEdge(const Primitive& pr, double wireRadius, double overA = 0.0,
                     double overB = 0.0);
bool primEdgeGrowsAnalytically(const Primitive& pr);

// Post-checks on assembled copper.
bool hasDegenerateSheetFace(const TopoDS_Shape& shape, double wireRadius);
TopoDS_Shape pruneDegenerateSolids(const TopoDS_Shape& shape, double wireRadius);

} // namespace mvb
