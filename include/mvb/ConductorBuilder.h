#pragma once

#include "MAS.hpp"
#include "mvb/NamedShape.h"
#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
#include <array>
#include <cstdlib>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace OpenMagnetics { class Coil; }

namespace mvb {

// Real-winding geometry: ONE continuous copper body per (winding, parallel), replacing
// TurnBuilder's independent closed loops when DrawConfig::useRealWindingGeometry is on.
//
// Contract (the hard constraint):
//   - Every MKF turn keeps its exact MAS position: each turn is an on-station arc swept at
//     exactly (coordinates[0], coordinates[1]); MVB++ adds ONLY connecting geometry, and
//     only inside the seam sector / terminal routes MKF leaves unspecified.
//   - Electrical order comes from turnsDescription vector order filtered by
//     (turn.winding, turn.parallel) — MKF's documented contract.
//   - The path planner collision-checks every primitive (capsule model) against all
//     conductors and the winding-window bounds, and THROWS on any overlap beyond contact;
//     nothing is ever moved to "make it fit".
//
// Scope: concentric ROUND, RECTANGULAR and OBLONG columns and TOROIDS, with round/litz
// wire. Rectangular/planar/foil wire throws (the lead cross-section orientation through
// the exit bends is undefined until specified).
class ConductorBuilder {
public:
    struct Options {
        Options() {
            if (std::getenv("MVB_TOROID_MITRE_CORNERS")) toroidMitreCorners = true;
            if (const char* v = std::getenv("MVB_MIN_BEND_RADIUS")) {
                const double m = std::atof(v);
                if (m > 0) minBendRadius = m;
            }
        }
        int  wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS;
        // true  -> conductor swept at the OUTER (insulation) footprint (visualisation)
        // false -> swept at the CONDUCTING (copper) footprint (FEM)
        bool paintCoating = true;
        // ABT #685 (Alf, 2026-08-15): the TWO standard corner constructions for lead/connection
        // polylines, selected here. false (DEFAULT) -> BISECTION MITRE: plain segments meeting at
        // the waypoint, the conformal assembler slicing both sides on the angle-bisector plane —
        // the same joint as the wrap elbows, identical for every parallel regardless of stub
        // height. true -> ROUNDED: exact tangent fillet arcs where a leg has room (a corner whose
        // leg is too short for a valid fillet radius still falls back to the mitre — a fillet at
        // the wire radius is a horn torus OCC cannot build). The conditional fallback is exactly
        // what made 06_llc's three parallels show three different exit joints, which is why the
        // deterministic mitre is the default.
        bool roundedLeadCorners = false;
        // TOROID turn corners: false (default) = the round tangent fillets (today's certified
        // construction); true = sharp MITRE corners -- five straight runs per turn joined on
        // 45-degree bisector planes (Alf, 2026-08-27). Mitres trade the fillets' osculating
        // tangent junctions (OCC's worst contact class) for transversal plane joints, and at
        // segments > 0 every piece becomes an exact boolean-free polygon prism. Env:
        // MVB_TOROID_MITRE_CORNERS=1.
        bool toroidMitreCorners = false;
        // MINIMUM BEND RADIUS for every drawn corner/fillet, in METRES (Alf, 2026-08-27:
        // "a setting that is the minimum radius in any turn"). 0 (default) keeps each site's
        // own policy radius (kRoundCornerBendFactor x wire radius, 1.02 x wireRadius for rect
        // stations, 1.5 x wireRadius for lead fillets). A positive value FLOORS every one of
        // those choices: effectiveBend() below is the single chokepoint, so no corner anywhere
        // can be drawn tighter than this. It can only widen a bend, never tighten one -- a
        // radius below the wire's own is a horn torus OCC cannot build, and that guard stays.
        // Configure via ConductorBuilder::Options, the MVB_MIN_BEND_RADIUS env (metres), or the
        // step generator's --min-bend-radius flag.
        double minBendRadius = 0.0;
        double effectiveBend(double policyBend) const {
            return policyBend < minBendRadius ? minBendRadius : policyBend;
        }
        // false (OM / drawing) -> fast per-run compound: overlapping swept pipes, perfect for the
        //   3D viewer, seconds to build.
        // true  (FEM export)   -> the slow one-piece machinery: single continuous body per parallel
        //   where a sweep closes (columns, sparse toroids), a conformal mitre-jointed compound for
        //   dense toroids, fuse for single-turns. Meshable (no self-overlap) but far slower.
        bool femReady = false;
        // Azimuth (radians, azPointC convention: az -> (x = r cos az, z = -r sin az)) the
        // terminal leads should EXIT through -- the core's winding-window OPENING. NAN =
        // keep the default (az = pi/2, the -Z direction). Cores whose lateral legs cover
        // -Z (PQ, EP, ...) MUST redirect the leads or they terminate inside the core plate
        // (measured on 03_buck_pq3230: the lead tip landed exactly on the oblique plate
        // face and tetgen could not recover the contact edge).
        double leadExitAzimuth = std::numeric_limits<double>::quiet_NaN();
        // The CORE solids, for aiming the leads at the true window opening. Column metadata
        // under-describes real cores (a PQ's plates wrap most of the perimeter, leaving only
        // narrow slots); classifying the actual solids around the lead-tip radius finds the
        // genuine free arc. Empty -> no aiming beyond leadExitAzimuth.
        std::vector<TopoDS_Shape> coreObstacles;
        // DIAGNOSTIC ONLY. Build the conductors even though two of them overlap, so the
        // overlap can be LOOKED AT in CAD instead of only read about in an exception. A
        // collision means two pieces of copper occupy the same space: the result is not a
        // manufacturable part and must never be exported as one, meshed, or handed to a
        // solver. Nothing in the library sets this — only the tools and tests that exist to
        // investigate a specific refusal (see the [realwinding][diagnostic] tests).
        bool diagnosticSkipCollisionCheck = false;
        // ABT #871 — MULTI-COLUMN PLACEMENT. Section name -> the core column that section's
        // turns wrap, for the sections that do NOT wrap the main column. Resolved by the
        // caller (MagneticBuilder::WoundColumnResolver reads the MAS placement chain
        // section -> windingWindow -> column), because the column geometry lives on the CORE
        // and this builder is only handed the coil and the bobbin.
        //
        // A conductor whose turns wrap a lateral leg is built in that leg's own frame — turn
        // radials measured from the leg axis, the racetrack laid on the leg's half-dims —
        // and the finished path is translated onto the leg. Sections absent from this map
        // wrap the main column at the origin, which is every single-window design, so an
        // empty map is byte-identical to the pre-#871 builder.
        std::map<std::string, TurnBuilder::WoundColumnSpec> woundColumnPerSection;
        // ABT #1173 (WP4): when non-null, a toroidal build writes the y of its common terminal
        // plane here (NaN when the build has no toroidal conductor). A toroid on a base draws the
        // base with its top face on that plane (BaseBuilder.h).
        double* toroidTerminalPlaneOut = nullptr;
    };

    // REAL-PATH POLYLINES: the fully-assembled, collision-checked conductor centrelines
    // (wraps, links, Z end-runs, terminal leads -- everything emitConductor would sweep),
    // sampled per primitive, WITHOUT building any solid. No booleans anywhere: this is the
    // input for implicit (signed-distance / level-set) winding meshing, which sidesteps the
    // OCC weld/fragment failure classes entirely (ABT #490/#491).
    struct PathPolyline {
        std::string name;                 // "<winding> parallel <k>"
        double wireRadius = 0.0;          // round capsule radius (= half min dim for rect)
        bool isRectangular = false;
        double wireWidth = 0.0, wireHeight = 0.0;
        // Per-primitive sampled points (metres). Distance = min over prims of the
        // point-to-polyline capsule distance; per-prim grouping avoids phantom bridges
        // between non-contiguous runs.
        std::vector<std::vector<std::array<double, 3>>> prims;
        std::array<double, 3> end0{}, end1{};   // free ends (terminal port centres)
        std::array<double, 3> dir0{}, dir1{};   // OUTWARD end tangents (port normals)
    };
    // TERMINAL-LEAD COPPER LENGTH (ABT #1215). What MVB++ actually draws beyond the turns, per
    // winding, measured on the SAME finished centreline the conductor solids are swept from (the
    // whole buildAll pipeline runs -- assembly, lead routing, terminal fillets, the toroid
    // drop-to-plane, the collision gate -- and stops before any solid is built).
    //
    // WHAT COUNTS AS LEAD (one rule, flags only, no label matching):
    //   - every primitive with Primitive::isLead: the toroid lead chain (axial stub out of the
    //     hole, its round corner, the radial leg -- side-stepped at the rim when the router had
    //     to -- the second round corner, the -Y drop to the terminal plane), the concentric exit
    //     polylines and their fan legs, the concentric "lead corner" spirals and runs, foil
    //     lead wires;
    //   - every primitive with Primitive::terminalFillet: the fillet arcs TerminalFillet puts
    //     where a concentric wrap's terminal stub meets its lead (they replace a piece of the
    //     stub and a piece of the lead, and exist only because the lead is there).
    // NOT lead: turn copper (wraps, terminal stubs, bumps), isConnection primitives (inter-layer
    // links and their smoothing arcs, U/Z returns).
    //
    // WHICH END: a primitive's `terminal` tag (0 entrance, 1 exit) where the router set one
    // (toroids); otherwise its position -- a lead primitive before the conductor's first
    // non-lead primitive belongs to the entrance, after its last one to the exit. A lead piece
    // between two turn pieces, or a tag that contradicts the position, THROWS. A conductor that
    // is nothing but lead (a foil's soldered lead wire) takes its end from its construction.
    //
    // Lengths are exact (primLength); a lead piece without a closed-form length throws. A
    // conductor whose emission re-shapes its lead corners after this measurement (the
    // rectangular-wire per-primitive emitter inserts elbows between two straight lead legs)
    // throws instead of reporting a number the solid does not have.
    struct TerminalLeadEnd {
        size_t parallel = 0;
        std::string end;              // "entrance" | "exit"
        double length_m = 0.0;
        // One centreline primitive of this end, in path (current-flow) order, for audit: its
        // label, kind ("SEG" | "ARC3" | "SPIRAL"), exact length, endpoints (metres, MVB++ frame),
        // and for an ARC3 its bend radius and sweep (radians); both 0 for a SEG.
        struct Piece {
            std::string label;
            std::string kind;
            double length_m = 0.0;
            std::array<double, 3> start{}, end{};
            double radius_m = 0.0, sweep_rad = 0.0;
        };
        std::vector<Piece> pieces;
    };
    struct TerminalLeadLength {
        double total_m = 0.0;         // summed over both ends of every parallel
        size_t parallels = 0;
        std::vector<TerminalLeadEnd> per_end;   // parallel ascending, entrance before exit
    };
    // Keyed by winding name. `opts` must be the options the conductors are (or were) built with:
    // paintCoating changes the wire radius every bend is sized from, femReady changes which
    // corners are filleted.
    static std::map<std::string, TerminalLeadLength> measureTerminalLeadLengths(
        const OpenMagnetics::Coil& coil, const MAS::CoreBobbinProcessedDescription& bobbin,
        bool isToroidal, const Options& opts = {});

    static std::vector<PathPolyline> buildAllPaths(const OpenMagnetics::Coil& coil,
                                                   const MAS::CoreBobbinProcessedDescription& bobbin,
                                                   bool isToroidal,
                                                   const Options& opts = {});

    // Builds one solid per (winding, parallel), named "<winding> parallel <k>" (MAS style).
    // `coil` must be an MKF-wound coil (turnsDescription present, positions final);
    // `bobbin` is the resolved+patched processed description (same object MagneticBuilder
    // hands to TurnBuilder).
    static std::vector<NamedShape> buildAll(const MAS::Coil& coil,
                                            const MAS::CoreBobbinProcessedDescription& bobbin,
                                            bool isToroidal,
                                            const Options& opts = {});
    static std::vector<NamedShape> buildAll(const OpenMagnetics::Coil& coil,
                                            const MAS::CoreBobbinProcessedDescription& bobbin,
                                            bool isToroidal,
                                            const Options& opts = {});
};

} // namespace mvb
