#pragma once
// THE TERMINAL CORNER IS A FILLET, NOT A MITRE (ABT #961 / #839 follow-up, 2026-09-02).
//
// The terminal stub used to be the helix's OSCULATING CIRCLE so that its corner with the lead
// could be mitred analytically (a spiral pipe end cannot be cut reliably by OCC). Every circle
// leaves the helix by the torsion cubic (-kappa*tau/6 * s^3), ~100 nm over a mm-scale stub, and on
// boost_inductor_complete that ate a sibling lead's +60 nm margin (47 nm certified overlap).
// Now the stub IS the helix (an exact SPIRAL primitive), and the corner between it and the lead
// is what a bent wire actually is: a tangent BIARC (two circular arcs, radii >= the bend radius),
// consuming equal tangent lengths from the helix's end and the lead's start. Every joint is
// tangent, so nothing is mitred and nothing is bridged.
#include <mvb/WireAssembler.h>

#include <gp_Vec.hxx>

#include <functional>
#include <map>
#include <optional>
#include <limits>
#include <string>
#include <vector>

namespace mvb {

// WHAT "EQUAL" MEANS FOR A LEG MKF MADE EXACTLY ONE TANGENT LONG: a picometre, three orders below the 1 nm grid MKF
// ships its coordinates on, and far above the round-off of the sums a leg's room comes out of. The corner tests used
// to allow 8 ulps of the larger number (3.3e-19 m on a 0.1835 mm leg); a room computed as leg - what the neighbour
// took missed that by 13 and 26 ulps (3.5e-19 and 7.0e-19 m) on realwinding_interleaved_full_section_e16's Secondary,
// and four corners MKF planned round (need 0.183500 mm, room 0.183500 mm) were drawn sharp and reported as findings.
constexpr double kLegFitSlack = 1e-12;

// Fillets every terminal corner in `prims` -- an exact-helix "(terminal stub)" SPIRAL adjacent to a
// "lead seg" SEG, in either order -- in place. Returns how many corners were filleted. Throws,
// naming the corner, when a fillet with radius >= minBend does not fit.
// entranceRoll / exitRoll: the ROLL of the corner about the helix's own tangent, in radians from
// the principal normal towards the binormal. NaN (the default) keeps the historical corner, which
// snaps to whichever of the two natural planes the lead lies nearer.
//
// The corner alone cannot choose it. A fillet leaves the helix by a few nanometres to a few
// microns, and MKF places conductors at EXACTLY one coated diameter, so whether that excursion
// matters depends entirely on where the neighbour is -- which only the fan can see. What the fan
// computes is not a search: at exact touch a wire may not bend TOWARDS a neighbour at all, so the
// admissible rolls are a half-circle (those pointing away from it) and the answer is the one
// nearest the direction the lead actually needs. One construction, every design.
// FILLETS ARE PROVEN, NOT PICKED BLIND (06_llc, stepvalidity). The first biarc that fits the bend
// radius is not necessarily one the neighbours leave room for: on 06_llc Secondary p1's entrance
// the first fit (L = 0.588 mm, 2.79 deg of helix) reached under the layer-1 wrap passing one coated
// diameter overhead, 3.34 um inside its envelope. With `clears` given, every fitting candidate is
// proven against the other conductors before it is taken, and the search walks on to the next one
// when it does not clear. When NO candidate of a corner clears: with `uncleared` given the corner is
// left as it was and `*uncleared` counts it (the caller re-rolls and, failing that, throws with
// `*witness`); without it, it throws with the witness here.
// Fillets the corner between the SEGs prims[i] and prims[i + 1], which share a vertex, with an arc of
// `radius` tangent to both (flags -- lead, connection, terminal, turn -- from `flags`). False, nothing
// changed, when either leg is too short for its tangent length.
bool filletSegCorner(std::vector<Primitive>& prims, size_t i, double radius, const Primitive& flags,
                     const std::string& label);

// THE SLANTED STEP (Alf, 2026-10-28, "Slanted step"). A step of about one wire OD between two
// parallel straights holds two 90 deg bends only when it is at least 2.1 wire radii long, and MKF's
// steps (a Z dragback's radial drop, a lead's axial run to its turn's row) are one OD. The step
// [short SEG prims[s0]][its corner ARC3][long SEG] (or the mirror, longAfter = false) is redrawn as
// one straight from the step's junction end V to a point Q on the long leg, with a round kink at Q:
// V keeps its junction (a wrap face, a helix), whose corner the caller rounds; `away` points from V
// into that junction's piece, and the junction bend takes junctionBend * tan(theta / 2) /
// junctionShare of the slant. Q is the nearest point that holds both bends, so the copper leaves
// MKF's step by as little as possible. nullopt, nothing changed, when no Q on the long leg does.
struct SlantedStep {
    double slant = 0.0;       // how far along the long leg Q stands from the step's corner (m)
    double departure = 0.0;   // the step's corner's distance from the slant (m)
};
std::optional<SlantedStep> slantOneOdStep(std::vector<Primitive>& prims, size_t s0, bool longAfter,
                                          const gp_Vec& away, double junctionBend, double junctionShare,
                                          double kinkBend);

using FilletClears = std::function<bool(const std::vector<Primitive>& arcs, std::string* witness)>;
// ABT #2274 item E (agent 48, Alf: "MKF plans what is built"): the terminal corner MKF planned for a ramped
// round-column lead (ConnectionRoute::terminalClimbTangentLength / terminalClimbPhiDegrees). A planned corner is
// drawn exactly as MKF planned it -- the natural-plane roll, MKF's tangent length and bend angle -- whatever roll the
// caller asks for and without the clearance search: its clearance is VALIDATED by the caller (a collision finding),
// never cured here by another corner. MVB++'s construction landing anywhere else than MKF's plan throws.
struct PlannedClimb {
    double tangentLength = 0.0;   // m
    double phiDegrees = 0.0;
};
size_t filletTerminalCorners(std::vector<Primitive>& prims, double minBend, double wireRadius,
                             const std::string& who,
                             double entranceRoll = std::numeric_limits<double>::quiet_NaN(),
                             double exitRoll = std::numeric_limits<double>::quiet_NaN(),
                             const FilletClears& clears = {}, size_t* uncleared = nullptr,
                             std::string* witness = nullptr,
                             const std::optional<PlannedClimb>& entranceClimb = std::nullopt,
                             const std::optional<PlannedClimb>& exitClimb = std::nullopt);

// The same corner where a wrap meets a Z dragback's straight (its "(dragback) seg"), in either order: a
// round dragback has round ends (Alf, 2026-10-28), bent at `minBend` -- the radius its own corners take.
// `allowSlant` false: an end with no room on its one-OD step is sharp at once, never slanted.
size_t filletDragbackEnds(std::vector<Primitive>& prims, double minBend, double wireRadius, const std::string& who,
                          const FilletClears& clears = {}, size_t* uncleared = nullptr,
                          std::string* witness = nullptr, bool allowSlant = true);

// ABT #2274 (Alf, 2026-11-13: "use round corners everywhere!"; agent 45): a Z dragback's INTERIOR corners, bent at
// MKF's dragbackBendRadius `bend` once the dragback's ENDS have taken their share (call after filletDragbackEnds and
// filletLeadJunctions). See TerminalFillet.cpp. Returns the corners bent; a corner whose tangent does not fit is sharp
// and reported.
size_t filletDragbackInteriorCorners(std::vector<Primitive>& prims, double bend, double wireRadius, const std::string& who);
// A lead's corners that appendFilletedPolyline deferred (Primitive::deferredCorners: no measured claim for its terminal
// fillet), bent at `bend` by the same allotted-space rule once the terminal fillet has taken its share. Call after the
// terminal fillets. Returns the corners bent; a corner whose tangent does not fit is sharp and reported.
size_t filletDeferredLeadCorners(std::vector<Primitive>& prims, double bend, double wireRadius, const std::string& who);

// Replaces each straight radial LAYER LINK by a biarc tangent to the two wraps it joins, taking
// the azimuth it needs from each (ABT #969, Alf's option 2). `clears` decides between the
// candidate transitions: the roundest one is the kindest to the wire, but a long transition
// sweeps azimuth and can reach a neighbouring parallel (14_dab: 11.7 um), so the caller -- which
// is the only place that can see the other conductors -- vetoes. Candidates run from the roundest
// down to the tightest that still respects the bend radius; without a predicate the roundest wins.
// `maxStepOf` caps a link's transition (by the link's label): the roundest candidate taken is the one at that step of
// the azimuth sequence or below; `stepsOf`, when given, receives per link label the step taken and the first step that
// fits (the tightest), so a caller can lower a cap and try again (ConductorBuilder: the parallels of one winding take
// their sibling links at the SAME step, the roundest at which every one of them clears the others).
size_t smoothLayerLinks(std::vector<Primitive>& prims, double minBend, const std::string& who,
                        const std::function<bool(const Primitive&, const Primitive&)>& clears = {},
                        const std::function<int(const std::string&)>& maxStepOf = {},
                        std::map<std::string, std::pair<int, int>>* stepsOf = nullptr);

}  // namespace mvb
