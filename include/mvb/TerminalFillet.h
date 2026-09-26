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

#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace mvb {

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
using FilletClears = std::function<bool(const std::vector<Primitive>& arcs, std::string* witness)>;
size_t filletTerminalCorners(std::vector<Primitive>& prims, double minBend, double wireRadius,
                             const std::string& who,
                             double entranceRoll = std::numeric_limits<double>::quiet_NaN(),
                             double exitRoll = std::numeric_limits<double>::quiet_NaN(),
                             const FilletClears& clears = {}, size_t* uncleared = nullptr,
                             std::string* witness = nullptr);

// Replaces each straight radial LAYER LINK by a biarc tangent to the two wraps it joins, taking
// the azimuth it needs from each (ABT #969, Alf's option 2). `clears` decides between the
// candidate transitions: the roundest one is the kindest to the wire, but a long transition
// sweeps azimuth and can reach a neighbouring parallel (14_dab: 11.7 um), so the caller -- which
// is the only place that can see the other conductors -- vetoes. Candidates run from the roundest
// down to the tightest that still respects the bend radius; without a predicate the roundest wins.
size_t smoothLayerLinks(std::vector<Primitive>& prims, double minBend, const std::string& who,
                        const std::function<bool(const Primitive&, const Primitive&)>& clears = {});
// Bend radius = kRoundCornerBendFactor * wire radius (the corner rule everywhere else).
inline size_t filletTerminalCorners(std::vector<Primitive>& prims, double minBend,
                                    const std::string& who,
                                    double entranceRoll = std::numeric_limits<double>::quiet_NaN(),
                                    double exitRoll = std::numeric_limits<double>::quiet_NaN()) {
    return filletTerminalCorners(prims, minBend, minBend / kRoundCornerBendFactor, who,
                                 entranceRoll, exitRoll);
}

}  // namespace mvb
