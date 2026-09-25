#pragma once
// THE PORT-PLANE RULE (ABT #1400) -- one definition, shared by MVB++ (which places the terminal
// lead tips and gates them) and OMFEM (which snaps its air box onto them and refuses a model that
// breaks it). Before this header the rule lived only in OMFEM's MasMesher as the literals 0.5 and
// 0.9, and MVB++ placed its concentric lead tips by a rule of its own (outermost turn + four wire
// ODs of the window's round wire): on the two-switch forward the Primary entrance Litz lead's tip
// ended 0.6995 mm past the outermost Secondary foil copper against the 0.7545 mm OMFEM demands, and
// the mesh was refused ("the lead drop is too short for this cap").
//
// WHAT OMFEM DOES WITH A PORT (MasMesher, port snap + clearance guard):
//   1. It takes the terminal CAP face (the planar end of a free lead) and its axis-aligned bounding
//      box, and calls half that box's diagonal the cap's size, h = portCapHalfDiagonal(box).
//      For a round cap of radius r whose normal is along an axis the box is 2r x 2r x 0, so
//      h = sqrt(2) r; a faceted (inscribed polygon) cap is a little smaller.
//   2. It pulls the air-box face PAST the cap, INTO the lead, by kPortPlaneInsetCaps * h, so the
//      boundary slices the conductor on a clean planar cross-section (this is the inset for a cap
//      whose normal is axis-aligned, |n_axis| > 0.99 -- every MVB++ lead tip; OMFEM uses 1.5 h for
//      an oblique cap).
//   3. Every conductor solid that does NOT cross that face must stay at least
//      kPortPlaneGuardCaps * h inside it; closer, the sizing field between the copper and the box
//      face collapses and gmsh's HXT kernel aborts.
// So, measured along the lead's outward direction, the tip must lie at least
//      portPlaneClearance(h) = (kPortPlaneInsetCaps + kPortPlaneGuardCaps) * h
// beyond the outermost copper of every conductor that does not reach the port plane.
#include <Bnd_Box.hxx>

#include <cmath>
#include <stdexcept>

namespace mvb {

// OMFEM's inset of the port face into an axis-aligned cap, in cap half-diagonals.
constexpr double kPortPlaneInsetCaps = 0.5;
// OMFEM's guard between the port face and any conductor that does not cross it, in cap
// half-diagonals.
constexpr double kPortPlaneGuardCaps = 0.9;

// Half the diagonal of a terminal cap's axis-aligned bounding box (metres in, metres out) -- the
// cap size h both sides of the rule use. Takes the box corners exactly as Bnd_Box::Get returns
// them, so OMFEM can pass what it already measures.
inline double portCapHalfDiagonal(double xmin, double ymin, double zmin, double xmax, double ymax,
                                  double zmax) {
    const double dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
    if (!(dx >= 0.0 && dy >= 0.0 && dz >= 0.0))
        throw std::invalid_argument("portCapHalfDiagonal: the cap's bounding box is inverted");
    return 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline double portCapHalfDiagonal(const Bnd_Box& capBox) {
    if (capBox.IsVoid())
        throw std::invalid_argument("portCapHalfDiagonal: the cap's bounding box is void");
    double x0, y0, z0, x1, y1, z1;
    capBox.Get(x0, y0, z0, x1, y1, z1);
    return portCapHalfDiagonal(x0, y0, z0, x1, y1, z1);
}

// How far, along the lead's outward direction, a port's tip must lie beyond the outermost copper of
// every conductor that does not cross the port plane: the inset plus the guard.
inline double portPlaneClearance(double capHalfDiagonal) {
    if (!(capHalfDiagonal > 0.0))
        throw std::invalid_argument("portPlaneClearance: a cap half-diagonal must be positive");
    return (kPortPlaneInsetCaps + kPortPlaneGuardCaps) * capHalfDiagonal;
}

}  // namespace mvb
