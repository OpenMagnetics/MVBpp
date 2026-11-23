#include <mvb/TerminalFillet.h>

#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <gp_Lin.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <cmath>
#include <functional>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <optional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace mvb {
namespace {


struct SpiralEval {
    gp_Pnt p;
    gp_Vec t;        // unit travel tangent
    double dsdaz;    // arc length per radian
};

// azPointC convention (WireAssembler): P(az) = (cx + r cos az, y, cz - r sin az); r and y linear
// in az for a non-blend spiral.
SpiralEval evalSpiral(const Spiral& sp, double az) {
    const double span = sp.az1 - sp.az0;
    const double f = std::abs(span) > 0 ? (az - sp.az0) / span : 0.0;
    const double r = sp.r0 + (sp.r1 - sp.r0) * f;
    const double y = sp.y0 + (sp.y1 - sp.y0) * f;
    const double rp = std::abs(span) > 0 ? (sp.r1 - sp.r0) / span : 0.0;
    const double yp = std::abs(span) > 0 ? (sp.y1 - sp.y0) / span : 0.0;
    SpiralEval e;
    e.p = gp_Pnt(sp.cx + r * std::cos(az), y, sp.cz - r * std::sin(az));
    gp_Vec d(rp * std::cos(az) - r * std::sin(az), yp, -rp * std::sin(az) - r * std::cos(az));
    e.dsdaz = d.Magnitude();
    if (span < 0) d.Reverse();
    if (d.Magnitude() > 0) d.Normalize();
    e.t = d;
    return e;
}

struct ArcOut {
    bool straight = false;
    Primitive pr;
    double radius = 0.0;
    gp_Vec tEnd;
};

// The circular arc that leaves P along the unit tangent t and reaches Q.
ArcOut arcFromTangent(const gp_Pnt& P, const gp_Vec& t, const gp_Pnt& Q, const Primitive& like,
                      const char* what) {
    ArcOut out;
    const gp_Vec q(P, Q);
    gp_Vec n = q - t * q.Dot(t);
    if (n.Magnitude() < 1e-15 * std::max(q.Magnitude(), 1e-12)) {
        out.straight = true;
        out.pr = like;
        out.pr.kind = Primitive::SEG;
        out.pr.seg = Seg{P, Q};
        out.pr.label = like.label + std::string(" ") + what;
        out.tEnd = t;
        out.radius = std::numeric_limits<double>::infinity();
        return out;
    }
    const double R = q.SquareMagnitude() / (2.0 * n.Magnitude());
    n.Normalize();
    const gp_Pnt C = P.Translated(n * R);
    gp_Vec ax = t.Crossed(n);
    ax.Normalize();
    const gp_Vec vP(C, P), vQ(C, Q);
    const double sweep = vP.Angle(vQ);
    out.pr = like;
    out.pr.kind = Primitive::ARC3;
    out.pr.arc.c = C;
    out.pr.arc.axis = ax.XYZ();
    out.pr.arc.v0 = vP.XYZ();
    out.pr.arc.sweep = sweep;
    out.pr.label = like.label + std::string(" ") + what;
    out.radius = R;
    gp_Vec tE = ax.Crossed(vQ);
    tE.Normalize();
    out.tEnd = tE;
    // Self-check: the arc must land on Q.
    const gp_Pnt landed = P.Rotated(gp_Ax1(C, gp_Dir(ax)), sweep);
    if (landed.Distance(Q) > 1e-9 * std::max(R, 1e-9)) {
        std::ostringstream m;
        m << "TerminalFillet: arc '" << out.pr.label << "' misses its landing by "
          << landed.Distance(Q) * 1e6 << " um";
        throw std::runtime_error(m.str());
    }
    return out;
}

struct Biarc {
    ArcOut a1, a2;
    double minRadius = 0.0;
};

// A fillet never loops: an arc sweeping more than 120 deg is the polygon construction picking
// the far solution (06_llc's two 178 deg arcs, an S-loop of one wire diameter).
constexpr double kMaxFilletSweep = 120.0 * kPi / 180.0;
bool loops(const ArcOut& a) { return !a.straight && a.pr.arc.sweep > kMaxFilletSweep; }

// Unit principal normal of the spiral at `az` (the direction its curvature points, i.e. where a
// bend "in the plane of the helix" goes): the tangent-orthogonal part of the second derivative.
gp_Vec spiralPrincipalNormal(const Spiral& sp, double az) {
    const double span = sp.az1 - sp.az0;
    const double h = 1e-4 * std::max(std::abs(span), 1e-6);
    const SpiralEval m = evalSpiral(sp, az), a = evalSpiral(sp, az - h), b = evalSpiral(sp, az + h);
    gp_Vec dd(a.p.X() + b.p.X() - 2.0 * m.p.X(), a.p.Y() + b.p.Y() - 2.0 * m.p.Y(),
              a.p.Z() + b.p.Z() - 2.0 * m.p.Z());
    dd -= m.t * dd.Dot(m.t);
    if (dd.Magnitude() < 1e-30) return gp_Vec(0, 0, 0);
    dd.Normalize();
    return dd;
}

// The biarc with UNEQUAL tangent lengths a1 (at P1) and a2 (at P2): control polygon P1, A, B,
// P2 with |B - A| = a1 + a2 (each arc tangent to two consecutive legs). a2 follows from a1.
std::optional<Biarc> biarcSplit(const gp_Pnt& P1, const gp_Vec& t1, const gp_Pnt& P2,
                                const gp_Vec& t2, double a1, const Primitive& like) {
    const gp_Pnt A = P1.Translated(t1 * a1);
    const gp_Vec d(A, P2);
    const double den = 2.0 * (d.Dot(t2) + a1);
    if (std::abs(den) < 1e-300) return std::nullopt;
    const double a2 = (d.SquareMagnitude() - a1 * a1) / den;
    if (!(a2 > 0.0) || !std::isfinite(a2)) return std::nullopt;
    const gp_Pnt B = P2.Translated(t2 * (-a2));
    gp_Vec AB(A, B);
    if (AB.Magnitude() < 1e-15) return std::nullopt;
    const gp_Pnt J = A.Translated(AB * (a1 / (a1 + a2)));
    AB.Normalize();
    Biarc out;
    out.a1 = arcFromTangent(P1, t1, J, like, "fillet arc 1");
    out.a2 = arcFromTangent(J, AB, P2, like, "fillet arc 2");
    if (out.a2.tEnd.Angle(t2) > 1e-9) return std::nullopt;
    out.minRadius = std::min(out.a1.radius, out.a2.radius);
    return out;
}

struct Fillet {
    std::vector<ArcOut> arcs;   // in travel order
    double minRadius = 0.0;
};

std::optional<Biarc> biarc(const gp_Pnt& P1, const gp_Vec& t1, const gp_Pnt& P2, const gp_Vec& t2,
                           const Primitive& like);

// THREE ARCS: the helix-side arc is a planar bend of radius `R` in the helix's osculating plane
// at the tangency (tangent + principal normal), turning by `phi` towards the lead's side; from
// its far end a biarc reaches the lead. A two-arc fillet cannot do this: its junction lies on
// the polygon leg between the lead's tangent point and the helix's, and the lead's point is off
// the osculating plane by the helix's sagitta (13 um on boost), so the helix-side arc always
// tilts (2.1 mrad measured) and grazes the axial sibling one pitch away by ~1 nm. A radial
// departure is first-order neutral and second-order clearing towards both siblings; the
// junction sits R(1-cos phi) off the helix, microns, before the biarc turns towards the lead.
std::optional<Fillet> triarcHelixPlane(const gp_Pnt& P1, const gp_Vec& t1, const gp_Pnt& P2,
                                       const gp_Vec& t2, bool helixAtP2, const gp_Vec& nHelix,
                                       double R, double phi, const Primitive& like) {
    const gp_Pnt PH = helixAtP2 ? P2 : P1;
    const gp_Vec tH = helixAtP2 ? t2 : t1;        // travel tangent at the helix point
    const gp_Pnt PL = helixAtP2 ? P1 : P2;
    if (nHelix.Magnitude() < 1e-30) return std::nullopt;
    gp_Vec binormal = tH.Crossed(nHelix);
    if (binormal.Magnitude() < 1e-30) return std::nullopt;
    binormal.Normalize();
    // bend towards the side of the osculating plane where the lead lies
    const double side = gp_Vec(PH, PL).Dot(nHelix) >= 0 ? 1.0 : -1.0;
    const gp_Vec n = nHelix * side;
    const gp_Pnt C = PH.Translated(n * R);
    Fillet out;
    if (helixAtP2) {
        // entrance: ... -> biarc -> [J -> P2 along the planar arc]. Walk BACK from P2 by phi.
        gp_Vec ax = tH.Crossed(n);   // rotation axis for travel J -> P2 (right-handed sweep)
        ax.Normalize();
        const gp_Pnt J = P2.Rotated(gp_Ax1(C, gp_Dir(ax)), -phi);
        gp_Vec tJ = tH.Rotated(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(ax)), -phi);
        tJ.Normalize();
        const auto b = biarc(P1, t1, J, tJ, like);
        if (!b) return std::nullopt;
        ArcOut a3 = arcFromTangent(J, tJ, P2, like, "fillet arc 3");
        if (a3.tEnd.Angle(t2) > 1e-9) return std::nullopt;
        out.arcs = {b->a1, b->a2, a3};
        out.arcs[0].pr.label = like.label + " fillet arc 1";
        out.arcs[1].pr.label = like.label + " fillet arc 2";
        out.minRadius = std::min({b->a1.radius, b->a2.radius, a3.radius});
    }
    else {
        // exit: [P1 -> J along the planar arc] -> biarc -> P2. Walk FORWARD from P1 by phi.
        gp_Vec ax = tH.Crossed(n);
        ax.Normalize();
        const gp_Pnt J = P1.Rotated(gp_Ax1(C, gp_Dir(ax)), phi);
        gp_Vec tJ = tH.Rotated(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(ax)), phi);
        tJ.Normalize();
        ArcOut a1 = arcFromTangent(P1, t1, J, like, "fillet arc 1");
        if (a1.tEnd.Angle(tJ) > 1e-9) return std::nullopt;
        const auto b = biarc(J, tJ, P2, t2, like);
        if (!b) return std::nullopt;
        out.arcs = {a1, b->a1, b->a2};
        out.arcs[1].pr.label = like.label + " fillet arc 2";
        out.arcs[2].pr.label = like.label + " fillet arc 3";
        out.minRadius = std::min({a1.radius, b->a1.radius, b->a2.radius});
    }
    return out;
}

// The biarc from (P1, t1) to (P2, t2) whose HELIX-side arc lies in the helix's osculating plane
// at the tangency (tangent + principal normal `nHelix`), `helixAtP2` saying which end is the
// helix. An arc in that plane leaves the helix purely radially/tangentially: first-order neutral
// and second-order clearing towards the axial siblings one pitch away, where the equal-tangent
// biarc (which absorbs the level-lead-vs-climbing-helix height change symmetrically) tilted
// 2.1 mrad and put the arc 1.3 nm inside the sibling (11_pushpull, measured). The residual goes
// into the lead-side arc, where the sibling is a tangent length away. The split a1 is the one
// unknown; the plane condition is one equation, root-found by bracketing on a1.
std::optional<Biarc> biarcInHelixPlane(const gp_Pnt& P1, const gp_Vec& t1, const gp_Pnt& P2,
                                       const gp_Vec& t2, bool helixAtP2, const gp_Vec& nHelix,
                                       const Primitive& like) {
    const gp_Vec tH = helixAtP2 ? t2 : t1;
    const gp_Pnt PH = helixAtP2 ? P2 : P1;
    gp_Vec binormal = tH.Crossed(nHelix);
    if (binormal.Magnitude() < 1e-30) return std::nullopt;
    binormal.Normalize();
    // residual: the helix-side arc's chord (from the tangency to the junction) must have no
    // binormal component
    auto residual = [&](double a1) -> std::optional<double> {
        const auto b = biarcSplit(P1, t1, P2, t2, a1, like);
        if (!b) return std::nullopt;
        const Primitive& helixArc = helixAtP2 ? b->a2.pr : b->a1.pr;
        // the junction J: end of arc 1 == start of arc 2; recover it from arc 1's landing
        gp_Pnt J;
        if (b->a1.straight) J = b->a1.pr.seg.b;
        else {
            const gp_Pnt c = b->a1.pr.arc.c;
            const gp_Vec v0(b->a1.pr.arc.v0);
            const gp_Pnt start = c.Translated(v0);
            J = start.Rotated(gp_Ax1(c, gp_Dir(b->a1.pr.arc.axis)), b->a1.pr.arc.sweep);
        }
        (void)helixArc;
        return gp_Vec(PH, J).Dot(binormal);
    };
    const double D = P1.Distance(P2);
    const int n = 64;
    std::optional<double> prevA, prevG;
    std::optional<double> rootA;
    for (int k = 1; k <= n; ++k) {
        const double a1 = D * k / (n + 1.0);
        const auto g = residual(a1);
        if (!g) { prevA.reset(); prevG.reset(); continue; }
        if (std::abs(*g) < 1e-15 * D) { rootA = a1; break; }
        if (prevA && prevG && ((*g < 0) != (*prevG < 0))) {
            double lo = *prevA, hi = a1, glo = *prevG;
            for (int it = 0; it < 200 && hi - lo > 1e-16 * D; ++it) {
                const double mid = 0.5 * (lo + hi);
                const auto gm = residual(mid);
                if (!gm) break;
                if ((*gm < 0) == (glo < 0)) { lo = mid; glo = *gm; }
                else { hi = mid; }
            }
            rootA = 0.5 * (lo + hi);
            break;
        }
        prevA = a1;
        prevG = *g;
    }
    if (!rootA) return std::nullopt;
    auto b = biarcSplit(P1, t1, P2, t2, *rootA, like);
    if (!b) return std::nullopt;
    const auto g = residual(*rootA);
    if (!g || std::abs(*g) > 1e-12 * D) return std::nullopt;
    return b;
}

// The biarc from (P1, t1) to (P2, t2) with equal tangent lengths: control polygon P1, A = P1 + a t1,
// B = P2 - a t2, P2 with |B - A| = 2a; junction J = (A + B) / 2.
std::optional<Biarc> biarc(const gp_Pnt& P1, const gp_Vec& t1, const gp_Pnt& P2, const gp_Vec& t2,
                           const Primitive& like) {
    const gp_Vec d(P1, P2);
    const gp_Vec s = t1 + t2;
    const double c = t1.Dot(t2);
    const double qa = 2.0 * c - 2.0, qb = -2.0 * d.Dot(s), qc = d.SquareMagnitude();
    double a = -1.0;
    if (std::abs(qa) < 1e-12) {
        if (std::abs(qb) < 1e-300) return std::nullopt;
        a = -qc / qb;
    }
    else {
        const double disc = qb * qb - 4.0 * qa * qc;
        if (disc < 0) return std::nullopt;
        const double r1 = (-qb - std::sqrt(disc)) / (2.0 * qa), r2 = (-qb + std::sqrt(disc)) / (2.0 * qa);
        a = std::numeric_limits<double>::max();
        for (double r : {r1, r2})
            if (r > 0 && r < a) a = r;
        if (a == std::numeric_limits<double>::max()) return std::nullopt;
    }
    const gp_Pnt A = P1.Translated(t1 * a), B = P2.Translated(t2 * (-a));
    const gp_Pnt J((A.X() + B.X()) / 2, (A.Y() + B.Y()) / 2, (A.Z() + B.Z()) / 2);
    Biarc out;
    out.a1 = arcFromTangent(P1, t1, J, like, "fillet arc 1");
    gp_Vec tJ(A, B);
    if (tJ.Magnitude() < 1e-15) return std::nullopt;
    tJ.Normalize();
    out.a2 = arcFromTangent(J, tJ, P2, like, "fillet arc 2");
    // Tangency at P2 must come out of the construction to the bit.
    if (out.a2.tEnd.Angle(t2) > 1e-9) return std::nullopt;
    if (loops(out.a1) || loops(out.a2)) return std::nullopt;
    out.minRadius = std::min(out.a1.radius, out.a2.radius);
    return out;
}

// Two SPIRAL pieces that are one and the same helix: same axis, same radial and axial pitch per
// radian, and `a` ends exactly where `b` starts. Then a fillet may run through their boundary.
bool sameHelix(const Primitive& a, const Primitive& b) {
    if (a.kind != Primitive::SPIRAL || b.kind != Primitive::SPIRAL || a.spiral.blend || b.spiral.blend ||
        a.spiral.levelOut || b.spiral.levelOut)
        return false;
    const Spiral& x = a.spiral;
    const Spiral& y = b.spiral;
    if (std::abs(x.cx - y.cx) > 1e-12 || std::abs(x.cz - y.cz) > 1e-12) return false;
    const double sx = x.az1 - x.az0, sy = y.az1 - y.az0;
    if (std::abs(sx) < 1e-15 || std::abs(sy) < 1e-15 || (sx > 0) != (sy > 0)) return false;
    const double rpx = (x.r1 - x.r0) / sx, rpy = (y.r1 - y.r0) / sy;
    const double ypx = (x.y1 - x.y0) / sx, ypy = (y.y1 - y.y0) / sy;
    if (std::abs(rpx - rpy) > 1e-9 * std::max(1.0, std::abs(rpx)) ||
        std::abs(ypx - ypy) > 1e-9 * std::max(1.0, std::abs(ypx)))
        return false;
    const SpiralEval ea = evalSpiral(x, x.az1), eb = evalSpiral(y, y.az0);
    return ea.p.Distance(eb.p) < 1e-12;
}

size_t segIdx0(size_t i, bool exitCorner) { return exitCorner ? i + 1 : i; }

bool isStub(const Primitive& p) {
    return p.kind == Primitive::SPIRAL && !p.spiral.blend && !p.spiral.levelOut &&
           p.label.find("(terminal stub)") != std::string::npos;
}

bool isLayerLink(const Primitive& p) {
    return p.kind == Primitive::SEG && p.label.find("(layer link)") != std::string::npos;
}

bool isWrapSpiral(const Primitive& p) {
    return p.kind == Primitive::SPIRAL && !p.spiral.blend && !p.spiral.levelOut &&
           !p.isConnection;
}

bool isLeadSeg(const Primitive& p) {
    if (p.kind != Primitive::SEG) return false;
    // The lead's own segments, and -- Alf, 2026-09-02 -- the RADIAL LAYER LINK, which meets the
    // wrap at the same kind of corner. The link's chunk is untouched (still the approved straight
    // segment at one azimuth); only the corner where it leaves the wrap is filleted, exactly as
    // the terminal corner is. At stations placed one coated OD apart the sharp corner's true
    // curve-to-curve minimum dips a nanometre or two under the touch value (06_llc 1.3 nm at
    // 13 um along the link, 23_illc 0.5-22 nm); a tangent departure cannot.
    return p.label.find("lead seg") != std::string::npos;
}

// The two corners this pass owns: a terminal stub meeting its lead, and a wrap meeting the
// radial layer link that leaves it. Both are a helix against a straight, and both were sharp.
bool filletable(const Primitive& helixSide, const Primitive& straightSide) {
    if (isStub(helixSide) && isLeadSeg(straightSide)) return true;
    // MVB_LINK_FILLET=1: fillet the wrap-to-layer-link corner too (Alf, 2026-09-02). OFF by
    // default because the corner does not admit one on the designs that need it: 06_llc's link
    // is 0.435 mm long against a 0.42 mm bend radius, so no arc of the wire's own radius fits
    // between the two wraps it joins -- REAL WIRE CANNOT MAKE THAT CORNER EITHER. Either the
    // whole link becomes one smooth arc (a redefinition of the approved U/Z chunk), or MKF
    // gives the radial step room. Recorded on ABT #969; the switch keeps the machinery live.
    if (std::getenv("MVB_LINK_FILLET") != nullptr && isWrapSpiral(helixSide) &&
        isLayerLink(straightSide))
        return true;
    return false;
}

}  // namespace

// THE LAYER LINK IS A SMOOTH TRANSITION (Alf, 2026-09-02, option 2 on ABT #969). The approved
// U/Z chunk drew the radial step as a STRAIGHT segment at one azimuth, so the wire turned 90 deg
// at each of its ends: on 06_llc that link is 0.435 mm long while two bends of the wire's own
// radius need 0.84 mm of tangent -- a corner real wire cannot make, and the sharp corner's true
// curve-to-curve minimum dips a nanometre or two under the sibling's exact-touch envelope
// (1.3 nm on 06, 0.5-22 nm on 23). A winder makes that step by leaning the wire out over a few
// degrees of azimuth. So the link becomes a BIARC tangent to both wraps, taking the azimuth it
// needs from each: stepping one coated OD with bends of the bend radius takes about 0.84 mm of
// travel, which at an 11.5 mm radius is ~4 deg. The straight chunk is gone, not reshaped.
// MVB_STRAIGHT_LINK=1 restores it (bisect switch).
size_t smoothLayerLinks(std::vector<Primitive>& prims, double minBend, const std::string& who,
                        const std::function<bool(const Primitive&, const Primitive&)>& clears,
                        const std::function<int(const std::string&)>& maxStepOf,
                        std::map<std::string, std::pair<int, int>>* stepsOf) {
    if (std::getenv("MVB_STRAIGHT_LINK") != nullptr) return 0;
    size_t done = 0;
    for (size_t i = 0; i + 2 < prims.size(); ++i) {
        if (!isWrapSpiral(prims[i]) || !isLayerLink(prims[i + 1]) || !isWrapSpiral(prims[i + 2]))
            continue;
        Spiral& a = prims[i].spiral;
        Spiral& b = prims[i + 2].spiral;
        const double sgnA = a.az1 >= a.az0 ? 1.0 : -1.0;
        const double sgnB = b.az1 >= b.az0 ? 1.0 : -1.0;
        const double spanA = std::abs(a.az1 - a.az0), spanB = std::abs(b.az1 - b.az0);
        const SpiralEval endA = evalSpiral(a, a.az1);
        std::optional<Biarc> best;
        std::vector<std::pair<Biarc, double>> candidates;   // (transition, azimuth taken per wrap)
        std::vector<int> candidateStep;                     // the step of the azimuth sequence each was found at
        const int maxStep = maxStepOf ? maxStepOf(prims[i + 1].label) : 39;
        double dA = 0.0, dB = 0.0;
        // THE GENTLEST TRANSITION THAT FITS, not the tightest. A winder does not bend at the
        // minimum radius unless the space forces it, and the tightest biarc leans hardest
        // towards the sibling: on 23_illc the first fit left arc 2 grazing the neighbouring
        // parallel by 0.89 nm at exact touch. So keep growing the azimuth taken from each wrap
        // while the arcs keep getting rounder, and take the last fit before the wraps run out.
        for (int it = 0; it < 40 && it <= maxStep; ++it) {
            const double L = minBend * std::tan(0.25 * kPi) * (0.5 + 0.15 * it);
            const double dTry = L / std::max(endA.dsdaz, 1e-18);
            if (dTry > 0.45 * spanA || dTry > 0.45 * spanB) break;
            const double azA = a.az1 - sgnA * dTry, azB = b.az0 + sgnB * dTry;
            const SpiralEval cutA = evalSpiral(a, azA), cutB = evalSpiral(b, azB);
            const auto fit = biarc(cutA.p, cutA.t, cutB.p, cutB.t, prims[i + 1]);
            if (!fit || fit->minRadius < minBend * (1.0 - 1e-9)) {
                if (best) break;   // past the useful range; keep the roundest one found
                continue;
            }
            if (best && fit->minRadius < best->minRadius) break;   // getting tighter again
            best = fit;
            dA = dB = dTry;
            candidates.push_back({*fit, dTry});
            candidateStep.push_back(it);
        }
        // Roundest first, then progressively tighter, until the caller stops vetoing.
        if (clears) {
            bool accepted = false;
            for (size_t ci = candidates.size(); ci-- > 0;) {
                if (!clears(candidates[ci].first.a1.pr, candidates[ci].first.a2.pr)) continue;
                best = candidates[ci].first;
                dA = dB = candidates[ci].second;
                accepted = true;
                break;
            }
            if (!accepted && !candidates.empty()) {
                best = candidates.front().first;   // nothing clears: the tightest, and the gate speaks
                dA = dB = candidates.front().second;
            }
        }
        if (stepsOf && best) {
            int taken = candidateStep.front();
            for (size_t ci = 0; ci < candidates.size(); ++ci)
                if (candidates[ci].second == dA) taken = candidateStep[ci];
            (*stepsOf)[prims[i + 1].label] = {taken, candidateStep.front()};
        }
        if (std::getenv("MVB_LEAD_DIAG") && best) {
            size_t chosen = candidates.size();
            for (size_t ci = 0; ci < candidates.size(); ++ci)
                if (candidates[ci].second == dA) chosen = ci;
            std::fprintf(stderr, "[link-smooth] %s '%s': %zu transitions fit (azimuth %.3f..%.3f deg, min radius %.4f..%.4f mm); "
                         "taken #%zu: azimuth %.3f deg, min radius %.4f mm%s\n", who.c_str(), prims[i + 1].label.c_str(),
                         candidates.size(), candidates.empty() ? 0.0 : candidates.front().second * 180.0 / kPi,
                         candidates.empty() ? 0.0 : candidates.back().second * 180.0 / kPi,
                         candidates.empty() ? 0.0 : candidates.front().first.minRadius * 1e3,
                         candidates.empty() ? 0.0 : candidates.back().first.minRadius * 1e3, chosen, dA * 180.0 / kPi,
                         best->minRadius * 1e3, chosen + 1 < candidates.size() ? " (rounder ones vetoed: they reach a neighbour)" : "");
        }
        if (!best) {
            std::ostringstream m;
            m.precision(9);
            m << "the layer link of " << who << " at '" << prims[i + 1].label
              << "' cannot be made smooth: no biarc of radius >= " << minBend * 1e3
              << " mm (MKF's wire bend) joins the two wraps within 45% of either's azimuth (" << spanA * 180.0 / kPi
              << " and " << spanB * 180.0 / kPi << " deg available)";
            // ABT #2274 (Alf, 2026-11-13: "use round corners everywhere! if not let me know!"): the link stays the
            // straight radial step -- sharp -- and is reported, never bent tighter than MKF's radius.
            std::cerr << "[ConductorBuilder] FINDING: " << m.str() << "; it is drawn sharp (straight radial step)\n";
            if (std::getenv("MVB_CORNER_LIST"))
                std::fprintf(stderr, "[corner-list]\tlayer-link\t%s\t%s\t-\t%.6f\t-\t%.6f\t%.6f\n", who.c_str(),
                             prims[i + 1].label.c_str(), minBend * 1e3, spanA * 180.0 / kPi, spanB * 180.0 / kPi);
            continue;
        }
        auto cut = [](Spiral& sp, double az, bool cutEnd) {
            const double span = sp.az1 - sp.az0;
            const double f = (az - sp.az0) / span;
            const double r = sp.r0 + (sp.r1 - sp.r0) * f;
            const double y = sp.y0 + (sp.y1 - sp.y0) * f;
            if (cutEnd) { sp.az1 = az; sp.r1 = r; sp.y1 = y; }
            else        { sp.az0 = az; sp.r0 = r; sp.y0 = y; }
        };
        cut(a, a.az1 - sgnA * dA, /*cutEnd=*/true);
        cut(b, b.az0 + sgnB * dB, /*cutEnd=*/false);
        Primitive a1 = best->a1.pr, a2 = best->a2.pr;
        a1.label = prims[i + 1].label + " arc 1";
        a2.label = prims[i + 1].label + " arc 2";
        a1.isConnection = a2.isConnection = true;
        a1.turnOrdinal = a2.turnOrdinal = prims[i + 1].turnOrdinal;
        if (std::getenv("MVB_FILLET_DIAG")) {
            std::fprintf(stderr,
                         "[link] '%s' smoothed: %.4f deg taken from each wrap, arc radii %.6f mm "
                         "(bend %.6f)\n",
                         prims[i + 1].label.c_str(), dA * 180.0 / kPi, best->minRadius * 1e3,
                         minBend * 1e3);
        }
        prims.erase(prims.begin() + static_cast<long>(i + 1));
        prims.insert(prims.begin() + static_cast<long>(i + 1), {a1, a2});
        i += 2;
        ++done;
    }
    return done;
}

bool filletSegCorner(std::vector<Primitive>& prims, size_t i, double radius, const Primitive& flags,
                     const std::string& label) {
    Primitive& a = prims[i];
    Primitive& b = prims[i + 1];
    const gp_XYZ va = a.seg.b.XYZ() - a.seg.a.XYZ();
    const gp_XYZ vb = b.seg.b.XYZ() - b.seg.a.XYZ();
    const double lenA = va.Modulus();
    const double lenB = vb.Modulus();
    if (lenA < 1e-12 || lenB < 1e-12) return false;
    const gp_XYZ da = va / lenA;
    const gp_XYZ db = vb / lenB;
    const double cosT = std::max(-1.0, std::min(1.0, da.Dot(db)));
    if (cosT > 1.0 - 1e-12) return true;   // collinear: nothing to round
    const double tanLen = radius * std::sqrt((1.0 - cosT) / (1.0 + cosT));
    if (tanLen >= lenA || tanLen >= lenB) return false;
    const gp_Pnt vertex = a.seg.b;
    const gp_Pnt tA(vertex.XYZ() - da * tanLen);
    const gp_Pnt tB(vertex.XYZ() + db * tanLen);
    gp_XYZ normal = db - da * cosT;
    normal /= normal.Modulus();
    gp_XYZ axis = da.Crossed(db);
    axis /= axis.Modulus();
    Primitive arc;
    arc.kind = Primitive::ARC3;
    arc.arc.c = gp_Pnt(tA.XYZ() + normal * radius);
    arc.arc.axis = axis;
    arc.arc.v0 = tA.XYZ() - arc.arc.c.XYZ();
    arc.arc.sweep = std::acos(cosT);
    arc.label = label;
    arc.turnOrdinal = flags.turnOrdinal;
    arc.isLead = flags.isLead;
    arc.isConnection = flags.isConnection;
    arc.terminal = flags.terminal;
    a.seg.b = tA;
    b.seg.a = tB;
    prims.insert(prims.begin() + static_cast<std::ptrdiff_t>(i) + 1, arc);
    return true;
}

std::optional<SlantedStep> slantOneOdStep(std::vector<Primitive>& prims, size_t s0, bool longAfter,
                                          const gp_Vec& away, double junctionBend, double junctionShare,
                                          double kinkBend) {
    const long k = longAfter ? static_cast<long>(s0) + 1 : static_cast<long>(s0) - 1;
    const long s1 = longAfter ? static_cast<long>(s0) + 2 : static_cast<long>(s0) - 2;
    if (k < 0 || s1 < 0 || s1 >= static_cast<long>(prims.size())) return std::nullopt;
    if (prims[s0].kind != Primitive::SEG || prims[static_cast<size_t>(k)].kind != Primitive::ARC3 ||
        prims[static_cast<size_t>(s1)].kind != Primitive::SEG)
        return std::nullopt;
    const Primitive& sh = prims[s0];
    const Primitive& lg = prims[static_cast<size_t>(s1)];
    const gp_Pnt V = longAfter ? sh.seg.a : sh.seg.b;
    // The long leg, oriented away from the step, and the step's own line, oriented from V to its corner.
    const gp_Pnt nearEnd = longAfter ? lg.seg.a : lg.seg.b;
    const gp_Pnt farEnd = longAfter ? lg.seg.b : lg.seg.a;
    gp_Vec dir1(nearEnd, farEnd);
    const double lenLong = dir1.Magnitude();
    gp_Vec dir0 = longAfter ? gp_Vec(sh.seg.a, sh.seg.b) : gp_Vec(sh.seg.b, sh.seg.a);
    if (lenLong < 1e-12 || dir0.Magnitude() < 1e-12 || away.Magnitude() < 1e-12) return std::nullopt;
    dir1.Normalize();
    dir0.Normalize();
    gp_Vec awayDir = away;
    awayDir.Normalize();
    // The step's corner: where the long leg's line passes nearest the step's line.
    const gp_Vec w(nearEnd, V);
    const double bb = dir0.Dot(dir1);
    const double den = 1.0 - bb * bb;
    if (den < 1e-12) return std::nullopt;   // the step runs along the long leg: no step
    const double tC = (dir1.Dot(w) - bb * dir0.Dot(w)) / den;
    const gp_Pnt C = nearEnd.Translated(dir1 * tC);
    const double reach = lenLong - tC;       // from the corner to the long leg's far end
    if (reach < 1e-12) return std::nullopt;
    // Along the slant from V to Q = C + d dir1: the junction bend's share and the kink's tangent
    // length must both fit; the long leg must keep the kink's tangent length beyond Q.
    auto slack = [&](double d, double& kinkTan) {
        const gp_Pnt Q = C.Translated(dir1 * d);
        gp_Vec u(V, Q);
        const double len = u.Magnitude();
        u.Normalize();
        const double thetaJ = M_PI - awayDir.Angle(u);
        const double thetaQ = u.Angle(dir1);
        kinkTan = kinkBend * std::tan(0.5 * thetaQ);
        return len - junctionBend * std::tan(0.5 * thetaJ) / junctionShare - kinkTan;
    };
    double kinkTan = 0.0;
    if (slack(reach, kinkTan) <= 0.0) return std::nullopt;
    double lo = 0.0, hi = reach;
    for (int it = 0; it < 80; ++it) {
        const double mid = 0.5 * (lo + hi);
        double kt = 0.0;
        (slack(mid, kt) > 0.0 ? hi : lo) = mid;
    }
    slack(hi, kinkTan);
    if (hi + kinkTan >= reach) return std::nullopt;
    const gp_Pnt Q = C.Translated(dir1 * hi);
    SlantedStep out;
    out.slant = hi;
    out.departure = gp_Lin(V, gp_Dir(gp_Vec(V, Q))).Distance(C);
    Primitive slant = sh;
    slant.seg.a = longAfter ? V : Q;
    slant.seg.b = longAfter ? Q : V;
    slant.label = sh.label + " (slanted step)";
    Primitive longFlags = lg;
    if (longAfter) prims[static_cast<size_t>(s1)].seg.a = Q;
    else           prims[static_cast<size_t>(s1)].seg.b = Q;
    prims[s0] = slant;
    prims.erase(prims.begin() + k);
    // [slant][long] or [long][slant]: the kink between them, at kinkBend.
    const size_t kinkAt = longAfter ? s0 : static_cast<size_t>(s1);
    if (!filletSegCorner(prims, kinkAt, kinkBend, longFlags, longFlags.label + " step kink"))
        throw std::logic_error("slantOneOdStep: the kink solved to fit does not fit");
    return out;
}

namespace {
// A Z dragback's own straight (Alf, 2026-10-28: a dragback's ends are round when its corners are).
bool isDragbackSeg(const Primitive& p) {
    return p.kind == Primitive::SEG && p.isConnection && p.label.find("(dragback) seg") != std::string::npos;
}
// Undo the round corner at prims[k] (an arc between two straights): both straights run on to the corner's
// vertex again, where the assembler mitres them. False when prims[k] is no such corner.
bool sharpenSegCorner(std::vector<Primitive>& prims, size_t k) {
    if (k == 0 || k + 1 >= prims.size() || prims[k].kind != Primitive::ARC3 || prims[k - 1].kind != Primitive::SEG ||
        prims[k + 1].kind != Primitive::SEG)
        return false;
    gp_Vec da(prims[k - 1].seg.a, prims[k - 1].seg.b);
    if (da.Magnitude() < 1e-12) return false;
    da.Normalize();
    const double tanLen = prims[k].arc.v0.Modulus() * std::tan(0.5 * std::abs(prims[k].arc.sweep));
    const gp_Pnt vertex = prims[k - 1].seg.b.Translated(da * tanLen);
    prims[k - 1].seg.b = vertex;
    prims[k + 1].seg.a = vertex;
    prims.erase(prims.begin() + static_cast<std::ptrdiff_t>(k));
    return true;
}
}  // namespace

static size_t filletHelixCorners(std::vector<Primitive>& prims, double minBend, double wireRadius,
                                 const std::string& who, double entranceRoll, double exitRoll,
                                 const FilletClears& clears, size_t* uncleared, std::string* witness,
                                 bool dragbackEnds, bool allowSlant,
                                 const std::optional<PlannedClimb>& entranceClimb = std::nullopt,
                                 const std::optional<PlannedClimb>& exitClimb = std::nullopt);

size_t filletTerminalCorners(std::vector<Primitive>& prims, double minBend, double wireRadius,
                             const std::string& who, double entranceRoll, double exitRoll,
                             const FilletClears& clears, size_t* uncleared, std::string* witness,
                             const std::optional<PlannedClimb>& entranceClimb,
                             const std::optional<PlannedClimb>& exitClimb) {
    return filletHelixCorners(prims, minBend, wireRadius, who, entranceRoll, exitRoll, clears, uncleared, witness,
                              /*dragbackEnds=*/false, /*allowSlant=*/true, entranceClimb, exitClimb);
}

namespace {
// The unit travel direction at a primitive's start (atEnd false) or end.
gp_Vec travelDirection(const Primitive& pr, bool atEnd) {
    switch (pr.kind) {
        case Primitive::SEG: return gp_Vec(pr.seg.a, pr.seg.b).Normalized();
        case Primitive::ARC3: {
            gp_Vec v(pr.arc.v0);
            if (atEnd) v.Rotate(gp_Ax1(gp::Origin(), gp_Dir(pr.arc.axis)), pr.arc.sweep);
            return gp_Vec(pr.arc.axis).Crossed(v).Normalized();
        }
        case Primitive::SPIRAL: {
            return evalSpiral(pr.spiral, atEnd ? pr.spiral.az1 : pr.spiral.az0).t;   // already the travel tangent
        }
        default: throw std::runtime_error("TerminalFillet: no travel direction for primitive '" + pr.label + "'");
    }
}
}  // namespace

// Marks a level wrap arc while filletDragbackEnds draws it as a zero-pitch SPIRAL (see there).
static constexpr const char* kLevelArcTag = "\x01level-arc";

size_t filletDragbackEnds(std::vector<Primitive>& prims, double minBend, double wireRadius, const std::string& who,
                          const FilletClears& clears, size_t* uncleared, std::string* witness, bool allowSlant) {
    // ABT #2274 (Alf, 2026-11-13: every corner round): a wrap that ends LEVEL at its dragback slot -- the flat arc
    // carrying it from the crossing plane to the slot, axis along the column -- is a helix of zero pitch. Drawn as an
    // ARC3 it was no wrap to the bend pass, and its corner onto the dragback's radial step stayed a silent 90-degree
    // mitre (06_llc: 10 of them; the coated gate found them 30 um into the enamel of the wrap riding over). It is
    // re-expressed as that exact zero-pitch SPIRAL, so the step bends round it like any other wrap.
    for (size_t i = 0; i < prims.size(); ++i) {
        if (!isDragbackSeg(prims[i])) continue;
        for (const size_t w : {i - 1, i + 1}) {
            if (w >= prims.size()) continue;   // i - 1 wraps round when i == 0
            Primitive& a = prims[w];
            if (a.kind != Primitive::ARC3 || a.isConnection || a.isLead) continue;
            const double axLen = std::sqrt(a.arc.axis.SquareModulus());
            if (!(axLen > 0.0) || std::abs(std::abs(a.arc.axis.Y()) / axLen - 1.0) > 1e-12 || std::abs(a.arc.v0.Y()) > 1e-12)
                continue;   // not a level arc round the column's axis
            Spiral sp;
            sp.cx = a.arc.c.X();
            sp.cz = a.arc.c.Z();
            sp.r0 = sp.r1 = std::hypot(a.arc.v0.X(), a.arc.v0.Z());
            sp.y0 = sp.y1 = a.arc.c.Y();
            sp.az0 = std::atan2(-a.arc.v0.Z(), a.arc.v0.X());   // pos = (cx + r cos az, y, cz - r sin az)
            sp.az1 = sp.az0 + (a.arc.axis.Y() > 0 ? 1.0 : -1.0) * a.arc.sweep;   // rotation about +y advances az
            a.kind = Primitive::SPIRAL;
            a.spiral = sp;
            a.label += kLevelArcTag;
        }
    }
    const size_t done = filletHelixCorners(prims, minBend, wireRadius, who, std::numeric_limits<double>::quiet_NaN(),
                                           std::numeric_limits<double>::quiet_NaN(), clears, uncleared, witness,
                                           /*dragbackEnds=*/true, allowSlant);
    // What the bend left of each level arc is drawn as the level ARC3 it was (the wrap's climb eases into level and
    // meets it tangentially; as a SPIRAL the assembler swept it differently and the two solids overlapped at the
    // crossing plane, "mitred corner did not abut"). The tag is dropped from every label it reached (the bend's arc
    // takes its wrap's label).
    for (auto& pr : prims) {
        const size_t at = pr.label.find(kLevelArcTag);
        if (at == std::string::npos) continue;
        pr.label.erase(at, std::string(kLevelArcTag).size());
        if (pr.kind != Primitive::SPIRAL) continue;
        const Spiral sp = pr.spiral;
        if (sp.y0 != sp.y1 || sp.r0 != sp.r1)
            throw std::logic_error("filletDragbackEnds: a level arc came back from the bend pass no longer level");
        const double sweep = sp.az1 - sp.az0;
        pr.kind = Primitive::ARC3;
        pr.arc.c = gp_Pnt(sp.cx, sp.y0, sp.cz);
        pr.arc.axis = gp_XYZ(0.0, sweep >= 0 ? 1.0 : -1.0, 0.0);
        pr.arc.v0 = gp_XYZ(sp.r0 * std::cos(sp.az0), 0.0, -sp.r0 * std::sin(sp.az0));
        pr.arc.sweep = std::abs(sweep);
    }
    // A dragback end left a mitre (a level wrap arc, or no round bend clear) grows wireRadius x tan(turn/2) into its
    // first or last straight; where the dragback's own round corner left less, that corner is sharp too (Alf,
    // 2026-10-30, "Fully sharp there"; flyback: 0.530 mm left of a 1.1155 mm step for a 0.558 mm mitre).
    for (size_t i = 0; i < prims.size(); ++i) {
        if (!isDragbackSeg(prims[i])) continue;
        for (const bool wrapBefore : {true, false}) {
            const size_t w = wrapBefore ? i - 1 : i + 1;
            const size_t corner = wrapBefore ? i + 1 : i - 1;
            if ((wrapBefore && i == 0) || (!wrapBefore && i + 1 >= prims.size())) continue;
            if (prims[w].isConnection || prims[w].isLead) continue;   // not a wrap
            if (corner >= prims.size() || prims[corner].kind != Primitive::ARC3 || !prims[corner].isConnection) continue;
            const gp_Vec tIn = travelDirection(wrapBefore ? prims[w] : prims[i], true);
            const gp_Vec tOut = travelDirection(wrapBefore ? prims[i] : prims[w], false);
            const double turn = tIn.Angle(tOut);
            if (turn < 1e-9) continue;   // tangent: no mitre
            const double mitre = wireRadius * std::tan(0.5 * turn);
            const double leg = prims[i].seg.a.Distance(prims[i].seg.b);
            if (leg >= mitre) continue;
            const std::string cornerLabel = prims[corner].label;
            if (!sharpenSegCorner(prims, corner)) continue;
            if (!wrapBefore) --i;
            std::cerr << "[TerminalFillet] FINDING: '" << prims[i].label << "' of " << who << " meets its wrap in a mitre "
                      << "that needs " << mitre * 1e3 << " mm of it, and its round corner '" << cornerLabel << "' left "
                      << leg * 1e3 << " mm; that corner is drawn sharp (bisection mitre).\n";
        }
    }
    return done;
}

// ABT #2274 (Alf, 2026-11-13: "use round corners everywhere!"; agent 45): a Z dragback's INTERIOR corners, bent at
// MKF's dragbackBendRadius once the dragback's ENDS have taken their share. appendFilletedPolyline leaves them sharp
// (deferDragbackCorners) because the end fillets (filletDragbackEnds on a helix, filletLeadJunctions on a straight
// face) are drawn later, and the old rule, not knowing their claim, split every leg between two corners in half and
// let a corner next to an end take the whole leg. Here each straight is what is really left: a round end has already
// trimmed it to its tangent point; a SHARP end (mitre) still needs wireRadius tan(turn/2) of it. A corner is round
// where its tangent length bend tan(turn/2) fits both legs' room (a leg's length minus what the corner or mitre at its
// other end takes); where two corners compete for one leg, the one needing less takes it first (ties: the earlier).
// Otherwise it is SHARP (bisection mitre) and reported.
namespace {
size_t filletInteriorCornersOf(std::vector<Primitive>& prims, double bend, double wireRadius, const std::string& who,
                               const std::function<bool(const Primitive&)>& inRun, const char* what, const char* joint);
}

size_t filletDragbackInteriorCorners(std::vector<Primitive>& prims, double bend, double wireRadius, const std::string& who) {
    return filletInteriorCornersOf(prims, bend, wireRadius, who, isDragbackSeg, "dragback", " bend ");
}

// The same pass for a lead whose corners appendFilletedPolyline deferred (Primitive::deferredCorners): run after the
// terminal fillet, so each straight is what that fillet left. A lead's joints are named " corner ", as
// appendFilletedPolyline names them.
size_t filletDeferredLeadCorners(std::vector<Primitive>& prims, double bend, double wireRadius, const std::string& who) {
    const size_t bent = filletInteriorCornersOf(
        prims, bend, wireRadius, who,
        [](const Primitive& p) { return p.kind == Primitive::SEG && p.isLead && p.deferredCorners; }, "lead", " corner ");
    for (auto& p : prims) p.deferredCorners = false;   // settled: round or reported sharp
    return bent;
}

namespace {
size_t filletInteriorCornersOf(std::vector<Primitive>& prims, double bend, double wireRadius, const std::string& who,
                               const std::function<bool(const Primitive&)>& inRun, const char* what, const char* joint) {
    auto dirOf = [](const Primitive& pr) {
        const gp_XYZ d = pr.seg.b.XYZ() - pr.seg.a.XYZ();
        return d / d.Modulus();
    };
    size_t bent = 0;
    for (size_t i0 = 0; i0 < prims.size(); ++i0) {
        if (!inRun(prims[i0])) continue;
        size_t i1 = i0;   // the run [i0, i1]: dragback straights meeting end to end
        while (i1 + 1 < prims.size() && inRun(prims[i1 + 1]) &&
               prims[i1].seg.b.Distance(prims[i1 + 1].seg.a) <= 1e-9)
            ++i1;
        const size_t m = i1 - i0 + 1;
        if (m < 2) { i0 = i1; continue; }
        std::vector<double> len(m), halfTan(m, 0.0);
        for (size_t k = 0; k < m; ++k) len[k] = prims[i0 + k].seg.a.Distance(prims[i0 + k].seg.b);
        for (size_t k = 1; k < m; ++k) {
            const double cosT = std::clamp(dirOf(prims[i0 + k - 1]).Dot(dirOf(prims[i0 + k])), -1.0, 1.0);
            if (cosT > 1.0 - 1e-9) continue;
            const double t = std::sqrt((1.0 - cosT) / (1.0 + cosT));
            if (bend * t < 1e-7) continue;   // OCCT's confusion tolerance: no arc it can build
            halfTan[k] = t;
        }
        // What a SHARP end takes of the run's first / last straight: the mitre against the wrap it meets.
        auto mitreAt = [&](bool atStart) {
            if ((atStart && i0 == 0) || (!atStart && i1 + 1 >= prims.size())) return 0.0;
            const Primitive& wrap = prims[atStart ? i0 - 1 : i1 + 1];
            // Agent 46: measured for EVERY neighbour, by the turn it really makes. A round end (filletDragbackEnds'
            // arc) meets the straight tangentially, so its turn is 0 and it takes nothing more; any other neighbour
            // -- a plain helix, or the level / blended wrap end a round column hands the dragback -- meets it in a
            // mitre that takes wireRadius tan(turn/2) of the straight. (Skipping every non-plain-helix neighbour as
            // "already trimmed" let 06/23's dragback bend 1 eat that mitre and overlap its own wrap, 0.0089 mm^3.)
            auto dirOf = [](const Primitive& pr, bool atEnd) {
                if (pr.kind == Primitive::ARC3)   // P(t) = c + rot(axis, t) v0, t from 0 to sweep (either sign)
                    return pr.arc.sweep < 0.0 ? travelDirection(pr, atEnd).Reversed() : travelDirection(pr, atEnd);
                return travelDirection(pr, atEnd);   // throws for a primitive with no defined travel direction
            };
            const gp_Vec tIn = atStart ? dirOf(wrap, true) : dirOf(prims[i1], true);
            const gp_Vec tOut = atStart ? dirOf(prims[i0], false) : dirOf(wrap, false);
            const double turn = tIn.Angle(tOut);
            return turn < 1e-9 ? 0.0 : wireRadius * std::tan(0.5 * turn);
        };
        std::vector<double> takenAtStart(m, 0.0), takenAtEnd(m, 0.0), tanLen(m, 0.0);
        std::vector<double> cornerRadius(m, bend);   // corner k's radius: MKF's bend, or that plus an end fillet's residue
        takenAtStart.front() = mitreAt(true);
        takenAtEnd.back() = mitreAt(false);
        std::vector<size_t> order;
        for (size_t k = 1; k < m; ++k)
            if (halfTan[k] > 0.0) order.push_back(k);
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return halfTan[a] < halfTan[b]; });
        // ABT #2274 (Alf, 2026-11-13: "I want ALL of them round when we say round"): a short straight between two
        // corners that turn the wire off its direction and back onto it (a dragback's one-OD radial step between two
        // axial legs) cannot hold two 90-degree bends of radius R when it is shorter than 2R. It is drawn as ONE
        // S of two arcs of that same radius R, each turning by alpha with 2R (1 - cos alpha) = the step's offset:
        // round all through, the step's own straight consumed, R sin(alpha) taken from each long leg.
        std::vector<double> sAlpha(m + 1, 0.0);   // sAlpha[k] > 0: corners k and k+1 (around leg k) are one S
        for (size_t k = 1; k + 1 < m; ++k) {
            if (halfTan[k] <= 0.0 || halfTan[k + 1] <= 0.0) continue;
            const gp_XYZ a = dirOf(prims[i0 + k - 1]), u = dirOf(prims[i0 + k]), b = dirOf(prims[i0 + k + 1]);
            if (a.Dot(b) < 1.0 - 1e-9 || std::abs(a.Dot(u)) > 1e-6) continue;   // not an offset step
            if (bend * (halfTan[k] + halfTan[k + 1]) <= len[k]) continue;        // the two bends fit as they are
            const double h = len[k];
            if (!(h < 2.0 * bend)) continue;
            const double alpha = std::acos(1.0 - h / (2.0 * bend));
            const double take = bend * std::sin(alpha);
            if (take > len[k - 1] - takenAtStart[k - 1] || take > len[k + 1] - takenAtEnd[k + 1]) continue;
            sAlpha[k] = alpha;
            tanLen[k] = take;          // what the S takes of the leg before it
            tanLen[k + 1] = take;      // and of the leg after it
            takenAtEnd[k - 1] = take;
            takenAtStart[k + 1] = take;
            halfTan[k] = halfTan[k + 1] = 0.0;   // settled: no plain bend at either corner
            ++k;
        }
        order.erase(std::remove_if(order.begin(), order.end(), [&](size_t k) { return halfTan[k] <= 0.0; }), order.end());
        for (size_t k : order) {
            const double need = bend * halfTan[k];
            const double roomBefore = len[k - 1] - takenAtStart[k - 1];
            const double roomAfter = len[k] - takenAtEnd[k];
            // Equal is a fit: MKF makes a dragback's radial leg exactly one bend long, and need == room then
            // differs only by floating-point round-off (17_cllc: 0.00098225 vs 0.00098225 m). A few ulps, no length.
            auto fits = [](double nd, double rm) { return nd <= rm + 8.0 * std::numeric_limits<double>::epsilon() * std::max(std::abs(nd), std::abs(rm)); };
            if (fits(need, roomBefore) && fits(need, roomAfter)) {
                // TWO BENDS THAT MEET BY PLAN MEET (ABT #2474). A radial step of one wire diameter, bent at half a
                // diameter at both ends, has no room for a straight; drawn at exactly MKF's radius its end fillet on
                // the wrap's CIRCLE stops a few micrometres short of where a flat face would put it, and that residue
                // stood as a 3-9 um straight between the two bends (06_llc: 10 micron-long solids per parallel).
                // filletHelixCorners solved the one radius at which the fillet and this corner meet; the corner is
                // drawn at it, starting where the fillet stops.
                double take = need;
                const std::optional<double> rBefore = k == 1 ? prims[i0].endBendRadius : std::optional<double>();
                const std::optional<double> rAfter = k + 1 == m ? prims[i1].startBendRadius : std::optional<double>();
                if (rBefore.has_value() != rAfter.has_value()) {
                    const double rest = rBefore ? roomBefore : roomAfter, other = rBefore ? roomAfter : roomBefore;
                    const double shared = rBefore ? *rBefore : *rAfter;
                    // The solve closed the leg to 1e-14 m; a corner that does not agree with it to 10 pm was solved
                    // on another leg than the one drawn here.
                    if (std::abs(rest - shared * halfTan[k]) > 1e-11) {
                        std::ostringstream msg;
                        msg << "TerminalFillet: the all-bend step at '" << prims[i0 + k].label << "' of " << who
                            << " does not close: its end fillet leaves " << rest * 1e6 << " um, the corner at the shared radius "
                            << shared * 1e3 << " mm needs " << shared * halfTan[k] * 1e6 << " um.";
                        throw std::runtime_error(msg.str());
                    }
                    if (!fits(rest, other)) {
                        std::ostringstream msg;
                        msg << "TerminalFillet: the all-bend step at '" << prims[i0 + k].label << "' of " << who
                            << " takes " << rest * 1e6 << " um of the next straight, which has only " << other * 1e6 << " um.";
                        throw std::runtime_error(msg.str());
                    }
                    take = rest;
                    cornerRadius[k] = take / halfTan[k];
                    if (std::getenv("MVB_FILLET_DIAG"))
                        std::fprintf(stderr, "[bend90] '%s' corner %zu meets its end fillet: radius %.6f -> %.6f mm (shared %.6f mm)\n",
                                     prims[i0 + k].label.c_str(), k, bend * 1e3, cornerRadius[k] * 1e3, shared * 1e3);
                }
                tanLen[k] = take;
                takenAtEnd[k - 1] = take;
                takenAtStart[k] = take;
                continue;
            }
            const gp_Pnt at = prims[i0 + k].seg.a;
            std::cerr << "[ConductorBuilder] FINDING: " << what << " corner " << k << " of '" << prims[i0 + k].label << "' of "
                      << who << " (at " << at.X() << ", " << at.Y() << ", " << at.Z() << " m) cannot hold MKF's bend of "
                      << bend << " m: its tangent length " << need << " m against legs " << roomBefore << " and "
                      << roomAfter << " m left beside what its neighbours take (" << takenAtStart[k - 1] << " and "
                      << takenAtEnd[k] << " m); it is drawn sharp (bisection mitre)\n";
            if (std::getenv("MVB_CORNER_LIST"))
                std::fprintf(stderr, "[corner-list]\t%s\t%s\t%s\t%.6f,%.6f,%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", what, who.c_str(),
                             prims[i0 + k].label.c_str(), at.X() * 1e3, at.Y() * 1e3, at.Z() * 1e3, bend * 1e3,
                             need * 1e3, roomBefore * 1e3, roomAfter * 1e3);
        }
        std::vector<Primitive> rebuilt;
        for (size_t k = 0; k < m; ++k) {
            if (k >= 1 && sAlpha[k] > 0.0) continue;   // the S's own short step: consumed by its two arcs
            Primitive sg = prims[i0 + k];
            const gp_XYZ d = dirOf(sg);
            const gp_Pnt b0 = sg.seg.b;
            sg.seg.a = gp_Pnt(sg.seg.a.XYZ() + d * tanLen[k]);
            sg.seg.b = gp_Pnt(b0.XYZ() - d * (k + 1 < m ? tanLen[k + 1] : 0.0));
            if (sg.seg.a.Distance(sg.seg.b) >= 1e-12) rebuilt.push_back(sg);
            if (k + 1 < m && sAlpha[k + 1] > 0.0) {
                // The S across leg k+1: from this leg's trimmed end, tangent to it, to the step's midpoint, then on
                // to the next leg's trimmed start -- two arcs of radius `bend` (arcFromTangent derives the radius
                // from the chord; 2R(1 - cos alpha) = the step makes it exactly `bend`).
                const Primitive& step = prims[i0 + k + 1];
                const gp_Pnt mid((step.seg.a.XYZ() + step.seg.b.XYZ()) * 0.5);
                const gp_XYZ dNext = dirOf(prims[i0 + k + 2]);
                const gp_Pnt endS(step.seg.b.XYZ() + dNext * tanLen[k + 2]);
                const std::string& lb = step.label;
                const std::string base = lb.substr(0, lb.rfind(" seg "));
                ArcOut a1 = arcFromTangent(sg.seg.b, gp_Vec(d), mid, step, "s1");
                ArcOut a2 = arcFromTangent(mid, a1.tEnd, endS, step, "s2");
                if (a2.tEnd.Angle(gp_Vec(dNext)) > 1e-7 || std::abs(a1.radius - bend) > 1e-9 * bend ||
                    std::abs(a2.radius - bend) > 1e-9 * bend)
                    throw std::logic_error("filletDragbackInteriorCorners: the S of " + who + " at '" + lb +
                                           "' does not close at radius " + std::to_string(bend));
                a1.pr.label = base + joint + std::to_string(k + 1);
                a2.pr.label = base + joint + std::to_string(k + 2);
                rebuilt.push_back(a1.pr);
                rebuilt.push_back(a2.pr);
                bent += 2;
                continue;
            }
            if (k + 1 >= m || tanLen[k + 1] <= 0.0) continue;
            const gp_XYZ d1 = dirOf(prims[i0 + k + 1]);
            const double cosT = std::clamp(d.Dot(d1), -1.0, 1.0);
            gp_XYZ normal = d1 - d * cosT;
            const gp_XYZ axis = d.Crossed(d1);
            if (normal.Modulus() < 1e-12 || axis.Modulus() < 1e-12)
                throw std::logic_error("filletDragbackInteriorCorners: corner " + std::to_string(k + 1) + " of " + who +
                                       " has a turn but no plane");
            normal /= normal.Modulus();
            const gp_Pnt tA(b0.XYZ() - d * tanLen[k + 1]);
            Primitive pr = prims[i0 + k];
            pr.kind = Primitive::ARC3;
            pr.arc.c = gp_Pnt(tA.XYZ() + normal * cornerRadius[k + 1]);
            pr.arc.axis = axis / axis.Modulus();
            pr.arc.v0 = tA.XYZ() - pr.arc.c.XYZ();
            pr.arc.sweep = std::acos(cosT);
            const std::string& lb = prims[i0 + k].label;
            // " bend " as appendFilletedPolyline names a connection's joint (not " corner ": 02_flyback reads those).
            pr.label = lb.substr(0, lb.rfind(" seg ")) + joint + std::to_string(k + 1);
            rebuilt.push_back(std::move(pr));
            ++bent;
        }
        if (rebuilt.empty())
            throw std::logic_error("filletDragbackInteriorCorners: the dragback run of " + who + " at '" + prims[i0].label +
                                   "' has no copper left");
        prims.erase(prims.begin() + static_cast<std::ptrdiff_t>(i0), prims.begin() + static_cast<std::ptrdiff_t>(i1) + 1);
        prims.insert(prims.begin() + static_cast<std::ptrdiff_t>(i0), rebuilt.begin(), rebuilt.end());
        i0 += rebuilt.size() - 1;
    }
    return bent;
}
}  // namespace

static size_t filletHelixCorners(std::vector<Primitive>& prims, double minBend, double wireRadius,
                                 const std::string& who, double entranceRoll, double exitRoll,
                                 const FilletClears& clears, size_t* uncleared, std::string* witness,
                                 bool dragbackEnds, bool allowSlant,
                                 const std::optional<PlannedClimb>& entranceClimb,
                                 const std::optional<PlannedClimb>& exitClimb) {
    // The corners this pass owns: a terminal stub and its lead, or (dragbackEnds) a wrap and the Z
    // dragback straight it meets -- the wrap it leaves (exit order) or the wrap it lands on (entrance).
    auto cornerOf = [&](const Primitive& helixSide, const Primitive& straightSide) {
        return dragbackEnds ? isWrapSpiral(helixSide) && isDragbackSeg(straightSide)
                            : filletable(helixSide, straightSide);
    };
    size_t done = 0;
    // ABT #2274 (Alf, 2026-11-13): every corner here is ONE planar 90-degree bend at MKF's radius, or sharp and
    // reported -- never a rolled biarc, so the fan's rolls, the slanted step and the per-candidate clearance search
    // have nothing left to choose. The bend is judged with the rest of the copper by the enamel gate.
    (void)entranceRoll; (void)exitRoll; (void)clears; (void)uncleared; (void)witness; (void)allowSlant;
    for (size_t i = 0; i + 1 < prims.size(); ++i) {
        const bool exitCorner = cornerOf(prims[i], prims[i + 1]);
        const bool entranceCorner = cornerOf(prims[i + 1], prims[i]);
        if (!exitCorner && !entranceCorner) continue;
        const size_t stubIdx = exitCorner ? i : i + 1;
        size_t leadIdx = exitCorner ? i + 1 : i;
        const std::string stubLabel = prims[stubIdx].label;
        // ABT #2274 item E (agent 48): the chain's TERMINAL corner (no helix before an entrance corner, none after an
        // exit corner) takes the climb MKF planned for it, when MKF planned one.
        std::optional<PlannedClimb> planned;
        if (!dragbackEnds) {
            bool helixBeyond = false;
            if (entranceCorner)
                for (size_t k = 0; k < i + 1 && !helixBeyond; ++k) helixBeyond = prims[k].kind == Primitive::SPIRAL;
            else
                for (size_t k = i + 1; k < prims.size() && !helixBeyond; ++k) helixBeyond = prims[k].kind == Primitive::SPIRAL;
            if (!helixBeyond) planned = entranceCorner ? entranceClimb : exitClimb;
        }
        Primitive& sp = prims[stubIdx];
        const double sgn = sp.spiral.az1 >= sp.spiral.az0 ? 1.0 : -1.0;
        const double azJoint = exitCorner ? sp.spiral.az1 : sp.spiral.az0;
        const SpiralEval atJoint = evalSpiral(sp.spiral, azJoint);

        // The lead side: the lead leg adjacent to the stub, or -- when that leg is one of MKF's
        // one-OD verticals too short to host both the terminal fillet and its own route corner
        // (boost, 11_pushpull: 0.9 mm legs, bend radius 0.48 mm) -- the leg BEFORE it, the short
        // leg being consumed: the wire makes one continuous bend from the helix onto that leg.
        auto legInfo = [&](size_t k, bool travelTowardJoint, gp_Vec& t, double& len, double& usable) {
            const Primitive& sg = prims[k];
            t = gp_Vec(sg.seg.a, sg.seg.b);
            len = t.Magnitude();
            if (len < 1e-12) { usable = 0.0; return; }
            t.Normalize();   // the stored a->b IS the travel direction (entrance legs run
                             // tip->attach, exit legs attach->tip)
            // The leg's far end is a route corner the assembler mitres (SEG against SEG): the
            // neighbour grows r*tan(phi/2) past the corner point. Only what lies beyond that
            // reach (plus one radius of tube beside it) is available to the fillet.
            usable = len;
            const long other = travelTowardJoint ? static_cast<long>(k) - 1 : static_cast<long>(k) + 1;
            if (other >= 0 && other < static_cast<long>(prims.size()) &&
                prims[static_cast<size_t>(other)].kind == Primitive::SEG) {
                const Primitive& o = prims[static_cast<size_t>(other)];
                gp_Vec tO(o.seg.a, o.seg.b);
                if (tO.Magnitude() > 1e-12) {
                    tO.Normalize();
                    gp_Vec tk(sg.seg.a, sg.seg.b);
                    tk.Normalize();
                    const double phi = tk.Angle(tO);
                    if (phi > 1e-9 && phi < kPi - 1e-9)
                        usable = std::max(0.0, len - (wireRadius * std::tan(0.5 * phi) + wireRadius));
                }
            }
        };
        // travel: entrance = lead then helix (lead travels TOWARD the joint); exit = helix then lead.
        gp_Vec tSeg;
        double segLen = 0.0, usable = 0.0;
        legInfo(leadIdx, /*travelTowardJoint=*/entranceCorner, tSeg, segLen, usable);
        if (segLen < 1e-12) continue;
        const gp_Vec tIn = exitCorner ? atJoint.t : tSeg;
        const gp_Vec tOut = exitCorner ? tSeg : atJoint.t;
        const double theta = tIn.Angle(tOut);
        if (theta < 1e-9) continue;   // already tangent
        if (theta > 170.0 * kPi / 180.0) {
            std::ostringstream m;
            m << "TerminalFillet: the terminal corner of " << who << " at '" << sp.label
              << "' reverses the wire by " << theta * 180.0 / kPi << " deg; no fillet exists.";
            throw std::runtime_error(m.str());
        }
        // ABT #2274 (Alf, 2026-11-13: "always do 90 degree in either radial or tangential direction"; "the inner
        // radius of the turns should be the coating thickness"): the helix-to-lead corner is ONE planar bend at MKF's
        // radius (minBend: the coated half size plus the coating), never a rolled biarc. The helix is cut back by the
        // bend's tangent length L = R tan(theta/2); its tangent line there and the lead's leg meet at a vertex V, so
        // the arc from the cut point to V + d L is a circle of radius R in the plane of the two lines -- no roll, no
        // search. When the stub is shorter than L the cut continues into the same helix's neighbouring piece. When
        // neither the helix nor the leg has the room, the corner stays sharp and says so.
        {
            double R = minBend;   // MKF's radius; an all-bend step's shared radius once that is solved (below)
            const bool nbOk = dragbackEnds ? false
                                           : exitCorner ? (i > 0 && sameHelix(prims[i - 1], prims[i]))
                                                        : (i + 2 < prims.size() && sameHelix(prims[i + 1], prims[i + 2]));
            const size_t nbI = exitCorner ? (i > 0 ? i - 1 : 0) : i + 2;
            const double span0 = std::abs(sp.spiral.az1 - sp.spiral.az0);
            const double nbSpan = nbOk ? std::abs(prims[nbI].spiral.az1 - prims[nbI].spiral.az0) : 0.0;
            Primitive& lg = prims[leadIdx];
            gp_Pnt farEnd = entranceCorner ? lg.seg.a : lg.seg.b;
            double L = R * std::tan(0.5 * theta);
            // TWO BENDS THAT MEET BY PLAN (ABT #2474). A dragback's radial step is one ride over a wire: one diameter,
            // bent at half a diameter at both ends, so on a flat face its two bends would meet. On the wrap's circle,
            // both at MKF's radius, they leave the circle's sagitta over one bend between them (R^2 / 2r: 3-4 um on
            // 06_llc), and a lane that is not the one the next layer's ride was sized for leaves that ride's excess too
            // (another 3-4 um). Such a step is drawn with ONE radius for both bends, the one at which they meet (solved
            // below): a little above MKF's, never below it. What counts as "such a step" is measured, not assumed: the
            // straight the leg would keep at MKF's radius is shorter than the wire's coating is thick, so the bend's
            // inner radius -- the coating, by MKF's rule -- grows by less than half of itself. A longer straight is a
            // straight, and stays.
            bool stepAllBend = false;
            gp_Vec stepNextDir;   // the straight beyond the leg's far corner, in travel direction
            const double coating = minBend - wireRadius;
            if (dragbackEnds && coating > 0.0) {
                const long kNext = exitCorner ? static_cast<long>(leadIdx) + 1 : static_cast<long>(leadIdx) - 1;
                if (kNext >= 0 && kNext < static_cast<long>(prims.size()) && isDragbackSeg(prims[static_cast<size_t>(kNext)])) {
                    const Primitive& nx = prims[static_cast<size_t>(kNext)];
                    gp_Vec tNext(nx.seg.a, nx.seg.b);
                    if (tNext.Magnitude() > 1e-12) {
                        const double turn = tSeg.Angle(tNext);
                        if (turn > 1e-9 && turn < 170.0 * kPi / 180.0) {
                            stepAllBend = true;
                            stepNextDir = tNext.Normalized();
                        }
                    }
                }
            }
            std::optional<double> stepBendRadius;
            double Lv = L;   // the tangent length V was placed at: Q uses the same, so the arc is exactly symmetric
            double azC = azJoint;
            bool intoNb = false, ok = true;
            gp_Pnt Pc, V, Q;
            gp_Vec tc, d;
            // EXACT ON THE LEG'S OWN LINE (ABT #2274): when the helix's tangent and the leg are coplanar (a level wrap
            // and a level step, 06_llc's dragback starts), the bend is solved with the leg unmoved -- the cut point
            // whose tangent line meets the leg's line exactly R tan(theta/2) away. Re-aiming the leg from its far end
            // (the general method below, for skew lines) tilted 06_llc's dragback steps by 0.15 mrad, which turned
            // their next bend 2 nm into the sibling lane.
            bool exact = false;
            auto solveExact = [&]() {
                exact = false;
                L = R * std::tan(0.5 * theta);
                auto cutPoint = [&](double Ltry, gp_Pnt& P, gp_Vec& t, double& az, bool& nb) {
                    const double dAz = Ltry / std::max(atJoint.dsdaz, 1e-18);
                    if (dAz <= 0.98 * span0) {
                        nb = false;
                        az = exitCorner ? azJoint - sgn * dAz : azJoint + sgn * dAz;
                    } else if (nbOk && dAz - span0 <= 0.98 * nbSpan) {
                        nb = true;
                        const Spiral& s2 = prims[nbI].spiral;
                        az = exitCorner ? s2.az1 - sgn * (dAz - span0) : s2.az0 + sgn * (dAz - span0);
                    } else {
                        return false;
                    }
                    const SpiralEval e = evalSpiral(nb ? prims[nbI].spiral : sp.spiral, az);
                    P = e.p;
                    t = e.t;
                    return true;
                };
                // g(L): how far the tangent line at the cut point L back meets the leg's line, less the tangent
                // length the bend there needs. Coplanarity is required (skew distance below a picometre).
                auto g = [&](double Ltry, double& out) {
                    gp_Pnt P; gp_Vec t; double az; bool nb;
                    if (!cutPoint(Ltry, P, t, az, nb)) return false;
                    const gp_Vec n = t.Crossed(tSeg);
                    if (n.SquareMagnitude() < 1e-24) return false;
                    const gp_Vec w(P, farEnd);
                    if (std::abs(w.Dot(n)) / n.Magnitude() > 1e-12) return false;   // skew: no exact planar bend
                    const double s = w.Crossed(tSeg).Dot(n) / n.SquareMagnitude();   // P + t s is on the leg's line
                    const double th = t.Angle(tSeg);
                    if (th > 170.0 * kPi / 180.0) return false;
                    out = (exitCorner ? s : -s) - R * std::tan(0.5 * th);
                    return true;
                };
                double x0 = L, x1 = L * 1.001, g0 = 0.0, g1 = 0.0;
                bool okSolve = g(x0, g0) && g(x1, g1);
                for (int it = 0; okSolve && it < 60 && std::abs(g1) > 1e-14; ++it) {
                    if (std::abs(g1 - g0) < 1e-300) { okSolve = false; break; }
                    const double x2 = x1 - g1 * (x1 - x0) / (g1 - g0);
                    x0 = x1; g0 = g1; x1 = x2;
                    if (!(x1 > 0.0) || !g(x1, g1)) okSolve = false;
                }
                if (okSolve && std::abs(g1) <= 1e-14) {
                    gp_Pnt P; gp_Vec t; double az; bool nb;
                    cutPoint(x1, P, t, az, nb);
                    const double th = t.Angle(tSeg);
                    const double Lx = R * std::tan(0.5 * th);
                    // the vertex: on the tangent line, Lx from the cut point (= on the leg's line, to |g| <= 1e-14)
                    const gp_Pnt Vx = P.Translated(t * (exitCorner ? Lx : -Lx));
                    const double legLen = Vx.Distance(farEnd);
                    if (Lx < legLen) {
                        Pc = P; tc = t; azC = az; intoNb = nb; L = Lx; Lv = Lx; d = tSeg; V = Vx;
                        exact = true;
                    }
                }
            };
            solveExact();
            // SKEW LINES (a pitched helix against a vertical lead): no planar bend keeps the leg's line, so the leg is
            // re-aimed. When its far end is the end of the lead's own round corner, it is re-aimed about that corner's
            // VERTEX and the corner is rebuilt at its new angle (same radius, same pieces): re-aiming about the arc's
            // end left a 0.2-0.7 degree kink between the corner and the leg (06_llc's entrance leads).
            const long kArc = exitCorner ? static_cast<long>(leadIdx) + 1 : static_cast<long>(leadIdx) - 1;
            const long kSeg = exitCorner ? static_cast<long>(leadIdx) + 2 : static_cast<long>(leadIdx) - 2;
            bool farCorner = false;
            gp_Vec farDir;      // the corner's other straight, travel direction
            double farRadius = 0.0;
            if (!exact && kArc >= 0 && kSeg >= 0 && kSeg < static_cast<long>(prims.size()) &&
                prims[static_cast<size_t>(kArc)].kind == Primitive::ARC3 &&
                prims[static_cast<size_t>(kSeg)].kind == Primitive::SEG) {
                const Primitive& ca = prims[static_cast<size_t>(kArc)];
                const Primitive& cs = prims[static_cast<size_t>(kSeg)];
                const auto ce = primEndpoints(ca);
                const bool joins = exitCorner ? (ce.first.Distance(farEnd) < 1e-12 && ce.second.Distance(cs.seg.a) < 1e-12)
                                              : (ce.second.Distance(farEnd) < 1e-12 && ce.first.Distance(cs.seg.b) < 1e-12);
                gp_Vec dv(cs.seg.a, cs.seg.b);
                if (joins && dv.Magnitude() > 1e-12) {
                    dv.Normalize();
                    farRadius = ca.arc.v0.Modulus();
                    const double tl = farRadius * std::tan(0.5 * std::abs(ca.arc.sweep));
                    farDir = dv;
                    farEnd = exitCorner ? cs.seg.a.Translated(dv * (-tl)) : cs.seg.b.Translated(dv * tl);   // the vertex
                    farCorner = true;
                }
            }
            // THE LEAD STAYS IN ITS LANE'S PLANE (Alf, 2026-11-09, decision 2: "every run stands on MKF's lane"). The
            // helix's tangent line never meets the leg's line, so the leg is re-aimed from the bend's vertex V to its far
            // end; with the helix cut exactly one tangent length back, V stood a few micrometres off the lane sideways
            // and the re-aimed leg -- and the lead's own corner rebuilt at its far end -- leaned out of the lane's plane
            // by about a milliradian: on the isolated buck-boost and the complete flyback that corner bulged 1.9 / 1.7 nm
            // into the neighbouring lead, which MKF lays touching. The plane is the lead's own: the leg and the straight
            // beyond its far corner, or the leg and the column axis where the leg is the lead's last. The cut is moved
            // along the helix until V lies in it; the leg then leans only within the plane, and its far corner stays there.
            std::optional<gp_Vec> lanePlaneNormal;
            if (!exact) {
                gp_Vec nrm = farCorner ? tSeg.Crossed(farDir) : tSeg.Crossed(gp_Vec(0.0, 1.0, 0.0));
                if (nrm.Magnitude() > 1e-9) {
                    nrm.Normalize();
                    lanePlaneNormal = nrm;
                }
            }
            auto solveSkew = [&]() {
              ok = true;
              L = R * std::tan(0.5 * theta);
              for (int pass = 0; pass < 4 && ok; ++pass) {
                auto cutAt = [&](double dAzTry, gp_Pnt& P, gp_Vec& t, double& az, bool& nb) {
                    if (dAzTry <= 0.98 * span0) {
                        nb = false;
                        az = exitCorner ? azJoint - sgn * dAzTry : azJoint + sgn * dAzTry;
                    } else if (nbOk && dAzTry - span0 <= 0.98 * nbSpan) {
                        nb = true;
                        const double rest = dAzTry - span0;
                        const Spiral& s2 = prims[nbI].spiral;
                        az = exitCorner ? s2.az1 - sgn * rest : s2.az0 + sgn * rest;
                    } else {
                        return false;
                    }
                    const SpiralEval e = evalSpiral(nb ? prims[nbI].spiral : sp.spiral, az);
                    P = e.p;
                    t = e.t;
                    return true;
                };
                double dAz = L / std::max(atJoint.dsdaz, 1e-18);
                if (lanePlaneNormal) {
                    // off(dAz): how far the vertex of the cut dAz back stands off the lane's plane.
                    auto off = [&](double dAzTry, double& out) {
                        gp_Pnt P; gp_Vec t; double az; bool nb;
                        if (!cutAt(dAzTry, P, t, az, nb)) return false;
                        const gp_Pnt vertex = exitCorner ? P.Translated(t * L) : P.Translated(t * (-L));
                        out = gp_Vec(farEnd, vertex).Dot(*lanePlaneNormal);
                        return true;
                    };
                    double x0 = dAz, x1 = dAz * (1.0 + 1e-3), h0 = 0.0, h1 = 0.0;
                    bool solved = off(x0, h0) && off(x1, h1);
                    for (int it = 0; solved && it < 40 && std::abs(h1) > 1e-15; ++it) {
                        if (std::abs(h1 - h0) < 1e-300) { solved = false; break; }
                        const double x2 = x1 - h1 * (x1 - x0) / (h1 - h0);
                        x0 = x1; h0 = h1; x1 = x2;
                        if (!(x1 > 0.0) || !off(x1, h1)) solved = false;
                    }
                    if (solved && std::abs(h1) <= 1e-15) dAz = x1;   // else: the cut one tangent length back, as before
                }
                if (!cutAt(dAz, Pc, tc, azC, intoNb)) {
                    ok = false;
                    break;
                }
                V = exitCorner ? Pc.Translated(tc * L) : Pc.Translated(tc * (-L));
                Lv = L;
                d = exitCorner ? gp_Vec(V, farEnd) : gp_Vec(farEnd, V);
                const double legLen = d.Magnitude();
                if (legLen < 1e-12) { ok = false; break; }
                d.Normalize();
                const double th = exitCorner ? tc.Angle(d) : d.Angle(tc);
                if (th > 170.0 * kPi / 180.0) { ok = false; break; }
                const double Lnew = R * std::tan(0.5 * th);
                // The leg must hold this bend's tangent and still reach its far end (whose own corner the
                // dragback/lead passes round with their own bends) -- and, re-aimed about a round corner's vertex,
                // that corner's tangent at its new angle too.
                double farTl = 0.0;
                if (farCorner) {
                    const double thc = exitCorner ? d.Angle(farDir) : farDir.Angle(d);
                    farTl = farRadius * std::tan(0.5 * thc);
                }
                if (Lnew + farTl >= legLen) { ok = false; break; }
                if (std::abs(Lnew - L) < 1e-12) { L = Lnew; break; }
                L = Lnew;
              }
            };
            const bool planar = exact;   // which of the two methods this corner is solved by
            auto solveBend = [&]() {
                if (planar) {
                    solveExact();
                    ok = exact;
                } else {
                    solveSkew();
                }
            };
            if (!planar) solveSkew();
            if (stepAllBend && ok) {
                // left(R): what the leg keeps between this fillet and its far corner's bend, both drawn at the radius
                // the fillet really has (the skew method's last step leaves it a hair off the R it was asked for).
                double drawn = 0.0;
                auto left = [&](double Rtry, double& out) {
                    R = Rtry;
                    solveBend();
                    if (!ok) return false;
                    const double th = exitCorner ? tc.Angle(d) : d.Angle(tc);
                    drawn = Lv / std::tan(0.5 * th);
                    out = V.Distance(farEnd) - Lv - drawn * std::tan(0.5 * d.Angle(stepNextDir));
                    return true;
                };
                double r0 = minBend, h0 = 0.0;
                const bool measured = left(r0, h0);
                if (std::getenv("MVB_FILLET_DIAG") && measured)
                    std::fprintf(stderr, "[bend90] step at '%s': %.4f um between its bends at MKF's radius %.6f mm; coating %.4f um\n",
                                 sp.label.c_str(), h0 * 1e6, minBend * 1e3, coating * 1e6);
                if (measured && h0 > 0.0 && h0 < coating) {
                    const double atMkf = h0;
                    double r1 = minBend + h0 / (std::tan(0.5 * d.Angle(stepNextDir)) + std::tan(0.5 * theta)), h1 = 0.0;
                    bool solved = left(r1, h1);
                    for (int it = 0; solved && it < 60 && std::abs(h1) > 1e-14; ++it) {
                        if (std::abs(h1 - h0) < 1e-300) { solved = false; break; }
                        const double r2 = r1 - h1 * (r1 - r0) / (h1 - h0);
                        r0 = r1; h0 = h1; r1 = r2;
                        if (!(r1 >= minBend) || !left(r1, h1)) solved = false;
                    }
                    if (!solved || std::abs(h1) > 1e-14) {
                        std::ostringstream msg;
                        msg << "TerminalFillet: the all-bend step of " << who << " at '" << sp.label << "' leaves "
                            << atMkf * 1e6 << " um between its two bends at MKF's radius " << minBend * 1e3
                            << " mm, and no shared radius at or above it closes them.";
                        throw std::runtime_error(msg.str());
                    }
                    stepBendRadius = drawn;
                    if (std::getenv("MVB_FILLET_DIAG"))
                        std::fprintf(stderr, "[bend90] all-bend step at '%s' (%s): %.4f um between its bends at MKF's radius "
                                     "%.6f mm -> both at %.6f mm (+%.3f %%), nothing between\n", sp.label.c_str(),
                                     planar ? "level" : "pitched", atMkf * 1e6, minBend * 1e3, drawn * 1e3,
                                     (drawn / minBend - 1.0) * 100.0);
                } else {
                    R = minBend;   // nothing to close (or the two bends do not fit at all): MKF's radius, as planned
                    solveBend();
                }
            }
            if (ok) {
                // Equal tangent lengths from V on both lines: the circle is tangent to the helix at Pc and to the
                // leg at Q exactly (its radius is Lv / tan(theta/2), R to the iteration's last step).
                Q = exitCorner ? V.Translated(d * Lv) : V.Translated(d * (-Lv));
                ArcOut a = exitCorner ? arcFromTangent(Pc, tc, Q, sp, "fillet arc 1")
                                      : arcFromTangent(Q, d, Pc, sp, "fillet arc 1");
                Primitive q = a.pr;
                q.isConnection = false;
                q.terminalFillet = isStub(sp);
                q.terminal = lg.terminal;
                if (exitCorner) lg.seg.a = Q; else lg.seg.b = Q;
                if (stepBendRadius) (exitCorner ? lg.endBendRadius : lg.startBendRadius) = stepBendRadius;
                if (farCorner && !exact) {
                    // Rebuild the lead's round corner at the re-aimed leg's angle, about the same vertex.
                    Primitive& ca = prims[static_cast<size_t>(kArc)];
                    Primitive& cs = prims[static_cast<size_t>(kSeg)];
                    const gp_XYZ dIn = exitCorner ? d.XYZ() : farDir.XYZ();     // travel into the corner
                    const gp_XYZ dOut = exitCorner ? farDir.XYZ() : d.XYZ();    // travel out of it
                    const double cosT = std::max(-1.0, std::min(1.0, dIn.Dot(dOut)));
                    const double tl = farRadius * std::sqrt((1.0 - cosT) / (1.0 + cosT));
                    const gp_Pnt tA(farEnd.XYZ() - dIn * tl), tB(farEnd.XYZ() + dOut * tl);
                    if (exitCorner) { lg.seg.b = tA; cs.seg.a = tB; } else { cs.seg.b = tA; lg.seg.a = tB; }
                    gp_XYZ nrm = dOut - dIn * cosT;
                    gp_XYZ ax = dIn.Crossed(dOut);
                    if (nrm.Modulus() < 1e-15 || ax.Modulus() < 1e-15)
                        throw std::logic_error("filletHelixCorners: the re-aimed lead corner of " + who + " is straight");
                    nrm /= nrm.Modulus();
                    ax /= ax.Modulus();
                    ca.arc.c = gp_Pnt(tA.XYZ() + nrm * farRadius);
                    ca.arc.axis = ax;
                    ca.arc.v0 = tA.XYZ() - ca.arc.c.XYZ();
                    ca.arc.sweep = std::acos(cosT);
                }
                auto cutAt = [&](Spiral& s2, double az, bool cutEnd) {
                    const double f = (az - s2.az0) / (s2.az1 - s2.az0);
                    const double r = s2.r0 + (s2.r1 - s2.r0) * f;
                    const double y = s2.y0 + (s2.y1 - s2.y0) * f;
                    if (cutEnd) { s2.az1 = az; s2.r1 = r; s2.y1 = y; }
                    else        { s2.az0 = az; s2.r0 = r; s2.y0 = y; }
                };
                cutAt(intoNb ? prims[nbI].spiral : sp.spiral, azC, /*cutEnd=*/exitCorner);
                size_t insertAt = i + 1;
                // What the bend leaves of the stub, when no longer than the wire is thick, is folded into the same
                // helix's neighbouring piece (exact: one helix, one law). A piece that short cannot be swept or
                // glued cleanly (ABT #1905), and faceted it has no room for the one-chord transition to the round
                // bend (06_llc at --segments 12: "spiral shorter than its round transition (129.4 um)").
                bool folded = false;
                if (!intoNb && nbOk && primLength(prims[stubIdx]) < 2.0 * wireRadius) {
                    Spiral& nb = prims[nbI].spiral;
                    const Spiral& st = prims[stubIdx].spiral;
                    if (exitCorner) { nb.az1 = st.az1; nb.r1 = st.r1; nb.y1 = st.y1; }
                    else            { nb.az0 = st.az0; nb.r0 = st.r0; nb.y0 = st.y0; }
                    folded = true;
                }
                if (intoNb || folded) {
                    prims.erase(prims.begin() + static_cast<long>(stubIdx));
                    if (stubIdx < insertAt) --insertAt;
                }
                prims.insert(prims.begin() + static_cast<long>(insertAt), q);
                i = insertAt;
                ++done;
                if (std::getenv("MVB_FILLET_DIAG"))
                    std::fprintf(stderr, "[bend90] '%s' %s R=%.6f mm L=%.6f mm radius drawn %.6f mm%s cut (%.12g,%.12g,%.12g) onLeg (%.12g,%.12g,%.12g) far (%.12g,%.12g,%.12g)\n",
                                 stubLabel.c_str(), exitCorner ? "EXIT" : "ENTRANCE", R * 1e3, L * 1e3,
                                 a.radius * 1e3, intoNb ? " (into the neighbouring piece)" : "", Pc.X(), Pc.Y(), Pc.Z(),
                                 Q.X(), Q.Y(), Q.Z(), farEnd.X(), farEnd.Y(), farEnd.Z());
                continue;
            }
            // The lead-slot planner's trial lead (who "... (fan)") is not the emitted chain: a corner it cannot bend is
            // that slot's verdict, not a corner of the part, so it is not reported (23_interleaved_llc printed a sharp
            // Secondary corner that the emitted leads, all round, do not have).
            if (who.find(" (fan") != std::string::npos) continue;
            static thread_local std::set<std::string> reportedSharp;
            if (reportedSharp.insert(who + "|" + stubLabel).second)
                std::cerr << "[TerminalFillet] FINDING: the " << (dragbackEnds ? "dragback end" : "terminal corner")
                          << " of " << who << " at '" << stubLabel << "' has no room for its "
                          << R * 1e3 << " mm bend (tangent " << L * 1e3 << " mm; stub "
                          << span0 * 180.0 / kPi << " deg, lead leg " << segLen * 1e3
                          << " mm); it is drawn sharp.\n";
            if (std::getenv("MVB_CORNER_LIST"))
                std::cerr << "[corner-list]\t" << who << "\t" << stubLabel << "\t"
                          << (dragbackEnds ? "dragback-end" : "terminal") << "\twant=" << R * 1e3
                          << "\tleg=" << segLen * 1e3 << "\n";
            continue;
        }
    }
    return done;
}

}  // namespace mvb
