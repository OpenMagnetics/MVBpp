#include <mvb/TerminalFillet.h>

#include <gp_Ax1.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <cmath>
#include <functional>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <optional>
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
    if (a.kind != Primitive::SPIRAL || b.kind != Primitive::SPIRAL || a.spiral.blend || b.spiral.blend)
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
    return p.kind == Primitive::SPIRAL && !p.spiral.blend &&
           p.label.find("(terminal stub)") != std::string::npos;
}

bool isLayerLink(const Primitive& p) {
    return p.kind == Primitive::SEG && p.label.find("(layer link)") != std::string::npos;
}

bool isWrapSpiral(const Primitive& p) {
    return p.kind == Primitive::SPIRAL && !p.spiral.blend && !p.isConnection;
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
                        const std::function<bool(const Primitive&, const Primitive&)>& clears) {
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
        double dA = 0.0, dB = 0.0;
        // THE GENTLEST TRANSITION THAT FITS, not the tightest. A winder does not bend at the
        // minimum radius unless the space forces it, and the tightest biarc leans hardest
        // towards the sibling: on 23_illc the first fit left arc 2 grazing the neighbouring
        // parallel by 0.89 nm at exact touch. So keep growing the azimuth taken from each wrap
        // while the arcs keep getting rounder, and take the last fit before the wraps run out.
        for (int it = 0; it < 40; ++it) {
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
        if (!best) {
            std::ostringstream m;
            m.precision(9);
            m << "TerminalFillet: the layer link of " << who << " at '" << prims[i + 1].label
              << "' cannot be made smooth: no biarc of radius >= " << minBend * 1e3
              << " mm joins the two wraps within 60% of either's azimuth (" << spanA * 180.0 / kPi
              << " and " << spanB * 180.0 / kPi << " deg available).";
            throw std::runtime_error(m.str());
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

size_t filletTerminalCorners(std::vector<Primitive>& prims, double minBend, double wireRadius,
                             const std::string& who, double entranceRoll, double exitRoll) {
    size_t done = 0;
    for (size_t i = 0; i + 1 < prims.size(); ++i) {
        const bool exitCorner = filletable(prims[i], prims[i + 1]);
        const bool entranceCorner = filletable(prims[i + 1], prims[i]);
        if (!exitCorner && !entranceCorner) continue;
        const size_t stubIdx = exitCorner ? i : i + 1;
        size_t leadIdx = exitCorner ? i + 1 : i;
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
        double L0 = std::max(minBend * std::tan(0.5 * theta), 1e-9);
        bool consumeLeg = false;
        if (0.9 * usable < L0) {
            const long prev = entranceCorner ? static_cast<long>(leadIdx) - 1 : static_cast<long>(leadIdx) + 1;
            if (prev >= 0 && prev < static_cast<long>(prims.size()) &&
                prims[static_cast<size_t>(prev)].kind == Primitive::SEG) {
                gp_Vec t2; double len2 = 0.0, usable2 = 0.0;
                legInfo(static_cast<size_t>(prev), entranceCorner, t2, len2, usable2);
                if (len2 > 1e-12) {
                    consumeLeg = true;
                    leadIdx = static_cast<size_t>(prev);
                    tSeg = t2; segLen = len2; usable = usable2;
                }
            }
        }
        Primitive& sg = prims[leadIdx];
        // the lead-side anchor: where the fillet's straight tangent leg starts/ends
        const gp_Pnt anchor = entranceCorner ? sg.seg.b : sg.seg.a;
        const gp_Vec tIn2 = exitCorner ? atJoint.t : tSeg;
        const gp_Vec tOut2 = exitCorner ? tSeg : atJoint.t;
        const double theta2 = tIn2.Angle(tOut2);
        double L = std::max(minBend * std::tan(0.5 * theta2), 1e-9);
        const bool neighbourOk = exitCorner ? (i > 0 && sameHelix(prims[i - 1], prims[i]))
                                            : (i + 2 < prims.size() && sameHelix(prims[i + 1], prims[i + 2]));
        const size_t nbIdx = exitCorner ? (i > 0 ? i - 1 : 0) : i + 2;
        const double spanAbs = std::abs(sp.spiral.az1 - sp.spiral.az0);
        const double nbSpanAbs = neighbourOk ? std::abs(prims[nbIdx].spiral.az1 - prims[nbIdx].spiral.az0) : 0.0;
        if (std::getenv("MVB_FILLET_DIAG")) {
            std::fprintf(stderr, "[fillet] '%s' %s theta=%.4f deg L0=%.6f mm dsdaz=%.6f m/rad span=%.4f deg (neighbour helix %s, %.4f deg) lead leg %.6f mm (%.6f usable)%s\n",
                         sp.label.c_str(), exitCorner ? "EXIT" : "ENTRANCE", theta2 * 180.0 / kPi, L * 1e3,
                         atJoint.dsdaz, spanAbs * 180.0 / kPi, neighbourOk ? "yes" : "no",
                         nbSpanAbs * 180.0 / kPi, segLen * 1e3, usable * 1e3,
                         consumeLeg ? " [short leg consumed, previous leg used]" : "");
        }
        std::optional<Fillet> best, gentle;
        gp_Pnt P1, P2, gentleP1, gentleP2;
        bool intoNeighbour = false, gentleNb = false;
        double azCut = azJoint, gentleCut = azJoint;
        for (int it = 0; it < 24; ++it, L *= 1.2) {
            const double dAz = L / std::max(atJoint.dsdaz, 1e-18);
            if (L > 0.9 * usable) break;
            bool useNb = false;
            if (dAz <= 0.98 * spanAbs) {
                azCut = exitCorner ? azJoint - sgn * dAz : azJoint + sgn * dAz;
            }
            else if (neighbourOk && dAz - spanAbs <= 0.9 * nbSpanAbs) {
                useNb = true;
                const double rest = dAz - spanAbs;
                const Spiral& nb = prims[nbIdx].spiral;
                azCut = exitCorner ? nb.az1 - sgn * rest : nb.az0 + sgn * rest;
            }
            else {
                break;
            }
            const SpiralEval atCut = evalSpiral(useNb ? prims[nbIdx].spiral : sp.spiral, azCut);
            gp_Vec t1, t2;
            if (exitCorner) {
                P1 = atCut.p; t1 = atCut.t;
                P2 = anchor.Translated(tSeg * L); t2 = tSeg;
            }
            else {
                P1 = anchor.Translated(tSeg * (-L)); t1 = tSeg;
                P2 = atCut.p; t2 = atCut.t;
            }
            {
                const Spiral& cutSp = useNb ? prims[nbIdx].spiral : sp.spiral;
                // The helix-side bend turns TOWARDS THE LEAD: in the plane of the helix tangent
                // and the direction from the helix point to the lead's tangent point. For a
                // radial leg that is the osculating plane (radial departure, neutral to the
                // axial siblings); for an axial leg it is the layer's tangent plane (no radial
                // excursion into the next layer -- a bend in the osculating plane there swung
                // 150 um outward on pushpull/06 and 20-100 um into the neighbouring layer).
                gp_Vec nH = spiralPrincipalNormal(cutSp, azCut);
                // THE ROLL: which way, about the helix's own tangent, the corner escapes. The
                // caller gives it when it can see the neighbours; NaN keeps the historical snap
                // to the natural plane nearest the lead.
                const double roll = exitCorner ? exitRoll : entranceRoll;
                {
                    // Snap the bend direction to one of the helix's two natural planes: the
                    // osculating plane (T, N) when the lead lies radially, the rectifying plane
                    // (T, B) when it lies axially. Using the raw direction to the lead point
                    // carried the leg's few-um axial offset into the plane -> a 2.25 mrad tilt
                    // of the helix-side arc -> 1.3 nm into the axial sibling (boost, measured).
                    // The residual belongs to the lead-side arcs of the biarc.
                    const gp_Pnt& PHx = exitCorner ? P1 : P2;
                    const gp_Pnt& PLx = exitCorner ? P2 : P1;
                    const gp_Vec& tHx = exitCorner ? t1 : t2;
                    gp_Vec toLead(PHx, PLx);
                    toLead -= tHx * toLead.Dot(tHx);
                    if (toLead.Magnitude() > 1e-12 && nH.Magnitude() > 1e-30) {
                        toLead.Normalize();
                        gp_Vec B = tHx.Crossed(nH);
                        B.Normalize();
                        if (std::isfinite(roll)) {
                            // Exactly the direction the caller solved for.
                            nH = nH * std::cos(roll) + B * std::sin(roll);
                        }
                        else {
                            const double cN = toLead.Dot(nH), cB = toLead.Dot(B);
                            nH = (std::abs(cN) >= std::abs(cB)) ? nH * (cN >= 0 ? 1.0 : -1.0)
                                                                : B * (cB >= 0 ? 1.0 : -1.0);
                        }
                    }
                }
                best.reset();
                // PLANAR CORNER: when the lead leaves in the helix's osculating plane (t1, t2 and
                // the chord coplanar), one circular arc of radius L/tan(theta/2) is the exact
                // fillet -- no tilt towards the axial siblings, nothing to absorb.
                {
                    const gp_Vec chord(P1, P2);
                    gp_Vec nPlane = t1.Crossed(t2);
                    if (nPlane.Magnitude() > 1e-12 && chord.Magnitude() > 1e-15) {
                        nPlane.Normalize();
                        const double off = std::abs(chord.Dot(nPlane));
                        const double th = t1.Angle(t2);
                        // tangent lengths must be equal for a single arc: |chord| = 2 L sin... check
                        // the arc through P1 tangent t1 lands on P2 with tangent t2
                        if (off < 1e-12 && th > 1e-9) {
                            // Corner point X of the two tangent lines; the single arc needs equal
                            // tangent lengths from X. The helix side is fixed (P2 on the helix);
                            // slide the lead-side point along its line to match.
                            const gp_Vec w(P2, P1);
                            const double c12 = t1.Dot(t2);
                            const double den = 1.0 - c12 * c12;
                            if (den > 1e-18) {
                                // P1 + s1 t1 = P2 + s2 t2 (least squares on the coplanar lines)
                                const double s1 = (-(w.Dot(t1)) + c12 * w.Dot(t2)) / den;
                                const double s2 = (w.Dot(t2) - c12 * w.Dot(t1)) / den;
                                const gp_Pnt X = P1.Translated(t1 * s1);
                                // a corner: P1 before X along t1, P2 after X along t2
                                if (s1 > 0 && s2 < 0) {
                                    gp_Pnt P1x = P1, P2x = P2;
                                    if (exitCorner) P2x = X.Translated(t2 * s1);      // helix at P1
                                    else            P1x = X.Translated(t1 * s2);      // helix at P2 (s2 < 0)
                                    const gp_Pnt& leadPt = exitCorner ? P2x : P1x;
                                    // the lead-side point must stay on the usable part of the leg
                                    const double along = exitCorner ? gp_Vec(anchor, leadPt).Dot(tSeg)
                                                                    : -gp_Vec(anchor, leadPt).Dot(tSeg);
                                    if (along > 0 && along <= 0.9 * usable) {
                                        ArcOut a = arcFromTangent(P1x, t1, P2x, sp, "fillet arc 1");
                                        if (!a.straight && a.tEnd.Angle(t2) < 1e-9 &&
                                            a.radius >= minBend * (1.0 - 1e-9)) {
                                            P1 = P1x;
                                            P2 = P2x;
                                            Fillet f;
                                            f.arcs = {a};
                                            f.minRadius = a.radius;
                                            best = std::move(f);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                // The planar helix-side bend, from a SHALLOW angle up: its excursion from the
                // helix is R(1-cos phi), so the smallest phi that closes the biarc is the one
                // that stays nearest the wire's own path -- and the neighbours, whichever side
                // they are on, sit one coated OD away with nothing to spare. (Measured the other
                // way round on 14_dab: forcing the GENTLEST fit instead of the tightest took the
                // intrusion from 4 nm to 551 nm.) 1 deg steps, so a 5 deg corner is not rounded
                // to 5 deg of excursion when 1 deg would do.
                const int wantVariant = 0;
                int seen = 0;
                std::optional<Fillet> lastFit;
                // Spread the candidates across the whole admissible range rather than in 1 deg
                // steps: six neighbouring angles are six almost identical fillets, and the fan
                // needs genuinely different ones to choose from (14_dab).
                // VARIANT 0 IS THE HISTORICAL FILLET: the first fit walking 5 degrees at a time,
                // which is what the corpus was verified on. The finer and wider angles are extra
                // CHOICES for the fan, never a change of default -- moving variant 0 to a 1 degree
                // first step silently re-cornered every terminal in the corpus and broke two
                // designs that had been clean (24_margin_interleaved_flyback, complete_flyback).
                static const double kPhiDeg[] = {5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 85};
                for (double phiDeg : kPhiDeg) {
                    const double phi = phiDeg * kPi / 180.0;
                    auto f = triarcHelixPlane(P1, t1, P2, t2, /*helixAtP2=*/!exitCorner, nH,
                                              minBend, phi, sp);
                    if (!(f && f->minRadius >= minBend * (1.0 - 1e-9))) continue;
                    lastFit = f;
                    if (seen++ == wantVariant) { best = std::move(f); break; }
                }
                if (!best && lastFit) best = std::move(lastFit);   // fewer variants than asked
            }
            if (std::getenv("MVB_FILLET_DIAG")) {
                std::fprintf(stderr, "[fillet] '%s' it=%d L=%.6f mm dAz=%.4f deg %s%s minR=%.6f mm (need %.6f)\n",
                             sp.label.c_str(), it, L * 1e3, dAz * 180.0 / kPi, best ? "biarc" : "NO-BIARC",
                             useNb ? " [into neighbour]" : "", best ? best->minRadius * 1e3 : -1.0, minBend * 1e3);
            }
            if (best && best->minRadius >= minBend * (1.0 - 1e-9)) {
                intoNeighbour = useNb;
                // MVB_FILLET_GENTLE=1: keep growing, and take the LAST fit rather than the
                // first. The tightest fillet leans hardest towards whatever the corner turns
                // past -- the same effect the layer link showed against an exact-touch sibling.
                if (std::getenv("MVB_FILLET_GENTLE") == nullptr) break;
                gentle = best;
                gentleCut = azCut;
                gentleP1 = P1;
                gentleP2 = P2;
                gentleNb = useNb;
                best.reset();
                continue;
            }
            best.reset();
        }
        if (!best && gentle) {   // the loop ran past the last fit: take it
            best = gentle;
            azCut = gentleCut;
            P1 = gentleP1;
            P2 = gentleP2;
            intoNeighbour = gentleNb;
        }
        if (!best) {
            std::ostringstream m;
            m.precision(9);
            m << "TerminalFillet: no tangent biarc with radius >= " << minBend * 1e3
              << " mm fits the terminal corner of " << who << " at '" << sp.label << "' ("
              << theta2 * 180.0 / kPi << " deg between the helix and the lead; lead leg "
              << segLen * 1e3 << " mm (" << usable * 1e3 << " mm beyond its own corner"
              << (consumeLeg ? ", after consuming the short leg" : "") << "), stub span "
              << spanAbs * 180.0 / kPi << " deg"
              << (neighbourOk ? ", same-helix neighbour available" : ", no same-helix neighbour") << ").";
            throw std::runtime_error(m.str());
        }
        auto cutSpiral = [&](Spiral& s2, double az, bool cutEnd) {
            const double span = s2.az1 - s2.az0;
            const double f = (az - s2.az0) / span;
            const double r = s2.r0 + (s2.r1 - s2.r0) * f;
            const double y = s2.y0 + (s2.y1 - s2.y0) * f;
            if (cutEnd) { s2.az1 = az; s2.r1 = r; s2.y1 = y; }
            else        { s2.az0 = az; s2.r0 = r; s2.y0 = y; }
        };
        std::vector<Primitive> arcs;
        // ABT #1215: a terminal corner's arcs are tagged for the lead-length measurement (and only
        // for it); a wrap-to-layer-link corner (MVB_LINK_FILLET) is not a terminal and stays untagged.
        const bool terminalCorner = isStub(sp);
        for (const auto& ao : best->arcs) {
            Primitive q = ao.pr;
            q.isConnection = false;
            q.terminalFillet = terminalCorner;
            arcs.push_back(q);
        }
        // Apply. Order of operations keeps indices valid: shorten pieces in place first, then
        // erase consumed pieces from the higher index down, then insert the arcs at the joint.
        if (exitCorner) sg.seg.a = P2; else sg.seg.b = P1;
        if (intoNeighbour) cutSpiral(prims[nbIdx].spiral, azCut, /*cutEnd=*/exitCorner);
        else               cutSpiral(sp.spiral, azCut, /*cutEnd=*/exitCorner);
        // pieces to erase: the stub when consumed, the short lead leg when consumed
        std::vector<size_t> erase;
        if (intoNeighbour) erase.push_back(stubIdx);
        if (consumeLeg) erase.push_back(exitCorner ? i + 1 : i);
        std::sort(erase.begin(), erase.end(), std::greater<size_t>());
        size_t insertAt = exitCorner ? i + 1 : i + 1;   // between the (remaining) helix and lead
        for (size_t k : erase) {
            prims.erase(prims.begin() + static_cast<long>(k));
            if (k < insertAt) --insertAt;
        }
        prims.insert(prims.begin() + static_cast<long>(insertAt), arcs.begin(), arcs.end());
        i = insertAt + arcs.size() - 1;
        ++done;
    }
    return done;
}

}  // namespace mvb
