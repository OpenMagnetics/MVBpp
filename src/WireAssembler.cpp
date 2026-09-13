#include "mvb/WireAssembler.h"
#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"
#include <array>
#include <functional>
#include <set>
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <BRep_Tool.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <Geom_BezierCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <cstdio>
#include <BRepAlgoAPI_Cut.hxx>
#include <ShapeFix_ShapeTolerance.hxx>
#include <Precision.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <ShapeUpgrade_ShapeDivideClosed.hxx>
#include <gp_Pln.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <gp_Elips.hxx>
#include <Precision.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepOffsetAPI_MakePipeShell.hxx>
#include <BRepLib.hxx>
#include <GCE2d_MakeSegment.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <gp_Vec.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Circ.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Dir.hxx>
#include <gp_XY.hxx>
#include <gp_XYZ.hxx>
#include <BRepAlgoAPI_Check.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BOPAlgo_GlueEnum.hxx>
#include <TopTools_ListOfShape.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Solid.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Wire.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <iostream>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace mvb {

int curveSampleCount(double radius, double azSpan, double wireRadius) {
    double maxSag = samplingSag(wireRadius);
    double stepAz = (radius > 1e-12)
                        ? 2.0 * std::acos(std::clamp(1.0 - maxSag / radius, 0.0, 1.0))
                        : std::max(azSpan, 1e-3);
    stepAz = std::clamp(stepAz, 1e-3, 0.2);
    return std::max(2, static_cast<int>(std::ceil(azSpan / stepAz)) + 1);
}

// A SPIRAL'S SAMPLING IS NOT A CIRCLE'S. ABT #373 (:818, Alf 2026-08-23). curveSampleCount
// sizes its step from the AZIMUTHAL geometry alone -- radius x azimuth -- which is exactly
// right for an arc or a helix and badly wrong for a spiral whose RADIUS runs while its azimuth
// barely turns. Measured on the FEM toroid's face crossing (realwinding_toroid, wrap
// 'turn 6' -> 'turn 7' bottom chord): 7.98 mm of copper travelled across 0.0377 rad, radius
// 19.978 -> 12.023 mm. The circle rule saw a 0.86 mm arc and returned TWO samples, so the
// "spiral" was emitted as a straight two-point chord -- the very shape faceSpiral exists to
// replace (ABT #685) -- and, worse, its swept end tangents were then the spiral's START
// tangent at BOTH ends while the assembler mitred the joint against the analytic END tangent.
// The junction reported 3.6e-11 deg (tangent, bridged, uncut) and was really a 4.21 deg corner:
// the chord and the following bottom-inner-corner arc interpenetrated by (2/3) r^3 tan(4.21 deg)
// = 0.0051 mm^3 (measured 0.00512 mm^3, 3050/216000 strict-IN grid points).
//
// The step therefore comes from the curve's OWN speed and curvature. For
// P(u) = (cx + r cos u, y, cz - r sin u) with r' = k and y' = m,
//     |P'| = sqrt(k^2 + r^2 + m^2)
//     kappa = sqrt(m^2 (4k^2 + r^2) + (2k^2 + r^2)^2) / (k^2 + r^2 + m^2)^(3/2)
// and a chord of length L sags by L^2 kappa / 8, so the admissible parameter step is
// sqrt(8 sag / kappa) / |P'|. For k = m = 0 that reduces to sqrt(8 sag / r), which is the
// small-angle form of curveSampleCount's own 2 acos(1 - sag/r): arcs and helices are sampled
// EXACTLY as before. The circle rule is kept as a floor so no primitive anywhere can come out
// coarser than it does today; the same [1e-3, 0.2] step clamp applies, for the same reason.
int spiralSampleCount(const Spiral& sp, double wireRadius) {
    const double dAz = sp.az1 - sp.az0;
    const double azSpan = std::abs(dAz);
    int n = curveSampleCount(std::max(sp.r0, sp.r1), azSpan, wireRadius);
    if (azSpan < 1e-12) return n;
    const double maxSag = samplingSag(wireRadius);
    // The blend's r and y advance with f'(t), so k and m are read LOCALLY at each station.
    constexpr int kStations = 9;
    double stepAz = std::numeric_limits<double>::max();
    for (int i = 0; i < kStations; ++i) {
        const double t = static_cast<double>(i) / (kStations - 1);
        const double f = sp.blend ? 0.5 * (1.0 - std::cos(kPi * t)) : t;
        const double fp = sp.blend ? 0.5 * kPi * std::sin(kPi * t) : 1.0;
        const double r = sp.r0 + (sp.r1 - sp.r0) * f;
        const double k = (sp.r1 - sp.r0) * fp / dAz;
        const double m = (sp.y1 - sp.y0) * fp / dAz;
        const double speed2 = k * k + r * r + m * m;
        if (speed2 < 1e-24) continue;
        const double quad = 2.0 * k * k + r * r;
        const double kappa =
            std::sqrt(m * m * (4.0 * k * k + r * r) + quad * quad) / (speed2 * std::sqrt(speed2));
        if (kappa < 1e-12) continue;   // locally straight: no sag to bound
        stepAz = std::min(stepAz, std::sqrt(8.0 * maxSag / kappa) / std::sqrt(speed2));
    }
    if (stepAz == std::numeric_limits<double>::max()) return n;
    stepAz = std::clamp(stepAz, 1e-3, 0.2);
    return std::max(n, static_cast<int>(std::ceil(azSpan / stepAz)) + 1);
}

std::vector<gp_Pnt> samplePrim(const Primitive& p, double wireRadius) {
    if (p.kind == Primitive::SEG) return {p.seg.a, p.seg.b};
    if (p.kind == Primitive::ARC3) {
        const Arc3& a = p.arc;
        double radius = a.v0.Modulus();
        int n = curveSampleCount(radius, std::abs(a.sweep), wireRadius);
        std::vector<gp_Pnt> pts;
        pts.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            double t = a.sweep * i / (n - 1);
            pts.push_back(gp_Pnt(a.c.XYZ() + rotateXYZ(a.v0, a.axis, t)));
        }
        return pts;
    }
    if (p.kind == Primitive::BLEND) {
        const Blend& bl = p.blendc;
        gp_XYZ d = bl.b.XYZ() - bl.a.XYZ();
        double L = d.Dot(bl.u);
        gp_XYZ perp = d - bl.u * L;
        double perpM = perp.Modulus();
        if (perpM < 1e-12) return {bl.a, bl.b};
        // Peak curvature of the cosine blend: kappa = perp * pi^2 / (2 L^2); sample with
        // the same sag rule as the arcs, on the tightest osculating radius.
        double rMin = 2.0 * L * L / (kPi * kPi * perpM);
        double maxSag = samplingSag(wireRadius);
        double step = 2.0 * std::sqrt(std::max(1e-18, 2.0 * rMin * maxSag));
        double arcLen = std::sqrt(L * L + perpM * perpM) * 1.1;
        int n = std::clamp(static_cast<int>(std::ceil(arcLen / step)) + 1, 8, 512);
        std::vector<gp_Pnt> pts;
        pts.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            double t = static_cast<double>(i) / (n - 1);
            pts.push_back(gp_Pnt(bl.a.XYZ() + bl.u * (L * t) +
                                 perp * (0.5 * (1.0 - std::cos(kPi * t)))));
        }
        return pts;
    }
    const Spiral& sp = p.spiral;
    int n = spiralSampleCount(sp, wireRadius);
    std::vector<gp_Pnt> pts;
    pts.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        double t = static_cast<double>(i) / (n - 1);
        double f = sp.blend ? 0.5 * (1.0 - std::cos(kPi * t)) : t;
        pts.push_back(azPointC(sp.cx, sp.cz, sp.r0 + (sp.r1 - sp.r0) * f,
                               sp.y0 + (sp.y1 - sp.y0) * f, sp.az0 + (sp.az1 - sp.az0) * t));
    }
    return pts;
}

std::pair<gp_Pnt, gp_Pnt> primEndpoints(const Primitive& p) {
    if (p.kind == Primitive::ARC3) {
        return {gp_Pnt(p.arc.c.XYZ() + p.arc.v0),
                gp_Pnt(p.arc.c.XYZ() + rotateXYZ(p.arc.v0, p.arc.axis, p.arc.sweep))};
    }
    if (p.kind == Primitive::SPIRAL) {
        return {azPointC(p.spiral.cx, p.spiral.cz, p.spiral.r0, p.spiral.y0, p.spiral.az0),
                azPointC(p.spiral.cx, p.spiral.cz, p.spiral.r1, p.spiral.y1, p.spiral.az1)};
    }
    if (p.kind == Primitive::BLEND) return {p.blendc.a, p.blendc.b};
    return {p.seg.a, p.seg.b};
}

double primLength(const Primitive& p) {
    switch (p.kind) {
        case Primitive::SEG:
            return p.seg.a.Distance(p.seg.b);
        case Primitive::ARC3: {
            // P(t) = c + rot(axis, t) v0: the point runs on a circle whose radius is v0's
            // component perpendicular to the (unit) axis, at unit angular rate.
            const double axisNorm = p.arc.axis.Modulus();
            if (!(axisNorm > 0.0))
                throw std::runtime_error("primLength: ARC3 '" + p.label + "' has a zero axis");
            const gp_XYZ axis = p.arc.axis / axisNorm;
            const gp_XYZ radial = p.arc.v0 - axis * axis.Dot(p.arc.v0);
            return radial.Modulus() * std::abs(p.arc.sweep);
        }
        case Primitive::SPIRAL: {
            const Spiral& sp = p.spiral;
            if (sp.blend)
                throw std::runtime_error(
                    "primLength: '" + p.label + "' is a BLENDED spiral, whose length has no closed "
                    "form; refusing to return an estimate");
            // r(t) = r0 + dr t, y(t) = y0 + dy t, az(t) = az0 + daz t, t in [0, 1]:
            // |P'(t)|^2 = dr^2 + dy^2 + (daz r(t))^2.
            const double dr = sp.r1 - sp.r0, dy = sp.y1 - sp.y0, daz = std::abs(sp.az1 - sp.az0);
            const double a = dr * dr + dy * dy;
            if (daz == 0.0) return std::sqrt(a);
            if (dr == 0.0) return std::sqrt(a + daz * daz * sp.r0 * sp.r0);
            // u = daz r, dt = du / (daz dr), and a > 0 because dr != 0:
            // integral sqrt(a + u^2) du = u/2 sqrt(a + u^2) + a/2 asinh(u / sqrt(a)).
            auto F = [&](double r) {
                const double u = daz * r;
                return 0.5 * u * std::sqrt(a + u * u) + 0.5 * a * std::asinh(u / std::sqrt(a));
            };
            return (F(sp.r1) - F(sp.r0)) / (daz * dr);
        }
        case Primitive::BLEND:
            throw std::runtime_error(
                "primLength: '" + p.label + "' is a BLEND (cosine S-curve), whose length has no "
                "closed form; refusing to return an estimate");
    }
    throw std::runtime_error("primLength: unknown primitive kind on '" + p.label + "'");
}

// --- geometry emission ------------------------------------------------------------------
// The section frame every round profile must share. gp_Ax2(center, normal) lets OCC pick the X
// direction, so two primitives meeting with the SAME tangent can still get sections rotated
// relative to each other; their abutting faces then only partially coincide. Deriving X from a
// fixed world reference makes the orientation a pure function of the tangent, so tangent
// neighbours share an EXACTLY coincident section -- which is what conformal abutment requires.
static gp_Ax2 deterministicSectionFrame(const gp_Pnt& center, const gp_Dir& normal) {
    gp_XYZ ref(0, 1, 0);
    if (std::abs(normal.XYZ().Dot(ref)) > 0.9) ref = gp_XYZ(1, 0, 0);
    const gp_XYZ xdir = ref - normal.XYZ() * normal.XYZ().Dot(ref);
    return gp_Ax2(center, normal, gp_Dir(xdir));
}

TopoDS_Wire wireProfileWire(const gp_Pnt& center, const gp_Dir& normal, double radius,
                            int segments) {
    // Deterministic frame (see deterministicSectionFrame): a mitre junction whose two sides
    // disagree on the section's X direction leaves faces that only partially coincide -- measured
    // as a 56-primitive rect-column conductor falling into 7 components with no endpoint gaps and
    // nothing dropped.
    gp_Ax2 plane = deterministicSectionFrame(center, normal);
    if (segments <= 0) {
        gp_Circ circ(plane, radius);
        return BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire();
    }
    BRepBuilderAPI_MakePolygon poly;
    gp_Dir dx = plane.XDirection();
    gp_Dir dy = plane.YDirection();
    double offset = kPi / segments;
    // INSCRIBED n-gon (vertices on the nominal circle): the faceted wire never protrudes past
    // the wire's real envelope, which is what the zero-tolerance pairwise-overlap contracts rely
    // on. Its flats consequently sit at the apothem, r*cos(pi/n) -- 0.981*r at n=16 -- so any
    // clearance or coverage check against a faceted section must use the apothem, not r.
    for (int i = 0; i < segments; ++i) {
        double ang = kTwoPi * i / segments + offset;
        gp_XYZ p = center.XYZ() + dx.XYZ() * (radius * std::cos(ang)) +
                   dy.XYZ() * (radius * std::sin(ang));
        poly.Add(gp_Pnt(p));
    }
    poly.Close();
    return poly.Wire();
}

TopoDS_Face wireProfile(const gp_Pnt& center, const gp_Dir& normal, double radius,
                        int segments) {
    return BRepBuilderAPI_MakeFace(wireProfileWire(center, normal, radius, segments)).Face();
}

// wireProfileWire with an explicit phase: the polygon rotated by `phase0` about the axis.
// Used to CONTINUE a neighbour's section phase through a junction: the sweep keeps its exact
// nominal on-plane profile (feeding an off-plane neighbour wire to MakePipeShell re-anchors
// the sweep and corrupts it -- measured: an invalid 121-face solid and -0.22 mm^3 of copper
// on 10_emi), while the vertex azimuths line up with the section the copper actually arrives
// with, which is what removes the phase-mismatch lens at the junction.
static TopoDS_Wire wireProfileWirePhased(const gp_Pnt& center, const gp_Dir& normal,
                                         double radius, int segments, double phase0) {
    gp_Ax2 plane = deterministicSectionFrame(center, normal);
    BRepBuilderAPI_MakePolygon poly;
    gp_Dir dx = plane.XDirection();
    gp_Dir dy = plane.YDirection();
    const double offset = kPi / segments;
    for (int i = 0; i < segments; ++i) {
        const double ang = kTwoPi * i / segments + offset + phase0;
        gp_XYZ p = center.XYZ() + dx.XYZ() * (radius * std::cos(ang)) +
                   dy.XYZ() * (radius * std::sin(ang));
        poly.Add(gp_Pnt(p));
    }
    poly.Close();
    return poly.Wire();
}

// The neighbour section's phase in this piece's own deterministic frame at (center, normal):
// the angle of its first vertex, minus the polygon's base offset. Returns false when the wire
// is unusable (wrong vertex count, vertex on the axis).
// Kill-switches for the section handoff, for bisecting a mesh regression against a build that
// predates it: MVB_NO_SECTION_HANDOFF disables both halves, MVB_NO_CAP_POSITIONS keeps the
// phase continuation but never adopts a neighbour's vertex POSITIONS.
static bool handoffDisabled() {
    static const bool off = std::getenv("MVB_NO_SECTION_HANDOFF") != nullptr;
    return off;
}
static bool capPositionsDisabled() {
    static const bool off = std::getenv("MVB_NO_CAP_POSITIONS") != nullptr ||
                            std::getenv("MVB_NO_SECTION_HANDOFF") != nullptr;
    return off;
}
// PHASE CONTINUATION IS OFF BY DEFAULT. The 2026-08-30 control matrix separated the handoff's
// two halves and they behave oppositely: adopting a neighbour's vertex POSITIONS is what makes
// 12_boost mesh (262/262; without it 1/262), while ROTATING a piece's own section to a
// neighbour's phase is what breaks 04_forward (5/418 with phase continuation in ANY gating,
// 418/418 with the whole handoff off). A position adoption is exact -- the two pieces then own
// literally the same cap polygon -- whereas a phase rotation re-generates a DIFFERENT polygon
// that only agrees with the neighbour at one vertex when the junction is oblique. So: adopt
// positions where they fit, and never fall back to a rotation. MVB_PHASE_CONTINUATION=1
// restores it for A/B measurement.
static bool phaseContinuationEnabled() {
    // ON by default: this is part of the configuration measured green on 10_emi (45/45),
    // 12_boost (262/262) and 02_flyback (588/588). The 2026-08-30 attempt to default it off
    // ("positions only") was validated on 04 only -- 04 stayed red regardless -- and 10's
    // riser junctions were never re-tested under it. MVB_NO_PHASE_CONTINUATION=1 for A/B.
    static const bool on = std::getenv("MVB_NO_PHASE_CONTINUATION") == nullptr &&
                           std::getenv("MVB_NO_SECTION_HANDOFF") == nullptr;
    return on;
}

static bool sectionPhaseOf(const TopoDS_Wire& w, const gp_Pnt& center, const gp_Dir& normal,
                           double r, int segments, double& phase0) {
    if (handoffDisabled() || !phaseContinuationEnabled()) return false;
    int n = 0;
    gp_Pnt v0;
    bool got = false;
    for (BRepTools_WireExplorer we(w); we.More(); we.Next()) {
        if (!got) { v0 = BRep_Tool::Pnt(we.CurrentVertex()); got = true; }
        ++n;
    }
    if (!got || n != segments) return false;
    const gp_Ax2 fr = deterministicSectionFrame(center, normal);
    const gp_XYZ d = v0.XYZ() - center.XYZ();
    const double x = d.Dot(fr.XDirection().XYZ());
    const double y = d.Dot(fr.YDirection().XYZ());
    if (x * x + y * y < 0.25 * r * r) return false;
    phase0 = std::atan2(y, x) - kPi / segments;
    return true;
}

// The same round profile, but as TWO half-circle edges instead of one closed circle.
// Sweeping a CLOSED circular edge produces a single PERIODIC B-spline surface carrying a seam
// curve, and gmsh cannot order the boundary loop of such a face -- it aborts the 2D stage with
// "the 1D mesh seems not to be forming a closed loop" (the offending face reported simply as
// "BSpline surface"). That is what made the real-winding conductor unmeshable: every other
// primitive is an analytic cylinder or torus, which mesh fine; only the swept axial transitions
// carried a seam. Split in two and the sweep yields ordinary patches with real boundaries -- the
// identical solid, now meshable. (ABT #332)
TopoDS_Wire wireProfileWireSplit(const gp_Pnt& center, const gp_Dir& normal, double radius,
                                 int pieces = 2) {
    // The SAME deterministic frame as wireProfileWire. This used to be a bare
    // gp_Ax2(center, normal), so a SEG built with the split profile and a tangent neighbour
    // built with the polygon/circle profile started their sections at different angles: the
    // abutting faces then coincide only partially, which is precisely the "co-located but not
    // bit-coincident" junction that survives OCCBooleanGlue and reaches gmsh as overlapping
    // facets (RC2 of the 2026-08-24 corpus audit).
    gp_Circ circ(deterministicSectionFrame(center, normal), radius);
    BRepBuilderAPI_MakeWire w;
    for (int i = 0; i < pieces; ++i)
        w.Add(BRepBuilderAPI_MakeEdge(circ, kTwoPi * i / pieces, kTwoPi * (i + 1) / pieces).Edge());
    return w.Wire();
}

// Oriented RECTANGLE profile for rectangular/planar wire, in the plane perpendicular to the wire's
// travel (tangent). The AXIAL in-plane axis is the column axis projected perpendicular to the
// tangent; the RADIAL axis is tangent x axial. width spans radial, height spans axial -- matching
// TurnBuilder::build_rect_profile so the real-winding section is consistent with the per-turn one.
TopoDS_Wire rectProfileWire(const gp_Pnt& center, const gp_Dir& tangent, const gp_Dir& axialAxis,
                            double width, double height) {
    gp_XYZ t = tangent.XYZ();
    gp_XYZ a = axialAxis.XYZ() - t * (axialAxis.XYZ().Dot(t));  // axial, ⊥ tangent
    if (a.Modulus() < 1e-9) {                                   // tangent ∥ axial: pick any ⊥ axis
        gp_XYZ ref = std::abs(t.X()) < 0.9 ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 0, 1);
        a = ref - t * (ref.Dot(t));
    }
    a.Normalize();
    gp_XYZ w = t.Crossed(a);                                    // radial
    w.Normalize();
    double hw = width / 2.0, hh = height / 2.0;
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(center.XYZ() - w * hw - a * hh));
    poly.Add(gp_Pnt(center.XYZ() + w * hw - a * hh));
    poly.Add(gp_Pnt(center.XYZ() + w * hw + a * hh));
    poly.Add(gp_Pnt(center.XYZ() - w * hw + a * hh));
    poly.Close();
    return poly.Wire();
}

// One analytic edge per primitive, so the whole conductor is ONE path:
//  - SEG: straight edge.
//  - ARC3: exact circle edge trimmed to [0, sweep].
//  - SPIRAL at constant radius: TRUE HELIX — a 2D line in the (U = azimuth, V = height)
//    parametric space of a cylinder about its vertical axis (the canonical OCCT
//    construction: Geom2d line on Geom_CylindricalSurface + BRepLib::BuildCurves3d).
//  - SPIRAL with varying radius (conical wraps, oblong cap ramps): BSpline through the
//    sampled centreline.
// A spiral drawn on an analytic surface can absorb its mitre growth INTO THE PCURVE, so the
// grown spine stays ONE edge on ONE surface. ABT #685 (Alf, 2026-08-18): the alternative --
// splicing a straight tangent edge onto each end, which is what every other primitive does --
// makes MakePipeShell emit a lateral face PER SPINE EDGE, and the mitre plane then falls exactly
// on the seam between two of them (the plane passes through the primitive's own endpoint, which
// is precisely where the splice sits). OCC's boolean cannot cut a pipe on a face boundary: every
// failure measured had faces=5 (3 lateral + 2 caps) and returned debris -- 06_llc lost 3.0 of
// 3.87 mm3 to a 0.59 mm3 knife. With the growth inside the pcurve there is no seam to cut on.
bool primEdgeGrowsAnalytically(const Primitive& pr) {
    if (pr.kind != Primitive::SPIRAL) return false;
    const Spiral& sp = pr.spiral;
    return std::abs(sp.r1 - sp.r0) < 1e-12 || std::abs(sp.y1 - sp.y0) > 1e-12;
}

TopoDS_Edge primEdge(const Primitive& pr, double wireRadius, double overA, double overB) {
    if (pr.kind == Primitive::SEG) {
        if (pr.seg.a.Distance(pr.seg.b) < 1e-12) return TopoDS_Edge();
        return BRepBuilderAPI_MakeEdge(pr.seg.a, pr.seg.b).Edge();
    }
    if (pr.kind == Primitive::ARC3) {
        double radius = pr.arc.v0.Modulus();
        if (radius < 1e-12 || std::abs(pr.arc.sweep) < 1e-12) return TopoDS_Edge();
        try {
            // gp_Circ parametrisation P(u) = c + r (cos u XDir + sin u YDir) with
            // YDir = axis x XDir matches rotateXYZ (axis is perpendicular to v0).
            gp_Ax2 ax2(pr.arc.c, gp_Dir(pr.arc.axis), gp_Dir(pr.arc.v0));
            gp_Circ circ(ax2, radius);
            return BRepBuilderAPI_MakeEdge(circ, 0.0, pr.arc.sweep).Edge();
        } catch (const Standard_Failure&) {
            return TopoDS_Edge();
        }
    }
    // THE SPINE'S 3D CURVE IS AN APPROXIMATION OF THE PCURVE, AND ITS TOLERANCE IS THE
    // COPPER'S (2026-09-02, flyback_transformer_complete ride-over). BRepLib::BuildCurves3d's
    // default tolerance is 1e-5 in the edge's own unit -- these edges are in METRES, so the
    // helix reaches MakePipeShell as a 5-pole B-spline up to 10 um off the true helix.
    // Measured on the flyback's 'Secondary parallel 0' raised wrap: spine 2.39 um off, swept
    // tube radius about the exact helix 0.5553..0.5601 mm for a 0.55771 mm wire. The
    // centreline model is exact (the ride-over sits 1.115499990 mm from the dragback run, one
    // coated OD to the bit), and the bare-to-bare gap that leaves is 80 nm -- so the
    // approximation alone puts the swept wrap 1.50 um INSIDE the exact cylinder of the run
    // (omfem_step_intersect: OVERLAP on 8 ride-over pairs; every helix-vs-straight pair in the
    // corpus is exposed to it, helix-vs-helix neighbours only hide it because they share the
    // error). MakePipeShell's own surface fit (mm frame, 1e-4 mm) is within 1 nm and is not
    // the deficit. 1e-8 m = 10 nm: the linear helix and cone fit to <= 1 nm with 8 and 6 poles
    // (5 and 4 before), the cosine blends to <= 7 nm with ~50 poles, and the sweep time is
    // unchanged (measured standalone, scratchpad/ride/helix_sweep_exp*.log).
    constexpr double kSpine3dTol = 1e-8;   // metres
    if (pr.kind == Primitive::SPIRAL) {
        // Every spiral lives on an ANALYTIC surface of revolution about its vertical
        // axis: a cylinder when the radius is constant, a cone otherwise (radius and
        // height share the blend factor, so the (r, y) trace is a straight meridian).
        // Building the edge as a 2D curve in the surface's (U = azimuth, V) parameter
        // space (+ BRepLib::BuildCurves3d) is the one construction OCCT's pipe-shell
        // sweeps reliably — free-form interpolated 3D splines, though positionally
        // tight, produce unbounded surface flares when swept.
        const Spiral& sp = pr.spiral;
        double dr = sp.r1 - sp.r0;
        double dy = sp.y1 - sp.y0;
        bool constantRadius = std::abs(dr) < 1e-12;
        if (constantRadius || std::abs(dy) > 1e-12) {
            try {
                // Frame with U measured from +X toward -Z — exactly the azPointC()
                // convention: P(U,V) matches (cx + r cos U, y, cz - r sin U).
                gp_Ax3 axis(gp_Pnt(sp.cx, sp.y0, sp.cz), gp_Dir(0, 1, 0),
                            gp_Dir(1, 0, 0));
                Handle(Geom_Surface) surf;
                double v1;   // V at the end (V = 0 at the start)
                if (constantRadius) {
                    axis.SetLocation(gp_Pnt(sp.cx, 0, sp.cz));
                    surf = new Geom_CylindricalSurface(axis, sp.r0);
                    // Cylinder V is plain height.
                    v1 = dy;
                } else {
                    double alpha = std::atan(dr / dy);
                    surf = new Geom_ConicalSurface(axis, alpha, sp.r0);
                    v1 = dy / std::cos(alpha);
                }
                double v0 = constantRadius ? sp.y0 : 0.0;
                // Mitre growth, in arc length, converted to azimuth on this surface: a step du
                // advances sqrt(radius^2 + (dv/du)^2) of arc. Extending the pcurve's range is an
                // EXACT continuation of the same helix/conical spiral for the linear case, and a
                // tangent (constant-V) continuation for the blend, whose end tangents are purely
                // azimuthal by construction. Either way it stays one edge -- see the note above.
                const double dAz = sp.az1 - sp.az0;
                const double sgn = (dAz >= 0.0) ? 1.0 : -1.0;
                const double dvdu = (std::abs(dAz) > 1e-12) ? v1 / dAz : 0.0;
                const double dsduA = std::sqrt(sp.r0 * sp.r0 + dvdu * dvdu);
                const double dsduB = std::sqrt(sp.r1 * sp.r1 + dvdu * dvdu);
                const double az0e =
                    sp.az0 - (overA > 0.0 && dsduA > 1e-12 ? sgn * overA / dsduA : 0.0);
                const double az1e =
                    sp.az1 + (overB > 0.0 && dsduB > 1e-12 ? sgn * overB / dsduB : 0.0);
                TopoDS_Edge e;
                if (!sp.blend) {
                    auto vAt = [&](double az) { return v0 + (az - sp.az0) * dvdu; };
                    Handle(Geom2d_TrimmedCurve) seg2d =
                        GCE2d_MakeSegment(gp_Pnt2d(az0e, vAt(az0e)),
                                          gp_Pnt2d(az1e, vAt(az1e))).Value();
                    e = BRepBuilderAPI_MakeEdge(seg2d, surf).Edge();
                } else {
                    // Cosine blend of V over U — end tangents purely azimuthal, so the
                    // junctions into the neighbouring on-station geometry are exact.
                    int n = curveSampleCount(std::max(sp.r0, sp.r1),
                                             std::abs(az1e - az0e), wireRadius);
                    Handle(TColgp_HArray1OfPnt2d) arr =
                        new TColgp_HArray1OfPnt2d(1, n);
                    for (int i = 0; i < n; ++i) {
                        double az = az0e + (az1e - az0e) * (static_cast<double>(i) / (n - 1));
                        // Clamped: outside the primitive's own range the blend holds its end
                        // value, which IS the tangent continuation (the cosine's slope is zero
                        // at both ends), so the imposed (1,0) end tangents stay exact.
                        const double t = std::clamp(
                            std::abs(dAz) > 1e-12 ? (az - sp.az0) / dAz : 0.0, 0.0, 1.0);
                        double f = 0.5 * (1.0 - std::cos(kPi * t));
                        arr->SetValue(i + 1, gp_Pnt2d(az, v0 + v1 * f));
                    }
                    Geom2dAPI_Interpolate interp(arr, Standard_False, 1e-12);
                    interp.Load(gp_Vec2d(1, 0), gp_Vec2d(1, 0), Standard_True);
                    interp.Perform();
                    if (!interp.IsDone() || interp.Curve().IsNull())
                        return TopoDS_Edge();
                    e = BRepBuilderAPI_MakeEdge(interp.Curve(), surf).Edge();
                }
                BRepLib::BuildCurves3d(e, kSpine3dTol);
                return e;
            } catch (const Standard_Failure&) {
                return TopoDS_Edge();
            }
        }
        // FLAT SPIRAL (radius varies at exactly constant height) -- blended OR linear.
        // A 2D curve in the horizontal plane's parameter space, interpolated with the EXACT
        // analytic end tangents.
        //
        // ABT #373 (:818): the LINEAR case used to fall through to a 3D BSpline through
        // samplePrim's points, and on a toroid's face crossing samplePrim returned TWO of them
        // (see spiralSampleCount) -- so the emitted piece was a straight chord swept by
        // MakePipeShell along a straight spine, which carries the START tangent's section to
        // BOTH ends. The assembler, meanwhile, mitres against spiralTangent's analytic END
        // tangent. On realwinding_toroid's 'turn 6' -> 'turn 7' bottom chord those differ by the
        // spiral's whole turning, 4.21 deg: the joint was read as tangent (3.6e-11 deg, bridged,
        // uncut) and the chord and the bottom-inner-corner arc interpenetrated by 0.00512 mm^3.
        // Interpolating in the plane with the analytic tangents IMPOSED makes the swept end
        // tangent equal the one the junction logic uses, by construction -- which is what
        // "neighbours meet on the plane bisecting their tangents" requires.
        {
            const double dAz = sp.az1 - sp.az0;
            // No azimuthal advance at all: the locus IS the straight radial segment between the
            // endpoints (r linear, y and az constant), and that is the exact edge for it --
            // faceSpiralTangent reports the same purely radial direction for this case.
            if (std::abs(dAz) < 1e-12) {
                const auto ends = primEndpoints(pr);
                if (ends.first.Distance(ends.second) < 1e-12) return TopoDS_Edge();
                return BRepBuilderAPI_MakeEdge(ends.first, ends.second).Edge();
            }
            try {
                gp_Ax3 frame(gp_Pnt(sp.cx, sp.y0, sp.cz), gp_Dir(0, 1, 0),
                             gp_Dir(1, 0, 0));
                Handle(Geom_Plane) plane = new Geom_Plane(frame);
                // Plane P(U,V) = O + U XDir + V YDir with YDir = -Z, so the azPointC
                // trace maps to (U, V) = (r cos az, r sin az).
                int n = spiralSampleCount(sp, wireRadius);
                Handle(TColgp_HArray1OfPnt2d) arr = new TColgp_HArray1OfPnt2d(1, n);
                for (int i = 0; i < n; ++i) {
                    double t = static_cast<double>(i) / (n - 1);
                    double f = sp.blend ? 0.5 * (1.0 - std::cos(kPi * t)) : t;
                    double az = sp.az0 + dAz * t;
                    double r = sp.r0 + dr * f;
                    arr->SetValue(i + 1, gp_Pnt2d(r * std::cos(az), r * std::sin(az)));
                }
                // d/dt (r cos az, r sin az) with az' = dAz and r' = dr f'(t). A cosine blend has
                // f'(0) = f'(1) = 0, so this reduces to dAz * r * (-sin az, cos az) -- the purely
                // azimuthal tangent the blend has always been given, now carrying dAz's SIGN as
                // well (a spiral running the other way round was previously handed a reversed
                // end tangent).
                auto tangentAt = [&](double az, double r, double rp) {
                    return gp_Vec2d(rp * std::cos(az) - r * dAz * std::sin(az),
                                    rp * std::sin(az) + r * dAz * std::cos(az));
                };
                const double rp = sp.blend ? 0.0 : dr;
                gp_Vec2d ta = tangentAt(sp.az0, sp.r0, rp);
                gp_Vec2d tb = tangentAt(sp.az1, sp.r1, rp);
                Geom2dAPI_Interpolate interp(arr, Standard_False, 1e-12);
                interp.Load(ta, tb, Standard_True);
                interp.Perform();
                if (!interp.IsDone() || interp.Curve().IsNull()) return TopoDS_Edge();
                TopoDS_Edge e = BRepBuilderAPI_MakeEdge(interp.Curve(), plane).Edge();
                BRepLib::BuildCurves3d(e, kSpine3dTol);
                return e;
            } catch (const Standard_Failure&) {
                return TopoDS_Edge();
            }
        }
    }
    if (pr.kind == Primitive::BLEND) {
        // The S-blend is PLANAR (the longitudinal direction and the offset span one
        // plane), so it too builds as a 2D cosine in an analytic surface's parameter
        // space — the pipe-reliable construction (see the spiral comment above).
        const Blend& bl = pr.blendc;
        gp_XYZ d = bl.b.XYZ() - bl.a.XYZ();
        double L = d.Dot(bl.u);
        gp_XYZ perp = d - bl.u * L;
        double perpM = perp.Modulus();
        if (L < 1e-12) return TopoDS_Edge();
        if (perpM < 1e-12) return BRepBuilderAPI_MakeEdge(bl.a, bl.b).Edge();
        try {
            gp_Ax3 frame(bl.a, gp_Dir(bl.u.Crossed(perp)), gp_Dir(bl.u));
            Handle(Geom_Plane) plane = new Geom_Plane(frame);
            // Plane P(U,V) = O + U*XDir + V*YDir with XDir = u, YDir = axis x XDir =
            // unit(perp).
            double stepL = 2.0 * std::sqrt(std::max(
                1e-18, 4.0 * L * L / (kPi * kPi * perpM) *
                           samplingSag(wireRadius)));
            int n = std::clamp(static_cast<int>(std::ceil(L / stepL)) + 1, 8, 512);
            Handle(TColgp_HArray1OfPnt2d) arr = new TColgp_HArray1OfPnt2d(1, n);
            for (int i = 0; i < n; ++i) {
                double t = static_cast<double>(i) / (n - 1);
                arr->SetValue(i + 1, gp_Pnt2d(L * t,
                                              perpM * 0.5 * (1.0 - std::cos(kPi * t))));
            }
            Geom2dAPI_Interpolate interp(arr, Standard_False, 1e-12);
            interp.Load(gp_Vec2d(1, 0), gp_Vec2d(1, 0), Standard_True);
            interp.Perform();
            if (!interp.IsDone() || interp.Curve().IsNull()) return TopoDS_Edge();
            TopoDS_Edge e = BRepBuilderAPI_MakeEdge(interp.Curve(), plane).Edge();
            BRepLib::BuildCurves3d(e, kSpine3dTol);
            return e;
        } catch (const Standard_Failure&) {
            return TopoDS_Edge();
        }
    }
    // Every kind returns above: SEG, ARC3, and all three spiral families (on a cylinder, on a
    // cone, and flat in a plane), plus BLEND. ABT #373: the flat LINEAR spiral used to land here
    // on a 3D BSpline through samplePrim's points -- the construction whose two-point degenerate
    // case produced the :818 mitre interpenetration -- and it now has its own analytic branch.
    // Nothing may fall through to a silent empty edge; a new Primitive kind must bring its own
    // construction.
    throw std::runtime_error(
        "WireAssembler: primEdge has no construction for primitive '" + pr.label +
        "' (unknown kind " + std::to_string(static_cast<int>(pr.kind)) + ")");
}

// Fuse all solids of a conductor's compound into as FEW bodies as possible, in ONE general
// boolean (BRepAlgoAPI_Fuse over the whole set via SetArguments/SetTools). A single BOP with
// OBB culling is both far faster than N sequential pairwise fuses (which is O(N) booleans on
// swept BSpline solids — minutes for a many-piece winding) AND more robust: OCCT resolves all
// the mutual intersections together, so overlapping chunk pipes + elbow/junction spheres merge
// into one solid while any piece that genuinely touches nothing stays its own solid. No glue:
// the pieces OVERLAP volumetrically (junction spheres straddle the pipe caps), which is a real
// intersection, not the coincident-face case glue is for. A gross-corruption volume guard
// rejects the known OCCT "ate a whole body" defect; on any failure the exact (unfused) compound
// is returned rather than losing copper — same guarantee as before, just reached far less often.
// A conductor face whose bounding extent is many wire-radii long but whose area is a sliver:
// the degenerate SEAM faces OCCT's boolean union leaves where two swept pipes (or a pipe and
// its junction sphere) merge -- the thin "line/cone" sheets seen on a multilayer winding. A
// genuine tube-wall face of the same extent has area ~ 2*pi*r*length, orders of magnitude
// larger; end caps and short faces (diag <= a few radii) are exempt.
bool hasDegenerateSheetFace(const TopoDS_Shape& shape, double wireRadius) {
    for (TopExp_Explorer fe(shape, TopAbs_FACE); fe.More(); fe.Next()) {
        Bnd_Box fb;
        BRepBndLib::Add(fe.Current(), fb);
        if (fb.IsVoid()) continue;
        double x0, y0, z0, x1, y1, z1;
        fb.Get(x0, y0, z0, x1, y1, z1);
        double diag = std::sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) +
                                (z1 - z0) * (z1 - z0));
        if (diag <= 4.0 * wireRadius) continue;
        GProp_GProps fp;
        BRepGProp::SurfaceProperties(fe.Current(), fp);
        if (fp.Mass() < 0.2 * wireRadius * diag) return true;
    }
    return false;
}

// Drop degenerate (near-zero-volume) solids from a conductor result. OCCT's fuse can leave
// thin sheet/shell solids behind when it welds a dense multilayer winding: they carry no
// copper, render as stray sheets, and break tetrahedral meshing. Any real conductor piece --
// even a lead-junction sphere -- is orders of magnitude above the 1e-12 m^3 (1e-3 mm^3) floor.
TopoDS_Shape pruneDegenerateSolids(const TopoDS_Shape& shape, double wireRadius) {
    // The threshold is RELATIVE TO THE WIRE, never an absolute volume. This exists to drop
    // boolean SLIVERS -- sheets with no thickness -- so the scale that matters is the wire's
    // own: a piece thinner than a hundredth of a wire radius along the wire is a sliver, any
    // longer piece is real copper. An absolute 1e-12 m^3 cutoff silently deleted every bump
    // riser of 13_current_sense's 0.049 mm-radius secondary (riser volume 8.4e-13 m^3, a
    // legitimate 0.112 mm step), leaving the bump arcs disconnected -- Alf, 2026-08-08:
    // "Secondary parallel 0 16 and 17 are not connected, they are missing the straight
    // segment of the bump".
    const double kMinSolidVolume =
        kPi * wireRadius * wireRadius * (wireRadius / 100.0);
    std::vector<TopoDS_Shape> kept;
    int total = 0;
    double droppedVolume = 0.0;
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        ++total;
        GProp_GProps props;
        BRepGProp::VolumeProperties(exp.Current(), props);
        if (std::abs(props.Mass()) >= kMinSolidVolume) kept.push_back(exp.Current());
        else droppedVolume += std::abs(props.Mass());
    }
    // Dropping copper is never routine: say so.
    if (static_cast<int>(kept.size()) != total)
        std::cerr << "[ConductorBuilder] pruneDegenerateSolids dropped "
                  << (total - static_cast<int>(kept.size())) << " of " << total
                  << " solids as slivers (" << droppedVolume * 1e9 << " mm3, threshold "
                  << kMinSolidVolume * 1e9 << " mm3)\n";
    if (kept.empty() || static_cast<int>(kept.size()) == total) return shape;
    if (kept.size() == 1) return kept.front();
    TopoDS_Compound out;
    BRep_Builder b;
    b.MakeCompound(out);
    for (const auto& s : kept) b.Add(out, s);
    return out;
}

// ---------------------------------------------------------------------------------------
// Conformal per-primitive toroid (mitre joints). A dense toroid's winding cannot become one
// solid (MakePipe self-intersects on the packed hole spine; a boolean fuse of the per-run pipes
// comes out invalid), and the plain per-primitive compound OVERLAPS at every joint (double
// material -> not a valid FEM assembly). Here each primitive is built as its OWN round solid,
// grown a hair past both endpoints, then sliced by the ANGLE-BISECTOR plane through each shared
// point. Two equal round tubes mitred on their bisector leave IDENTICAL elliptical faces (same
// centre, tilt and radius), so neighbouring solids ABUT on a coincident face instead of
// overlapping: a conformal, mesh-shareable assembly. Centrelines/crossings never move.
// ANALYTIC end tangents, not sampled chords. Chord directions from samplePrim run up to half
// a sample step's turn away from the true tangent (a 16-sample wrap chord is ~5-6 deg off),
// which made emitToroidConformal misread every TANGENT wrap junction as a corner and slice it
// on a chord-tilted bisector: the mitre faces then no longer coincide and the 'conformal'
// assembly falls apart into disconnected pieces (measured: PQ33 round column, 18 prims,
// 34 spurious cuts, 14 components).
static gp_Dir spiralTangent(const Spiral& sp, bool atStart) {
    const double daz = sp.az1 - sp.az0;
    const double az = atStart ? sp.az0 : sp.az1;
    double rp = 0.0, yp = 0.0;
    const double r = atStart ? sp.r0 : sp.r1;
    if (std::fabs(daz) > 1e-12 && !sp.blend) {
        rp = (sp.r1 - sp.r0) / daz;
        yp = (sp.y1 - sp.y0) / daz;
    }   // cosine-blend spirals have purely azimuthal end tangents (r' = y' = 0)
    // P(az) = (cx + r cos az, y, cz - r sin az)  [azPointC convention]
    gp_XYZ t(rp * std::cos(az) - r * std::sin(az), yp,
             -rp * std::sin(az) - r * std::cos(az));
    if (daz < 0) t *= -1.0;
    if (t.Modulus() < 1e-15) return gp_Dir(1, 0, 0);
    return gp_Dir(t);
}

gp_Dir primFwdStart(const Primitive& p, double r) {
    switch (p.kind) {
        case Primitive::SEG: {
            gp_Vec v(p.seg.a, p.seg.b);
            if (v.Magnitude() > 1e-12) return gp_Dir(v);
            break;
        }
        case Primitive::ARC3: {
            gp_XYZ t = p.arc.axis.Crossed(p.arc.v0);
            if (p.arc.sweep < 0) t *= -1.0;
            if (t.Modulus() > 1e-15) return gp_Dir(t);
            break;
        }
        case Primitive::SPIRAL:
            return spiralTangent(p.spiral, /*atStart=*/true);
        case Primitive::BLEND:
            if (p.blendc.u.Modulus() > 1e-15) return gp_Dir(p.blendc.u);
            break;
    }
    auto pts = samplePrim(p, r);
    for (size_t i = 1; i < pts.size(); ++i) {
        gp_Vec v(pts.front(), pts[i]);
        if (v.Magnitude() > 1e-12) return gp_Dir(v);
    }
    return gp_Dir(1, 0, 0);
}

gp_Dir primFwdEnd(const Primitive& p, double r) {
    switch (p.kind) {
        case Primitive::SEG: {
            gp_Vec v(p.seg.a, p.seg.b);
            if (v.Magnitude() > 1e-12) return gp_Dir(v);
            break;
        }
        case Primitive::ARC3: {
            gp_XYZ ve = rotateXYZ(p.arc.v0, p.arc.axis, p.arc.sweep);
            gp_XYZ t = p.arc.axis.Crossed(ve);
            if (p.arc.sweep < 0) t *= -1.0;
            if (t.Modulus() > 1e-15) return gp_Dir(t);
            break;
        }
        case Primitive::SPIRAL:
            return spiralTangent(p.spiral, /*atStart=*/false);
        case Primitive::BLEND:
            if (p.blendc.u.Modulus() > 1e-15) return gp_Dir(p.blendc.u);
            break;
    }
    auto pts = samplePrim(p, r);
    for (size_t i = pts.size(); i-- > 1;) {
        gp_Vec v(pts[i - 1], pts.back());
        if (v.Magnitude() > 1e-12) return gp_Dir(v);
    }
    return gp_Dir(1, 0, 0);
}


// ---- CHUNK 7: THE JOINT SPHERE (Alf, 2026-08-25) -----------------------------------------
// A wire IS the set of points within r of its centreline -- nothing else is copper. Two tubes
// meeting at a corner cover all of that except the notch on the OUTSIDE of the bend, and the
// exact shape filling that notch is a SPHERE of radius r at the corner point. That is why this
// is the physically right fill and the axial "bridge growth" it replaces was not: growing a tube
// by r*tan(theta) past the joint puts copper up to r*tan(theta)*sin(theta) OUTSIDE the envelope,
// copper the real wire does not have. Measured on 14_dab: an ~8.6 deg bridged junction poked out
// 9.67 um and interpenetrated the sibling parallel by 9.6 um, while MKF's stations were correctly
// spaced (0.855247 mm against a 0.855 mm insulated envelope) and the analytic centrelines cleared.
// A sphere of radius r is within r of the centreline everywhere, so it CANNOT reach a foreign
// conductor whose centreline is rA + rB away. No margin, no tolerance, no chosen constant.
//
// FACETED, and INSCRIBED like wireProfileWire, for two reasons: the drawn copper must never
// protrude past the real envelope, and a quadric sphere is closed in BOTH parameters, so it stays
// periodic even after the NURBS pass and would land straight back on gmsh's fragile periodic
// mesher (OCCFace.cpp:142). Planar facets are never periodic.
static TopoDS_Shape jointSphere(const gp_Pnt& c, double r, int segments) {
    if (r <= 0.0) return {};
    if (segments <= 0) return BRepPrimAPI_MakeSphere(c, r).Shape();
    const int nu = std::max(4, segments);          // divisions around the equator
    const int nv = std::max(2, segments / 2);      // divisions pole to pole
    const gp_Pnt north(c.X(), c.Y(), c.Z() + r);
    const gp_Pnt south(c.X(), c.Y(), c.Z() - r);
    std::vector<std::vector<gp_Pnt>> ring(nv);     // latitudes j = 1 .. nv-1 (index j-1)
    for (int j = 1; j < nv; ++j) {
        const double phi = kPi * j / nv;
        const double z = r * std::cos(phi), rr = r * std::sin(phi);
        std::vector<gp_Pnt> row;
        row.reserve(nu);
        for (int i = 0; i < nu; ++i) {
            const double th = kTwoPi * i / nu;
            row.emplace_back(c.X() + rr * std::cos(th), c.Y() + rr * std::sin(th), c.Z() + z);
        }
        ring[j - 1] = std::move(row);
    }
    BRepBuilderAPI_Sewing sew(1e-9);
    auto addFace = [&](std::initializer_list<gp_Pnt> pts) {
        BRepBuilderAPI_MakePolygon poly;
        for (const auto& p : pts) poly.Add(p);
        poly.Close();
        if (!poly.IsDone()) return;
        BRepBuilderAPI_MakeFace f(poly.Wire());
        if (f.IsDone()) sew.Add(f.Face());
    };
    for (int i = 0; i < nu; ++i) {
        const int i2 = (i + 1) % nu;
        addFace({north, ring[0][i], ring[0][i2]});
        addFace({south, ring[nv - 2][i2], ring[nv - 2][i]});
    }
    for (int j = 0; j + 1 < nv - 1; ++j)
        for (int i = 0; i < nu; ++i) {
            const int i2 = (i + 1) % nu;
            addFace({ring[j][i], ring[j][i2], ring[j + 1][i2], ring[j + 1][i]});
        }
    sew.Perform();
    const TopoDS_Shape shell = sew.SewedShape();
    if (shell.IsNull()) return {};
    for (TopExp_Explorer ex(shell, TopAbs_SHELL); ex.More(); ex.Next()) {
        BRepBuilderAPI_MakeSolid mk(TopoDS::Shell(ex.Current()));
        if (mk.IsDone() && !mk.Shape().IsNull()) return mk.Shape();
    }
    return {};
}

// A round solid for one primitive, grown by overA/overB past its start/end so a tilted bisector
// plane fully crosses the tube there (SEG -> longer cylinder, ARC3 -> wider revolve, else pipe).
// A tangent junction passes 0 -> the tube ends flush on a perpendicular cap, no growth, no cut.
// segments <= 0 keeps the EXACT round section (analytic cylinders/tori); segments > 0 emits the
// n-gon prism/revolve. This is not cosmetic: an exact round sweep produces PERIODIC surfaces
// (seam-carrying cylinders and tori), and gmsh's "Impossible to mesh periodic surface" is a
// known failure cluster on them -- the faceted section has no seam at all, which is why the
// pipeline has always meshed with segments=12. The conformal assembler must honour the same
// knob the retired sweep path did.
// startProfileOverride / endCapOut: SECTION-WIRE HANDOFF ACROSS TANGENT JUNCTIONS. At
// segments > 0 a prism's polygon phase comes from deterministicSectionFrame while a pipe's
// end phase is whatever MakePipeShell's transport produced -- two 16-gons in the same plane
// with different phases, which is the measured source of the ride-chain junction lenses
// (7.2e-5 mm^3) and of the imprint shards their trims leave (ABT #935: a 30 um grazing band
// whose 2D mesh folds). Handing the ACTUAL end-section wire to the next piece as its start
// profile makes every tangent junction bit-identical by construction. The override is only
// adopted after a geometric sanity check (vertices in the start plane, at the wire radius);
// anything else falls back to the piece's own profile.
static bool sectionWireFits(const TopoDS_Wire& w, const gp_Pnt& a, const gp_Dir& dir, double r,
                            int segments) {
    int n = 0;
    for (TopExp_Explorer vx(w, TopAbs_VERTEX); vx.More(); vx.Next()) ++n;
    if (segments > 0 && n != 2 * segments) return false;  // closed wire: each vertex twice
    // Same plane tolerance as the prism-side adoption (0.25 r): a neighbour's real end
    // section sits where its sweep put it, measured 0.13 r off the nominal plane on 10_emi,
    // while genuine endpoint mismatches are 1.5 r. Refusing the real cap left an
    // audit-visible 3.9e-5 mm^3 lens at exactly the riser->pipe junctions.
    const double planeTol = std::max(2e-6, 0.25 * r);
    for (TopExp_Explorer vx(w, TopAbs_VERTEX); vx.More(); vx.Next()) {
        const gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(vx.Current()));
        const gp_Vec d(a, p);
        if (std::abs(d.Dot(gp_Vec(dir))) > planeTol) return false;  // out of the start plane
        const double rad = d.Magnitude();
        if (rad < 0.8 * r || rad > 1.6 * r) return false;           // not this wire's section
    }
    return true;
}

static TopoDS_Shape rawGrownSolid(const Primitive& pr, double r, double overA, double overB,
                                  int segments,
                                  const TopoDS_Wire* startProfileOverride = nullptr,
                                  TopoDS_Wire* endCapOut = nullptr,
                                  bool* startPhasedOut = nullptr) {
    if (pr.kind == Primitive::SEG) {
        gp_Vec d(pr.seg.a, pr.seg.b);
        double len = d.Magnitude();
        if (len < 1e-12) return {};
        gp_Dir dir(d);
        gp_Pnt a = pr.seg.a.Translated(gp_Vec(dir) * (-overA));
        const double total = len + overA + overB;
        // Never a quadric MakeCylinder: its side is ONE closed periodic face, which gmsh sends
        // to its fragile periodic mesher (OCCFace.cpp:142). A prism over the SPLIT profile is the
        // identical solid built from ordinary patches -- the ABT #332 remedy, applied here too.
        TopoDS_Wire segProfWire;
        double segPhase = 0.0;
        if (startProfileOverride && !startProfileOverride->IsNull() && segments > 0 &&
            overA <= 0.0 && sectionPhaseOf(*startProfileOverride, a, dir, r, segments, segPhase)) {
            segProfWire = wireProfileWirePhased(a, dir, r, segments, segPhase);
            if (startPhasedOut) *startPhasedOut = true;
        } else if (segments > 0) {
            segProfWire = wireProfileWire(a, dir, r, segments);
        }
        TopoDS_Face prof =
            !segProfWire.IsNull()
                ? BRepBuilderAPI_MakeFace(segProfWire).Face()
                : BRepBuilderAPI_MakeFace(wireProfileWireSplit(a, dir, r, 2)).Face();
        if (endCapOut && !segProfWire.IsNull()) {
            gp_Trsf tr;
            tr.SetTranslation(gp_Vec(dir) * total);
            *endCapOut =
                TopoDS::Wire(BRepBuilderAPI_Transform(segProfWire, tr, Standard_True).Shape());
        }
        return BRepPrimAPI_MakePrism(prof, gp_Vec(dir) * total).Shape();
    }
    if (pr.kind == Primitive::ARC3) {
        double radius = pr.arc.v0.Modulus();
        if (radius < 1e-12 || std::abs(pr.arc.sweep) < 1e-12) return {};
        // TIGHT-ARC REPORT (diagnostic only -- changes no geometry). This ARC3 is built by
        // REVOLVING the wire profile about an axis at distance `radius`. That revolve
        // self-intersects when radius <= r, because the profile then sweeps THROUGH the axis;
        // it is marginal for radius only slightly above r, where the inner facets converge on
        // the axis. Nothing here checked it, and the result is copper that is BRepCheck-VALID,
        // watertight and overlap-clean while being BOPAlgo SELF-INTERSECTING (ABT #958:
        // 02_flyback ships 42 such chunks, exactly one per turn of Primary parallel 0).
        // MVB_TIGHTARC_DIAG=<ratio> widens the report to every arc below that bend ratio
        // (any non-numeric value keeps the historical 2.0).
        const double tightReportRatio = [] {
            const char* v = std::getenv("MVB_TIGHTARC_DIAG");
            if (!v) return 0.0;
            const double f = std::atof(v);
            return f > 0.0 ? f : 2.0;
        }();
        if (tightReportRatio > 0.0 && radius <= r * tightReportRatio) {
            std::fprintf(stderr,
                "[tight-arc] '%s' radius=%.6f mm wireR=%.6f mm ratio=%.3f sweep=%.4f deg "
                "arclen=%.6f mm%s\n",
                pr.label.c_str(), radius * 1e3, r * 1e3, radius / r,
                pr.arc.sweep * 180.0 / kPi, std::abs(pr.arc.sweep) * radius * 1e3,
                radius <= r ? "  <-- REVOLVE CROSSES THE AXIS" : "");
        }
        double sgn = pr.arc.sweep < 0 ? -1.0 : 1.0;
        double ddA = (overA / radius) * sgn, ddB = (overB / radius) * sgn;
        double total = pr.arc.sweep + ddA + ddB;
        if (std::getenv("MVB_GROW_DIAG") && (overA > 0 || overB > 0)) {
            gp_XYZ vs = rotateXYZ(pr.arc.v0, pr.arc.axis, -ddA);
            gp_XYZ ve = rotateXYZ(pr.arc.v0, pr.arc.axis, -ddA + total);
            std::fprintf(stderr,
                "[grow-arc] '%s' overA=%.4f overB=%.4f grownStart=(%.4f,%.4f,%.4f) grownEnd=(%.4f,%.4f,%.4f)\n",
                pr.label.c_str(), overA * 1e3, overB * 1e3,
                (pr.arc.c.X() + vs.X()) * 1e3, (pr.arc.c.Y() + vs.Y()) * 1e3,
                (pr.arc.c.Z() + vs.Z()) * 1e3, (pr.arc.c.X() + ve.X()) * 1e3,
                (pr.arc.c.Y() + ve.Y()) * 1e3, (pr.arc.c.Z() + ve.Z()) * 1e3);
        }
        if (std::abs(total) > 1.9 * kPi) return {};  // a near-full revolve would self-close: caller retries
        gp_XYZ v0e = rotateXYZ(pr.arc.v0, pr.arc.axis, -ddA);
        gp_Pnt start(pr.arc.c.XYZ() + v0e);
        gp_XYZ tangent = pr.arc.axis.Crossed(v0e);
        // Same for the arc: revolving a CLOSED circle makes a torus face that wraps in the
        // profile direction. The split profile revolves to ordinary patches instead.
        TopoDS_Wire profWire;
        double arcPhase = 0.0;
        // TRIED AND REVERTED (2026-08-30): phase-adopting GROWN arcs at mitred junctions
        // (reading the phase at the ungrown junction and rotating the profile back by -ddA),
        // on the theory that canonical arcs between phase-adopted prisms create the measured
        // 1.6 um rim mismatches. Measured: the rotation SPREAD into the corner pieces (run
        // rotation doubled to +-3.2 um, corners left canonical positions entirely) and 04
        // stayed 5/418. Arcs adopt only at ungrown tangent junctions, as before.
        if (startProfileOverride && !startProfileOverride->IsNull() && segments > 0 &&
            overA <= 0.0 &&
            sectionPhaseOf(*startProfileOverride, start, gp_Dir(tangent), r, segments,
                           arcPhase)) {
            profWire = wireProfileWirePhased(start, gp_Dir(tangent), r, segments, arcPhase);
            if (startPhasedOut) *startPhasedOut = true;
        } else if (segments > 0) {
            profWire = wireProfileWire(start, gp_Dir(tangent), r, segments);
        }
        // Tight arcs are revolved exactly even in faceted mode -- see kFacetTightArcExactRatio.
        // The exact piece carries no polygon phase, so a neighbour cannot adopt one from it.
        if (!profWire.IsNull() && radius < kFacetTightArcExactRatio * r) {
            if (std::getenv("MVB_MITRE_DIAG"))
                std::cerr << "[tight-arc] '" << pr.label << "' bend " << radius * 1e6 << " um < "
                          << kFacetTightArcExactRatio << " x r (" << r * 1e6
                          << " um): exact round revolve\n";
            profWire.Nullify();
            if (startPhasedOut) *startPhasedOut = false;
        }
        TopoDS_Face prof =
            !profWire.IsNull()
                ? BRepBuilderAPI_MakeFace(profWire).Face()
                : BRepBuilderAPI_MakeFace(wireProfileWireSplit(start, gp_Dir(tangent), r, 2)).Face();
        BRepPrimAPI_MakeRevol rev(prof, gp_Ax1(pr.arc.c, gp_Dir(pr.arc.axis)), total);
        if (rev.IsDone() && std::getenv("MVB_CAP_DEBUG")) {
            // Where the revolve's PLANAR caps land against the arc's analytic ends (the same
            // question [pipe-cap] asks of a sweep): centre offset and normal tilt at each end.
            const gp_Pnt endPt(pr.arc.c.XYZ() + rotateXYZ(pr.arc.v0, pr.arc.axis, pr.arc.sweep));
            const gp_Dir tStart = primFwdStart(pr, r), tEnd = primFwdEnd(pr, r);
            for (TopExp_Explorer fe(rev.Shape(), TopAbs_FACE); fe.More(); fe.Next()) {
                BRepAdaptor_Surface sf(TopoDS::Face(fe.Current()));
                if (sf.GetType() != GeomAbs_Plane) continue;
                GProp_GProps fg;
                BRepGProp::SurfaceProperties(fe.Current(), fg);
                const gp_Pnt c = fg.CentreOfMass();
                const gp_Dir n = sf.Plane().Axis().Direction();
                const double dS = c.Distance(start), dE = c.Distance(endPt);
                const bool isStart = dS < dE;
                const gp_Dir want = isStart ? tStart : tEnd;
                std::cerr << "[revolve-cap] '" << pr.label << "' " << (isStart ? "start" : "end")
                          << " cap: centre off nominal " << (isStart ? dS : dE) * 1e6
                          << " um, normal off analytic tangent "
                          << std::min(n.Angle(want), n.Angle(want.Reversed())) * 180.0 / kPi
                          << " deg, area " << fg.Mass() * 1e6 << " mm2 vs pi r^2 "
                          << kPi * r * r * 1e6 << "\n";
            }
        }
        if (rev.IsDone() && endCapOut && !profWire.IsNull()) {
            // A revolve's end section is EXACTLY the start profile rotated by the sweep.
            gp_Trsf rot;
            rot.SetRotation(gp_Ax1(pr.arc.c, gp_Dir(pr.arc.axis)), total);
            *endCapOut = TopoDS::Wire(
                BRepBuilderAPI_Transform(profWire, rot, Standard_True).Shape());
        }
        return rev.IsDone() ? rev.Shape() : TopoDS_Shape();
    }
    // BLEND / SPIRAL: swept as exact pipes. A pipe has no revolve/extrude axis to overhang
    // along, so mitre growth EXTENDS THE SPINE instead: a short straight run along the analytic
    // end tangent (G1-continuous, so MakePipeShell sweeps it as one smooth body) which the
    // bisector cut then trims back -- identical overhang semantics to the SEG/ARC3 branches.
    // Without this, a corner cut at a spiral/blend end carved a wedge with nothing to fill it
    // (measured: every lead<->wrap corner on PQ33 left a gap -> 5 disconnected components).
    auto pts = samplePrim(pr, r);
    if (pts.size() < 2) return {};
    // SHORT-SPINE REPORT (diagnostic only -- changes no geometry). MakePipeShell sweeps the
    // profile along this spine; the only guard above is "at least two points". A spine whose
    // LENGTH is at or below the profile DIAMETER cannot be swept cleanly -- the inner laterals
    // converge and cross, giving a solid that is BRepCheck-VALID, watertight and overlap-clean
    // while being BOPAlgo SELF-INTERSECTING. That is ABT #958: 02_flyback ships 42 such chunks,
    // exactly one per turn of Primary parallel 0, each 18 faces [bspline:16, plane:2] and each
    // ~0.42 mm across against a 0.3923 mm wire. Measured here rather than assumed: the ARC3
    // revolve was cleared first (zero tight-arc reports on 02), so the fold is on THIS path.
    if (std::getenv("MVB_SHORTSPINE_DIAG")) {
        double spineLen = 0.0;
        for (size_t i = 1; i < pts.size(); ++i) spineLen += pts[i - 1].Distance(pts[i]);
        if (spineLen <= 4.0 * r) {
            std::fprintf(stderr,
                "[short-spine] '%s' kind=%s spineLen=%.6f mm wireR=%.6f mm ratio(len/diam)=%.3f "
                "pts=%zu%s\n",
                pr.label.c_str(),
                pr.kind == Primitive::SPIRAL ? "SPIRAL" :
                pr.kind == Primitive::BLEND  ? "BLEND"  : "other",
                spineLen * 1e3, r * 1e3, spineLen / (2.0 * r), pts.size(),
                spineLen <= 2.0 * r ? "  <-- SHORTER THAN THE WIRE IS THICK" : "");
        }
    }
    // Where the growth can live INSIDE the curve (spirals on an analytic surface), put it there:
    // the spine is then a single edge and the pipe a single lateral face, with no seam for the
    // mitre plane to land on. See primEdgeGrowsAnalytically.
    const bool analyticGrowth = primEdgeGrowsAnalytically(pr);
    TopoDS_Edge e = analyticGrowth ? primEdge(pr, r, overA, overB) : primEdge(pr, r);
    if (e.IsNull()) return {};
    try {
        const auto ends = primEndpoints(pr);
        const gp_Dir tA = primFwdStart(pr, r), tB = primFwdEnd(pr, r);
        BRepBuilderAPI_MakeWire mw;
        gp_Pnt spineStart = ends.first;
        gp_Dir spineDir = tA;
        if (analyticGrowth) {
            mw.Add(e);
            // The profile must sit square to the GROWN spine's start, which the extended pcurve
            // has moved; read both straight off the edge rather than re-deriving them.
            BRepAdaptor_Curve bac(e);
            gp_Pnt p0;
            gp_Vec d0;
            bac.D1(bac.FirstParameter(), p0, d0);
            if (d0.Magnitude() > 1e-12) {
                spineStart = p0;
                spineDir = gp_Dir(d0);
                if (spineDir.Dot(tA) < 0.0) spineDir.Reverse();
            }
        } else {
            if (overA > 0.0) {
                spineStart = gp_Pnt(ends.first.XYZ() - tA.XYZ() * overA);
                mw.Add(BRepBuilderAPI_MakeEdge(spineStart, ends.first).Edge());
            }
            mw.Add(e);
            if (overB > 0.0)
                mw.Add(BRepBuilderAPI_MakeEdge(
                           ends.second, gp_Pnt(ends.second.XYZ() + tB.XYZ() * overB)).Edge());
        }
        if (!mw.IsDone()) return {};
        TopoDS_Wire spine = mw.Wire();
        // ONE attempt, at the REQUESTED faceting -- no round-profile retry (Alf: no
        // fallbacks). A failed faceted sweep returns null and the caller throws naming the
        // primitive, so the geometry that produced it gets fixed instead of silently mixing
        // a round (periodic-surface, gmsh-hostile) piece into a faceted conductor.
        //
        // The sweep runs in a MILLIMETRE frame: OCC's MakePipeShell mis-frames
        // sub-millimetre faceted profiles in METRES (measured standalone: the identical
        // helix/profile that fails at x1 -- PipeNotDone on a 0.1 mm 12-gon, 10_emi/
        // 13_current_sense/20_iso/24_margin -- sweeps cleanly at x1000). Pure scaling is an
        // EXACT affine map both ways, so this is a numerical frame choice, not an
        // approximation. TODO: migrate the other MakePipeShell sites to the same frame as
        // they are touched.
        // A CLOSED circular profile sweeps into ONE periodic B-spline carrying a seam, which
        // gmsh refuses (ABT #332, and confirmed 2026-08-25: gmsh dispatches on the SURFACE's
        // periodicity, OCCFace.cpp:142). Measured on 14_dab at segments=0: 246 of its 453 faces
        // were full-period wraps from exactly this sweep. The split profile yields ordinary
        // patches with real boundaries -- the IDENTICAL solid, no extra boolean, no geometry
        // change -- and the NURBS pass then clears the residual periodicity. Faceted profiles
        // (segments > 0) are already non-periodic, so they keep the polygon.
        // PHASE-ONLY continuation (see wireProfileWirePhased): never hand a neighbour's raw
        // wire to MakePipeShell -- take its section phase, keep the exact nominal profile.
        TopoDS_Wire prof;
        double phase0 = 0.0;
        if (startProfileOverride && !startProfileOverride->IsNull() && segments > 0 &&
            overA <= 0.0 &&
            sectionPhaseOf(*startProfileOverride, spineStart, spineDir, r, segments, phase0)) {
            // EXACT phase continuation (the v18-green behaviour for 10_emi's pipes). The
            // half-facet stagger belongs only at tangent junctions with a STRAIGHT source --
            // see the stagger rule in the assembler loop; pipes are curved receivers and a
            // staggered chord cuts their tightly curved strips (measured: 10_emi fell to 5/50
            // with the stagger applied here).
            prof = wireProfileWirePhased(spineStart, spineDir, r, segments, phase0);
            if (startPhasedOut) *startPhasedOut = true;
            if (std::getenv("MVB_CAP_DEBUG"))
                std::cerr << "[cap-debug] pipe start phase continued (" << phase0 << " rad)\n";
        } else {
            prof = segments > 0 ? wireProfileWire(spineStart, spineDir, r, segments)
                                : wireProfileWireSplit(spineStart, spineDir, r, 2);
        }
        gp_Trsf up, down;
        up.SetScale(gp_Pnt(0, 0, 0), 1000.0);
        down.SetScale(gp_Pnt(0, 0, 0), 1.0 / 1000.0);
        TopoDS_Wire spineMm =
            TopoDS::Wire(BRepBuilderAPI_Transform(spine, up, Standard_True).Shape());
        TopoDS_Wire profMm =
            TopoDS::Wire(BRepBuilderAPI_Transform(prof, up, Standard_True).Shape());
        BRepOffsetAPI_MakePipeShell ps(spineMm);
        // MVB_PIPE_MODE: a MEASUREMENT harness for ABT #958, not a tuning knob. The sweep has
        // always run on OCC defaults, and the default transition is BRepBuilderAPI_Transformed,
        // which the OCC docs describe as "assuming the result is a self-intersected shell" --
        // i.e. the default tolerates exactly the defect we are chasing. The corner that folds is
        // geometrically VALID but marginal: R/r = 1.073, so the swept tube's inner surface
        // passes 0.0147 mm from the corner axis on a 0.4 mm wire. Sampling is NOT the deficit
        // (chord sag there is 0.0011 mm, a 14x margin), so the fold comes from the sweep's own
        // approximation. Each variant below preserves the SHAPE -- none moves the spine or the
        // profile -- so whichever removes the self-intersections does so without altering the
        // part, which is the requirement. Unset = today's behaviour, so this is inert by default.
        if (const char* pm = std::getenv("MVB_PIPE_MODE")) {
            const std::string m = pm;
            if (m == "rightcorner")      ps.SetTransitionMode(BRepBuilderAPI_RightCorner);
            else if (m == "roundcorner") ps.SetTransitionMode(BRepBuilderAPI_RoundCorner);
            else if (m == "frenet")      ps.SetMode(Standard_True);   // true Frenet
            else if (m == "cfrenet")     ps.SetMode(Standard_False);  // corrected Frenet
            else if (m == "approxc1")    ps.SetForceApproxC1(Standard_True);
            else if (m == "tight")       ps.SetTolerance(1e-7, 1e-7, 1e-4);
        }
        ps.Add(profMm);
        ps.Build();
        // MakeSolid() is not idempotent -- a second call returns false -- so it runs exactly
        // once, here, and every later test reads the result.
        const bool pipeMade = ps.IsDone() && ps.MakeSolid();
        if (pipeMade && std::getenv("MVB_CAP_DEBUG")) {
            // Where the sweep's PLANAR caps actually are, against the analytic ends: a pipe whose
            // cap is tilted or displaced from its nominal section is what leaves wedges at
            // junctions the assembler believes are tangent.
            const gp_Pnt endNominal = ends.second;
            for (TopExp_Explorer fe(ps.Shape(), TopAbs_FACE); fe.More(); fe.Next()) {
                BRepAdaptor_Surface sf(TopoDS::Face(fe.Current()));
                if (sf.GetType() != GeomAbs_Plane) continue;
                GProp_GProps fg;
                BRepGProp::SurfaceProperties(fe.Current(), fg);
                gp_Pnt c = fg.CentreOfMass().Transformed(down);
                gp_Dir n = sf.Plane().Axis().Direction();
                const double dStart = c.Distance(spineStart), dEnd = c.Distance(endNominal);
                const bool isStart = dStart < dEnd;
                const gp_Dir want = isStart ? spineDir : tB;
                std::cerr << "[pipe-cap] '" << pr.label << "' " << (isStart ? "start" : "end")
                          << " cap: centre off nominal " << (isStart ? dStart : dEnd) * 1e6
                          << " um, normal off analytic tangent "
                          << std::min(n.Angle(want), n.Angle(want.Reversed())) * 180.0 / kPi
                          << " deg\n";
            }
        }
        if (pipeMade) {
            if (endCapOut) {
                // The sweep's own final section, downscaled -- the truth of what the next
                // piece must start from.
                try {
                    const TopoDS_Shape last = ps.LastShape();
                    // LastShape hands back the swept END SECTION but not always as a bare
                    // wire (measured: 90x "not a wire" on 10_emi) -- unwrap whatever it is.
                    TopoDS_Wire lastWire;
                    if (!last.IsNull()) {
                        if (last.ShapeType() == TopAbs_WIRE) {
                            lastWire = TopoDS::Wire(last);
                        } else {
                            TopExp_Explorer wx(last, TopAbs_WIRE);
                            if (wx.More()) lastWire = TopoDS::Wire(wx.Current());
                        }
                    }
                    if (!lastWire.IsNull())
                        *endCapOut = TopoDS::Wire(
                            BRepBuilderAPI_Transform(lastWire, down, Standard_True).Shape());
                    else if (std::getenv("MVB_CAP_DEBUG"))
                        std::cerr << "[cap-debug] pipe LastShape "
                                  << (last.IsNull() ? "NULL" : "carries no wire") << " for '"
                                  << pr.label << "'\n";
                } catch (const Standard_Failure& f) {
                    if (std::getenv("MVB_CAP_DEBUG"))
                        std::cerr << "[cap-debug] pipe LastShape threw ("
                                  << (f.GetMessageString() ? f.GetMessageString() : "?")
                                  << ") for '" << pr.label << "'\n";
                }
            }
            // VALIDATE THE SWEEP. MakePipeShell reports IsDone/MakeSolid without guaranteeing
            // a valid B-Rep, and this branch used to ship whatever it produced: measured on
            // 06_llc, 'Primary parallel 0 [solid 25]' (3.65 mm3, 80 faces) reached the STEP
            // BRepCheck-invalid, which is what made the whole design CAD DEFECTIVE. Repair it
            // with OCC's own ShapeFix and re-check; a sweep that stays invalid returns null so
            // the caller THROWS naming the primitive, rather than shipping bad copper.
            TopoDS_Shape swept = BRepBuilderAPI_Transform(ps.Shape(), down, Standard_True).Shape();
            if (!swept.IsNull() && !BRepCheck_Analyzer(swept).IsValid()) {
                GProp_GProps gs0;
                BRepGProp::VolumeProperties(swept, gs0);
                try {
                    Handle(ShapeFix_Shape) fix = new ShapeFix_Shape(swept);
                    fix->Perform();
                    const TopoDS_Shape fixed = fix->Shape();
                    if (!fixed.IsNull() && BRepCheck_Analyzer(fixed).IsValid()) {
                        GProp_GProps gs1;
                        BRepGProp::VolumeProperties(fixed, gs1);
                        // The repair may only clean the representation, never move copper.
                        if (gs0.Mass() > 0 &&
                            std::abs(gs1.Mass() - gs0.Mass()) <= 1e-3 * gs0.Mass()) {
                            std::cerr << "[sweep-fix] '" << pr.label
                                      << "' swept solid was invalid; ShapeFix repaired it"
                                      << std::endl;
                            return fixed;
                        }
                    }
                } catch (const Standard_Failure&) {
                }
                std::cerr << "[sweep-fail] '" << pr.label
                          << "' swept solid is an invalid B-Rep and could not be repaired"
                          << std::endl;
                return {};
            }
            return swept;
        }
        if (std::getenv("MVB_DIAG"))
            std::cerr << "[sweep-fail] '" << pr.label << "' MakePipeShell status=" << (int)ps.GetStatus()
                      << " isDone=" << ps.IsDone() << " segments=" << segments
                      << " r=" << r << "\n";
    } catch (const Standard_Failure& f) {
        if (std::getenv("MVB_DIAG"))
            std::cerr << "[sweep-fail] '" << pr.label << "' OCC exception: "
                      << (f.GetMessageString() ? f.GetMessageString() : "(null)") << "\n";
    }
    return {};
}

// SCALING A SHAPE SCALES ITS TOLERANCE TOO -- AND A 1000x TOLERANCE BLINDS THE BOOLEAN.
// ABT #860 (2026-08-23). Every mitre boolean runs in the MILLIMETRE frame because OCC mis-frames
// sub-millimetre profiles in metres. BRepBuilderAPI_Transform carries the shape's tolerances
// through the same scale, so a solid holding OCC's own confusion (1e-7) in metres arrives in the
// millimetre frame claiming 1e-4 -- a tenth of a micron of slop on geometry that is analytically
// exact. At that tolerance OCC 7.9.3 declares a TOROIDAL face and the knife's PLANE
// non-intersecting: Cut returns its argument UNCHANGED and Common returns EMPTY, both reporting
// IsDone with no error and no warning. Measured on 24_margin's entrance corner (--segments 0):
// 0.0581 mm3 of the terminal stub stands inside the knife box by point classification, and both
// booleans see nothing; with the tolerance reset to 1e-7 the same Cut removes 0.0583 mm3. That is
// the whole '--segments 0 mitred corner did not abut' failure -- the faceted profile never hit it
// because a plane/plane intersection survives the inflated tolerance that a plane/torus one does
// not.
// A shape does not become less accurate by being measured in millimetres, so the millimetre copy
// is given the millimetre frame's confusion -- and only if it is still VALID with it, because a
// swept BSpline pipe may genuinely need the slop it carries.
static const gp_Trsf& mitreUpTrsf() {
    static const gp_Trsf t = [] { gp_Trsf x; x.SetScale(gp_Pnt(0, 0, 0), 1000.0); return x; }();
    return t;
}
static const gp_Trsf& mitreDownTrsf() {
    static const gp_Trsf t = [] { gp_Trsf x; x.SetScale(gp_Pnt(0, 0, 0), 1.0 / 1000.0); return x; }();
    return t;
}
static TopoDS_Shape toMitreFrame(const TopoDS_Shape& metres) {
    TopoDS_Shape mm = BRepBuilderAPI_Transform(metres, mitreUpTrsf(), Standard_True).Shape();
    if (mm.IsNull()) return mm;
    TopoDS_Shape tight = BRepBuilderAPI_Copy(mm).Shape();
    ShapeFix_ShapeTolerance().SetTolerance(tight, Precision::Confusion());
    if (!tight.IsNull() && BRepCheck_Analyzer(tight).IsValid()) return tight;
    return mm;   // the shape really does need its slop: leave it, the probe below still judges it
}
// Back to metres. The inverse scale divides the tolerances too, which would leave the result
// claiming 1e-10 -- below the kernel's own confusion, which OCC reports as an invalid tolerance.
// Raise the too-small ones back to confusion and leave anything larger alone.
static TopoDS_Shape fromMitreFrame(const TopoDS_Shape& mm) {
    TopoDS_Shape metres = BRepBuilderAPI_Transform(mm, mitreDownTrsf(), Standard_True).Shape();
    if (!metres.IsNull()) ShapeFix_ShapeTolerance().LimitTolerance(metres, Precision::Confusion());
    return metres;
}

// DID THE MITRE ACTUALLY HAPPEN? ABT #860 (2026-08-23). Every previous verdict on a mitre cut came
// from the SAME boolean kernel that performs it -- 'removed' from the Cut, 'shouldRemove' from a
// Common with the knife -- so when the kernel goes blind (above) it reports success AND reports
// nothing left to remove, and the untrimmed overhang ships. The suggested cure of probing with a
// large half-space is blind in exactly the same way: measured on the same corner, Common against a
// 4691 mm3 box returned 0 mm3 while point classification found 0.058 mm3 of the solid inside it.
// So the check must not be a boolean at all. It asks the trimmed solid a question with a
// yes/no answer: is a point on the piece's OWN CENTRELINE, past the bisector plane, still inside
// it? The overhang runs from the piece's un-grown end `stubEnd` along `stubDir` for `grow`; a
// correct mitre leaves none of it. Only points standing at least kProbeClear past the plane are
// asked, so the answer never depends on the tolerance of the mitre face itself; an overhang
// thinner than that is thinner than the modelling tolerance and there is nothing to detect.
// Returns the largest stand-off (in metres) at which copper was found beyond the plane, or -1
// when the piece is clean.
static double overhangBeyondPlane(const TopoDS_Shape& trimmed, const gp_Pnt& stubEnd,
                                  const gp_Dir& stubDir, const gp_Pnt& P, const gp_XYZ& w,
                                  double grow) {
    if (trimmed.IsNull() || grow <= 0.0) return -1.0;
    const double rate = stubDir.XYZ().Dot(w);          // how fast the overhang leaves the plane
    if (rate <= 1e-9) return -1.0;                     // it does not head into the discard side
    const double s0 = (stubEnd.XYZ() - P.XYZ()).Dot(w);
    const double kProbeClear = 1e-6;                   // 1 um: 10x the shapes' own tolerance
    const double tCross = (kProbeClear - s0) / rate;   // first t standing clear of the plane
    if (tCross >= grow) return -1.0;                   // the grown end never reaches past it
    const double tFrom = std::max(0.0, tCross);
    double worst = -1.0;
    for (double f : {0.30, 0.55, 0.80, 0.97}) {
        const double t = tFrom + f * (grow - tFrom);
        const gp_Pnt probe(stubEnd.XYZ() + stubDir.XYZ() * t);
        const double stand = (probe.XYZ() - P.XYZ()).Dot(w);
        if (stand < kProbeClear) continue;
        // Per SOLID: a trim that fragmented the piece is a compound, and the classifier
        // silently misclassifies one of those.
        for (TopExp_Explorer ex(trimmed, TopAbs_SOLID); ex.More(); ex.Next()) {
            BRepClass3d_SolidClassifier c(TopoDS::Solid(ex.Current()), probe, Precision::Confusion());
            if (c.State() == TopAbs_IN) { worst = std::max(worst, stand); break; }
        }
    }
    return worst;
}

// Trim the grown overhang of `solid` that pokes past the mitre plane (through P, material side =
// keepDir): SUBTRACT a small knife box sitting on the discard side, localized to the junction.
// Neither of the two 'obvious' cuts is correct here:
//  - a global half-space keep is WRONG for near-closed wraps -- the wrap legitimately curves back
//    across the junction plane far from the joint, so the far side would be sliced off;
//  - the previous 'tight beam' knife (8r-wide box, depth spanning the solid, kept via Common) was
//    NOT a half-space: it discarded everything outside the beam laterally. Measured on PQ33: it
//    kept 17% of the first wrap and split it into 3 solids -> disconnected slivers.
// Locality bounds: the grown overhang lies within `grow`+r of P along -keepDir and within
// ~2.5r+`grow` of P laterally (bisector tilt <= ~50 deg), so a (3r+grow)-lateral x (grow+2r)-deep
// knife covers it with margin. A volume guard rejects any trim that removed more than a stub's
// worth of material (i.e. the knife reached a neighbouring pass of the same solid).
// ABT #685 (Alf, 2026-08-18): a failed mitre THROWS. "Don't allow stuff to return stuff when
// things fail, make them throw!" -- an untrimmed mitre is not a degraded result, it is copper
// standing proud of the corner, and in FEM that is fatal.
//
// MVB_MITRE_KEEP=1 is the DEBUGGING escape: the same complaint goes to stderr, in full, and the
// untrimmed solid is handed back so the geometry reaches disk and the fault can be LOOKED AT.
// Never a quiet fallback -- it shouts either way, and the loud path is opt-in.

// MITRED FACETED PRISM, built exactly -- no boolean. A faceted SEG tube with bisector-cut ends
// does not need the knife at all: every vertex of the section polygon travels along the tube
// axis, so its start/end are the LINE-PLANE intersections with the two end planes. Caps are
// planar polygons IN those planes; sides are planar quads ON the original facet planes. Sew,
// solid, done. This exists because the knife route is structurally unsafe exactly here: a
// 16-gon's facet normals include the knife-box axes, so knife faces land COPLANAR with facet
// planes and OCC's cut leaves self-intersecting trims -- measured on 05_pfc seg=16: the two
// straight tubes of EVERY turn (~40 pieces) carried 2-3 self-intersections each and 9/399
// meshed, while the same cuts on seg=0 arc-panel tubes were clean. Boolean-free construction
// is the codebase's own doctrine (conformal-by-construction); this applies it to the one
// primitive kind where the knife is provably the wrong tool.
// startCapShared / endCapOut: BIT-IDENTICAL JUNCTION CAPS. Adjacent prisms each project their
// own polygon onto the shared bisector plane, and the two projections agree only to fp
// precision -- the fragment then imprints nm slivers and micro-edges at every such junction
// (measured on 12_boost: the raw STEP has 0 sliver faces and 0 micro-edges, the fragmented
// model has 2 and 32, and one cap face degenerates outright -> "overlapping facets", 167/206
// volumes lost). When the previous piece hands its exact end-cap vertices in and they agree
// with this piece's own projection to 100 nm, ADOPT them verbatim: the two solids then share
// bit-identical vertices and the glue merges the caps exactly. A mismatch beyond 100 nm means
// the frames' phases differ at this corner -- adopting would twist the prism, so fall back to
// the own projection (today's behaviour).
static TopoDS_Shape mitredFacetPrism(const gp_Pnt& a, const gp_Dir& dir, double r, int segments,
                                     const gp_Pnt& Ps, const gp_XYZ& ns,
                                     const gp_Pnt& Pe, const gp_XYZ& ne,
                                     const std::vector<gp_Pnt>* startCapShared = nullptr,
                                     std::vector<gp_Pnt>* endCapOut = nullptr,
                                     bool sharedTrusted = false,
                                     bool* startAdoptedOut = nullptr,
                                     bool staggerAdopt = false,
                                     std::string* why = nullptr) {
    auto refuse = [&](const std::string& r) { if (why) *why = r; return TopoDS_Shape(); };
    if (segments <= 0) return refuse("segments <= 0");
    const double dS = dir.XYZ().Dot(ns), dE = dir.XYZ().Dot(ne);
    // The upstream bisector guarantee is tilt <= ~50 deg; below 0.2 the plane is close to
    // parallel to the axis and the projection blows up -- refuse and let the knife path run.
    if (std::abs(dS) < 0.2 || std::abs(dE) < 0.2) return refuse("cut plane too oblique (|d.n| < 0.2)");
    gp_Ax2 plane = deterministicSectionFrame(a, dir);
    const gp_Dir dx = plane.XDirection(), dy = plane.YDirection();
    const double offset = kPi / segments;
    std::vector<gp_Pnt> vs(segments);
    for (int i = 0; i < segments; ++i) {
        const double ang = kTwoPi * i / segments + offset;
        const gp_XYZ v = a.XYZ() + dx.XYZ() * (r * std::cos(ang)) + dy.XYZ() * (r * std::sin(ang));
        const double ts = (Ps.XYZ() - v).Dot(ns) / dS;
        vs[i] = gp_Pnt(v + dir.XYZ() * ts);
    }
    if (startCapShared && !handoffDisabled() && (int)startCapShared->size() == segments) {
        // The neighbour's ACTUAL end section, at ANY junction. The two pieces' deterministic
        // frames differ per axis, so their projections onto the shared plane are ROTATED
        // polygons -- measured on 12_boost: mitre-corner cap mismatches of 18 um..1.4 mm, only
        // ~50 of ~240 junctions coincidentally aligned. Adoption gates on GEOMETRY only:
        // vertices in the start plane, at the section's projected radius (up to r/cos(tilt),
        // tilt <= ~50 deg upstream => 1.6 r). The wire's facet phase then rotates continuously
        // through the chain -- physically arbitrary for copper, and every junction is exact.
        (void)sharedTrusted;
        {
            bool fits = true;
            double worstPlane = 0.0, worstRad = 0.0;
            const gp_Dir nsd(ns);
            // Plane tolerance 25% of the wire radius, from the two measured populations on
            // 10_emi (r = 0.1 mm): construction offsets -- where the sweep ACTUALLY ended --
            // are 12-13 um (0.13 r), genuine endpoint mismatches are 0.14-0.15 mm (1.5 r).
            // 0.25 r sits between them with a >5x margin each way. The cap IS the copper's
            // real end face; the next piece must start there, not at the nominal plane.
            const double planeTol = std::max(2e-6, 0.25 * r);
            // THE CAP MUST BELONG TO THIS PIECE'S OWN CYLINDER. Project every adopted vertex
            // onto the plane PERPENDICULAR to this piece's axis and require its distance from
            // the axis to be exactly the wire radius. This is analytic, not a fit: at a clean
            // symmetric corner the mitre plane's normal is n ~ dA + dB, and reflection through
            // it maps dA -> -dB, so it carries the neighbour's regular section onto a regular
            // section of THIS piece -- every projected vertex lands at exactly r. It fails in
            // exactly one case, which is the one that broke 04_forward: where an ENDPOINT
            // MISMATCH moved the plane off the true bisector (Ps = midpoint of the two
            // endpoints), the corner is no longer symmetric, the reflection no longer maps one
            // cylinder onto the other, and sweeping the adopted cap bulges the piece OUTSIDE
            // the nominal wire -- measured 24 um proud (solid 8 reached x=5.994 where the run's
            // surface is at 6.018), which is what grazed the neighbouring run's side facet and
            // cost 413 of 418 volumes. The old 0.8r..1.6r band could not see it: an oblique cap
            // legitimately spans r..r/cos(tilt), so a 1.10 r bulge passed.
            // BASELINE GATE -- the configuration measured green on 10/12/02. TRIED AND
            // REVERTED (2026-08-30): a Newell-normal tilt test (1e-3 rad) plus a strict
            // projected-radius shape test (|rad - r| <= 0.01 r perpendicular to this axis).
            // Both passed every one of 04's 464 adoptions -- the caps ARE geometrically exact
            // -- so they discriminate nothing there, and they were never validated on 10_emi.
            for (const auto& p : *startCapShared) {
                const gp_Vec d(Ps, p);
                worstPlane = std::max(worstPlane, std::abs(d.Dot(gp_Vec(nsd))));
                const gp_XYZ axisFoot = a.XYZ() + dir.XYZ() * gp_Vec(a, p).Dot(gp_Vec(dir));
                const double rad = (p.XYZ() - axisFoot).Modulus();
                worstRad = std::max(worstRad, std::abs(rad - r));
                if (std::abs(d.Dot(gp_Vec(nsd))) > planeTol || rad < 0.8 * r || rad > 1.6 * r)
                    fits = false;
            }
            if (std::getenv("MVB_CAP_DEBUG"))
                std::cerr << "[cap-debug] neighbour cap "
                          << (fits ? "ADOPTED" : "phase-only") << " planeOff=" << worstPlane
                          << " radOff=" << worstRad << "\n";
            if (fits && !capPositionsDisabled() && staggerAdopt) {
                // PIPE/REVOLVE -> PRISM: STAGGERED CONTINUATION, NOT VERBATIM ADOPTION
                // (2026-08-30, the 04_forward osculation class). Verbatim adoption aligns this
                // prism's facets with the neighbouring sweep's strips, and at a G1 fillet
                // junction the two surfaces then OSCULATE -- gmsh meshes both faces
                // independently into overlapping triangles (04: turn-0 entrance junction, PLC
                // "two facets intersect"). Keep what exactness is FOR -- the junction PLANE
                // (the mean axial offset mu, which is what killed 10_emi's 3.9e-5 mm^3 volume
                // lenses) -- and rotate the polygon a half facet, so the caps are two coplanar
                // 16-gons whose rim wedges are FLAT (zero volume, audit-clean) and the side
                // facets cross transversally, the configuration every control meshes.
                double mu = 0.0;
                for (const auto& p : *startCapShared) mu += gp_Vec(Ps, p).Dot(gp_Vec(nsd));
                mu /= static_cast<double>(startCapShared->size());
                const gp_XYZ d0 = startCapShared->front().XYZ() -
                                  (a.XYZ() + dir.XYZ() *
                                             gp_Vec(a, startCapShared->front()).Dot(gp_Vec(dir)));
                const double px = d0.Dot(dx.XYZ()), py = d0.Dot(dy.XYZ());
                const double phi = std::atan2(py, px) - offset + kPi / segments;
                for (int i = 0; i < segments; ++i) {
                    const double ang = kTwoPi * i / segments + offset + phi;
                    const gp_XYZ v =
                        a.XYZ() + dx.XYZ() * (r * std::cos(ang)) + dy.XYZ() * (r * std::sin(ang));
                    const double ts = (Ps.XYZ() + nsd.XYZ() * mu - v).Dot(ns) / dS;
                    vs[i] = gp_Pnt(v + dir.XYZ() * ts);
                }
                if (std::getenv("MVB_CAP_DEBUG"))
                    std::cerr << "[cap-debug]   staggered continuation (mu=" << mu
                              << " phi=" << phi << ")\n";
            } else if (fits && !capPositionsDisabled()) {
                vs = *startCapShared;
                // PRISM -> PRISM (mitred or collinear): verbatim adoption stays -- both sides
                // are planar-faceted and meet at a dihedral or as a straight continuation, so
                // there is no osculation, and 12_boost REQUIRES the exactly shared face
                // (262/262 with it, 1/262 without).
                if (startAdoptedOut) *startAdoptedOut = true;
            } else {
                // POSITIONS REJECTED. Phase may continue ONLY at a TANGENT junction, where the
                // neighbour's section projects into this piece's start plane as a rigid
                // rotation and its vertex azimuths are the copper's real facet azimuths. At a
                // MITRED junction the projection is oblique (the section becomes an ellipse in
                // the bisector plane), so a phase read off one vertex aligns that vertex and
                // leaves the rest off -- measured on 04_forward: 413 of 418 volumes lost to a
                // 50x95 um sliver facet and a 30 um-wide split side strip at a corner
                // (surfaces 220/237). There the deterministic frame is already shared by both
                // sides, so it stays the right answer.
                // TRIED AND REVERTED (2026-08-31, 20_iso_buckboost_epc25_3c91): continuing the
                // phase here through the BISECTOR REFLECTION of the neighbour's cap (an
                // isometry that maps its end tangent onto this piece's axis, so its azimuths
                // ARE this cylinder's facet azimuths -- checked by requiring every reflected
                // vertex to land at radius r, measured 0.016-0.023 r on 20's six mitred pipe
                // corners). It did what it claimed -- the junction 16-gons' 0.0516 rad relative
                // rotation, 9.0 um at the rim, collapsed to a 3.8 um offset -- and 20 STAYED
                // RED, because the phase was never the defect: the WORKING mirror-image
                // entrance corner has a 14.5 um worst vertex mismatch, larger than either. The
                // defect was the neighbour's knifed face being NON-PLANAR (see curveGrow).
                const bool tangentStart = std::abs(dS) > 1.0 - 1e-6 && phaseContinuationEnabled();
                const gp_XYZ d0 = startCapShared->front().XYZ() -
                                  (a.XYZ() + dir.XYZ() *
                                             gp_Vec(a, startCapShared->front()).Dot(gp_Vec(dir)));
                const double px = d0.Dot(dx.XYZ()), py = d0.Dot(dy.XYZ());
                if (tangentStart && px * px + py * py > 0.25 * r * r) {
                    const double phi = std::atan2(py, px) - offset;
                    for (int i = 0; i < segments; ++i) {
                        const double ang = kTwoPi * i / segments + offset + phi;
                        const gp_XYZ v = a.XYZ() + dx.XYZ() * (r * std::cos(ang)) +
                                         dy.XYZ() * (r * std::sin(ang));
                        const double ts = (Ps.XYZ() - v).Dot(ns) / dS;
                        vs[i] = gp_Pnt(v + dir.XYZ() * ts);
                    }
                    if (std::getenv("MVB_CAP_DEBUG"))
                        std::cerr << "[cap-debug]   phase continued (" << phi << " rad)\n";
                }
            }
        }
    } else if (startCapShared && std::getenv("MVB_CAP_DEBUG")) {
        std::cerr << "[cap-debug] shared start cap SIZE mismatch (" << startCapShared->size()
                  << " vs " << segments << ")\n";
    }
    // The end cap is the FINAL start polygon projected to the end plane -- computed from vs
    // so an adopted (phase-shifted) start cap yields a straight, untwisted prism.
    std::vector<gp_Pnt> ve(segments);
    for (int i = 0; i < segments; ++i) {
        const double te = (Pe.XYZ() - vs[i].XYZ()).Dot(ne) / dE;
        // The span must exceed what the B-Rep can RESOLVE, not merely zero. A cap adoption can
        // rotate the polygon so one vertex lands within nanometres of where the two end planes
        // cross; a 1 nm-tall side quad then collapses to coincident mesh nodes and gmsh emits
        // one exactly-degenerate tet (measured on 04_forward: SICN -2.2e-14 at
        // (6.0034, -8.2317, 5.9955) mm INSIDE the piece, the whole design's NOT_MESHABLE).
        // Refusing here falls back to the knife path, whose trimmed piece is transversal at
        // that corner -- exactly the configuration the handoff-off control meshes green.
        if (!(te > Precision::Confusion())) return {};
        ve[i] = gp_Pnt(vs[i].XYZ() + dir.XYZ() * te);
    }
    // ThruSections between the two cap polygons, RULED: OCC's own constructor for exactly this
    // ruled prism, with face orientations it manages itself. The first implementation hand-sewed
    // caps and side quads from MakeFace defaults, which enforces nothing about a consistently
    // outward shell -- and every prism came back FAULTY from BOPAlgo_ArgumentAnalyzer (88
    // "self-intersections" on 05_pfc that were really orientation conflicts). Build the two
    // wires in the SAME vertex order; ruled lofting connects vertex i to vertex i.
    try {
        BRepBuilderAPI_MakePolygon ws, we;
        for (const auto& q : vs) ws.Add(q);
        ws.Close();
        for (const auto& q : ve) we.Add(q);
        we.Close();
        if (!ws.IsDone() || !we.IsDone()) return {};
        BRepOffsetAPI_ThruSections loft(Standard_True /*solid*/, Standard_True /*ruled*/);
        loft.AddWire(ws.Wire());
        loft.AddWire(we.Wire());
        loft.Build();
        if (!loft.IsDone() || loft.Shape().IsNull()) return refuse("loft failed");
        const TopoDS_Shape out = loft.Shape();
        if (!BRepCheck_Analyzer(out).IsValid()) return refuse("loft invalid (BRepCheck)");
        // SUB-RESOLUTION SLIVER GUARD. "> 0" let through prisms whose two end planes almost
        // coincide: positive volume below what the B-Rep can represent (measured on the
        // mitre-corner toroid: 18-face solids at 0.000000 mm3, flagged DEFECTIVE by stage A).
        // Both bounds are DERIVED: the mean plane-to-plane span along the axis must exceed the
        // model's own length resolution, and the volume must be at least half of
        // section-area x span -- a true prism is exactly area x span, so half catches any
        // degenerate or folded loft without a tuned constant.
        double spanSum = 0, spanMin = std::numeric_limits<double>::max();
        for (int i = 0; i < segments; ++i) {
            const double sp = (ve[i].XYZ() - vs[i].XYZ()).Dot(dir.XYZ());
            spanSum += sp;
            spanMin = std::min(spanMin, sp);
        }
        const double span = spanSum / segments;
        if (span <= Precision::Confusion()) return refuse("non-positive mean span " + std::to_string(span * 1e6) + " um");
        // FOLDED LOFT GUARD (2026-09-05, 13_current_sense at --segments 12). When the piece is
        // shorter than its two mitre bevels the end planes CROSS inside the section: some polygon
        // vertices end up with a NEGATIVE axial span, the side quads between them are bow-ties,
        // and the loft is a self-intersecting solid that BRepCheck accepts and the volume gate
        // (mean span, total volume) cannot see -- two such pieces (15 um and 105 um long on
        // 97 um and 490 um wires) reached the STEP as 3-self-intersection crumbs. The mean span
        // says nothing about this; every vertex must clear its own plane.
        if (spanMin <= Precision::Confusion()) return refuse("folded loft: min vertex span " + std::to_string(spanMin * 1e6) + " um (mean " + std::to_string(span * 1e6) + " um)");
        // Inscribed polygon area: n/2 * r^2 * sin(2 pi / n).
        const double area = 0.5 * segments * r * r * std::sin(kTwoPi / segments);
        GProp_GProps g;
        BRepGProp::VolumeProperties(out, g);
        if (g.Mass() < 0.5 * area * span) return refuse("volume " + std::to_string(g.Mass() * 1e9) + " mm3 < half of area x span " + std::to_string(0.5 * area * span * 1e9) + " mm3");
        if (endCapOut) *endCapOut = ve;
        return out;
    } catch (const Standard_Failure& e) {
        return refuse(std::string("OCCT exception: ") + (e.GetMessageString() ? e.GetMessageString() : "?"));
    }
}

// MITRED ROUND PRISM, on analytic geometry only (ABT #961, the round twin of the facet prism
// above). A round SEG with a mitred end used to be a SPLIT-PROFILE cylinder sliced by a box
// knife, and the knife is structurally unsafe there: the bisector plane's intersection ellipse
// crosses the profile's seam, and OCC's cut then leaves a sliver -- measured on 11_pushpull's
// exit lead seg 1: a 13.7 um edge whose two vertices carry 33 um and 15 um tolerances, which
// BOPAlgo (correctly) reports as a self-intersecting solid; 14_dab's entrance lead seg 1 on both
// parallels is the same defect. There is nothing to approximate here: the piece is an analytic
// cylinder between two planes. So it is built as exactly that -- the split-profile prism (two
// exact half-cylinder faces, never a periodic surface) sliced by HALF-SPACES. Plane against an
// analytic cylinder is the well-conditioned boolean -- every knife failure in this file is a
// BOX against a swept B-spline, or a box face landing on a seam -- and a half-space has no box
// faces to land anywhere.
// TRIED FIRST (2026-09-02): a ruled ThruSections loft between the two elliptical sections. Its
// rulings are the cylinder's generators only if the two sections share their parametrisation,
// and OCC re-approximates the conics on the way in: the loft came out 2 % short in volume on a
// 45-degree mitre. The analytic route has nothing to re-approximate.
//   Ps/ns, Pe/ne: the start/end cut planes (a point and a normal each). A flush end passes its
//   own perpendicular plane.
// Every check is exact-or-refuse: an invalid result, or a volume off the analytic pi r^2 L by
// more than the boolean's own precision, returns null and the caller falls back to the knife,
// saying so.
static TopoDS_Shape mitredRoundPrism(const gp_Pnt& a, const gp_Dir& dir, double r,
                                     const gp_Pnt& Ps, const gp_XYZ& ns,
                                     const gp_Pnt& Pe, const gp_XYZ& ne,
                                     std::string* why = nullptr) {
    auto refuse = [&](const std::string& reason) {
        if (why) *why = reason;
        return TopoDS_Shape();
    };
    const gp_XYZ d = dir.XYZ();
    const double dS = d.Dot(ns), dE = d.Dot(ne);
    if (std::abs(dS) < 0.2 || std::abs(dE) < 0.2) return refuse("cut plane too oblique");
    // Axial stations of the two cut planes on the axis, and the true truncated-cylinder volume.
    const double tS = (Ps.XYZ() - a.XYZ()).Dot(ns) / dS;
    const double tE = (Pe.XYZ() - a.XYZ()).Dot(ne) / dE;
    const double span = tE - tS;
    if (span <= Precision::Confusion()) return refuse("non-positive axial span");
    const double expect = kPi * r * r * span;
    try {
        // Work in the millimetre frame, like every boolean in this file.
        gp_Trsf up, down;
        up.SetScale(gp_Pnt(0, 0, 0), 1000.0);
        down.SetScale(gp_Pnt(0, 0, 0), 1.0 / 1000.0);
        // A cylinder long enough to contain both cut ellipses (each reaches r*tan(tilt) past its
        // plane's axial station; tilt <= ~78 deg from the 0.2 guard, so 5 r covers it).
        const double reach = 5.0 * r;
        const gp_Pnt base(a.XYZ() + d * (tS - reach));
        // The SPLIT-profile prism the flush segments already use: two exact half-cylinder faces,
        // no periodic surface, and -- unlike the box knife -- a half-space plane meets the seam
        // lines transversally, so there is nothing for the cut to graze.
        TopoDS_Face profile =
            BRepBuilderAPI_MakeFace(wireProfileWireSplit(base, dir, r, 2)).Face();
        TopoDS_Shape cyl = BRepPrimAPI_MakePrism(profile, gp_Vec(d * (span + 2.0 * reach))).Shape();
        cyl = BRepBuilderAPI_Transform(cyl, up, Standard_True).Shape();
        // Keep the side of each plane that contains the axis midpoint between the planes.
        const gp_Pnt mid(a.XYZ() + d * (0.5 * (tS + tE)));
        auto slice = [&](const TopoDS_Shape& body, const gp_Pnt& P, const gp_XYZ& n) {
            gp_Pln pln(P, gp_Dir(n));
            TopoDS_Face f = BRepBuilderAPI_MakeFace(pln).Face();
            f = TopoDS::Face(BRepBuilderAPI_Transform(f, up, Standard_True).Shape());
            const gp_Pnt midMm = mid.Transformed(up);
            // MakeHalfSpace keeps the side containing the reference point, so this is the
            // half-space to REMOVE: build it about the point mirrored through the plane.
            const gp_XYZ nn = gp_Dir(n).XYZ();
            const double h = (midMm.XYZ() - P.Transformed(up).XYZ()).Dot(nn);
            const gp_Pnt away(midMm.XYZ() - nn * (2.0 * h));
            TopoDS_Solid half = BRepPrimAPI_MakeHalfSpace(f, away).Solid();
            BRepAlgoAPI_Cut cut(body, half);
            if (!cut.IsDone() || cut.Shape().IsNull()) return TopoDS_Shape();
            return cut.Shape();
        };
        TopoDS_Shape s1 = slice(cyl, Ps, ns);
        if (s1.IsNull()) return refuse("start-plane cut failed");
        TopoDS_Shape s2 = slice(s1, Pe, ne);
        if (s2.IsNull()) return refuse("end-plane cut failed");
        // One solid expected.
        TopoDS_Shape out;
        int nSolids = 0;
        for (TopExp_Explorer ex(s2, TopAbs_SOLID); ex.More(); ex.Next(), ++nSolids) out = ex.Current();
        if (nSolids != 1) return refuse("cut produced " + std::to_string(nSolids) + " solids");
        if (!BRepCheck_Analyzer(out).IsValid()) return refuse("cut result invalid (BRepCheck, mm frame)");
        out = BRepBuilderAPI_Transform(out, down, Standard_True).Shape();
        if (!BRepCheck_Analyzer(out).IsValid()) return refuse("result invalid after downscale (BRepCheck)");
        GProp_GProps g;
        BRepGProp::VolumeProperties(out, g);
        // The half-space cut works in the mm frame at OCC's default tolerance: measured 1.2e-4
        // relative on 10_emi's 0.007 mm3 risers (nanometres of length), so the gate sits at 1e-3
        // -- against a collapsed or doubled result, not against the boolean's own precision.
        if (std::abs(g.Mass() - expect) > 1e-3 * expect)
            return refuse("volume " + std::to_string(g.Mass() * 1e9) + " mm3 vs exact " +
                          std::to_string(expect * 1e9) + " mm3");
        return out;
    } catch (const Standard_Failure&) {
        return refuse("OCCT exception");
    }
}

// Ordered polygon vertices of a closed wire (edge order), and back. The section-wire handoff
// moves junction sections between the prism path (points) and the pipe/revolve path (wires).
static std::vector<gp_Pnt> wireToPoints(const TopoDS_Wire& w) {
    std::vector<gp_Pnt> pts;
    for (BRepTools_WireExplorer we(w); we.More(); we.Next())
        pts.push_back(BRep_Tool::Pnt(we.CurrentVertex()));
    return pts;
}
static TopoDS_Wire pointsToWire(const std::vector<gp_Pnt>& pts) {
    if (pts.size() < 3) return {};
    BRepBuilderAPI_MakePolygon poly;
    for (const auto& p : pts) poly.Add(p);
    poly.Close();
    if (!poly.IsDone()) return {};
    return BRepBuilderAPI_MakeWire(poly.Wire()).Wire();
}

static TopoDS_Shape refuseTrim(const std::string& why, const TopoDS_Shape& untrimmed) {
    if (std::getenv("MVB_MITRE_KEEP")) {
        std::cerr << "\n*** MITRE REFUSED (MVB_MITRE_KEEP set -- geometry is NOT valid) ***\n"
                  << why << "\n*** the piece below keeps its growth and overshoots the joint ***\n"
                  << std::endl;
        return untrimmed;
    }
    throw std::runtime_error(why + " [set MVB_MITRE_KEEP=1 to export it anyway and look at it]");
}

// ABT #685 (Alf, 2026-08-18): " A rejected trim used to return the solid UNCHANGED and silently, so the piece kept the
// growth it was given for the mitre and visibly ran past the corner -- "the turn continues a bit
// further than the corner", on Primary parallel 0's turn 30 against its dragback. An untrimmed
// mitre is not a degraded result, it is copper in the wrong place, and in FEM that is fatal. It
// throws, naming the primitive and the junction, so the geometry that caused it gets fixed.
// A SPIRAL that closes a FULL REVOLUTION is the one piece whose own copper comes back beside
// its corner: its other end sits one pitch away along the column axis, and on a coil wound with
// touching turns that is exactly one wire radius from this end's surface. No knife can sever such
// a corner without biting the piece itself (measured: PQ33 pitch 0.959133 mm against 2r 0.959 mm,
// 0.13 um to spare; the single-switch forward is identical at 0.06 um). ABT #685 (Alf,
// 2026-08-18): the adjacency is an artefact of it being ONE SOLID -- a knife only ever subtracts
// from the solid it is applied to, so a pass belonging to a different primitive is untouchable by
// construction. Sweeping the revolution as two halves puts each half's other end 180 degrees
// away, where no knife reaches, and makes the pass beside the corner somebody else's solid.
// The halves overlap by half a wire radius of arc at a tangent seam and are fused back into ONE
// solid after both corners are trimmed, so a wrap still ships as a single named part.
static TopoDS_Shape localMitreTrim(const TopoDS_Shape& solid, const gp_Pnt& P, const gp_Dir& keepDir,
                                   const gp_Dir& stubDir, const gp_Pnt& stubEnd, double r,
                                   double grow, const std::string& what, double neighbourStep,
                                   const std::string& diagnosis) {
    // HOW FAR MAY EACH KNIFE AXIS TRAVEL? ABT #685 (Alf, 2026-08-18). The same wire's next pass
    // sits `neighbourStep` away ALONG THE COLUMN AXIS (one pitch per revolution -- and a
    // full-revolution wrap has its OWN other end exactly there), so its copper begins at
    // `neighbourStep - r`. A knife face may go that far and no further.
    //
    // The bound is PER AXIS, not one scalar. A single `reach` applied to every axis is wrong in
    // both directions: it clamped the DEPTH, which here runs almost perpendicular to the column
    // axis and never approaches the neighbour, while leaving the one axis that does point along
    // it (v) free to run to 1.6r. On the flyback's Primary that put the knife 0.511 mm into a gap
    // whose neighbouring copper starts at 0.334 mm: it gouged the turn one pitch below (Alf: "a
    // dent below the mitre, as if the mitre had cut into the turn below") and sheared off a
    // fragment that shipped as its own solid ("[2/2] ... shouldn't exist").
    const gp_XYZ kColumnAxis(0, 1, 0);
    const double neighbourGap = neighbourStep > 0.0 ? std::max(0.0, neighbourStep - r) : 0.0;
    auto axisLimit = [&](const gp_XYZ& e) {
        if (neighbourGap <= 0.0) return std::numeric_limits<double>::max();
        const double c = std::abs(e.Dot(kColumnAxis));
        return c < 1e-9 ? std::numeric_limits<double>::max() : neighbourGap / c;
    };
    // Oriented knife: the overhang stub extends along `stubDir` beyond the plane, so the box
    // follows its lateral DRIFT direction and stays tube-radius-thin on the perpendicular
    // (axial/pitch) axis -- a symmetric (3r+grow)-wide box nicked the neighbouring wrap pass
    // one pitch (~3r) away and detached a 0.37 mm fragment (measured, PQ33 junction 9->10).
    const gp_XYZ w = keepDir.XYZ() * -1.0;                              // into the discard side
    const double cosT = std::max(0.2, std::abs(stubDir.XYZ().Dot(w)));  // clamp: tilt <= ~78 deg
    const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    // Knife extent. It must cover the overhang and NOTHING else: a CURVED piece bends back
    // into the far part of the box and gets a rectangular notch gouged out of its keep side
    // (visible as square voids beside every dragback). Depth is therefore just the overhang plus
    // the tube radius -- what the tilted section actually needs -- and the perpendicular half is
    // barely wider than the tube.
    const double ell = r / cosT + 0.25 * r;   // mitre ellipse semi-extent + margin
    const double drift = (grow + r) * sinT;   // stub tip's lateral reach
    // LOCALITY BY SIZE, NOT BY INTERSECTION. ABT #685 (Alf, 2026-08-18): the knife used to be a
    // box intersected with a sphere centred on the joint, to stop its straight sides biting a
    // curved piece further along. The sphere was a disaster as a *cutting tool*: its face meets
    // the tube almost tangentially, and OCC's boolean returns a shape it then reports as INVALID
    // whose volume is garbage -- on 10_emi the cut of a 0.421 mm3 wrap measured 0.006 mm3, so the
    // guard read a 0.415 mm3 "removal" from a knife that was only 0.0165 mm3, refused the trim,
    // and left the piece standing proud. That refusal was the reported "mitre overshoot".
    // The box alone is local as long as its EXTENTS are: no face of it lies further than `reach`
    // from the joint, which is the same bound the sphere enforced, with flat faces the boolean can
    // actually handle.
    const double depth = std::min(grow + r, axisLimit(w));
    // The lateral faces must CLEAR the tube, never grave it. ABT #685 (Alf, 2026-08-18): at a
    // near-tangent joint the cut normal runs along the wire, so u and v are both lateral and the
    // knife becomes a sleeve around the tube -- and a half-extent of 1.05r puts its side face
    // 5.5 um outside a 0.11 mm cylinder. A plane grazing a cylinder is the classic boolean
    // degeneracy: OCC returned an INVALID cut whose volume read 0.006 mm3 for a 0.42 mm3 wrap.
    // A mitre cuts the WHOLE cross-section, so the box must contain the tube laterally with
    // margin; `reach` still caps it away from the neighbouring pass.
    // The knife must COVER THE WHOLE STUB in the lateral directions. ABT #685 (Alf, 2026-08-18):
    // clamping these by `reach` let the box's far face stop INSIDE the tube, so the cut severed
    // only part of the section and left a sliver -- OCC then returned an invalid shape (measured
    // on the push-pull: the stub needed 3.114 mm of drift, the clamp gave 2.529 mm, and a 418 mm3
    // wrap came back as 7.46 mm3). `reach` guards the direction the knife DIGS (depth, below),
    // which is where a neighbouring pass of the same wire sits; it has no business trimming the
    // knife across its own cut.
    const double lateral = 1.6 * r;
    gp_XYZ u = stubDir.XYZ() - w * stubDir.XYZ().Dot(w);
    if (u.Modulus() > 1e-9) u.Normalize();
    else u = gp_Ax2(P, gp_Dir(w)).XDirection().XYZ();   // stub ~normal to plane: no drift
    const gp_XYZ v = w.Crossed(u);
    const double uPos = std::min(std::max(ell + drift, lateral), axisLimit(u));
    const double uNeg = std::min(std::max(ell, lateral), axisLimit(u));
    const double vHalf = std::min(lateral, axisLimit(v));
    // A mitre SEVERS the section: the knife must span the tube on every lateral axis. If the
    // neighbouring pass is so close that it cannot, there is no knife that both cuts this corner
    // and spares the turn beside it -- say so rather than emit a partial cut (which leaves a
    // sliver and an invalid solid) or a gouged neighbour.
    // The knife needs real clearance outside the tube, not a tangent touch: a face lying ON the
    // cylinder is the grazing degeneracy that makes OCC return debris. Measured: the flyback's
    // Primary has 4.6% of a radius to spare (pitch 0.6536, 2r 0.639) and cuts cleanly, while the
    // single-switch forward's Primary is wound TOUCHING (pitch 0.383553 against 2r 0.3835), which
    // put the knife face 0.06 um outside the tube and returned an invalid shape that removed
    // nothing at all (faces 3 -> 9, volume unchanged). 2% of the radius separates the two.
    const double kSeverClearance = 1.02;
    if (uNeg < kSeverClearance * r || vHalf < kSeverClearance * r) {
        return refuseTrim(
            "ConductorBuilder: at " + what + " the same wire's next pass is only " +
                std::to_string(neighbourStep * 1e3) + " mm away (its copper starts at " +
                std::to_string(neighbourGap * 1e3) + " mm) against a " + std::to_string(r * 1e3) +
                " mm wire radius, so no knife can sever this corner without cutting into it "
                "(the turns are wound touching, and this wrap is a full revolution, so the pass "
                "beside the corner is this very piece's own other end); " +
                diagnosis,
            solid);
    }
    gp_Pnt corner(P.XYZ() - u * uNeg - v * vHalf);
    TopoDS_Shape box = BRepPrimAPI_MakeBox(gp_Ax2(corner, gp_Dir(w), gp_Dir(u)),
                                           uNeg + uPos, 2.0 * vHalf, depth).Shape();
    if (std::getenv("MVB_KNIFE_DIAG"))
        std::cerr << "[knife] " << what << " r=" << r * 1e3 << " grow=" << grow * 1e3
                  << " cosT=" << cosT << " ell=" << ell * 1e3 << " drift=" << drift * 1e3
                  << " uNeg=" << uNeg * 1e3 << " uPos=" << uPos * 1e3
                  << " vHalf=" << vHalf * 1e3 << " depth=" << depth * 1e3
                  << " nstep=" << neighbourStep * 1e3 << " gap=" << neighbourGap * 1e3
                  << " limU=" << axisLimit(u) * 1e3 << " limV=" << axisLimit(v) * 1e3
                  << " limW=" << axisLimit(w) * 1e3 << " mm" << std::endl;
    // ZERO REMOVAL = FAILED CUT (ABT #685, 2026-08-20). A mitred end always carries growth
    // (overS/overE >= mitreGrow > 0), so its knife must always remove material; OCC's Cut on
    // one steep spiral pipe (14_dab's steep-exit stub) returned the solid unchanged while
    // reporting Done+valid -- removed -9e-12 mm3 with the overhang tip verified inside the
    // box -- and the untrimmed wedge shipped overlapping the exit lead by 0.148 mm3. When the
    // local box removes nothing and the piece has no neighbouring pass to protect (the only
    // reason the knife is local), retry once with a half-space-scale knife; a cut that still
    // removes nothing is refused loudly.
    const bool canRetryHalfSpace = neighbourStep <= 0.0;
    TopoDS_Shape bigBox;
    if (canRetryHalfSpace) {
        const double big = 20.0 * (r + grow);
        gp_Pnt bigCorner(P.XYZ() - u * big - v * big);
        bigBox = BRepPrimAPI_MakeBox(gp_Ax2(bigCorner, gp_Dir(w), gp_Dir(u)),
                                     2.0 * big, 2.0 * big, big).Shape();
    }
    try {
        // AN INVALID INPUT CANNOT BE CUT: OCC's booleans return an invalid argument unchanged
        // while reporting Done (measured, 14_dab's steep-exit stub: a 2.5-degree steep spiral
        // pipe whose sweep comes out invalid -- the local knife AND a half-space both "removed"
        // -9e-12 mm3). Repair the piece first; if it will not repair, the refusal below names it.
        TopoDS_Shape solidIn = solid;
        if (!BRepCheck_Analyzer(solidIn).IsValid()) {
            if (std::getenv("MVB_MITRE_DIAG"))
                std::cerr << "[trim-input-invalid] " << what << " -- repairing before the cut"
                          << std::endl;
            try {
                ShapeFix_Shape fixIn(solidIn);
                fixIn.Perform();
                if (!fixIn.Shape().IsNull() && BRepCheck_Analyzer(fixIn.Shape()).IsValid())
                    solidIn = fixIn.Shape();
            } catch (const Standard_Failure&) {
            }
        }
        GProp_GProps gpBefore, gpAfter;
        BRepGProp::VolumeProperties(solidIn, gpBefore);
        // CUT IN THE MILLIMETRE FRAME. ABT #685 (Alf, 2026-08-18): rawGrownSolid already sweeps
        // spiral/blend pipes at x1000 because OCC mis-frames sub-millimetre profiles in metres --
        // the boolean has exactly the same problem one step later. Measured on the push-pull's
        // 90-degree layer-link corner, in METRES: Cut(wrap, knife) returned a 7.46 mm3 fragment of
        // a 418 mm3 wrap, and Common(wrap, knife) returned ZERO -- OCC could not even see that the
        // two solids touch, though the joint sits exactly on the wrap's centreline. Scaling is an
        // EXACT affine map both ways, so this is a numerical frame choice, not an approximation.
        const TopoDS_Shape solidMm = toMitreFrame(solidIn);
        const TopoDS_Shape boxMm = toMitreFrame(box);
        BRepAlgoAPI_Cut cut(solidMm, boxMm);
        if (cut.IsDone() && !cut.Shape().IsNull()) {
            TopoDS_Shape trimmed = fromMitreFrame(cut.Shape());
            // REPAIR BEFORE JUDGING. ABT #685 (Alf, 2026-08-18): cutting a swept-pipe wrap leaves a
            // shell OCC reports as invalid -- and the volume of an invalid shell is meaningless
            // (the divergence integral needs a closed one). Measured on the push-pull: a 418 mm3
            // wrap cut by a 31 mm3 knife came back reading 7.46 mm3, so the guard below saw a
            // 411 mm3 "removal" and refused a trim that had actually done the right thing (faces
            // 5 -> 6: the mitre face WAS added). ShapeFix is the same repair the assembler already
            // runs on every trimmed piece for the STEP round-trip; running it here, before the
            // volume is read, means the guard judges the repaired solid instead of the debris.
            const bool rawValid = BRepCheck_Analyzer(trimmed).IsValid();
            double rawVol = -1;
            { GProp_GProps g; BRepGProp::VolumeProperties(trimmed, g); rawVol = g.Mass(); }
            if (!BRepCheck_Analyzer(trimmed).IsValid()) {
                try {
                    ShapeFix_Shape fix(trimmed);
                    fix.Perform();
                    if (!fix.Shape().IsNull() && BRepCheck_Analyzer(fix.Shape()).IsValid())
                        trimmed = fix.Shape();
                } catch (const Standard_Failure&) {
                }
            }
            // A shape OCC reports as invalid has a MEANINGLESS volume, so the guard below cannot
            // be trusted on one -- that is exactly how a good trim got read as a 98.5% removal.
            if (!BRepCheck_Analyzer(trimmed).IsValid()) {
                if (std::getenv("MVB_MITRE_DIAG")) {
                    GProp_GProps gk, ga;
                    BRepGProp::VolumeProperties(box, gk);
                    BRepGProp::VolumeProperties(trimmed, ga);
                    int nf = 0, nsol = 0;
                    for (TopExp_Explorer e(trimmed, TopAbs_FACE); e.More(); e.Next()) ++nf;
                    for (TopExp_Explorer e(trimmed, TopAbs_SOLID); e.More(); e.Next()) ++nsol;
                    int nfIn = 0;
                    for (TopExp_Explorer e(solid, TopAbs_FACE); e.More(); e.Next()) ++nfIn;
                    std::cerr << "[mitre-invalid] " << what << "\n   r=" << r * 1e3
                              << " grow=" << grow * 1e3 << " nstep=" << neighbourStep * 1e3
                              << " depth=" << depth * 1e3 << " uNeg=" << uNeg * 1e3
                              << " uPos=" << uPos * 1e3 << " vHalf=" << vHalf * 1e3
                              << " mm\n   knife=" << gk.Mass() * 1e9 << " before="
                              << gpBefore.Mass() * 1e9 << " after=" << ga.Mass() * 1e9
                              << " mm3, faces " << nfIn << "->" << nf << ", solids " << nsol
                              << "\n   cosT=" << cosT << " sinT=" << sinT << " drift="
                              << drift * 1e3 << " ell=" << ell * 1e3 << " mm" << std::endl;
                }
                return refuseTrim("ConductorBuilder: the mitre cut at " + what +
                                      " produced a shape OCC reports as INVALID, so neither the"
                                      " geometry nor its volume can be trusted; " + diagnosis,
                                  solid);
            }
            BRepGProp::VolumeProperties(trimmed, gpAfter);
            double removed = gpBefore.Mass() - gpAfter.Mass();
            if (removed < 1e-15 && canRetryHalfSpace && !bigBox.IsNull()) {
                // The local knife cut nothing off a grown end: retry with the half-space knife.
                const TopoDS_Shape bigMm = toMitreFrame(bigBox);
                BRepAlgoAPI_Cut cut2(solidMm, bigMm);
                if (cut2.IsDone() && !cut2.Shape().IsNull()) {
                    TopoDS_Shape trimmed2 = fromMitreFrame(cut2.Shape());
                    if (!BRepCheck_Analyzer(trimmed2).IsValid()) {
                        try {
                            ShapeFix_Shape fix2(trimmed2);
                            fix2.Perform();
                            if (!fix2.Shape().IsNull() &&
                                BRepCheck_Analyzer(fix2.Shape()).IsValid())
                                trimmed2 = fix2.Shape();
                        } catch (const Standard_Failure&) {
                        }
                    }
                    if (BRepCheck_Analyzer(trimmed2).IsValid()) {
                        GProp_GProps g2;
                        BRepGProp::VolumeProperties(trimmed2, g2);
                        const double removed2 = gpBefore.Mass() - g2.Mass();
                        if (removed2 > 1e-15) {
                            if (std::getenv("MVB_MITRE_DIAG"))
                                std::cerr << "[trim-retry] " << what
                                          << " half-space knife removed "
                                          << removed2 * 1e9 << " mm3 after the local box "
                                          << "removed nothing" << std::endl;
                            trimmed = trimmed2;
                            gpAfter = g2;
                            removed = removed2;
                        }
                    }
                }
            }
            // THE ONLY VERDICT THAT IS NOT THE BOOLEAN'S OWN. ABT #860 (2026-08-23): both of the
            // old verdicts came from the kernel that does the cutting -- `removed` from the Cut and
            // a Common against the knife -- so a kernel that goes blind passes both while the
            // overhang still stands (it did, on every --segments 0 terminal stub). Ask the trimmed
            // solid directly whether copper of its own is still standing past the bisector plane.
            const double leftBeyond =
                overhangBeyondPlane(trimmed, stubEnd, stubDir, P, w, grow);
            if (leftBeyond > 0.0) {
                std::ostringstream why2;
                why2 << "ConductorBuilder: the mitre cut at " << what << " left the overhang"
                     << " standing -- a point on the piece's own centreline " << leftBeyond * 1e3
                     << " mm past the bisector plane is still INSIDE the trimmed solid (the cut"
                     << " reported success and removed " << removed * 1e9 << " mm3 of "
                     << gpBefore.Mass() * 1e9 << " mm3; growth " << grow * 1e3 << " mm, wire r "
                     << r * 1e3 << " mm, half-space retry "
                     << (canRetryHalfSpace ? "available" : "unavailable") << "); " << diagnosis;
                return refuseTrim(why2.str(), solid);
            }
            if (removed < 1e-15 && std::getenv("MVB_MITRE_DIAG"))
                std::cerr << "[trim-zero-ok] " << what
                          << " grown end reaches but does not cross the plane" << std::endl;
            // stub bound: the overhang is at most a full-radius tube of length grow+2r, doubled
            // for the tilted-ellipse wedge. More than that = the knife ate distant material.
            const double stubMax = 2.0 * kPi * r * r * (grow + 2.0 * r);
            if (removed <= stubMax && gpAfter.Mass() > 0.0) {
                if (std::getenv("MVB_MITRE_DIAG")) {
                    std::cerr << "[trim] " << what << " removed=" << removed * 1e9
                              << " mm3 (grow=" << grow * 1e3 << " mm)" << std::endl;
                    if (removed < 1e-15) {
                        Bnd_Box bs;
                        BRepBndLib::Add(solid, bs);
                        double sx0, sy0, sz0, sx1, sy1, sz1;
                        bs.Get(sx0, sy0, sz0, sx1, sy1, sz1);
                        std::cerr << "[trim-miss] P=(" << P.X() * 1e3 << "," << P.Y() * 1e3
                                  << "," << P.Z() * 1e3 << ") w=(" << w.X() << "," << w.Y()
                                  << "," << w.Z() << ") stubDir=(" << stubDir.X() << ","
                                  << stubDir.Y() << "," << stubDir.Z() << ") depth="
                                  << depth * 1e3 << " uPos=" << uPos * 1e3 << " uNeg="
                                  << uNeg * 1e3 << " vHalf=" << vHalf * 1e3
                                  << " solidBB=[" << sx0 * 1e3 << "," << sy0 * 1e3 << ","
                                  << sz0 * 1e3 << "]..[" << sx1 * 1e3 << "," << sy1 * 1e3
                                  << "," << sz1 * 1e3 << "]" << std::endl;
                    }
                    int ns = 0;
                    for (TopExp_Explorer e(trimmed, TopAbs_SOLID); e.More(); e.Next()) ++ns;
                    if (ns != 1)
                        std::cerr << "[mitre-split] " << what << " -> " << ns << " solids"
                                  << ", removed " << removed * 1e9 << " mm3 of "
                                  << gpBefore.Mass() * 1e9 << ", r=" << r * 1e3
                                  << " grow=" << grow * 1e3 << " nstep=" << neighbourStep * 1e3
                                  << " depth=" << depth * 1e3 << " uNeg=" << uNeg * 1e3
                                  << " uPos=" << uPos * 1e3 << " vHalf=" << vHalf * 1e3
                                  << " mm; " << diagnosis << std::endl;
                }
                return trimmed;
            }
            if (std::getenv("MVB_MITRE_DIAG")) {
                std::cerr << "[stages] raw cut valid=" << rawValid << " vol=" << rawVol * 1e9
                          << " mm3 -> after fix vol=" << gpAfter.Mass() * 1e9
                          << " mm3 (input " << gpBefore.Mass() * 1e9 << " mm3)" << std::endl;
            }
            std::ostringstream why;
            why << "ConductorBuilder: the mitre knife at " << what << " REMOVED "
                << removed * 1e9 << " mm3 against a " << stubMax * 1e9
                << " mm3 bound for a corner stub"
                << " (before " << gpBefore.Mass() * 1e9 << " mm3, after "
                << gpAfter.Mass() * 1e9 << " mm3)"
                << "; knife at (" << P.X() * 1e3 << "," << P.Y() * 1e3 << "," << P.Z() * 1e3
                << ") mm, keep dir (" << keepDir.X() << "," << keepDir.Y() << "," << keepDir.Z()
                << "), wire r=" << r * 1e3 << " mm, growth=" << grow * 1e3 << " mm, neighbour step="
                << neighbourStep * 1e3 << " mm, " << diagnosis << "."
                << " The knife reached material away from the joint, so the trim is refused and"
                << " the piece would stand proud of the corner by its growth.";
            return refuseTrim(why.str(), solid);
        }
    } catch (const Standard_Failure& f) {
        return refuseTrim("ConductorBuilder: the mitre knife at " + what + " threw in OCC (" +
                              std::string(f.GetMessageString() ? f.GetMessageString() : "(null)") +
                              "); knife at (" + std::to_string(P.X() * 1e3) + "," +
                              std::to_string(P.Y() * 1e3) + "," + std::to_string(P.Z() * 1e3) +
                              ") mm, wire r=" + std::to_string(r * 1e3) + " mm, growth=" +
                              std::to_string(grow * 1e3) + " mm",
                          solid);
    }
    return refuseTrim("ConductorBuilder: the mitre cut at " + what + " did not complete (OCC "
                      "reported not-done or an empty result), so the piece would stand proud of "
                      "the corner by its " + std::to_string(grow * 1e3) + " mm growth",
                      solid);
}

// EXACT-ROUND SHARED CAP (segments <= 0). ABT #1186. At a TANGENT junction the two pieces' end
// caps are two INDEPENDENTLY constructed circles on the same plane -- one disc geometrically,
// two faces topologically -- and that is the worst possible input for BOPAlgo_PaveFiller's
// edge/face stage. Measured on 13_current_sense_er95_n87 at segments = 0: the abutment gate's
// two Common booleans cost a mean 6.42 s on each of the 134 ARC3->ARC3 biarc-riser joints, 64%
// of an 833 s build, for a lens that cannot exist.
//
// The faceted path already exempts exactly this case (capsProvablyShared below), and its clause
// is gated on segments > 0 -- correctly, because ITS proof is the polygon handoff, and
// startProfileOverride is a PHASE continuation (sectionPhaseOf / wireProfileWirePhased) that
// only means something on an n-gon. The guard is right about the mechanism and silent about the
// conclusion: at segments <= 0 the section is wireProfileWireSplit(point, dir, r, 2), a circle
// fully determined by the junction point, the tangent and the wire radius. There is nothing to
// hand over because two pieces meeting at one point with one tangent already build THE SAME
// DISC -- and therefore nothing for the boolean to find.
//
// That is not assumed here, it is MEASURED on the solids as built. capDiscAt reads the planar
// boundary faces a solid actually presents in the junction plane and requires them to add up to
// the whole nominal disc: same plane to 0.1 nm, total area pi*r^2 to 1e-6 relative, area-
// weighted centroid on the junction point to 1e-6 r. An exact revolve or prism passes; a sweep
// whose cap landed where MakePipeShell transported it does not (the [pipe-cap] probe has
// measured 0.13 r off nominal), and keeps its booleans. So the exemption follows the geometry,
// never the primitive kind.
//
// The residual is bounded, not hand-waved: two caps whose planes are within 1e-10 m can enclose
// at most pi*r^2 * 2e-10 m^3 -- 6e-7 mm^3 even for a 1 mm wire radius, below the gate's own
// 1e-6 mm^3 floor, and 1.6e-9 mm^3 on this design's 50 um wire. What is NOT re-proved is a
// piece looping back onto its immediate neighbour far from the shared cap; the faceted
// exemption does not test for that either, and the per-solid interpenetration sweep
// (MVB_MITRE_INTERPEN) is where that question lives.
static bool capDiscAt(const TopoDS_Shape& s, const gp_Pnt& j, const gp_Dir& n, double r,
                      double* areaOut, double* planeOffOut) {
    if (areaOut) *areaOut = 0.0;
    if (planeOffOut) *planeOffOut = -1.0;
    if (s.IsNull() || !(r > 0.0)) return false;
    const double kCapPlaneTol = 1e-10;   // metres between the two cap planes
    const double kCapNormalTol = 1e-9;   // radians
    double area = 0.0, worstOff = 0.0;
    gp_XYZ moment(0.0, 0.0, 0.0);
    int nf = 0;
    for (TopExp_Explorer fx(s, TopAbs_FACE); fx.More(); fx.Next()) {
        const TopoDS_Face f = TopoDS::Face(fx.Current());
        BRepAdaptor_Surface sf(f);
        if (sf.GetType() != GeomAbs_Plane) continue;
        const gp_Pln pln = sf.Plane();
        const gp_Dir pn = pln.Axis().Direction();
        if (std::min(pn.Angle(n), pn.Angle(n.Reversed())) > kCapNormalTol) continue;
        const double off = pln.Distance(j);
        if (off > kCapPlaneTol) continue;
        GProp_GProps g;
        BRepGProp::SurfaceProperties(f, g);
        const double a = g.Mass();
        if (!(a > 0.0)) continue;
        area += a;
        moment += g.CentreOfMass().XYZ() * a;
        worstOff = std::max(worstOff, off);
        ++nf;
    }
    if (areaOut) *areaOut = area;
    if (planeOffOut) *planeOffOut = nf ? worstOff : -1.0;
    if (nf == 0) return false;
    const double want = kPi * r * r;
    if (std::abs(area - want) > 1e-6 * want) return false;   // not the whole section
    return gp_Pnt(moment / area).Distance(j) <= 1e-6 * r;    // and centred on the junction
}

TopoDS_Shape assembleWire(const std::vector<const Primitive*>& ptrs, double wireRadius,
                          int segments, CornerStyle corners,
                          std::vector<size_t>* primIndexPerSolid) {
    // WHEN IS A JUNCTION A CORNER? ABT #685 (Alf, 2026-08-18). Not "below 3 degrees", which was
    // another chosen number. A junction that is NOT mitred is BRIDGED: the earlier piece grows
    // flush past the joint until it fills the wedge the direction change opens on the outer side
    // of the bend, r*tan(theta), and the union swallows it. The only cost of bridging is that the
    // bridge stub pokes out of the neighbouring tube, by r*tan(theta)*sin(theta).
    // So mitre exactly when that poke-out is something this model can even REPRESENT -- i.e. when
    // it exceeds the sweep's own chordal sag, the error the faceted tube already carries
    // (samplingSag). Below that the mitre is beneath the model's resolution, and bridging is
    // strictly better: no boolean to fail, and no wedge gap either.
    // This matters because the boolean DOES fail. On 10_emi a 3.47 deg riser/wrap joint (poke-out
    // 0.4 um against a 2.2 um sag) was mitred, and OCC could not cut that spiral pipe surface at
    // ANY slice thickness -- 143 um or 3.3 um both came back INVALID, with the volume collapsing
    // from 0.411 mm3 to 0.0008 mm3 and a face LOST. The trim was refused, the piece stood proud,
    // and that was the reported "mitre overshoot".
    const double sag = samplingSag(wireRadius);
    // Growth that fills the wedge at a bridged (un-mitred) joint. Clamped well inside the range
    // where bridging is ever chosen, so tan() cannot run away.
    auto bridgeGrow = [&](double ang) { return wireRadius * std::tan(std::min(ang, 0.5)); };
    // A bridged joint fills its wedge by growing the earlier piece PAST the joint, so it pokes
    // into the neighbour: real copper inside copper. Alf, 2026-08-25: "no overlap must be
    // allowed in copper copper". The sag test alone is no longer sufficient -- it asks only
    // whether the mitre is BELOW THE MODEL'S RESOLUTION, not whether the overlap it leaves is
    // acceptable, and it is not (measured 0.0031 mm3, 1.26% of the smaller piece, at a 1.44 deg
    // SEG->SEG lead joint on 16_coupled).
    //
    // Mitre whenever the joint actually bends, EXCEPT where the cut is known to be unsafe. The
    // documented failure (10_emi, 3.47 deg) was OCC refusing to slice a SPIRAL PIPE surface --
    // a swept BSpline -- at a shallow angle. Analytic SEG/ARC3 solids slice cleanly at the same
    // angle, so the exemption belongs to the PIPE KINDS, not to the angle.
    // MVB_MITRE_PIPES=1: treat SPIRAL/BLEND pipe ends as cuttable too, i.e. mitre EVERY bend.
    // The pipe exemption rests on a 2026-08 measurement (10_emi, 3.47 deg) taken before the
    // polygon-profile rebuild, the ShapeFix pass and the flush-tube repair fallback existed, so
    // it is worth re-measuring rather than assuming. Any slice that still fails is caught by the
    // existing validity check and rebuilt as a flush tube -- loudly, never silently.
    const bool mitrePipes = std::getenv("MVB_MITRE_PIPES") != nullptr;
    auto cuttableKind = [mitrePipes](int k) {
        return mitrePipes || k == Primitive::SEG || k == Primitive::ARC3;
    };
    auto worthMitring = [&](double ang, int ka, int kb) {
        if (ang <= 1e-9) return false;                       // truly tangent: caps already coincide
        if (cuttableKind(ka) && cuttableKind(kb)) return true;
        // PIPE END (ABT #961). The sag rule measured the WRONG quantity. It compared the bridged
        // stub's lateral poke-out, r*tan(ang)*sin(ang) -- quadratic in the angle -- against the
        // model's chordal sag, and so called a 3.47 deg riser/wrap joint "beneath resolution"
        // (0.37 um on 10_emi) and bridged it. But a bridge grows the earlier piece r*tan(ang)
        // = 6 um straight into its neighbour: a full-section overlap of 0.00019 mm3, forty times
        // per design, which the weld then has to hide and cannot when the fused body
        // self-intersects. What a flush, un-mitred joint actually leaves is a WEDGE: each
        // piece's perpendicular cap misses the shared bisector by r*tan(ang/2) at its rim,
        // LINEAR in the angle -- 3.0 um at 3.47 deg, above the 2.2 um sag, so the wedge is real
        // and the joint must be mitred (measured: NO OVERLAPS, watertight, same copper to
        // 0.004 mm3). At 2.8e-5 deg (00_debug's closing wrap) the same quantity is 2.5e-8 mm,
        // nothing the model can represent, and mitring THERE is what let a half-space knife
        // bite a closed revolution's other end (REMOVED 89.2 mm3 against a 27.2 mm3 bound) --
        // so those joints stay bridged, with a growth of a fraction of a nanometre.
        return wireRadius * std::tan(0.5 * ang) > sag;
    };
    // |unit + unit| = 2 cos(angle/2): 0.2 admits joints up to ~168 degrees and rejects the folds
    // beyond, whose bisector carries no usable direction (ABT #685).
    const double kMinBisector = 0.2;
    // GROWTH PER MITRED END -- derived, not chosen. ABT #685 (Alf, 2026-08-18: "why 1.3? that
    // looks a magic number" -- it was). A bisector plane cuts a tube of radius r whose axis meets
    // the plane NORMAL at half the joint angle, so the ellipse it carves reaches exactly
    // r*tan(theta/2) past the joint along the axis. Grow that much and the plane is guaranteed to
    // find material everywhere it cuts, with nothing left over to slice away.
    // The old fixed 1.3*wireRadius is tan(52 deg): sized for a ~105 deg worst case and then
    // applied to every junction. It was wrong at BOTH ends -- on 10_emi's 3.47 deg riser/wrap
    // joint the mitre needs 3.3 um and it grew 143 um, a 43x stub that OCC then could not slice
    // off a spiral pipe surface (the cut came back INVALID, so the trim was refused and the piece
    // stood proud: the reported "mitre overshoot"); and past ~105 deg it grew LESS than the plane
    // cuts, leaving the mitre face notched.
    // The cap is not a new constant -- it is the SAME conditioning limit as kMinBisector:
    // |a+b| = 2cos(theta/2) > 0.2 admits theta < 168.5 deg, so the growth is bounded by ~10r.
    const double kMaxMitreAngle = 2.0 * std::acos(0.5 * kMinBisector);
    auto mitreGrow = [&](double ang) {
        return wireRadius * std::tan(0.5 * std::min(ang, kMaxMitreAngle));
    };
    // A PIECE GROWN ALONG ITS OWN CURVE FALLS SHORT OF ITS OWN TANGENT (2026-08-31,
    // 20_iso_buckboost_epc25_3c91). r*tan(theta/2) is the exact overhang a STRAIGHT piece needs
    // for the extreme point of its tilted section to cross the bisector plane. But an ARC3
    // (revolved: ddA/ddB = over/radius) and an analytically grown SPIRAL travel that distance
    // ON THE ARC, not on the tangent, so the grown tip lands R*(1 - cos(g/R)) short of where
    // the straight formula puts it. The mitre plane then cuts a face that is part plane and
    // part TUBE -- non-planar, with an extra vertex where the two surfaces meet -- and OCC's
    // glue can no longer identify it with the neighbouring prism's planar mitre face, so the
    // fragment FUSES the two bodies instead of imprinting a shared one. Measured on 20's exit
    // lead corner: the stub's knifed face carried a 4.54 um dent on two of its 17 vertices
    // (against 0.00 um on all 16 of the mirror-image entrance corner, whose knife extents are
    // bit-identical: r=0.1775 grow=0.1775 ell=0.295398 uNeg=0.295398 uPos=0.546421 vHalf=0.284
    // depth=0.355 mm, nothing clamped on either side). The fragment then returned ONE solid of
    // 3.42834 mm3 instead of two (3.42262 + 0.29981), and gmsh reported "3D Meshing 2 volumes
    // with 1 connected component / Found volume 26" -> "No elements in volume 27". The
    // entrance corner's fragment produced the shared PLANE face (0.1351735 mm2 on both solids)
    // and meshed.
    // The correction is the arc's own departure from its tangent, R*(1 - cos(g/R)) -- exact for
    // the ARC3 revolve, and for a SPIRAL the wrap radius at the grown end (a helix's radius of
    // curvature is R + b^2/R >= R, so using R over-grows rather than under-grows). It is never
    // shipped copper: everything past the bisector is what the knife removes.
    auto curveGrow = [&](const Primitive& pr, double g, bool atStart) {
        if (!(g > 0.0)) return 0.0;
        double R = 0.0;
        if (pr.kind == Primitive::ARC3) R = pr.arc.v0.Modulus();
        else if (pr.kind == Primitive::SPIRAL && primEdgeGrowsAnalytically(pr))
            R = atStart ? pr.spiral.r0 : pr.spiral.r1;
        if (!(R > 1e-12)) return 0.0;   // straight, or grown along a straight spine extension
        return R * (1.0 - std::cos(std::min(g / R, kPi)));
    };
    const bool diag = std::getenv("MVB_MITRE_DIAG") != nullptr;
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    size_t n = ptrs.size();
    if (n == 0) return compound;
    std::vector<gp_Dir> fs(n), fe(n);
    for (size_t i = 0; i < n; ++i) {
        fs[i] = primFwdStart(*ptrs[i], wireRadius);
        fe[i] = primFwdEnd(*ptrs[i], wireRadius);
    }
    int nCut = 0, nRepaired = 0, nInvalid = 0;
    std::vector<TopoDS_Shape> built;
    // MITRED-JOINT ABUTMENT CHECK (ABT #685, Alf 2026-08-20: "there is no real bisection, just
    // overlapping solids ... learn how to detect those and add that to the checkings"). A
    // mitred corner slices BOTH sides on one shared bisector plane, so the two solids must
    // ABUT: their boolean COMMON is empty. Any residual volume is copper counted twice -- not
    // a valid FEM assembly and not a real winding. Checked for every CORNER joint as it is
    // built; bridged (near-tangent) joints are exempt, their wedge fill overlaps by design.
    // The floor is a (10 um)^3 cube in the mm frame: orders of magnitude above the boolean's
    // face-coincidence noise and orders below any real un-bisected corner.
    TopoDS_Shape prevBuilt;
    // Pieces grouped per primitive + whether each primitive is BRIDGED onto its predecessor.
    std::vector<std::vector<TopoDS_Shape>> perPrimSolids;
    std::vector<bool> bridgedToPrev;
    // A joint sphere spans the junction it fills: half of it lies in the piece it is fused into
    // and half in the NEXT piece. Those two must therefore end up as ONE solid, or the sphere
    // simply becomes copper overlapping its own neighbour -- which is exactly what chunk 7
    // measured before this (28 self-overlaps against 19 without it, with weldfail=0, i.e. the
    // weld never even attempted those junctions). Set when a sphere is emitted at the junction
    // leaving piece i; consumed when piece i+1 decides whether it is welded onto its predecessor.
    bool sphereAtPrevJunction = false;
    std::vector<std::string> overlappingJoints;
    // ABT #685: a silent check is worth nothing -- "no overlap reported" must not be able to mean
    // "the boolean never answered". Count the joints that got a real verdict against the ones that
    // did not, so a clean corpus can be told apart from a blind one (MVB_MITRE_DIAG prints both).
    size_t abutChecked = 0, abutNoVerdict = 0;
    // MEASURE ONLY WHAT THE KERNEL CAN VALIDLY PRODUCE, AND TAKE THE WORST OF IT. ABT #860
    // (2026-08-23): this gate used to read the volume of whatever BRepAlgoAPI_Common handed back,
    // in the millimetre frame, without asking whether the result was a valid solid -- the one thing
    // the rest of this file already knows never to do ("the volume of an invalid shell is
    // meaningless"). Measured on 13_current_sense at --segments 0: the Common of two VALID solids
    // came back INVALID reading 0.000973882 mm3 -- 87% of the 0.001118 mm3 dragback stub, i.e. the
    // stub reported as very nearly buried inside its neighbour -- and the design was refused for a
    // corner that does not overlap at all. The SAME two solids intersect in exactly 0, valid, when
    // the boolean runs in metres (and after a BRep round trip, and at any fuzzy value >= 1e-5 mm);
    // no point of the stub classifies inside the neighbour; and the faceted build of that very
    // design reports the joint clean. So: run the boolean in BOTH frames, discard any answer OCC
    // reports as invalid, and take the LARGEST of the trustworthy ones. That is conservative in the
    // direction a gate must be conservative -- a valid answer that shows overlap is never
    // discarded, and the frame that goes blind (also measured, see toMitreFrame) cannot argue a
    // corner clean on its own. If NEITHER frame gives a valid answer the joint is counted as
    // unchecked, loudly, and never as clean.
    auto checkMitredAbutment = [&](const TopoDS_Shape& a, const TopoDS_Shape& b,
                                   const std::string& joint) {
        if (a.IsNull() || b.IsNull()) { ++abutNoVerdict; return; }
        double worst = -1.0;                 // mm^3, over the frames that answered validly
        Bnd_Box worstBox;
        // scale: how many mm^3 one cubic unit of that frame's result is worth
        const std::array<double, 2> frameScale{1.0, 1e9};
        for (int frame = 0; frame < 2; ++frame) {
            try {
                const TopoDS_Shape x = frame == 0 ? toMitreFrame(a) : a;
                const TopoDS_Shape y = frame == 0 ? toMitreFrame(b) : b;
                BRepAlgoAPI_Common common(x, y);
                if (!common.IsDone() || common.Shape().IsNull()) continue;
                if (!BRepCheck_Analyzer(common.Shape()).IsValid()) continue;   // debris, not a volume
                GProp_GProps gcp;
                BRepGProp::VolumeProperties(common.Shape(), gcp);
                const double mm3 = gcp.Mass() * frameScale[frame];
                if (mm3 > worst) {
                    worst = mm3;
                    Bnd_Box cb;
                    BRepBndLib::Add(common.Shape(), cb);
                    if (frame == 1 && !cb.IsVoid()) {   // report every box in mm
                        double bx0, by0, bz0, bx1, by1, bz1;
                        cb.Get(bx0, by0, bz0, bx1, by1, bz1);
                        Bnd_Box scaled;
                        scaled.Update(bx0 * 1e3, by0 * 1e3, bz0 * 1e3,
                                      bx1 * 1e3, by1 * 1e3, bz1 * 1e3);
                        cb = scaled;
                    }
                    worstBox = cb;
                }
            } catch (const Standard_Failure&) {
                // The boolean itself failing is a different defect; the invalid-cut paths
                // already throw where the trim happens. Let the other frame answer.
            }
        }
        if (worst < 0.0) { ++abutNoVerdict; return; }   // no frame produced a valid answer
        ++abutChecked;
        if (worst > 1e-6) {
            std::ostringstream where;
            where << joint << " overlaps by " << worst << " mm^3";
            if (!worstBox.IsVoid()) {
                double cx0, cy0, cz0, cx1, cy1, cz1;
                worstBox.Get(cx0, cy0, cz0, cx1, cy1, cz1);
                where << ", common bbox mm [" << cx0 << "," << cy0 << "," << cz0 << "]..["
                      << cx1 << "," << cy1 << "," << cz1 << "]";
            }
            overlappingJoints.push_back(where.str());
        }
    };
    // The previous SEG prism's exact end-cap vertices, handed to the next prism so mitre
    // junction caps are BIT-IDENTICAL (see mitredFacetPrism). Valid only across a directly
    // consecutive prism pair; anything else (sweep piece, refused prism) breaks the chain.
    std::vector<gp_Pnt> prevEndCap;
    bool prevEndCapValid = false;
    // Whether prevEndCap came from a pipe/revolve's ACTUAL end section (trusted geometric
    // adoption) rather than a prism's deterministic-frame projection (identity adoption).
    bool prevEndCapTrusted = false;
    // Whether the piece that PRODUCED prevEndCap is straight at its end (a SEG). A straight
    // source's facets stay coplanar with an aligned receiver indefinitely (osculation, the
    // 04_forward killer); a curved source separates quadratically within a facet length and
    // NEEDS the exact shared face (10_emi's ride stack). See the stagger rule in
    // mitredFacetPrism.
    bool prevEndCapSourceStraight = false;
    // The straight source's LENGTH decides the stagger (two measured populations, a decade
    // apart): 04_forward's 12 mm runs (48 r) produce coplanar overlap bands no local mesh can
    // resolve -- stagger required (verbatim kept it at 5/418; full stagger 418/418) -- while
    // 10_emi's 0.9 mm risers (3.6 r) need the verbatim shared face (stagger broke them to
    // 5/50). The boundary is one section circumference, 2*pi*r: a coplanar band longer than
    // the full perimeter cannot cross transversally anywhere.
    double prevSourceLen = 0.0;
    for (size_t i = 0; i < n; ++i) {
        // Bend at each end: a tangent joint (wrap SEG<->ARC) needs no growth and no boolean; only a
        // real corner (the entrance/exit leads) is grown and sliced on its angle-bisector plane.
        // Junctions may carry an ENDPOINT MISMATCH (layouts hand out lead ends offset from the wrap
        // start by up to a wire radius -- measured 0.36 mm on PQ33): the mitre plane goes through
        // the MIDPOINT of the two endpoints (same plane for both sides) and the growth is widened
        // by the mismatch so the overhang always reaches past the shared plane.
        const double dpS = (i > 0)
            ? primEndpoints(*ptrs[i - 1]).second.Distance(primEndpoints(*ptrs[i]).first) : 0.0;
        const double dpE = (i + 1 < n)
            ? primEndpoints(*ptrs[i]).second.Distance(primEndpoints(*ptrs[i + 1]).first) : 0.0;
        bool bentS = false;
        gp_Dir nS = fs[i];
        if (diag && i > 0) {
            static const char* kn[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
            std::cerr << "[slope]   " << ptrs[i - 1]->label << "  dy/ds=" << fe[i - 1].Y()
                      << "   |   " << ptrs[i]->label << "  dy/ds=" << fs[i].Y() << "\n";
            std::cerr << "[mitre]   junction " << i - 1 << "->" << i << " " << kn[ptrs[i - 1]->kind]
                      << "->" << kn[ptrs[i]->kind] << " angle=" << fe[i - 1].Angle(fs[i]) * 180.0 / kPi
                      << " deg" << (worthMitring(fe[i - 1].Angle(fs[i]), ptrs[i - 1]->kind, ptrs[i]->kind) ? " CORNER" : " bridged")
                      << " ['" << ptrs[i - 1]->label << "' -> '" << ptrs[i]->label << "']\n";
        }
        // ABT #685 (Alf, 2026-08-14): a lead/turn corner is mitred exactly like the elbow between
        // the two lead runs themselves (the vertical bridge into the radial run out) — one bisector
        // plane, both sides grown to it and sliced, so on a 90 degree corner the cut is at 45 and
        // the two pieces share the face. Tried and rejected on the way here: slicing on the LEAD's
        // face, and leaving both flush with a wire-radius overshoot — the first is a long slanted
        // ellipse, the second leaves the corner's outer side notched.
        if (i > 0 && worthMitring(fe[i - 1].Angle(fs[i]), ptrs[i - 1]->kind, ptrs[i]->kind)) {
            gp_Vec s(fe[i - 1].XYZ());
            s += gp_Vec(fs[i].XYZ());
            // ABT #685 (Alf, 2026-08-18): the bisector must be WELL CONDITIONED. Two unit
            // directions meeting at nearly 180 degrees sum to nearly nothing, and normalising
            // that gives an essentially arbitrary plane -- the knife then keeps the wrong side
            // and eats the piece. Measured on 10_emi: a riser/wrap joint whose keep direction
            // came out (0, 0.03, -0.9995), pointing along the wrap rather than across it, and the
            // cut removed 0.415 mm3 of a 0.421 mm3 half-revolution. |sum| = 2 cos(angle/2), so
            // requiring 0.2 stops mitring anything past ~168 degrees; such a joint is a fold, not
            // a corner, and is left tangent -- ungrown and uncut, which cannot overshoot.
            if (s.Magnitude() > kMinBisector) { nS = gp_Dir(s); bentS = true; }
        }
        bool bentE = false;
        gp_Dir nE = fe[i];
        if (i + 1 < n && worthMitring(fe[i].Angle(fs[i + 1]), ptrs[i]->kind, ptrs[i + 1]->kind)) {
            gp_Vec s(fe[i].XYZ());
            s += gp_Vec(fs[i + 1].XYZ());
            if (s.Magnitude() > kMinBisector) { nE = gp_Dir(s); bentE = true; }   // see above
        }
        // Growth per end: corners grow by over+mismatch (then trimmed on the mitre plane);
        // a TANGENT junction with an endpoint gap is bridged by growing the EARLIER prim's end
        // flush forward (no cut) -- one side only, so the bridge is never doubled.
        const double angS = (i > 0) ? fe[i - 1].Angle(fs[i]) : 0.0;
        const double angE = (i + 1 < n) ? fe[i].Angle(fs[i + 1]) : 0.0;
        const double growS0 = mitreGrow(angS) + dpS;
        const double overS = bentS ? growS0 + curveGrow(*ptrs[i], growS0, true) : 0.0;
        // A bridged end grows by the WEDGE the bend opens (r*tan(theta)) as well as any endpoint
        // mismatch -- growing only by dpE left a wedge gap of up to r*tan(theta) on the outer side
        // of every near-tangent joint, open copper-to-copper. One side only, so it is never doubled.
        // A bridged end no longer grows by the wedge r*tan(theta): that growth left copper
        // OUTSIDE the wire's envelope (see jointSphere). It still closes any endpoint MISMATCH
        // (dpE), which is real centreline, and the joint sphere fills the wedge.
        // The wedge is filled by EITHER the joint sphere (physically exact) or, while that is
        // gated off, the old axial growth. Never neither: dropping both would leave open
        // copper-to-copper wedges on every near-tangent joint.
        const bool sphereFill = std::getenv("MVB_JOINT_SPHERE") != nullptr;
        const double growE0 = mitreGrow(angE) + dpE;
        const double overE = bentE ? growE0 + curveGrow(*ptrs[i], growE0, false)
                                   : dpE + ((!sphereFill && angE > 1e-12) ? bridgeGrow(angE) : 0.0);
        // SPLIT A CLOSED REVOLUTION. ABT #685 (Alf, 2026-08-18): a SPIRAL that closes a full
        // turn has its own other end one pitch away ALONG THE COLUMN AXIS -- and a coil wound
        // with touching turns puts that end's copper exactly one wire radius from this one's
        // surface (measured: PQ33 pitch 0.959133 against 2r 0.959000, 0.13 um to spare; the
        // single-switch forward is flush to 0.06 um). A mitre knife must span the whole section
        // to sever it, so on those coils NO box can cut the corner without biting the piece
        // itself: too wide and it gouges the turn below and shears off a fragment (the flyback's
        // "dent below the mitre" and its stray [2/2] solid), exactly wide enough and its face
        // lies ON the cylinder, which is the grazing degeneracy that makes the boolean return
        // debris.
        // Sweeping the revolution as two HALVES dissolves the conflict: each half's other end is
        // 180 degrees away, far outside any knife, and the pass one pitch along the axis now
        // belongs to a DIFFERENT solid -- which a knife may overlap freely, since it only ever
        // subtracts from its own. The halves overlap at a tangent seam and are fused back into
        // ONE solid after both corners are trimmed, so the part list is unchanged.
        auto revolutionHalves = [&](const Primitive& pr) {
            std::vector<Primitive> out;
            const Spiral& sp = pr.spiral;
            const bool closed = pr.kind == Primitive::SPIRAL && !sp.blend &&
                                std::abs(sp.r1 - sp.r0) < 1e-12 &&
                                std::abs(sp.az1 - sp.az0) > 0.97 * kTwoPi;
            if (!closed) {
                out.push_back(pr);
                return out;
            }
            const double azM = 0.5 * (sp.az0 + sp.az1);
            const double yM = 0.5 * (sp.y0 + sp.y1);
            // The halves ABUT EXACTLY at the mid-azimuth: same centreline point, same tangent,
            // meeting on one flush perpendicular disc -- the assembler's own junction doctrine.
            // The previous seam OVERLAPPED the halves along the same analytic curve, which made
            // their tube surfaces exactly coincident over the window: the one input class OCC's
            // fuse cannot handle. On the pushpull's steep wrap it returned "done", one valid
            // solid -- half A alone, half B silently gone (A=147.242, B=147.294, fused=147.240
            // mm3). Flush abutment plus glue mode below is the coincidence OCC does support.
            Primitive a = pr, b = pr;
            a.spiral.az1 = azM;
            a.spiral.y1 = yM;
            b.spiral.az0 = azM;
            b.spiral.y0 = yM;
            out.push_back(a);
            out.push_back(b);
            return out;
        };
        // Set by whichever build path starts this piece on its predecessor's EXACT section
        // (prism cap adoption, or pipe/arc phase continuation). An exactly shared cap is an
        // internal face -- see the weld gate below.
        bool startCapAdopted = false;
        bool staggerHereForAudit = false;
        std::vector<Primitive> pieces = revolutionHalves(*ptrs[i]);
        // The section-wire handoff (see rawGrownSolid): the previous piece's ACTUAL end
        // section enters this piece as its start profile when the junction is tangent and
        // exact, and this piece's own end section is captured for the next.
        //
        // TRIED AND REVERTED (2026-08-30): restricting the continuation to ANALYTICALLY tangent
        // junctions (angS <= 1e-12) on the theory that a bend makes the two sections coplanar
        // and their side facets graze. Measured: 12_boost fell from 262/262 to 1/262 -- its
        // toroid chain NEEDS the section continued across its mitred corners -- and 04_forward
        // did not recover. The junction gate is not what 04 is about.
        // THE GATE IS ENDPOINT SHARING, NOT TANGENCY (2026-08-30, derived then measured).
        // At a corner with directions dA, dB the mitre plane normal is n proportional to
        // dA + dB, and reflection through that plane maps dA -> -dB while FIXING every point
        // lying in the plane. So a neighbour's cap polygon is a valid section of THIS piece's
        // cylinder -- adoption is exact and undistorted -- precisely when the plane is the true
        // angular bisector through the shared corner point. That holds however sharp the bend
        // is, which is why gating on tangency destroyed 12_boost (its toroid corners are
        // mitred but endpoint-exact: 262/262 -> 1/262).
        // It FAILS when the two pieces do not share an endpoint: this code then puts the plane
        // through the MIDPOINT of the two endpoints (see Ps below), the corner is no longer
        // symmetric, the reflection no longer maps one cylinder onto the other, and the adopted
        // cap is not a section of the receiving piece -- it sweeps into a solid that bulges
        // outside the nominal wire (measured 24 um on 04_forward's lead corner) and grazes its
        // neighbour's side facet, which is what cost 04 413 of its 418 volumes.
        // Endpoint mismatches are a real datum from the layout (up to a wire radius, 0.36 mm
        // measured on PQ33), so this gate is physical, not a tuned threshold.
        // A SHARED SECTION NEEDS A SHARED CUTTING PLANE. There are exactly three kinds of
        // junction, and only two of them have one:
        //   MITRED (bentS): both pieces are cut by the SAME bisector plane, so a shared cap is
        //     exact and both pieces' facets terminate on that plane -- 12_boost's toroid
        //     corners, which need this (without it: 262/262 -> 1/262).
        //   COLLINEAR (angS ~ 0): one axis, so the shared cap is a plain continuation.
        //   BRIDGED WITH A BEND (not mitred, angS > 0): NO common plane -- each piece keeps its
        //     own perpendicular cap, tilted from the other by the bend angle. Sharing the
        //     section here makes the two pieces' side facets NEARLY parallel instead of skew,
        //     so they graze over a strip rather than crossing at a line, and gmsh meshes the
        //     two faces independently into overlapping triangles. That is 04_forward's lead
        //     junction: 413 of 418 volumes lost to surfaces 220/237, unchanged by every gate
        //     that varied HOW MUCH is adopted (the caps themselves are geometrically exact --
        //     parallel planes, projected radius exactly r -- so no cap-quality test can see it).
        // Endpoint sharing is still required: without it the plane is a fabricated midpoint,
        // not a bisector (see the reflection argument below).
        // TRIED AND REVERTED (2026-08-30): gating the handoff offer on endpoint sharing
        // (dpS <= 1e-9) and on shared-cutting-plane junction kinds (bentS || collinear). Both
        // kept 12_boost green but left 04 red -- and neither was ever validated on 10_emi,
        // whose green state was measured with the offer UNCONDITIONAL. The receiving piece's
        // own geometric fit check remains the only gate.
        const bool junctionEndpointExact = i > 0;
        // Did THIS piece receive its predecessor's exact end cap? Read before the piece is built
        // (building it overwrites prevEndCapValid for the next one). Used by the abutment gate.
        const bool receivedCap = junctionEndpointExact && prevEndCapValid;

        // Pipes/revolves receiving from a LONG STRAIGHT source get the staggered offer (the
        // cap rotated half a facet about the source axis): 04_forward's corner pipes receive
        // from its 12 mm runs and were the remaining osculation sites (prism-side-only stagger
        // left 04 at 5/418; all-pipe stagger fixed 04 but broke 10's dragbacks, whose source
        // is the 0.9 mm riser -- under the 2*pi*r length boundary those stay verbatim).
        TopoDS_Wire startCapWire;
        if (junctionEndpointExact && prevEndCapValid) {
            const bool staggerPipeOffer =
                i > 0 && !bentS && dpS <= 1e-9 && ptrs[i]->kind != Primitive::SEG &&
                ptrs[i - 1]->kind == Primitive::SEG &&
                ptrs[i - 1]->seg.a.Distance(ptrs[i - 1]->seg.b) > kTwoPi * wireRadius;
            if (staggerPipeOffer) {
                gp_Vec pd(ptrs[i - 1]->seg.a, ptrs[i - 1]->seg.b);
                if (pd.Magnitude() > 1e-12) {
                    gp_Trsf rot;
                    rot.SetRotation(gp_Ax1(primEndpoints(*ptrs[i - 1]).second, gp_Dir(pd)),
                                    kPi / segments);
                    std::vector<gp_Pnt> rotated = prevEndCap;
                    for (auto& p : rotated) p = p.Transformed(rot);
                    startCapWire = pointsToWire(rotated);
                    if (std::getenv("MVB_CAP_DEBUG"))
                        std::cerr << "[cap-debug] staggered pipe offer at '" << ptrs[i]->label
                                  << "'\n";
                }
            }
            if (startCapWire.IsNull()) startCapWire = pointsToWire(prevEndCap);
        }
        TopoDS_Wire pieceEndCapWire;
        auto sweepPieces = [&](double a0, double b0) {
            std::vector<TopoDS_Shape> got;
            for (size_t q = 0; q < pieces.size(); ++q) {
                TopoDS_Wire endW;
                TopoDS_Shape sp = rawGrownSolid(
                    pieces[q], wireRadius, q == 0 ? a0 : 0.0,
                    q + 1 == pieces.size() ? b0 : 0.0, segments,
                    (q == 0 && a0 <= 0.0 && !startCapWire.IsNull()) ? &startCapWire : nullptr,
                    q + 1 == pieces.size() ? &endW : nullptr,
                    q == 0 ? &startCapAdopted : nullptr);
                if (sp.IsNull()) return std::vector<TopoDS_Shape>{};
                got.push_back(sp);
                if (q + 1 == pieces.size()) pieceEndCapWire = endW;
            }
            return got;
        };
        // FACETED SEG: build the mitred prism EXACTLY and skip both the growth and the knife
        // (see mitredFacetPrism -- coplanar knife faces on polygon facets are what left every
        // toroid tube self-intersecting at seg=16). End planes: the junction bisector where the
        // end is mitred, the perpendicular cap (grown by the endpoint mismatch) where tangent.
        bool prismDone = false;
        TopoDS_Shape solid;
        // ROUND SEG WITH A MITRED END: the exact elliptical-section prism, never the knife
        // (see mitredRoundPrism). Flush-ended round segments keep today's plain cylinder.
        if (ptrs[i]->kind == Primitive::SEG && segments <= 0 && (bentS || bentE)) {
            const gp_Pnt& A = ptrs[i]->seg.a;
            const gp_Pnt& B = ptrs[i]->seg.b;
            gp_Vec dv(A, B);
            if (dv.Magnitude() > 1e-12) {
                const gp_Dir dir(dv);
                gp_Pnt Ps = A; gp_XYZ nsx = dir.XYZ();
                if (bentS) {
                    Ps = gp_Pnt(0.5 * (primEndpoints(*ptrs[i - 1]).second.XYZ() + A.XYZ()));
                    nsx = nS.XYZ();
                } else if (dpS > 0) {
                    Ps = A.Translated(gp_Vec(dir) * (-dpS));
                }
                gp_Pnt Pe = B; gp_XYZ nex = dir.XYZ();
                if (bentE) {
                    Pe = gp_Pnt(0.5 * (B.XYZ() + primEndpoints(*ptrs[i + 1]).first.XYZ()));
                    nex = nE.XYZ();
                } else if (dpE > 0) {
                    Pe = B.Translated(gp_Vec(dir) * dpE);
                }
                std::string why;
                TopoDS_Shape prism = mitredRoundPrism(A, dir, wireRadius, Ps, nsx, Pe, nex, &why);
                if (!prism.IsNull()) {
                    solid = prism;
                    prismDone = true;
                    if (bentS) ++nCut;
                    if (bentE) ++nCut;
                    // A round cap is an ellipse, not a polygon: nothing to hand on.
                    prevEndCapValid = false;
                    prevEndCapTrusted = false;
                    prevEndCapSourceStraight = true;
                    prevSourceLen = dv.Magnitude();
                } else if (std::getenv("MVB_MITRE_DIAG")) {
                    std::cerr << "[round-prism] '" << ptrs[i]->label
                              << "' exact mitred cylinder refused (" << why
                              << "); falling back to the knife\n";
                }
            }
        }
        if (!prismDone && ptrs[i]->kind == Primitive::SEG && segments > 0) {
            const gp_Pnt& A = ptrs[i]->seg.a;
            const gp_Pnt& B = ptrs[i]->seg.b;
            gp_Vec dv(A, B);
            if (dv.Magnitude() > 1e-12) {
                const gp_Dir dir(dv);
                gp_Pnt Ps = A; gp_XYZ nsx = dir.XYZ();
                if (bentS) {
                    Ps = gp_Pnt(0.5 * (primEndpoints(*ptrs[i - 1]).second.XYZ() + A.XYZ()));
                    nsx = nS.XYZ();
                } else if (dpS > 0) {
                    Ps = A.Translated(gp_Vec(dir) * (-dpS));
                }
                gp_Pnt Pe = B; gp_XYZ nex = dir.XYZ();
                if (bentE) {
                    Pe = gp_Pnt(0.5 * (B.XYZ() + primEndpoints(*ptrs[i + 1]).first.XYZ()));
                    nex = nE.XYZ();
                } else if (dpE > 0) {
                    Pe = B.Translated(gp_Vec(dir) * dpE);
                }
                // Offer the previous prism's exact end cap across EVERY junction, mitred or
                // tangent -- a tangent junction's two caps are the same polygon projected onto
                // essentially the same plane, and their fp-level mismatch is exactly what makes
                // the glue-fragment prefer one welded body there (2026-08-27 measurement). The
                // 100 nm adoption gate in mitredFacetPrism rejects any junction where the caps
                // are not the same face (endpoint mismatches, phase changes).
                std::vector<gp_Pnt> endCap;
                // Trusted mode when the previous piece was a pipe/revolve at a tangent
                // junction: its cap's phase is legitimately different from the deterministic
                // frame (see mitredFacetPrism).
                // THE PRISM TAKES THE NEIGHBOUR'S SECTION AT EVERY JUNCTION, mitred corners
                // included -- gated only by the geometric plane/radius check inside
                // mitredFacetPrism. TRIED AND REVERTED (2026-08-30): gating this on
                // junctionTangentIn "for symmetry with the pipe path" took 12_boost from
                // 262/262 to 1/262. A mitred corner's two prisms flank one bisector plane, and
                // continuing the section is exactly what makes their shared faces coincide.
                // ...but only where the two pieces actually SHARE the endpoint, so the mitre
                // plane is the true bisector and the neighbour's cap is a real section of this
                // cylinder (see junctionEndpointExact above).
                // STAGGER RULE (2026-08-30, measured on the 04/10/12/02 quartet): a tangent
                // junction whose SOURCE piece is straight (SEG) gets the staggered
                // continuation -- verbatim adoption there leaves the two pieces' coplanar
                // facets overlapping over the axial mismatch band, the osculation class that
                // kept 04_forward at 5/418 through every other rule (full stagger fixed it,
                // 418/418). Everything else -- mitred corners (12_boost requires the exact
                // shared face) and curved sources (10_emi's ride stack, where a staggered
                // chord cuts the neighbour's tightly curved strip) -- stays verbatim.
                const bool tangentIn = i > 0 && !bentS && angS <= 1e-12;
                // Qualify by the RECEIVER's length: the osculation band lives on the receiving
                // prism's own facet (04's 12 mm runs, 48 r, adopting from corner pipes -- red
                // verbatim at 5/418, green staggered at 418/418), while a short receiver has
                // nothing to osculate along and NEEDS the exact face (10_emi's 0.9 mm risers,
                // 3.6 r: green verbatim, 5/50 staggered). Boundary: one section circumference.
                const bool staggerHere =
                    tangentIn && dv.Magnitude() > kTwoPi * wireRadius;
                staggerHereForAudit = staggerHere;   // a staggered adoption is a half-facet rotation: NOT the same face
                std::string facetWhy;
                TopoDS_Shape prism = mitredFacetPrism(
                    A, dir, wireRadius, segments, Ps, nsx, Pe, nex,
                    (junctionEndpointExact && prevEndCapValid) ? &prevEndCap : nullptr, &endCap,
                    /*sharedTrusted=*/prevEndCapTrusted, &startCapAdopted,
                    /*staggerAdopt=*/staggerHere, &facetWhy);
                if (prism.IsNull() && std::getenv("MVB_MITRE_DIAG")) {
                    std::cerr << "[facet-prism] '" << ptrs[i]->label << "' len " << dv.Magnitude() * 1e6
                              << " um r " << wireRadius * 1e6 << " um bentS=" << bentS << " bentE=" << bentE
                              << " refused (" << facetWhy << "); falling back to the swept piece + knife\n";
                }
                if (!prism.IsNull()) {
                    solid = prism;
                    prismDone = true;
                    if (bentS) ++nCut;
                    if (bentE) ++nCut;
                    prevEndCap = std::move(endCap);
                    prevEndCapValid = true;
                    prevEndCapTrusted = false;   // a prism cap follows the deterministic frame
                    prevEndCapSourceStraight = true;   // a SEG prism is straight by definition
                    prevSourceLen = ptrs[i]->seg.a.Distance(ptrs[i]->seg.b);
                }
            }
        }
        if (!prismDone) prevEndCapValid = false;
        std::vector<TopoDS_Shape> parts;
        if (!prismDone) parts = sweepPieces(overS, overE);
        if (!prismDone && parts.empty()) {  // ARC clamp (near-full revolve): flush, uncut tube
            parts = sweepPieces(0.0, 0.0);
            bentS = bentE = false;
        }
        if (!prismDone) solid = parts.empty() ? TopoDS_Shape() : parts.front();
        if (!prismDone && !pieceEndCapWire.IsNull()) {
            // A pipe/revolve delivered its true end section: hand it on (trusted) so a
            // following riser prism starts from the section the copper actually ends with.
            auto capPts = wireToPoints(pieceEndCapWire);
            if ((int)capPts.size() == segments) {
                prevEndCap = std::move(capPts);
                prevEndCapValid = true;
                prevEndCapTrusted = true;
                prevEndCapSourceStraight = false;   // pipes/revolves are curved at their ends
                prevSourceLen = 0.0;
            } else if (std::getenv("MVB_CAP_DEBUG")) {
                std::cerr << "[cap-debug] end cap wire has " << capPts.size()
                          << " vertices (want " << segments << ") for '" << ptrs[i]->label
                          << "'\n";
            }
        } else if (!prismDone && std::getenv("MVB_CAP_DEBUG")) {
            std::cerr << "[cap-debug] no end cap delivered by '" << ptrs[i]->label << "'\n";
        }
        if (solid.IsNull()) {
            // Never skip silently: a missing primitive is missing copper, and its neighbours
            // may still touch each other so even the connectivity contract can stay satisfied
            // (measured: faceted blend sweeps returned null and every rect-column turn quietly
            // detached from the next).
            static const char* kindName[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
            const auto ends = primEndpoints(*ptrs[i]);
            throw std::runtime_error(
                "ConductorBuilder: the conformal assembler could not sweep primitive " +
                std::to_string(i) + " (" + kindName[ptrs[i]->kind] + ") '" + ptrs[i]->label +
                "' (" + std::to_string(ends.first.X()) + "," + std::to_string(ends.first.Y()) +
                "," + std::to_string(ends.first.Z()) + ")->(" + std::to_string(ends.second.X()) +
                "," + std::to_string(ends.second.Y()) + "," + std::to_string(ends.second.Z()) +
                "). Refusing to emit a conductor with missing copper.");
        }
        // How close the SAME wire comes to this joint on a neighbouring pass: for a spiral that
        // is its axial advance per revolution. The knife may not reach that far, or it cuts the
        // wire's other pass instead of the corner. Zero = unbounded (no neighbouring pass).
        auto junctionDiag = [](const gp_Dir& a, const gp_Dir& b, const Primitive& pr) {
            static const char* kn[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
            std::ostringstream d;
            d << "joint angle " << a.Angle(b) * 180.0 / kPi << " deg (in ("
              << a.X() << "," << a.Y() << "," << a.Z() << ") out (" << b.X() << "," << b.Y()
              << "," << b.Z() << ")), cut piece is " << kn[pr.kind];
            if (pr.kind == Primitive::SPIRAL) {
                d << " az span " << (pr.spiral.az1 - pr.spiral.az0) * 180.0 / kPi << " deg r "
                  << pr.spiral.r0 * 1e3 << " y " << pr.spiral.y0 * 1e3 << "->"
                  << pr.spiral.y1 * 1e3 << " mm";
            }
            return d.str();
        };
        auto reachLimitOf = [&](const Primitive& pr) {
            if (pr.kind != Primitive::SPIRAL) return 0.0;
            const double sweep = pr.spiral.az1 - pr.spiral.az0;
            if (std::abs(sweep) < 1e-9) return 0.0;
            // ONLY A CLOSED REVOLUTION HAS COPPER OF ITS OWN BESIDE THE CORNER. ABT #685 (Alf,
            // 2026-08-18): the knife subtracts from THIS solid and nothing else, so a neighbouring
            // pass that belongs to a different primitive is out of its reach by construction --
            // bounding against it only shrank the knife below the tube it had to sever. Applying
            // the bound to every spiral refused ordinary windings: PQ33's half-revolution wrap
            // 8->9 was rejected for a pass 0.959133 mm away that is not part of it at all, and
            // that design had been building for weeks. A full turn IS its own neighbour -- and is
            // split into halves above precisely so that neither half is.
            if (std::abs(sweep) < 0.97 * kTwoPi) return 0.0;
            const double pitch = std::abs((pr.spiral.y1 - pr.spiral.y0) * kTwoPi / sweep);
            if (pitch < 1e-12) return 0.0;
            return pitch;   // the neighbouring pass sits exactly one pitch along the axis
        };
        const double reachHere = reachLimitOf(*ptrs[i]);
        // Once split, a half's own other end is half a revolution away, so nothing of this piece
        // sits near either knife: the neighbour bound only applies to a still-closed revolution.
        const double neighbourHere = pieces.size() > 1 ? 0.0 : reachHere;
        if (!prismDone && bentS) {
            const gp_Pnt J(0.5 * (primEndpoints(*ptrs[i - 1]).second.XYZ() +
                                  primEndpoints(*ptrs[i]).first.XYZ()));
            parts.front() =
                localMitreTrim(parts.front(), J, nS, gp_Dir(fs[i].XYZ() * -1.0),
                               primEndpoints(*ptrs[i]).first, wireRadius, overS,
                               "'" + ptrs[i - 1]->label + "' -> '" + ptrs[i]->label + "'",
                               neighbourHere, junctionDiag(fe[i - 1], fs[i], *ptrs[i]));
            ++nCut;
        }
        if (!prismDone && bentE) {
            const gp_Pnt J(0.5 * (primEndpoints(*ptrs[i]).second.XYZ() +
                                  primEndpoints(*ptrs[i + 1]).first.XYZ()));
            parts.back() =
                localMitreTrim(parts.back(), J, gp_Dir(nE.XYZ() * -1.0), fe[i],
                               primEndpoints(*ptrs[i]).second, wireRadius, overE,
                               "'" + ptrs[i]->label + "' -> '" + ptrs[i + 1]->label + "'",
                               neighbourHere, junctionDiag(fe[i], fs[i + 1], *ptrs[i]));
            ++nCut;
        }
        if (!prismDone) solid = parts.empty() ? TopoDS_Shape() : parts.front();
        if (parts.size() > 1) {
            // Fuse the halves back: both corners are cut, so the closed revolution is safe again.
            gp_Trsf up, down;
            up.SetScale(gp_Pnt(0, 0, 0), 1000.0);
            down.SetScale(gp_Pnt(0, 0, 0), 1.0 / 1000.0);
            try {
                // Glue mode: the halves meet on one exactly-shared flush disc, and
                // BOPAlgo_GlueFull is OCC's declared contract for coincident faces -- the
                // general boolean treats them as an intersection problem and loses a half.
                BRepAlgoAPI_Fuse fu;
                TopTools_ListOfShape args, tools;
                args.Append(BRepBuilderAPI_Transform(parts[0], up, Standard_True).Shape());
                tools.Append(BRepBuilderAPI_Transform(parts[1], up, Standard_True).Shape());
                fu.SetArguments(args);
                fu.SetTools(tools);
                fu.SetGlue(BOPAlgo_GlueShift);
                fu.Build();
                if (!fu.IsDone() || fu.Shape().IsNull())
                    throw std::runtime_error("fuse of the two half-revolutions did not complete");
                solid = BRepBuilderAPI_Transform(fu.Shape(), down, Standard_True).Shape();
                // THE FUSE MUST ACCOUNT FOR EVERY CUBIC MICRON. ABT #685 (Alf, 2026-08-19,
                // "why is turn 2 -> turn 2_ending not connecting with the layer link?"): on the
                // pushpull's steep single-turn wrap, BRepAlgoAPI_Fuse reported IsDone, returned
                // one VALID solid -- and that solid was half A alone (A=147.242, B=147.294,
                // fused=147.240 mm3). Half B, the half that reaches the layer link, vanished
                // without a diagnostic; the conductor shipped visibly disconnected. IsDone and
                // IsValid cannot see this, only conservation can: the halves coincide exactly
                // over the seam (same analytic curve), so the union's volume must be
                // A + B - seamOverlap, where seamOverlap is the tube over the shared arc --
                // computable, not estimated. If the fuse comes up short, keep BOTH halves as
                // flush-overlapping solids (the weld-lens pattern the rect path already uses;
                // the consumer's fragment welds them) and say so.
                {
                    GProp_GProps gA, gB, gF;
                    BRepGProp::VolumeProperties(parts[0], gA);
                    BRepGProp::VolumeProperties(parts[1], gB);
                    BRepGProp::VolumeProperties(solid, gF);
                    // Flush abutment: no shared volume, the union must carry every cubic
                    // micron of both halves. The failure mode loses an entire half (~50%); one
                    // percent is pure measurement headroom, not a configuration in between.
                    const double expected = gA.Mass() + gB.Mass();
                    if (gF.Mass() < expected * 0.99) {
                        std::cerr << "[ConductorBuilder] fuse of '" << ptrs[i]->label
                                  << "' lost copper (A=" << gA.Mass() * 1e9
                                  << " B=" << gB.Mass() * 1e9 << " fused=" << gF.Mass() * 1e9
                                  << " mm3, expected >= " << expected * 0.99 * 1e9
                                  << "); keeping both halves as flush-abutting solids for the "
                                     "consumer's weld."
                                  << std::endl;
                        TopoDS_Compound both;
                        BRep_Builder bb2;
                        bb2.MakeCompound(both);
                        bb2.Add(both, parts[0]);
                        bb2.Add(both, parts[1]);
                        solid = both;
                    }
                    if (std::getenv("MVB_FUSE_DIAG")) {
                        int ns = 0;
                        for (TopExp_Explorer e(solid, TopAbs_SOLID); e.More(); e.Next()) ++ns;
                        std::cerr << "[fuse] '" << ptrs[i]->label << "' A=" << gA.Mass() * 1e9
                                  << " B=" << gB.Mass() * 1e9 << " fused=" << gF.Mass() * 1e9
                                  << " mm3, solids=" << ns << std::endl;
                        int si = 0;
                        for (TopExp_Explorer e(solid, TopAbs_SOLID); e.More(); e.Next(), ++si) {
                            Bnd_Box bx;
                            BRepBndLib::Add(e.Current(), bx);
                            double x0, y0z, z0, x1, y1z, z1;
                            bx.Get(x0, y0z, z0, x1, y1z, z1);
                            std::cerr << "[fuse]   solid " << si << " bbox x[" << x0 * 1e3 << ","
                                      << x1 * 1e3 << "] y[" << y0z * 1e3 << "," << y1z * 1e3
                                      << "] z[" << z0 * 1e3 << "," << z1 * 1e3 << "] mm"
                                      << std::endl;
                        }
                    }
                }
            } catch (const Standard_Failure& f) {
                throw std::runtime_error(
                    "ConductorBuilder: could not rejoin the two halves of the closed revolution '" +
                    ptrs[i]->label + "' (" +
                    std::string(f.GetMessageString() ? f.GetMessageString() : "(null)") + ")");
            }
        }
        // STEP round-trip robustness: the boolean's cut curve on a torus can carry an edge tolerance
        // that a STEP export degrades into an invalid solid (valid in memory, invalid on reload).
        // ShapeFix tightens edges/tolerances in place so the written solid survives the round-trip.
        if (!solid.IsNull() && (bentS || bentE)) {
            try {
                ShapeFix_Shape fix(solid);
                fix.Perform();
                if (!fix.Shape().IsNull()) {
                    if (std::getenv("MVB_FUSE_DIAG")) {
                        GProp_GProps g0, g1;
                        BRepGProp::VolumeProperties(solid, g0);
                        BRepGProp::VolumeProperties(fix.Shape(), g1);
                        if (std::abs(g1.Mass() - g0.Mass()) > 0.01 * std::abs(g0.Mass()))
                            std::cerr << "[shapefix] '" << ptrs[i]->label << "' volume "
                                      << g0.Mass() * 1e9 << " -> " << g1.Mass() * 1e9 << " mm3"
                                      << std::endl;
                    }
                    solid = fix.Shape();
                }
            } catch (const Standard_Failure&) {
            }
        }
        // A slice that self-intersected or emptied the tube (very short/sharp primitive) is caught
        // here and rebuilt as a plain flush-capped tube -- always valid, at worst a hair of overlap
        // at that single joint. Never emit an invalid solid.
        if (solid.IsNull() || !BRepCheck_Analyzer(solid).IsValid()) {
            TopoDS_Shape flush = rawGrownSolid(*ptrs[i], wireRadius, 0.0, 0.0, segments);
            if (!flush.IsNull() && BRepCheck_Analyzer(flush).IsValid()) { solid = flush; ++nRepaired; }
            else {
                // NEVER drop a primitive silently: the copper it carries simply disappears from
                // the conductor, and downstream checks can miss it (a short corner's neighbours
                // may still touch, so even the connectivity contract stays satisfied -- that is
                // exactly how a missing dragback corner reached a STEP). Fail loudly instead,
                // naming the primitive so the geometry that produced it can be fixed.
                static const char* kindName[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
                const auto ends = primEndpoints(*ptrs[i]);
                std::ostringstream detail;
                detail << "ConductorBuilder: the conformal assembler cannot build a valid solid for "
                       << "primitive " << i << " (" << kindName[ptrs[i]->kind] << ") '"
                       << ptrs[i]->label << "' (" << ends.first.X() << "," << ends.first.Y() << ","
                       << ends.first.Z() << ")->(" << ends.second.X() << "," << ends.second.Y()
                       << "," << ends.second.Z() << ")";
                if (ptrs[i]->kind == Primitive::ARC3) {
                    detail << " sweep=" << ptrs[i]->arc.sweep
                           << " centrelineRadius=" << ptrs[i]->arc.v0.Modulus()
                           << " wireRadius=" << wireRadius;
                    if (ptrs[i]->arc.v0.Modulus() <= wireRadius * (1.0 + 1e-6)) {
                        detail << " -- the centreline radius does not exceed the wire radius, so the "
                                  "swept tube is a HORN TORUS touching its own axis of revolution "
                                  "(no valid solid exists; give the corner a larger bend radius)";
                    }
                }
                detail << ". Refusing to emit a conductor with missing copper.";
                throw std::runtime_error(detail.str());
            }
        }
        // ABT #685: add SOLIDS, never a sub-compound. A mitre trim can fragment one primitive
        // into two solids, and adding that compound as one shape made the assembly's component
        // count differ from its solid count -- 21 solids, 20 components on the pushpull's
        // Primary 1. Every consumer that walks solids and then indexes components (the STEP
        // exporter's per-solid naming, above all) is off by one from that point on, which is how
        // 'turn 15 -> turn 16' ended up naming the wrong piece of copper. One solid, one
        // component, one name.
        // A TANGENT joint (no bend, no bridging growth, no gap) must abut exactly as well: the
        // two pieces share one end cap, so any common volume is a construction error in one of
        // the caps, not a bridged stub. Bridged joints (angle > 1e-12, grown on purpose) stay
        // exempt.
        // A FACETED PRISM THAT ADOPTED ITS PREDECESSOR'S CAP VERBATIM SHARES THAT FACE BY
        // CONSTRUCTION (same vertices, same plane): there is no lens to measure. Running the
        // Common anyway is what hung single_switch at --segments 12 (28 min inside
        // BOPAlgo_PaveFiller::PerformEF <- BRepAlgoAPI_Common, from this very call: twelve
        // coincident facets per joint, two frames per joint, every joint) and what crashed
        // 14_dab (the same chain, one frame deeper, in Extrema_GenExtPS::BuildGrid). The gate
        // stays on every joint whose caps are NOT provably the same face -- refused prisms,
        // knife-trimmed pieces, curved sources -- which is where a lens can exist.
        // ...and the same holds for a faceted PIPE/REVOLVE at a TANGENT joint that received the
        // cap: the sweep starts on the neighbour's polygon (startProfileOverride), so the two
        // caps are one face by construction. single_switch's last line before the 30-minute
        // hang was exactly such a joint: "junction 49->50 SPIRAL->SPIRAL angle=1.4e-14 deg".
        // Mitred (bent) joints keep the gate whatever the piece kind: their caps are cut, not
        // handed over, and a lens there is a real defect.
        // ...and at segments <= 0 the two exact round caps ARE one disc, verified on the built
        // solids rather than inferred from the handoff that does not exist there (see capDiscAt,
        // ABT #1186). Only a tangent, endpoint-exact, ungrown junction qualifies: a mitred joint's
        // caps are cut, not shared, and a lens there is a real defect.
        bool roundCapsShared = false;
        if (segments <= 0 && i > 0 && !bentS && angS <= 1e-12 && dpS <= 1e-9 &&
            !sphereAtPrevJunction && !prevBuilt.IsNull() && !solid.IsNull()) {
            const gp_Pnt jPt = primEndpoints(*ptrs[i]).first;
            double aPrev = 0.0, aCur = 0.0, offPrev = 0.0, offCur = 0.0;
            roundCapsShared = capDiscAt(prevBuilt, jPt, fs[i], wireRadius, &aPrev, &offPrev) &&
                              capDiscAt(solid, jPt, fs[i], wireRadius, &aCur, &offCur);
            if (diag)
                std::cerr << "[abut]   round-cap exemption "
                          << (roundCapsShared ? "TAKEN   " : "declined") << " '"
                          << ptrs[i - 1]->label << "' -> '" << ptrs[i]->label << "' capArea prev="
                          << aPrev * 1e6 << " cur=" << aCur * 1e6 << " mm2 (want "
                          << kPi * wireRadius * wireRadius * 1e6 << "), planeOff prev="
                          << offPrev * 1e9 << " cur=" << offCur * 1e9 << " nm\n";
        }
        const bool capsProvablyShared = (prismDone && startCapAdopted && !staggerHereForAudit)
                                     || (segments > 0 && !bentS && receivedCap)
                                     || roundCapsShared;
        if (i > 0 && (bentS || (angS <= 1e-12 && dpS <= 1e-9))) {
            if (capsProvablyShared) ++abutChecked;   // verdict: shared face, zero lens, by construction
            else checkMitredAbutment(prevBuilt, solid,
                                std::string(bentS ? "" : "tangent ") + "'" + ptrs[i - 1]->label +
                                    "' -> '" + ptrs[i]->label + "'");
        }
        prevBuilt = solid;
        // Collect per PRIMITIVE instead of emitting straight into the compound: bridged
        // junctions have to be fused after the fact (see the run-fuse below), and that needs
        // the pieces grouped, not already flattened.
        // CHUNK 7. This primitive's END junction is BRIDGED (not mitred), so the notch on the
        // outside of the bend is filled by a SPHERE at the corner rather than by growing this tube
        // past it (see jointSphere). Fused into THIS piece immediately: emitting it as its own
        // solid multiplies the assembly's solid count by the number of joints, and the bobbin cut
        // then drowns -- OCC segfaulted in BOPAlgo_PaveFiller::PerformEF via cut_bobbin on 14_dab.
        // The sphere overlaps the tube end substantially, so this is a well-conditioned union, not
        // the near-tangent kind that fails.
        // GATED OFF BY DEFAULT, and not because the geometry is wrong -- it is the physically
        // exact fill. OCC cannot carry it at this scale on this geometry: measured on 14_dab at
        // segments=16 the build went 37 s -> 815 s (a boolean per joint on a 112-face sphere) and
        // then SEGFAULTED downstream of the weld, with 3 spheres refusing to fuse at all. Emitting
        // the spheres as separate solids instead segfaulted OCC in cut_bobbin
        // (BOPAlgo_PaveFiller::PerformEF). Turn on with MVB_JOINT_SPHERE=1; until the cost and the
        // crash are resolved the bridged wedge is left as it was.
        if (std::getenv("MVB_JOINT_SPHERE") && i + 1 < n && angE > 1e-12 && !bentE) {
            const gp_Pnt jA = primEndpoints(*ptrs[i]).second;
            const gp_Pnt jB = primEndpoints(*ptrs[i + 1]).first;
            std::vector<gp_Pnt> corners{jA};
            if (jA.Distance(jB) > 1e-9) corners.push_back(jB);
            for (const auto& cpt : corners) {
                TopoDS_Shape sph = jointSphere(cpt, wireRadius, segments);
                if (sph.IsNull() || !BRepCheck_Analyzer(sph).IsValid()) {
                    std::cerr << "[joint-sphere] '" << ptrs[i]->label
                              << "': corner fill not built; the wedge stays open" << std::endl;
                    continue;
                }
                bool ok = false;
                try {
                    BRepAlgoAPI_Fuse fu(solid, sph);
                    fu.Build();
                    if (fu.IsDone() && !fu.Shape().IsNull() &&
                        BRepCheck_Analyzer(fu.Shape()).IsValid()) {
                        GProp_GProps g0, g1;
                        BRepGProp::VolumeProperties(solid, g0);
                        BRepGProp::VolumeProperties(fu.Shape(), g1);
                        // A union only adds; it must not lose the tube it started from.
                        if (g1.Mass() >= g0.Mass() * (1.0 - 1e-9)) { solid = fu.Shape(); ok = true; }
                    }
                } catch (const Standard_Failure&) {
                }
                if (!ok)
                    std::cerr << "[joint-sphere] '" << ptrs[i]->label
                              << "': corner fill would not fuse; the wedge stays open" << std::endl;
                else
                    sphereAtPrevJunction = true;   // piece i+1 must weld onto this one
            }
        }
        perPrimSolids.emplace_back();
        for (TopExp_Explorer solidExplorer(solid, TopAbs_SOLID); solidExplorer.More();
             solidExplorer.Next()) {
            perPrimSolids.back().push_back(solidExplorer.Current());
        }
        // Junction (i-1 -> i) needs welding only if it was BRIDGED **and** the bridge actually
        // grew: the wedge fill is r*tan(theta), so a junction at theta ~ 0 adds no copper and has
        // nothing to weld away. Nearly every wrap junction is tangent to ~1e-14 deg, so welding
        // "every bridged junction" ran a fuse per turn against an ever-growing accumulator for no
        // benefit -- measured on 01_etd34 at segments=16: 665 s with the weld, 27 s without, with
        // NURBS accounting for none of it. Gating on the same angle that creates the overlap is
        // what makes the weld proportional to the defect instead of to the turn count.
        // MVB_WELD_ALL=1 restores welding at EVERY bridged junction, tangent ones included.
        // The angS gate cut 01_etd34's seg=16 build from 665 s to 29 s, but the resulting
        // artifact regressed from meshable (32/32 volumes) to boundary-recovery failure: a
        // tangent junction's two pieces abut on bit-identical 16-gon faces, and gmsh's
        // fragment+recovery handles ONE welded body better than a glued pair there. The A/B
        // knob exists so speed vs meshability is a measurement, not an argument.
        const bool weldAll = std::getenv("MVB_WELD_ALL") != nullptr;
        // Mitre corners (bentS) weld ONLY when this piece was built on its predecessor's EXACT
        // cap polygon (startCapAdopted). Welding every mitre was tried on 2026-08-28 and made
        // the geometry worse -- 4 non-manifold solids, 0/24 meshed -- but that was before the
        // section handoff, when the two sides' bisector-plane faces were merely co-planar and
        // not coincident, which is the one input class OCC's fuse cannot do. With an adopted
        // cap the fuse is two solids sharing one identical planar face. It is still only an
        // OFFER: the weld's own acceptance test (valid + volume within the union bounds)
        // rejects a bad fuse, and the pieces then stay abutting exactly as before, so this can
        // never ship a fused solid that is worse than the pair.
        // Why it matters (04_forward): two mitred prisms flanking one bisector plane also have
        // near-parallel SIDE facets, and with the sections aligned those facets touch over a
        // strip instead of meeting along an edge -- COMMON = 0.000000 mm^3 over an 80 x 278 x
        // 87 um patch between 'Primary parallel 0' solids 7 and 8, which the audit passes and
        // gmsh then triangulates twice ("Invalid boundary mesh (overlapping facets)", 413 of
        // 418 volumes lost). Fusing the pair makes that contact INTERNAL, which is what it
        // physically is: one wire, not two touching bodies.
        // TRIED AND REVERTED (2026-08-30): offering the weld at every junction whose cap was
        // adopted (mitred ones included) to make 04_forward's grazing contact internal. It
        // welded hard -- 418 solids down to 186 -- and did NOT fix it, because the graze was
        // never a junction: it was a piece BULGING outside the nominal wire from a bad adopted
        // cap (see the projected-radius test above). Fix the cap, not the topology.
        bridgedToPrev.push_back((i > 0 && !bentS && (weldAll || angS > 1e-12)) ||
                                sphereAtPrevJunction);
        sphereAtPrevJunction = false;
        if (std::getenv("MVB_CAP_DEBUG"))
            std::cerr << "[cap-debug] junction " << i << " angS=" << angS << " bentS=" << bentS
                      << " capAdopted=" << startCapAdopted << " welded=" << bridgedToPrev.back()
                      << "  '" << ptrs[i]->label << "'\n";
        if (diag) {
            built.push_back(solid);
            GProp_GProps gpr;
            BRepGProp::VolumeProperties(solid, gpr);
            Bnd_Box bb;
            BRepBndLib::Add(solid, bb);
            double x0, y0, z0, x1, y1, z1;
            bb.Get(x0, y0, z0, x1, y1, z1);
            const double len = primEndpoints(*ptrs[i]).first.Distance(primEndpoints(*ptrs[i]).second);
            const double nominal = kPi * wireRadius * wireRadius * len;
            std::cerr << "[prim] " << i << " '" << ptrs[i]->label << "' vol=" << gpr.Mass()
                      << " nominal=" << nominal
                      << " ratio=" << (nominal > 0 ? gpr.Mass() / nominal : 0.0)
                      << " y=[" << y0 << "," << y1 << "]\n";
        }
    }
    if (abutNoVerdict > 0 || diag) {
        std::cerr << "[abut] mitred joints with a boolean verdict: " << abutChecked
                  << ", without one: " << abutNoVerdict
                  << (abutNoVerdict > 0 ? "  <-- those corners are UNCHECKED, not clean\n" : "\n");
    }
    if (!overlappingJoints.empty()) {
        std::string msg =
            "ConductorBuilder: " + std::to_string(overlappingJoints.size()) +
            " mitred corner(s) did not abut -- the two sides overlap instead of sharing the "
            "bisector face:";
        for (const auto& j : overlappingJoints) msg += "\n  " + j;
        // ABT #685: the same escape the certified enamel gate carries (MVB_ALLOW_ENAMEL) -- OFF by
        // default, so anything that ships is a design whose corners abut. Set it only to EXPORT a
        // known-faulty design for visual review; the fault is still printed in full.
        if (std::getenv("MVB_ALLOW_ABUTMENT")) {
            std::cerr << "[corners] " << msg
                      << "\n  (MVB_ALLOW_ABUTMENT set: reported, not refused)\n";
        } else {
            throw std::runtime_error(msg);
        }
    }
    // Per-solid analysis (union-find connectivity, volumes, interpenetration) is O(n^2) EXACT
    // BRepExtrema/classification -- ~15k queries on a 170-prim toroid, tens of minutes. Keep it
    // behind its own flag so plain MVB_MITRE_DIAG (junction angles, DROPPED prims, the counters)
    // stays usable during normal iteration.
    if (diag && std::getenv("MVB_MITRE_INTERPEN")) {
        // component structure over ALL pairs (mirrors the test's union-find), with volumes
        std::vector<int> uf(built.size());
        for (size_t i = 0; i < uf.size(); ++i) uf[i] = (int)i;
        std::function<int(int)> find = [&](int a) {
            while (uf[a] != a) a = uf[a] = uf[uf[a]];
            return a;
        };
        for (size_t i = 0; i < built.size(); ++i)
            for (size_t j = i + 1; j < built.size(); ++j) {
                try {
                    BRepExtrema_DistShapeShape d(built[i], built[j]);
                    if (d.IsDone() && d.Value() < 1e-6) uf[find((int)i)] = find((int)j);
                } catch (const Standard_Failure&) {
                }
            }
        for (size_t i = 0; i < built.size(); ++i) {
            GProp_GProps gp;
            BRepGProp::VolumeProperties(built[i], gp);
            int nsol = 0;
            std::string sub;
            for (TopExp_Explorer ex(built[i], TopAbs_SOLID); ex.More(); ex.Next()) {
                ++nsol;
                GProp_GProps gs;
                BRepGProp::VolumeProperties(ex.Current(), gs);
                char vb[32]; std::snprintf(vb, sizeof(vb), " %.4g", gs.Mass()); sub += vb;
            }
            std::cerr << "[mitre]   solid " << i << " comp=" << find((int)i) << " vol=" << gp.Mass()
                      << " nsolids=" << nsol << (nsol > 1 ? " subvols:" + sub : "") << "\n";
        }
        for (size_t i = 1; i < built.size(); ++i) {
            const double dp = primEndpoints(*ptrs[i - 1]).second.Distance(primEndpoints(*ptrs[i]).first);
            double ds = -1.0;
            try {
                BRepExtrema_DistShapeShape ext(built[i - 1], built[i]);
                if (ext.IsDone()) ds = ext.Value();
            } catch (const Standard_Failure&) {
            }
            if (dp > 1e-9 || ds > 1e-9)
                std::cerr << "[mitre]   junction " << i - 1 << "->" << i << " endpointGap=" << dp
                          << " solidGap=" << ds << "\n";
        }
        // Interpenetration sweep over ALL bbox-touching pairs (not just chain neighbours), the
        // same junction-grid classifier as the test battery, with LABELS -- post-prune solid
        // indices in test output cannot be mapped back to primitives. O(n^2) solid
        // classification dominates the run (tens of minutes on a 170-prim toroid), so it needs
        // its OWN flag: plain MVB_MITRE_DIAG stays fast enough for routine use.
        {
            std::vector<Bnd_Box> bxs(built.size());
            for (size_t i = 0; i < built.size(); ++i) BRepBndLib::Add(built[i], bxs[i]);
            for (size_t i = 0; i < built.size(); ++i)
                for (size_t j = i + 1; j < built.size(); ++j) {
                    if (bxs[i].Distance(bxs[j]) > 1e-9) continue;
                    double x0, y0, z0, x1, y1, z1, u0, v0, w0, u1, v1, w1;
                    bxs[i].Get(x0, y0, z0, x1, y1, z1);
                    bxs[j].Get(u0, v0, w0, u1, v1, w1);
                    const double ax = std::max(x0, u0), bx = std::min(x1, u1);
                    const double ay = std::max(y0, v0), by = std::min(y1, v1);
                    const double az = std::max(z0, w0), bz = std::min(z1, w1);
                    if (ax >= bx || ay >= by || az >= bz) continue;
                    // classify per SOLID: a built[] entry can be a multi-solid compound (a trim
                    // that fragmented), and BRepClass3d_SolidClassifier on a compound silently
                    // misclassifies -- exploding is what makes this sweep agree with the test
                    // battery's flat-solid enumeration.
                    auto insideAny = [](const TopoDS_Shape& sh, const gp_Pnt& pp) {
                        for (TopExp_Explorer ex(sh, TopAbs_SOLID); ex.More(); ex.Next()) {
                            BRepClass3d_SolidClassifier c(TopoDS::Solid(ex.Current()), pp, 1e-9);
                            if (c.State() == TopAbs_IN) return true;
                        }
                        return false;
                    };
                    int hits = 0;
                    gp_Pnt hit;
                    constexpr int N = 6;
                    for (int gx = 0; gx < N && hits == 0; ++gx)
                        for (int gy = 0; gy < N && hits == 0; ++gy)
                            for (int gz = 0; gz < N && hits == 0; ++gz) {
                                gp_Pnt pp(ax + (bx - ax) * (gx + 0.5) / N,
                                          ay + (by - ay) * (gy + 0.5) / N,
                                          az + (bz - az) * (gz + 0.5) / N);
                                if (insideAny(built[i], pp) && insideAny(built[j], pp)) {
                                    ++hits;
                                    hit = pp;
                                }
                            }
                    if (hits > 0)
                        std::cerr << "[mitre]   INTERPENETRATION prims " << i << "<->" << j
                                  << " ['" << ptrs[i]->label << "' vs '" << ptrs[j]->label
                                  << "'] at (" << hit.X() << "," << hit.Y() << "," << hit.Z()
                                  << ")\n";
                }
        }
    }
    // ---- WELD THE BRIDGED NEIGHBOURS ---------------------------------------------------------
    // Mitred neighbours abut on a shared bisector face and stay separate solids. A BRIDGED
    // neighbour is different: it is left un-mitred ON PURPOSE, because a SPIRAL/BLEND pipe end
    // cannot be sliced (re-measured 2026-08-25 with MVB_MITRE_PIPES=1: on 10_emi the 3.47 deg
    // cut reported success, removed 0.000125 of 0.205 mm3 and left the overhang standing; on
    // 06_llc four mitred corners came out overlapping instead of abutting). So the wedge fill
    // POKES INTO the neighbour -- real copper inside copper, 523 instances across the corpus.
    //
    // Alf, 2026-08-25: "no overlap must be allowed in copper copper". Two pieces of ONE
    // conductor that overlap by construction are one piece of copper, so weld them and the
    // overlap ceases to exist. INCREMENTAL AND PAIRWISE, never one N-way union: a single
    // BRepAlgo fuse over a whole run is the most fragile step in this file (measured here: a
    // 135-piece run came back with 0 mm3). Each step is guarded on validity and volume, and a
    // step that fails simply ends the accumulation -- pieces are kept, never dropped.
    // MVB_ALLOW_WELD_LENS=1 restores the old overlapping pieces.
    {
        const bool keepLens = std::getenv("MVB_ALLOW_WELD_LENS") != nullptr;
        TopoDS_Shape acc;            // accumulated weld
        double accVol = 0.0;
        int accOwner = 0;
        // The UNWELDED pieces that went into `acc`, kept so a bad accumulation can be undone.
        std::vector<TopoDS_Shape> accPieces;
        auto flush = [&]() {
            if (acc.IsNull()) return;
            // SELF-INTERSECTION GATE, ONCE PER ACCUMULATION (2026-08-31).
            // A fuse can be BRepCheck-VALID yet BOPAlgo-SELF-INTERSECTING -- a check BRepCheck
            // simply does not perform. Such a solid is marginal: it survives in memory and
            // tips to invalid when the STEP round-trip reconstructs it, which is exactly how
            // 06_llc shipped 'Primary parallel 0 [solid 25]' as CAD DEFECTIVE (measured:
            // valid pre-scale, valid post-scale, invalid on read-back; no writer setting
            // changes it; ShapeFix and UnifySameDomain only "repair" it by deleting copper,
            // 3.652 -> 3.633/3.536 mm3, which the volume guards correctly reject).
            // Checking every fuse is correct but prohibitively slow (the analyzer re-runs on a
            // growing accumulator: ~6 turns in 25 min). Checking ONCE per flush costs one
            // analysis per conductor run instead of hundreds, and the remedy needs no repair:
            // fall back to the pieces we already hold, unwelded -- exactly the shape every
            // green design's weld-refusal path produces.
            if (accPieces.size() > 1) {
                bool selfInt = false;
                try {
                    BOPAlgo_ArgumentAnalyzer an;
                    an.SetShape1(acc);
                    an.ArgumentTypeMode() = Standard_True;
                    an.SelfInterMode() = Standard_True;
                    an.Perform();
                    selfInt = an.HasFaulty();
                } catch (const Standard_Failure&) {
                }
                if (selfInt) {
                    std::cerr << "[weld-selfint] '"
                              << (accOwner < (int)ptrs.size() ? ptrs[accOwner]->label
                                                              : std::string("?"))
                              << "' welded body self-intersects (BOPAlgo); emitting its "
                              << accPieces.size() << " pieces unwelded" << std::endl;
                    for (const auto& pc : accPieces) {
                        for (TopExp_Explorer px(pc, TopAbs_SOLID); px.More(); px.Next()) {
                            builder.Add(compound, px.Current());
                            if (primIndexPerSolid != nullptr)
                                primIndexPerSolid->push_back(accOwner);
                        }
                    }
                    acc.Nullify();
                    accVol = 0.0;
                    accPieces.clear();
                    return;
                }
            }
            // FINAL GATE: nothing invalid leaves the assembler. Every producer upstream
            // validates its own result, yet 06_llc still shipped one BRepCheck-invalid solid
            // ('Primary parallel 0 [solid 25]', 3.65 mm3, 80 faces) that was valid when swept
            // -- so some composition step degrades it. Repair here with ShapeFix (accepted
            // only if the volume is unchanged to 0.1%: representation may be cleaned, copper
            // may not move) and SAY SO if it cannot be repaired, rather than emitting a solid
            // the consumer will choke on.
            for (TopExp_Explorer ex(acc, TopAbs_SOLID); ex.More(); ex.Next()) {
                TopoDS_Shape out = ex.Current();
                if (!BRepCheck_Analyzer(out).IsValid()) {
                    GProp_GProps g0;
                    BRepGProp::VolumeProperties(out, g0);
                    bool repaired = false;
                    try {
                        Handle(ShapeFix_Shape) fix = new ShapeFix_Shape(out);
                        fix->Perform();
                        const TopoDS_Shape fixed = fix->Shape();
                        if (!fixed.IsNull() && BRepCheck_Analyzer(fixed).IsValid()) {
                            GProp_GProps g1;
                            BRepGProp::VolumeProperties(fixed, g1);
                            if (g0.Mass() > 0 &&
                                std::abs(g1.Mass() - g0.Mass()) <= 1e-3 * g0.Mass()) {
                                out = fixed;
                                repaired = true;
                                std::cerr << "[assembler-fix] repaired an invalid solid of '"
                                          << (accOwner < (int)ptrs.size() ? ptrs[accOwner]->label
                                                                          : std::string("?"))
                                          << "' (" << g0.Mass() * 1e9 << " mm^3)" << std::endl;
                            }
                        }
                    } catch (const Standard_Failure&) {
                    }
                    if (!repaired)
                        std::cerr << "[assembler-invalid] '"
                                  << (accOwner < (int)ptrs.size() ? ptrs[accOwner]->label
                                                                  : std::string("?"))
                                  << "' emits an INVALID solid (" << g0.Mass() * 1e9
                                  << " mm^3) that ShapeFix could not repair" << std::endl;
                }
                builder.Add(compound, out);
                if (primIndexPerSolid != nullptr) primIndexPerSolid->push_back(accOwner);
            }
            acc.Nullify();
            accVol = 0.0;
            accPieces.clear();
        };
        for (size_t k = 0; k < perPrimSolids.size(); ++k) {
            // One primitive may have been fragmented into several solids by its own mitre trim;
            // treat that group as a unit.
            TopoDS_Shape group;
            if (perPrimSolids[k].size() == 1) {
                group = perPrimSolids[k].front();
            } else {
                TopoDS_Compound c; BRep_Builder gb; gb.MakeCompound(c);
                for (const auto& sh : perPrimSolids[k]) gb.Add(c, sh);
                group = c;
            }
            double groupVol = 0.0;
            { GProp_GProps g; BRepGProp::VolumeProperties(group, g); groupVol = g.Mass(); }

            const bool weldToPrev = !keepLens && !acc.IsNull() &&
                                    k < bridgedToPrev.size() && bridgedToPrev[k];
            if (!weldToPrev) {
                flush();
                acc = group; accVol = groupVol; accOwner = (int)k;
                accPieces.assign(1, group);
                continue;
            }
            bool welded = false;
            // TOLERANT-IMPRINT LADDER. Standard assembly-meshing practice (Sandia/Coreform Cubit,
            // "Resolving Problems with Conformal Assemblies") is imprint-and-merge, falling back
            // to a TOLERANT imprint when the exact one fails on sub-tolerance discrepancies --
            // "normal imprinting when possible, tolerant imprinting only when normal imprinting
            // fails". A fuzzy boolean value is OCC's tolerant imprint. Exact-first, then the same
            // 1e-7 m ConductorBuilder already uses to weld its chain primitives, then 1e-6.
            // Each fuzzy rung is tried in BOTH operand orders. Measured on 10_emi (EP 13):
            // fusing a 4%-of-chain bump riser as the TOOL onto the accumulated chain drops the
            // accumulator at fuzzy=0 (result = the riser alone) and invents +0.29% volume at
            // fuzzy=1e-7 -- an operand-order pathology in the boolean, not a geometry fault.
            for (bool swapped : {false, true}) {
            for (double fuzzy : {0.0, 1e-7, 1e-6}) {
            if (welded) break;
            try {
                BRepAlgoAPI_Fuse fu;
                TopTools_ListOfShape fargs, ftools;
                if (swapped) { fargs.Append(group); ftools.Append(acc); }
                else         { fargs.Append(acc); ftools.Append(group); }
                fu.SetArguments(fargs); fu.SetTools(ftools);
                if (fuzzy > 0) fu.SetFuzzyValue(fuzzy);
                fu.Build();
                if (std::getenv("MVB_WELD_DEBUG")) {
                    std::cerr << "[weld-debug] '" << ptrs[k]->label << "' swapped=" << swapped
                              << " fuzzy=" << fuzzy
                              << " done=" << fu.IsDone()
                              << " null=" << (fu.IsDone() ? fu.Shape().IsNull() : true);
                    if (fu.IsDone() && !fu.Shape().IsNull()) {
                        const bool valid = BRepCheck_Analyzer(fu.Shape()).IsValid();
                        GProp_GProps gd; BRepGProp::VolumeProperties(fu.Shape(), gd);
                        std::cerr << " valid=" << valid << " vol=" << gd.Mass()
                                  << " acc=" << accVol << " group=" << groupVol
                                  << " sum=" << (accVol + groupVol);
                    }
                    std::cerr << "\n";
                }
                if (fu.IsDone() && !fu.Shape().IsNull() &&
                    BRepCheck_Analyzer(fu.Shape()).IsValid()) {
                    GProp_GProps gf; BRepGProp::VolumeProperties(fu.Shape(), gf);
                    const double sum = accVol + groupVol;
                    // ONLY THE TWO BOUNDS THAT ARE EXACTLY TRUE OF A UNION: it cannot be smaller
                    // than its largest input (that would mean OCC dropped an operand -- the failure
                    // this guard exists for) and cannot exceed their sum (that would be invented
                    // material). There used to be a third clause, >= 0.90 * sum, which was a guess
                    // at how much overlap is plausible. It is wrong for the joint sphere: half a
                    // ball sits inside each adjoining tube, so against a SMALL neighbour the true
                    // union legitimately falls below 90% of the sum and a perfectly good fuse was
                    // rejected -- leaving the sphere overlapping the next piece instead of welded
                    // into it, which is why chunk 7 measured WORSE than the growth it replaces
                    // (28 self-overlaps against 19 on 14_dab).
                    // Tolerance 1e-3, MEASURED on 10_emi (EP 13): legitimate welds carry up to
                    // +3.9e-5 relative volume noise at fuzzy=0 (quadrature on bspline junction
                    // faces) and +0.29% at fuzzy=1e-7 (the tolerant imprint displaces the
                    // junction boundary by <=0.1 um) -- the old 1e-6 bound refused them all,
                    // leaving real 1.8e-4 mm^3 lens overlaps that fail the audit AND poison
                    // tetgen (04_forward's PLC error sat exactly on a refused junction). The
                    // failures this bound exists for -- OCC dropping an operand -- are >=2%
                    // errors (a dropped bump riser is 4% of its chain), an order of magnitude
                    // beyond the noise. 1e-3 sits in the gap.
                    if (gf.Mass() <= sum * (1.0 + 1e-3) &&
                        gf.Mass() >= std::max(accVol, groupVol) * (1.0 - 1e-3)) {
                        // SELF-INTERSECTION ACCEPTANCE (MVB_WELD_SELFINT_CHECK, 2026-08-31).
                        // BRepCheck_Analyzer passes a fused solid that BOPAlgo's own argument
                        // analyzer calls self-intersecting, and such a solid is MARGINAL: it
                        // survives in memory but the STEP round-trip tips it into
                        // BRepCheck-invalid (measured on 06_llc: 'Primary parallel 0
                        // [solid 25]' reads back with SelfIntersectingWire + UnorientableShape
                        // + BadOrientationOfSubshape, while PRE- and POST-scale in memory it
                        // is valid; no writer setting and no post-hoc repair fixes it --
                        // ShapeFix and UnifySameDomain only "repair" it by deleting copper).
                        // Rejecting the fuse here sends the junction to the cut/abut fallback,
                        // which is the configuration every green design already uses.
                        // Env-gated until measured against the green corpus: the check is not
                        // free (BOPAlgo on a growing accumulator).
                        static const bool selfIntCheck =
                            std::getenv("MVB_WELD_SELFINT_CHECK") != nullptr;
                        bool fusedSelfIntersects = false;
                        if (selfIntCheck) {
                            try {
                                BOPAlgo_ArgumentAnalyzer an;
                                an.SetShape1(fu.Shape());
                                an.ArgumentTypeMode() = Standard_True;
                                an.SelfInterMode() = Standard_True;
                                an.Perform();
                                fusedSelfIntersects = an.HasFaulty();
                            } catch (const Standard_Failure&) {
                                fusedSelfIntersects = true;   // cannot prove it clean: refuse
                            }
                            if (fusedSelfIntersects) {
                                std::cerr << "[weld-selfint] '" << ptrs[k]->label
                                          << "' fuse is self-intersecting; refusing it"
                                          << std::endl;
                            }
                        }
                        if (fusedSelfIntersects) break;   // fall through to the cut/abut path
                        acc = fu.Shape(); accVol = gf.Mass(); welded = true;
                        accPieces.push_back(group);   // undo material for the flush gate
                    }
                }
            } catch (const Standard_Failure&) {
            }
            }
            if (welded) break;
            }
            if (!welded) {
                // FUSE-REFUSED FALLBACK: cut the junction lens off the incoming piece so it
                // ABUTS the chain exactly. Bump-riser junctions are boolean-sick for Fuse in
                // BOTH operand orders and at every fuzzy rung (measured on 10_emi: fuzzy=0
                // drops the accumulator, fuzzy=1e-7 invents 0.29%, swapped collapses to 0) --
                // but Cut(group, acc) only has to trim a <=e-4 mm^3 lens off one small piece.
                // Exact abutment of separate solids is a shape the pipeline already proves
                // meshable (the toroid chain is exactly that), and the audit sees no overlap.
                // Acceptance: the cut may only REMOVE the lens -- result within [90%, 100%] of
                // the piece (losing >10% means the junction is not a lens; refuse loudly).
                bool trimmed = false;
                // The lens involves only the immediately preceding primitive, so cut against
                // THAT, not the whole accumulated chain -- a 200-solid accumulator makes the
                // boolean as sick as the fuse it is standing in for (measured: cut-vs-acc
                // failed on all five 10_emi bump risers, cut-vs-predecessor is two small
                // convex-ish operands).
                try {
                    TopoDS_Shape prevTool;
                    if (k > 0 && !perPrimSolids[k - 1].empty()) {
                        if (perPrimSolids[k - 1].size() == 1) {
                            prevTool = perPrimSolids[k - 1].front();
                        } else {
                            TopoDS_Compound pc; BRep_Builder pb; pb.MakeCompound(pc);
                            for (const auto& sh : perPrimSolids[k - 1]) pb.Add(pc, sh);
                            prevTool = pc;
                        }
                    }
                    // Measure the LENS directly with a Common: its volume integrates on the
                    // lens itself, so the error is relative to the lens -- inferring it from
                    // groupVol - cutVol drowns a real 2.8e-5 mm^3 lens in the group's own
                    // +-1.6e-4 mm^3 quadrature noise (measured, 10_emi turn-30).
                    double lensVol = -1.0;
                    if (!prevTool.IsNull()) {
                        try {
                            BRepAlgoAPI_Common comOp;
                            TopTools_ListOfShape mArgs, mTools;
                            mArgs.Append(group); mTools.Append(prevTool);
                            comOp.SetArguments(mArgs); comOp.SetTools(mTools);
                            comOp.Build();
                            if (comOp.IsDone() && !comOp.Shape().IsNull()) {
                                GProp_GProps gl; BRepGProp::VolumeProperties(comOp.Shape(), gl);
                                lensVol = std::max(0.0, gl.Mass());
                            }
                        } catch (const Standard_Failure&) {
                        }
                    }
                    // The Common's testimony is only trusted when it produced SOLIDS or the cut
                    // agrees: on the boolean-sick riser junctions Common returns EMPTY for a
                    // piece that is ENTIRELY inside its predecessor (measured on 10_emi -- the
                    // lens-floor version shipped the risers fully overlapping). So: run the CUT
                    // ladder first; absorb on an empty remainder; and only skip the trim when
                    // Common POSITIVELY measured a sub-tolerance lens.
                    if (!prevTool.IsNull()) {
                        for (double cfuzzy : {0.0, 1e-7, 1e-6}) {
                        if (trimmed) break;
                        BRepAlgoAPI_Cut cutOp;
                        TopTools_ListOfShape cargs, ctools;
                        cargs.Append(group); ctools.Append(prevTool);
                        cutOp.SetArguments(cargs); cutOp.SetTools(ctools);
                        if (cfuzzy > 0) cutOp.SetFuzzyValue(cfuzzy);
                        cutOp.Build();
                        if (std::getenv("MVB_WELD_DEBUG")) {
                            std::cerr << "[weld-debug] cut '" << ptrs[k]->label
                                      << "' cfuzzy=" << cfuzzy << " done=" << cutOp.IsDone();
                            if (cutOp.IsDone() && !cutOp.Shape().IsNull()) {
                                GProp_GProps gd; BRepGProp::VolumeProperties(cutOp.Shape(), gd);
                                std::cerr << " valid=" << BRepCheck_Analyzer(cutOp.Shape()).IsValid()
                                          << " vol=" << gd.Mass() << " group=" << groupVol;
                            }
                            std::cerr << "\n";
                        }
                        if (!cutOp.IsDone() || cutOp.Shape().IsNull()) continue;
                        GProp_GProps gc; BRepGProp::VolumeProperties(cutOp.Shape(), gc);
                        // Piece ENTIRELY inside its predecessor (measured: 4 of 5 10_emi bump
                        // risers cut to ~0 -- the dragback piece already covers the riser's
                        // span). A ∪ B = A when B ⊆ A: absorbing the redundant piece IS the
                        // weld. No validity demand on the empty remainder.
                        if (gc.Mass() <= groupVol * 0.01) {
                            std::cerr << "[weld-cut] '" << ptrs[k]->label
                                      << "' lies entirely inside its predecessor; absorbed"
                                      << std::endl;
                            trimmed = true;   // acc unchanged: the union is the predecessor
                            break;
                        }
                        // LENS FLOOR: when the Common POSITIVELY measured the lens below the
                        // audit tolerance (1e-6 mm^3 = 1e-15 m^3), the pieces already abut --
                        // do NOT cut. Measured on 04_forward: cutting e-14 mm^3 noise lenses
                        // imprinted micro-edges (a 40 um boundary pocket) whose curves CROSS
                        // after the fragment, and gmsh's 1D-intersection recovery then
                        // retry-loops past any budget. Requires lensVol >= 0 (Common produced
                        // a result) AND the cut to be sane (remainder ~ the whole piece), so a
                        // lying empty Common cannot bless a real overlap.
                        if (lensVol >= 0.0 && lensVol < 1e-15 &&
                            gc.Mass() >= groupVol * (1.0 - 1e-3)) {
                            std::cerr << "[weld-cut] '" << ptrs[k]->label
                                      << "' junction lens below audit tolerance ("
                                      << lensVol * 1e9 << " mm^3); left abutting, uncut"
                                      << std::endl;
                            flush();
                            acc = group; accVol = groupVol; accOwner = (int)k;
                            accPieces.assign(1, group);
                            trimmed = true;
                            break;
                        }
                        // Upper bound 1e-3, not 1e-9: a CUT cannot add material, so any
                        // apparent excess is quadrature noise (measured +5e-4 relative on the
                        // 10_emi turn-30 wrap junction, where the true lens is only -0.9e-4 --
                        // the noise swamps the lens and a 1e-9 bound rejects a correct cut).
                        TopoDS_Shape cutRes = cutOp.Shape();
                        bool cutValid = BRepCheck_Analyzer(cutRes).IsValid();
                        if (!cutValid) {
                            // REPAIR, THEN RE-JUDGE. On 06_llc's dragback junctions the cut
                            // completes but returns a solid BRepCheck rejects at every fuzzy
                            // rung (measured: 7 of its 27 junctions), so the lens survived into
                            // the STEP -- 1 invalid solid, ~100 sliver faces, CAD DEFECTIVE.
                            // ShapeFix is OCC's own repair for exactly this (boolean debris:
                            // wire order, tiny faces, tolerances); the result is then held to
                            // the SAME volume bounds and validity test as an unrepaired cut, so
                            // nothing unverified can ship.
                            try {
                                Handle(ShapeFix_Shape) fix = new ShapeFix_Shape(cutRes);
                                fix->Perform();
                                const TopoDS_Shape fixed = fix->Shape();
                                if (!fixed.IsNull() && BRepCheck_Analyzer(fixed).IsValid()) {
                                    GProp_GProps gf2;
                                    BRepGProp::VolumeProperties(fixed, gf2);
                                    if (gf2.Mass() <= groupVol * (1.0 + 1e-3) &&
                                        gf2.Mass() >= groupVol * 0.90) {
                                        cutRes = fixed;
                                        gc = gf2;
                                        cutValid = true;
                                        std::cerr << "[weld-cut] '" << ptrs[k]->label
                                                  << "' cut was invalid; ShapeFix repaired it"
                                                  << std::endl;
                                    }
                                }
                            } catch (const Standard_Failure&) {
                            }
                        }
                        if (cutValid &&
                            gc.Mass() <= groupVol * (1.0 + 1e-3) &&
                            gc.Mass() >= groupVol * 0.90) {
                            std::cerr << "[weld-cut] '" << ptrs[k]->label
                                      << "' would not fuse; trimmed its junction lens ("
                                      << (groupVol - gc.Mass()) * 1e9
                                      << " mm^3) so it abuts its predecessor exactly"
                                      << std::endl;
                            flush();
                            acc = cutRes; accVol = gc.Mass(); accOwner = (int)k;
                            accPieces.assign(1, cutRes);
                            trimmed = true;
                        }
                        }
                    }
                } catch (const Standard_Failure&) {
                }
                if (!trimmed) {
                    // Keep both, unwelded, and SAY SO: this is a copper-copper overlap that
                    // will reach the audit, not something to swallow quietly.
                    std::cerr << "[weld-lens] '" << ptrs[k]->label
                              << "' would not weld onto its bridged predecessor; the junction "
                                 "overlap remains" << std::endl;
                    flush();
                    acc = group; accVol = groupVol; accOwner = (int)k;
                    accPieces.assign(1, group);
                }
            }
        }
        flush();
    }

    if (diag) {
        std::cerr << "[mitre] prims=" << n << " boolean-cuts=" << nCut << " repaired=" << nRepaired
                  << " dropped-invalid=" << nInvalid << "\n";
    }
    return compound;
}

} // namespace mvb
