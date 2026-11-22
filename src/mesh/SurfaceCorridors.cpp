#include "mvb/mesh/SurfaceCorridors.h"

#include "mvb/mesh/SizeField.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <Bnd_Box.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <IntCurvesFace_Intersector.hxx>
#include <Precision.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Dir.hxx>
#include <gp_Lin.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>

namespace mvb::mesh {

// ---- CORRIDORS OF RECTANGULAR CONDUCTORS, MEASURED ON THEIR BUILT SURFACES (ABT #2460) ------------
// The chord bound is the same for every conductor: a facet of size h on a surface whose radius of curvature is rho
// sags h^2/8rho off it, so a gap d needs h <= 0.8 sqrt(4 rho d). What differs is where rho and d come from. A round
// wire carries them in its centreline and its radius (chord_points_3d in Mesh3d.cpp). A rectangular wire or a foil
// does not: its section has no radius, its broad faces are flat, and the only curvature it has is the one its PATH
// gives its faces
// (the wrap radius about the column, a bend at a column corner, a lead's own round section, a compressed bundle's
// corner fillet). So for these conductors both are read off the solids that will be meshed:
//   (i)  every CURVED face of a rectangular conductor is sampled; at each sample the radius is the face's own
//        (the larger principal curvature), and the gap is what the outward normal meets first on ANY conductor
//        solid, its own included (the next turn of the same foil). A flat face is never sampled: it has no sag.
//   (ii) every centreline sample of a ROUND wire is measured against the faces of the rectangular conductors near
//        it; the radius is the wire's.
// A pair with a curved rectangular face on one side and a round wire on the other is seen from both sides and keeps
// the finer bound. What (i) does not see: two rectangular conductors that only approach corner to corner, with no
// face of either under the other's normal.
SurfaceCorridors surface_corridor_points(const std::vector<ConductorSolid>& solids,
                                         const std::vector<RoundCentreline>& roundCentrelines) {
    SurfaceCorridors out;
    const auto t0 = std::chrono::steady_clock::now();
    for (const auto& rc : roundCentrelines)
        if (!(rc.radius > 0.0) || !(rc.target > 0.0))
            throw std::runtime_error("surface corridors: round conductor '" + rc.name + "' has no positive radius and target");
    double tmax = 0.0;
    for (const auto& cs : solids) {
        if (!(cs.target > 0.0))
            throw std::runtime_error("surface corridors: conductor solid of winding '" + cs.winding + "' has no positive target");
        tmax = std::max(tmax, cs.target);
        if (cs.rectangular) ++out.rectangularSolids;
    }
    if (out.rectangularSolids == 0) return out;
    const double touching = 1e-7;   // a gap this small is contact: the layout check owns it (as in chord_points_3d)

    // every face of every conductor solid, with its box; the faces of rectangular solids are also what (ii) measures to
    struct FaceRec {
        TopoDS_Face face; int solid = -1; double lo[3], hi[3];
        std::unique_ptr<IntCurvesFace_Intersector> ray;
        std::unique_ptr<BRepAdaptor_Surface> surf;
    };
    std::vector<FaceRec> faces;
    double glo[3] = {1e300, 1e300, 1e300}, ghi[3] = {-1e300, -1e300, -1e300};
    for (size_t si = 0; si < solids.size(); ++si)
        for (TopExp_Explorer ex(solids[si].shape, TopAbs_FACE); ex.More(); ex.Next()) {
            FaceRec fr; fr.face = TopoDS::Face(ex.Current()); fr.solid = static_cast<int>(si);
            Bnd_Box b; BRepBndLib::Add(fr.face, b, /*useTriangulation=*/false);
            if (b.IsVoid()) continue;
            b.SetGap(0.0);
            b.Get(fr.lo[0], fr.lo[1], fr.lo[2], fr.hi[0], fr.hi[1], fr.hi[2]);
            for (int d = 0; d < 3; ++d) { glo[d] = std::min(glo[d], fr.lo[d]); ghi[d] = std::max(ghi[d], fr.hi[d]); }
            faces.push_back(std::move(fr));
        }
    if (faces.empty()) throw std::runtime_error("surface corridors: the conductor solids have no faces");
    auto surfOf = [&](FaceRec& fr) -> BRepAdaptor_Surface& {
        if (!fr.surf) fr.surf = std::make_unique<BRepAdaptor_Surface>(fr.face, Standard_True);
        return *fr.surf;
    };
    // outward normal and radius of curvature of a face at (u, v); false where the surface has neither
    auto surfaceAt = [&](FaceRec& fr, double u, double v, gp_Pnt& pnt, gp_Dir& normal, double& curvature) {
        BRepLProp_SLProps props(surfOf(fr), u, v, 2, Precision::Confusion());
        if (!props.IsNormalDefined()) return false;
        pnt = props.Value();
        normal = props.Normal();
        if (fr.face.Orientation() == TopAbs_REVERSED) normal.Reverse();
        curvature = props.IsCurvatureDefined() ? std::max(std::fabs(props.MaxCurvature()), std::fabs(props.MinCurvature())) : 0.0;
        return true;
    };

    // how far anything is ever looked for from a point of a rectangular conductor's surface: the widest binding gap of
    // any radius under the coarsest target (the two limits of widest_binding_gap meet at rho = sqrt(0.08) t / 0.8)
    const double rhoWidest = std::sqrt(0.08) * tmax / 0.8;
    const double reach = widest_binding_gap(rhoWidest, tmax);
    // ... and from a round wire's centreline: its radius plus its own widest binding gap
    double reachRound = 0.0;
    for (const auto& rc : roundCentrelines)
        reachRound = std::max(reachRound, rc.radius + widest_binding_gap(rc.radius, std::max(rc.target, tmax)));
    // WHICH FACES CAN BE NEAR A POINT: a cloud of points over every face of every conductor solid, one per cell of a
    // parameter grid no coarser than half the solid's own target (cells are not trimmed to the face: the cloud is a
    // superset). A surface point within d of some point P has its cell's centre within d + that cell's half-diagonal of
    // P, so the faces worth an exact measurement are those with a cloud point that close. Boxes cannot do this job: the
    // box of one turn's helical face is the whole ring, and every sample of 03_buck was ray-tested against every turn
    // (213 s, for no hit at all).
    struct CloudPoint { double x, y, z; int face; double u, v; };
    struct FaceGrid { double du = 0.0, dv = 0.0, lenU = 0.0, lenV = 0.0; };   // a cell's parameter and physical sizes
    std::vector<FaceGrid> faceGrid(faces.size());
    std::vector<CloudPoint> cloud;
    std::vector<double> cellRadius(faces.size(), 0.0);   // each face's cell half-diagonal [m]
    std::vector<double> tightest(faces.size(), 0.0);     // each face's largest probed curvature [1/m]
    double cellRadiusMax = 0.0;
    for (size_t fi = 0; fi < faces.size(); ++fi) {
        BRepAdaptor_Surface& S = surfOf(faces[fi]);
        const double u0 = S.FirstUParameter(), u1 = S.LastUParameter(), v0 = S.FirstVParameter(), v1 = S.LastVParameter();
        if (!(u1 > u0) || !(v1 > v0) || !std::isfinite(u0) || !std::isfinite(u1) || !std::isfinite(v0) || !std::isfinite(v1))
            throw std::runtime_error("surface corridors: a conductor face of winding '" + solids[faces[fi].solid].winding +
                                     "' has no finite parameter range");
        double lenU = 0.0, lenV = 0.0, kmax = 0.0;
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
            const double u = u0 + (a + 0.5) / 3.0 * (u1 - u0), v = v0 + (b + 0.5) / 3.0 * (v1 - v0);
            gp_Pnt p; gp_Vec du, dv; S.D1(u, v, p, du, dv);
            lenU = std::max(lenU, du.Magnitude() * (u1 - u0)); lenV = std::max(lenV, dv.Magnitude() * (v1 - v0));
            gp_Pnt q; gp_Dir n(0, 0, 1); double k = 0.0;
            if (S.GetType() != GeomAbs_Plane && surfaceAt(faces[fi], u, v, q, n, k)) kmax = std::max(kmax, k);
        }
        tightest[fi] = kmax;
        // the probe is nine points: a quarter more length than it saw, so a cell is not larger than it is taken for
        lenU *= 1.25; lenV *= 1.25;
        const double step = 0.5 * solids[faces[fi].solid].target;
        const int nu = std::max(1, static_cast<int>(std::ceil(lenU / step))), nv = std::max(1, static_cast<int>(std::ceil(lenV / step)));
        faceGrid[fi] = {(u1 - u0) / nu, (v1 - v0) / nv, lenU / nu, lenV / nv};
        cellRadius[fi] = 0.5 * std::hypot(lenU / nu, lenV / nv);
        cellRadiusMax = std::max(cellRadiusMax, cellRadius[fi]);
        for (int a = 0; a < nu; ++a) for (int b = 0; b < nv; ++b) {
            const double u = u0 + (a + 0.5) / nu * (u1 - u0), v = v0 + (b + 0.5) / nv * (v1 - v0);
            const gp_Pnt p = S.Value(u, v);
            cloud.push_back({p.X(), p.Y(), p.Z(), static_cast<int>(fi), u, v});
        }
    }
    const double cell = cellRadiusMax + std::max(reach, reachRound);
    auto cellKey = [](long long i, long long j, long long k) {
        return (i & 0x1FFFFF) | ((j & 0x1FFFFF) << 21) | ((k & 0x1FFFFF) << 42);
    };
    auto cellIndex = [&](double x, int d) { return static_cast<long long>(std::floor((x - glo[d]) / cell)); };
    std::unordered_map<long long, std::vector<int>> grid;
    for (size_t ci = 0; ci < cloud.size(); ++ci)
        grid[cellKey(cellIndex(cloud[ci].x, 0), cellIndex(cloud[ci].y, 1), cellIndex(cloud[ci].z, 2))].push_back(static_cast<int>(ci));
    // Whether the part of a face under one cloud cell can come within `within` of p: the cell is halved, and halved
    // again, a quarter dropped as soon as its centre is further than `within` plus its own half-diagonal, until a piece
    // smaller than `within` is still in reach (then the exact measurement decides). This is what keeps the exact tests
    // for the faces that matter: on 03_buck the neighbouring turn is 0.1 mm away, inside a cloud cell's own radius but
    // twenty times the 5 um its wrap can bind, and 93,644 exact ray tests found nothing in 34 s.
    std::function<bool(int, double, double, double, double, const gp_Pnt&, double)> cellReaches =
        [&](int fi, double u, double v, double du, double dv, const gp_Pnt& p, double within) -> bool {
        const FaceGrid& fg = faceGrid[fi];
        const double radius = 0.5 * std::hypot(fg.lenU * du / fg.du, fg.lenV * dv / fg.dv);
        const double d = surfOf(faces[fi]).Value(u, v).Distance(p);
        if (d > within + radius) return false;
        if (radius < within || d <= within) return true;
        for (int a = -1; a <= 1; a += 2) for (int b = -1; b <= 1; b += 2)
            if (cellReaches(fi, u + 0.25 * a * du, v + 0.25 * b * dv, 0.5 * du, 0.5 * dv, p, within)) return true;
        return false;
    };
    // the faces that can come within `within` of p
    std::vector<int> stamp(faces.size(), -1), near;
    int query = 0;
    auto facesNear = [&](const gp_Pnt& p, double within) -> const std::vector<int>& {
        near.clear(); ++query;
        const double R = within + cellRadiusMax;
        long long lo[3], hi[3];
        const double q[3] = {p.X(), p.Y(), p.Z()};
        for (int d = 0; d < 3; ++d) { lo[d] = cellIndex(q[d] - R, d); hi[d] = cellIndex(q[d] + R, d); }
        for (long long i = lo[0]; i <= hi[0]; ++i) for (long long j = lo[1]; j <= hi[1]; ++j) for (long long k = lo[2]; k <= hi[2]; ++k) {
            const auto it = grid.find(cellKey(i, j, k));
            if (it == grid.end()) continue;
            for (int ci : it->second) {
                const CloudPoint& c = cloud[ci];
                if (stamp[c.face] == query) continue;
                const double reachHere = within + cellRadius[c.face];
                const double dx = c.x - q[0], dy = c.y - q[1], dz = c.z - q[2];
                if (dx * dx + dy * dy + dz * dz > reachHere * reachHere) continue;
                if (!cellReaches(c.face, c.u, c.v, faceGrid[c.face].du, faceGrid[c.face].dv, p, within)) continue;
                stamp[c.face] = query; near.push_back(c.face);
            }
        }
        return near;
    };
    auto boxDistance = [](const FaceRec& fr, const gp_Pnt& p) {
        const double q[3] = {p.X(), p.Y(), p.Z()};
        double s2 = 0.0;
        for (int d = 0; d < 3; ++d) { const double e = std::max({fr.lo[d] - q[d], 0.0, q[d] - fr.hi[d]}); s2 += e * e; }
        return std::sqrt(s2);
    };

    // (i) the curved faces of the rectangular conductors
    const size_t nfaces = faces.size();
    for (size_t fi = 0; fi < nfaces; ++fi) {
        if (!solids[faces[fi].solid].rectangular) continue;
        BRepAdaptor_Surface& S = surfOf(faces[fi]);
        if (S.GetType() == GeomAbs_Plane) continue;
        const double tOwn = solids[faces[fi].solid].target;
        const double u0 = S.FirstUParameter(), u1 = S.LastUParameter(), v0 = S.FirstVParameter(), v1 = S.LastVParameter();
        if (!(u1 > u0) || !(v1 > v0) || !std::isfinite(u0) || !std::isfinite(u1) || !std::isfinite(v0) || !std::isfinite(v1)) continue;
        // sample spacing: half the smaller of the winding's target and the face's own radius (the centreline model's
        // half-a-radius), from the probe of the face
        const double kmax = tightest[fi];
        double lenU = 0.0, lenV = 0.0;
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
            const double u = u0 + (a + 0.5) / 3.0 * (u1 - u0), v = v0 + (b + 0.5) / 3.0 * (v1 - v0);
            gp_Pnt p; gp_Vec du, dv; S.D1(u, v, p, du, dv);
            lenU = std::max(lenU, du.Magnitude() * (u1 - u0)); lenV = std::max(lenV, dv.Magnitude() * (v1 - v0));
        }
        if (!(kmax > 0.0)) continue;   // flat wherever it was probed
        ++out.curvedFaces;
        const double spacing = 0.5 * std::min(tOwn, 1.0 / kmax);
        const int nu = std::max(1, static_cast<int>(std::ceil(lenU / spacing))), nv = std::max(1, static_cast<int>(std::ceil(lenV / spacing)));
        BRepTopAdaptor_FClass2d inFace(faces[fi].face, Precision::PConfusion());
        for (int a = 0; a < nu; ++a) for (int b = 0; b < nv; ++b) {
            const double u = u0 + (a + 0.5) / nu * (u1 - u0), v = v0 + (b + 0.5) / nv * (v1 - v0);
            const TopAbs_State st = inFace.Perform(gp_Pnt2d(u, v));
            if (st != TopAbs_IN && st != TopAbs_ON) continue;
            gp_Pnt P; gp_Dir n(0, 0, 1); double k = 0.0;
            if (!surfaceAt(faces[fi], u, v, P, n, k) || !(k > 0.0)) continue;
            const double rho = 1.0 / k;
            const double dmax = widest_binding_gap(rho, tmax);
            if (!(dmax > touching)) { ++out.unbound; continue; }   // too gently curved to bind any gap wider than contact
            ++out.samples;
            const gp_Lin line(P, n);
            double best = dmax; int bestFace = -1; double bu = 0.0, bv = 0.0;
            for (int ci : facesNear(P, dmax)) {
                // its own face again: only where the face curves back on itself within the reach (an outward normal
                // leaves a surface of radius rho for at least 2 rho; half of that, the probe being nine points)
                if (ci == static_cast<int>(fi) && dmax * kmax < 1.0) continue;
                FaceRec& g = faces[ci];
                if (boxDistance(g, P) > best) continue;
                if (!g.ray) g.ray = std::make_unique<IntCurvesFace_Intersector>(g.face, Precision::Confusion());
                g.ray->Perform(line, touching, best);
                if (!g.ray->IsDone()) continue;
                for (int m = 1; m <= g.ray->NbPnt(); ++m) {
                    if (g.ray->State(m) != TopAbs_IN && g.ray->State(m) != TopAbs_ON) continue;
                    const double w = g.ray->WParameter(m);
                    if (w > touching && w < best) { best = w; bestFace = ci; bu = g.ray->UParameter(m); bv = g.ray->VParameter(m); }
                }
            }
            if (bestFace < 0) continue;
            ++out.hits;
            const double tPair = std::max(tOwn, solids[faces[bestFace].solid].target);
            gp_Pnt Q; gp_Dir nq(0, 0, 1); double kq = 0.0;
            if (!surfaceAt(faces[bestFace], bu, bv, Q, nq, kq)) continue;
            if (!surfaces_face(n.Dot(nq), rho, tPair)) { ++out.inclined; continue; }
            if (best >= widest_binding_gap(rho, tPair)) continue;
            const double h = chord_bound(rho, best);
            if (!corridor_binds(h, tOwn, solids[faces[bestFace].solid].target)) continue;
            if (faces[bestFace].solid == faces[fi].solid) ++out.ownHits;
            CorridorCandidate c; c.h = h; c.gap = best;
            c.p[0] = P.X() + 0.5 * best * n.X(); c.p[1] = P.Y() + 0.5 * best * n.Y(); c.p[2] = P.Z() + 0.5 * best * n.Z();
            out.points.push_back(c);
        }
    }

    // (ii) round centrelines against the faces of the rectangular conductors
    auto partOf = [](const RoundCentreline& rc, const ConductorSolid& cs) {
        return !cs.conductor.empty() && (rc.name == cs.conductor || rc.name.rfind(cs.conductor + " ", 0) == 0);
    };
    for (const auto& pl : roundCentrelines) {
        const double tRound = pl.target;
        const double r = pl.radius;
        // the solids this centreline is the axis of a part of (a foil's own lead): its distance to them is its radius
        std::vector<char> own(solids.size(), 0);
        for (size_t si = 0; si < solids.size(); ++si) own[si] = partOf(pl, solids[si]) ? 1 : 0;
        for (const auto& prim : pl.pieces)
            for (size_t s = 0; s + 1 < prim.size(); ++s) {
                const gp_Pnt A(prim[s][0], prim[s][1], prim[s][2]), B(prim[s + 1][0], prim[s + 1][1], prim[s + 1][2]);
                const double len = A.Distance(B);
                if (!(len > 0.0)) continue;
                const int n = std::max(1, static_cast<int>(std::ceil(len / (0.5 * r))));
                for (int m = 0; m < n; ++m) {
                    const double t = (m + 0.5) / n;
                    const gp_Pnt C(A.X() + t * (B.X() - A.X()), A.Y() + t * (B.Y() - A.Y()), A.Z() + t * (B.Z() - A.Z()));
                    double bestGap = 1e300, bestH = 1e300; gp_Pnt bestQ; bool sampled = false;
                    for (int ci : facesNear(C, r + widest_binding_gap(r, std::max(tRound, tmax)))) {
                        const FaceRec& g = faces[ci];
                        if (!solids[g.solid].rectangular || own[g.solid]) continue;
                        const double tPair = std::max(tRound, solids[g.solid].target);
                        const double dmax = widest_binding_gap(r, tPair);
                        if (boxDistance(g, C) > r + dmax) continue;
                        if (!sampled) { ++out.roundSamples; sampled = true; }
                        BRepExtrema_DistShapeShape dss(BRepBuilderAPI_MakeVertex(C).Vertex(), g.face);
                        if (!dss.IsDone() || dss.NbSolution() < 1) continue;
                        const double gap = dss.Value() - r;
                        if (gap <= touching || gap >= dmax) continue;
                        const double h = chord_bound(r, gap);
                        if (!corridor_binds(h, tRound, solids[g.solid].target) || h >= bestH) continue;
                        bestH = h; bestGap = gap; bestQ = dss.PointOnShape2(1);
                    }
                    if (bestH >= 1e299) continue;
                    ++out.roundHits;
                    const double D = C.Distance(bestQ), f = (r + 0.5 * bestGap) / D;
                    CorridorCandidate c; c.h = bestH; c.gap = bestGap;
                    c.p[0] = C.X() + f * (bestQ.X() - C.X()); c.p[1] = C.Y() + f * (bestQ.Y() - C.Y()); c.p[2] = C.Z() + f * (bestQ.Z() - C.Z());
                    out.points.push_back(c);
                }
            }
    }
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

}  // namespace mvb::mesh
