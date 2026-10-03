#include "mvb/mesh/Mesh3d.h"
// Moved verbatim from OMFEM src/meshing/MasMesher.cpp (ABT #1588, step 2): namespace omfem -> mvb::mesh,
// MasMeshOptions -> MeshOptions, file-local helpers shared with OMFEM's 2D mesher made external. No logic changed.

#include "mvb/mesh/MeshSupport.h"
#include "mvb/mesh/SizeField.h"
#include <cctype>
#include <cerrno>
#include <cstring>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <unordered_map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <dlfcn.h>
#include <sys/stat.h>
#include <gmsh.h>
#include <mmg/mmg3d/libmmg3d.h>
#include "mvb/MagneticBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/SectionBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Symmetry.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/Resistivity.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/StrayCapacitance.h"   // get_wire_insulation_relative_permittivity
#include "physical_models/Temperature.h"          // lumped thermal network (temperature from losses)
#include "processors/MagneticSimulator.h"         // MKF full loss simulation (for the MKF reference)
#include "support/Utils.h"                         // inputs_autocomplete
#include "support/Settings.h"                      // accessories OUT: pin routing off (WP3)
#include "Definitions.h"                          // resolve_dimensional_values
#include "processors/Inputs.h"
#include "support/Utils.h"
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Bnd_Box.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Tool.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <cfloat>

namespace mvb::mesh {

using nlohmann::json;

// CAD EDGES SIDECAR (<out>.cadedges, 2026-09-29). Every mesh edge that lies on a CAD CURVE bounding a
// face between two different regions (or on the model boundary): the facet creases of the swept
// conductors, the mitres where MVB++ primitives join, the core's and the box's edges, the port faces'
// rims. omfem_skinlayer hands exactly these to MMG as ridges, with MMG's angle-based ridge detection
// off: the base mesh's own triangulation of a curved face (20_iso: ~4 deg between 0.4 mm triangles on a
// 5.6 mm turn) and a real junction mitre (1.4-10 deg) overlap in angle, so no angle can tell them
// apart -- only the CAD can. Format: "cadedges 1", the number of edges, then "tagA tagB" per edge
// (node tags of the written .msh).
static void write_cad_edges(const std::string& out) {
    std::map<int, int> volGroup;                                   // volume entity -> physical group tag
    gmsh::vectorpair pgs; gmsh::model::getPhysicalGroups(pgs, 3);
    for (const auto& pg : pgs) {
        std::vector<int> ents; gmsh::model::getEntitiesForPhysicalGroup(3, pg.second, ents);
        for (int e : ents) volGroup[e] = pg.second;
    }
    std::set<int> curves;
    gmsh::vectorpair faces; gmsh::model::getEntities(faces, 2);
    for (const auto& f : faces) {
        std::vector<int> up, down; gmsh::model::getAdjacencies(2, f.second, up, down);
        std::set<int> groups; for (int v : up) { auto it = volGroup.find(v); groups.insert(it == volGroup.end() ? -v : it->second); }
        if (up.size() >= 2 && groups.size() == 1) continue;       // inside one region: not a surface of the model
        for (int c : down) curves.insert(std::abs(c));
    }
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    for (int c : curves) {
        std::vector<int> et; std::vector<std::vector<std::size_t>> tg, nd;
        gmsh::model::mesh::getElements(et, tg, nd, 1, c);
        for (size_t t = 0; t < et.size(); ++t) {
            if (et[t] != 1) throw std::runtime_error("mesh3d: CAD curve " + std::to_string(c) + " carries 1D elements of type " + std::to_string(et[t]) + ", not 2-node lines");
            for (size_t k = 0; k + 1 < nd[t].size(); k += 2) edges.push_back({nd[t][k], nd[t][k + 1]});
        }
    }
    const std::string part = out + ".cadedges.part";
    std::FILE* f = std::fopen(part.c_str(), "w");
    if (!f) throw std::runtime_error("mesh3d: cannot write " + part);
    std::fprintf(f, "cadedges 1\n%zu\n", edges.size());
    for (const auto& e : edges) std::fprintf(f, "%zu %zu\n", e.first, e.second);
    std::fclose(f);
    if (std::rename(part.c_str(), (out + ".cadedges").c_str()) != 0)
        throw std::runtime_error("mesh3d: could not rename " + part + ": " + std::strerror(errno));
    std::fprintf(stderr, "[mesh3d] CAD edges: %zu mesh edges on %zu CAD curves -> %s.cadedges\n", edges.size(), curves.size(), out.c_str());
}

namespace {

// ---- LOCAL CHORD-BOUND REFINEMENT FROM THE BUILT 3D CENTRELINES -----------------------------------
// A tight corridor is where two conductor surfaces come closer than the design-wide conductor size can
// resolve: a facet of size h on a wire of radius r sags h^2/8r into the gap, so the gap d needs
// h <= 0.8 sqrt(4 r d) (the same bound, same factor, as the layout corridors). The layout cross-section
// cannot say WHERE such a corridor is in 3D: MVB++ winds helices, stadiums around rectangular columns and
// irregular EFD paths, not flat circles (measured 2026-09-23: a flat circle of the layout radius overlaps
// the real turn 57 % on ETD49, 15 % on the E55 stadium). So the corridors are found on the real
// centrelines MVB++ built: every sample of every conductor against every segment of every OTHER
// conductor, and of its own conductor further than half a wrap away along the wire.
struct ChordPoints {
    std::vector<double> x, y, z, h, reach;
    double cell = 0.0;                                       // query-grid cell = the largest influence
    std::unordered_map<long long, std::vector<int>> grid;
    size_t samples = 0, close = 0;
    double hmin = 1e30, dmin = 1e30; std::array<double, 3> at{0, 0, 0};
    static long long key(long long i, long long j, long long k) {
        return (i & 0x1FFFFF) | ((j & 0x1FFFFF) << 21) | ((k & 0x1FFFFF) << 42);
    }
    // the size the corridors ask for at a point: h inside the reach, graded to the conductor target over 4 h
    double size_at(double px, double py, double pz, double condT) const {
        double best = 1e30;
        const long long ci = static_cast<long long>(std::floor(px / cell)), cj = static_cast<long long>(std::floor(py / cell)),
                        ck = static_cast<long long>(std::floor(pz / cell));
        for (long long a = -1; a <= 1; ++a) for (long long b = -1; b <= 1; ++b) for (long long c = -1; c <= 1; ++c) {
            const auto it = grid.find(key(ci + a, cj + b, ck + c));
            if (it == grid.end()) continue;
            for (int i : it->second) {
                const double dist = std::sqrt((px - x[i]) * (px - x[i]) + (py - y[i]) * (py - y[i]) + (pz - z[i]) * (pz - z[i]));
                if (dist >= reach[i] + 4.0 * h[i]) continue;
                const double sz = dist <= reach[i] ? h[i] : h[i] + (dist - reach[i]) * (condT - h[i]) / (4.0 * h[i]);
                best = std::min(best, sz);
            }
        }
        return best;
    }
};

template <class Paths>
ChordPoints chord_points_3d(const Paths& paths, double condT) {
    struct Seg { double a[3], b[3]; int c; double s0, s1, r; };
    struct Smp { double p[3]; int c; double s, r; };
    std::vector<Seg> segs; std::vector<Smp> smp;
    double rmin = 1e30, rmax = 0.0, segMax = 0.0;
    int ci = 0;
    for (const auto& pl : paths) {
        const double r = pl.wireRadius;
        if (!(r > 0.0)) throw std::runtime_error("chord bound: conductor '" + pl.name + "' has no positive wire radius");
        rmin = std::min(rmin, r); rmax = std::max(rmax, r);
        double sAcc = 0.0;
        for (const auto& prim : pl.prims)
            for (size_t k = 0; k + 1 < prim.size(); ++k) {
                Seg g; for (int d = 0; d < 3; ++d) { g.a[d] = prim[k][d]; g.b[d] = prim[k + 1][d]; }
                const double len = std::sqrt((g.b[0]-g.a[0])*(g.b[0]-g.a[0]) + (g.b[1]-g.a[1])*(g.b[1]-g.a[1]) + (g.b[2]-g.a[2])*(g.b[2]-g.a[2]));
                if (!(len > 0.0)) continue;
                g.c = ci; g.s0 = sAcc; g.s1 = sAcc + len; g.r = r;
                segs.push_back(g); segMax = std::max(segMax, len);
                const int n = std::max(1, static_cast<int>(std::ceil(len / (0.5 * r))));
                for (int m = 0; m < n; ++m) {
                    const double t = (m + 0.5) / n;
                    Smp q; for (int d = 0; d < 3; ++d) q.p[d] = g.a[d] + t * (g.b[d] - g.a[d]);
                    q.c = ci; q.s = sAcc + t * len; q.r = r; smp.push_back(q);
                }
                sAcc += len;
            }
        ++ci;
    }
    ChordPoints out; out.samples = smp.size();
    if (segs.empty()) throw std::runtime_error("chord bound: the real-winding centrelines have no segments");
    // the widest gap whose bound can still bind the target, for the thinnest wire
    const double dmax = (condT / 0.8) * (condT / 0.8) / (4.0 * rmin);
    const double R = 2.0 * rmax + dmax + 0.5 * segMax;         // candidate segment midpoints within R
    std::unordered_map<long long, std::vector<int>> sg;
    auto cellOf = [&](const double* p, double c, long long& i, long long& j, long long& k) {
        i = static_cast<long long>(std::floor(p[0] / c)); j = static_cast<long long>(std::floor(p[1] / c)); k = static_cast<long long>(std::floor(p[2] / c)); };
    for (size_t i = 0; i < segs.size(); ++i) {
        double m[3]; for (int d = 0; d < 3; ++d) m[d] = 0.5 * (segs[i].a[d] + segs[i].b[d]);
        long long a, b, c; cellOf(m, R, a, b, c); sg[ChordPoints::key(a, b, c)].push_back(static_cast<int>(i));
    }
    // dedupe the corridor points on a grid of half the target: one point (the finest) per cell
    const double dcell = 0.5 * condT;
    std::unordered_map<long long, int> kept;
    for (const auto& q : smp) {
        long long a, b, c; cellOf(q.p, R, a, b, c);
        double bestH = 1e30, bestD = 0.0, mid[3] = {0, 0, 0};
        for (long long u = -1; u <= 1; ++u) for (long long v = -1; v <= 1; ++v) for (long long w = -1; w <= 1; ++w) {
            const auto it = sg.find(ChordPoints::key(a + u, b + v, c + w));
            if (it == sg.end()) continue;
            for (int si : it->second) {
                const Seg& g = segs[si];
                double ab[3], ap[3]; for (int d = 0; d < 3; ++d) { ab[d] = g.b[d] - g.a[d]; ap[d] = q.p[d] - g.a[d]; }
                const double L2 = ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2];
                const double t = std::clamp((ap[0]*ab[0] + ap[1]*ab[1] + ap[2]*ab[2]) / L2, 0.0, 1.0);
                if (g.c == q.c && std::fabs(g.s0 + t * (g.s1 - g.s0) - q.s) < M_PI * (q.r + g.r + dmax)) continue;  // own wire, nearby
                double cq[3], dv[3]; for (int d = 0; d < 3; ++d) { cq[d] = g.a[d] + t * ab[d]; dv[d] = cq[d] - q.p[d]; }
                const double c2c = std::sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
                const double gap = c2c - q.r - g.r;
                if (gap <= 1e-7 || !(c2c > 0.0)) continue;   // touching / overlapping: the layout check owns those
                const double h = 0.8 * std::sqrt(4.0 * std::min(q.r, g.r) * gap);
                if (h >= condT || h >= bestH) continue;
                bestH = h; bestD = gap;
                const double f = (q.r + 0.5 * gap) / c2c;
                for (int d = 0; d < 3; ++d) mid[d] = q.p[d] + f * dv[d];
            }
        }
        if (bestH >= 1e29) continue;
        ++out.close;
        if (bestH < out.hmin) { out.hmin = bestH; out.dmin = bestD; out.at = {mid[0], mid[1], mid[2]}; }
        long long a2, b2, c2; cellOf(mid, dcell, a2, b2, c2);
        const long long kk = ChordPoints::key(a2, b2, c2);
        const auto it = kept.find(kk);
        if (it != kept.end() && out.h[it->second] <= bestH) continue;
        if (it == kept.end()) {
            kept[kk] = static_cast<int>(out.h.size());
            out.x.push_back(mid[0]); out.y.push_back(mid[1]); out.z.push_back(mid[2]); out.h.push_back(bestH); out.reach.push_back(1.5 * bestH);
        } else {
            const int i = it->second;
            out.x[i] = mid[0]; out.y[i] = mid[1]; out.z[i] = mid[2]; out.h[i] = bestH; out.reach[i] = 1.5 * bestH;
        }
    }
    double infl = 0.0; for (size_t i = 0; i < out.h.size(); ++i) infl = std::max(infl, out.reach[i] + 4.0 * out.h[i]);
    out.cell = infl > 0.0 ? infl : condT;
    for (size_t i = 0; i < out.h.size(); ++i) {
        double p[3] = {out.x[i], out.y[i], out.z[i]}; long long a, b, c; cellOf(p, out.cell, a, b, c);
        out.grid[ChordPoints::key(a, b, c)].push_back(static_cast<int>(i));
    }
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// 3D tetrahedral mesh from a MAS magnetic. Same MVB++ -> STEP -> gmsh idiom as
// mesh_from_mas, but keeps the FULL 3D solids (real turns; no cut2DFaces) and
// tets them. Volume groups: core / winding / air. Boundary surface: outer.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------
// MAPPED COPPER (OMFEM_SKIN_MAPPED, ABT #1156, Alf 2026-09-11). The skin-effect mesh of a
// rectangular conductor is a MAPPED hex grid graded toward all four faces (Flux: "mapped or
// extruded structured elements"; COMSOL: boundary layers, first cell delta/2..delta, coarse
// along the perimeter). MVB++ sweeps the wire as one solid per parallel whose lateral faces are
// split per centreline primitive, so the solid is cut at its OWN junction rings (the four
// straight section edges between consecutive primitives) into 6-face blocks, each block meshed
// transfinite: a calibrated two-sided ("Bump") distribution across width and thickness, uniform
// along the wire. Validated on the 3 x 0.5 mm bar against 1.72M-tet references: 1 MHz -1.2 %
// (3780 hexes, 16 s), 100 kHz -1.3 %. The air/core stay tets; gmsh puts pyramids on the quads.
// ---------------------------------------------------------------------------------------------
namespace {
struct MappedCurve { double len = 0.0; bool straight = false; };
MappedCurve mapped_curve_info(int c) {
    MappedCurve mc; gmsh::model::occ::getMass(1, c, mc.len);
    gmsh::vectorpair pts; gmsh::model::getBoundary({{1, c}}, pts, false, false, false);
    if (pts.size() == 2 && mc.len > 0.0) {
        std::vector<double> a, b;
        gmsh::model::getValue(0, std::abs(pts[0].second), {}, a);
        gmsh::model::getValue(0, std::abs(pts[1].second), {}, b);
        const double chord = std::sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));
        mc.straight = chord / mc.len > 0.9995;      // STEP export turns lines into B-splines: judge by chord/length
    }
    return mc;
}
bool mapped_is_width(const MappedCurve& mc, double W)  { return mc.straight && std::fabs(mc.len - W) < 0.02 * W; }
bool mapped_is_thick(const MappedCurve& mc, double T)  { return mc.straight && std::fabs(mc.len - T) < 0.05 * T; }

// Split one swept conductor solid at its junction rings. Returns the block tags (the input solid
// is consumed by the fragment). Throws when a block is not a 6-face hexahedral topology.
std::vector<int> mapped_split_at_rings(int solid, double W, double T, const std::string& region) {
    gmsh::vectorpair fdt; gmsh::model::getBoundary({{3, solid}}, fdt, false, false, false);
    std::vector<int> faces; for (auto& f : fdt) faces.push_back(std::abs(f.second));
    std::map<int, std::set<int>> faceCurves; std::set<int> sect, seenCurves;
    std::map<int, std::vector<int>> ends; std::map<int, std::set<int>> adj;
    for (int f : faces) {
        gmsh::vectorpair cdt; gmsh::model::getBoundary({{2, f}}, cdt, false, false, false);
        for (auto& c : cdt) {
            const int ct = std::abs(c.second); faceCurves[f].insert(ct);
            if (!seenCurves.insert(ct).second) continue;
            const MappedCurve mc = mapped_curve_info(ct);
            if (!(mapped_is_width(mc, W) || mapped_is_thick(mc, T))) continue;
            sect.insert(ct);
            gmsh::vectorpair pdt; gmsh::model::getBoundary({{1, ct}}, pdt, false, false, false);
            for (auto& p : pdt) { ends[ct].push_back(std::abs(p.second)); adj[std::abs(p.second)].insert(ct); }
        }
    }
    // rings = connected components of the section curves over shared end points
    std::set<int> seen; std::vector<std::vector<int>> rings;
    for (int c : sect) {
        if (seen.count(c)) continue;
        std::set<int> comp; std::vector<int> stack{c};
        while (!stack.empty()) { int x = stack.back(); stack.pop_back(); if (comp.count(x)) continue; comp.insert(x);
            for (int p : ends[x]) for (int y : adj[p]) if (!comp.count(y)) stack.push_back(y); }
        for (int x : comp) seen.insert(x);
        rings.emplace_back(comp.begin(), comp.end());
    }
    gmsh::vectorpair tools; int nIncomplete = 0, nCap = 0;
    for (auto& r : rings) {
        if (r.size() != 4) { ++nIncomplete; continue; }
        bool cap = false; const std::set<int> rs(r.begin(), r.end());
        for (int f : faces) if (faceCurves[f] == rs) { cap = true; break; }
        if (cap) { ++nCap; continue; }
        const int loop = gmsh::model::occ::addCurveLoop(r);
        tools.push_back({2, gmsh::model::occ::addPlaneSurface({loop})});
    }
    gmsh::model::occ::synchronize();
    gmsh::vectorpair out; std::vector<gmsh::vectorpair> outMap;
    gmsh::model::occ::fragment({{3, solid}}, tools, out, outMap, -1, true, true);
    gmsh::model::occ::synchronize();
    std::vector<int> blocks; std::map<int,int> facesPer;
    for (auto& o : out) if (o.first == 3) {
        blocks.push_back(o.second);
        gmsh::vectorpair bf; gmsh::model::getBoundary({{3, o.second}}, bf, false, false, false);
        facesPer[(int)bf.size()]++;
    }
    std::fprintf(stderr, "[mesh3d] mapped copper: %s solid %d: %zu section rings (%d caps, %d incomplete) -> %zu cut(s) -> %zu block(s)",
                 region.c_str(), solid, rings.size(), nCap, nIncomplete, tools.size(), blocks.size());
    for (auto& kv : facesPer) std::fprintf(stderr, " [%d faces x%d]", kv.first, kv.second);
    std::fprintf(stderr, "\n");
    if (facesPer.size() != 1 || facesPer.begin()->first != 6) {
        std::ostringstream m; m << "mesh3d_from_mas: mapped copper: " << region << " solid " << solid
          << " does not split into 6-face blocks (rings " << rings.size() << ", incomplete " << nIncomplete
          << "; faces per block:"; for (auto& kv : facesPer) m << " " << kv.first << "x" << kv.second;
        m << "). A junction whose section edges are not W/T straight lines (mitre?) or a face imprinted by a neighbour breaks the block topology.";
        throw std::runtime_error(m.str());
    }
    return blocks;
}

// gmsh's two-sided "Bump" distribution has no closed form for its first cell: calibrate (n, coef)
// on a 1D line so the first cell is <= hs and the wall growth <= growth, with the fewest cells.
struct BumpFit { int n = 0; double coef = 0.0, first = 0.0, grow = 0.0, maxCell = 0.0; };
BumpFit mapped_calibrate_bump(double length, double hs, double growth) {
    std::string cur; gmsh::model::getCurrent(cur);
    auto probe = [&](int n, double coef, double& first, double& grow, double& maxc) {
        gmsh::model::add("__omfem_bump_cal");
        const int p0 = gmsh::model::geo::addPoint(0, 0, 0), p1 = gmsh::model::geo::addPoint(length, 0, 0);
        const int c = gmsh::model::geo::addLine(p0, p1); gmsh::model::geo::synchronize();
        gmsh::model::mesh::setTransfiniteCurve(c, n, "Bump", coef);
        gmsh::model::mesh::generate(1);
        std::vector<std::size_t> nt; std::vector<double> xyz, par;
        gmsh::model::mesh::getNodes(nt, xyz, par, 1, c, true, false);
        std::vector<double> xs; for (size_t i = 0; i < xyz.size(); i += 3) xs.push_back(xyz[i]);
        std::sort(xs.begin(), xs.end());
        first = xs[1] - xs[0]; grow = (xs[2] - xs[1]) / first; maxc = 0.0;
        for (size_t i = 1; i < xs.size(); ++i) maxc = std::max(maxc, xs[i] - xs[i-1]);
        gmsh::model::remove(); gmsh::model::setCurrent(cur);
    };
    BumpFit best;
    for (int k = 1; k <= 49; ++k) {
        const double coef = 0.02 * k;
        int lo = 4, hi = 400; double first, grow, maxc;
        probe(hi, coef, first, grow, maxc); if (first > 1.05 * hs) continue;
        while (hi - lo > 1) { const int mid = (lo + hi) / 2; probe(mid, coef, first, grow, maxc); (first <= 1.05 * hs ? hi : lo) = mid; }
        probe(hi, coef, first, grow, maxc);
        if (grow > 1.05 * growth) continue;
        if (best.n == 0 || hi < best.n || (hi == best.n && std::fabs(first - hs) < std::fabs(best.first - hs)))
            best = {hi, coef, first, grow, maxc};
    }
    if (best.n == 0) {
        std::ostringstream m; m << "mesh3d_from_mas: mapped copper: no Bump distribution on a " << length*1e3
          << " mm edge reaches a first cell of " << hs*1e3 << " mm with growth <= " << growth;
        throw std::runtime_error(m.str());
    }
    return best;
}
} // namespace

std::string mesh3d_from_mas(json magnetic, const MeshOptions& opt) {

    auto step = [](const char* s){ std::fprintf(stderr, "[mesh3d] %s\n", s); };
    // REAL-WINDING geometry is the DEFAULT: ONE continuous copper body per (winding, parallel)
    // (MVB++ useRealWindingGeometry+femReady) with the actual helical pitch and lead-in/out, and NO
    // symmetry (the leads break every mirror plane). The magnetostatic solver's stranded J=j0*e_phi
    // coil model is geometry-agnostic, so it drives the continuous winding region directly; drive
    // full amp-turns with OMFEM_SYMMULT=1.
    //
    // OMFEM_CLOSED_TURNS restores the old per-turn CLOSED-LOOP model -- independent idealised
    // rings, no pitch, no leads -- which is TO BE DEPRECATED for FEM use: it is a geometry
    // approximation (rings carry no net axial transport, terminals must be faked from symmetry cut
    // planes, and the SERIES bordered machinery exists solely to compensate for it). It remains
    // available for now because it admits mirror symmetry (quarter/eighth meshes the continuous
    // winding cannot) and the thermal path can place per-turn losses on its turn_<w>_<i> regions.
    // Do not build new work on it. (The per-turn compound remains fully supported for DRAWING in
    // MVB++ -- this deprecation is about FEM meshing only.) OMFEM_CONTINUOUS is accepted for
    // backward compatibility and overrides OMFEM_CLOSED_TURNS.
    const bool continuous = std::getenv("OMFEM_CONTINUOUS") != nullptr ||
                            std::getenv("OMFEM_CLOSED_TURNS") == nullptr;
    if (!continuous)
        std::fprintf(stderr,
            "[mesh3d] WARNING: OMFEM_CLOSED_TURNS selects the per-turn closed-loop winding model "
            "(idealised rings, no pitch, no leads). This model is TO BE DEPRECATED for FEM use; "
            "the real/continuous winding is the default and the supported path.\n");
    else if (opt.symmetry != "none" && !opt.symmetry.empty())
        std::fprintf(stderr,
            "[mesh3d] note: symmetry '%s' ignored -- the real/continuous winding's helical "
            "lead-in/out breaks every mirror plane (use OMFEM_CLOSED_TURNS for the deprecated "
            "symmetric per-turn model)\n", opt.symmetry.c_str());
    // ACCESSORIES OUT (Alf, 2026-09-21): "make sure that all work for sleeves, bobbins, pins, etc is
    // flagged out for OMFEM for now". OMFEM simulates core + copper; the manufacturing accessories
    // MVB++ has started drawing (WP0-WP7: bobbin, pins, dividers, toroid bases, spacers, shunts,
    // sleeves) stay out until they are validated as FEM inputs. Two halves:
    //   (a) WP3 routes terminal leads to the bobbin's pins when MKF's coil_connect_leads_to_pins is
    //       on. MKF defaults it OFF, but a default is not a decision: set it here, so the lead
    //       geometry OMFEM meshes cannot change because some other caller flipped a global.
    //   (b) the accessory solids themselves are dropped by ROLE after buildAllNamed, below.
    // OMFEM_ACCESSORIES=1 restores both. Every dropped solid is named on stderr -- never silent.
    const bool accessories = std::getenv("OMFEM_ACCESSORIES") != nullptr;
    if (!accessories) {
        auto& mkfs = OpenMagnetics::Settings::GetInstance();
        mkfs.set_coil_connect_leads_to_pins(false);
        mkfs.set_coil_quick_bobbin_generate_pins(false);
        std::fprintf(stderr, "[mesh3d] accessories OUT: leads not routed to bobbin pins, no pins generated "
                             "(OMFEM_ACCESSORIES=1 restores)\n");
    }
    step(continuous ? "enrich (real winding)" : "enrich");
    OpenMagnetics::Magnetic enriched =
        mvb::magnetic_autocomplete_safe(magnetic, /*useRealWindingGeometry=*/continuous);
    // TERMINALS ON ONE PLANE (ABT #1155): the FEM port is a planar box face, so every lead tip
    // of a winding's parallels must lie on it -- MVB++ pins the concentric fan and (2026-09-11)
    // the toroidal lead tips to that plane under this switch. Same mode as the Ansys batches.
    setenv("MVB_FAN_TERMINALS_ON_PLANE", "1", 0);
    // CENTRELINE SIDECAR (<out>.path, metres): the conductors' real centrelines, one polyline per
    // primitive, for omfem_skinlayer's wire tangent (round wire has no perpendicular faces to
    // derive it from). Written for the continuous (real-winding) build only.
    // the same centrelines place the local chord-bound refinement (chord_points_3d), so they are kept
    using RealPaths = decltype(std::declval<mvb::MagneticBuilder&>().buildRealWindingPaths(enriched));
    std::optional<RealPaths> realPaths;
    std::string realPathsError = "not a real-winding build";
    if (continuous) {
        try {
            mvb::MagneticBuilder pb; auto paths = pb.buildRealWindingPaths(enriched);
            realPaths = paths;
            std::FILE* pf = std::fopen((opt.out_msh + ".path").c_str(), "w");
            if (pf) {
                std::fprintf(pf, "%zu\n", paths.size());
                for (const auto& pl : paths) {
                    std::fprintf(pf, "%s|%zu|%.6g\n", pl.name.c_str(), pl.prims.size(), pl.wireRadius);
                    for (const auto& seg : pl.prims) { std::fprintf(pf, "%zu", seg.size()); for (const auto& q : seg) std::fprintf(pf, " %.9g %.9g %.9g", q[0], q[1], q[2]); std::fprintf(pf, "\n"); }
                }
                std::fclose(pf);
                std::fprintf(stderr, "[mesh3d] centreline sidecar: %s (%zu conductor(s))\n", (opt.out_msh + ".path").c_str(), paths.size());
            }
        } catch (const std::exception& ex) { realPathsError = ex.what(); std::fprintf(stderr, "[mesh3d] centreline sidecar skipped: %s\n", ex.what()); }
    }
    step("buildAllNamed (3D solids)");
    mvb::MagneticBuilder builder;
    // Turn-part fusing (MVB++ global): OFF by default here to preserve the existing solids-path
    // behaviour (thermal FEA meshes the unfused compound fine). OMFEM_FUSE_TURNS=1 fuses each turn
    // into ONE solid -- needed for the winding-loss eddy mesh (overlapping compound sub-solids crash
    // tet meshing), at the cost of a boolean per turn. NB: fusing alone is not yet sufficient for a
    // clean eddy mesh -- the turn sub-solids overlap by design, so the fuse leaves artifact faces;
    // the real fix is a single swept-solid turn with proper corners (see 3D_SERIES_WINDING doc).
    mvb::TurnBuilder::setFuseTurnParts(std::getenv("OMFEM_FUSE_TURNS") != nullptr);
    // Real winding has no usable symmetry plane -> force full geometry.
    const int symPlanes = continuous ? 0 : mvb::parse_symmetry_spec(opt.symmetry);
    // Include the bobbin (default) so the 3D thermal FEA has a winding->bobbin->core conduction
    // bridge -- without it the winding is thermally isolated (sheds all loss by convection alone,
    // over-hot, and the core stays at ambient). Opt out with OMFEM_NO_BOBBIN.
    const bool include_bobbin = !std::getenv("OMFEM_NO_BOBBIN");
    // OMFEM_COATING: emit each turn as bare-copper core + outer insulated shell ("<turn> coating"),
    // so the thermal FEA resolves the low-k wire enamel between turns (the mesh fragment splits the
    // overlap into copper + coating annulus). Forces copper turns (paint_coating ignored for turns).
    const bool emit_coating = std::getenv("OMFEM_COATING") != nullptr;
    // OMFEM_INSULATION: also build inter-layer/inter-section insulation solids ("insulation_layer_*")
    // when they carry real thickness, tagged as a low-k conduction barrier between windings.
    const bool incl_insul = std::getenv("OMFEM_INSULATION") != nullptr;
    // Core coating (epoxy/parylene) as a conformal shell solid around the core. The thickness MUST
    // come from the MAS core's EXPLICIT coating: MKF shrinks the winding window by exactly that
    // thickness (Core.cpp, TOROIDAL/CLOSED_SHAPE) only for an explicit catalogue coating, so the
    // turns are placed on the COATED surface and the shell fills the reserved gap -- no collision.
    // Forcing a thickness on an uncoated core (where MKF reserved no space) makes the shell collide
    // with the turns, so OMFEM_CORE_COATING is a manual override for testing geometry only.
    double core_coat_t = 0.0;
    if (enriched.get_core().get_functional_description().get_coating()) {
        try { core_coat_t = enriched.get_core().get_coating_thickness(); } catch (...) { core_coat_t = 0.0; }
    }
    if (const char* e = std::getenv("OMFEM_CORE_COATING")) core_coat_t = std::atof(e);  // testing override
    // Enriched json for adaptive gap sizing (concentric core axis is Y). Cheap (no re-enrich).
    json enr3; OpenMagnetics::to_json(enr3, enriched);
    const std::vector<CoreGap> core_gaps3 = extract_core_gaps(enr3);
    const double post_half3 = central_column_half_width(enr3);

    // AUTO CONDUCTOR TARGET (used when opt.conductor_target == 0; < 0 disables conductor
    // refinement entirely). The viable surface-mesh size on a winding is bounded by CLEARANCE,
    // not wire size: two nearly-parallel curved surfaces with clearance g, meshed with chord
    // sagitta h^2/(8 r) each, CROSS when 2*h^2/(8 r) > g -- gmsh then dies with "PLC Error: a
    // segment and a facet intersect" or "overlapping facets". The bound h < sqrt(4 r g) is
    // computable EXACTLY from the enriched turnsDescription (true pairwise turn clearances) and
    // the wire dimensions. This replaces per-design hand-tuned targets, which caused every
    // earlier sweep to measure its own guess instead of the geometry (ABT #332).
    // LAZY -- computed only when the meshing path actually seeds a conductor target, so the
    // ANISO and STEP_ONLY paths (which never consume it) cannot be aborted by it.
    // ---- LOCAL CHORD-BOUND REFINEMENT (2026-09-20) ------------------------------------------
    // A CORRIDOR is the narrow air channel between two nearly-touching EMITTED surfaces. Its
    // width d bounds the chord sagitta of the facets that face into it -- and it bounds it THERE,
    // nowhere else. The previous rule reduced every corridor in the layout to a single minimum
    // and applied that size to every conductor in the design: on 01_simple_inductor_etd34_n87
    // that meant 2,314,933 tets instead of 383,187 (6.0x) and 4253 s instead of 1770 s, for a
    // 0.2% change in the loss ratio; on 17_cllc_xfmr_e5528_3c92a two sub-micron pairs out of 780
    // drove the whole winding to 0.06 mm and the mesher to 36 GB. The per-pair geometry needed to
    // place the refinement was already being computed in the double loop below and thrown away.
    struct TurnCorridor {
        double x, y;        // narrowest point of the channel, in LAYOUT cross-section coords [m]
        double clearance;   // measured surface-to-surface clearance there [m]
        double size;        // element size the chord bound demands there [m]
        double reach;       // how far along/around the channel that size must carry [m]
        std::string turnA, turnB;   // the two turns that bound it (named in every refusal)
    };
    struct LayoutSizing {
        double rsmin = std::numeric_limits<double>::max();  // smallest emitted surface radius [m]
        double gmin  = std::numeric_limits<double>::max();  // smallest measured pair clearance [m]
        double coat_corridor = std::numeric_limits<double>::max();  // ro - rc [m]
        int touching = 0, overlapping = 0, touchingAcrossLayers = 0;   // touching: across layers per MKF + same-layer at d <= 1e-7
        double worst_overlap = 0.0;
        std::vector<TurnCorridor> corridors;      // measurable corridors, tightest first
        std::vector<std::array<double,2>> turn_xy;  // every turn centre, layout coords [m]
    };
    auto compute_layout_sizing = [&]() -> LayoutSizing {
        // Radii per winding: rc = copper, rs = the EMITTED surface (the coating shell when
        // coating solids are emitted, bare copper otherwise) -- the clearance that bounds the
        // mesh is between the surfaces gmsh actually triangulates.
        struct WDim { double rc, ro, rs; };
        std::map<std::string, WDim> wdims;
        OpenMagnetics::Coil coilW = enriched.get_coil();  // mutable copy: resolve_wire is non-const
        const auto& fdW = coilW.get_functional_description();
        // MKF's resolver (nominal -> (min+max)/2 -> max -> min; throws if none), not a
        // hand-read .nominal -- same rule as every other dimension in this file.
        auto res = [](const auto& optDim) -> double {
            return OpenMagnetics::resolve_dimensional_values(optDim.value());
        };
        for (size_t i = 0; i < fdW.size(); ++i) {
            const std::string wname = fdW[i].get_name();
            auto wire = coilW.resolve_wire(i);   // handles catalogue-name wires ("Round 0.4 ...")
            const auto cd = wire.get_conducting_diameter();
            const auto od = wire.get_outer_diameter();
            const auto cw = wire.get_conducting_width();
            const auto ch = wire.get_conducting_height();
            double rc = 0.0, ro = 0.0;
            if (cd) {
                rc = 0.5 * res(cd);
                ro = od ? 0.5 * res(od) : rc;
                // the layout's copper radius must be the one MVB++ builds (see the litz branch)
                if (realPaths)
                    for (const auto& pl : *realPaths)
                        if (pl.name.rfind(wname + " parallel ", 0) == 0 && std::fabs(pl.wireRadius - rc) > 1e-9)
                            throw std::runtime_error("mesh3d auto-sizing: winding '" + wname + "' is built at radius " +
                                std::to_string(pl.wireRadius * 1e3) + " mm but its conducting radius is " +
                                std::to_string(rc * 1e3) + " mm");
            } else if (cw && ch) {
                rc = 0.5 * std::min(res(cw), res(ch));
                ro = rc;
                const auto ow = wire.get_outer_width();
                const auto oh = wire.get_outer_height();
                if (ow && oh) ro = 0.5 * std::min(res(ow), res(oh));
            } else if (od) {
                // Litz: the emitted conductor is ONE bundle cylinder. Its radius is what MVB++ BUILDS,
                // not the MAS outer diameter: 17_cllc's litz (outer 1.880..2.049 mm -> 0.98225 mm) is
                // emitted at 0.96598 mm (Primary solid 5743.2 mm3 over 2051.5 mm of centreline = a
                // 12-gon of that circumradius), so sizing on the outer radius made every litz gap
                // 33 um too small and invented "tight corridors" the model does not have. A real-winding
                // build reads the emitted radius from the centrelines; otherwise the outer radius stands.
                ro = 0.5 * res(od);
                rc = ro;
                if (realPaths) {
                    std::optional<double> pr;
                    for (const auto& pl : *realPaths)
                        if (pl.name.rfind(wname + " parallel ", 0) == 0) pr = pl.wireRadius;
                    if (!pr) throw std::runtime_error("mesh3d auto-sizing: no built centreline for litz winding '" + wname + "'");
                    rc = *pr;
                }
            } else {
                throw std::runtime_error("mesh3d auto-sizing: winding '" + wname + "' wire has "
                                         "neither conductingDiameter nor conductingWidth/Height "
                                         "nor outerDiameter -- set OMFEM_COND_TARGET explicitly");
            }
            const double rs = emit_coating ? ro : rc;
            if (!(rs > 0.0))
                throw std::runtime_error("mesh3d auto-sizing: winding '" + wname +
                                         "' resolved to a non-positive wire radius");
            wdims[wname] = {rc, ro, rs};
        }
        // When only the COPPER is emitted, the tightest real air corridor around it is one
        // coating thickness: the wire's outer surface is what rests on the bobbin and on
        // neighbouring turns, so the bare-copper surface sits ro - rc away from those contacts.
        // This is RESOLVED wire data, not a fabricated floor -- a wire with no declared coating
        // (ro == rc) contributes nothing. Measured on e138: the binding clearance is exactly
        // this 19.5 um coating corridor, not the 0.714 mm turn-centre spacing.
        double gmin0 = std::numeric_limits<double>::max();
        if (!emit_coating)
            for (const auto& [wn, wd] : wdims) {
                (void)wn;
                if (wd.ro > wd.rc + 1e-9) gmin0 = std::min(gmin0, wd.ro - wd.rc);
            }
        const auto& tdW = coilW.get_turns_description();
        if (!tdW || tdW->empty())
            throw std::runtime_error("mesh3d auto-sizing: enriched MAS has no turnsDescription -- "
                                     "set OMFEM_COND_TARGET explicitly");
        struct TP { double x, y, rs; std::string name, layer; };
        std::vector<TP> tps;
        for (const auto& t : *tdW) {
            const auto it = wdims.find(t.get_winding());
            if (it == wdims.end()) continue;
            const auto& c = t.get_coordinates();
            tps.push_back({c.at(0), c.at(1), it->second.rs, t.get_name(), t.get_layer() ? *t.get_layer() : std::string()});
        }
        // TOUCHING ACROSS LAYERS, FROM MKF's LAYER DATA (Alf 2026-09-25, ABT #1402). Two turns in DIFFERENT layers with
        // nothing between those layers -- no layer at all, or only zero-thickness (virtual) ones -- are placed by MKF to
        // touch at their radius sum. Where the layers' pitches differ (flyback_complete's primary: 0.654 vs 1.154 mm)
        // their helices cross, so the 2D section finds them anywhere from touching to a fraction of a micron apart
        // (0.164 um there), and the old fixed d <= 1e-7 cut turned that into a corridor demanding 11.6 um elements. MKF
        // always respects this, so the pair is TOUCHING, whatever the sampled d; nothing between them is meshed. A layer
        // of non-zero thickness between them (insulation, or another conduction layer) keeps the pair a real corridor.
        // Same-layer pairs keep the gap rule below.
        struct LayerBox { double c[2]; double w[2]; int orient; };
        std::map<std::string, LayerBox> layerBox;
        if (const auto& lds = coilW.get_layers_description()) {
            for (const auto& l : *lds) {
                const auto& lc = l.get_coordinates();
                const auto& ld = l.get_dimensions();
                if (lc.size() < 2 || ld.size() < 2)
                    throw std::runtime_error("mesh3d auto-sizing: layer '" + l.get_name() + "' has no 2D coordinates/dimensions");
                LayerBox b;
                b.c[0] = lc[0]; b.c[1] = lc[1]; b.w[0] = ld[0]; b.w[1] = ld[1];
                b.orient = static_cast<int>(l.get_orientation());
                layerBox[l.get_name()] = b;
            }
        }
        // true when layers A and B differ, stack along one axis, and every layer strictly between them along it is
        // zero-thickness
        auto layers_touch = [&](const std::string& A, const std::string& B) -> bool {
            if (A.empty() || B.empty() || A == B) return false;
            const auto ia = layerBox.find(A), ib = layerBox.find(B);
            if (ia == layerBox.end() || ib == layerBox.end())
                throw std::runtime_error("mesh3d auto-sizing: turn layer '" + (ia == layerBox.end() ? A : B) +
                                         "' is not in the coil's layersDescription");
            const LayerBox &a = ia->second, &b = ib->second;
            if (a.orient != b.orient) return false;
            const int ax = (a.orient == static_cast<int>(MAS::WindingOrientation::CONTIGUOUS)) ? 1 : 0;   // overlapping layers stack along axis 0
            const double lo = std::min(a.c[ax], b.c[ax]), hi = std::max(a.c[ax], b.c[ax]);
            const int ot = 1 - ax;
            for (const auto& [n, k] : layerBox) {
                if (n == A || n == B) continue;
                if (!(k.c[ax] > lo && k.c[ax] < hi)) continue;
                // it must lie across the two layers, not beside them
                const double kLo = k.c[ot] - 0.5 * k.w[ot], kHi = k.c[ot] + 0.5 * k.w[ot];
                const double aLo = std::max(a.c[ot] - 0.5 * a.w[ot], b.c[ot] - 0.5 * b.w[ot]);
                const double aHi = std::min(a.c[ot] + 0.5 * a.w[ot], b.c[ot] + 0.5 * b.w[ot]);
                if (kHi <= aLo || kLo >= aHi) continue;
                if (k.w[ax] > 0.0) return false;
            }
            return true;
        };
        LayoutSizing ls;
        ls.coat_corridor = gmin0;
        for (const auto& [wn, wd] : wdims) { (void)wn; ls.rsmin = std::min(ls.rsmin, wd.rs); }
        // A pair only ever matters if its bound could bind SOME target we would plausibly ask
        // for; 1.2 * rsmin (wire-scale elements) is the coarsest conductor target this mesher
        // ever uses, so a pair whose bound is above it constrains nothing and is not stored.
        // That keeps the corridor list a NEAR-NEIGHBOUR list (O(turns)), not the O(turns^2) pair
        // list: a bound below 1.2 rs means a clearance below 0.56 rs.
        const double keep_below = 1.2 * ls.rsmin;
        for (size_t i = 0; i < tps.size(); ++i)
            for (size_t j = i + 1; j < tps.size(); ++j) {
                const double dx = tps[j].x - tps[i].x, dy = tps[j].y - tps[i].y;
                const double c2c = std::hypot(dx, dy);
                const double d = c2c - tps[i].rs - tps[j].rs;
                // only a pair tight enough to become a corridor (see keep_below below) needs the layer verdict
                if (d > -1e-9 && 0.8 * std::sqrt(4.0 * std::min(tps[i].rs, tps[j].rs) * std::max(d, 0.0)) < keep_below &&
                    layers_touch(tps[i].layer, tps[j].layer)) { ++ls.touching; ++ls.touchingAcrossLayers; continue; }
                if (d <= 1e-7) {
                    if (d > -1e-9) ++ls.touching;
                    else { ++ls.overlapping; ls.worst_overlap = std::max(ls.worst_overlap, -d); }
                    continue;
                }
                ls.gmin = std::min(ls.gmin, d);
                // The chord bound takes the SMALLER radius: the tighter curvature is the one
                // whose facets cut into the corridor first.
                const double reff = std::min(tps[i].rs, tps[j].rs);
                const double h = 0.8 * std::sqrt(4.0 * reff * d);
                if (h >= keep_below || !(c2c > 0.0)) continue;
                TurnCorridor tc;
                const double f = (tps[i].rs + 0.5 * d) / c2c;   // the channel's narrowest point
                tc.x = tps[i].x + f * dx;
                tc.y = tps[i].y + f * dy;
                tc.clearance = d;
                tc.size = h;
                // Reach. Offset s along the channel from its narrowest point, the two surfaces
                // separate as d + s^2/reff (each curving away by s^2/2reff), so at s = 1.5 h =
                // 2.4 sqrt(reff d) they are ~6.8 d apart and the LOCAL bound there is already
                // 2.6x looser. 1.5 h is therefore where the fine size stops being needed.
                tc.reach = 1.5 * h;
                tc.turnA = tps[i].name; tc.turnB = tps[j].name;
                ls.corridors.push_back(tc);
            }
        // Interpenetrating emitted surfaces are a GEOMETRY defect (they will break the OCC
        // fragment or the boundary recovery downstream) -- surface it here, where the diagnosis
        // lands on the layout, instead of letting the mesher blame the element size.
        if (ls.overlapping > 0)
            std::fprintf(stderr, "[mesh3d] WARNING: %d turn pair(s) INTERPENETRATE (worst "
                         "%.4g mm) in the enriched layout -- meshing will likely fail on "
                         "geometry, not sizing\n", ls.overlapping, ls.worst_overlap * 1e3);
        std::sort(ls.corridors.begin(), ls.corridors.end(),
                  [](const TurnCorridor& a, const TurnCorridor& b) { return a.size < b.size; });
        // Merge corridors that are the same channel seen from two pairs: a later (coarser) one
        // whose narrowest point already lies inside a kept one's tube adds nothing.
        std::vector<TurnCorridor> merged;
        for (const auto& c : ls.corridors) {
            bool covered = false;
            for (const auto& k : merged)
                if (std::hypot(c.x - k.x, c.y - k.y) <= k.reach && c.size >= k.size)
                    { covered = true; break; }
            if (!covered) merged.push_back(c);
        }
        ls.corridors.swap(merged);
        for (const auto& t : tps) ls.turn_xy.push_back({t.x, t.y});
        return ls;
    };
    // Memoised: the layout does not change, and both the pre-CAD prediction and the size field
    // want the same measurement.
    std::optional<LayoutSizing> layout_sizing_cache;
    auto layout_sizing = [&]() -> const LayoutSizing& {
        if (!layout_sizing_cache) layout_sizing_cache = compute_layout_sizing();
        return *layout_sizing_cache;
    };
    // The GLOBAL conductor target is now a PHYSICS choice (resolve the wire), not a corridor
    // one: the corridors are refined locally by add_ring_bands below. `contact_floor` is the one
    // case where a corridor cannot be placed from the layout: when a bobbin / insulation solid
    // is emitted, the wire rests ON it and its bare-copper surface sits one coating thickness
    // (ro - rc) from that solid, everywhere the winding touches it. That contact is not a turn
    // PAIR, so it has no place in the layout cross-section this code can locate, and it stays a
    // global bound. When nothing is emitted there (no bobbin, no insulation, no coating shell)
    // the ro - rc corridor exists in NO surface gmsh will triangulate -- it was a phantom, and
    // it is what drove 7 of the 10 corpus designs (01, 05, 12, 15, 18, buck, CMC: every one of
    // them reported "min clearance" exactly equal to ro - rc with ZERO touching pairs, while
    // their measured copper-to-copper clearance was 2x-70x larger).
    auto compute_auto_cond_target = [&](bool contact_floor) -> double {
        const LayoutSizing& ls = layout_sizing();
        double t = 1.2 * ls.rsmin;                    // wire-scale elements
        if (contact_floor && ls.coat_corridor < 1e30)
            t = std::min(t, 0.8 * std::sqrt(4.0 * ls.rsmin * ls.coat_corridor));
        if (!(t >= 2e-5))
            throw std::runtime_error(
                "mesh3d auto-sizing: the emitted contact clearance (" +
                std::to_string(ls.coat_corridor * 1e6) + " um) demands elements under 20 um -- "
                "an infeasible mesh. Fix the layout, or set OMFEM_COND_TARGET explicitly.");
        std::fprintf(stderr, "[mesh3d] auto conductor target: %.4g mm  (min surface radius %.4g "
                     "mm, min measured pair clearance %.4g mm, %zu local corridor(s), touching "
                     "pairs %d (%d across layers with nothing between them), contact floor %s)\n", t * 1e3, ls.rsmin * 1e3,
                     ls.gmin < 1e30 ? ls.gmin * 1e3 : -1.0, ls.corridors.size(), ls.touching, ls.touchingAcrossLayers,
                     contact_floor ? "ON" : "off");
        return t;
    };
    // SLIVER RISK, PREDICTED BEFORE ANY GEOMETRY IS BUILT (2026-09-19, ABT #1257). Micron-scale
    // slivers come from surfaces that nearly touch: where the clearance g between emitted surfaces
    // of radius ~r is much smaller than the element size, the corridor triangulates into degenerate
    // tets and the AMS solve stops converging -- or worse, returns a spurious answer
    // (10_emi_filter: 318 slivers, 99 W instead of 0.2 W). This runs on the enriched layout alone,
    // so it costs milliseconds and precedes the hours of CAD + meshing. OMFEM_PREDICT_ONLY=1
    // prints it and stops.
    //
    // WHAT THE REFUSAL TESTS NOW (2026-09-20). It used to test the requested target against the
    // GLOBAL minimum bound -- the right test only while the mesher had a single size for all
    // conductors. With the corridors refined locally that test answers the wrong question: a
    // target above the tightest corridor's bound is no longer degenerate, because that corridor
    // gets its own size. What still cannot be meshed is a corridor whose OWN bound is below the
    // feasibility floor -- no local field can rescue a channel that demands sub-20-um elements --
    // so THAT is what is refused here, per corridor, naming the corridor. The global-vs-target
    // comparison survives only under OMFEM_CHORD_GLOBAL=1, where the operator has explicitly
    // asked for the old single-size behaviour and there is no local field to save it.
    if (opt.conductor_target > 0.0 && !std::getenv("OMFEM_NO_SLIVER_PREDICT")) {
        // The measurement itself may legitimately fail (a layout it cannot measure) -- that is
        // caught and reported. The REFUSALS must not be: they are thrown AFTER the try block,
        // because a throw inside it was swallowed by this very catch and printed as "not
        // predictable", which read like a harmless diagnostic while the mesher carried on into
        // the degenerate mesh. Verified 2026-09-19: the guard fired, was caught here, and meshing
        // continued as if nothing had happened.
        const bool global_mode = std::getenv("OMFEM_CHORD_GLOBAL") != nullptr;
        bool infeasible = false, above_bound = false;
        double worst = 0.0, worst_x = 0.0, worst_y = 0.0, worst_g = 0.0;
        std::string worst_pair;
        double bound = 0.0, ratio = 0.0;
        size_t n_local = 0;
        try {
            const LayoutSizing& ls = layout_sizing();
            for (const auto& c : ls.corridors) {
                if (c.size >= opt.conductor_target) continue;   // the global target already fits
                ++n_local;
                if (n_local == 1 || c.size < worst)
                    { worst = c.size; worst_x = c.x; worst_y = c.y; worst_g = c.clearance; worst_pair = "'" + c.turnA + "' / '" + c.turnB + "'"; }
            }
            infeasible = (n_local > 0 && worst < 2e-5);
            std::fprintf(stderr, "[mesh3d] chord bound: requested conductor target %.4g mm; %zu of "
                         "%zu measured corridor(s) are tighter and will be refined LOCALLY"
                         "%s\n", opt.conductor_target * 1e3, n_local, ls.corridors.size(),
                         n_local ? "" : " (none: the requested target already fits every corridor)");
            if (n_local)
                std::fprintf(stderr, "[mesh3d] tightest corridor: clearance %.4g mm at layout "
                             "(%.4g, %.4g) mm between %s -> local size %.4g mm\n",
                             worst_g * 1e3, worst_x * 1e3, worst_y * 1e3, worst_pair.c_str(), worst * 1e3);
            if (global_mode) {
                bound = ls.corridors.empty() ? 1.2 * ls.rsmin : ls.corridors.front().size;
                ratio = opt.conductor_target / bound;
                above_bound = (ratio > 1.0);
                std::fprintf(stderr, "[mesh3d] OMFEM_CHORD_GLOBAL: one size for every conductor -- "
                             "target %.4g mm vs global bound %.4g mm -> ratio %.2f%s\n",
                             opt.conductor_target * 1e3, bound * 1e3, ratio,
                             above_bound ? "  *** ABOVE THE BOUND ***" : "");
            }
        } catch (const std::exception& e) {
            std::fprintf(stderr, "[mesh3d] chord bound: not predictable (%s)\n", e.what());
        }
        // REFUSE rather than mesh known-degenerate input (Alf 2026-09-19). A corridor under the
        // feasibility floor triangulates into degenerate tets whatever size field is applied;
        // that is what gmsh's HXT kernel turns into a heap corruption and what leaves slivers
        // behind when it does not crash. The mesher does NOT substitute a size of its own -- it
        // names the corridor and stops, so the choice is visible in the corpus row.
        if (infeasible && !std::getenv("OMFEM_ALLOW_ABOVE_CHORD_BOUND"))
            throw std::runtime_error(
                "mesh3d: the corridor at layout (" + std::to_string(worst_x * 1e3) + ", " +
                std::to_string(worst_y * 1e3) + ") mm between " + worst_pair + " has a clearance of " +
                std::to_string(worst_g * 1e6) + " um and demands " + std::to_string(worst * 1e3) +
                " mm elements -- under the 20 um feasibility floor, so no local refinement can "
                "mesh it. Fix the layout, or set OMFEM_ALLOW_ABOVE_CHORD_BOUND=1 to mesh it "
                "anyway. (No CHORD_BOUND_MM is emitted: this is not a target the caller can "
                "lower into -- re-meshing at 0.95x a sub-20-um bound is not a mesh.)");
        if (above_bound && !std::getenv("OMFEM_ALLOW_ABOVE_CHORD_BOUND"))
            throw std::runtime_error(
                "mesh3d: conductor target " + std::to_string(opt.conductor_target * 1e3) +
                " mm is ABOVE the layout's global chord bound " + std::to_string(bound * 1e3) +
                " mm (ratio " + std::to_string(ratio) + ") and OMFEM_CHORD_GLOBAL disabled the "
                "local refinement that would otherwise handle it. Choose a target below the "
                "bound, drop OMFEM_CHORD_GLOBAL, or set OMFEM_ALLOW_ABOVE_CHORD_BOUND=1. "
                "CHORD_BOUND_MM=" + std::to_string(bound * 1e3));
        if (std::getenv("OMFEM_PREDICT_ONLY")) { std::fprintf(stderr, "[mesh3d] OMFEM_PREDICT_ONLY: stopping after the prediction\n"); return opt.out_msh; }
    }
    // What MVB++'s enamel gate proved about THIS build (NotRun unless the build completed). Only
    // Certified is a proof that no two conductors share volume; the same-region weld below relies
    // on it and on nothing weaker.
    mvb::EnamelGateVerdict gateVerdict = mvb::EnamelGateVerdict::NotRun;
    std::vector<mvb::NamedShape> named =
        builder.buildAllNamed(enriched, include_bobbin, symPlanes,
                              opt.polygon_segments, opt.polygon_segments, opt.paint_coating,
                              emit_coating, incl_insul, core_coat_t,
                              /*useRealWindingGeometry=*/continuous, /*femReady=*/continuous,
                              /*skipGeometryChecks=*/false, &gateVerdict);
    if (named.empty()) throw std::runtime_error("mesh3d_from_mas: MVB++ produced no solids");
    if (!accessories) {
        // (b) of ACCESSORIES OUT above. Solder and FR4 are NOT accessories in this sense: solder is the
        // copper path of a foil terminal and FR4 is a planar coil's own substrate, so both stay.
        std::vector<mvb::NamedShape> kept; kept.reserve(named.size());
        for (auto& ns : named) {
            switch (ns.role) {
                case mvb::Role::Bobbin: case mvb::Role::Pin: case mvb::Role::Divider: case mvb::Role::Base:
                case mvb::Role::Spacer: case mvb::Role::Shunt: case mvb::Role::Sleeve:
                    std::fprintf(stderr, "[mesh3d] accessory OUT: '%s' (%s)\n", ns.name.c_str(), mvb::role_name(ns.role));
                    continue;
                default:
                    kept.push_back(std::move(ns));
            }
        }
        named.swap(kept);
        if (named.empty()) throw std::runtime_error("mesh3d_from_mas: nothing left after dropping the accessory solids");
    }
    std::fprintf(stderr, "[mesh3d] %zu solids (segments=%d coating=%d symmetry=%s/%d planes)\n",
                 named.size(), opt.polygon_segments, (int)opt.paint_coating,
                 opt.symmetry.c_str(), symPlanes);
    // TERMINAL LEAD COPPER (ABT #1215, 2026-09-13): the 3D copper includes the leads MVB++ draws, MKF's
    // R_dc counts turns only (Alf: keep the leads in 3D; MKF will compute its own lead length, MVB++
    // reports what it drew as the cross-check). Write MVB++'s measurement next to the mesh as
    // <out.msh>.leads.json so the corpus runner can state the lead share of the copper in its row.
    // Real-winding path only; a failed measurement fails the run (no silent absence).
    if (continuous) {
      // The lead measurement is a CROSS-CHECK for MKF's own lead length (ABT #1215/#1217), not an input to the mesh or
      // the solve. When it cannot be made (00_debug 2026-09-17: "'Primary parallel 0' has no turn copper"), record the
      // failure in the sidecar and on stderr and carry on meshing -- refusing the mesh over a diagnostic is the wrong
      // trade. Nothing is swallowed: the sidecar says why it is empty.
      try {
        const auto leads = builder.measureTerminalLeadLengths(
              enriched, emit_coating ? false : opt.paint_coating, /*femReady=*/continuous,
              opt.polygon_segments, opt.polygon_segments, core_coat_t);
          const nlohmann::json lj = mvb::MagneticBuilder::terminalLeadLengthsToJson(leads);
          const std::string lp = opt.out_msh + ".leads.json";
          std::ofstream lf(lp);
          if (!lf) throw std::runtime_error("mesh3d_from_mas: cannot write " + lp);
          lf << lj.dump(2) << "\n";
          for (const auto& kv : leads)
              std::fprintf(stderr, "[mesh3d] terminal leads %s: %.3f mm total over %d parallel(s)\n",
                           kv.first.c_str(), kv.second.total_m * 1e3, (int)kv.second.parallels);
        std::fprintf(stderr, "[mesh3d] lead lengths: %s\n", lp.c_str());
      } catch (const std::exception& e) {
        const std::string lp = opt.out_msh + ".leads.json";
        std::ofstream lf(lp);
        if (lf) lf << nlohmann::json{{"error", e.what()}}.dump(2) << "\n";
        // "has no turn copper" is not a measurement problem: it is MVB++ telling us the conductor
        // has NO TURNS -- only its two lead stubs were emitted (a one-turn conductor on a ROUND
        // column, MVB++ ABT #1260). Meshing on is pointless and actively misleading: 00_debug went
        // on to produce a mesh whose $PhysicalNames carry only air, core and outer, no winding and
        // no ports, which reads downstream as a meshing failure and sends the search to the wrong
        // stage. Making this non-fatal was my own shortcut; MVB++ threw and we continued anyway.
        // Other measurement failures stay non-fatal -- a lead length we cannot measure costs us a
        // row note, not the physics.
        if (std::string(e.what()).find("no turn copper") != std::string::npos)
            throw std::runtime_error(
                std::string("mesh3d: MVB++ emitted NO TURN COPPER for this winding -- the conductor "
                            "is lead stubs only, so any mesh built from it has no winding region and "
                            "no ports (MVB++ ABT #1260: a one-turn conductor on a ROUND column gets "
                            "no ring). Refusing to mesh it. Original: ") + e.what());
        std::fprintf(stderr, "[mesh3d] WARNING terminal-lead measurement failed (%s); recorded in %s, meshing continues\n",
                     e.what(), lp.c_str());
      }
    }
    // OMFEM_MESH_DEBUG: what MVB++ actually handed us, BEFORE any STEP round-trip / import /
    // region matching. Without this the geometry stage is a black box -- a winding that arrives
    // as a lead stub instead of a coil looks identical to a solver bug downstream.
    if (std::getenv("OMFEM_MESH_DEBUG")) {
        for (const auto& ns : named) {
            GProp_GProps gp; BRepGProp::VolumeProperties(ns.shape, gp);
            Bnd_Box bb; BRepBndLib::Add(ns.shape, bb);
            double x0,y0,z0,x1,y1,z1;
            if (bb.IsVoid()) { x0=y0=z0=x1=y1=z1=0.0; } else bb.Get(x0,y0,z0,x1,y1,z1);
            int nsol = 0; for (TopExp_Explorer ex(ns.shape, TopAbs_SOLID); ex.More(); ex.Next()) ++nsol;
            std::fprintf(stderr, "[mesh3d.dbg] named '%s': vol=%.4gmm3 solids=%d "
                         "bbox x[%.3f,%.3f] y[%.3f,%.3f] z[%.3f,%.3f] mm\n",
                         ns.name.c_str(), gp.Mass()*1e9, nsol,
                         x0*1e3,x1*1e3, y0*1e3,y1*1e3, z0*1e3,z1*1e3);
        }
    }
    // OMFEM_VERIFY_CAD=1: per-solid WATERTIGHTNESS battery on everything MVB++ handed over,
    // BEFORE any meshing -- BRepCheck validity, every shell topologically CLOSED, positive
    // volume. One "[verify-cad]" line per named shape plus a final verdict line the sweep
    // runner can gate on. Reports rather than throws: the runner decides, and a failed check
    // must not mask the mesh-stage diagnosis that would otherwise follow.
    if (std::getenv("OMFEM_VERIFY_CAD")) {
        int bad = 0;
        for (const auto& ns : named) {
            if (ns.name.find(" terminal ") != std::string::npos) continue;   // port cap faces
            int nsol = 0, openShells = 0;
            for (TopExp_Explorer ex(ns.shape, TopAbs_SOLID); ex.More(); ex.Next()) {
                ++nsol;
                for (TopExp_Explorer sh(ex.Current(), TopAbs_SHELL); sh.More(); sh.Next())
                    if (!TopoDS::Shell(sh.Current()).Closed() &&
                        !BRep_Tool::IsClosed(sh.Current()))
                        ++openShells;
            }
            GProp_GProps gp; BRepGProp::VolumeProperties(ns.shape, gp);
            const bool valid = BRepCheck_Analyzer(ns.shape).IsValid();
            const bool ok = valid && openShells == 0 && nsol >= 1 && gp.Mass() > 1e-15;
            if (!ok) ++bad;
            std::fprintf(stderr, "[verify-cad] %s '%s': valid=%d openShells=%d solids=%d "
                         "vol=%.6gmm3\n", ok ? "OK  " : "FAIL", ns.name.c_str(), (int)valid,
                         openShells, nsol, gp.Mass() * 1e9);
        }
        std::fprintf(stderr, "[verify-cad] VERDICT: %s (%d bad shape(s))\n",
                     bad == 0 ? "ALL WATERTIGHT" : "FAILURES", bad);
    }

    // OMFEM_STEP_ONLY: dump the real-winding 3D solids to STEP for VISUAL review and stop (no mesh).
    // This is the exact geometry buildAllNamed produced (core + bobbin + one continuous conductor per
    // winding-parallel), and it sidesteps the swept-winding tet mesher entirely. Terminal cap faces
    // are dropped (degenerate sheets that clutter a CAD viewer); the winding solids carry them anyway.
    if (std::getenv("OMFEM_STEP_ONLY")) {
        std::vector<mvb::NamedShape> solids;
        solids.reserve(named.size());
        for (const auto& ns : named)
            if (ns.name.find(" terminal ") == std::string::npos) solids.push_back(ns);
        std::ostringstream sink; std::streambuf* old = std::cout.rdbuf(sink.rdbuf());
        // The FEM product's periodic-face -> B-spline pass runs inside the exporter, on the mm
        // geometry (MVB++ StepExportOptions; ABT #1111).
        mvb::StepExportOptions xo; xo.nurbsPeriodicSolids = continuous;
        const bool ok = mvb::exportSTEP(solids, opt.out_msh, xo);
        std::cout.rdbuf(old);
        if (!ok) throw std::runtime_error("mesh3d_from_mas(STEP_ONLY): exportSTEP failed");
        std::fprintf(stderr, "[mesh3d] STEP_ONLY: wrote %zu solids -> %s\n", solids.size(), opt.out_msh.c_str());
        return opt.out_msh;
    }

    // Per-solid VOLUME + lumped role (core / winding). We match imported volumes by
    // VOLUME, not centroid: 3D turns are rings centred on the axis, so every turn (and
    // the core) has its centroid ~on the axis -> nearest-centroid collapses them. Volume
    // is distinctive (core pieces are large, turns small; ties among turns share a role).
    // Per-solid signature: VOLUME + bounding-box centre/extent. Turns get an INDIVIDUAL region
    // "turn_<winding>_<index>" (not a single lumped "winding") so the thermal FEA can place the
    // real per-turn loss and resolve the hot-spot turn. We match imported volumes by bbox, not by
    // volume alone: identical wire turns share a volume, but each sits at a distinct axial/radial
    // position so their bboxes differ. (OMFEM_LUMP_WINDING falls back to one "winding" region.)
    const bool lump_winding = std::getenv("OMFEM_LUMP_WINDING") != nullptr;
    const NameContext nctx3d = name_context(enriched);   // ABT #1169
    struct NV { double cx,cy,cz, ex,ey,ez; std::string region; int group; };
    std::vector<NV> nvs;
    // Terminal PORT faces from femReady single-body windings ("<winding> parallel <p> terminal <k>"):
    // planar caps at the two free ends of the continuous conductor. Collect their centres so the meshed
    // winding-boundary face at each can be tagged term<2*port+k> -> the full-wave A-V solver's per-
    // winding current-injection surfaces. Port index is per WINDING (all its parallels share the net).
    struct TermFace { double cx,cy,cz; int idx; double rad; double nx,ny,nz; std::string name; };   // rad = cap bounding radius; n = cap plane normal
    std::vector<TermFace> termFaces;
    std::map<std::string,int> termPort;
    struct TermCand { std::string name, winding, ownBody; int side = 0; TopoDS_Shape shape; };
    std::vector<TermCand> termCands;
    const std::string coatingSuffix = " coating";
    int group = -1;
    for (const auto& ns : named) {
        ++group;   // index of the NAMED SHAPE this signature came from
        std::string name = ns.name;
        if (auto tp = name.find(" terminal "); tp != std::string::npos) {   // a port cap face, not a volume
            // Collected here, decided below: MVB++ caps BOTH path ends of every conductor path it emits, and a foil
            // parallel is three paths -- the sheet and its entrance / exit leads -- so it carries six caps, of which
            // only the two lead TIPS are free ends. The others sit in the solder joints (two_switch_forward: 16 sheet-
            // end caps of 12.12 mm were tagged as ports, ABT #1400).
            TermCand tc;
            tc.name = name;
            tc.winding = name.substr(0, name.find(" parallel "));
            tc.ownBody = name.substr(0, tp);
            tc.side = std::atoi(name.c_str() + tp + 10);
            static const std::regex leadRe(R"(^.+ parallel \d+ (entrance|exit) lead$)");
            std::smatch lm;
            if (std::regex_match(tc.ownBody, lm, leadRe)) tc.side = (lm[1].str() == "entrance") ? 0 : 1;
            tc.shape = ns.shape;
            termCands.push_back(tc);
            continue;
        }
        if (name.rfind("insulation", 0) == 0) {   // insulation_layer_<i> -> own low-k region
            Bnd_Box bb; BRepBndLib::Add(ns.shape, bb);
            double x0,y0,z0,x1,y1,z1; bb.Get(x0,y0,z0,x1,y1,z1);
            nvs.push_back({0.5*(x0+x1),0.5*(y0+y1),0.5*(z0+z1), x1-x0,y1-y0,z1-z0, name, group});
            continue;
        }
        // "<turn> coating" = the outer insulated shell (enamel/insulation/serving/tape -- the type
        // only changes its thickness, which is baked into the outer footprint). Tag it separately.
        bool is_coating = name.size() > coatingSuffix.size() &&
            name.compare(name.size() - coatingSuffix.size(), coatingSuffix.size(), coatingSuffix) == 0;
        if (is_coating) name = name.substr(0, name.size() - coatingSuffix.size());
        std::string winding; int index = 0;
        std::string role = classify(name, winding, index, nctx3d);
        if (name.rfind("Bobbin", 0) == 0) role = "bobbin";  // classify() ignores it; keep for thermal
        if (role.empty()) continue;
        std::string region =
              is_coating ? (role == "core" ? ("coating_core_" + std::to_string(index))
                                            : ("coating_" + winding + "_" + std::to_string(index)))
            : (role == "core") ? "core" : (role == "bobbin") ? "bobbin"
            // ABT #1169: accessory solids. A pin is copper and a base is plastic, but neither is
            // magnetic and neither carries a driven current, so both share the bobbin's bucket
            // until they get their own conductivity; the rest are their own regions so a
            // per-region material (D4/WP0) can be attached without touching this tagger again.
            : (role == "pin")     ? "bobbin"
            : (role == "spacer")  ? ("spacer_"  + std::to_string(index))
            : (role == "shunt")   ? ("shunt_"   + std::to_string(index))
            : (role == "divider") ? ("divider_" + std::to_string(index))
            : (role == "sleeve")  ? ("sleeve_"  + winding + "_" + std::to_string(index))
            : (role == "solder")  ? ("solder_"  + winding + "_" + std::to_string(index))
            // Continuous conductor (real winding): all parallels of a winding lump into ONE coil
            // region "winding_<name>" so the magnetostatic solver drives it as a single stranded
            // coil (its "winding*" prefix scan). Distinct windings -> distinct coils.
            : (role == "winding") ? ("winding_" + winding)
            : lump_winding ? "winding" : ("turn_" + winding + "_" + std::to_string(index));
        // One signature PER SOLID, not per named shape. Imported volumes are classified by nearest
        // bbox centre+extent, so a multi-solid shape (the bobbin's 3 pieces; a real winding that
        // comes back as a chain of per-primitive conductor solids) must contribute one signature
        // each -- otherwise its pieces are matched against the bbox of the WHOLE assembly, which is
        // nowhere near any of them, and they get classified as whatever happens to be nearest.
        // (Measured before this: the bobbin's own pieces matched at L1 = 9-13 mm.)
        int nsolids = 0;
        for (TopExp_Explorer ex(ns.shape, TopAbs_SOLID); ex.More(); ex.Next()) {
            Bnd_Box sb; BRepBndLib::Add(ex.Current(), sb);
            if (sb.IsVoid()) continue;
            double a0,b0,c0,a1,b1,c1; sb.Get(a0,b0,c0,a1,b1,c1);
            nvs.push_back({0.5*(a0+a1),0.5*(b0+b1),0.5*(c0+c1), a1-a0,b1-b0,c1-c0, region, group});
            ++nsolids;
        }
        if (nsolids == 0) {   // sheet/shell shape with no solid: keep the whole-shape signature
            Bnd_Box bb; BRepBndLib::Add(ns.shape, bb);
            if (bb.IsVoid()) continue;
            double x0,y0,z0,x1,y1,z1; bb.Get(x0,y0,z0,x1,y1,z1);
            nvs.push_back({0.5*(x0+x1),0.5*(y0+y1),0.5*(z0+z1), x1-x0,y1-y0,z1-z0, region, group});
        }
    }
    if (nvs.empty()) throw std::runtime_error("mesh3d_from_mas: no core/turn solids classified");
    // PORTS = FREE conductor ends only (ABT #1400). A terminal cap that touches another conductor or solder body is a
    // JOINT (a foil sheet end in its solder film, a lead's inner end on the sheet), not a place where current enters
    // the model. Decided by exact contact (BRepExtrema) with a tolerance of 1e-3 of the cap's own radius -- a round
    // wire's lead tip touches nothing and stays a port. Port numbering is per winding, entrance/start 0, exit/end 1.
    {
        std::vector<std::pair<std::string, const TopoDS_Shape*>> bodies;   // conductor + solder bodies
        for (const auto& ns : named) {
            if (ns.name.find(" terminal ") != std::string::npos) continue;
            std::string w; int ix = 0;
            const std::string role = classify(ns.name, w, ix, nctx3d);
            if (role == "winding" || role == "solder" || role == "turn") bodies.push_back({ns.name, &ns.shape});
        }
        // A parallel that has its own LEAD bodies ("<w> parallel <p> entrance|exit lead", the foil construction) takes its
        // ports at the lead tips; the sheet's own two caps are its free EDGES -- dead ends the leads are soldered beside,
        // where no current enters or leaves -- so they are not ports either (two_switch: 12.12 mm caps at the sheet ends).
        std::set<std::string> parallelsWithLeads;
        {
            static const std::regex leadBodyRe(R"(^(.+ parallel \d+) (?:entrance|exit) lead$)");
            for (const auto& tc : termCands) {
                std::smatch lm;
                if (std::regex_match(tc.ownBody, lm, leadBodyRe)) parallelsWithLeads.insert(lm[1].str());
            }
        }
        std::map<std::string, std::array<int, 2>> portsPerWinding;
        for (const auto& tc : termCands) {
            if (parallelsWithLeads.count(tc.ownBody)) {
                std::fprintf(stderr, "[mesh3d] terminal cap '%s' is a sheet edge of a parallel with lead bodies -> not a port\n", tc.name.c_str());
                continue;
            }
            Bnd_Box bb; BRepBndLib::Add(tc.shape, bb);
            double x0,y0,z0,x1,y1,z1; bb.Get(x0,y0,z0,x1,y1,z1);
            const double rad = 0.5*std::sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
            const double tol = 1e-3 * rad;
            std::string joinedTo;
            for (const auto& [bn, bs] : bodies) {
                if (bn == tc.ownBody) continue;
                Bnd_Box ob; BRepBndLib::Add(*bs, ob); ob.Enlarge(tol);
                if (ob.IsOut(bb)) continue;
                BRepExtrema_DistShapeShape dist(tc.shape, *bs, Extrema_ExtFlag_MIN, Extrema_ExtAlgo_Tree);
                if (!dist.IsDone())
                    throw std::runtime_error("mesh3d_from_mas: cannot measure the contact of terminal cap '" + tc.name + "' with '" + bn + "'");
                if (dist.Value() <= tol) { joinedTo = bn; break; }
            }
            if (!joinedTo.empty()) {
                std::fprintf(stderr, "[mesh3d] terminal cap '%s' joins '%s' -> a joint, not a port\n", tc.name.c_str(), joinedTo.c_str());
                continue;
            }
            const int port = termPort.emplace(tc.winding, (int)termPort.size()).first->second;
            double nx = 0, ny = 0, nz = 0;
            for (TopExp_Explorer ex(tc.shape, TopAbs_FACE); ex.More(); ex.Next()) {
                BRepAdaptor_Surface surf(TopoDS::Face(ex.Current()));
                if (surf.GetType() != GeomAbs_Plane) continue;
                const gp_Dir d = surf.Plane().Axis().Direction(); nx = d.X(); ny = d.Y(); nz = d.Z(); break;
            }
            if (nx == 0 && ny == 0 && nz == 0)
                throw std::runtime_error("mesh3d_from_mas: terminal cap '" + tc.name + "' has no planar face -- cannot orient its port");
            termFaces.push_back({0.5*(x0+x1),0.5*(y0+y1),0.5*(z0+z1), 2*port + tc.side, rad, nx, ny, nz, tc.name});
            portsPerWinding[tc.winding][tc.side]++;
        }
        for (const auto& [w, n] : portsPerWinding)
            if (n[0] == 0 || n[1] == 0)
                throw std::runtime_error("mesh3d_from_mas: winding '" + w + "' has " + std::to_string(n[0]) + " entrance and " +
                                         std::to_string(n[1]) + " exit port cap(s) after removing the joints -- it needs both ends");
    }

    // STEP round-trip (heals geometry for gmsh's kernel, as in the 2D path). Terminal cap faces are
    // NOT exported — they are already faces of the winding SOLID; a standalone duplicate would confuse
    // the import/fragment. We tag them post-mesh from the collected centres (termFaces) instead.
    const std::string step_path = opt.out_msh + ".step";
    {
        std::vector<mvb::NamedShape> forStep;
        forStep.reserve(named.size());
        for (const auto& ns : named)
            if (ns.name.find(" terminal ") == std::string::npos) forStep.push_back(ns);
        std::ostringstream sink;
        std::streambuf* old = std::cout.rdbuf(sink.rdbuf());
        // The FEM product's periodic-face -> B-spline pass runs inside the exporter, on the mm
        // geometry (MVB++ StepExportOptions; ABT #1111).
        mvb::StepExportOptions xo; xo.nurbsPeriodicSolids = continuous;
        const bool ok = mvb::exportSTEP(forStep, step_path, xo);
        std::cout.rdbuf(old);
        if (!ok) throw std::runtime_error("mesh3d_from_mas: exportSTEP failed");
    }

    bool own_init = !gmsh::isInitialized();
    if (own_init) gmsh::initialize(); apply_gmsh_threads();
    // OMFEM_GMSH_VERBOSE: let gmsh name the entity it choked on. Without it a mesh failure is a
    // bare "the 1D mesh seems not to be forming a closed loop" with no surface/curve tag.
    const bool gmsh_verbose = std::getenv("OMFEM_GMSH_VERBOSE") != nullptr;
    gmsh::option::setNumber("General.Terminal", gmsh_verbose ? 1 : 0);
    if (gmsh_verbose) gmsh::option::setNumber("General.Verbosity", 99);
    // NO SILENT ALGORITHM SWITCH. When the 2D surface mesher fails on a surface, gmsh by default
    // falls back to MeshAdapt -- whose BDS edge-collapse SEGFAULTS on these conductor surfaces
    // (reproduced on PQ5050: BDS_Mesh::add_edge <- collapse_edge_parametric, killing the process
    // with no error at all; ETD34/ETD49/ETD3920 died the same way). A crash is strictly worse than
    // a failure: it cannot be caught, classified or retried. With the switch off the surface fails
    // CLEANLY and the retry ladder gets to do its job. OMFEM_ALGO_SWITCH=1 restores the old
    // behaviour for comparison.
    if (!std::getenv("OMFEM_ALGO_SWITCH"))
        gmsh::option::setNumber("Mesh.AlgorithmSwitchOnFailure", 0);
    // FACET-OVERLAP ANGLE TOLERANCE. gmsh flags two facets as "overlapping" when they are within
    // this dihedral angle (default 0.1 deg) AND project onto each other -- a heuristic aimed at
    // genuinely folded surfaces. A dense winding VIOLATES its premise by construction: adjacent
    // turns present near-PARALLEL sheets of the same lateral surface a few tens of microns apart,
    // so the heuristic reports "overlapping facets on surface N surface N" for geometry that is
    // provably non-intersecting (measured: the flagged pair on ETD34 met at 0.0077 deg, and the
    // CAD clearance there is 27.5 um against a 3 um chord error). Tighten it for the whole
    // pipeline so the DELAUNAY stage judges real intersections rather than parallel proximity.
    // OMFEM_FACET_ANGLE_TOL overrides (degrees).
    double facet_tol = 1e-4;
    if (const char* ft = std::getenv("OMFEM_FACET_ANGLE_TOL")) {
        char* ftEnd = nullptr;
        facet_tol = std::strtod(ft, &ftEnd);
        if (ftEnd == ft || !(facet_tol > 0.0))
            throw std::runtime_error("mesh3d_from_mas: OMFEM_FACET_ANGLE_TOL must be a positive "
                                     "number in degrees (got '" + std::string(ft) + "')");
    }
    gmsh::option::setNumber("Mesh.AngleToleranceFacetOverlap", facet_tol);
    // The option is PROCESS-GLOBAL gmsh state; restore the stock default on every exit
    // (success or throw) so it cannot leak into a caller's later, winding-less meshes. A session
    // this function opened is finalized before the guard runs, and takes the option with it.
    struct FacetTolRestore {
        ~FacetTolRestore() {
            if (gmsh::isInitialized()) gmsh::option::setNumber("Mesh.AngleToleranceFacetOverlap", 0.1);
        }
    } facet_tol_restore;
    // Optional geometry healing (env OMFEM_HEAL=1) for faceted turns. NOTE: OCCSewFaces
    // /OCCMakeSolids are too aggressive here -- they merge/destroy the separate turn solids
    // (-> "no volumes classified"), so healing is OFF by default. Kept env-gated for
    // experiments; the real fix for faceting is sizing, not healing.
    if (std::getenv("OMFEM_HEAL")) {
        gmsh::option::setNumber("Geometry.OCCFixSmallEdges", 1);
        gmsh::option::setNumber("Geometry.OCCFixSmallFaces", 1);
        gmsh::option::setNumber("Geometry.Tolerance", 1e-6);
    }
    gmsh::model::add("omfem_mas3d");
    gmsh::vectorpair imported;
    // MVB++ exportSTEP owns the metres->millimetres conversion (ABT #317), so the file on disk
    // declares MILLIMETRES while everything in OMFEM -- mesh size targets, sigma, mu0, the
    // solvers -- is SI metres. Without this the round-trip silently returns geometry 1000x too
    // big (measured: a core bbox of x[-6.325,6.325] mm came back as x[-6325,6325]), which also
    // wrecks the bbox-based region matching. Ask OCC to convert to metres on import so the unit
    // declared in the file is honoured whatever the exporter does.
    gmsh::option::setString("Geometry.OCCTargetUnit", "M");
    gmsh::model::occ::importShapes(step_path, imported);
    gmsh::model::occ::synchronize();

    // Classify each imported volume by nearest recorded centroid (geometry, not labels).
    // label = the STEP PRODUCT the volume came from ("Shapes/Primary parallel 1"), read back through gmsh's
    // label import. It is the exact identity of the MVB++ named shape; the bbox match below cannot separate
    // interleaved parallels whose boxes agree to 0.05 mm (23_llc: both primary parallels matched one shape).
    struct Imp { int tag; std::string region; double vol; int group; std::string label; };
    std::vector<Imp> imps;
    {
        gmsh::vectorpair vols; gmsh::model::occ::getEntities(vols, 3);
        for (auto& dt : vols) {
            double m = 0.0; gmsh::model::occ::getMass(3, dt.second, m);
            if (m <= 1e-15) continue;
            double x0,y0,z0,x1,y1,z1; gmsh::model::occ::getBoundingBox(3, dt.second, x0,y0,z0,x1,y1,z1);
            const double cx=0.5*(x0+x1),cy=0.5*(y0+y1),cz=0.5*(z0+z1), ex=x1-x0,ey=y1-y0,ez=z1-z0;
            double best = DBL_MAX; const std::string* region = nullptr;
            for (const auto& nv : nvs) {   // nearest by bbox centre+extent -> distinguishes individual turns
                const double d = std::fabs(cx-nv.cx)+std::fabs(cy-nv.cy)+std::fabs(cz-nv.cz)
                               + std::fabs(ex-nv.ex)+std::fabs(ey-nv.ey)+std::fabs(ez-nv.ez);
                if (d < best) { best = d; region = &nv.region; group = nv.group; }
            }
            std::string label; gmsh::model::getEntityName(3, dt.second, label);
            if (region) imps.push_back({dt.second, *region, m, group, label});
            if (std::getenv("OMFEM_MESH_DEBUG"))
                std::fprintf(stderr, "[mesh3d.dbg] imported vol tag=%d '%s' mass=%.4gmm3 -> region '%s' group %d "
                             "(match L1=%.4gmm) bbox x[%.3f,%.3f] y[%.3f,%.3f] z[%.3f,%.3f] mm\n",
                             dt.second, label.c_str(), m*1e9, region?region->c_str():"<none>", group, best*1e3,
                             x0*1e3,x1*1e3, y0*1e3,y1*1e3, z0*1e3,z1*1e3);
        }
    }
    if (imps.empty()) { if (own_init) gmsh::finalize(); throw std::runtime_error("mesh3d_from_mas: no volumes classified"); }
    std::fprintf(stderr, "[mesh3d] imported %zu solids\n", imps.size());

    // Union OVERLAPPING same-region solids. A winding's parallel strands (allowed to overlap now that
    // same-net envelope overlaps are permitted, ConductorBuilder) arrive as separate solids whose
    // volumes intersect -> the mesher sees overlapping facets and fails. They are one electrical net,
    // so fuse them into one solid per region (real volume overlap unions robustly). If a region's fuse
    // fails (e.g. a toroidal winding's coincident mitre faces, not a true overlap), keep its pieces.
    // Single-solid regions (a clean single-body winding) pass straight through.
    //
    // Only pieces of ONE conductor are candidates when MVB++ certified the build. A winding's
    // parallels are distinct conductors, and the enamel gate proved every pair of those at or
    // beyond its coated envelope: they share no volume and CANNOT fuse into one solid, so trying
    // is pure cost. Across the corpus the parallel weld never merged anything (0 of 1,110 steps,
    // every one nsolid=2 with the masses adding up), and on the all-B-spline windings of the LLC
    // transformers the first such fuse ran past 44 min without finishing (23_llc, 2026-09-26).
    // Pieces of the SAME named conductor (one MVB++ left multi-body) are still welded. Identity is
    // the imported STEP product label, never the bbox match; a solid without a label, or any verdict
    // other than Certified, keeps the old region-wide weld.
    const bool gateCertified = gateVerdict == mvb::EnamelGateVerdict::Certified;
    {
        static const char* const kVerdict[] = {"NotRun", "Certified", "SampledOnly", "ReportedAllowed", "Skipped"};
        const int vi = static_cast<int>(gateVerdict);
        std::fprintf(stderr, "[mesh3d] MVB++ enamel gate verdict: %s\n",
                     (vi >= 0 && vi < 5) ? kVerdict[vi] : "<unknown>");
    }
    {
        std::map<std::pair<std::string,std::string>,std::vector<int>> byreg;
        // Conductors only: the certificate is about copper. Core halves, bobbin pieces and the rest keep the
        // region-wide weld (an ungapped core's halves DO fuse into one body).
        for (auto& im : imps) {
            const bool conductor = im.region.rfind("winding", 0) == 0 || im.region.rfind("turn_", 0) == 0;
            byreg[{im.region, (gateCertified && conductor) ? im.label : std::string()}].push_back(im.tag);
        }
        bool anymulti=false; for (auto& kv : byreg) if (kv.second.size()>1) { anymulti=true; break; }
        if (gateCertified) {
            std::map<std::string,int> perRegion;
            for (auto& kv : byreg) ++perRegion[kv.first.first];
            int skipped = 0; for (auto& kv : perRegion) if (kv.second > 1) skipped += kv.second;
            if (skipped)
                std::fprintf(stderr, "[mesh3d] enamel gate CERTIFIED: %d distinct conductor(s) "
                             "sharing a winding region are not welded to each other (proven disjoint)\n", skipped);
        }
        if (anymulti) {
            std::vector<Imp> fused;
            for (auto& kv : byreg) {
                const std::string& reg = kv.first.first;
                auto& tags = kv.second;
                if (tags.size()==1) { double m; gmsh::model::occ::getMass(3,tags[0],m); fused.push_back({tags[0],reg,m,0}); continue; }
                std::vector<Imp> got; bool ok=false;
                // ITERATIVE pairwise fusion, biggest piece first, each step COVERAGE-GUARDED.
                // One N-way BRepAlgo fuse can silently DROP an input (measured on ETD34: the
                // 5-piece winding fused to a body missing the whole exit lead -- 2.75 mm3 of
                // copper gone, one FEM port unusable -- while the total mass moved only 0.54%,
                // invisible to any relative guard). Fusing one piece at a time with the inputs
                // kept alive (removeObject/removeTool=false) lets a bad step be REJECTED: the
                // piece that will not weld stays a separate solid for the meshing fragment,
                // and every piece that welds cleanly still does. Coverage per step: the merged
                // body must retain ~the accumulated volume plus the new piece.
                try {
                    // Fuzzy boolean tolerance for the welds: the lead prims are deliberately
                    // EXTENDED into each other (coaxial same-radius overlap closes the mitre),
                    // and coincident cylindrical surfaces are the classic BRepAlgo degenerate
                    // case -- at the default tolerance the fuse "succeeds" while dropping the
                    // piece. Same remedy as MVB++'s own fuseAllSolids (MVB_FUSE_FUZZY): a
                    // micron-scale fuzzy value. Restored right after the region loop.
                    gmsh::option::setNumber("Geometry.ToleranceBoolean", 1e-6);
                    std::vector<std::pair<double,int>> order;
                    for (int t : tags) { double m; gmsh::model::occ::getMass(3,t,m); order.push_back({m,t}); }
                    std::sort(order.rbegin(), order.rend());
                    int acc = order[0].second; double accMass = order[0].first;
                    std::vector<int> loose;
                    // Every tag this weld ever touched (inputs, intermediate accumulators,
                    // rejected fuse outputs). gmsh booleans MAY alias an input tag in their
                    // output, so nothing is removed mid-flight -- one guarded sweep at the end
                    // deletes everything not kept, each tag individually try-caught.
                    std::set<int> everSeen(tags.begin(), tags.end());
                    auto weldOne = [&](int t, double mt) -> bool {
                        bool stepOk = false; int newAcc = acc; double newMass = 0.0;
                        try {
                            gmsh::vectorpair outd; std::vector<gmsh::vectorpair> outm;
                            gmsh::model::occ::fuse({{3,acc}}, {{3,t}}, outd, outm, -1, false, false);
                            gmsh::model::occ::synchronize();
                            int nsolid = 0;
                            for (auto& dt : outd) if (dt.first==3) { double m;
                                gmsh::model::occ::getMass(3,dt.second,m);
                                everSeen.insert(dt.second);
                                if (m>1e-15) { newAcc = dt.second; newMass += m; ++nsolid; } }
                            stepOk = nsolid == 1 && newMass >= 0.98 * accMass + 0.5 * mt
                                     && newMass >= 0.90 * (accMass + mt);
                            if (!stepOk)
                                std::fprintf(stderr, "[mesh3d]   weld diag: nsolid=%d "
                                             "newMass=%.6g accMass=%.6g mt=%.6g mm3\n",
                                             nsolid, newMass*1e9, accMass*1e9, mt*1e9);
                        } catch (...) { stepOk = false; }
                        if (stepOk) { acc = newAcc; accMass = newMass; }
                        return stepOk;
                    };
                    for (size_t ti = 1; ti < order.size(); ++ti)
                        if (!weldOne(order[ti].second, order[ti].first))
                            loose.push_back(order[ti].second);
                    // SECOND PASS: a piece that refused to weld against the early accumulator can
                    // weld once the neighbours it coincides with are in.
                    for (auto it = loose.begin(); it != loose.end();) {
                        double mt; gmsh::model::occ::getMass(3, *it, mt);
                        if (weldOne(*it, mt)) it = loose.erase(it);
                        else {
                            std::fprintf(stderr, "[mesh3d] region '%s': piece tag %d (%.4g mm3) "
                                         "will not weld -- kept separate for fragment welding\n",
                                         reg.c_str(), *it, mt*1e9);
                            ++it;
                        }
                    }
                    std::set<int> keepTags(loose.begin(), loose.end());
                    keepTags.insert(acc);
                    for (int t : everSeen)
                        if (!keepTags.count(t)) {
                            try { gmsh::model::occ::remove({{3,t}}, /*recursive=*/true); }
                            catch (...) {}
                        }
                    gmsh::model::occ::synchronize();
                    { double m; gmsh::model::occ::getMass(3,acc,m); got.push_back({acc,reg,m,0}); }
                    for (int t : loose) { double m; gmsh::model::occ::getMass(3,t,m); got.push_back({t,reg,m,0}); }
                    ok = true;
                } catch (...) { ok=false; got.clear(); }
                gmsh::option::setNumber("Geometry.ToleranceBoolean", 0.0);
                if (ok) for (auto& g : got) fused.push_back(g);
                else    for (int t : tags) { double m; gmsh::model::occ::getMass(3,t,m); fused.push_back({t,reg,m,0}); }
            }
            if (fused.size() != imps.size()) {
                imps.swap(fused);
                std::fprintf(stderr, "[mesh3d] fused overlapping same-region solids -> %zu solids\n", imps.size());
            }
        }
    }

    // MAPPED COPPER: cut every rectangular-wire conductor at its junction rings into 6-face blocks
    // (see the helpers above). Wire dimensions come from MKF's resolved wire, never from the geometry.
    const bool skinMapped = std::getenv("OMFEM_SKIN_MAPPED") != nullptr;
    std::map<std::string, std::pair<double,double>> mappedDims;   // region -> (W, T) in metres
    if (skinMapped) {
        OpenMagnetics::Coil coilM = enriched.get_coil();
        const auto& fdM = coilM.get_functional_description();
        for (size_t i = 0; i < fdM.size(); ++i) {
            auto wire = coilM.resolve_wire(i);
            const auto cw = wire.get_conducting_width(), ch = wire.get_conducting_height();
            if (!(cw && ch))
                throw std::runtime_error("mesh3d_from_mas: OMFEM_SKIN_MAPPED needs a rectangular wire (conductingWidth/Height); winding '"
                                         + fdM[i].get_name() + "' has none -- round-wire O-grids are not implemented");
            const double W = OpenMagnetics::resolve_dimensional_values(cw.value());
            const double T = OpenMagnetics::resolve_dimensional_values(ch.value());
            mappedDims["winding_" + fdM[i].get_name()] = {std::max(W, T), std::min(W, T)};
        }
        std::vector<Imp> split;
        for (auto& im : imps) {
            auto it = mappedDims.find(im.region);
            if (it == mappedDims.end()) { split.push_back(im); continue; }
            for (int b : mapped_split_at_rings(im.tag, it->second.first, it->second.second, im.region)) {
                double m = 0.0; gmsh::model::occ::getMass(3, b, m);
                split.push_back({b, im.region, m, im.group});
            }
        }
        imps.swap(split);
        std::fprintf(stderr, "[mesh3d] mapped copper: %zu solids after the ring split\n", imps.size());
    }

    // Far-field air box around all solids.
    double xmin=1e30,ymin=1e30,zmin=1e30,xmax=-1e30,ymax=-1e30,zmax=-1e30;
    { gmsh::vectorpair ents; gmsh::model::occ::getEntities(ents, 3);
      for (auto& dt : ents) { double a,b,c,d,e,f; gmsh::model::occ::getBoundingBox(3,dt.second,a,b,c,d,e,f);
        xmin=std::min(xmin,a); ymin=std::min(ymin,b); zmin=std::min(zmin,c);
        xmax=std::max(xmax,d); ymax=std::max(ymax,e); zmax=std::max(zmax,f); } }
    const double cx=0.5*(xmin+xmax), cy=0.5*(ymin+ymax), cz=0.5*(zmin+zmax);
    const double hx=0.5*(xmax-xmin)*opt.air_margin, hy=0.5*(ymax-ymin)*opt.air_margin, hz=0.5*(zmax-zmin)*opt.air_margin;
    // SYMMETRY CUTS (2026-09-08): with a reduced (half/quarter/eighth) geometry the air box must
    // STOP AT the cut planes so they become box faces -- 'outer' (a x n = 0) for the two planes
    // containing the winding axis (x = 0, z = 0: B lies in them), NOT outer (natural n x H = 0)
    // for the axial mid-plane (y = 0: B is normal to it). Padding the box past the cuts, as this
    // path did, left the quarter coil floating in open air with its electrodes as interior faces:
    // the 03_buck quarter mesh spanned x[-4,20] z[-3.3,16.6] mm -- no symmetry at all. Detected
    // from the reduced geometry (it does not extend to -x / -z / -y), as the level-set path does.
    const bool cutX = symPlanes > 0 && xmin > -1e-4, cutZ = symPlanes > 0 && zmin > -1e-4,
               cutY = symPlanes > 0 && ymin > -1e-4;
    double bx0 = cutX ? 0.0 : cx-hx, by0 = cutY ? 0.0 : cy-hy, bz0 = cutZ ? 0.0 : cz-hz;
    double bx1 = cx+hx, by1 = cy+hy, bz1 = cz+hz;
    // CONDUCTION PATH MUST REACH THE BOUNDARY (2026-09-08). A real winding is NOT a closed loop: it
    // has two free lead ends. The imposed current must be divergence-free over the whole solution
    // space, which (Ansys Maxwell, Eddy Current Excitations) leaves exactly two legal cases -- the
    // conduction path is closed, or it "must begin and end at the boundaries". Our lead tips stopped
    // INSIDE the air box (03_buck: caps at z = -23.22 mm, box to -33.45), so the current had nowhere
    // to go and the bordered system was INCONSISTENT: FGMRES stalled with the field rows at 0.399
    // (their RHS is zero) and the current rows 16% short, and UMFPACK called the gauged matrix
    // singular (ABT #1150). Closed rings never showed it because a ring satisfies the first case.
    // Fix without inventing geometry: pull the air-box face back onto the terminal plane, so each
    // lead END IS a boundary face. OMFEM_NO_PORT_SNAP disables.
    std::vector<int> snapAxis;
    struct SnapPlane { int axis; int side; double plane; double rad; int idx; };
    std::vector<SnapPlane> snaps;
    if (!termFaces.empty() && !std::getenv("OMFEM_NO_PORT_SNAP")) {
        for (const auto& tf : termFaces) {
            struct Cand { double d; double* bound; double val; };
            const Cand c[6] = {{std::fabs(tf.cx-bx0), &bx0, tf.cx}, {std::fabs(tf.cx-bx1), &bx1, tf.cx},
                               {std::fabs(tf.cy-by0), &by0, tf.cy}, {std::fabs(tf.cy-by1), &by1, tf.cy},
                               {std::fabs(tf.cz-bz0), &bz0, tf.cz}, {std::fabs(tf.cz-bz1), &bz1, tf.cz}};
            // FACE BY CAP NORMAL (2026-09-12, ABT #1159): the nearest bound is the wrong criterion
            // for a toroid whose leads exit radially -- a +z-pointing lead sitting near the x+ bound
            // got the x+ face, and a face pulled onto one cap then sliced another lead off ("clipping
            // a lead to the air box removed it entirely" on common_mode_choke_complete). The port must
            // be a cross-section perpendicular to the lead, so the axis is the cap normal's dominant
            // component; only the SIDE on that axis is chosen by distance.
            const double an[3] = {std::fabs(tf.nx), std::fabs(tf.ny), std::fabs(tf.nz)};
            const int axis = (an[0] >= an[1] && an[0] >= an[2]) ? 0 : (an[1] >= an[2] ? 1 : 2);
            const int best = c[2*axis].d <= c[2*axis+1].d ? 2*axis : 2*axis + 1;
            snapAxis.push_back(best);                       // 0/1=x-,x+ 2/3=y-,y+ 4/5=z-,z+
            { static const char* fn[6] = {"x-","x+","y-","y+","z-","z+"};
              std::fprintf(stderr, "[mesh3d] port cap term%d '%s': centre (%.2f,%.2f,%.2f) mm normal (%.4f,%.4f,%.4f) rad %.2f mm -> face %s (pre-snap box x[%.2f,%.2f] y[%.2f,%.2f] z[%.2f,%.2f] mm)\n",
                           tf.idx, tf.name.c_str(), tf.cx*1e3, tf.cy*1e3, tf.cz*1e3, tf.nx, tf.ny, tf.nz, tf.rad*1e3, fn[best],
                           bx0*1e3, bx1*1e3, by0*1e3, by1*1e3, bz0*1e3, bz1*1e3); }
            // Snap PAST the cap, INTO the lead, by one cap radius: the box must SLICE the conductor
            // so the port is a clean planar cross-section lying in the boundary face. Snapping onto
            // the cap itself is not enough when the cap is OBLIQUE (03_buck's lead caps span 12 um
            // in z): the cap then touches the boundary along a line while the rest of the conduction
            // path stays inside, and the system is still ill-posed -- the straight-bar case, whose
            // caps are flush and planar, solves to 5e-13 while 03 blew up to 1e8 V.
            // OMFEM_PORT_INSET (default 1.5): the inset in cap radii; a probe knob for the CMC HXT crash
            // (the snapped face passes 1.5 rad below the caps, i.e. ~1.2 mm from the wrap copper).
            // INSET RULE (2026-09-12, ABT #1159 with mvb-ab): a cap that is already an axis-aligned
            // planar cross-section (|n_axis| > 0.99 -- every MVB++ drop lead) only needs the face to
            // pass 0.5 cap radii into the lead; the 1.5 radii were for OBLIQUE caps (03_buck's 12 um
            // tilt), and on common_mode_choke 1.5 radii put the box face 1.2 mm from the wrap copper,
            // where gmsh's HXT kernel aborted (inset 0.5 and no-snap both mesh). OMFEM_PORT_INSET
            // overrides both.
            const double nAxis = an[axis];
            const double insetDefault = nAxis > 0.99 ? 0.5 : 1.5;
            const double insetFactor = std::getenv("OMFEM_PORT_INSET") ? std::atof(std::getenv("OMFEM_PORT_INSET")) : insetDefault;
            const double inset = insetFactor * tf.rad;
            *c[best].bound = (best % 2 == 0) ? c[best].val + inset : c[best].val - inset;
            snaps.push_back({axis, best % 2, *c[best].bound, tf.rad, tf.idx});
        }
        std::fprintf(stderr, "[mesh3d] port snap: air box pulled onto the terminal planes -> "
                     "x[%.2f,%.2f] y[%.2f,%.2f] z[%.2f,%.2f] mm (leads now END on the boundary)\n",
                     bx0*1e3, bx1*1e3, by0*1e3, by1*1e3, bz0*1e3, bz1*1e3);
        // CLEARANCE GUARD (2026-09-12, with mvb-ab, ABT #1159): a snapped face may only cut the
        // leads it was snapped for. Any other conductor solid that stops short of the plane must
        // stay at least 0.9 cap radii away from it (MVB++'s drop rule guarantees >= 1 radius);
        // closer than that the sizing field between wrap copper and the box face collapses and
        // gmsh's HXT kernel aborts (common_mode_choke at inset 1.5). Refuse loudly with a witness
        // instead of meshing a model that dies later or slices a wrap.
        for (const auto& sp : snaps) {
            for (const auto& im : imps) {
                if (im.region.rfind("winding", 0) != 0 && im.region.rfind("conductor", 0) != 0 && im.region.rfind("turn", 0) != 0) continue;
                double a,b,cc,d,e,f; gmsh::model::occ::getBoundingBox(3, im.tag, a,b,cc,d,e,f);
                const double lo = sp.axis == 0 ? a : sp.axis == 1 ? b : cc, hi = sp.axis == 0 ? d : sp.axis == 1 ? e : f;
                if (lo <= sp.plane && hi >= sp.plane) continue;            // crosses the plane: a lead being clipped
                const double gap = sp.side == 0 ? lo - sp.plane : sp.plane - hi;   // inside distance to the face
                if (gap < 0.9 * sp.rad)
                    throw std::runtime_error("mesh3d_from_mas: port face for term" + std::to_string(sp.idx) + " (" +
                        std::string(sp.axis == 0 ? "x" : sp.axis == 1 ? "y" : "z") + (sp.side ? "+" : "-") + " = " +
                        std::to_string(sp.plane*1e3) + " mm) passes " + std::to_string(gap*1e3) + " mm from conductor solid " +
                        std::to_string(im.tag) + " ('" + im.region + "'), less than 0.9 cap radii (" +
                        std::to_string(0.9*sp.rad*1e3) + " mm) -- the lead drop is too short for this cap");
            }
        }
    }
    if (symPlanes > 0) std::fprintf(stderr, "[mesh3d] symmetry box: cutX=%d cutZ=%d cutY=%d (box x[%.2f,%.2f] y[%.2f,%.2f] z[%.2f,%.2f] mm)\n",
                                    (int)cutX, (int)cutZ, (int)cutY, bx0*1e3, bx1*1e3, by0*1e3, by1*1e3, bz0*1e3, bz1*1e3);
    const int box_tag = gmsh::model::occ::addBox(bx0, by0, bz0, bx1-bx0, by1-by0, bz1-bz0);
    gmsh::model::occ::synchronize();

    // Fragment the box against every solid -> one conforming model.
    // TOLERANCE LADDER: BOPAlgo can reject a fragment outright ("Fragments failed -
    // BOPAlgo_AlertBuilderFailed") on assemblies with many coincident faces -- the toroidal
    // conductor is a chain of mitre-jointed solids ABUTTING on shared elliptical faces, exactly
    // that case (seen on 05_pfc_inductor_t4020 and 12_boost_inductor_t5026). A small fuzzy value
    // lets the boolean absorb the sub-micron seams; retry with it rather than failing the whole
    // pipeline. Loud on retry, throws with the ladder when all tolerances fail.
    // CLIP THE PROTRUDING LEADS TO THE BOX (2026-09-08). After the port snap the box face cuts
    // through each lead, but the fragment runs NON-DESTRUCTIVE + GLUE (needed for the abutting
    // conductor chain), so it does not slice the stub off: 03_buck kept a lead reaching z = -23.22
    // while the box started at -20.94, i.e. the conduction path pierced the boundary instead of
    // ENDING on it, and the solve stayed ill-posed. Intersect each solid that pokes out with the
    // box first, so every lead ends in a clean planar cross-section lying in the boundary face.
    if (!termFaces.empty() && !std::getenv("OMFEM_NO_PORT_SNAP")) {
        int nclip = 0; std::vector<int> clippedTags;
        std::vector<std::pair<int,std::string>> extraSolids;   // (tag, region) from multi-piece clips
        for (auto& im : imps) {
            double a,b,c,d,e,f;
            try { gmsh::model::occ::getBoundingBox(3, im.tag, a,b,c,d,e,f); } catch (...) { continue; }
            const double eps = 1e-9;
            if (a >= bx0-eps && d <= bx1+eps && b >= by0-eps && e <= by1+eps &&
                c >= bz0-eps && f <= bz1+eps) continue;                   // wholly inside
            gmsh::vectorpair res; std::vector<gmsh::vectorpair> resMap;
            try {
                gmsh::model::occ::intersect({{3, im.tag}}, {{3, box_tag}}, res, resMap,
                                            -1, /*removeObject=*/true, /*removeTool=*/false);
            } catch (const std::exception& ex) {
                throw std::runtime_error(std::string("mesh3d_from_mas: clipping a lead to the air box "
                    "failed (the conduction path must END on the boundary): ") + ex.what());
            }
            if (res.empty()) throw std::runtime_error("mesh3d_from_mas: clipping a lead to the air box "
                                                      "removed it entirely");
            // An intersection can return SEVERAL solids (a lead that leaves and re-enters the box, a
            // winding cut in more than one place). Keeping only the first silently threw the rest of
            // the conductor away -- and with it the other port face (ETD34: "no port face found").
            im.tag = res.front().second;
            for (auto& r : res) clippedTags.push_back(r.second);
            for (size_t j = 1; j < res.size(); ++j) extraSolids.emplace_back(res[j].second, im.region);
            ++nclip;
        }
        // Extra pieces from a multi-solid clip must stay in the region list, or they lose their
        // material and silently become air.
        for (auto& e : extraSolids) imps.push_back({e.first, e.second, 0.0});
        if (nclip) { gmsh::model::occ::synchronize();
            std::fprintf(stderr, "[mesh3d] port clip: %d solid(s) trimmed to the air box so each lead "
                         "ENDS in a planar cross-section on the boundary\n", nclip);
            // The caps MOVED to the clip plane, so the post-mesh port matcher (which works from the
            // collected centres) must follow them: replace each terminal centre by the centre of the
            // clipped solid's face that lies IN that box plane and is nearest the old cap.
            const double bnd[6] = {bx0, bx1, by0, by1, bz0, bz1};
            // ONLY the clipped conductors' own faces: the big box face is also flat in this plane
            // and would win the nearest-centre test (it did: both ports collapsed onto (0,0,plane)).
            gmsh::vectorpair allf;
            for (int ct : clippedTags) {
                gmsh::vectorpair bd; gmsh::model::getBoundary({{3, ct}}, bd, false, false, false);
                for (auto& dt : bd) allf.push_back({2, std::abs(dt.second)});
            }
            for (size_t t = 0; t < termFaces.size() && t < snapAxis.size(); ++t) {
                const int ax = snapAxis[t] / 2, side = snapAxis[t];
                const double plane = bnd[side];
                double best = 1e30; double nx=0, ny=0, nz=0; bool found=false;
                double bestPlaneDev = 1e30;                 // how flat the flattest candidate is
                for (auto& dt : allf) {
                    double a,b,c2,d,e,f; gmsh::model::occ::getBoundingBox(2, dt.second, a,b,c2,d,e,f);
                    const double lo[3]={a,b,c2}, hi[3]={d,e,f};
                    bestPlaneDev = std::min(bestPlaneDev,
                        std::max(std::fabs(lo[ax]-plane), std::fabs(hi[ax]-plane)));
                    // GEOMETRY-RELATIVE flatness, not an absolute epsilon: the boolean leaves the cut
                    // face a few microns off the nominal plane (ETD34: 5 um), which a 1 um test threw
                    // away. A quarter of the cap's own radius is far below any other face's distance
                    // from this plane, so it cannot match the wrong face.
                    const double planeTol = 0.25 * termFaces[t].rad;
                    if (std::fabs(lo[ax]-plane) > planeTol || std::fabs(hi[ax]-plane) > planeTol) continue;
                    const double fx=0.5*(a+d), fy=0.5*(b+e), fz=0.5*(c2+f);
                    const double dd = std::fabs(fx-termFaces[t].cx)+std::fabs(fy-termFaces[t].cy)
                                    + std::fabs(fz-termFaces[t].cz);
                    if (dd < best) { best = dd; nx=fx; ny=fy; nz=fz; found=true; }
                }
                if (found) { termFaces[t].cx=nx; termFaces[t].cy=ny; termFaces[t].cz=nz; }
                else {
                    char msg[512];
                    std::snprintf(msg, sizeof msg,
                        "mesh3d_from_mas: lead %zu was clipped to the air box but no port face lies in its "
                        "clip plane (axis %d, side %d, plane %.4f mm; cap was at %.3f,%.3f,%.3f mm; "
                        "%zu candidate face(s) on %zu clipped solid(s); flattest candidate is %.4g mm "
                        "off the plane)",
                        t, ax, side, plane*1e3, termFaces[t].cx*1e3, termFaces[t].cy*1e3,
                        termFaces[t].cz*1e3, allf.size(), clippedTags.size(), bestPlaneDev*1e3);
                    throw std::runtime_error(msg);
                }
            }
        }
    }
    std::vector<gmsh::vectorpair> outMap;
    {
        gmsh::vectorpair obj{{3, box_tag}}, tools;
        for (auto& im : imps) tools.push_back({3, im.tag});
        gmsh::vectorpair outdt; std::vector<gmsh::vectorpair> outMap0;
        // OMFEM_FRAG_TOL starts the ladder at a given fuzzy value: a butt-joined conductor
        // CHAIN can fragment "successfully" at exact tolerance yet leave a face with a broken
        // edge loop that only surfaces at 2D meshing ("the 1D mesh seems not to be forming a
        // closed loop") -- a fuzzy fragment glues those coincident discs properly.
        const double fragTol0 = std::getenv("OMFEM_FRAG_TOL")
                              ? std::atof(std::getenv("OMFEM_FRAG_TOL")) : 0.0;
        // gmsh >= 4.13: BOPAlgo GLUE (shift mode) + NON-DESTRUCTIVE for the fragment. The
        // conformal winding abuts on EXACTLY coincident faces with no true intersections --
        // precisely the input class SetGlue is specified for (2026-08 research: OCCT boolean
        // spec + gmsh 4.13 CHANGELOG); non-destructive stops tolerance creep on the inputs.
        // OMFEM_NO_GLUE=1 restores the old behaviour for comparison.
        if (!std::getenv("OMFEM_NO_GLUE")) {
            gmsh::option::setNumber("Geometry.OCCBooleanGlue", 1);
            gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", 1);
        }
        const double tolLadder[] = {fragTol0, 1e-9, 1e-7};
        std::string ferrs;
        bool done = false;
        for (double tol : tolLadder) {
            try {
                if (tol > 0.0) {
                    std::fprintf(stderr, "[mesh3d] RETRY fragment with boolean tolerance %g\n", tol);
                    gmsh::option::setNumber("Geometry.ToleranceBoolean", tol);
                }
                gmsh::vectorpair o = obj, t = tools;
                gmsh::model::occ::fragment(o, t, outdt, outMap0);
                done = true;
                break;
            } catch (const std::exception& e) {
                char tolbuf[32]; std::snprintf(tolbuf, sizeof tolbuf, "%g", tol);
                ferrs += (ferrs.empty() ? "" : "; ") + std::string("tol=") + tolbuf +
                         ": " + e.what();
            }
        }
        gmsh::option::setNumber("Geometry.ToleranceBoolean", 0.0);   // never leak into meshing
        gmsh::option::setNumber("Geometry.OCCBooleanGlue", 0);
        gmsh::option::setNumber("Geometry.OCCBooleanNonDestructive", 0);
        if (!done)
            throw std::runtime_error("mesh3d_from_mas: fragment failed at every boolean "
                                     "tolerance [" + ferrs + "]");
        outMap = outMap0;
    }
    gmsh::model::occ::synchronize();

    // TRIM TO THE BOX (2026-09-08). The fragment keeps every piece, including the part of a lead
    // that now lies OUTSIDE the air box after the port snap. Leaving it there means the conduction
    // path still does not terminate on the boundary -- it just pokes through it -- and the system
    // stays ill-posed (03_buck: 1e8 V, while the flush straight bar solves to 5e-13). Delete the
    // volumes whose centre of mass is outside the box, so each lead ENDS in a clean planar
    // cross-section lying in the boundary face, which is what the port needs.
    if (!termFaces.empty() && !std::getenv("OMFEM_NO_PORT_SNAP")) {
        gmsh::vectorpair allv; gmsh::model::occ::getEntities(allv, 3);
        gmsh::vectorpair outside;
        const double tol = 1e-9;
        for (auto& dt : allv) {
            double cmx, cmy, cmz;
            try { gmsh::model::occ::getCenterOfMass(3, dt.second, cmx, cmy, cmz); } catch (...) { continue; }
            if (cmx < bx0-tol || cmx > bx1+tol || cmy < by0-tol || cmy > by1+tol ||
                cmz < bz0-tol || cmz > bz1+tol) outside.push_back(dt);
        }
        std::fprintf(stderr, "[mesh3d] port trim: %zu volume(s) after fragment, %zu outside the box "
                     "z[%.3f,%.3f] mm\n", allv.size(), outside.size(), bz0*1e3, bz1*1e3);
        if (!outside.empty()) {
            gmsh::model::occ::remove(outside, true);
            gmsh::model::occ::synchronize();
            std::fprintf(stderr, "[mesh3d] port trim: removed %zu volume(s) (lead stubs beyond the port plane)\n",
                         outside.size());
        }
    }

    // Assign region per result volume via the fragment map (outMap[0]=box/air,
    // outMap[i+1]=children of tool i). Descending volume so specifics win on overlap.
    std::map<int,std::string> vol_region;
    std::vector<size_t> order(imps.size());
    for (size_t i=0;i<imps.size();++i) order[i]=i;
    std::sort(order.begin(), order.end(), [&](size_t a,size_t b){ return imps[a].vol>imps[b].vol; });
    for (size_t k=0;k<order.size();++k) {
        size_t i = order[k];
        for (auto& dt : outMap[i+1]) if (dt.first==3) vol_region[dt.second]=imps[i].region;
    }

    std::map<std::string,std::vector<int>> groups;
    gmsh::vectorpair allv; gmsh::model::getEntities(allv, 3);
    for (auto& dt : allv) {
        auto it = vol_region.find(dt.second);
        groups[it==vol_region.end() ? "air" : it->second].push_back(dt.second);
    }
    int pg = 1;
    for (auto& [name, vs] : groups) {
        if (vs.empty()) continue;
        gmsh::model::addPhysicalGroup(3, vs, pg); gmsh::model::setPhysicalName(3, pg, name); ++pg;
        std::fprintf(stderr, "[mesh3d] group %-8s %zu vols\n", name.c_str(), vs.size());
    }

    // OMFEM_MESH_DEBUG: report DUPLICATE faces left by the fragment. Two coincident faces that
    // were not unified give the 3D stage identical triangles from two different surfaces
    // ("Invalid boundary mesh (overlapping facets) on surface A surface B") and no mesh at all;
    // the error names only tags, so locate and size them here.
    if (std::getenv("OMFEM_MESH_DEBUG")) {
        struct FInfo { int tag; double cx, cy, cz, area; };
        std::vector<FInfo> fi;
        gmsh::vectorpair fs; gmsh::model::getEntities(fs, 2);
        for (auto& dt : fs) {
            double a,b,c,d,e2,f2; gmsh::model::occ::getBoundingBox(2,dt.second,a,b,c,d,e2,f2);
            double m = 0.0; gmsh::model::occ::getMass(2, dt.second, m);
            fi.push_back({dt.second, 0.5*(a+d), 0.5*(b+e2), 0.5*(c+f2), m});
        }
        int ndup = 0;
        for (size_t i = 0; i < fi.size(); ++i)
            for (size_t j = i+1; j < fi.size(); ++j) {
                const double tol = 1e-9;
                if (std::fabs(fi[i].cx-fi[j].cx) > tol || std::fabs(fi[i].cy-fi[j].cy) > tol ||
                    std::fabs(fi[i].cz-fi[j].cz) > tol) continue;
                if (std::fabs(fi[i].area-fi[j].area) > tol*std::max(fi[i].area, fi[j].area) + 1e-15)
                    continue;
                std::fprintf(stderr, "[mesh3d.dbg] DUPLICATE faces %d and %d: centre "
                             "(%.4f,%.4f,%.4f)mm area=%.4gmm2\n", fi[i].tag, fi[j].tag,
                             fi[i].cx*1e3, fi[i].cy*1e3, fi[i].cz*1e3, fi[i].area*1e6);
                ++ndup;
            }
        std::fprintf(stderr, "[mesh3d.dbg] %zu faces, %d duplicate pair(s)\n", fi.size(), ndup);
        if (std::atoi(std::getenv("OMFEM_MESH_DEBUG")) >= 2)
            for (const auto& f : fi) {
                double a,b,c,d,e2,f2; gmsh::model::occ::getBoundingBox(2,f.tag,a,b,c,d,e2,f2);
                std::fprintf(stderr, "[mesh3d.dbg] face %d area=%.4gmm2 bbox x[%.3f,%.3f] "
                             "y[%.3f,%.3f] z[%.3f,%.3f] mm\n", f.tag, f.area*1e6,
                             a*1e3,d*1e3, b*1e3,e2*1e3, c*1e3,f2*1e3);
            }
    }

    // 'outer' = the box faces (far-field A x n = 0), INCLUDING the x = 0 / z = 0 symmetry cuts
    // (azimuthal planes, B tangential -> a x n = 0) and EXCLUDING an axial mid-plane cut at y = 0
    // (B normal -> natural n x H = 0).
    std::vector<int> outer;
    const double e = 1e-4 * std::max({2*hx, 2*hy, 2*hz});
    // A terminal cap that now LIES IN a box face (port snap, above) must NOT join 'outer': the
    // far-field condition there is a x n = 0, which forbids exactly the current we inject through
    // it. Maxwell's rule is the same -- a conductor touching the boundary uses its cross-section as
    // the excitation surface, not the field boundary condition. Skip any face centred on a cap.
    const double ptol = 0.5 * std::max({2*hx, 2*hy, 2*hz}) * 5e-3 + 1e-4;
    // CAP-SIZED as well as cap-centred (00_debug, 2026-09-23): a one-turn winding at x = y = 0 puts its lead
    // cap at the centre of the z- box plane, and the box face -- an annulus around that cap, symmetric in x
    // and y -- has its bounding-box centre on the cap too. Centre-only matching skipped it as a terminal and
    // the plane had no outer face. A cap spans ~2 rad in the plane; the box face spans the box.
    auto isTerminalFace = [&](double x0, double y0, double z0, double x1, double y1, double z1) {
        const double fx = 0.5*(x0+x1), fy = 0.5*(y0+y1), fz = 0.5*(z0+z1);
        const double span = std::max({x1-x0, y1-y0, z1-z0});
        for (const auto& tf : termFaces)
            if (std::fabs(fx-tf.cx) + std::fabs(fy-tf.cy) + std::fabs(fz-tf.cz) < ptol && span <= 2.5 * tf.rad + ptol) return true;
        return false;
    };
    // ON-PLANE TEST (2026-09-12): OCC inflates a face's bounding box by its shape tolerance, and
    // the box face that was CUT by the port caps comes out of the fragment with a tolerance of
    // several um (15_gan: z[-4.572,-4.562] mm for the plane z = -4.5675). Testing both bbox ends
    // against the plane with e (1.9 um there) silently dropped that face from 'outer': the cap
    // rims lost their a x n = 0 edges and the A-V system stalled with the current rows flat at
    // 0.9999 (buck_inductor_complete, 15_gan, ETD34 on the CPU lane). A planar face lies on the
    // plane when its bbox is THIN (its mid-plane within e + half its thickness of the bound);
    // the thickness cap is scaled to the box so no solid face can pass.
    const double thinCap = 1e-2 * std::max({2*hx, 2*hy, 2*hz});
    auto onPlane = [&](double lo_, double hi_, double bound) {
        const double thick = hi_ - lo_;
        return thick < thinCap && std::fabs(0.5*(lo_+hi_) - bound) < e + 0.5*thick;
    };
    int planeHits[6] = {0,0,0,0,0,0};
    gmsh::vectorpair faces; gmsh::model::getEntities(faces, 2);
    for (auto& dt : faces) {
        double a,b,c,d,ee,f; gmsh::model::occ::getBoundingBox(2,dt.second,a,b,c,d,ee,f);
        if (isTerminalFace(a, b, c, d, ee, f)) continue;
        int hit = -1;
        if (onPlane(a,d,bx0)) hit = 0; else if (onPlane(a,d,bx1)) hit = 1;
        else if (!cutY && onPlane(b,ee,by0)) hit = 2; else if (onPlane(b,ee,by1)) hit = 3;
        else if (onPlane(c,f,bz0)) hit = 4; else if (onPlane(c,f,bz1)) hit = 5;
        if (hit >= 0) { outer.push_back(dt.second); planeHits[hit]++; }
    }
    // Every box plane must carry at least one outer face (the y- plane is the natural-BC cut when
    // cutY). A plane without one means the far-field condition is missing on that side -- refuse.
    {
        static const char* pn[6] = {"x-","x+","y-","y+","z-","z+"};
        std::string missing;
        for (int i = 0; i < 6; i++) if (!(i == 2 && cutY) && planeHits[i] == 0) missing += std::string(missing.empty() ? "" : ",") + pn[i];
        std::fprintf(stderr, "[mesh3d] outer: %zu faces on the box planes (x- %d x+ %d y- %d y+ %d z- %d z+ %d)\n",
                     outer.size(), planeHits[0], planeHits[1], planeHits[2], planeHits[3], planeHits[4], planeHits[5]);
        if (!missing.empty())
            throw std::runtime_error("mesh3d: air-box plane(s) " + missing + " have no 'outer' face -- the far-field "
                                     "boundary is incomplete (a fragmented box face failed the on-plane test?)");
    }
    if (!outer.empty()) { gmsh::model::addPhysicalGroup(2, outer, pg); gmsh::model::setPhysicalName(2, pg, "outer"); ++pg; }

    // ---- PORT TERMINAL faces: tag the meshed winding-boundary face nearest each collected terminal
    // cap centre as 'term<2*port+k>' (the A-V solver's current-injection surfaces). Match by centre-of-
    // mass; a terminal cap is a small planar face whose centroid is far from the winding's lateral
    // surface, so the nearest-centre match is unambiguous within a wire-radius tolerance.
    if (!termFaces.empty()) {
        const double ttol = 0.5 * std::max({2*hx, 2*hy, 2*hz}) * 5e-3 + 1e-4;  // generous, feature-relative
        std::map<int, std::vector<int>> termGroups;
        gmsh::vectorpair f2; gmsh::model::getEntities(f2, 2);
        for (auto& dt : f2) {
            double a,b,c,d,ee,f; gmsh::model::occ::getBoundingBox(2,dt.second,a,b,c,d,ee,f);
            const double fx=0.5*(a+d), fy=0.5*(b+ee), fz=0.5*(c+f);
            double best=1e30; int bi=-1;
            for (const auto& tf : termFaces) {
                const double dd=std::fabs(fx-tf.cx)+std::fabs(fy-tf.cy)+std::fabs(fz-tf.cz);
                if (dd<best){ best=dd; bi=tf.idx; }
            }
            if (bi>=0 && best<ttol) termGroups[bi].push_back(dt.second);
        }
        // How many CAPS feed each port, so "2 face(s)" on a 3-parallel winding is visibly wrong.
        std::map<int,int> capsPerIdx;
        for (const auto& tf : termFaces) capsPerIdx[tf.idx]++;
        for (auto& [idx, fs] : termGroups) {
            gmsh::model::addPhysicalGroup(2, fs, pg);
            gmsh::model::setPhysicalName(2, pg, "term" + std::to_string(idx)); ++pg;
            std::fprintf(stderr, "[mesh3d] port term%d: %zu face(s) from %d cap(s)\n",
                         idx, fs.size(), capsPerIdx[idx]);
        }
        // PER-CAP check, not per-port. All parallels of a winding share one port index, so a port
        // whose OTHER parallels matched looks healthy while one strand is left with no terminal at
        // all -- electrically floating, Gd = 0, |C| = 0, and the SERIES solve aborts on a
        // "mesh-degenerate turn" hours later with no hint of the cause (buck_inductor_complete:
        // 2 of 3 strands, ABT #1155). Report every cap that found no face of its own.
        std::map<int,int> faceOwner;   // meshed face -> how many caps call it their nearest
        for (const auto& tf : termFaces) {
            double best = 1e30; int bf = -1;
            for (auto& dt : f2) {
                double a,b,c,d,ee,f; gmsh::model::occ::getBoundingBox(2,dt.second,a,b,c,d,ee,f);
                const double dd = std::fabs(0.5*(a+d)-tf.cx) + std::fabs(0.5*(b+ee)-tf.cy)
                                + std::fabs(0.5*(c+f)-tf.cz);
                if (dd < best) { best = dd; bf = dt.second; }
            }
            std::fprintf(stderr, "[mesh3d]   cap term%d at (%.4f, %.4f, %.4f) mm -> face %d at "
                         "%.4g mm%s\n", tf.idx, tf.cx*1e3, tf.cy*1e3, tf.cz*1e3, bf, best*1e3,
                         best >= ttol ? "  *** BEYOND TOL, no terminal for this parallel ***" : "");
            if (best < ttol) faceOwner[bf]++;
        }
        for (const auto& [fid, n] : faceOwner)
            if (n > 1)
                std::fprintf(stderr, "[mesh3d] WARN meshed face %d is the nearest face for %d caps -- "
                             "those parallels collapse onto ONE terminal, so all but one are left "
                             "with no conducting path (Gd = 0)\n", fid, n);
        // A port with NO matched face is fatal for the eddy/full-wave solvers -- say WHERE the
        // cap was expected and how close the nearest meshed face came, or the failure surfaces
        // hours later as an inscrutable femready NOT_READY.
        for (const auto& tf : termFaces)
            if (termGroups.find(tf.idx) == termGroups.end()) {
                std::vector<std::pair<double,int>> near;
                for (auto& dt : f2) {
                    double a,b,c,d,ee,f; gmsh::model::occ::getBoundingBox(2,dt.second,a,b,c,d,ee,f);
                    near.push_back({std::fabs(0.5*(a+d)-tf.cx)+std::fabs(0.5*(b+ee)-tf.cy)
                                    +std::fabs(0.5*(c+f)-tf.cz), dt.second});
                }
                std::sort(near.begin(), near.end());
                double best = near.empty() ? 1e30 : near[0].first;
                for (size_t q = 0; q < near.size() && q < 3; ++q) {
                    double a,b,c,d,ee,f; gmsh::model::occ::getBoundingBox(2,near[q].second,a,b,c,d,ee,f);
                    std::fprintf(stderr, "[mesh3d]   near-face %d: centre (%.3f, %.3f, %.3f) mm "
                                 "bbox y[%.3f,%.3f] z[%.3f,%.3f]\n", near[q].second,
                                 0.5*(a+d)*1e3, 0.5*(b+ee)*1e3, 0.5*(c+f)*1e3,
                                 b*1e3, ee*1e3, c*1e3, f*1e3);
                }
                std::fprintf(stderr, "[mesh3d] WARN port term%d UNMATCHED: cap centre "
                             "(%.4f, %.4f, %.4f) mm, nearest meshed face %.4g mm away "
                             "(tol %.4g mm)\n", tf.idx, tf.cx*1e3, tf.cy*1e3, tf.cz*1e3,
                             best*1e3, ttol*1e3);
            }
    }

    // ---- ANISOTROPIC winding mesh (OMFEM_ANISO): fine ACROSS each wire, coarse ALONG it.
    // The wire runs azimuthally about the column axis (MVB++ revolves round-column turns
    // about Y), so the metric is the SAME analytic frame everywhere: fine in radial+axial,
    // coarse in azimuthal. We build it as a MathEvalAniso metric (1/h^2 convention),
    // Restrict it to the winding volumes (coarse isotropic elsewhere), and mesh with BAMG
    // (2D anisotropic surface) + MMG3D (3D anisotropic) -- the only combo gmsh honours.
    const auto t_mesh0 = std::chrono::steady_clock::now();
    const bool aniso = std::getenv("OMFEM_ANISO");
    if (aniso) {
        const char axis = std::getenv("OMFEM_AXIS") ? std::getenv("OMFEM_AXIS")[0] : 'Y';
        const double h_across = std::getenv("OMFEM_HACROSS") ? std::atof(std::getenv("OMFEM_HACROSS"))
                                : std::max(opt.conductor_target, 2.0e-4);
        const double h_along  = std::getenv("OMFEM_HALONG") ? std::atof(std::getenv("OMFEM_HALONG")) : 2.0e-3;
        const double a = 1.0/(h_across*h_across);
        const double b = 1.0/(h_along*h_along) - a;     // < 0: coarsen the azimuthal direction
        auto F = [](double v){ char s[64]; std::snprintf(s, sizeof s, "%.10g", v); return std::string(s); };
        // metric = a*I + b*(e_phi (x) e_phi); e_phi azimuthal about `axis`.
        std::string r2, m11, m22, m33, m12="0", m13="0", m23="0";
        if (axis=='Y')      { r2="(x^2+z^2+1e-30)"; m11=F(a)+"+("+F(b)+")*z^2/"+r2; m22=F(a); m33=F(a)+"+("+F(b)+")*x^2/"+r2; m13="("+F(b)+")*(-x*z)/"+r2; }
        else if (axis=='Z') { r2="(x^2+y^2+1e-30)"; m11=F(a)+"+("+F(b)+")*y^2/"+r2; m22=F(a)+"+("+F(b)+")*x^2/"+r2; m33=F(a); m12="("+F(b)+")*(-x*y)/"+r2; }
        else                { r2="(y^2+z^2+1e-30)"; m11=F(a); m22=F(a)+"+("+F(b)+")*z^2/"+r2; m33=F(a)+"+("+F(b)+")*y^2/"+r2; m23="("+F(b)+")*(-y*z)/"+r2; }
        const int fa = gmsh::model::mesh::field::add("MathEvalAniso");
        gmsh::model::mesh::field::setString(fa,"m11",m11); gmsh::model::mesh::field::setString(fa,"m22",m22);
        gmsh::model::mesh::field::setString(fa,"m33",m33); gmsh::model::mesh::field::setString(fa,"m12",m12);
        gmsh::model::mesh::field::setString(fa,"m13",m13); gmsh::model::mesh::field::setString(fa,"m23",m23);
        // Region names: "winding" (level-set), "winding_<name>" (REAL/continuous conductor,
        // OMFEM_CONTINUOUS) and "turn_<w>_<i>" (per-turn solids). An exact lookup of "winding"
        // matches only the first, so on the continuous and per-turn paths the Restrict field got an
        // EMPTY volume list: the anisotropic metric then applied NOWHERE, leaving BAMG to mesh
        // curved conductor CAD under the bare air target -- which is why it failed outright
        // ("BAMG failed") instead of producing the ~10-100x cheaper winding mesh it exists for.
        // Same defect, and the same fix, as the isotropic branch below.
        std::vector<double> wv;
        for (const auto& [gname, vs] : groups) {
            if (gname.rfind("winding", 0) != 0 && gname.rfind("turn_", 0) != 0) continue;
            for (int v : vs) wv.push_back(v);
        }
        if (wv.empty())
            throw std::runtime_error("mesh3d_from_mas(ANISO): no winding/turn volume groups to "
                                     "restrict the anisotropic metric to -- it would apply nowhere");
        const int fr = gmsh::model::mesh::field::add("Restrict");
        gmsh::model::mesh::field::setNumber(fr,"InField",fa);
        gmsh::model::mesh::field::setNumbers(fr,"VolumesList",wv);
        const int fc = gmsh::model::mesh::field::add("MathEval");
        gmsh::model::mesh::field::setString(fc,"F",F(opt.air_target));
        const int fm = gmsh::model::mesh::field::add("Min");
        gmsh::model::mesh::field::setNumbers(fm,"FieldsList",{(double)fr,(double)fc});
        gmsh::model::mesh::field::setAsBackgroundMesh(fm);
        gmsh::option::setNumber("Mesh.MeshSizeFromPoints",0);
        gmsh::option::setNumber("Mesh.MeshSizeFromCurvature",0);
        gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary",0);
        // 2D surface algorithm: BAMG(7) is the anisotropic one but fragile on curved CAD;
        // MeshAdapt(1)/Frontal(6) are robust isotropic and let MMG3D do the 3D anisotropy.
        int algo2d = std::getenv("OMFEM_ALGO2D") ? std::atoi(std::getenv("OMFEM_ALGO2D")) : 7;
        gmsh::option::setNumber("Mesh.Algorithm", algo2d);
        gmsh::option::setNumber("Mesh.Algorithm3D",7);   // MMG3D: 3D anisotropic
        std::fprintf(stderr,"[mesh3d] ANISO axis=%c across=%.3gmm along=%.3gmm winding_vols=%zu (BAMG+MMG3D)\n",
                     axis, h_across*1e3, h_along*1e3, wv.size());
        gmsh::model::mesh::generate(3);
    } else {
        // Adaptive sizing (shared SizeField): conductor skin refinement + air ceiling composed into
        // one Min field. Then a ROBUST algorithm chain (Delaunay=1 first, HXT=10 fallback) -- the
        // old default Mesh.Algorithm3D=7 (MMG3D) is an adaptive backend that ignores general size
        // fields and frequently fails the first-pass write; it is never the generator here.
        // ---- RETRY LADDER -------------------------------------------------------------------
        // A dense real winding is bounded on BOTH sides: too coarse and the surface triangles of
        // nearly-touching curved surfaces CROSS (adjacent turns sit one coating apart -- ~0.04-0.09
        // mm on ETD34 -- so gmsh dies with "PLC Error: a segment and a facet intersect",
        // "Invalid boundary mesh (overlapping facets)" or "Unable to recover the edge"); too fine
        // and the mesh is unaffordable. The viable window is design-dependent (chord sagitta
        // h^2/8r vs the clearance), so no fixed target suits every design. When the failure is one
        // of those clearance signatures, retry with the conductor target scaled by 0.65 -- LOUDLY,
        // and rethrowing after the ladder is exhausted. This is adaptive search, not a silent
        // fallback: every retry is printed and the final failure carries the whole ladder.
        const int max_attempts = 1 + (std::getenv("OMFEM_MESH_RETRIES")
                                      ? std::atoi(std::getenv("OMFEM_MESH_RETRIES")) : 3);
        // conductor_target contract: > 0 = use as-is; 0 (default) = derive from the layout's
        // real clearances (compute_auto_cond_target); < 0 = NO conductor refinement at all.
        // CONTACT FLOOR: is there an emitted solid the wire RESTS on? When a bobbin wall, an
        // insulation layer or the coating shell is in the model, the bare-copper surface sits
        // one coating thickness (ro - rc) from it everywhere the winding touches it -- a real
        // corridor between surfaces gmsh must triangulate, but NOT a turn pair, so the layout
        // cross-section cannot place it and it stays a global bound. With none of them emitted
        // (the corpus default: groups are air / core / winding_<name> only) that corridor exists
        // in no emitted surface at all and must not bound anything.
        bool contact_floor = false;
        for (const auto& [gname, gvs] : groups) {
            (void)gvs;
            if (gname.rfind("bobbin", 0) == 0 || gname.rfind("insulation", 0) == 0 ||
                gname.rfind("coating", 0) == 0) { contact_floor = true; break; }
        }
        const double cond_seed = opt.conductor_target > 0.0 ? opt.conductor_target
                               : opt.conductor_target < 0.0 ? 0.0
                                                            : compute_auto_cond_target(contact_floor);
        double cond_t = cond_seed;
        double ladder_seed = 0.0;   // fixed reference for the retry ladder (set on first failure)
        std::string ladder;
        for (int attempt = 0; ; ++attempt) {
            try {
                if (attempt > 0) gmsh::model::mesh::clear();
                gmsh::model::mesh::removeSizeCallback();   // an earlier attempt's corridor sizes do not carry over
                mvb::mesh::SizeFieldBuilder sizer(3, opt.air_target);
                sizer.add_ceiling(opt.air_target);
                double cmin = opt.core_target;
                if (cond_t > 0.0) {
                    std::vector<double> cond_faces;
                    // PER-WINDING TARGETS (OMFEM_COND_TARGET_WINDINGS "winding_<name>=<m>;..."): each
                    // winding's copper refined at its own target, from its own wire. The design-wide
                    // cond_t is the finest of them; a retry that lowers cond_t scales all of them by the
                    // same factor, so the ladder keeps its meaning. A named winding group missing from the
                    // map is an error, not a default.
                    std::map<std::string, double> perW;
                    if (const char* pw = std::getenv("OMFEM_COND_TARGET_WINDINGS"); pw && *pw) {
                        std::stringstream ss(pw); std::string item;
                        while (std::getline(ss, item, ';')) {
                            const auto eq = item.rfind('=');
                            if (eq == std::string::npos || eq == 0)
                                throw std::runtime_error("OMFEM_COND_TARGET_WINDINGS: bad entry '" + item + "'");
                            const double t = std::stod(item.substr(eq + 1));
                            if (!(t > 0.0)) throw std::runtime_error("OMFEM_COND_TARGET_WINDINGS: non-positive target in '" + item + "'");
                            perW[item.substr(0, eq)] = t;
                        }
                    }
                    double perWmin = 0.0;
                    for (const auto& [g, t] : perW) perWmin = perWmin > 0.0 ? std::min(perWmin, t) : t;
                    const double perWscale = perW.empty() ? 1.0 : cond_t / perWmin;
                    std::map<double, std::vector<double>> faces_by_target;   // per-winding targets only
                    for (auto& [name, vs] : groups) {
                        // Region names: "winding" (level-set), "turn_<w>_<i>" (per-turn solids) and
                        // "winding_<name>" (REAL/continuous conductor, OMFEM_CONTINUOUS). The last one
                        // matched NEITHER test, so the continuous path silently skipped conductor
                        // refinement entirely — e138 came out at 5k tets with an unresolved 0.4 mm wire,
                        // i.e. no skin/proximity physics at all. Match the "winding" prefix too.
                        if (name.rfind("winding", 0) != 0 && name.rfind("turn_", 0) != 0) continue;
                        std::vector<double>* into = &cond_faces;
                        if (!perW.empty() && name.rfind("winding_", 0) == 0) {
                            const auto it = perW.find(name);
                            if (it == perW.end())
                                throw std::runtime_error("mesh3d: conductor group '" + name + "' has no entry in "
                                                         "OMFEM_COND_TARGET_WINDINGS");
                            into = &faces_by_target[it->second * perWscale];
                        }
                        for (int v : vs) {
                            gmsh::vectorpair bd; gmsh::model::getBoundary({{3, v}}, bd, false, false, false);
                            for (auto& dt : bd) if (dt.first == 2) into->push_back(std::abs(dt.second));
                        }
                    }
                    // Halo: how far the conductor-target refinement reaches into the surrounding air,
                    // in multiples of conductor_target. 8x is generous — on a densely wound design the
                    // halos of neighbouring turns merge and fill the whole winding window at conductor
                    // resolution, which is what actually sets the DOF count (PQ50, 38 turns, 0.16 mm:
                    // coarsening air 4->8 mm and core 3->6 mm changed the mesh by only 10%, because the
                    // window is halo, not far field). Tighten it when the direct solve must fit in RAM;
                    // the wire surfaces stay resolved either way.
                    const double halo = std::getenv("OMFEM_COND_HALO")
                                      ? std::atof(std::getenv("OMFEM_COND_HALO")) : 8.0;
                    // DistMin must cover the conductor's own INTERIOR, not just a skin around it.
                    // It used to be conductor_target (0.08 mm) while a 0.4 mm wire has a 0.2 mm radius,
                    // so the wire centre sat in the growth band and got ~0.5 mm elements — coarser than
                    // the wire itself, i.e. no skin/proximity resolution at all (e138 continuous: 7.6k
                    // tets). Fine out to the full halo, then grow to the air ceiling beyond it.
                    // FINE-BAND RADIUS: a PHYSICAL distance, not a multiple of the element size.
                    // fine_radius = halo * cond_t couples the air refinement to a sizing DECISION:
                    // raise the conductor target and the fine band around every conductor widens with
                    // it. Measured on ETD34 (2026-09-20): the target moved 0.1594 -> 0.2 mm, the band
                    // 31.9 -> 40 um, and the base mesh's COPPER tets FELL 1.8% while air+core ROSE
                    // 15.6% and the copper-boundary triangle count rose 11.8% -- which then fed MMG
                    // ~6% more junction metric nodes. A coarser target producing a bigger mesh is the
                    // surprise this coupling creates.
                    // The comment above says what the band is FOR: covering the wire's own interior.
                    // That is a property of the WIRE, so OMFEM_COND_HALO_WIRE expresses it in wire
                    // radii and is independent of the target. OMFEM_COND_HALO (x target) still works
                    // and still wins when set, so nothing outside the corpus changes behaviour.
                    double fine_radius = halo * cond_t;
                    if (const char* hw = std::getenv("OMFEM_COND_HALO_WIRE")) {
                        const double rw = layout_sizing().rsmin;
                        if (!(rw > 0.0))
                            throw std::runtime_error("OMFEM_COND_HALO_WIRE: the layout has no positive "
                                                     "wire surface radius to scale the fine band by");
                        fine_radius = std::atof(hw) * rw;
                        std::fprintf(stderr, "[mesh3d] fine band from the WIRE: %.4g x r_wire %.4g mm "
                                     "= %.4g mm (target-independent)\n", std::atof(hw), rw*1e3,
                                     fine_radius*1e3);
                    }
                    std::fprintf(stderr, "[mesh3d] conductor refinement: %zu faces target=%.4gmm fine_r=%.4gmm\n",
                                 cond_faces.size(), cond_t*1e3, fine_radius*1e3);
                    if (!cond_faces.empty()) {
                        sizer.add_distance_refinement(cond_faces, cond_t, fine_radius,
                                                      2.0 * fine_radius, 100);
                        cmin = std::min(cmin, cond_t);
                    }
                    for (const auto& [t, faces] : faces_by_target) {
                        // the fine band is the WIRE's (OMFEM_COND_HALO_WIRE) when set; else it scales
                        // with this winding's own target, as the design-wide band does with cond_t
                        const double fr = std::getenv("OMFEM_COND_HALO_WIRE") ? fine_radius : halo * t;
                        std::fprintf(stderr, "[mesh3d]   per-winding refinement: %zu faces target=%.4gmm fine_r=%.4gmm\n",
                                     faces.size(), t * 1e3, fr * 1e3);
                        sizer.add_distance_refinement(faces, t, fr, 2.0 * fr, 100);
                        cmin = std::min(cmin, t);
                    }
                    // ---- LOCAL CHORD-BOUND REFINEMENT --------------------------------------
                    // Everything above is the design-wide conductor size. THIS is where the tight
                    // corridors -- and only they -- get the finer size their clearance demands.
                    //
                    // THE COORDINATE FRAME. turnsDescription coordinates are the 2D winding-window
                    // cross-section: x is the distance from the core's central axis, y the axial
                    // position. MVB++ builds a round-column winding by carrying each turn around
                    // that axis (the axis is Y in the model; the ANISO path above relies on the
                    // same frame), so a point of the cross-section is not a point of the model but
                    // a CIRCLE: radius x about the axis, at height y. A corridor of the layout is
                    // therefore a thin torus, which is exactly what SizeFieldBuilder::RingBand is.
                    // This mapping holds for a winding revolved about an axis -- the corpus's
                    // E/ETD/PQ round-column designs. It does NOT hold for a winding whose path is
                    // a stadium around a rectangular column (the layout x is then a distance from
                    // the column FACE, not from the axis) nor for a toroid (whose turns wrap the
                    // ring and whose layout is the ring plane, not a revolved cross-section). That
                    // is not asserted, it is VERIFIED below against the built geometry, and the
                    // mesher refuses rather than guess when the verification fails.
                    std::vector<mvb::mesh::SizeFieldBuilder::RingBand> bands;
                    if (!std::getenv("OMFEM_CHORD_GLOBAL")) {
                        const LayoutSizing& ls = layout_sizing();
                        std::vector<TurnCorridor> tight;
                        for (const auto& c : ls.corridors)
                            if (c.size < cond_t) tight.push_back(c);
                        if (!tight.empty() && continuous) {
                            // REAL WINDING: the corridors are located on the built 3D centrelines
                            // (chord_points_3d), not by revolving the layout -- MVB++'s turns are
                            // helices/stadiums, and the circle mapping below refused 17_cllc, both
                            // pushpulls and isolated_buck ("cannot locate layout").
                            if (!realPaths)
                                throw std::runtime_error("mesh3d: " + std::to_string(tight.size()) + " layout corridor(s) "
                                    "are tighter than the conductor target but the real-winding centrelines that "
                                    "place their refinement could not be built: " + realPathsError);
                            auto cp = std::make_shared<const ChordPoints>(chord_points_3d(*realPaths, cond_t));
                            if (cp->h.empty())
                                throw std::runtime_error("mesh3d: the layout has " + std::to_string(tight.size()) +
                                    " corridor(s) tighter than the conductor target (tightest " +
                                    std::to_string(tight.front().size * 1e3) + " mm) but the built centrelines "
                                    "have none: layout and geometry disagree");
                            // VERIFICATION -- a corridor point is AIR. One inside copper means the
                            // centrelines or radii disagree with the solids that will be meshed.
                            std::vector<int> cvols;
                            for (auto& [gname, gvs] : groups)
                                if (gname.rfind("winding", 0) == 0 || gname.rfind("turn_", 0) == 0)
                                    for (int v : gvs) cvols.push_back(v);
                            const size_t nck = std::min<size_t>(cp->h.size(), 256), stride = std::max<size_t>(1, cp->h.size() / nck);
                            for (size_t i = 0; i < cp->h.size(); i += stride)
                                for (int v : cvols)
                                    if (gmsh::model::isInside(3, v, {cp->x[i], cp->y[i], cp->z[i]}) > 0)
                                        throw std::runtime_error("mesh3d: chord-bound corridor point (" +
                                            std::to_string(cp->x[i] * 1e3) + ", " + std::to_string(cp->y[i] * 1e3) + ", " +
                                            std::to_string(cp->z[i] * 1e3) + ") mm lies inside conductor volume " +
                                            std::to_string(v) + ": the centrelines disagree with the solids");
                            const double ct = cond_t;
                            gmsh::model::mesh::setSizeCallback([cp, ct](int, int, double x, double y, double z, double lc) {
                                return std::min(lc, cp->size_at(x, y, z, ct));
                            });
                            cmin = std::min(cmin, cp->hmin);
                            std::fprintf(stderr, "[mesh3d] local chord-bound refinement (3D centrelines): %zu corridor "
                                         "point(s) from %zu of %zu centreline samples; tightest %.4g mm (gap %.4g um) "
                                         "at (%.3f, %.3f, %.3f) mm; design target %.4g mm; %zu checked in air\n",
                                         cp->h.size(), cp->close, cp->samples, cp->hmin * 1e3, cp->dmin * 1e6,
                                         cp->at[0] * 1e3, cp->at[1] * 1e3, cp->at[2] * 1e3, cond_t * 1e3,
                                         (cp->h.size() + stride - 1) / stride);
                        } else if (!tight.empty()) {
                            const char axis = std::getenv("OMFEM_AXIS")
                                            ? std::getenv("OMFEM_AXIS")[0] : 'Y';
                            std::vector<int> cond_vols;
                            for (auto& [gname, gvs] : groups) {
                                if (gname.rfind("winding", 0) != 0 && gname.rfind("turn_", 0) != 0)
                                    continue;
                                for (int v : gvs) cond_vols.push_back(v);
                            }
                            if (cond_vols.empty())
                                throw std::runtime_error(
                                    "mesh3d: local chord-bound refinement has corridors to place "
                                    "but the model has no conductor volume to verify them against");
                            // Map a layout point onto its circle and test one azimuth.
                            auto at = [&](double lx, double ly, double phi, double p[3]) {
                                const double c = std::cos(phi), sn = std::sin(phi);
                                if (axis == 'X' || axis == 'x') { p[0] = ly; p[1] = lx * c; p[2] = lx * sn; }
                                else if (axis == 'Z' || axis == 'z') { p[0] = lx * c; p[1] = lx * sn; p[2] = ly; }
                                else { p[0] = lx * c; p[1] = ly; p[2] = lx * sn; }
                            };
                            auto in_copper = [&](const double p[3]) {
                                for (int v : cond_vols)
                                    if (gmsh::model::isInside(3, v, {p[0], p[1], p[2]}) > 0) return true;
                                return false;
                            };
                            // VERIFICATION 1 -- the mapping puts the TURN CENTRES inside copper.
                            // Not "at one azimuth": a stadium path around a rectangular column
                            // crosses any given radius near its corner arcs, so a single hit is
                            // no evidence at all (measured on 17_cllc_xfmr_e5528_3c92a, a
                            // rectangular-centre-leg E core: 42 of its 44 turn centres were
                            // inside copper at SOME azimuth while the path is not a circle).
                            // A genuinely revolved turn is inside at essentially EVERY azimuth --
                            // it is interrupted only by a crossing, a layer step or a lead
                            // cut-out -- so the discriminating quantity is the FRACTION, averaged
                            // over the winding so that a legitimately partial first/last turn
                            // cannot veto a correct mapping.
                            const int NAZ = 16;
                            const double MIN_MEAN = 0.70;
                            size_t hits_total = 0, dead = 0; double fx = 0, fy = 0;
                            for (const auto& cc : ls.turn_xy) {
                                int hits = 0;
                                for (int k = 0; k < NAZ; ++k) {
                                    double p[3]; at(cc[0], cc[1], 2.0 * M_PI * k / NAZ, p);
                                    if (in_copper(p)) ++hits;
                                }
                                hits_total += static_cast<size_t>(hits);
                                if (hits == 0) { ++dead; if (dead == 1) { fx = cc[0]; fy = cc[1]; } }
                            }
                            const size_t verified = ls.turn_xy.size();
                            const double mean_in = verified
                                ? static_cast<double>(hits_total) / (verified * NAZ) : 0.0;
                            if (dead || mean_in < MIN_MEAN)
                                throw std::runtime_error(
                                    "mesh3d: local chord-bound refinement cannot locate the "
                                    "layout in the model -- mapped as a circle of radius x about "
                                    "the " + std::string(1, axis) + " axis, the turn centres lie "
                                    "inside copper at only " + std::to_string(mean_in * 100.0) +
                                    "% of azimuths (need " + std::to_string(MIN_MEAN * 100.0) +
                                    "%), and " + std::to_string(dead) + " of " +
                                    std::to_string(verified) + " are inside at none (first at "
                                    "layout (" + std::to_string(fx * 1e3) + ", " +
                                    std::to_string(fy * 1e3) + ") mm). That mapping holds for a "
                                    "winding revolved about the column axis; it does not hold for "
                                    "a stadium/rectangular-column path or a toroid. Refusing to "
                                    "place the refinement by guesswork -- set OMFEM_CHORD_GLOBAL=1 "
                                    "to mesh this design with one size for every conductor.");
                            // VERIFICATION 2 -- a corridor is AIR. A mid-channel point that is
                            // inside copper at every azimuth is not a corridor; the mapping (or
                            // the pairing) is wrong and the refinement would land in the metal.
                            for (const auto& c : tight) {
                                bool air = false;
                                for (int k = 0; k < NAZ && !air; ++k) {
                                    double p[3]; at(c.x, c.y, 2.0 * M_PI * k / NAZ, p);
                                    air = !in_copper(p);
                                }
                                if (!air)
                                    throw std::runtime_error(
                                        "mesh3d: the corridor at layout (" + std::to_string(c.x * 1e3) +
                                        ", " + std::to_string(c.y * 1e3) + ") mm maps inside solid "
                                        "copper at every azimuth -- it is not an air channel, so "
                                        "the layout-to-model mapping is wrong. Refusing to place "
                                        "the refinement there.");
                                mvb::mesh::SizeFieldBuilder::RingBand b;
                                b.size_in = c.size;
                                b.radius = c.x;
                                b.axis_pos = c.y;
                                b.tube = c.reach;
                                b.thickness = 4.0 * c.size;   // grade back to the ceiling
                                bands.push_back(b);
                                cmin = std::min(cmin, c.size);
                            }
                            std::fprintf(stderr, "[mesh3d] local chord-bound refinement: %zu "
                                         "corridor ring(s) about the %c axis (tightest %.4g mm at "
                                         "r=%.4g mm, y=%.4g mm; design target %.4g mm), %zu turn "
                                         "centres map inside copper at %.0f%% of azimuths\n",
                                         bands.size(), axis,
                                         tight.front().size * 1e3, tight.front().x * 1e3,
                                         tight.front().y * 1e3, cond_t * 1e3, verified,
                                         mean_in * 100.0);
                            sizer.add_ring_bands(axis, bands);
                        }
                    }
                }
                // GAP REFINEMENT (opt-in, OMFEM_GAP_DIV_3D). Off by default because isotropic 3D
                // refinement of the wide thin gap disk is expensive -- a 5 mm-radius disk at 0.08 mm is
                // O(1e5) tets -- and the level-set path resolves gap fringing far more cheaply with an
                // anisotropic MMG metric. The gap being explicit conforming geometry does NOT by itself
                // make the magnetic circuit right, though: measured on e138 (0.8 mm subtractive gap), the
                // elements spanning the gap were 3.55 mm on the meshes small enough to solve directly and
                // still 1.14 mm on the 2.05M-tet mesh -- so one element bridges the gap and the flux
                // compression that DOMINATES this inductor's reluctance is not represented. That matters
                // for winding loss, because the gap fringing field sets the proximity loss in the nearby
                // turns. Symptom when it bites: the answer stops responding to the core at all (measured:
                // mu_r 1 -> 1109 moved R_ac by 5e-5 relative on the 91k mesh while moving it 13% on the
                // 39k one -- the magnetic circuit was not resolved on either, so the sensitivity was
                // numerical, not physical). Set OMFEM_GAP_DIV_3D=N to put >= N elements across each gap.
                if (const char* gd3 = std::getenv("OMFEM_GAP_DIV_3D")) {
                    const double gapdiv = std::atof(gd3);
                    if (!(gapdiv > 0.0))
                        throw std::runtime_error("mesh3d_from_mas: OMFEM_GAP_DIV_3D must be positive");
                    // Radial extent of the refined slab: the gap RETURN FLUX runs gap -> inner air ->
                    // winding window, so it must cover at least the centre post plus the fringing that
                    // reaches the innermost turn. A box tight to the post alone under-resolves that path
                    // (the level-set path records L collapsing 42% on PQ50 when it was tightened).
                    const double rgap = std::getenv("OMFEM_GAP_R_3D")
                                      ? std::atof(std::getenv("OMFEM_GAP_R_3D"))
                                      : (post_half3 > 0.0 ? 2.0 * post_half3 : 0.0);
                    for (const auto& g : core_gaps3) {
                        const double gt = g.len / gapdiv;
                        const double yb = 0.5 * g.len + 4.0 * gt;    // gap + a fringing band each side
                        if (!(rgap > 0.0)) {
                            std::fprintf(stderr, "[mesh3d] gap %.4gmm: no central-column half-width and no "
                                         "OMFEM_GAP_R_3D -- cannot size the refinement box, skipping\n",
                                         g.len*1e3);
                            continue;
                        }
                        std::fprintf(stderr, "[mesh3d] gap refinement: len=%.4gmm yc=%.4gmm -> target "
                                     "%.4gmm over r<%.4gmm\n", g.len*1e3, g.yc*1e3, gt*1e3, rgap*1e3);
                        sizer.add_box(gt, -rgap, g.yc - yb, -rgap, rgap, g.yc + yb, rgap, 6.0 * gt);
                        cmin = std::min(cmin, gt);
                    }
                }
                (void)post_half3;
                sizer.finalize(cmin);
                // ---- SKIN-DEPTH BOUNDARY LAYERS (OMFEM_SKIN_DELTA, 2026-09-11) ------------------
                // What Ansys Maxwell's skin-depth refinement and gmsh's advancing-layer method do:
                // mesh the copper SURFACE, extrude that surface mesh INWARD in layers whose thickness
                // grows geometrically from the skin depth (COMSOL: first layer delta/2 .. delta, 2-8
                // layers, coarse along the perimeter), drop the original copper volume and rebuild the
                // core from the closed shell of layer tops. Measured on the 000_debug bar at 1 MHz
                // (t/delta = 7): 4284 tets against 1.72M for the isotropic mesh at the same first-layer
                // size, gamma median 0.67, zero non-manifold faces. Every copper face is extruded --
                // caps included -- so the tops CLOSE; extruding the lateral faces alone left the old
                // volume under the layers and reproduced ABT #1154's "interior triangular face".
                // Parameters: OMFEM_SKIN_DELTA (m, required), OMFEM_SKIN_CELLS (first layer =
                // delta/cells, default 2), OMFEM_SKIN_GROWTH (1.3), OMFEM_SKIN_LAYERS (4),
                // OMFEM_SKIN_DMIN (m; the band is capped at 0.45*dmin so the layers never meet).
                if (skinMapped) {
                    // MAPPED COPPER: transfinite hex blocks (see the helpers above mesh3d_from_mas).
                    // OMFEM_SKIN_DELTA (m, required), OMFEM_SKIN_CELLS (first cell = delta/cells, 2),
                    // OMFEM_SKIN_GROWTH (1.3), OMFEM_SKIN_HLONG (m, along the wire, 1e-3).
                    const char* sd = std::getenv("OMFEM_SKIN_DELTA");
                    if (!sd) throw std::runtime_error("mesh3d_from_mas: OMFEM_SKIN_MAPPED needs OMFEM_SKIN_DELTA (skin depth, metres)");
                    const double delta  = std::atof(sd);
                    if (!(delta > 0.0)) throw std::runtime_error("mesh3d_from_mas: OMFEM_SKIN_DELTA must be > 0 (metres)");
                    const double cells  = std::getenv("OMFEM_SKIN_CELLS")  ? std::atof(std::getenv("OMFEM_SKIN_CELLS"))  : 2.0;
                    const double growth = std::getenv("OMFEM_SKIN_GROWTH") ? std::atof(std::getenv("OMFEM_SKIN_GROWTH")) : 1.3;
                    const double hlong  = std::getenv("OMFEM_SKIN_HLONG")  ? std::atof(std::getenv("OMFEM_SKIN_HLONG"))  : 1.0e-3;
                    const double hs = delta / cells;
                    std::map<double, BumpFit> fits;
                    auto fitFor = [&](double L) -> const BumpFit& {
                        auto it = fits.find(L);
                        if (it == fits.end()) it = fits.emplace(L, mapped_calibrate_bump(L, hs, growth)).first;
                        return it->second;
                    };
                    size_t nblocks = 0; std::string errs;
                    for (auto& [gname, gvols] : groups) {
                        auto dit = mappedDims.find(gname);
                        if (dit == mappedDims.end()) continue;
                        const double W = dit->second.first, T = dit->second.second;
                        const BumpFit& fw = fitFor(W); const BumpFit& ft = fitFor(T);
                        std::fprintf(stderr, "[mesh3d] mapped copper: %s W=%.4g mm -> %d cells (first %.4g mm x%.2f, max %.4g mm); "
                                     "T=%.4g mm -> %d cells (first %.4g mm x%.2f, max %.4g mm); along the wire %.4g mm\n",
                                     gname.c_str(), W*1e3, fw.n-1, fw.first*1e3, fw.grow, fw.maxCell*1e3,
                                     T*1e3, ft.n-1, ft.first*1e3, ft.grow, ft.maxCell*1e3, hlong*1e3);
                        for (int v : gvols) {
                            gmsh::vectorpair fdt; gmsh::model::getBoundary({{3, v}}, fdt, false, false, false);
                            if (fdt.size() != 6) { errs += " vol " + std::to_string(v) + " has " + std::to_string(fdt.size()) + " faces;"; continue; }
                            std::set<int> along; bool ok = true;
                            for (auto& f : fdt) {
                                gmsh::vectorpair cdt; gmsh::model::getBoundary({{2, std::abs(f.second)}}, cdt, false, false, false);
                                if (cdt.size() != 4) { ok = false; errs += " vol " + std::to_string(v) + " face " + std::to_string(std::abs(f.second)) + " has " + std::to_string(cdt.size()) + " edges;"; break; }
                                for (auto& c : cdt) {
                                    const int ct = std::abs(c.second); const MappedCurve mc = mapped_curve_info(ct);
                                    if (mapped_is_width(mc, W))      gmsh::model::mesh::setTransfiniteCurve(ct, fw.n, "Bump", fw.coef);
                                    else if (mapped_is_thick(mc, T)) gmsh::model::mesh::setTransfiniteCurve(ct, ft.n, "Bump", ft.coef);
                                    else along.insert(ct);
                                }
                            }
                            if (!ok) continue;
                            double lsum = 0.0; for (int c : along) { double l; gmsh::model::occ::getMass(1, c, l); lsum += l; }
                            const int nz = std::max(2, (int)std::lround(lsum / std::max<size_t>(1, along.size()) / hlong) + 1);
                            for (int c : along) gmsh::model::mesh::setTransfiniteCurve(c, nz);
                            for (auto& f : fdt) { gmsh::model::mesh::setTransfiniteSurface(std::abs(f.second)); gmsh::model::mesh::setRecombine(2, std::abs(f.second)); }
                            gmsh::model::mesh::setTransfiniteVolume(v); gmsh::model::mesh::setRecombine(3, v);
                            ++nblocks;
                        }
                    }
                    if (!errs.empty()) throw std::runtime_error("mesh3d_from_mas: mapped copper: blocks that are not hexahedral:" + errs);
                    std::fprintf(stderr, "[mesh3d] mapped copper: %zu transfinite hex blocks\n", nblocks);
                }

                else if (const char* sd = std::getenv("OMFEM_SKIN_DELTA")) {
                    const double delta  = std::atof(sd);
                    if (!(delta > 0.0)) throw std::runtime_error("mesh3d_from_mas: OMFEM_SKIN_DELTA must be > 0 (metres)");
                    const double cells  = std::getenv("OMFEM_SKIN_CELLS")  ? std::atof(std::getenv("OMFEM_SKIN_CELLS"))  : 2.0;
                    const double growth = std::getenv("OMFEM_SKIN_GROWTH") ? std::atof(std::getenv("OMFEM_SKIN_GROWTH")) : 1.3;
                    const int    nlay   = std::getenv("OMFEM_SKIN_LAYERS") ? std::atoi(std::getenv("OMFEM_SKIN_LAYERS")) : 4;
                    const double dmin   = std::getenv("OMFEM_SKIN_DMIN")   ? std::atof(std::getenv("OMFEM_SKIN_DMIN"))   : 0.0;
                    std::vector<double> heights; double h = delta / cells, acc = 0.0;
                    // BAND CAP. Inward layers from opposite faces of a thin conductor, and from the
                    // caps and the lateral faces at a lead end, cross each other when the band is
                    // too deep: on a faceted 0.8 mm round wire a 0.2175 mm band (0.54 R) gives
                    // "PLC Error: a segment and a facet intersect" at the cap edge while 0.14 mm
                    // (0.35 R) meshes cleanly with 1, 2 or 3 layers. 0.45*dmin was far too deep;
                    // 0.18*dmin (~0.36 R) is the measured safe side. OMFEM_SKIN_BANDFRAC overrides.
                    const double bandfrac = std::getenv("OMFEM_SKIN_BANDFRAC") ? std::atof(std::getenv("OMFEM_SKIN_BANDFRAC")) : 0.18;
                    for (int i = 0; i < nlay; ++i) {
                        if (dmin > 0.0 && acc + h > bandfrac * dmin) break;
                        acc += h; heights.push_back(-acc); h *= growth;
                    }
                    if (heights.empty()) throw std::runtime_error("mesh3d_from_mas: skin layers: first layer "
                        "exceeds 0.45*dmin -- the conductor is thinner than the layer it asked for");
                    std::string ht; for (double v : heights) ht += (ht.empty() ? "" : ", ") + std::to_string(-v * 1e3);
                    std::fprintf(stderr, "[mesh3d] skin layers: delta=%.4gmm first=%.4gmm x%.2f -> %zu layer(s), band %.4gmm "
                                 "(%.2f delta)  [cumulative mm: %s]\n", delta*1e3, delta/cells*1e3, growth, heights.size(),
                                 -heights.back()*1e3, -heights.back()/delta, ht.c_str());
                    // surface mesh first (the size field is active), then extrude it
                    mvb::mesh::generate_robust(2, mvb::mesh::algorithm_chain_2d(false));
                    std::vector<int> numEl(heights.size(), 1);
                    for (auto& [gname, gvols] : groups) {
                        if (gname.rfind("winding", 0) != 0 && gname.rfind("turn_", 0) != 0) continue;
                        std::vector<int> newvols;
                        for (int v : gvols) {
                            gmsh::vectorpair faces; gmsh::model::getBoundary({{3, v}}, faces, false, false, false);
                            for (auto& f : faces) f.second = std::abs(f.second);
                            gmsh::vectorpair ext;
                            // PRISMS, not split tets. gmsh refuses the tet subdivision on a real
                            // swept winding ("Unable to subdivide extruded mesh: change surface mesh
                            // or recombine extrusion instead" -- 03_buck, 150 layer volumes), and the
                            // split is what made the layers sliver on the bar anyway. Prisms need no
                            // pyramids: their triangular tops meet the core tets and their triangular
                            // bottoms meet the air tets directly. MFEM 4.10 solves them.
                            const bool prisms = !std::getenv("OMFEM_SKIN_SPLIT_TETS");
                            gmsh::model::geo::extrudeBoundaryLayer(faces, ext, numEl, heights, prisms, false, -1);
                            std::vector<int> tops;
                            for (size_t i = 1; i < ext.size(); ++i)
                                if (ext[i].first == 3) {
                                    newvols.push_back(ext[i].second);
                                    if (ext[i-1].first == 2) tops.push_back(ext[i-1].second);
                                }
                            gmsh::model::removeEntities({{3, v}}, false);
                            const int sl = gmsh::model::geo::addSurfaceLoop(tops);
                            newvols.push_back(gmsh::model::geo::addVolume({sl}));
                        }
                        gmsh::model::geo::synchronize();
                        // re-point the physical group at the layers + core
                        gmsh::vectorpair pgs; gmsh::model::getPhysicalGroups(pgs, 3);
                        for (auto& pgt : pgs) {
                            std::string nm; gmsh::model::getPhysicalName(3, pgt.second, nm);
                            if (nm == gname) { gmsh::model::removePhysicalGroups({pgt});
                                gmsh::model::addPhysicalGroup(3, newvols, pgt.second);
                                gmsh::model::setPhysicalName(3, pgt.second, gname); }
                        }
                        std::fprintf(stderr, "[mesh3d] skin layers: %s -> %zu volume(s) (%zu conductor(s) x %zu layers + core)\n",
                                     gname.c_str(), newvols.size(), gvols.size(), heights.size());
                        gvols = newvols;
                    }
                }
                // Mapped copper: Delaunay (1) only -- HXT refuses quad boundary faces -- and no
                // optimiser pass (it moves tet nodes, which the hex blocks share).
                // gmsh's built-in optimiser OFF for the hybrid Delaunay pass (measured on 03_buck:
                // with it, 389 triangles shared by 4 tets -- duplicate tets -- in the far air; without,
                // a clean mesh). OMFEM_MAPPED_OPT=1 turns it back on.
                if (skinMapped) {
                    gmsh::option::setNumber("Mesh.Optimize", std::getenv("OMFEM_MAPPED_OPT") ? 1 : 0);
                    // PYRAMID APEX CAP (patched gmsh, Mesh.PyramidApexMaxHeight). gmsh creates each
                    // pyramid apex inside the tet region at 0.4 x the quad's short edge; between two
                    // facing copper faces (stacked turns: pitch - thickness) or copper and core the
                    // fans must not cross the gap. The clearance is known here: the coating
                    // thickness on each face (outer - conducting)/2, the inter-turn gap being two of
                    // them and the copper-to-post clearance one. OMFEM_SKIN_APEX_MAX (m) overrides.
                    double gapMin = 1e30;
                    { OpenMagnetics::Coil coilG = enriched.get_coil();
                      const auto& fdG = coilG.get_functional_description();
                      for (size_t i = 0; i < fdG.size(); ++i) {
                          auto wire = coilG.resolve_wire(i);
                          const auto cw = wire.get_conducting_width(), ch = wire.get_conducting_height();
                          const auto ow = wire.get_outer_width(), oh = wire.get_outer_height();
                          if (cw && ow) gapMin = std::min(gapMin, 0.5 * (OpenMagnetics::resolve_dimensional_values(ow.value()) - OpenMagnetics::resolve_dimensional_values(cw.value())));
                          if (ch && oh) gapMin = std::min(gapMin, 0.5 * (OpenMagnetics::resolve_dimensional_values(oh.value()) - OpenMagnetics::resolve_dimensional_values(ch.value())));
                      } }
                    const double apexMax = std::getenv("OMFEM_SKIN_APEX_MAX") ? std::atof(std::getenv("OMFEM_SKIN_APEX_MAX"))
                                         : (gapMin < 1e29 && gapMin > 0.0 ? 0.35 * gapMin : 0.0);
                    if (!(apexMax > 0.0)) throw std::runtime_error("mesh3d_from_mas: mapped copper: no coating clearance to cap the pyramid apex height (wire has no outer dimensions); set OMFEM_SKIN_APEX_MAX (m)");
                    gmsh::option::setNumber("Mesh.PyramidApexMaxHeight", apexMax);
                    std::fprintf(stderr, "[mesh3d] mapped copper: pyramid apex height capped at %.4g mm (clearance %.4g mm)\n", apexMax*1e3, gapMin*1e3);
                }
                mvb::mesh::generate_robust(3, skinMapped ? std::vector<int>{1} : mvb::mesh::algorithm_chain_3d());
                if (!skinMapped) mvb::mesh::optimize_mesh(3);
                mvb::mesh::validate_regions_nonempty(3);   // throw if any volume region is empty
                // MANIFOLD CHECK (ABT #1154). gmsh can return a mesh it considers finished whose
                // topology is invalid: common_mode_choke_complete produced 1060597 tets with a
                // triangular face shared by THREE tetrahedra, and nothing upstream noticed. The CAD
                // is clean -- 3 solids, all watertight, and omfem_stepcheck --audit reports NO
                // OVERLAPS -- so this is a mesher defect, not a geometry one. MFEM finds it hours
                // later and aborts on load ("Invalid mesh topology. Interior triangular face..."),
                // which reads as a solver failure. Catch it here, where the retry ladder can try a
                // different size instead.
                {
                    std::vector<int> et; std::vector<std::vector<std::size_t>> en, ev;
                    gmsh::model::mesh::getElements(et, en, ev, 3);
                    std::unordered_map<std::string,int> faceUse;
                    std::size_t bad = 0; std::string witness;
                    for (size_t i = 0; i < et.size(); ++i) {
                        if (et[i] != 4) continue;                       // 4 = 4-node tetrahedron
                        const auto& nodes = ev[i];
                        for (size_t t = 0; t + 3 < nodes.size() + 0; t += 4) {
                            static const int fidx[4][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
                            for (auto& fi : fidx) {
                                std::size_t a = nodes[t+fi[0]], b = nodes[t+fi[1]], c = nodes[t+fi[2]];
                                if (a > b) std::swap(a, b);
                                if (b > c) std::swap(b, c);
                                if (a > b) std::swap(a, b);
                                std::string key = std::to_string(a) + "," + std::to_string(b) + ","
                                                + std::to_string(c);
                                if (++faceUse[key] == 3) {
                                    ++bad;
                                    if (witness.empty()) witness = key;
                                }
                            }
                        }
                    }
                    if (bad)
                        throw std::runtime_error(
                            "Invalid mesh topology: " + std::to_string(bad) + " triangular face(s) "
                            "shared by three or more tetrahedra (first: nodes " + witness + ") -- "
                            "gmsh returned a non-manifold volume mesh" +
                            ([&]{ if (!std::getenv("OMFEM_MESH_DEBUG")) return std::string();
                                  const std::string bad = opt.out_msh + ".bad.msh";
                                  try { gmsh::option::setNumber("Mesh.MshFileVersion", 2.2); gmsh::write(bad); } catch (...) {}
                                  return " [dumped to " + bad + "]"; }()));
                }
                break;
            } catch (const std::exception& e) {
                const std::string msg = e.what();
                // "Impossible to mesh periodic surface N" -> say WHAT surface N is (type + bbox +
                // owning region), else the number is undebuggable without a GUI session.
                {
                    static const std::regex periodic_re(
                        "(?:periodic surface|on surface) ([0-9]+)");
                    std::smatch pm;
                    if (std::regex_search(msg, pm, periodic_re)) {
                        const int stag = std::stoi(pm[1]);
                        try {
                            std::string stype;
                            gmsh::model::getType(2, stag, stype);
                            double bx0, by0, bz0, bx1, by1, bz1;
                            gmsh::model::getBoundingBox(2, stag, bx0, by0, bz0, bx1, by1, bz1);
                            std::fprintf(stderr, "[mesh3d] problem surface %d: type=%s bbox "
                                         "x[%.3f,%.3f] y[%.3f,%.3f] z[%.3f,%.3f] mm\n",
                                         stag, stype.c_str(), bx0*1e3, bx1*1e3, by0*1e3, by1*1e3,
                                         bz0*1e3, bz1*1e3);
                        } catch (...) {}
                    }
                }
                // NB: there used to be a special-case retry here that RAISED
                // Mesh.AngleToleranceFacetOverlap to 0.002 on a "surface N surface N"
                // self-overlap. It was removed: the baseline is already 1e-4 deg, so 0.002 flags
                // strictly MORE facet pairs than the attempt that just failed (the knob's meaning
                // is "overlapping when the dihedral is SMALLER than this"), its unanchored
                // backreference regex also matched "surface 12 surface 123", and the `continue`
                // consumed a retry without recording the error in the ladder. A folding facet at
                // micron clearance is a SIZE problem, and the ladder below already refines.
                const bool clearance_class =
                    msg.find("overlapping facets") != std::string::npos ||
                    msg.find("PLC Error") != std::string::npos ||
                    msg.find("Unable to recover the edge") != std::string::npos ||
                    msg.find("Invalid boundary mesh") != std::string::npos ||
                    msg.find("No elements in volume") != std::string::npos ||
                    msg.find("Invalid mesh topology") != std::string::npos ||
                    msg.find("empty mesh") != std::string::npos;
                ladder += (ladder.empty() ? "" : ", ") +
                          std::to_string(cond_t * 1e3) + "mm: " +
                          msg.substr(0, msg.find('\n'));
                if (std::getenv("OMFEM_SKIN_DELTA"))   // geometry already rewritten by the layer extrusion
                    throw std::runtime_error("mesh3d_from_mas: meshing with skin layers failed (no ladder retry "
                                             "-- the copper volumes were replaced in place): " + msg.substr(0, msg.find('\n')));
                if (!clearance_class || attempt + 1 >= max_attempts)
                    throw std::runtime_error(
                        "mesh3d_from_mas: meshing failed after " + std::to_string(attempt + 1) +
                        " attempt(s) [conductor-target ladder: " + ladder + "]");
                // No conductor target set -> seed the ladder from the core target so the retry
                // has something to refine (the clearance failures are at the winding surfaces).
                if (ladder_seed <= 0.0) ladder_seed = cond_t > 0.0 ? cond_t : opt.core_target;
                // ALTERNATE around the seed instead of only refining. The window is bounded on both
                // sides (see above), so a one-directional ladder can only ever find one of its two
                // edges -- 03_buck failed at 0.1111 mm with a PLC error, the ladder went finer to
                // 0.0722 mm and failed again, while 0.12 mm (COARSER) meshes this same design fine.
                // Multipliers apply to the seed, not compounded, so the search stays near the size
                // the physics asked for rather than running away.
                static const double kLadder[] = {1.15, 0.65, 1.35, 0.50};
                const int li = attempt < 4 ? attempt : 3;
                cond_t = kLadder[li] * ladder_seed;
                std::fprintf(stderr, "[mesh3d] RETRY %d/%d: %s -> conductor target %.4g mm\n",
                             attempt + 1, max_attempts - 1,
                             msg.substr(0, 90).c_str(), cond_t * 1e3);
            }
        }
    }
    const double mesh_secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_mesh0).count();
    std::size_t ntet = 0;
    { std::vector<int> et; std::vector<std::vector<std::size_t>> en, ev;
      gmsh::model::mesh::getElements(et, en, ev, 3);
      for (size_t i=0;i<et.size();++i) if (et[i]==4) ntet = en[i].size(); }
    std::fprintf(stderr, "[mesh3d] %zu tets  (mesh generation = %.1fs)\n", ntet, mesh_secs);
    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
    write_cad_edges(opt.out_msh);
    write_mesh_atomically(opt.out_msh);
    gmsh::model::remove();
    if (own_init) gmsh::finalize();
    if (!std::getenv("OMFEM_KEEP_STEP")) std::remove(step_path.c_str());   // keep it for review/prototyping
    return opt.out_msh;
}

}  // namespace mvb::mesh
