// mvbpp_meshheal <in.msh> <out.msh> [vmin_m3=1e-18] [passes=3]
// Moved verbatim from OMFEM tools/omfem_meshheal.cpp at 9482f50 (ABT #1588, step 7); only the tool name changed.
//
// HEAL the micron-scale sliver tets that break the hypre-AMS solve (ABT #1257): tets whose volume is below vmin carry
// stiffness entries ~edge^2/V and destroy the preconditioner -- 0-7 of them solve, 42+ stagnate, and 318 (10_emi_filter)
// returned 99 W where 0.2 W is right. They come from near-tangency / seam geometry, not from the size field, so a
// re-mesh gains or loses them by chance; healing keeps the mesh and costs a minute instead of an hour.
//
// Method: gmsh's own optimisers, which move interior nodes without touching the boundary (Netgen's optimiser then the
// default Gmsh one, repeated while the sliver count still falls). Counts before and after, region by region, and NEVER
// claims success it did not achieve: the exit code is non-zero if slivers remain, so the caller decides what to do.
#include <gmsh.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <array>
#include <string>
#include <vector>

static long count_slivers(double vmin, std::map<int,long>& perRegion) {
    perRegion.clear();
    std::vector<std::size_t> nt; std::vector<double> co, pa;
    gmsh::model::mesh::getNodes(nt, co, pa, -1, -1, false, false);
    std::map<std::size_t,std::size_t> idx;
    for (std::size_t i = 0; i < nt.size(); ++i) idx[nt[i]] = i;
    gmsh::vectorpair vols; gmsh::model::getEntities(vols, 3);
    long n = 0;
    for (const auto& dt : vols) {
        std::vector<int> et; std::vector<std::vector<std::size_t>> en, ev;
        gmsh::model::mesh::getElements(et, en, ev, 3, dt.second);
        for (std::size_t k = 0; k < et.size(); ++k) {
            if (et[k] != 4) continue;                       // linear tets only
            for (std::size_t e = 0; e + 3 < ev[k].size(); e += 4) {
                const double* p[4];
                for (int j = 0; j < 4; ++j) p[j] = &co[3 * idx[ev[k][e + j]]];
                const double a[3] = {p[1][0]-p[0][0], p[1][1]-p[0][1], p[1][2]-p[0][2]};
                const double b[3] = {p[2][0]-p[0][0], p[2][1]-p[0][1], p[2][2]-p[0][2]};
                const double c[3] = {p[3][0]-p[0][0], p[3][1]-p[0][1], p[3][2]-p[0][2]};
                const double V = std::fabs(a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6.0;
                if (V < vmin) { ++n; perRegion[dt.second]++; }
            }
        }
    }
    return n;
}

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: %s in.msh out.msh [vmin_m3=1e-18] [passes=3]\n", argv[0]); return 2; }
    const double vmin = (argc > 3) ? std::atof(argv[3]) : 1e-18;
    const int passes  = (argc > 4) ? std::atoi(argv[4]) : 3;
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", std::getenv("OMFEM_GMSH_VERBOSE") ? 1 : 0);
    gmsh::open(argv[1]);
    std::map<int,long> per;
    const long before = count_slivers(vmin, per);
    std::printf("[meshheal] %s: %ld tets below %.3g m3 before healing\n", argv[1], before, vmin);
    if (before == 0) { gmsh::write(argv[2]); gmsh::finalize(); std::printf("[meshheal] nothing to heal; copied\n"); return 0; }
    // PERTURBATION HEALING (2026-09-18). gmsh's optimisers segfault on our meshes (only the port/outer surface groups
    // are kept, so they have no complete boundary; the MMG-written skin meshes have none at all), and MMG output has
    // twice been unsolvable here. The literature's other family -- perturbing the sliver's vertices (Perturbing Slivers
    // in 3D Delaunay Meshes; the volume floor of the Tempered FEM, arXiv 2606.14301) -- is exactly what is free on OUR
    // geometry: every sliver we have measured sits in AIR, where moving a node by a micron changes no physics and no
    // boundary. So: move ONE free interior node of each sliver along the normal of its opposite face until the tet has
    // a real volume, and accept the move only if every tet touching that node stays positive and above the floor.
    // Nodes on any surface element, and nodes shared by more than one region, are never moved.
    std::vector<std::size_t> nt; std::vector<double> co, pa;
    gmsh::model::mesh::getNodes(nt, co, pa, -1, -1, false, false);
    std::map<std::size_t,std::size_t> idx;
    for (std::size_t i = 0; i < nt.size(); ++i) idx[nt[i]] = i;
    std::vector<char> frozen(nt.size(), 0);
    { gmsh::vectorpair fs; gmsh::model::getEntities(fs, 2);          // every node of a surface element is frozen
      for (const auto& dt : fs) { std::vector<int> et; std::vector<std::vector<std::size_t>> en, ev;
          gmsh::model::mesh::getElements(et, en, ev, 2, dt.second);
          for (const auto& v : ev) for (std::size_t t : v) { auto it = idx.find(t); if (it != idx.end()) frozen[it->second] = 1; } } }
    // region membership per node: a node touched by two different volumes is an interface node -> frozen
    std::vector<int> firstReg(nt.size(), 0);
    std::vector<std::array<std::size_t,4>> tets; std::vector<int> tetReg;
    { gmsh::vectorpair vols; gmsh::model::getEntities(vols, 3);
      for (const auto& dt : vols) { std::vector<int> et; std::vector<std::vector<std::size_t>> en, ev;
          gmsh::model::mesh::getElements(et, en, ev, 3, dt.second);
          for (std::size_t k = 0; k < et.size(); ++k) { if (et[k] != 4) continue;
              for (std::size_t e = 0; e + 3 < ev[k].size(); e += 4) {
                  std::array<std::size_t,4> T{};
                  for (int j = 0; j < 4; ++j) { T[j] = idx[ev[k][e+j]];
                      if (!firstReg[T[j]]) firstReg[T[j]] = dt.second;
                      else if (firstReg[T[j]] != dt.second) frozen[T[j]] = 1; }
                  tets.push_back(T); tetReg.push_back(dt.second); } } } }
    std::vector<std::vector<std::size_t>> incident(nt.size());
    for (std::size_t t = 0; t < tets.size(); ++t) for (int j = 0; j < 4; ++j) incident[tets[t][j]].push_back(t);
    auto vol = [&](const std::array<std::size_t,4>& T) {
        const double* p0 = &co[3*T[0]]; const double* p1 = &co[3*T[1]]; const double* p2 = &co[3*T[2]]; const double* p3 = &co[3*T[3]];
        const double a[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
        const double b[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
        const double c[3] = {p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2]};
        return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6.0; };
    const double vtarget = vmin * (std::getenv("OMFEM_HEAL_TARGET") ? std::atof(std::getenv("OMFEM_HEAL_TARGET")) : 20.0);
    long fixed_ = 0, stuck = 0, nofree = 0;
    for (int p = 0; p < passes; ++p) {
        long fixedThisPass = 0;
        for (std::size_t t = 0; t < tets.size(); ++t) {
            if (std::fabs(vol(tets[t])) >= vmin) continue;
            bool done = false, anyFree = false;
            for (int j = 0; j < 4 && !done; ++j) {
                const std::size_t n = tets[t][j];
                if (frozen[n]) continue;
                anyFree = true;
                // face opposite to j, its normal, and the height needed for vtarget
                const std::size_t f0 = tets[t][(j+1)%4], f1 = tets[t][(j+2)%4], f2 = tets[t][(j+3)%4];
                const double u[3] = {co[3*f1]-co[3*f0], co[3*f1+1]-co[3*f0+1], co[3*f1+2]-co[3*f0+2]};
                const double w[3] = {co[3*f2]-co[3*f0], co[3*f2+1]-co[3*f0+1], co[3*f2+2]-co[3*f0+2]};
                double nrm[3] = {u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0]};
                const double area2 = std::sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]+nrm[2]*nrm[2]);
                if (area2 <= 0) continue;
                for (int k = 0; k < 3; ++k) nrm[k] /= area2;
                const double need = 6.0 * vtarget / area2;              // height for the target volume
                const double keep[3] = {co[3*n], co[3*n+1], co[3*n+2]};
                // Move ALONG the direction that grows this tet's own signed volume (moving the other way would flip it),
                // and line-search the step down until every incident tet keeps its orientation and clears the floor.
                const double vc = vol(tets[t]);
                const double dir = (vc >= 0 ? 1.0 : -1.0);
                std::vector<double> before_(incident[n].size());
                for (std::size_t q = 0; q < incident[n].size(); ++q) before_[q] = vol(tets[incident[n][q]]);
                for (double f : {1.0, 0.6, 0.35, 0.2, 0.1, 0.05}) {
                    for (int k = 0; k < 3; ++k) co[3*n+k] = keep[k] + dir * f * need * nrm[k];
                    bool ok = std::fabs(vol(tets[t])) >= vmin;             // the sliver itself must be fixed
                    if (ok) for (std::size_t q = 0; q < incident[n].size(); ++q) {
                        const double v = vol(tets[incident[n][q]]);
                        // never flip an element, and never push another one below the floor it was above
                        if ((v > 0) != (before_[q] > 0) || (std::fabs(v) < vmin && std::fabs(before_[q]) >= vmin)) { ok = false; break; }
                    }
                    if (ok) { done = true; ++fixed_; ++fixedThisPass; break; }
                    for (int k = 0; k < 3; ++k) co[3*n+k] = keep[k];
                }
            }
            if (!done) { if (anyFree) ++stuck; else ++nofree; }
        }
        std::printf("[meshheal] pass %d: %ld slivers moved\n", p + 1, fixedThisPass);
        if (!fixedThisPass) break;
    }
    for (std::size_t i = 0; i < nt.size(); ++i)
        gmsh::model::mesh::setNode(nt[i], {co[3*i], co[3*i+1], co[3*i+2]}, {});
    const long now = count_slivers(vmin, per);
    std::printf("[meshheal] moved %ld, stuck %ld, no free node %ld -> %ld slivers left (was %ld)\n", fixed_, stuck, nofree, now, before);
    gmsh::write(argv[2]);
    gmsh::finalize();
    std::printf("[meshheal] wrote %s\n", argv[2]);
    return now > 0 ? 1 : 0;                                   // non-zero while any remain: the caller decides
}
