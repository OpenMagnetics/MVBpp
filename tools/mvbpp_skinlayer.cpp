// mvbpp_skinlayer <in.msh> <out.msh> [skin_depth_mm cells h_tan_mm] [--recipe r.json]
//   env: OMFEM_SKIN_DEEP (m), OMFEM_MMG_HGRAD -- required, no defaults; needs <in.msh>.path and <in.msh>.recipe.json
//
// Moved verbatim from OMFEM tools/omfem_skinlayer.cpp at 59cc050 (ABT #1588, step 5): the recipe
// functions are mvb::mesh's (moved in step 1), the tool and its provenance are named mvbpp_skinlayer.
// No logic changed.
//
// SKIN-LAYER REMESH of a CAD-conformal tet mesh (2026-09-08). What TRAFOLO's "dense boundary
// layers on massive wire" and Maxwell's "skin depth based" mesh operation do: take the valid
// isotropic conformal mesh (omfem_mesh3d, copper at ~1 cell per skin depth), and hand MMG3D an
// anisotropic metric that is FINE NORMAL TO THE COPPER SURFACE inside the conductors --
// h_n = delta/cells at the surface, growing geometrically with depth -- and unchanged along it.
// Outside the copper and away from it the local size is kept (measured per node from the input
// mesh), so the core/air/bobbin resolution the conformal mesher chose is preserved. Every volume
// and boundary physical group is carried through (interfaces between different volume refs are
// MMG boundaries; nosurf keeps them where they are).
//
// Distance-to-copper: exact point-triangle distance to the copper surface triangles (faces of
// copper tets not shared with another copper tet), bucketed on a uniform grid; the direction to
// the nearest surface point is the layer normal (at strip corners it is the corner direction,
// which is the right stretching there).
#include <gmsh.h>
#include "mvb/mesh/MeshRecipe.h"
#include <nlohmann/json.hpp>
#include <unistd.h>
extern "C" {
#include <mmg/mmg3d/libmmg3d.h>
}
#include <sys/resource.h>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <chrono>
#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

using clk = std::chrono::steady_clock;
static double since(clk::time_point t){ return std::chrono::duration<double>(clk::now()-t).count(); }

struct V3 { double x,y,z; };
static V3 sub(V3 a, V3 b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static V3 add(V3 a, V3 b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static V3 mul(V3 a, double s){ return {a.x*s,a.y*s,a.z*s}; }
static double dot(V3 a, V3 b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static V3 cross(V3 a, V3 b){ return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }
static double norm(V3 a){ return std::sqrt(dot(a,a)); }

// Closest point on triangle abc to p (Ericson, Real-Time Collision Detection 5.1.5).
static V3 closestOnTri(V3 p, V3 a, V3 b, V3 c) {
    V3 ab=sub(b,a), ac=sub(c,a), ap=sub(p,a);
    double d1=dot(ab,ap), d2=dot(ac,ap);
    if (d1<=0 && d2<=0) return a;
    V3 bp=sub(p,b); double d3=dot(ab,bp), d4=dot(ac,bp);
    if (d3>=0 && d4<=d3) return b;
    double vc=d1*d4-d3*d2;
    if (vc<=0 && d1>=0 && d3<=0){ double v=d1/(d1-d3); return add(a,mul(ab,v)); }
    V3 cp=sub(p,c); double d5=dot(ab,cp), d6=dot(ac,cp);
    if (d6>=0 && d5<=d6) return c;
    double vb=d5*d2-d1*d6;
    if (vb<=0 && d2>=0 && d6<=0){ double w=d2/(d2-d6); return add(a,mul(ac,w)); }
    double va=d3*d6-d5*d4;
    if (va<=0 && (d4-d3)>=0 && (d5-d6)>=0){ double w=(d4-d3)/((d4-d3)+(d5-d6)); return add(b,mul(sub(c,b),w)); }
    double denom=1.0/(va+vb+vc); double v=vb*denom, w=vc*denom;
    return add(a, add(mul(ab,v), mul(ac,w)));
}

int main(int argc_in, char** argv_in) {
    // MESH RECIPE (2026-09-13): --recipe <file.json> supplies the "layers" section (delta_m, cells, band_k,
    // h_tan_m, deep_m, growth, ...) through the OMFEM_LAYERS_* / OMFEM_SKIN_* environment; positional
    // arguments given explicitly still win. The effective layer recipe is written to <out.msh>.recipe.json.
    std::vector<char*> args; std::string recipePath;
    for (int i = 0; i < argc_in; ++i) {
        const std::string a = argv_in[i];
        if (a == "--recipe") { if (i + 1 >= argc_in) { std::fprintf(stderr, "--recipe needs a file\n"); return 2; } recipePath = argv_in[++i]; continue; }
        args.push_back(argv_in[i]);
    }
    const int argc = (int)args.size(); char** argv = args.data();
    if (!recipePath.empty()) {
        try {
            const auto over = mvb::mesh::apply_mesh_recipe(recipePath);
            std::printf("mesh recipe %s applied", recipePath.c_str());
            if (!over.empty()) { std::printf(" (overrides the environment:"); for (const auto& o : over) std::printf(" %s", o.c_str()); std::printf(")"); }
            std::printf("\n");
        } catch (const std::exception& e) { std::fprintf(stderr, "ERROR: %s\n", e.what()); return 2; }
    }
    // NO DEFAULTS (Alf 2026-09-28, "no magic numbers"). Every quantity below is either an INPUT the
    // caller states (the design plan records why it chose it) or DERIVED from the inputs and the mesh.
    // A missing input is refused, never filled with a constant.
    //   inputs:  skin depth delta, cells per delta at the wall, h_tan (size along the wire),
    //            h_deep (size cap through the thickness, OMFEM_SKIN_DEEP, metres), MMG gradation
    //            (OMFEM_MMG_HGRAD), the facet count of the conductors (from <in>.recipe.json).
    //   derived: layer count, growth, stack depth, search reach, air-side layer, feature angle,
    //            hmin, copper hausd, far-field hausd, hmax, MMG memory.
    auto envReq = [](const char* k, const char* what) -> double {
        const char* v = std::getenv(k);
        if (!v) throw std::runtime_error(std::string(k) + " is not set: " + what + " is an input with no default");
        char* end = nullptr; const double x = std::strtod(v, &end);
        if (end == v || *end != '\0') throw std::runtime_error(std::string(k) + "='" + v + "' is not a number");
        return x;
    };
    if (argc < 3) { std::fprintf(stderr,"usage: %s in.msh out.msh [skin_depth_mm cells h_tan_mm] [--recipe recipe.json]\n"
                                        "  env (required): OMFEM_SKIN_DEEP (m), OMFEM_MMG_HGRAD; the positional values may come from the recipe instead\n",argv[0]); return 2; }
    if (argc != 3 && argc != 6) { std::fprintf(stderr,"give all three of skin_depth_mm cells h_tan_mm, or none (recipe)\n"); return 2; }
    if (std::getenv("OMFEM_SKIN_GROWTH") || std::getenv("OMFEM_LAYERS_BAND")) {
        std::fprintf(stderr,"OMFEM_SKIN_GROWTH / OMFEM_LAYERS_BAND are gone: the growth and the band are derived from the layer count, h_surface and h_deep\n"); return 2; }
    double delta = 0, cells = 0, hTan = 0, hDeep = 0, hgrad = 0;
    try {
        delta = argc == 6 ? std::atof(argv[3])*1e-3 : envReq("OMFEM_LAYERS_DELTA", "the skin depth");
        cells = argc == 6 ? std::atof(argv[4])      : envReq("OMFEM_LAYERS_CELLS", "cells per skin depth at the wall");
        hTan  = argc == 6 ? std::atof(argv[5])*1e-3 : envReq("OMFEM_LAYERS_HTAN", "the element size along the wire");
        // Cap on the layer size in the NEAREST face's direction (through the thickness of a strip, or
        // to the centre of a round wire), in metres like the recipe's deep_m.
        hDeep = envReq("OMFEM_SKIN_DEEP", "the layer size cap through the conductor thickness (h_deep)");
        hgrad = envReq("OMFEM_MMG_HGRAD", "MMG's gradation");
    } catch (const std::exception& e) { std::fprintf(stderr, "ERROR: %s\n", e.what()); return 2; }
    if (!(delta>0) || !(cells>0) || !(hTan>0) || !(hDeep>0)) { std::fprintf(stderr,"skin depth, cells, h_tan and h_deep must be positive\n"); return 2; }
    if (!(hgrad > 1.0)) { std::fprintf(stderr, "OMFEM_MMG_HGRAD must be > 1 (got %.4g): with no gradation there is no growth to derive the layers from\n", hgrad); return 2; }
    const double hs = delta/cells;                 // normal size in the fine layer

    // NUMBER OF LAYERS (2026-09-20). The stack is defined by three quantities of which only two are
    // free, exactly as in every boundary-layer mesher: the first layer thickness hs = delta/cells,
    // the growth ratio g, and the number of layers N between the wall and the size at which the
    // layers hand over to the bulk (hDeep, the through-thickness cap). Asking for N derives g from
    // the other two -- hs*g^(N-1) = hDeep -- instead of leaving the layer count to fall out of a
    // growth ratio nobody chose for this geometry. OMFEM_SKIN_LAYERS=N; N and OMFEM_SKIN_GROWTH
    // together are contradictory and are refused rather than silently ranked.
    // PREDICTING N (2026-09-21). The layer count does not have to be guessed. The stack exists to
    // carry the element size from the wall (hs) out to h_deep, where the current is past e^-2 and
    // the layers have done their job; and MMG's gradation hgrad is by definition the largest ratio
    // it will tolerate between adjacent element sizes. A stack that grows FASTER than hgrad cannot
    // survive -- MMG smooths it back. A stack that grows SLOWER adds layers the gradation would
    // never have demanded, which is mesh nobody asked for. So the growth that costs least while
    // still reaching h_deep is exactly hgrad, and
    //       N = 1 + ceil( ln(h_deep/hs) / ln(hgrad) ).
    // MEASURED on ETD34 (cells 1, hgrad 3, h_deep 2 delta -> N = 2), six frequencies each:
    //       N=2  859,568 tets  ratio 1.127      N=4  1,177,706 tets  ratio 1.126
    //       N=3  1,021,885 tets ratio 1.127      N=6  1,402,394 tets  (running)
    // Every N gives the same answer to three digits and the cost rises monotonically, so the
    // smallest N that reaches h_deep is the right one -- which is what this formula returns.
    // Unset or 'auto' uses it; auto+K asks for K layers beyond it (a deliberate check); a plain N
    // asks for exactly N (the convergence runs above).
    if (!(hDeep > hs)) { std::fprintf(stderr, "h_deep (%.4g mm) is not coarser than h_surface (%.4g mm): nothing to grow across\n", hDeep*1e3, hs*1e3); return 2; }
    const char* envLayers = std::getenv("OMFEM_SKIN_LAYERS");
    const std::string ls = envLayers ? envLayers : "auto";
    int nLayers = 0;
    if (ls.rfind("auto", 0) == 0) {
        if (ls.size() > 4 && ls[4] != '+') { std::fprintf(stderr, "OMFEM_SKIN_LAYERS: expected 'auto', 'auto+K' or N, got '%s'\n", ls.c_str()); return 2; }
        const int plus = ls.size() > 4 ? std::atoi(ls.c_str() + 5) : 0;
        // hDeep > hs, so the log ratio is positive and derived >= 2: a wall layer and the deep layer.
        const int derived = 1 + (int)std::ceil(std::log(hDeep/hs) / std::log(hgrad));
        nLayers = derived + plus;
        std::printf("layers: N=%d derived from h_deep/h_surface=%.4g and hgrad=%.4g (1 + ceil(ln %.4g / ln %.4g) = %d)%s\n",
                    nLayers, hDeep/hs, hgrad, hDeep/hs, hgrad, derived,
                    plus ? (" plus " + std::to_string(plus) + " requested").c_str() : "");
    } else {
        nLayers = std::atoi(ls.c_str());
        if (nLayers < 2) { std::fprintf(stderr,"OMFEM_SKIN_LAYERS must be >= 2 (one layer is not a stack)\n"); return 2; }
    }
    const double growthK = std::pow(hDeep/hs, 1.0/(nLayers - 1));   // > 1 because hDeep > hs
    // Total thickness of the N-layer stack, hs*(g^N-1)/(g-1): the depth the layers actually occupy.
    const double band = hs*(std::pow(growthK, nLayers) - 1.0)/(growthK - 1.0);

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::open(argv[1]);

    std::vector<std::size_t> nodeTags; std::vector<double> coord, par;
    gmsh::model::mesh::getNodes(nodeTags, coord, par, -1, -1, false, false);
    const int np = (int)nodeTags.size();
    std::map<std::size_t,int> id; for (int i=0;i<np;i++) id[nodeTags[i]] = i+1;
    std::vector<V3> P(np);
    for (int i=0;i<np;i++) P[i] = {coord[3*i], coord[3*i+1], coord[3*i+2]};

    struct Tet { int v[4]; int ref; };
    std::vector<Tet> tets; std::vector<char> copperNode(np+1, 0); std::vector<char> solidNode(np+1, 0); std::vector<char> copperRef;
    std::vector<std::string> volNames;
    gmsh::vectorpair pgs; gmsh::model::getPhysicalGroups(pgs, 3);
    for (auto& pg : pgs) {
        std::string nm; gmsh::model::getPhysicalName(3, pg.second, nm);
        volNames.push_back(nm);
        const int ref = (int)volNames.size();
        const bool isCu = nm.rfind("winding",0)==0 || nm.rfind("turn",0)==0;
        copperRef.push_back(isCu);
        std::vector<int> ents; gmsh::model::getEntitiesForPhysicalGroup(3, pg.second, ents);
        for (int e : ents) {
            std::vector<int> ety; std::vector<std::vector<std::size_t>> etag, enod;
            gmsh::model::mesh::getElements(ety, etag, enod, 3, e);
            for (size_t t=0;t<ety.size();t++) if (ety[t]==4) {
                const auto& nn = enod[t];
                for (size_t k=0;k+3<nn.size();k+=4) {
                    Tet tt; for (int j=0;j<4;j++) tt.v[j]=id[nn[k+j]]; tt.ref=ref; tets.push_back(tt);
                    if (isCu) for (int j=0;j<4;j++) copperNode[tt.v[j]]=1;
                    // A node of a non-conducting SOLID (core, bobbin): the air-side layer metric
                    // must not be imposed there -- see the OMFEM_SKIN_OUTSIDE_SOLID note below.
                    if (!isCu && nm.rfind("air",0) != 0) for (int j=0;j<4;j++) solidNode[tt.v[j]]=1;
                }
            }
        }
    }
    // Boundary refs follow the volume refs (1..nvol). MMG gives a boundary triangle it CREATES either a
    // volume ref, 0, or the ref of the face it came from (mmg3d hash_3d.c MMG5_bdryTria), so every
    // triangle with a ref above nvol is one of ours and the rest are dropped on output.
    const int kBndBase = (int)volNames.size();
    struct Tri { int v[3]; int ref; };
    std::vector<Tri> tris; std::vector<std::string> bndNames;
    gmsh::vectorpair pgs2; gmsh::model::getPhysicalGroups(pgs2, 2);
    for (auto& pg : pgs2) {
        std::string nm; gmsh::model::getPhysicalName(2, pg.second, nm);
        bndNames.push_back(nm);
        const int ref = kBndBase + (int)bndNames.size();
        std::vector<int> ents; gmsh::model::getEntitiesForPhysicalGroup(2, pg.second, ents);
        for (int e : ents) {
            std::vector<int> ety; std::vector<std::vector<std::size_t>> etag, enod;
            gmsh::model::mesh::getElements(ety, etag, enod, 2, e);
            for (size_t t=0;t<ety.size();t++) if (ety[t]==2)
                for (size_t k=0;k+2<enod[t].size();k+=3) {
                    Tri tr; for (int j=0;j<3;j++) tr.v[j]=id[enod[t][k+j]]; tr.ref=ref; tris.push_back(tr);
                }
        }
    }
    gmsh::finalize();
    int nCuTet=0; for (auto& t : tets) if (copperRef[t.ref-1]) ++nCuTet;
    std::printf("loaded: %d nodes, %zu tets (%d copper), %zu boundary tris; delta=%.4g mm -> h_surface=%.4g mm, band=%.3g mm, h_deep=%.3g mm, h_tan<=%.3g mm\n",
                np, tets.size(), nCuTet, tris.size(), delta*1e3, hs*1e3, band*1e3, hDeep*1e3, hTan*1e3);
    std::printf("layers: %d from the wall to h_deep, growth %.4f derived (stack %.4g mm = band)\n", nLayers, growthK, band*1e3);

    // FEATURE ANGLE, from the geometry (2026-09-28). The swept conductors are n-gon prisms
    // (MVB++ --segments n, recorded by omfem_mesh3d in <in>.recipe.json); two neighbouring facets
    // meet at a crease of 360/n degrees, which is the faceting of a SMOOTH surface. The real edges
    // of a conductor -- a rectangular section's corners -- are 90 degrees. The threshold that tells
    // the two apart with the most margin on both sides is the midpoint, (360/n + 90)/2 (n = 12:
    // 60 deg; a true cylinder, n = 0, has no crease: 45 deg). Normals closer than this belong to one
    // face direction of the layer metric: that grouping is about the FIELD, which sees a round wire.
    int segments = -1;
    {
        const std::string rp = std::string(argv[1]) + ".recipe.json";
        std::ifstream rin(rp);
        if (!rin) { std::fprintf(stderr, "ERROR: %s not found: the conductor facet count (segments) is needed for the feature angle\n", rp.c_str()); return 2; }
        try { nlohmann::json r; rin >> r; segments = r.at("provenance").at("segments").get<int>(); }
        catch (const std::exception& e) { std::fprintf(stderr, "ERROR: %s has no provenance.segments: %s\n", rp.c_str(), e.what()); return 2; }
        if (segments < 0 || segments == 1 || segments == 2) { std::fprintf(stderr, "ERROR: %s: segments=%d is not a prism\n", rp.c_str(), segments); return 2; }
    }
    const double creaseDeg  = segments > 0 ? 360.0 / segments : 0.0;
    const double featureDeg = 0.5 * (creaseDeg + 90.0);
    if (!(creaseDeg < 90.0)) { std::fprintf(stderr, "ERROR: segments=%d gives %.4g-deg facet creases, not below the 90-deg real edges\n", segments, creaseDeg); return 2; }
    std::printf("feature angle: %.4g deg = midpoint of the %d-segment facet crease (%.4g deg) and a 90-deg edge\n", featureDeg, segments, creaseDeg);

    // --- centreline sidecar (<in>.path): wire tangents for ROUND wire (2026-09-12) ---
    // A cylinder has one normal direction; the two-perpendicular-faces rule that gives a strip's
    // tangent finds nothing, the size stays conformal in every direction and the copper count
    // explodes. mesh3d writes the real centrelines; the tangent at a node is the direction of the
    // nearest centreline segment.
    // REQUIRED (2026-09-28): the centreline is the ONLY source of the wire direction. The former
    // second source -- n1 x n2 of two faces judged "perpendicular" past a chosen 60 degrees -- and the
    // air-side copy from the nearest copper node are gone with their thresholds; the drawn centreline
    // is the exact geometry, for round wire and strips alike.
    struct CSeg { V3 a, b, u; double len; };
    std::vector<CSeg> csegs;
    {
        const std::string pp = std::string(argv[1]) + ".path";
        std::ifstream pf(pp);
        if (!pf) { std::fprintf(stderr, "ERROR: %s not found: the wire direction comes from the real-winding centrelines (omfem_mesh3d writes them)\n", pp.c_str()); return 2; }
        // format: <npaths>, then per path "name|nprims|wireRadius" and nprims lines "npt x y z ..."
        size_t npaths = 0; std::string line; bool ok = static_cast<bool>(pf >> npaths) && static_cast<bool>(std::getline(pf, line));
        for (size_t k = 0; ok && k < npaths; k++) {
            ok = static_cast<bool>(std::getline(pf, line));
            const size_t b1 = line.find('|'), b2 = ok ? line.find('|', b1 + 1) : std::string::npos;
            if (!ok || b1 == std::string::npos || b2 == std::string::npos) { ok = false; break; }
            const size_t nprims = std::stoul(line.substr(b1 + 1, b2 - b1 - 1));
            for (size_t q = 0; ok && q < nprims; q++) {
                ok = static_cast<bool>(std::getline(pf, line)); if (!ok) break;
                std::istringstream ss(line); size_t npt = 0; ok = static_cast<bool>(ss >> npt);
                V3 prev{0,0,0}; bool have = false;
                for (size_t j = 0; ok && j < npt; j++) { V3 c; ok = static_cast<bool>(ss >> c.x >> c.y >> c.z); if (!ok) break;
                    // a repeated point has no direction; it is not a segment
                    if (have) { V3 d = sub(c, prev); const double L = norm(d);
                        if (L > 0) csegs.push_back({prev, c, mul(d, 1.0/L), L}); }
                    prev = c; have = true; }
            }
        }
        if (!ok || csegs.empty()) { std::fprintf(stderr, "ERROR: %s is truncated or holds no segments\n", pp.c_str()); return 2; }
        std::printf("centreline sidecar: %zu segments from %s\n", csegs.size(), pp.c_str());
    }
    // RIDGES = THE CAD'S EDGES, exactly (2026-09-29). MMG keeps a ridge sharp and rebuilds a smooth surface
    // across every other surface edge. An angle threshold cannot say which is which here: the CAD's real
    // edges run from 1.4 deg (a junction mitre between MVB++ primitives) to 90 deg, and the base mesh's
    // own triangulation of a curved face bends by ~4 deg (0.4 mm triangles on a 5.6 mm turn) -- the two
    // overlap. Measured on 20_iso_buckboost: MMG's 45 deg (and a 60 deg midpoint) rounded the 30-deg facet
    // creases and moved the copper up to 26 um off the CAD (1,410 slivers); 15 deg kept those but still
    // rounded the mitres (589 slivers, flat tets lying on the copper); 0.69 deg froze the triangulation
    // itself and MMG never finished (> 4 h). omfem_mesh3d writes the mesh edges that lie on CAD curves
    // (<in>.cadedges); they are handed to MMG as ridges and its angle detection is switched off.
    std::vector<std::array<int,2>> cadEdges;
    {
        const std::string ep = std::string(argv[1]) + ".cadedges";
        std::ifstream ef(ep);
        if (!ef) { std::fprintf(stderr, "ERROR: %s not found: the CAD's edges come from omfem_mesh3d (re-mesh with a build that writes them)\n", ep.c_str()); return 2; }
        std::string magic; int ver = 0; size_t ne = 0;
        if (!(ef >> magic >> ver >> ne) || magic != "cadedges" || ver != 1) { std::fprintf(stderr, "ERROR: %s is not a 'cadedges 1' file\n", ep.c_str()); return 2; }
        cadEdges.reserve(ne);
        for (size_t k = 0; k < ne; ++k) {
            std::size_t a = 0, b = 0;
            if (!(ef >> a >> b)) { std::fprintf(stderr, "ERROR: %s is truncated at edge %zu of %zu\n", ep.c_str(), k, ne); return 2; }
            auto ia = id.find(a), ib = id.find(b);
            if (ia == id.end() || ib == id.end()) { std::fprintf(stderr, "ERROR: %s edge %zu names node %zu/%zu that %s does not have\n", ep.c_str(), k, a, b, argv[1]); return 2; }
            cadEdges.push_back({ia->second, ib->second});
        }
        std::printf("CAD edges: %zu mesh edges from %s -> MMG ridges (angle detection off)\n", cadEdges.size(), ep.c_str());
    }
    // Segment grid. The cell is the mean segment length (only speed depends on it); the search
    // below is EXACT for any cell: it widens ring by ring until nothing outside can be closer.
    double cseg_cell = 0.0; for (const auto& c : csegs) cseg_cell += c.len; cseg_cell /= (double)csegs.size();
    std::unordered_map<long long, std::vector<int>> cgrid;
    auto ckey = [&](int i, int j, int k){ return ((long long)i * 73856093LL) ^ ((long long)j * 19349663LL) ^ ((long long)k * 83492791LL); };
    auto ccell = [&](const V3& q, int& i, int& j, int& k){ i = (int)std::floor(q.x / cseg_cell); j = (int)std::floor(q.y / cseg_cell); k = (int)std::floor(q.z / cseg_cell); };
    int cgLo[3] = {std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
    int cgHi[3] = {std::numeric_limits<int>::lowest(), std::numeric_limits<int>::lowest(), std::numeric_limits<int>::lowest()};
    for (size_t s2 = 0; s2 < csegs.size(); s2++) {
        int i0,j0,k0,i1,j1,k1; ccell(csegs[s2].a, i0,j0,k0); ccell(csegs[s2].b, i1,j1,k1);
        const int lo3[3] = {std::min(i0,i1), std::min(j0,j1), std::min(k0,k1)}, hi3[3] = {std::max(i0,i1), std::max(j0,j1), std::max(k0,k1)};
        for (int a = 0; a < 3; ++a) { cgLo[a] = std::min(cgLo[a], lo3[a]); cgHi[a] = std::max(cgHi[a], hi3[a]); }
        for (int i = lo3[0]; i <= hi3[0]; i++) for (int j = lo3[1]; j <= hi3[1]; j++) for (int k = lo3[2]; k <= hi3[2]; k++) cgrid[ckey(i,j,k)].push_back((int)s2);
    }
    // Nearest centreline segment to q. After scanning every cell within Chebyshev ring r of q's
    // cell, any segment not yet seen is at least r*cell away; stop once the best is within that.
    auto nearestTangent = [&](const V3& q, V3& t, double& dist) {
        int ci, cj, ck; ccell(q, ci, cj, ck);
        double best = std::numeric_limits<double>::infinity(); int bi = -1;
        const int rMax = std::max({std::abs(ci - cgLo[0]), std::abs(ci - cgHi[0]), std::abs(cj - cgLo[1]), std::abs(cj - cgHi[1]),
                                   std::abs(ck - cgLo[2]), std::abs(ck - cgHi[2])});
        for (int r = 0; r <= rMax; ++r) {
            for (int di=-r; di<=r; di++) for (int dj=-r; dj<=r; dj++) for (int dk=-r; dk<=r; dk++) {
                if (std::max({std::abs(di), std::abs(dj), std::abs(dk)}) != r) continue;   // this ring's shell only
                auto it = cgrid.find(ckey(ci+di, cj+dj, ck+dk)); if (it == cgrid.end()) continue;
                for (int s2 : it->second) { const CSeg& c = csegs[s2]; const V3 ap = sub(q, c.a); double tt = dot(ap, c.u); tt = std::max(0.0, std::min(c.len, tt));
                    const V3 cp = add(c.a, mul(c.u, tt)); const double d2 = dot(sub(q,cp), sub(q,cp)); if (d2 < best) { best = d2; bi = s2; } }
            }
            if (bi >= 0 && std::sqrt(best) <= r * cseg_cell) break;
        }
        t = csegs[bi].u; dist = std::sqrt(best);   // csegs is not empty, so bi >= 0 after the last ring
    };
    double cuCentrelineMax = 0.0;   // farthest copper node from its centreline (reported: a coverage check)

    // --- local size per node from the input mesh (mean edge length of the incident tets) ---
    std::vector<double> hloc(np+1, 0.0), hcnt(np+1, 0.0);
    for (auto& t : tets) {
        double s=0; int c=0;
        for (int i=0;i<4;i++) for (int j=i+1;j<4;j++){ s+=norm(sub(P[t.v[i]-1],P[t.v[j]-1])); ++c; }
        const double h = s/c;
        for (int j=0;j<4;j++){ hloc[t.v[j]]+=h; hcnt[t.v[j]]+=1; }
    }
    // A node no tet uses has no size to keep; MMG would be handed one from nowhere. Refuse the mesh.
    { long orphan = 0; for (int i=1;i<=np;i++) { if (hcnt[i]>0) hloc[i]/=hcnt[i]; else ++orphan; }
      if (orphan) { std::fprintf(stderr, "ERROR: %s has %ld node(s) that belong to no tetrahedron\n", argv[1], orphan); return 1; } }

    // --- per-PHYSICAL-VOLUME size measured from the input mesh (for MMG's local parameters) ---
    // The coarsest element the conformal mesher put in a region IS the size that region asked for
    // (air 6 mm, core 2 mm in the corpus recipe); a skin-layer pass must not make it finer.
    std::vector<double> hRefMax(volNames.size() + 1, 0.0);
    std::vector<long>   nRefTet(volNames.size() + 1, 0);
    for (auto& t : tets) {
        double s=0; int c=0;
        for (int i=0;i<4;i++) for (int j=i+1;j<4;j++){ s+=norm(sub(P[t.v[i]-1],P[t.v[j]-1])); ++c; }
        hRefMax[t.ref] = std::max(hRefMax[t.ref], s/c);
        ++nRefTet[t.ref];
    }

    // --- copper surface triangles: faces of copper tets not shared with another copper tet ---
    struct FaceKey { int a,b,c; bool operator==(const FaceKey& o) const { return a==o.a&&b==o.b&&c==o.c; } };
    struct FaceHash { size_t operator()(const FaceKey& k) const { return ((size_t)k.a*73856093u) ^ ((size_t)k.b*19349663u) ^ ((size_t)k.c*83492791u); } };
    std::unordered_map<FaceKey,int,FaceHash> fcount; fcount.reserve(4*nCuTet);
    static const int fidx[4][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
    for (auto& t : tets) { if (!copperRef[t.ref-1]) continue;
        for (int f=0;f<4;f++){ int v[3]={t.v[fidx[f][0]],t.v[fidx[f][1]],t.v[fidx[f][2]]}; std::sort(v,v+3); fcount[{v[0],v[1],v[2]}]++; } }
    std::vector<std::array<int,3>> sf; sf.reserve(fcount.size()/2);
    for (auto& kv : fcount) if (kv.second==1) sf.push_back({kv.first.a,kv.first.b,kv.first.c});
    std::printf("copper surface: %zu triangles\n", sf.size());

    // --- grid of surface triangles (by centroid) ---
    V3 lo{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    V3 hi{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    for (auto& p : P){ lo.x=std::min(lo.x,p.x); lo.y=std::min(lo.y,p.y); lo.z=std::min(lo.z,p.z); hi.x=std::max(hi.x,p.x); hi.y=std::max(hi.y,p.y); hi.z=std::max(hi.z,p.z); }
    double hCuIn = 0.0; for (int i=1;i<=np;i++) if (copperNode[i]) hCuIn = std::max(hCuIn, hloc[i]);
    // SEARCH REACH, derived (2026-09-28; was max(2 mm, min(6 band, ...)) with a 0.05 growth floor).
    // A face direction at distance d asks for the normal size hs + (g-1)*d (air side: hgrad*hs +
    // (g-1)*d, below), and it changes the metric only while that is finer than the size already
    // there, which is at most h_tan. So nothing past (h_tan - hs)/(g-1) can contribute: that is the
    // reach, exactly. The one exception is the copper node's NEAREST face, capped at h_deep and so
    // active at any depth; it is found by an exact search below (nearestSurface).
    const double reach = (hTan - hs) / (growthK - 1.0);
    if (!(reach > 0)) { std::fprintf(stderr, "h_tan (%.4g mm) is not coarser than h_surface (%.4g mm): no layer can form\n", hTan*1e3, hs*1e3); return 2; }
    // The grid is searched one cell around the node and holds each triangle at its centroid; a
    // triangle's centroid lies within rho (its farthest vertex from the centroid) of every point of
    // it, so a cell of reach + rho_max finds every triangle within reach.
    double rhoMax = 0.0;
    for (const auto& t3 : sf) { const V3 c = mul(add(add(P[t3[0]-1],P[t3[1]-1]),P[t3[2]-1]),1.0/3.0);
        for (int j=0;j<3;j++) rhoMax = std::max(rhoMax, norm(sub(P[t3[j]-1], c))); }
    const double cell = reach + rhoMax;
    std::printf("search reach %.3f mm = (h_tan - h_surface)/(growth - 1); grid cell %.3f mm (+ largest surface triangle radius); coarsest copper input size %.3f mm\n",
                reach*1e3, cell*1e3, hCuIn*1e3);
    const int nx=(int)((hi.x-lo.x)/cell)+1, ny=(int)((hi.y-lo.y)/cell)+1, nz=(int)((hi.z-lo.z)/cell)+1;
    auto cidx=[&](V3 p, int& i, int& j, int& k){ i=std::min(nx-1,std::max(0,(int)((p.x-lo.x)/cell))); j=std::min(ny-1,std::max(0,(int)((p.y-lo.y)/cell))); k=std::min(nz-1,std::max(0,(int)((p.z-lo.z)/cell))); };
    std::unordered_map<long long, std::vector<int>> grid; grid.reserve(sf.size());
    auto key=[&](int i,int j,int k){ return ((long long)i*ny + j)*nz + k; };
    for (size_t s=0;s<sf.size();s++){ V3 c = mul(add(add(P[sf[s][0]-1],P[sf[s][1]-1]),P[sf[s][2]-1]),1.0/3.0); int i,j,k; cidx(c,i,j,k); grid[key(i,j,k)].push_back((int)s); }

    // Nearest copper-surface triangle at ANY distance (exact ring search over the same grid). A face
    // beyond reach is irrelevant for every direction but one: inside copper the nearest face's
    // direction is capped at h_deep, so it refines at every depth -- a copper node deeper than reach
    // (the core of a thick strip) still needs it. After scanning ring r, an unseen triangle's
    // centroid is at least r*cell away, so the triangle itself at least r*cell - rho_max.
    auto nearestSurface = [&](const V3& p, V3& n, double& dist) {
        int ci,cj,ck; cidx(p,ci,cj,ck);
        const int rMax = std::max({ci, nx-1-ci, cj, ny-1-cj, ck, nz-1-ck});
        double best = std::numeric_limits<double>::infinity(); int bs = -1;
        for (int r = 0; r <= rMax; ++r) {
            for (int di=-r;di<=r;di++) for (int dj=-r;dj<=r;dj++) for (int dk=-r;dk<=r;dk++) {
                if (std::max({std::abs(di), std::abs(dj), std::abs(dk)}) != r) continue;
                const int i2 = ci+di, j2 = cj+dj, k2 = ck+dk;
                if (i2 < 0 || j2 < 0 || k2 < 0 || i2 >= nx || j2 >= ny || k2 >= nz) continue;
                auto it = grid.find(key(i2,j2,k2)); if (it==grid.end()) continue;
                for (int s : it->second) {
                    const double dd = norm(sub(p, closestOnTri(p, P[sf[s][0]-1], P[sf[s][1]-1], P[sf[s][2]-1])));
                    if (dd < best) { best = dd; bs = s; }
                }
            }
            if (bs >= 0 && best <= r*cell - rhoMax) break;
        }
        if (bs < 0) return false;
        const V3 a=P[sf[bs][0]-1], b=P[sf[bs][1]-1], c=P[sf[bs][2]-1];
        V3 nn = cross(sub(b,a), sub(c,a)); const double L = norm(nn);
        if (!(L > 0)) return false;
        n = mul(nn, 1.0/L); dist = best; return true;
    };

    // --- metric ---
    auto t0 = clk::now();
    MMG5_pMesh mmg=nullptr; MMG5_pSol sol=nullptr;
    MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh,&mmg, MMG5_ARG_ppMet,&sol, MMG5_ARG_end);
    MMG3D_Set_meshSize(mmg, np, (int)tets.size(), 0, (int)tris.size(), 0, (int)cadEdges.size());
    for (int i=0;i<np;i++) MMG3D_Set_vertex(mmg, P[i].x,P[i].y,P[i].z, 0, i+1);
    for (size_t k=0;k<tets.size();k++) MMG3D_Set_tetrahedron(mmg, tets[k].v[0],tets[k].v[1],tets[k].v[2],tets[k].v[3], tets[k].ref, (int)k+1);
    // NON-CONDUCTING SOLIDS ARE FROZEN (2026-10-01, ABT #1564). The core and bobbin carry no skin layer and
    // keep the conformal mesher's size (see OUTSIDE THE COPPER below), but leaving their tets to MMG let it
    // remesh their surface anyway, driven by the fine air-side metric across a thin gap: buck_inductor's
    // core grew 1,094 -> ~33k tets, and MMG split the core's straight vertical edges into collinear ridge
    // points and left flat tets on them -- 160 tets of EXACTLY zero volume (78 air + 78 core sharing a
    // zero-area triangle of three collinear points on a core corner line; smd_drum 44, same class). The base
    // mesh has none. Required tets keep their vertices, edges and faces as given, so the core/air interface
    // stays the base mesh's (already sliver-free) triangulation and the air builds its layers against it.
    // OMFEM_SKIN_FREE_SOLID=1 lets MMG remesh the solids again (A/B only).
    const bool freeSolid = std::getenv("OMFEM_SKIN_FREE_SOLID") && std::atoi(std::getenv("OMFEM_SKIN_FREE_SOLID")) != 0;
    long nFrozen = 0;
    if (!freeSolid)
        for (size_t k=0;k<tets.size();k++) {
            const std::string& nm = volNames[tets[k].ref-1];
            if (copperRef[tets[k].ref-1] || nm.rfind("air",0) == 0) continue;
            if (!MMG3D_Set_requiredTetrahedron(mmg, (int)k+1)) { std::fprintf(stderr, "ERROR: MMG3D refused to freeze solid tet %zu\n", k+1); return 1; }
            ++nFrozen;
        }
    std::printf("solids: %ld tet(s) of non-conducting solids frozen%s\n", nFrozen, freeSolid ? " (OMFEM_SKIN_FREE_SOLID: none, A/B)" : "");
    for (size_t k=0;k<tris.size();k++) MMG3D_Set_triangle(mmg, tris[k].v[0],tris[k].v[1],tris[k].v[2], tris[k].ref, (int)k+1);
    for (size_t k=0;k<cadEdges.size();k++) {
        if (!MMG3D_Set_edge(mmg, cadEdges[k][0], cadEdges[k][1], 0, (int)k+1) || !MMG3D_Set_ridge(mmg, (int)k+1)) {
            std::fprintf(stderr, "ERROR: MMG3D refused CAD edge %zu\n", k+1); return 1; } }
    MMG3D_Set_solSize(mmg, sol, MMG5_Vertex, np, MMG5_Tensor);
    long nLayer=0, nOut=0, nFar=0;
    // AIR SIDE (2026-09-11). Left isotropic at the input size, the air next to the conductor CAPS
    // the tangential size of the copper surface through gradation: the bar asked 1 mm along the
    // wire and got 0.15-0.18 (= the air size) whatever hgrad/hmax/hausd/input coarseness was tried.
    // So the air near copper gets an anisotropic layer metric too: coarse along the wire, and a
    // normal size that grows away from the copper.
    // DERIVED (2026-09-28; were OMFEM_SKIN_OUTSIDE = 2 and OMFEM_SKIN_OUTSIDE_BAND = 4 bands). The
    // air node nearest the wall neighbours a copper surface node of normal size hs, and MMG's
    // gradation lets it be at most hgrad times that; so the air normal size starts at hgrad*hs and
    // grows like the copper layers, hn = hgrad*hs + (g-1)*d, capped at the node's own size. Its
    // extent is where that stays finer than what is already there, which the reach bounds exactly
    // -- no band count. The inter-turn gaps stay resolved by the normal direction (both facing
    // copper faces are one merged direction).
    const bool outsideSolid = std::getenv("OMFEM_SKIN_OUTSIDE_SOLID") && std::atoi(std::getenv("OMFEM_SKIN_OUTSIDE_SOLID")) != 0;
    // CORNER-AWARE METRIC (2026-09-08). A rectangular conductor crowds its current into the four
    // CORNERS (the edge effect that distinguishes it from round wire), and a corner is the
    // INTERSECTION of two boundary layers: both face normals must be fine there. Taking only the
    // nearest surface point gives one fine direction along the diagonal and leaves both face
    // normals under-resolved -- wrong exactly where the current density peaks. So collect the
    // nearby surface triangles, CLUSTER them by normal direction, and intersect the per-face layer
    // metrics: M = (1/ht^2) I + sum_k [1/hn(d_k)^2 - 1/ht^2] n_k (x) n_k. On a flat face this
    // collapses to the single-normal metric (one cluster); on a rectangular edge it gives two fine
    // directions, at a corner three. COMSOL calls the equivalent step "corner refinement".
    // Two normals are one direction when they are closer than the FEATURE ANGLE (derived above
    // from the facet count; was a chosen 0.7 = 45.6 deg): 20 deg split the 30-deg facets of a
    // 12-segment sweep into separate directions and the metric went near-isotropic (ABT #1156).
    // Every distinct direction counts (the former cap of 3 per node is gone).
    const double cosMerge = std::cos(featureDeg * M_PI / 180.0);
    long nCorner = 0, nTan = 0, nTanAir = 0;
    // OMFEM_SKIN_STATS=1: where do the multi-direction nodes come from? Histogram of the angle
    // between the two FINEST directions and of the second direction's distance (ABT #1156 lever).
    // The bins are only how the diagnostic prints.
    const bool stats = std::getenv("OMFEM_SKIN_STATS") != nullptr;
    long angHist[10] = {0}, d2Hist[8] = {0}, nfineHist[4] = {0};
    // WIRE DIRECTION (2026-09-11). Coarse ALONG the wire, conformal-sized ACROSS the section: the
    // cut across a strip must still look like a boundary-layer mesh (thin layers, moderate width),
    // only the direction in which nothing varies may be long. The direction is the nearest
    // centreline segment's (sidecar above). Letting h_tan loose in-plane drew 1 mm slivers across
    // the width in the cut (bar_v5).
    double hMetricMin = std::numeric_limits<double>::infinity();   // finest size handed to MMG (-> hmin)
    for (int i=0;i<np;i++) {
        const bool inCu = copperNode[i+1];
        const V3 p = P[i];
        int ci,cj,ck; cidx(p,ci,cj,ck);
        // per-normal-cluster nearest distance, over the faces within reach
        struct Dir { V3 n; double d; };
        std::vector<Dir> dirs;
        for (int di=-1;di<=1;di++) for (int dj=-1;dj<=1;dj++) for (int dk=-1;dk<=1;dk++) {
            auto it = grid.find(key(ci+di,cj+dj,ck+dk)); if (it==grid.end()) continue;
            for (int s : it->second) {
                const V3 a=P[sf[s][0]-1], b=P[sf[s][1]-1], c=P[sf[s][2]-1];
                const V3 cp = closestOnTri(p, a, b, c);
                const double dd = norm(sub(p,cp));
                if (dd > reach) continue;                          // cannot change the metric (see reach)
                V3 nn = cross(sub(b,a), sub(c,a)); const double L = norm(nn);
                if (!(L > 0)) continue;                            // a zero-area triangle has no normal
                nn = mul(nn, 1.0/L);
                bool merged = false;
                for (auto& D : dirs)
                    if (std::fabs(dot(D.n, nn)) > cosMerge) { if (dd < D.d) { D.d = dd; D.n = nn; } merged = true; break; }
                if (!merged) dirs.push_back({nn, dd});
            }
        }
        if (inCu && dirs.empty()) {                               // deeper than reach: the nearest face still counts
            V3 n; double dn = 0.0;
            if (!nearestSurface(p, n, dn)) { std::fprintf(stderr, "ERROR: copper node %d has no copper surface triangle with a normal\n", i+1); return 1; }
            dirs.push_back({n, dn});
        }
        const double hl = hloc[i+1];
        double m11,m12,m13,m22,m23,m33;
        // OUTSIDE THE COPPER, ONLY IN THE AIR (2026-09-20). The air-side layer metric was added so
        // that the AIR next to a conductor does not cap the copper's tangential size through
        // gradation; it was applied to every non-copper node, so on a toroid -- where the core is a
        // fraction of a millimetre from the winding over its whole surface -- the CORE got the
        // conductor's layer sizes. Measured on common_mode_choke (corpus v7 base, local parameters
        // on): core 1,917 -> 195,953 tets, 102x, in a conductor skin-layer pass. A node that
        // belongs to a non-conducting SOLID (core, bobbin) keeps that solid's own size; the air
        // nodes in the gap still carry the layer. OMFEM_SKIN_OUTSIDE_SOLID=1 restores the old
        // behaviour (A/B only).
        const bool outsideOK = inCu || !solidNode[i+1] || outsideSolid;
        if (!dirs.empty() && outsideOK) {
            std::sort(dirs.begin(), dirs.end(), [](const Dir& a, const Dir& b){ return a.d < b.d; });
            V3 t; double dc = 0.0; nearestTangent(p, t, dc);
            if (inCu) { cuCentrelineMax = std::max(cuCentrelineMax, dc); ++nTan; } else ++nTanAir;
            // TANGENTIAL sizes: across the section the conformal size, capped at h_tan in copper
            // (air keeps its own size); along the wire h_tan when that is coarser.
            const double hIn = inCu ? std::min(hl, hTan) : hl, aIn = 1.0/(hIn*hIn);
            hMetricMin = std::min(hMetricMin, hIn);
            m11=m22=m33=aIn; m12=m13=m23=0.0;
            const bool alongCoarser = hTan > hIn;
            if (alongCoarser) {
                const double bt = 1.0/(hTan*hTan) - aIn;          // negative: coarser along t
                m11 += bt*t.x*t.x; m22 += bt*t.y*t.y; m33 += bt*t.z*t.z;
                m12 += bt*t.x*t.y; m13 += bt*t.x*t.z; m23 += bt*t.y*t.z;
            }
            int nfine = 0;
            const double dNear = dirs[0].d;
            for (const auto& D : dirs) {
                // GROWTH FROM THE SURFACE. Layer k has thickness hs*growth^k and starts at depth
                // hs*(growth^k-1)/(growth-1), i.e. the local layer thickness at depth d is
                // hs + (growth-1)*d. Caps: the nearest face's direction at hDeep (through the
                // thickness); any other face's direction -- across the width of a strip, from its
                // side faces -- at the in-plane size, so the middle of the width keeps the
                // conformal size and the edges get graded layers (ABT #1156).
                double hn = hs + (growthK - 1.0) * D.d;
                if (inCu) hn = std::min(hn, D.d <= dNear ? hDeep : std::max(hIn, hDeep));
                else      hn = std::min(hIn, hgrad*hs + (growthK - 1.0) * D.d);
                // the normal direction: remove whatever the tangential terms put there, set 1/hn^2
                const double nt = dot(D.n, t);
                const double cur = aIn + (alongCoarser ? (1.0/(hTan*hTan) - aIn) * nt*nt : 0.0);
                const double bb = 1.0/(hn*hn) - cur;
                if (bb <= 0) continue;                            // this face is too far to refine toward
                ++nfine;
                hMetricMin = std::min(hMetricMin, hn);
                m11 += bb*D.n.x*D.n.x; m22 += bb*D.n.y*D.n.y; m33 += bb*D.n.z*D.n.z;
                m12 += bb*D.n.x*D.n.y; m13 += bb*D.n.x*D.n.z; m23 += bb*D.n.y*D.n.z;
            }
            if (nfine > 1) ++nCorner;
            if (stats) {
                ++nfineHist[std::min(nfine, 3)];
                if (dirs.size() > 1) {
                    const double cs = std::fabs(dot(dirs[0].n, dirs[1].n));
                    const double ang = std::acos(std::min(1.0, cs)) * 180.0 / M_PI;   // 0..90
                    ++angHist[std::min(9, (int)(ang / 10.0))];
                    ++d2Hist[std::min(7, (int)(dirs[1].d / (0.5*band)))];
                }
            }
            if (inCu) ++nLayer; else ++nOut;
        } else { const double a = 1.0/(hl*hl); m11=m22=m33=a; m12=m13=m23=0; ++nFar; hMetricMin = std::min(hMetricMin, hl); }
        MMG3D_Set_tensorSol(sol, m11,m12,m13,m22,m23,m33, i+1);
    }
    std::printf("metric: wire tangent from the centreline at %ld copper and %ld air nodes; farthest copper node %.4g mm from its centreline\n", nTan, nTanAir, cuCentrelineMax*1e3);
    std::printf("metric: %ld corner/edge nodes got >1 fine direction (rectangular edge effect)\n", nCorner);
    std::printf("metric: %ld layer nodes, %ld air-side, %ld untouched (%.1fs)\n", nLayer, nOut, nFar, since(t0));
    if (stats) {
        std::printf("stats: fine directions per layer node: 1:%ld 2:%ld 3+:%ld\n", nfineHist[1], nfineHist[2], nfineHist[3]);
        std::printf("stats: angle between the two nearest fine directions (deg, 10-deg bins 0..90):");
        for (int k=0;k<10;k++) std::printf(" %ld", angHist[k]);
        std::printf("\nstats: distance of the 2nd direction in half-bands (%.3f mm):", 0.5*band*1e3);
        for (int k=0;k<8;k++) std::printf(" %ld", d2Hist[k]);
        std::printf("\n");
        if (std::getenv("OMFEM_SKIN_METRIC_ONLY")) { std::printf("metric only, exiting\n"); return 0; }
    }
    MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_verbose, std::getenv("OMFEM_MMG_VERBOSE")?5:-1);
    // MMG MEMORY, derived (2026-09-28; was 16000 MB): what this process can still get -- the machine's
    // MemAvailable, and under an address-space limit (ulimit -v) that limit minus what we already hold.
    // OMFEM_MMG_MEM (MB) overrides.
    int mmgMemMB = 0;
    if (const char* mm = std::getenv("OMFEM_MMG_MEM")) mmgMemMB = std::atoi(mm);
    else {
        double availB = -1.0;
        { std::ifstream mi("/proc/meminfo"); std::string k; double v; std::string unit;
          while (mi >> k >> v >> unit) if (k == "MemAvailable:") { availB = v * 1024.0; break; } }
        if (!(availB > 0)) { std::fprintf(stderr, "ERROR: no MemAvailable in /proc/meminfo; set OMFEM_MMG_MEM (MB)\n"); return 2; }
        rlimit rl{};
        if (getrlimit(RLIMIT_AS, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY) {
            long pages = 0; { std::ifstream sm("/proc/self/statm"); sm >> pages; }
            availB = std::min(availB, (double)rl.rlim_cur - (double)pages * (double)sysconf(_SC_PAGESIZE));
        }
        mmgMemMB = (int)(availB / (1024.0 * 1024.0));
        if (mmgMemMB <= 0) { std::fprintf(stderr, "ERROR: no memory left for MMG (%.0f MB)\n", availB / (1024.0*1024.0)); return 1; }
    }
    std::printf("mmg memory cap: %d MB%s\n", mmgMemMB, std::getenv("OMFEM_MMG_MEM") ? " (OMFEM_MMG_MEM)" : " (available to this process)");
    MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_mem, mmgMemMB);
    // nosurf by default: the copper/air/bobbin/core interfaces and the outer/cut faces stay exactly
    // where the conformal mesher put them (the tangential size there is already the conformal one);
    // only the volume is layered. OMFEM_SKIN_SURF=1 lets MMG remesh surfaces too.
    // SURFACE REMESHING ON by default (2026-09-11). With the surface frozen (nosurf) MMG cannot
    // build the layer: keeping our sizes at the required surface vertices (nosizreq) then refines
    // isotropically everywhere (0.013 mm in all directions, 540k copper tets on the bar), while
    // letting MMG remesh the copper/air interface gives the skin profile that was asked for --
    // measured on the bar for a 0.02 mm first layer: normal extent 0.021 -> 0.033 -> 0.052 ->
    // 0.089 mm with depth, tangential 0.05-0.07 mm, along the wire 0.15-0.18 mm, 126k copper tets.
    // Interfaces move at most hausd off the CAD facets. OMFEM_SKIN_NOSURF=1 freezes them.
    if (std::getenv("OMFEM_SKIN_NOSURF")) MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_nosurf, 1);
    MMG3D_Set_dparameter(mmg, sol, MMG3D_DPARAM_hgrad, hgrad);
    // NO GRADATION FROM THE FROZEN SOLIDS (2026-10-01, ABT #1564). MMG grades the size away from REQUIRED
    // entities with its own hgradreq, independently of hgrad and of the metric we hand it. With the solids
    // frozen that spread their faces' sizes through the whole air and copper: smd_drum copper 153k -> 237k
    // tets, buck adapt 610 -> 1200 s. Off (-1), the frozen core costs nothing: smd_drum 0 slivers, copper
    // 153,909 tets (unfrozen 153,489), core 674 = base. The size near a frozen face is still bounded by
    // hgrad through our metric. OMFEM_MMG_HGRADREQ overrides (A/B only).
    if (const char* gr = std::getenv("OMFEM_MMG_HGRADREQ")) MMG3D_Set_dparameter(mmg, sol, MMG3D_DPARAM_hgradreq, std::atof(gr));
    else if (nFrozen > 0) MMG3D_Set_dparameter(mmg, sol, MMG3D_DPARAM_hgradreq, -1.0);
    // hmax was never set; MMG then derives one from the bounding box and refined the core and air
    // it was told to leave alone (03_buck: 662k -> 1.5M non-copper tets). Cap at the largest local
    // size measured from the input mesh, so nothing outside the copper gets finer than it was.
    { double hmaxIn = 0.0; for (int i=1;i<=np;i++) hmaxIn = std::max(hmaxIn, hloc[i]);
      MMG3D_Set_dparameter(mmg, sol, MMG3D_DPARAM_hmax, std::getenv("OMFEM_MMG_HMAX")?std::atof(std::getenv("OMFEM_MMG_HMAX")):std::max(hmaxIn, hTan)); }
    // COPPER HAUSD, derived (2026-09-28; was 0.01 mm): the copper surface may deviate from its smooth
    // reconstruction by at most the first layer's thickness hs. Tighter asks for a geometric fidelity
    // the layers cannot represent; looser lets the wall move by more than the layer that resolves it.
    const double hausdCu = hs;
    MMG3D_Set_dparameter(mmg, sol, MMG3D_DPARAM_hausd, hausdCu);
    // HMIN, derived (2026-09-28; was 0.5*hs): the finest size anywhere in the metric we hand MMG --
    // the wall layer, or an input element finer than it (a thin gap). MMG's own curvature metric may
    // not go below what was asked, and nothing the conformal mesher resolved is forced coarser.
    // NOT the shortest CAD-edge segment (tried 2026-09-29): on 20_iso that is 0.89 um, on the
    // 1-4 um rings MVB++ leaves in wound turns (ABT #1528), and hmin there gave 49 slivers
    // instead of 33 (long flat tets appeared). Those rings are a geometry defect, not a size to honour.
    const double hmin = hMetricMin;
    MMG3D_Set_dparameter(mmg, sol, MMG3D_DPARAM_hmin, hmin);
    // RIDGES: only the CAD's edges handed over above; no angle detection.
    MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_angle, 0);
    std::printf("mmg: hgrad %.4g, hmin %.4g mm (finest metric size), copper hausd %.4g mm (= h_surface), ridges = %zu CAD edges\n",
                hgrad, hmin*1e3, hausdCu*1e3, cadEdges.size());
    // LOCAL REFINEMENT (2026-09-20). MMG3D_defsiz_ani computes a CURVATURE metric at every
    // boundary point of the model from the Hausdorff tolerance and INTERSECTS it with the metric
    // we hand it (MMG3D_intextmet); the intersection keeps the finer of the two. A single global
    // hausd therefore refines every curved surface in the model -- the core's round centre post,
    // the bobbin, the outer box -- to within hausd of its Bezier reconstruction, whether or not a
    // conductor is anywhere near. Measured on ETD34 (corpus v7): the CORE grew 1,287 -> 12,967
    // tets (10.1x) in a pass whose only job is copper skin layers.
    // MMG's own mechanism for keeping a tolerance local is a per-reference local parameter, and
    // it accepts one per TETRAHEDRON reference (MMG5_Tetrahedron): MMG5_defmetreg /
    // MMG5_defmetvol / mmg3d1 take hausd, hmin and hmax from the tets in the ball of the point,
    // so a boundary point of the core/air interface never sees the copper tolerance.
    // The far-field hausd is DERIVED, not chosen: MMG sets the surface eigenvalue to
    // 1/h^2 = (2/9)*kappa/hausd (MMG5_solveDefmetregSys), i.e. h = sqrt(4.5*hausd*R). Asking that
    // this never fall below the size the region already has, h >= h_reg, for any radius the input
    // mesh could resolve (R >= h_reg), gives hausd >= h_reg/4.5. That is the threshold at which
    // the curvature rule stops adding refinement to a region, and it is what we use.
    // OMFEM_SKIN_LOCALPAR=0 restores the single global tolerance (A/B only).
    const bool localPar = !std::getenv("OMFEM_SKIN_LOCALPAR") || std::atoi(std::getenv("OMFEM_SKIN_LOCALPAR")) != 0;
    if (localPar) {
        const double farK = std::getenv("OMFEM_SKIN_HAUSD_FAR_K") ? std::atof(std::getenv("OMFEM_SKIN_HAUSD_FAR_K")) : 4.5;
        if (!(farK > 0)) { std::fprintf(stderr, "OMFEM_SKIN_HAUSD_FAR_K must be positive\n"); return 2; }
        if (!MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_numberOfLocalParam, (int)volNames.size())) {
            std::fprintf(stderr, "ERROR: MMG3D refused %zu local parameters\n", volNames.size()); return 1; }
        for (size_t i = 0; i < volNames.size(); ++i) {
            const int ref = (int)i + 1;
            if (nRefTet[ref] == 0 || !(hRefMax[ref] > 0)) {     // no silent default: a region with no
                std::fprintf(stderr, "ERROR: physical volume \"%s\" (ref %d) has %ld tets and measured size %g m in %s -- cannot set its local parameters\n",
                             volNames[i].c_str(), ref, nRefTet[ref], hRefMax[ref], argv[1]);
                return 1;                                       // measurable size cannot be given one
            }
            const bool isCu = copperRef[i];
            const double hmaxR  = isCu ? std::max(hCuIn, hTan) : hRefMax[ref];
            const double hausdR = isCu ? hausdCu : hRefMax[ref] / farK;
            if (!MMG3D_Set_localParameter(mmg, sol, MMG5_Tetrahedron, ref, hmin, hmaxR, hausdR)) {
                std::fprintf(stderr, "ERROR: MMG3D_Set_localParameter failed for \"%s\" (ref %d)\n", volNames[i].c_str(), ref); return 1; }
            std::printf("local param: %-14s ref %d  %8ld tets  h_in<=%.3f mm -> hmax %.3f mm, hausd %.4f mm%s\n",
                        volNames[i].c_str(), ref, nRefTet[ref], hRefMax[ref]*1e3, hmaxR*1e3, hausdR*1e3, isCu ? "  (copper)" : "");
        }
    }
    // MMG honoured an ISOTROPIC request exactly (asked 0.03 mm, got 0.029/0.029/0.030) but returned
    // the coarse size in every direction for an anisotropic one (asked 0.02 normal x 0.15
    // tangential, got 0.134/0.146/0.150): the tensor's directions were being flattened. anisosize
    // is documented as "anisotropic metric creation when no metric is provided", but it is the
    // only switch that selects MMG's anisotropic code path. OMFEM_MMG_ANISO=0 restores the old call.
    if (!std::getenv("OMFEM_MMG_ANISO") || std::atoi(std::getenv("OMFEM_MMG_ANISO")) != 0)
        MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_anisosize, 1);
    // With nosurf the copper surface vertices are REQUIRED, and MMG by default overwrites the
    // metric at required vertices with the size they already have (0.15 mm, isotropic). The thin
    // normal size we prescribe AT the surface is therefore thrown away exactly where the layer
    // starts, and gradation from those coarse surface metrics flattens the whole band: measured
    // 0.079 mm normal in the first 0.03 mm and 0.13-0.15 beyond, for a 0.02 mm request. nosizreq
    // keeps our sizes at required vertices. OMFEM_MMG_NOSIZREQ=0 restores the default.
    if (!std::getenv("OMFEM_MMG_NOSIZREQ") || std::atoi(std::getenv("OMFEM_MMG_NOSIZREQ")) != 0)
        MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_nosizreq, 1);
    const int ier = MMG3D_mmg3dlib(mmg, sol);
    const double adapt_s = since(t0);
    if (ier == MMG5_STRONGFAILURE) { std::fprintf(stderr,"MMG3D STRONG FAILURE\n"); return 1; }
    if (ier != MMG5_SUCCESS) std::fprintf(stderr,"MMG3D returned %d (lowfailure, continuing)\n", ier);
    int np2=0,ne2=0,nt2=0; MMG3D_Get_meshSize(mmg, &np2,&ne2,nullptr,&nt2,nullptr,nullptr);
    std::printf("MMG3D: %d->%d nodes, %zu->%d tets, %zu->%d tris  (adapt %.1fs, ier=%d)\n", np, np2, tets.size(), ne2, tris.size(), nt2, adapt_s, ier);

    std::FILE* f = std::fopen(argv[2], "w");
    std::fprintf(f, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
    std::fprintf(f, "$PhysicalNames\n%zu\n", volNames.size() + bndNames.size());
    for (size_t i=0;i<bndNames.size();++i) std::fprintf(f, "2 %d \"%s\"\n", kBndBase + (int)i + 1, bndNames[i].c_str());
    for (size_t i=0;i<volNames.size();++i) std::fprintf(f, "3 %d \"%s\"\n", (int)i + 1, volNames[i].c_str());
    std::fprintf(f, "$EndPhysicalNames\n$Nodes\n%d\n", np2);
    for (int i=1;i<=np2;i++){ double x,y,z; int ref; MMG3D_Get_vertex(mmg,&x,&y,&z,&ref,nullptr,nullptr); std::fprintf(f,"%d %.10g %.10g %.10g\n",i,x,y,z); }
    std::vector<std::array<int,4>> outTris;
    for (int k=1;k<=nt2;k++){ int v0,v1,v2,ref; MMG3D_Get_triangle(mmg,&v0,&v1,&v2,&ref,nullptr);
        if (ref > kBndBase && ref <= kBndBase + (int)bndNames.size()) outTris.push_back({v0,v1,v2,ref}); }
    std::fprintf(f, "$EndNodes\n$Elements\n%zu\n", (size_t)ne2 + outTris.size());
    int eid = 0; std::vector<long> cnt(volNames.size()+1, 0);
    for (auto& t : outTris) std::fprintf(f,"%d 2 2 %d %d %d %d %d\n", ++eid, t[3], t[3], t[0], t[1], t[2]);
    for (int k=1;k<=ne2;k++){ int v0,v1,v2,v3,ref; MMG3D_Get_tetrahedron(mmg,&v0,&v1,&v2,&v3,&ref,nullptr);
        if (ref>=1 && ref<=(int)volNames.size()) cnt[ref]++;
        std::fprintf(f,"%d 4 2 %d %d %d %d %d %d\n", ++eid, ref, ref, v0,v1,v2,v3); }
    std::fprintf(f, "$EndElements\n"); std::fclose(f);
    long cu=0, other=0; for (size_t i=0;i<volNames.size();i++) (copperRef[i]?cu:other) += cnt[i+1];
    {   // effective layer recipe next to the output (the values this run used, positional or recipe)
        setenv("OMFEM_LAYERS_ENABLED", "1", 1);
        char b[64]; std::snprintf(b, sizeof b, "%.12g", delta); setenv("OMFEM_LAYERS_DELTA", b, 1);
        std::snprintf(b, sizeof b, "%.12g", cells); setenv("OMFEM_LAYERS_CELLS", b, 1);
        // the layer count this run used (derived or asked for); h_deep and hgrad are inputs already in the environment
        std::snprintf(b, sizeof b, "%d", nLayers); setenv("OMFEM_SKIN_LAYERS", b, 1);
        std::snprintf(b, sizeof b, "%.12g", hTan);  setenv("OMFEM_LAYERS_HTAN", b, 1);
        nlohmann::json prov; prov["tool"] = "mvbpp_skinlayer"; prov["input"] = argv[1]; prov["mesh"] = argv[2];
        prov["recipe_in"] = recipePath.empty() ? nlohmann::json(nullptr) : nlohmann::json(recipePath);
        prov["cwd"] = std::filesystem::current_path().string();
        prov["derived"] = {{"n_layers", nLayers}, {"growth", growthK}, {"band_m", band}, {"reach_m", reach},
                           {"feature_angle_deg", featureDeg}, {"cad_edges", cadEdges.size()}, {"segments", segments}, {"hmin_m", hmin}, {"hausd_copper_m", hausdCu},
                           {"mmg_memory_mb", mmgMemMB}};
        // the recipe is how this mesh is reproduced: failing to write it fails the run
        try { mvb::mesh::write_mesh_recipe(std::string(argv[2]) + ".recipe.json", mvb::mesh::effective_mesh_recipe({"layers"}, prov)); }
        catch (const std::exception& e) { std::fprintf(stderr, "ERROR: %s\n", e.what()); return 1; }
    }
    std::printf("wrote %s: copper %ld tets, other %ld tets; boundary tris kept: ", argv[2], cu, other);
    for (size_t i=0;i<bndNames.size();++i){ size_t c=0; for(auto&t:outTris) if(t[3]==kBndBase+(int)i+1)++c; std::printf("%s=%zu ", bndNames[i].c_str(), c); }
    std::printf("\n");
    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh,&mmg, MMG5_ARG_ppMet,&sol, MMG5_ARG_end);
    return 0;
}
