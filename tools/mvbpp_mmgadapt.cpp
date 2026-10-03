// mvbpp_mmgadapt <in.msh> <out.msh> [axis Y|Z|X] [across_mm] [along_mm] [coarse_mm]
// Moved verbatim from OMFEM tools/omfem_mmgadapt.cpp at b536a9a (ABT #1588, step 7); only the tool name changed.
//
// Standalone MMG3D anisotropic REMESH (the proper MMG workflow gmsh's integrated path
// can't do): take a VALID isotropic tet mesh, attach an anisotropic tensor metric -- fine
// ACROSS each wire, coarse ALONG it (azimuthal about the column axis), in the winding;
// coarse isotropic elsewhere -- and call MMG3D_mmg3dlib to coarsen along the wires.
// Reports element counts + runtime. Pure C++: gmsh reads/writes, libmmg adapts.
#include <gmsh.h>
extern "C" {
#include <mmg/mmg3d/libmmg3d.h>
}
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>

using clk = std::chrono::steady_clock;
static double since(clk::time_point t){ return std::chrono::duration<double>(clk::now()-t).count(); }

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr,"usage: %s in.msh out.msh [axis] [across_mm] [along_mm] [coarse_mm]\n",argv[0]); return 2; }
    const char  axis    = (argc>3)? argv[3][0] : 'Y';
    const double h_across = (argc>4? std::atof(argv[4]) : 0.20) * 1e-3;
    const double h_along  = (argc>5? std::atof(argv[5]) : 2.0)  * 1e-3;
    const double h_coarse = (argc>6? std::atof(argv[6]) : 4.0)  * 1e-3;

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::open(argv[1]);

    // --- nodes (contiguous 1-based for MMG) ---
    std::vector<std::size_t> nodeTags; std::vector<double> coord, par;
    gmsh::model::mesh::getNodes(nodeTags, coord, par, -1, -1, false, false);
    const int np = (int)nodeTags.size();
    std::map<std::size_t,int> id; for (int i=0;i<np;i++) id[nodeTags[i]] = i+1;
    std::vector<double> X(np), Y(np), Z(np);
    for (int i=0;i<np;i++){ X[i]=coord[3*i]; Y[i]=coord[3*i+1]; Z[i]=coord[3*i+2]; }

    // --- tets per volume physical group. EVERY named group keeps its own ref, so the adapted
    // mesh carries core/bobbin/winding_*/turn_*/air etc. through unchanged -- the original
    // core=1/winding=2/air=3 collapse silently relabelled the bobbin as air (wrong thermal
    // conductivity downstream) and merged distinct windings.
    struct Tet { int v[4]; int ref; };
    std::vector<Tet> tets; std::set<int> windNode;
    std::vector<std::string> volNames;                    // ref-1 -> name
    gmsh::vectorpair pgs; gmsh::model::getPhysicalGroups(pgs, 3);
    for (auto& pg : pgs) {
        std::string nm; gmsh::model::getPhysicalName(3, pg.second, nm);
        volNames.push_back(nm);
        const int ref = (int)volNames.size();
        const bool isWind = nm.rfind("winding",0)==0 || nm.rfind("turn",0)==0;
        std::vector<int> ents; gmsh::model::getEntitiesForPhysicalGroup(3, pg.second, ents);
        for (int e : ents) {
            std::vector<int> ety; std::vector<std::vector<std::size_t>> etag, enod;
            gmsh::model::mesh::getElements(ety, etag, enod, 3, e);
            for (size_t t=0;t<ety.size();t++) if (ety[t]==4) {
                const auto& nn = enod[t];
                for (size_t k=0;k+3<nn.size();k+=4) {
                    Tet tt; for (int j=0;j<4;j++) tt.v[j]=id[nn[k+j]]; tt.ref=ref; tets.push_back(tt);
                    if (isWind) for (int j=0;j<4;j++) windNode.insert(tt.v[j]);
                }
            }
        }
    }
    // --- boundary triangles from PHYSICAL surface groups (term0/term1/outer/...). Refs live in
    // a dedicated band (kBndBase+idx) so they can never collide with the volume refs; MMG keeps
    // the ref of every triangle patch it preserves/rebuilds, which is exactly what the solvers
    // need (term0/term1 are the eddy/RF drive electrodes, outer is the far-field BC).
    constexpr int kBndBase = 1000;
    struct Tri { int v[3]; int ref; };
    std::vector<Tri> tris;
    std::vector<std::string> bndNames;                    // idx -> name
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
    std::printf("loaded: %d nodes, %zu tets (%zu vol groups), %zu boundary tris (%zu bnd groups), %zu winding nodes\n",
                np, tets.size(), volNames.size(), tris.size(), bndNames.size(), windNode.size());

    // --- build MMG mesh + anisotropic tensor metric ---
    auto t0 = clk::now();
    MMG5_pMesh mmg=nullptr; MMG5_pSol sol=nullptr;
    MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh,&mmg, MMG5_ARG_ppMet,&sol, MMG5_ARG_end);
    MMG3D_Set_meshSize(mmg, np, (int)tets.size(), 0, (int)tris.size(), 0, 0);
    for (int i=0;i<np;i++) MMG3D_Set_vertex(mmg, X[i],Y[i],Z[i], 0, i+1);
    for (size_t k=0;k<tets.size();k++) MMG3D_Set_tetrahedron(mmg, tets[k].v[0],tets[k].v[1],tets[k].v[2],tets[k].v[3], tets[k].ref, (int)k+1);
    for (size_t k=0;k<tris.size();k++) MMG3D_Set_triangle(mmg, tris[k].v[0],tris[k].v[1],tris[k].v[2], tris[k].ref, (int)k+1);

    MMG3D_Set_solSize(mmg, sol, MMG5_Vertex, np, MMG5_Tensor);
    const double a = 1.0/(h_across*h_across);
    const double b = 1.0/(h_along*h_along) - a;      // <0: coarsen azimuthal
    const double c = 1.0/(h_coarse*h_coarse);        // isotropic outside winding
    for (int i=0;i<np;i++) {
        double m11,m12,m13,m22,m23,m33;
        if (windNode.count(i+1)) {
            double u,w;                              // coords in the plane perpendicular to axis
            if      (axis=='Y'){ u=X[i]; w=Z[i]; }   // e_phi about Y -> (-z,0,x)
            else if (axis=='Z'){ u=X[i]; w=Y[i]; }
            else               { u=Y[i]; w=Z[i]; }
            double r2 = u*u + w*w + 1e-30;
            // M = a*I + b*(ephi (x) ephi), ephi in (u,w) plane = (-w,u)/r
            double exx = b*w*w/r2, eww = b*u*u/r2, exw = -b*u*w/r2;
            if (axis=='Y'){ m11=a+exx; m33=a+eww; m13=exw; m22=a; m12=0; m23=0; }
            else if (axis=='Z'){ m11=a+exx; m22=a+eww; m12=exw; m33=a; m13=0; m23=0; }
            else { m22=a+exx; m33=a+eww; m23=exw; m11=a; m12=0; m13=0; }
        } else { m11=m22=m33=c; m12=m13=m23=0; }
        MMG3D_Set_tensorSol(sol, m11,m12,m13,m22,m23,m33, i+1);
    }
    // options: be quiet, allow size gradation jumps, preserve geometry tightly.
    MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_verbose, std::getenv("OMFEM_MMG_VERBOSE")?5:-1);
    MMG3D_Set_dparameter(mmg, sol, MMG3D_DPARAM_hgrad, std::getenv("OMFEM_MMG_HGRAD")?std::atof(std::getenv("OMFEM_MMG_HGRAD")):3.0);
    MMG3D_Set_dparameter(mmg, sol, MMG3D_DPARAM_hausd, std::getenv("OMFEM_MMG_HAUSD")?std::atof(std::getenv("OMFEM_MMG_HAUSD")):0.05e-3);

    const int ier = MMG3D_mmg3dlib(mmg, sol);
    const double adapt_s = since(t0);
    if (ier == MMG5_STRONGFAILURE) { std::fprintf(stderr,"MMG3D STRONG FAILURE\n"); return 1; }
    if (ier != MMG5_SUCCESS) std::fprintf(stderr,"MMG3D returned %d (lowfailure, continuing)\n", ier);

    // --- read back ---
    int np2=0,ne2=0,nt2=0; MMG3D_Get_meshSize(mmg, &np2,&ne2,nullptr,&nt2,nullptr,nullptr);
    std::printf("MMG3D: %d->%d nodes, %zu->%d tets, %zu->%d tris  (adapt %.1fs, ier=%d)\n",
                np, np2, tets.size(), ne2, tris.size(), nt2, adapt_s, ier);

    // --- write gmsh v2.2 with EVERY volume and boundary physical group preserved ---
    std::FILE* f = std::fopen(argv[2], "w");
    std::fprintf(f, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
    std::fprintf(f, "$PhysicalNames\n%zu\n", volNames.size() + bndNames.size());
    for (size_t i=0;i<bndNames.size();++i)
        std::fprintf(f, "2 %d \"%s\"\n", kBndBase + (int)i + 1, bndNames[i].c_str());
    for (size_t i=0;i<volNames.size();++i)
        std::fprintf(f, "3 %d \"%s\"\n", (int)i + 1, volNames[i].c_str());
    std::fprintf(f, "$EndPhysicalNames\n");
    std::fprintf(f, "$Nodes\n%d\n", np2);
    for (int i=1;i<=np2;i++){ double x,y,z; int ref; MMG3D_Get_vertex(mmg,&x,&y,&z,&ref,nullptr,nullptr); std::fprintf(f,"%d %.10g %.10g %.10g\n",i,x,y,z); }
    // Count the output triangles that carry a boundary-band ref (MMG can also emit interface
    // triangles it created itself with tet-derived refs; those are not ours and are skipped).
    std::vector<std::array<int,4>> outTris;
    for (int k=1;k<=nt2;k++){ int v0,v1,v2,ref; MMG3D_Get_triangle(mmg,&v0,&v1,&v2,&ref,nullptr);
        if (ref > kBndBase && ref <= kBndBase + (int)bndNames.size())
            outTris.push_back({v0,v1,v2,ref}); }
    std::fprintf(f, "$EndNodes\n$Elements\n%zu\n", (size_t)ne2 + outTris.size());
    int eid = 0;
    for (auto& t : outTris)
        std::fprintf(f,"%d 2 2 %d %d %d %d %d\n", ++eid, t[3], t[3], t[0], t[1], t[2]);
    for (int k=1;k<=ne2;k++){ int v0,v1,v2,v3,ref; MMG3D_Get_tetrahedron(mmg,&v0,&v1,&v2,&v3,&ref,nullptr);
        std::fprintf(f,"%d 4 2 %d %d %d %d %d %d\n", ++eid, ref, ref, v0,v1,v2,v3); }
    std::fprintf(f, "$EndElements\n");
    std::fclose(f);
    std::printf("wrote %s  (%zu boundary tris kept: ", argv[2], outTris.size());
    for (size_t i=0;i<bndNames.size();++i){ size_t c=0; for(auto&t:outTris) if(t[3]==kBndBase+(int)i+1)++c;
        std::printf("%s=%zu ", bndNames[i].c_str(), c); }
    std::printf(")\n");

    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh,&mmg, MMG5_ARG_ppMet,&sol, MMG5_ARG_end);
    gmsh::finalize();
    return 0;
}
