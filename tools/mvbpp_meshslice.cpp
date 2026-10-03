// mvbpp_meshslice <mesh.msh> <out.svg> [plane=z0|y0|x0] [coord_mm]
// Moved verbatim from OMFEM tools/omfem_meshslice.cpp at 070ae86 (ABT #1588, step 7); only the tool name changed.
// Cuts a 3D tet mesh by a coordinate plane and writes an SVG cross-section, colored by
// physical region (core grey / winding copper-red / air faint). Pure C++ (gmsh reads the
// mesh; we compute the plane-tet intersection polygons + emit SVG). No Python, no graphics.
#include <gmsh.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr,"usage: %s mesh.msh out.svg [z0|y0|x0] [coord_mm]\n",argv[0]); return 2; }
    const std::string plane = (argc>3)? argv[3] : "z0";
    const double cut = (argc>4? std::atof(argv[4]) : 0.0) * 1e-3;
    // axis index of the cut normal, and the two in-plane axes (we project those to SVG x,y).
    int n = plane[0]=='x'?0 : plane[0]=='y'?1 : 2;
    int u = (n==0)?1:0, v = (n==2)?1:2;

    gmsh::initialize(); gmsh::option::setNumber("General.Terminal",0); gmsh::open(argv[1]);
    std::vector<std::size_t> nt; std::vector<double> nc, par;
    gmsh::model::mesh::getNodes(nt, nc, par, -1,-1,false,false);
    std::map<std::size_t,int> id; for (size_t i=0;i<nt.size();++i) id[nt[i]]=(int)i;
    auto C=[&](int j,int a){ return nc[3*j+a]; };

    struct Poly { std::vector<std::pair<double,double>> p; int reg; };
    std::vector<Poly> polys;
    double U0=1e30,U1=-1e30,V0=1e30,V1=-1e30;
    const char* names[3]={"core","winding","air"};
    gmsh::vectorpair pgs; gmsh::model::getPhysicalGroups(pgs,3);
    for (auto& pg : pgs) {
        std::string nm; gmsh::model::getPhysicalName(3,pg.second,nm);
        int reg = nm.rfind("core",0)==0?0 : nm.rfind("winding",0)==0?1 : 2;
        std::vector<int> ents; gmsh::model::getEntitiesForPhysicalGroup(3,pg.second,ents);
        for (int e : ents) {
            std::vector<int> ty; std::vector<std::vector<std::size_t>> tg, no;
            gmsh::model::mesh::getElements(ty,tg,no,3,e);
            for (size_t t=0;t<ty.size();++t) if (ty[t]==4)
                for (size_t k=0;k+3<no[t].size();k+=4) {
                    int vtx[4]; double sd[4];
                    for (int j=0;j<4;j++){ vtx[j]=id[no[t][k+j]]; sd[j]=C(vtx[j],n)-cut; }
                    // edge crossings of the cut plane
                    static const int E[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
                    std::vector<std::pair<double,double>> pts;
                    for (auto& ed : E) {
                        double a=sd[ed[0]], b=sd[ed[1]];
                        if ((a<0)!=(b<0) && a!=b) { double s=a/(a-b);
                            int A=vtx[ed[0]],B=vtx[ed[1]];
                            double pu=C(A,u)+s*(C(B,u)-C(A,u)), pv=C(A,v)+s*(C(B,v)-C(A,v));
                            pts.push_back({pu,pv}); }
                    }
                    if (pts.size()<3) continue;
                    // order the (3 or 4) points by angle about their centroid
                    double cu=0,cv=0; for(auto&q:pts){cu+=q.first;cv+=q.second;} cu/=pts.size();cv/=pts.size();
                    std::sort(pts.begin(),pts.end(),[&](auto&A,auto&B){return std::atan2(A.second-cv,A.first-cu)<std::atan2(B.second-cv,B.first-cu);});
                    for(auto&q:pts){U0=std::min(U0,q.first);U1=std::max(U1,q.first);V0=std::min(V0,q.second);V1=std::max(V1,q.second);}
                    polys.push_back({pts,reg});
                }
        }
    }
    gmsh::finalize();
    if (polys.empty()) { std::fprintf(stderr,"no intersections on plane %s @ %.3gmm\n",plane.c_str(),cut*1e3); return 1; }

    // emit SVG (flip V so +up; scale to ~1000px wide)
    const double W=U1-U0, H=V1-V0, sc=1000.0/W, pad=10;
    std::FILE* f=std::fopen(argv[2],"w");
    std::fprintf(f,"<svg xmlns='http://www.w3.org/2000/svg' width='%.0f' height='%.0f' style='background:white'>\n",W*sc+2*pad,H*sc+2*pad);
    const char* fill[3]={"#9aa0a6","#c0392b","#eef4fb"};   // core grey, winding copper-red, air faint
    const char* strk[3]={"#5f6368","#7b241c","#cdd9e8"};
    for (int pass=2; pass>=0; --pass)  // draw air first, then core, then winding on top
      for (auto& po : polys) { if (po.reg!=pass) continue;
        std::fprintf(f,"<polygon points='");
        for (auto& q : po.p) std::fprintf(f,"%.2f,%.2f ",(q.first-U0)*sc+pad,(V1-q.second)*sc+pad);
        std::fprintf(f,"' fill='%s' stroke='%s' stroke-width='0.25'/>\n",fill[po.reg],strk[po.reg]);
      }
    std::fprintf(f,"</svg>\n"); std::fclose(f);
    std::printf("wrote %s : plane %s @ %.3gmm, %zu cross-section polygons (%.0fx%.0f mm)\n",
                argv[2],plane.c_str(),cut*1e3,polys.size(),W*1e3,H*1e3);
    return 0;
}
