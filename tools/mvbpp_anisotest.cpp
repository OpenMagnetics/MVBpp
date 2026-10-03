// mvbpp_anisotest  --  capability probe: can this gmsh do ANISOTROPIC 3D meshing?
// Moved verbatim from OMFEM tools/omfem_anisotest.cpp at df21432 (ABT #1588, step 8); renamed. No logic changed.
// Meshes a 10mm box with a metric asking for fine (hx) in x, coarse (hy=hz) in y,z.
// If anisotropic works: ~ (L/hx)*(L/hy)*(L/hz) elements (few hundred, stretched).
// If it falls back to isotropic min size: ~ (L/hx)^3 (tens of thousands, regular).
// This decides whether the MAS-metric winding-mesh plan is viable on this gmsh.
#include <gmsh.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

int main(int argc, char** argv) {
    const double L  = 10e-3;
    const double hx = (argc > 1 ? std::atof(argv[1]) : 0.3) * 1e-3;   // fine
    const double hy = (argc > 2 ? std::atof(argv[2]) : 3.0) * 1e-3;   // coarse
    const bool scalar = std::getenv("OMFEM_SCALAR");
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("probe");
    gmsh::model::occ::addBox(0,0,0, L,L,L);
    gmsh::model::occ::synchronize();

    if (scalar) {
        // known-good scalar background field: does the mechanism refine at all?
        const int fs = gmsh::model::mesh::field::add("MathEval");
        char e[64]; std::snprintf(e, sizeof e, "%.8g", hx);
        gmsh::model::mesh::field::setString(fs, "F", e);
        gmsh::model::mesh::field::setAsBackgroundMesh(fs);
        gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
        gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
        gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
        gmsh::model::mesh::generate(3);
        std::vector<int> et; std::vector<std::vector<std::size_t>> en, ev;
        gmsh::model::mesh::getElements(et, en, ev, 3);
        std::size_t nt=0; for (size_t i=0;i<et.size();++i) if (et[i]==4) nt=en[i].size();
        std::printf("SCALAR MathEval h=%.2fmm -> tets=%zu (expect ~%.0f if honored)\n",
                    hx*1e3, nt, std::pow(L/hx,3)*6);
        gmsh::finalize(); return 0;
    }

    // MathEvalAniso metric. Convention test: OMFEM_METRIC=size -> diag(hx,hy,hz) [size tensor],
    // else diag(1/hx^2,1/hy^2,1/hz^2) [inverse-square metric].
    const bool sizeConv = std::getenv("OMFEM_METRIC") && std::string(std::getenv("OMFEM_METRIC"))=="size";
    const int f = gmsh::model::mesh::field::add("MathEvalAniso");
    char m11[64], m22[64], m33[64];
    std::snprintf(m11, sizeof m11, "%.8g", sizeConv ? hx : 1.0/(hx*hx));
    std::snprintf(m22, sizeof m22, "%.8g", sizeConv ? hy : 1.0/(hy*hy));
    std::snprintf(m33, sizeof m33, "%.8g", sizeConv ? hy : 1.0/(hy*hy));
    gmsh::model::mesh::field::setString(f, "m11", m11);
    gmsh::model::mesh::field::setString(f, "m22", m22);
    gmsh::model::mesh::field::setString(f, "m33", m33);
    gmsh::model::mesh::field::setString(f, "m12", "0");
    gmsh::model::mesh::field::setString(f, "m13", "0");
    gmsh::model::mesh::field::setString(f, "m23", "0");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);
    if (const char* a = std::getenv("OMFEM_ALGO3D")) gmsh::option::setNumber("Mesh.Algorithm3D", std::atoi(a));
    if (std::getenv("OMFEM_ANISOMAX")) gmsh::option::setNumber("Mesh.AnisoMax", 1e30);
    gmsh::option::setNumber("Mesh.Algorithm", 7);   // BAMG: 2D ANISOTROPIC surface mesher (needed to feed MMG3D)
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);

    std::size_t ntet = 0; const char* st = "OK";
    try {
        gmsh::model::mesh::generate(3);
        std::vector<int> et; std::vector<std::vector<std::size_t>> en, ev;
        gmsh::model::mesh::getElements(et, en, ev, 3);
        for (size_t i=0;i<et.size();++i) if (et[i]==4) ntet = en[i].size();
    } catch (const std::exception& e) { st = "FAIL"; std::fprintf(stderr, "%s\n", e.what()); }

    const double iso  = std::pow(L/hx, 3) * 6;          // rough isotropic-fine count
    const double aniso= (L/hx)*(L/hy)*(L/hy) * 6;        // rough anisotropic count
    std::printf("hx=%.2fmm hy=hz=%.2fmm -> %s tets=%zu  (iso~%.0f, aniso~%.0f)\n",
                hx*1e3, hy*1e3, st, ntet, iso, aniso);
    std::printf(ntet > 0 && (double)ntet < 0.3*iso ? "ANISOTROPIC WORKS\n" : "looks isotropic / failed\n");
    gmsh::finalize();
    return 0;
}
