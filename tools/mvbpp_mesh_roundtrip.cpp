// mvbpp_mesh_roundtrip <mesh.msh>...: the round-trip gate on real meshes (ABT #1588).
// Each file (an msh2 that gmsh or OMFEM's skin layer wrote) is read into the mesh model, written
// back in every format, read again and compared: msh2 exactly, the BDF in metres and millimetres
// to its 16-character precision (11-13 significant digits, so 1e-11 relative), with every tet's volume keeping its sign and no new
// |V| < 1e-18 m^3. A material is not known from an msh2, so every region is stated as having
// none (the BDF then carries properties and no MAT cards). Exits 1 on the first failure.
#include "mvb/mesh/Mesher.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

using namespace mvb::mesh;

namespace {

std::string slurp(const char* path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error(std::string("cannot read ") + path);
    std::stringstream s;
    s << f.rdbuf();
    return s.str();
}

std::vector<double> tetVolumes(const Mesh& m) {
    std::map<std::int64_t, std::size_t> at;
    for (std::size_t i = 0; i < m.nodeIds.size(); ++i) at[m.nodeIds[i]] = i;
    std::vector<double> v;
    for (const auto& b : m.blocks) {
        if (b.type != ElementType::Tet4) continue;
        for (std::size_t e = 0; e < b.ids.size(); ++e) {
            const double* p[4];
            for (int k = 0; k < 4; ++k) p[k] = &m.xyz[3 * at.at(b.nodes[4 * e + k])];
            double u[3], w[3], z[3];
            for (int k = 0; k < 3; ++k) { u[k] = p[1][k] - p[0][k]; w[k] = p[2][k] - p[0][k]; z[k] = p[3][k] - p[0][k]; }
            v.push_back((u[0] * (w[1] * z[2] - w[2] * z[1]) - u[1] * (w[0] * z[2] - w[2] * z[0]) +
                         u[2] * (w[0] * z[1] - w[1] * z[0])) / 6.0);
        }
    }
    return v;
}

// A BDF has no elementary entity: drop it and join the blocks that differed only by it, as the
// reader does, so the rest of the mesh compares exactly.
Mesh withoutEntities(const Mesh& m) {
    Mesh out = m;
    out.blocks.clear();
    for (auto b : m.blocks) {
        b.entity = 0;
        if (!out.blocks.empty() && out.blocks.back().type == b.type && out.blocks.back().region == b.region) {
            auto& last = out.blocks.back();
            last.ids.insert(last.ids.end(), b.ids.begin(), b.ids.end());
            last.nodes.insert(last.nodes.end(), b.nodes.begin(), b.nodes.end());
        } else {
            out.blocks.push_back(b);
        }
    }
    return out;
}

std::string checkBdf(const Mesh& m, const std::string& unit) {
    const Mesh back = importMesh(exportMesh(m, "bdf", unit), "bdf");
    if (back.nodeIds != m.nodeIds) return "node ids differ";
    double worst = 0;
    for (std::size_t i = 0; i < m.xyz.size(); ++i) {
        const double rel = std::abs(m.xyz[i] - back.xyz[i]) / std::max(std::abs(m.xyz[i]), 1e-9);
        worst = std::max(worst, rel);
    }
    if (worst > 1e-11) return "coordinate relative error " + std::to_string(worst);
    Mesh a = withoutEntities(m), b = back;
    b.xyz = a.xyz;
    const auto d = diffMeshes(a, b);
    if (!d.empty()) return d;
    const auto v0 = tetVolumes(m), v1 = tetVolumes(back);
    std::size_t flipped = 0, newTiny = 0;
    for (std::size_t i = 0; i < v0.size(); ++i) {
        if ((v0[i] > 0) != (v1[i] > 0)) ++flipped;
        if (std::abs(v1[i]) < 1e-18 && std::abs(v0[i]) >= 1e-18) ++newTiny;
    }
    if (flipped || newTiny)
        return std::to_string(flipped) + " tets changed sign, " + std::to_string(newTiny) + " became |V| < 1e-18";
    std::ostringstream ok;
    ok << "ok (worst relative coordinate error " << worst << ", " << v0.size() << " tets)";
    return ok.str();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: mvbpp_mesh_roundtrip <mesh.msh>...\n";
        return 2;
    }
    for (int i = 1; i < argc; ++i) {
        try {
            Mesh m = importMesh(slurp(argv[i]), "msh2");
            for (auto& r : m.regions) r.material = {MaterialKind::None, {}};
            std::size_t elements = 0;
            for (const auto& b : m.blocks) elements += b.ids.size();
            std::cout << argv[i] << ": " << m.nodeIds.size() << " nodes, " << elements << " elements, "
                      << m.regions.size() << " regions, " << m.blocks.size() << " blocks\n";
            const auto d = diffMeshes(m, importMesh(exportMesh(m, "msh2", "m"), "msh2"));
            std::cout << "  msh2: " << (d.empty() ? "identical" : "DIFF " + d) << "\n";
            if (!d.empty()) return 1;
            for (const char* unit : {"m", "mm"}) {
                const auto r = checkBdf(m, unit);
                std::cout << "  bdf " << unit << ": " << r << "\n";
                if (r.rfind("ok", 0) != 0) return 1;
            }
        } catch (const std::exception& e) {
            std::cout << argv[i] << ": FAILED " << e.what() << "\n";
            return 1;
        }
    }
    return 0;
}
