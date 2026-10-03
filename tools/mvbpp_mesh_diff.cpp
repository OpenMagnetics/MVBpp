// mvbpp_mesh_diff <a.msh> <b.msh>: the identity gate as a command (ABT #1588). Both files are
// read as msh2 (gmsh's, the skin layer's or exportMesh's) and compared with diffMeshes: node ids
// and exact coordinates, regions, and every block's type, region, entity, element ids and node
// lists, in order. Prints "identical" (exit 0) or the first difference (exit 1).
#include "mvb/mesh/Mesher.h"

#include <fstream>
#include <iostream>
#include <sstream>

static std::string slurp(const char* path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error(std::string("cannot read ") + path);
    std::stringstream s;
    s << f.rdbuf();
    return s.str();
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: mvbpp_mesh_diff <a.msh> <b.msh>\n";
        return 2;
    }
    try {
        const auto d = mvb::mesh::diffMeshes(mvb::mesh::importMesh(slurp(argv[1]), "msh2"),
                                             mvb::mesh::importMesh(slurp(argv[2]), "msh2"));
        std::cout << (d.empty() ? "identical" : "DIFF " + d) << "\n";
        return d.empty() ? 0 : 1;
    } catch (const std::exception& e) {
        std::cout << "FAILED " << e.what() << "\n";
        return 2;
    }
}
