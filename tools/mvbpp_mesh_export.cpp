// mvbpp_mesh_export <mas.json> <out-base> [--format bdf|msh2|inp|vtk] [--unit m|mm] [--parts-only]
//                   [--recipe <recipe.json>] [--step]
//
// Meshes a MAS file's magnetic (mvb::mesh::meshMagnetic, ABT #1588) and writes it in one format:
// <out-base>.<format> ("bdf" is Nastran bulk data, written as .nas). --parts-only leaves the air
// out (mvb::mesh::partsOnly: the parts by name, their terminals and skins). --step also writes
// <out-base>.step, the geometry exactly as it was meshed (millimetres). The recipe defaults to
// {} (every knob at its documented default); the effective one is written as
// <out-base>.recipe.json. A BDF takes its reference temperature from the MAS operating points'
// ambient, and comes with <out-base>.materials.json, its MAS material companion.
#include "mvb/mesh/Mesher.h"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace mvb::mesh;

namespace {

std::string slurp(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read " + path);
    std::stringstream s;
    s << f.rdbuf();
    return s.str();
}

void spit(const std::string& path, const std::string& bytes) {
    std::ofstream f(path, std::ios::binary);
    if (!(f << bytes)) throw std::runtime_error("cannot write " + path);
    std::cout << "wrote " << path << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: mvbpp_mesh_export <mas.json> <out-base> [--format bdf|msh2|inp|vtk] [--unit m|mm] "
                     "[--parts-only] [--recipe <recipe.json>] [--step]\n";
        return 2;
    }
    try {
        const std::string base = argv[2];
        std::string format = "bdf", unit = "mm";
        bool parts = false, step = false;
        nlohmann::json recipe = nlohmann::json::object();
        for (int i = 3; i < argc; ++i) {
            const std::string a = argv[i];
            auto value = [&]() -> std::string {
                if (i + 1 >= argc) throw std::runtime_error(a + " needs a value");
                return argv[++i];
            };
            if (a == "--format") format = value();
            else if (a == "--unit") unit = value();
            else if (a == "--parts-only") parts = true;
            else if (a == "--step") step = true;
            else if (a == "--recipe") recipe = nlohmann::json::parse(slurp(value()));
            else throw std::runtime_error("unknown argument " + a);
        }
        const auto mas = nlohmann::json::parse(slurp(argv[1]));
        const auto r = meshMagnetic(mas.at("magnetic"), recipe, step);
        const Mesh mesh = parts ? partsOnly(r.mesh) : r.mesh;
        std::size_t elements = 0;
        for (const auto& b : mesh.blocks) elements += b.ids.size();
        std::cout << mesh.nodeIds.size() << " nodes, " << elements << " elements; regions:";
        for (const auto& reg : mesh.regions) std::cout << " " << reg.name;
        std::cout << "\n";
        nlohmann::json options = nlohmann::json::object();
        if (format == "bdf") options["ambientTemperature"] = ambientTemperature(mas);
        const std::string ext = format == "bdf" ? "nas" : format == "msh2" ? "msh" : format;
        spit(base + "." + ext, exportMesh(mesh, format, unit, options));
        if (format == "bdf") spit(base + ".materials.json", exportMaterials(mesh).dump(2) + "\n");
        spit(base + ".recipe.json", r.sidecars.at(".recipe.json"));
        if (step) spit(base + ".step", r.sidecars.at(".step"));
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAILED " << e.what() << "\n";
        return 1;
    }
}
