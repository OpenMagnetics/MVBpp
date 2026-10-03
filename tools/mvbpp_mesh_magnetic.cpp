// mvbpp_mesh_magnetic <mas.json> <reference.msh> [--write <out.msh>]
//
// The identity gate for meshMagnetic (ABT #1588): meshes the MAS file's magnetic with the recipe
// the reference run recorded (<reference.msh>.recipe.json, its provenance and documentation keys
// dropped) and compares with what that run wrote: the mesh parsed and exact (diffMeshes), the
// .cadedges/.path/.leads.json sidecars byte for byte, and the effective recipe apart from its
// provenance. Exits 1 on any difference. The process environment must not set OMFEM_*/MVB_*:
// every knob comes from the recipe.
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

// An effective recipe back into an input recipe: values only (null = default), no units,
// defaults or provenance.
nlohmann::json inputRecipe(const nlohmann::json& effective) {
    nlohmann::json in = nlohmann::json::object();
    for (auto it = effective.begin(); it != effective.end(); ++it) {
        if (it.key() == "provenance" || it.key() == "recipe_version") continue;
        if (!it.value().is_object()) throw std::runtime_error("recipe section " + it.key() + " is not an object");
        for (auto jt = it.value().begin(); jt != it.value().end(); ++jt) {
            const std::string k = jt.key();
            if (k.rfind("_unit_", 0) == 0) continue;
            if (k.size() > 8 && k.compare(k.size() - 8, 8, "_default") == 0) continue;
            in[it.key()][k] = jt.value();
        }
    }
    return in;
}

nlohmann::json withoutProvenance(nlohmann::json r) {
    r.erase("provenance");
    return r;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3 && !(argc == 5 && std::string(argv[3]) == "--write")) {
        std::cerr << "usage: mvbpp_mesh_magnetic <mas.json> <reference.msh> [--write <out.msh>]\n";
        return 2;
    }
    try {
        const std::string ref = argv[2];
        const auto mas = nlohmann::json::parse(slurp(argv[1]));
        const auto refRecipe = nlohmann::json::parse(slurp(ref + ".recipe.json"));
        const auto result = meshMagnetic(mas.at("magnetic"), inputRecipe(refRecipe));
        if (argc == 5) {
            std::ofstream(argv[4], std::ios::binary) << exportMesh(result.mesh, "msh2", "m");
            for (const auto& [suffix, bytes] : result.sidecars)
                std::ofstream(std::string(argv[4]) + suffix, std::ios::binary) << bytes;
        }
        int failures = 0;
        const auto d = diffMeshes(importMesh(slurp(ref), "msh2"), result.mesh);
        std::cout << "mesh: " << (d.empty() ? "identical" : "DIFF " + d) << "\n";
        failures += !d.empty();
        for (const char* suffix : {".cadedges", ".path", ".leads.json"}) {
            const bool same = slurp(ref + suffix) == result.sidecars.at(suffix);
            std::cout << suffix << ": " << (same ? "identical" : "DIFF") << "\n";
            failures += !same;
        }
        const bool recipeSame = withoutProvenance(refRecipe) == withoutProvenance(result.recipe);
        std::cout << "recipe (minus provenance): " << (recipeSame ? "identical" : "DIFF") << "\n";
        failures += !recipeSame;
        return failures ? 1 : 0;
    } catch (const std::exception& e) {
        std::cout << "FAILED " << e.what() << "\n";
        return 1;
    }
}
