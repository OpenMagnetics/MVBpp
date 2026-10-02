#include "mvb/mesh/Mesher.h"

#include <stdexcept>

// ABT #1588 skeleton: the target links gmsh + MMG; the meshing code moves in from OMFEM in
// gated steps. Until then every entry point refuses loudly -- no caller may get an empty mesh.

namespace mvb::mesh {

MeshResult meshMagnetic(const nlohmann::json&, const nlohmann::json&) {
    throw std::runtime_error("mvb::mesh::meshMagnetic: not implemented yet (ABT #1588 move in progress)");
}

nlohmann::json recipeTemplate() {
    throw std::runtime_error("mvb::mesh::recipeTemplate: not implemented yet (ABT #1588 move in progress)");
}

std::string exportMesh(const Mesh&, const std::string&, const std::string&) {
    throw std::runtime_error("mvb::mesh::exportMesh: not implemented yet (ABT #1588 move in progress)");
}

}  // namespace mvb::mesh
