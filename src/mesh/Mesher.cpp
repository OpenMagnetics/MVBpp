#include "mvb/mesh/Mesher.h"
#include "mvb/mesh/BuildRev.h"

#include <gmsh.h>
#include <mmg/common/mmgversion.h>

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

nlohmann::json buildRevision() {
    // gmsh::option needs an initialised gmsh; initialise only if nobody has (a caller in the
    // middle of a model must not have it torn down), and leave it as found.
    const bool mine = !gmsh::isInitialized();
    if (mine) gmsh::initialize();
    std::string gmshVersion, gmshBuildOptions;
    gmsh::option::getString("General.Version", gmshVersion);
    // The modules compiled in (Metis or not, OCC, MMG...): two gmsh builds of one version can
    // mesh differently, so the version alone does not say which gmsh made a mesh.
    gmsh::option::getString("General.BuildOptions", gmshBuildOptions);
    if (mine) gmsh::finalize();
    return {{"mvbpp", MVBPP_BUILD_REV_MVBPP},
            {"mkf", MVBPP_BUILD_REV_MKF},
            {"mas", MVBPP_BUILD_REV_MAS},
            {"gmsh", gmshVersion},
            {"gmshBuildOptions", gmshBuildOptions},
            {"mmg", MMG_VERSION_RELEASE}};
}


}  // namespace mvb::mesh
