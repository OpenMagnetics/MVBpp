#include "mvb/mesh/Mesher.h"
#include "mvb/mesh/BuildRev.h"
#include "mvb/mesh/Mesh3d.h"
#include "mvb/mesh/MeshOptions.h"
#include "mvb/mesh/MeshRecipe.h"
#include "mvb/Utils.h"

#include "constructive_models/Core.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"

#include <gmsh.h>
#include <mmg/common/mmgversion.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

extern char** environ;

namespace mvb::mesh {

namespace {

// The sections omfem_mesh3d records in its effective recipe (the skin layer's "layers" is not
// part of a mesh3d run; "geometry" holds the MVB++ drawing knobs a corpus plan sets).
const std::vector<std::string> kMesh3dSections{"wire", "conductor", "core", "air", "gap", "skin_mapped", "mesher", "geometry"};

// While the move is under way the knobs travel as OMFEM_* environment variables, read inside the
// moved code, and MVB++'s geometry reads MVB_* ones. A library call must not take hidden inputs
// from its caller's environment, nor leave its own behind: any OMFEM_* or MVB_* already set is
// refused (most of them are not recipe knobs, so a mesh made with them could not say so), and
// every knob the recipe set is unset on exit, throw included.
struct RecipeEnvironment {
    RecipeEnvironment() {
        std::string set;
        for (char** e = environ; e && *e; ++e)
            if (std::string(*e).rfind("OMFEM_", 0) == 0 || std::string(*e).rfind("MVB_", 0) == 0)
                set += std::string(set.empty() ? "" : ", ") + std::string(*e).substr(0, std::string(*e).find('='));
        if (!set.empty())
            throw std::runtime_error("meshMagnetic: takes every meshing knob from the recipe, but the process "
                                     "environment already sets " + set + "; unset them (a knob the recipe lacks has to join the recipe table first)");
    }
    ~RecipeEnvironment() {
        for (const auto& k : mesh_recipe_knobs()) unsetenv(k.env);
        unsetenv("OMFEM_KEEP_STEP");
    }
    RecipeEnvironment(const RecipeEnvironment&) = delete;
    RecipeEnvironment& operator=(const RecipeEnvironment&) = delete;
};

// A private directory for the moved code's files, removed on exit.
struct ScratchDir {
    std::filesystem::path path;
    ScratchDir() {
        std::string tmpl = (std::filesystem::temp_directory_path() / "mvbpp_mesh_XXXXXX").string();
        if (!mkdtemp(tmpl.data()))
            throw std::runtime_error("meshMagnetic: cannot create a scratch directory under " +
                                     std::filesystem::temp_directory_path().string());
        path = tmpl;
    }
    ~ScratchDir() {
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
    }
    ScratchDir(const ScratchDir&) = delete;
    ScratchDir& operator=(const ScratchDir&) = delete;
};

std::string slurp(const std::filesystem::path& p, const char* what) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error(std::string("meshMagnetic: the mesher wrote no ") + what + " (" + p.string() + ")");
    std::stringstream s;
    s << f.rdbuf();
    return s.str();
}

double envNumber(const char* name) {
    const char* v = std::getenv(name);
    return v ? std::atof(v) : 0.0;
}

// Each volume region's MAS material, from the enriched magnetic. "core" and "winding_<name>" are
// the regions mesh3d_from_mas names; air and every face set have none. A winding drawn at its
// outer (insulated) diameter is copper and enamel in one region, so it gets no material: an
// exporter that needs one then refuses it.
void assignMaterials(Mesh& mesh, const nlohmann::json& magnetic, bool outerWireDiameter) {
    auto enriched = mvb::magnetic_autocomplete_safe(magnetic);
    for (auto& r : mesh.regions) {
        if (r.dimension != 3 || r.name == "air") {
            r.material = {MaterialKind::None, {}};
            continue;
        }
        if (r.name == "core") {
            auto materials = enriched.get_mutable_core().resolve_materials();
            if (materials.size() != 1)
                throw std::runtime_error("meshMagnetic: the core has " + std::to_string(materials.size()) +
                                         " materials but the mesh one 'core' region");
            nlohmann::json j;
            MAS::to_json(j, materials.front());
            r.material = {MaterialKind::Core, j};
            continue;
        }
        if (r.name.rfind("winding_", 0) == 0) {
            if (outerWireDiameter) {
                r.material = {};
                continue;
            }
            const std::string winding = r.name.substr(8);
            const OpenMagnetics::Winding* fd = nullptr;
            for (const auto& w : enriched.get_coil().get_functional_description())
                if (w.get_name() == winding) fd = &w;
            if (!fd)
                throw std::runtime_error("meshMagnetic: region '" + r.name + "' names a winding the coil does not have");
            const auto* wire = std::get_if<OpenMagnetics::Wire>(&fd->get_wire());
            if (!wire)
                throw std::runtime_error("meshMagnetic: winding '" + winding + "' has an unresolved wire after enrichment");
            nlohmann::json j;
            MAS::to_json(j, OpenMagnetics::Wire::resolve_material(*wire));
            r.material = {MaterialKind::Wire, j};
            continue;
        }
        throw std::runtime_error("meshMagnetic: volume region '" + r.name + "' has no known MAS material role");
    }
}

}  // namespace

MeshResult meshMagnetic(const nlohmann::json& magnetic, const nlohmann::json& recipe, bool withStep) {
    if (!magnetic.is_object() || !magnetic.contains("core") || !magnetic.contains("coil"))
        throw std::runtime_error("meshMagnetic: expects a MAS magnetic (an object with core and coil)");
    if (!recipe.is_object()) throw std::runtime_error("meshMagnetic: the recipe must be a JSON object");
    RecipeEnvironment env;
    ScratchDir dir;
    {
        std::ofstream r(dir.path / "recipe.json");
        r << recipe.dump();
    }
    apply_mesh_recipe((dir.path / "recipe.json").string());
    // The moved mesher deletes its STEP unless told to keep it (unset again on exit, above).
    if (withStep) setenv("OMFEM_KEEP_STEP", "1", 1);

    // The options exactly as omfem_mesh3d builds them (its positional defaults, then the
    // recipe's environment), in the same order, so one recipe makes one mesh.
    MeshOptions opt;
    opt.out_msh = (dir.path / "real.msh").string();
    opt.polygon_segments = 12;
    if (const char* sg = std::getenv("OMFEM_MESH_SEGMENTS")) opt.polygon_segments = std::atoi(sg);
    opt.paint_coating = false;
    opt.symmetry = "none";
    if (const char* wd = std::getenv("OMFEM_WIRE_DIM")) {
        const std::string v(wd);
        if (v == "conducting") opt.paint_coating = false;
        else if (v == "outer") opt.paint_coating = true;
        else throw std::runtime_error("meshMagnetic: wire.dimension must be 'conducting' or 'outer', got '" + v + "'");
    }
    if (std::getenv("OMFEM_AIR_TARGET")) opt.air_target = envNumber("OMFEM_AIR_TARGET");
    if (std::getenv("OMFEM_CORE_TARGET")) opt.core_target = envNumber("OMFEM_CORE_TARGET");
    if (std::getenv("OMFEM_COND_TARGET")) opt.conductor_target = envNumber("OMFEM_COND_TARGET");
    if (std::getenv("OMFEM_AIR_MARGIN")) {
        opt.air_margin = envNumber("OMFEM_AIR_MARGIN");
        if (!(opt.air_margin > 1.0)) throw std::runtime_error("meshMagnetic: air.margin must be > 1 (box = bbox * margin)");
    }

    const std::string out = mesh3d_from_mas(magnetic, opt);

    MeshResult result;
    result.mesh = importMesh(slurp(out, "mesh"), "msh2");
    for (const char* suffix : {".cadedges", ".path", ".leads.json"})
        result.sidecars[suffix] = slurp(out + suffix, suffix);
    if (withStep)
        result.sidecars[".step"] = slurp(out + ".step", "STEP");
    else if (std::filesystem::exists(out + ".step"))
        throw std::runtime_error("meshMagnetic: the mesher left its STEP behind (" + out + ".step)");
    nlohmann::json provenance = {{"tool", "mvb::mesh::meshMagnetic"},
                                 {"segments", opt.polygon_segments},
                                 {"wire_dimension", opt.paint_coating ? "outer" : "conducting"},
                                 {"symmetry", opt.symmetry}};
    result.recipe = effective_mesh_recipe(kMesh3dSections, provenance);
    result.recipe["provenance"]["build"] = buildRevision();
    result.sidecars[".recipe.json"] = result.recipe.dump(2) + "\n";
    assignMaterials(result.mesh, magnetic, opt.paint_coating);
    return result;
}

nlohmann::json recipeTemplate() {
    return mesh_recipe_template();
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
