#include "mvb/mesh/MeshRecipe.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <ctime>

namespace mvb::mesh {

const std::vector<RecipeKnob>& mesh_recipe_knobs() {
    static const std::vector<RecipeKnob> K = {
        // ---- what is drawn -------------------------------------------------------------------
        {"wire", "dimension", "OMFEM_WIRE_DIM", "string", "-", "conducting",
         "conducting = bare copper cross-section (winding-loss FEM); outer = insulated envelope"},
        {"wire", "segments", "OMFEM_MESH_SEGMENTS", "number", "-", "12",
         "MVB++ curve faceting of wires and core (0 = exact round surfaces)"},
        // ---- region sizes --------------------------------------------------------------------
        {"conductor", "target_m", "OMFEM_COND_TARGET", "number", "m", "auto",
         "element size on and inside the conductors; auto = chord bound from the real turn clearances; corpus rule max(delta/2, dmin/4)"},
        {"conductor", "halo", "OMFEM_COND_HALO", "number", "x target", "8",
         "reach of the conductor-size refinement into the air, in multiples of conductor.target_m; grows to the air ceiling by 2x that radius"},
        {"core", "target_m", "OMFEM_CORE_TARGET", "number", "m", "3e-4",
         "element size in the core (the core is a few % of the tets; it only gets fine where it touches the halo)"},
        {"air", "target_m", "OMFEM_AIR_TARGET", "number", "m", "1e-3", "element size ceiling in the far air"},
        {"air", "box_margin", "OMFEM_AIR_MARGIN", "number", "x bbox", "1.5",
         "air box = assembly bounding box scaled about its centre (E4: 2.0 moved P_ac by -0.85 %)"},
        {"gap", "divisions", "OMFEM_GAP_DIV_3D", "number", "-", "off",
         "opt-in isotropic refinement across each core gap: elements per gap length"},
        {"gap", "radius_m", "OMFEM_GAP_R_3D", "number", "m", "column radius", "lateral reach of the gap refinement"},
        // ---- rectangular wire: mapped hex skin blocks (built inside mesh3d) -------------------
        {"skin_mapped", "enabled", "OMFEM_SKIN_MAPPED", "flag", "-", "off",
         "rectangular/foil wire: graded hex blocks toward all four faces instead of tets"},
        {"skin_mapped", "delta_m", "OMFEM_SKIN_DELTA", "number", "m", "required when enabled", "skin depth the blocks are graded to"},
        {"skin_mapped", "cells", "OMFEM_SKIN_CELLS", "number", "cells/delta", "2", "cells per skin depth at the surface (E5: 1 = +0.9 % P_ac, -16 % time)"},
        {"skin_mapped", "growth", "OMFEM_SKIN_GROWTH", "number", "-", "1.3", "geometric growth of the block layers toward the wire centre"},
        {"skin_mapped", "h_along_m", "OMFEM_SKIN_HLONG", "number", "m", "1e-3", "block length along the wire"},
        // ---- mesher --------------------------------------------------------------------------
        {"mesher", "algo3d", "OMFEM_ALGO3D", "number", "gmsh id", "auto (Delaunay 1 -> HXT 10 chain)",
         "gmsh Mesh.Algorithm3D override: 1 Delaunay, 10 HXT (parallel; heap crashes on some surfaces), 4 Frontal"},
        {"mesher", "threads", "OMFEM_GMSH_THREADS", "number", "-", "1 (gmsh default)",
         "gmsh threads for surface + HXT volume meshing (ETD34 base: 12 threads 175 s vs 337 s; mesh differs, re-validate before default)"},
        {"mesher", "retries", "OMFEM_MESH_RETRIES", "number", "-", "3", "extra attempts of the algorithm/size ladder after a failed 3D pass"},
        {"mesher", "quality_floor", "OMFEM_MESH_QMIN", "number", "-", "1e-6", "minimum accepted tet quality (gamma) before the mesh is refused"},
        {"mesher", "fragment_tolerance_m", "OMFEM_FRAG_TOL", "number", "m", "0", "OCC fuzzy value for the solid fragment (0 = exact)"},
        {"mesher", "port_inset", "OMFEM_PORT_INSET", "number", "x cap radius", "0.5 axis-aligned caps, 1.5 oblique",
         "how far a lead-end port face is pulled inside the air box, in cap radii"},
        // ---- omfem_skinlayer: MMG anisotropic layers on the conformal mesh (round wire) --------
        {"layers", "enabled", "OMFEM_LAYERS_ENABLED", "flag", "-", "off", "run omfem_skinlayer after the conformal mesh (round/litz wire)"},
        {"layers", "delta_m", "OMFEM_LAYERS_DELTA", "number", "m", "required when enabled",
         "layer depth = skin depth of the highest harmonic MKF keeps at a 10 % threshold (corpus rule)"},
        // Inputs have no default (the tool refuses a missing one; the design plan states them and why).
        // Everything else the tool derives and records under provenance.derived of its recipe.
        {"layers", "cells", "OMFEM_LAYERS_CELLS", "number", "cells/delta", "required", "cells across one skin depth at the copper surface"},
        {"layers", "h_tan_m", "OMFEM_LAYERS_HTAN", "number", "m", "required", "element size along the copper surface / wire"},
        {"layers", "deep_m", "OMFEM_SKIN_DEEP", "number", "m", "required", "layer size cap through the conductor thickness (h_deep)"},
        {"layers", "mmg_hgrad", "OMFEM_MMG_HGRAD", "number", "-", "required", "MMG gradation (max size ratio between neighbours); the layer count derives from it"},
        {"layers", "n_layers", "OMFEM_SKIN_LAYERS", "number", "-", "derived: 1 + ceil(ln(h_deep/h_surface)/ln(hgrad))", "layers between the wall and deep_m; the growth derives from it"},
        {"layers", "mmg_memory_mb", "OMFEM_MMG_MEM", "number", "MB", "derived: memory available to the process", "MMG memory cap"},
    };
    return K;
}

static std::string num_to_env(const nlohmann::json& v) {
    if (v.is_number_integer()) return std::to_string(v.get<long long>());
    char b[64]; std::snprintf(b, sizeof b, "%.12g", v.get<double>()); return b;
}

std::vector<std::string> apply_mesh_recipe(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("mesh recipe: cannot open " + path);
    nlohmann::json r; in >> r;
    if (!r.is_object()) throw std::runtime_error("mesh recipe: top level must be an object");
    const auto& K = mesh_recipe_knobs();
    std::vector<std::string> overridden, unknown;
    for (auto it = r.begin(); it != r.end(); ++it) {
        const std::string sec = it.key();
        if (sec == "recipe_version" || sec == "provenance" || sec == "_doc" || sec == "name" || sec == "notes") continue;
        if (!it.value().is_object()) { unknown.push_back(sec); continue; }
        for (auto jt = it.value().begin(); jt != it.value().end(); ++jt) {
            const std::string key = jt.key();
            if (key == "_doc" || key == "default") continue;
            const RecipeKnob* k = nullptr;
            for (const auto& kk : K) if (sec == kk.section && key == kk.key) { k = &kk; break; }
            if (!k) { unknown.push_back(sec + "." + key); continue; }
            const auto& v = jt.value();
            const bool had = std::getenv(k->env) != nullptr;
            if (v.is_null()) { if (had) { unsetenv(k->env); overridden.push_back(std::string(k->env) + " (unset by null)"); } continue; }
            std::string sv;
            if (std::string(k->kind) == "flag") {
                if (!v.is_boolean()) throw std::runtime_error("mesh recipe: " + sec + "." + key + " must be true/false");
                if (!v.get<bool>()) { if (had) { unsetenv(k->env); overridden.push_back(k->env); } continue; }
                sv = "1";
            } else if (std::string(k->kind) == "number") {
                if (!v.is_number()) throw std::runtime_error("mesh recipe: " + sec + "." + key + " must be a number (" + k->unit + ")");
                sv = num_to_env(v);
            } else {
                if (!v.is_string()) throw std::runtime_error("mesh recipe: " + sec + "." + key + " must be a string");
                sv = v.get<std::string>();
            }
            if (had && sv != std::getenv(k->env)) overridden.push_back(std::string(k->env) + "=" + std::getenv(k->env) + " -> " + sv);
            setenv(k->env, sv.c_str(), 1);
        }
    }
    if (!unknown.empty()) {
        std::string m = "mesh recipe " + path + ": unknown key(s):";
        for (const auto& u : unknown) m += " " + u;
        throw std::runtime_error(m + " (see omfem_mesh3d --recipe-template for the accepted set)");
    }
    return overridden;
}

nlohmann::json effective_mesh_recipe(const std::vector<std::string>& sections, const nlohmann::json& provenance) {
    nlohmann::json r; r["recipe_version"] = 1;
    for (const auto& k : mesh_recipe_knobs()) {
        if (!sections.empty()) { bool keep = false; for (const auto& s : sections) if (s == k.section) keep = true; if (!keep) continue; }
        nlohmann::json e;
        const char* v = std::getenv(k.env);
        if (std::string(k.kind) == "flag") e = (v != nullptr);
        else if (!v) e = nullptr;
        else if (std::string(k.kind) == "number") { char* end = nullptr; double d = std::strtod(v, &end); e = (end && *end == 0) ? nlohmann::json(d) : nlohmann::json(std::string(v)); }
        else e = std::string(v);
        r[k.section][k.key] = e;
        if (!v) r[k.section][std::string(k.key) + "_default"] = k.def;
        r[k.section]["_unit_" + std::string(k.key)] = k.unit;
    }
    r["provenance"] = provenance;
    std::time_t t = std::time(nullptr); char buf[32]; std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%S", std::localtime(&t));
    r["provenance"]["written"] = buf;
    return r;
}

nlohmann::json mesh_recipe_template() {
    nlohmann::json r; r["recipe_version"] = 1;
    r["_doc"] = "OMFEM mesh recipe: null = tool default; lengths in metres unless the unit says otherwise; unknown keys are refused";
    for (const auto& k : mesh_recipe_knobs()) {
        r[k.section][k.key] = nullptr;
        r[k.section]["_doc_" + std::string(k.key)] = std::string(k.doc) + " [" + k.unit + "; default " + k.def + "; env " + k.env + "]";
    }
    return r;
}

void write_mesh_recipe(const std::string& path, const nlohmann::json& recipe) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("mesh recipe: cannot write " + path);
    out << recipe.dump(2) << "\n";
}

} // namespace mvb::mesh
