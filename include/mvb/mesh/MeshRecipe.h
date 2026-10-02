// MESH RECIPE (2026-09-13, Alf: "I want all these knobs later accessible ... a software where these
// could be controlled by a human").
//
// One JSON document holds every meshing knob of omfem_mesh3d and omfem_skinlayer, with units and
// documentation, so a person (or a GUI) edits a document instead of environment variables, and
// every mesh records how it was made: the tools write the EFFECTIVE recipe next to their output
// (<out.msh>.recipe.json). Internally the knobs still travel as the OMFEM_* environment variables
// the meshing code reads; this table is the single place that names them.
//
// Rules: unknown keys throw (no silent typos); a null value means "the tool's default"; a value in
// the recipe overrides an inherited environment variable (the document is the source of truth, the
// override is reported); lengths are in metres.
#pragma once
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace mvb::mesh {

struct RecipeKnob {
    const char* section;   // JSON object
    const char* key;       // JSON key inside the section
    const char* env;       // environment variable the meshing code reads
    const char* kind;      // "number" | "string" | "flag" (present/absent) 
    const char* unit;      // "m", "-", "cells", ...
    const char* def;       // the tool's default, as text (documentation, also written when unset)
    const char* doc;       // one line
};

// The complete knob table (mesh3d sections + the "layers" section of omfem_skinlayer).
const std::vector<RecipeKnob>& mesh_recipe_knobs();

// Read a recipe file and export its values into the environment (setenv, overwriting). Throws on
// unknown sections/keys or wrong value kinds. Returns the names of environment variables that were
// already set and got overridden, for the caller to report.
std::vector<std::string> apply_mesh_recipe(const std::string& path);

// The effective recipe: every knob with the value currently in the environment, or its documented
// default when unset ("default": true marks those), plus a provenance block. `sections` limits the
// output to some sections (empty = all).
nlohmann::json effective_mesh_recipe(const std::vector<std::string>& sections,
                                     const nlohmann::json& provenance);

// A template with every knob at null and a "_doc" line per knob (for `--recipe-template`).
nlohmann::json mesh_recipe_template();

void write_mesh_recipe(const std::string& path, const nlohmann::json& recipe);

} // namespace mvb::mesh
