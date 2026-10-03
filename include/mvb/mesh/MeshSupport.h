#pragma once
// Moved verbatim from OMFEM src/meshing/MasMesher.cpp (ABT #1588, step 2): namespace omfem -> mvb::mesh,
// MasMeshOptions -> MeshOptions, file-local helpers shared with OMFEM's 2D mesher made external. No logic changed.

#include <optional>
#include <regex>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "constructive_models/Magnetic.h"

namespace mvb::mesh {

// ---- ABT #1169 (WP0): strict solid-name classification -----------------------------------
//
// MVB++ names the geometry; OMFEM decides which FEM region each name becomes. Exposed so the
// mapping can be unit-tested directly — including the fact that an unrecognised name THROWS
// instead of being simulated as ferrite.
struct SolidNameContext {
    std::string core_name;     // MVB++ names core pieces "<core_name>" / "<core_name>_<i>"
    std::string bobbin_name;   // MVB++ names the former after the MAS bobbin, else "Bobbin"
};
using NameContext = SolidNameContext;

// gmsh threads from OMFEM_GMSH_THREADS (unset = gmsh default, 1).
void apply_gmsh_threads();
// gmsh::write to "<out>.part.msh", then rename onto <out>.
void write_mesh_atomically(const std::string& out);

// Turn label "<winding> parallel <par> turn <idx>".
extern const std::regex kTurnRe;

NameContext name_context(const OpenMagnetics::Magnetic& enriched);
bool ends_with(const std::string& s, const std::string& suffix);
std::optional<std::string> classify_base(const std::string& name, const NameContext& ctx,
                                         std::string& winding, int& index);
std::string classify(const std::string& raw, std::string& winding, int& index,
                     const NameContext& ctx);
std::string region_for_role(const std::string& role, const std::string& winding, int index,
                            const std::string& turnRegion);

// Role of a named solid: "core", "turn", "winding", "bobbin", "terminal", "sleeve", "solder",
// "insulation", "spacer", "shunt", "divider", "pin"; "" for a shape to ignore entirely.
// `winding` and `index` are filled where the name carries them.
// Throws std::runtime_error on any name that matches no family.
std::string classify_solid_name(const std::string& name, const SolidNameContext& ctx,
                                std::string& winding, int& index);

// The FEM region string a classified solid becomes. `turnRegion` is used for role "turn" only.
std::string fem_region_for_role(const std::string& role, const std::string& winding, int index,
                                const std::string& turnRegion);

struct CoreGap { double yc; double len; double xc; };
std::vector<CoreGap> extract_core_gaps(const nlohmann::json& enr);
double central_column_half_width(const nlohmann::json& enr);

}  // namespace mvb::mesh
