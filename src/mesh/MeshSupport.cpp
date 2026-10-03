#include "mvb/mesh/MeshSupport.h"
// Moved verbatim from OMFEM src/meshing/MasMesher.cpp (ABT #1588, step 2): namespace omfem -> mvb::mesh,
// MasMeshOptions -> MeshOptions, file-local helpers shared with OMFEM's 2D mesher made external. No logic changed.

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <variant>

#include <gmsh.h>

namespace mvb::mesh {

using nlohmann::json;

// GMSH THREADS (2026-09-13): the mesher never set General.NumThreads, so every volume mesh of the corpus ran on ONE
// thread (gmsh's default) -- the "12" the corpus runner passes is the MVB++ faceting, not a thread count. Measured on
// the ETD34 base mesh: 12 threads = mesh generation 175 s vs 337 s (HXT Delaunay parallel; surface meshing, boundary
// recovery and optimisation stay largely serial), tet count 709k vs 635k (HXT's parallel insertion order gives a
// different, equivalent mesh -- re-validate P_total on ETD34/03 before making it the default). Knob: OMFEM_GMSH_THREADS
// (recipe mesher.threads); unset = gmsh default (1) so running corpora are reproducible.
void apply_gmsh_threads() {
    if (const char* t = std::getenv("OMFEM_GMSH_THREADS")) {
        const int n = std::atoi(t);
        if (n < 1) throw std::runtime_error("OMFEM_GMSH_THREADS must be a positive integer");
        gmsh::option::setNumber("General.NumThreads", n);
        gmsh::option::setNumber("Mesh.MaxNumThreads2D", n);
        gmsh::option::setNumber("Mesh.MaxNumThreads3D", n);
        std::fprintf(stderr, "[mesh3d] gmsh threads: %d (OMFEM_GMSH_THREADS)\n", n);
    }
}

// Write the mesh to "<out>.part" and RENAME it into place, so the final path appears only when
// the file is complete. gmsh::write goes straight to the destination, so the mesh existed from its
// first byte and grew for minutes (41 MB here, 182 MB for a skin mesh) -- a second runner testing
// existence could pick up a half-written file, and a mesher that died mid-write left a truncated
// mesh on the final path that would be "reused" forever. Same directory, so the rename is atomic.
void write_mesh_atomically(const std::string& out) {
    // The temp name must keep the FINAL extension: gmsh picks its output format from the
    // extension, so writing to "<out>.part" made it refuse with "Unknown output file format"
    // AFTER meshing -- 2.1M tets and 735 s thrown away per attempt, on every design that meshed.
    // "<out>.part.msh" keeps the .msh that gmsh needs and still renames onto <out>.
    const std::string part = out + ".part.msh";
    gmsh::write(part);
    if (std::rename(part.c_str(), out.c_str()) != 0) {
        std::remove(part.c_str());
        throw std::runtime_error("mesh3d: could not rename " + part + " onto " + out + ": "
                                 + std::strerror(errno));
    }
}

// Turn label "<winding> parallel <par> turn <idx>"; core "<base>[_<i>]"; Bobbin/FR4Board.
const std::regex kTurnRe(R"(^(.+) parallel (\d+) turn (\d+)$)");
const std::regex kSuffixRe(R"(^(.+)_(\d+)$)");
// Real-winding continuous conductor "<winding> parallel <par>" (ONE solid per winding-parallel,
// no per-turn "turn <idx>" suffix) -- the useRealWindingGeometry/femReady geometry.
const std::regex kParallelRe(R"(^(.+) parallel (\d+)$)");
// Port cap face at a free end of a continuous conductor: "<winding> parallel <p> terminal <k>".
const std::regex kTerminalRe(R"(^(.+) parallel (\d+) terminal (\d+)$)");
// Lead sleeve (proposal D2 / WP7): "<winding> parallel <p> <entrance|exit> sleeve".
const std::regex kSleeveRe(R"(^(.+) parallel (\d+) (?:entrance|exit) sleeve$)");
// Foil lead wire: "<winding> parallel <p> <entrance|exit> lead" -- the round copper lead MVB++ solders onto a foil
// sheet (the sheet itself is "<winding> parallel <p>"). Copper of that winding, so the same "winding" role and region
// as the sheet, joined to it through the solder body below (two_switch_forward_transformer_complete, 2026-09-23).
const std::regex kFoilLeadRe(R"(^(.+) parallel (\d+) (?:entrance|exit) lead$)");
// Foil terminal solder body: "<winding> parallel <p> <entrance|exit> lead solder joint".
const std::regex kSolderRe(R"(^(.+) parallel (\d+) .* solder joint$)");
// Bobbin chamber wall (WP2): "<bobbin> divider <i>".  Bobbin pin (WP4): "<bobbin> pin <name>".
const std::regex kDividerRe(R"(^.+ divider (\d+)$)");
const std::regex kPinRe(R"(^(.+) pin (.+)$)");
// Core spacer (WP1) and magnetic shunt (WP6): "Spacer_<i>", "Shunt_<i>".
const std::regex kSpacerRe(R"(^Spacer_(\d+)$)");
const std::regex kShuntRe(R"(^Shunt_(\d+)$)");

// ABT #1169 (WP0). What classify() needs to know about THIS magnetic before it can be strict.
//
// MVB++ names a core piece after the MAS core ("<core name>", or "<core name>_<i>" when the
// core has several pieces) and the former after the MAS bobbin. Those names are arbitrary
// catalogue strings, so "anything I do not recognise must be a core piece" used to be the only
// way to spot them -- and that is exactly the fallback that would silently simulate a spacer,
// a pin or a shunt as ferrite. Handing classify() the two names the magnetic actually carries
// turns core/bobbin recognition POSITIVE, which is what lets everything else throw.
// The public SolidNameContext (declared in MasMesher.hpp) IS this context; alias it so the
// helpers below read naturally and the exported entry point needs no conversion.
using NameContext = SolidNameContext;

NameContext name_context(const OpenMagnetics::Magnetic& enriched) {
    NameContext ctx;
    ctx.core_name = enriched.get_core().get_name().value_or("Core");
    ctx.bobbin_name = "Bobbin";
    const auto& bobVar = enriched.get_coil().get_bobbin();
    if (const auto* s = std::get_if<std::string>(&bobVar)) { if (!s->empty()) ctx.bobbin_name = *s; }
    else if (const auto* b = std::get_if<OpenMagnetics::Bobbin>(&bobVar)) {
        const auto nm = b->get_name(); if (nm && !nm->empty()) ctx.bobbin_name = *nm;
    }
    return ctx;
}

bool ends_with(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// One classification attempt on an already-normalised name. Returns the role, "" for a shape
// that must be IGNORED entirely, or nullopt when the name is not recognised at all (the caller
// then peels one "_<i>" child suffix and tries again, and finally throws).
std::optional<std::string> classify_base(const std::string& name, const NameContext& ctx,
                                         std::string& winding, int& index) {
    std::smatch m;
    if (std::regex_match(name, m, kTurnRe)) {
        winding = m[1].str(); index = std::stoi(m[3].str()); return std::string("turn");
    }
    if (std::regex_match(name, m, kTerminalRe)) {
        winding = m[1].str(); index = std::stoi(m[3].str()); return std::string("terminal");
    }
    if (std::regex_match(name, m, kSleeveRe)) {
        winding = m[1].str(); index = std::stoi(m[2].str()); return std::string("sleeve");
    }
    if (std::regex_match(name, m, kSolderRe)) {
        winding = m[1].str(); index = std::stoi(m[2].str()); return std::string("solder");
    }
    if (std::regex_match(name, m, kFoilLeadRe)) {   // after the solder body, whose name extends this one
        winding = m[1].str(); index = std::stoi(m[2].str()); return std::string("winding");
    }
    // Real-winding continuous conductor "<winding> parallel <par>": ONE solid per winding-parallel.
    // Roled "winding" (distinct from per-turn "turn") so the 3D tagging lumps it into a coil region.
    if (std::regex_match(name, m, kParallelRe)) {
        winding = m[1].str(); index = std::stoi(m[2].str()); return std::string("winding");
    }
    if (name.rfind("insulation_layer", 0) == 0) return std::string("insulation");
    if (std::regex_match(name, m, kSpacerRe))  { index = std::stoi(m[1].str()); return std::string("spacer"); }
    if (std::regex_match(name, m, kShuntRe))   { index = std::stoi(m[1].str()); return std::string("shunt"); }
    if (std::regex_match(name, m, kDividerRe)) { index = std::stoi(m[1].str()); return std::string("divider"); }
    // A toroid base plays the former's role: standoff + pin carrier. Plastic, like the bobbin.
    if (!ctx.bobbin_name.empty() && name == ctx.bobbin_name + " base") return std::string("bobbin");
    if (ends_with(name, " base") && name.rfind("Bobbin", 0) == 0)      return std::string("bobbin");
    if (std::regex_match(name, m, kPinRe))     { winding = m[2].str(); return std::string("pin"); }
    if (name == "Bobbin" || name == "FR4Board") return std::string("");   // ignore
    // Positive core / bobbin recognition from the magnetic's own names (see NameContext).
    auto matches_base = [&](const std::string& base) {
        if (base.empty()) return false;
        if (name == base) return true;
        if (name.size() > base.size() + 1 && name.compare(0, base.size(), base) == 0 &&
            name[base.size()] == '_') {
            const std::string tail = name.substr(base.size() + 1);
            if (tail.find_first_not_of("0123456789") == std::string::npos) {
                index = std::stoi(tail);
                return true;
            }
        }
        return false;
    };
    if (matches_base(ctx.core_name))   return std::string("core");
    if (matches_base(ctx.bobbin_name)) return std::string("bobbin");
    return std::nullopt;
}

// Classify an MVB++ shape NAME into an OMFEM FEM REGION role. Regions (the FEM
// subdomain tagging) are OMFEM's responsibility — MVB++ only names the geometry, OMFEM
// interprets those names into simulation regions. Handles STEP-label prefixes/suffixes
// and "_<int>" children from section/symmetry cuts. winding/index set for turns.
//
// Roles: "core", "turn", "winding", "bobbin", "terminal", "sleeve", "solder", "insulation",
// "spacer", "shunt", "divider", "pin"; "" for a shape to ignore.
//
// ABT #1169: there is NO fallback. A name this function does not recognise THROWS. It used to
// return "core" for anything unknown, which meant the first spacer, pin, divider, shunt or
// sleeve MVB++ ever drew would have been meshed and simulated as ferrite — silently, with a
// plausible-looking answer. Adding a solid to MVB++ now forces a branch here.
std::string classify(const std::string& raw, std::string& winding, int& index,
                     const NameContext& ctx) {
    static const std::string kCoating = " coating";
    // One candidate spelling: peel "_<i>" child suffixes (from section/symmetry cuts) until the
    // name is recognised. Returns the role, or nullopt when no peeling of THIS spelling matches.
    auto attempt = [&](std::string name) -> std::optional<std::string> {
        // "<x> coating" is the outer insulated shell of <x>; it classifies as whatever <x> is
        // (the 3D tagger strips the suffix itself before calling, and tags the shell separately).
        if (ends_with(name, kCoating)) name = name.substr(0, name.size() - kCoating.size());
        index = 0;
        int childIndex = -1;
        for (;;) {
            if (auto role = classify_base(name, ctx, winding, index)) {
                // A core piece cut into children by the section/symmetry keeps the CHILD index
                // as its region index, so siblings stay distinct regions (previous behaviour).
                if (*role == "core" && childIndex >= 0) index = childIndex;
                return role;
            }
            std::smatch ms;
            if (!std::regex_match(name, ms, kSuffixRe)) return std::nullopt;
            if (childIndex < 0) { try { childIndex = std::stoi(ms[2].str()); } catch (...) { return std::nullopt; } }
            name = ms[1].str();
        }
    };

    std::string name = raw;
    if (name.rfind("Shapes/", 0) == 0) name = name.substr(7);
    if (auto role = attempt(name)) return *role;
    // STEP label paths ("Shapes/<product>/<child>") are only split on '/' as a LAST resort:
    // core shape names contain slashes themselves ("ETD 44/22/15"), and splitting first turned
    // every such core into the unrecognised fragment "ETD 44".
    if (auto slash = name.find('/'); slash != std::string::npos)
        if (auto role = attempt(name.substr(0, slash))) return *role;
    throw std::runtime_error(
        "mesh: unclassified solid name '" + raw + "' (core '" + ctx.core_name + "', bobbin '" +
        ctx.bobbin_name + "'). Every solid MVB++ draws must map to an OMFEM region — add a "
        "branch to classify() rather than letting a new part be simulated as core (ABT #1169).");
}

// ABT #1169: the FEM region string for a classified solid. One place, so the 2D section tagger
// and the 3D volume tagger cannot drift into naming the same physics differently. `turnRegion`
// is what the caller wants a per-turn conductor called (the two paths differ: the 2D tagger
// numbers turns, the 3D tagger may lump them), and is used for role "turn" only.
std::string region_for_role(const std::string& role, const std::string& winding, int index,
                            const std::string& turnRegion) {
    if (role == "core")       return "core_" + std::to_string(index);
    if (role == "spacer")     return "spacer_" + std::to_string(index);
    if (role == "shunt")      return "shunt_" + std::to_string(index);
    if (role == "divider")    return "divider_" + std::to_string(index);
    if (role == "insulation") return "insulation_" + std::to_string(index);
    if (role == "solder")     return "solder_" + winding + "_" + std::to_string(index);
    if (role == "sleeve")     return "sleeve_" + winding + "_" + std::to_string(index);
    if (role == "terminal")   return "terminal_" + winding + "_" + std::to_string(index);
    if (role == "winding")    return "winding_" + winding;
    if (role == "turn")       return turnRegion;
    // A pin is copper and a base is plastic, but neither carries a driven current and neither is
    // magnetic: both join the bobbin's thermal bucket until they are given their own conductivity
    // (proposal WP0). They are named apart so that day is a rename, not a re-classification.
    if (role == "pin")        return "bobbin_pin_" + winding;
    if (role == "bobbin")     return "bobbin_" + std::to_string(index);
    throw std::runtime_error("mesh: no FEM region defined for role '" + role + "'");
}

// A functional core gap: axial centre yc [m], length len [m], and the x of the column it sits in,
// xc [m]. Extracted from the ENRICHED core gapping (subtractive/additive only; residual was retyped
// to subtractive above). Used to drive adaptive gap mesh refinement (>= N_gap elements across each
// gap + fringing band).
//
// ABT #1305 (2026-09-21): xc was DROPPED here, so the 2D refinement box sat on the centre post for
// every gap and an outer-leg gap was never refined -- on PQ 20/16 the two outer-leg 5 um residual
// gaps got 4 elements of ~52 um while the centre gap was at 0.8 um, and an Ampere loop in the empty
// strip beside the outer leg circulated +1.82 A where it must be 0 (mkf-57). The 3D level-set path
// already carried xc; only this extractor missed it. A gap MKF has not placed is now an error: with
// several gaps, a missing coordinate is not "the centre post", it is an unknown.
std::vector<CoreGap> extract_core_gaps(const json& enr) {
    std::vector<CoreGap> gaps;
    if (!enr.contains("core") || !enr["core"].contains("functionalDescription")) return gaps;
    const auto& cfd = enr["core"]["functionalDescription"];
    if (!cfd.contains("gapping")) return gaps;
    for (const auto& g : cfd["gapping"]) {
        const std::string ty = g.value("type", std::string());
        if (ty != "subtractive" && ty != "additive") continue;
        const double len = g.value("length", 0.0);
        if (!(len > 0.0)) continue;
        if (!(g.contains("coordinates") && g["coordinates"].is_array() && g["coordinates"].size() >= 2 &&
              !g["coordinates"][0].is_null() && !g["coordinates"][1].is_null()))
            throw std::runtime_error("extract_core_gaps: a " + ty + " gap of " + std::to_string(len * 1e3) +
                                     " mm has no (x, y) coordinates in the enriched core -- MKF did not place it, "
                                     "so there is no column to refine");
        gaps.push_back({g["coordinates"][1].get<double>(), len, g["coordinates"][0].get<double>()});
    }
    return gaps;
}

// Central-column radial half-extent [m] from the enriched processedDescription (0 if unavailable).
// Bounds the gap refinement region radially (the gap + its fringing live around the centre leg).
double central_column_half_width(const json& enr) {
    double r = 0.0;
    if (!enr.contains("core") || !enr["core"].contains("processedDescription")) return r;
    const auto& pd = enr["core"]["processedDescription"];
    if (!pd.contains("columns")) return r;
    for (const auto& col : pd["columns"]) {
        if (col.value("type", std::string()) != "central") continue;
        const double w = col.value("width", 0.0), d = col.value("depth", 0.0);
        r = std::max({r, 0.5 * w, 0.5 * d});
    }
    return r;
}

// ABT #1169: the public entry points (declared in MasMesher.hpp) into the classifier above.
std::string classify_solid_name(const std::string& name, const SolidNameContext& ctx,
                                std::string& winding, int& index) {
    return classify(name, winding, index, ctx);
}

std::string fem_region_for_role(const std::string& role, const std::string& winding, int index,
                                const std::string& turnRegion) {
    return region_for_role(role, winding, index, turnRegion);
}

}  // namespace mvb::mesh
