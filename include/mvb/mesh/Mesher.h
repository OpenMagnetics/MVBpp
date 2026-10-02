#pragma once
// mvbpp_mesh: MAS -> FEM mesh, for any solver (ABT #1588).
//
// The meshing stage moves here from OMFEM (MasMesher, SizeField, MeshRecipe, StepAudit, the MMG
// skin-layer pass) so MVB++ meshes what it draws and OMFEM keeps only the MFEM simulation. It is
// built on libmvb++ (geometry) with gmsh and MMG, in every MVB++ build: native, pybind and WASM.
//
// THE CONTRACT. A caller hands over a MAS magnetic and a MESH RECIPE and gets back a finished
// mesh plus its sidecars. Everything physical (skin depth, conductor target size, layer count and
// thickness, drive) arrives in the recipe: the solver decides those. Everything geometric (gap
// and clearance refinement, halos, CAD edges, ports, region names) MVB++ derives itself from the
// geometry it built. Unknown recipe keys throw; a knob left null takes its documented default,
// and the effective recipe is returned with every mesh, so a mesh always records how it was made.
//
// WHILE THE MOVE IS UNDER WAY the recipe knobs are the ones OMFEM's MeshRecipe table names, with
// the same keys, units and defaults, and the moved code must reproduce OMFEM's reference meshes
// node-for-node and tet-for-tet (gmsh fork om-4.15.2, threads as pinned in the recipe). The gate
// compares PARSED meshes, not file bytes (OMFEM has two msh2 writers today, gmsh's and the skin
// layer's own, formatted differently): node ids with exact coordinates, element ids, types, node
// lists, physical and elementary tags, all in order. Sidecars stay byte-identical apart from the
// recipe's provenance (paths, cwd, timings). Any change of output is its own, separately gated step.
//
// This header is the API sketch for review; nothing implements it yet.

#include <nlohmann/json.hpp>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace mvb::mesh {

// ---- The finished mesh, independent of the mesher that made it -------------------------------
// One neutral model for gmsh's output AND for the MMG skin-layer pass (which is not gmsh's
// in-memory model). Every exporter reads only this, so every format sees the same mesh.
enum class ElementType : std::uint8_t { Point, Line2, Tri3, Quad4, Tet4, Pyramid5, Prism6, Hex8 };

struct ElementBlock {
    ElementType type;
    int region = 0;                          // physical tag (Region::tag)
    int entity = 0;                          // elementary tag: the geometric entity meshed
    std::vector<std::int64_t> ids;           // element ids, one per element
    std::vector<std::int64_t> nodes;         // node ids, nodesPerElement(type) per element
};

struct Region {
    int tag = 0;
    int dimension = 0;                       // 3 volume, 2 boundary/port, 1 CAD edge
    std::string name;                        // "core", "winding_<w>", "air", "outer", "port_..."
};

struct Mesh {
    std::string lengthUnit = "m";            // node coordinates are in metres, always
    std::vector<std::int64_t> nodeIds;       // gmsh numbering, kept so the identity gate can diff
    std::vector<double> xyz;                 // 3 per node, metres
    std::vector<Region> regions;
    std::vector<ElementBlock> blocks;
};

// ---- Meshing ---------------------------------------------------------------------------------
struct MeshResult {
    Mesh mesh;
    nlohmann::json recipe;                   // the EFFECTIVE recipe (defaults filled in, provenance)
    // Sidecars, by the file suffix OMFEM writes them under today (".path", ".cadedges",
    // ".leads.json", ".recipe.json"), byte-for-byte what is written to disk.
    std::map<std::string, std::string> sidecars;
};

// Throws on any malformed MAS, unknown recipe key, or a mesh with an inverted or zero-volume
// element -- never returns a degraded mesh.
MeshResult meshMagnetic(const nlohmann::json& magnetic, const nlohmann::json& recipe);

// The recipe schema: every knob with section, key, unit, default and one documentation line.
nlohmann::json recipeTemplate();

// What THIS library was built from, for a mesh's provenance (ABT #1592): the MVB++, MKF and MAS
// commits compiled in (a "-dirty" suffix when the checkout had uncommitted changes), and the
// gmsh and MMG versions linked. The gmsh version is asked of the linked gmsh at run time, never
// read from a checkout. Before/after builds of a move step must report identical values.
nlohmann::json buildRevision();

// ---- Export ----------------------------------------------------------------------------------
// format: "msh2" (gmsh 2.2, what OMFEM writes today: the identity-gate reference), "msh4", "bdf" (Nastran bulk data: large-field
// GRID* only, PSOLID/MAT1 or PSHELL per region with the region name, throws on an element with no
// region), and the other writers as they gain a round-trip gate. `unit` scales the coordinates
// written ("m" or "mm") and is stated in the file. Every writer is gated by a round trip: read
// back, node and element counts per region, coordinates within the format's precision, tet
// volumes recomputed with the same sign and no new |V| < 1e-18 m^3. Text formats write
// coordinates round-trip exact (%.17g), never fewer digits.
std::string exportMesh(const Mesh& mesh, const std::string& format, const std::string& unit);

}  // namespace mvb::mesh
