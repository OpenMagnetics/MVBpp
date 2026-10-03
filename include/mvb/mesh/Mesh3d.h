#pragma once
// Moved verbatim from OMFEM src/meshing/MasMesher.cpp (ABT #1588, step 2): namespace omfem -> mvb::mesh,
// MasMeshOptions -> MeshOptions, file-local helpers shared with OMFEM's 2D mesher made external. No logic changed.

#include <string>

#include <nlohmann/json.hpp>

#include "mvb/mesh/MeshOptions.h"

namespace mvb::mesh {

// Build the 3D tetrahedral mesh for one MAS magnetic: MVB++ buildAllNamed (real
// turns, honouring opt.polygon_segments / opt.paint_coating / opt.symmetry) ->
// gmsh C++ (air box, fragment, region tags core/winding/air + 'outer') -> .msh
// (v2.2). Pure C++; no Python. Returns the written path. Throws on failure.
std::string mesh3d_from_mas(nlohmann::json magnetic, const MeshOptions& opt);

}  // namespace mvb::mesh
