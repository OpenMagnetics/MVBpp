#pragma once
// Moved verbatim from OMFEM src/meshing/MasMesher.cpp at 59cc050 (ABT #1588, step 6): namespace omfem ->
// mvb::mesh, MasMeshOptions -> MeshOptions, mesh_from_mas -> mesh2d_from_mas. No logic changed.

#include <string>

#include <nlohmann/json.hpp>

#include "mvb/mesh/MeshOptions.h"

namespace mvb::mesh {

// Build the 2D mesh for one MAS magnetic: the MVB++ section (axisymmetric (r,z) +X half, or the
// full planar XY cut of an E-type core, as opt says) -> gmsh (air, region tags, gap refinement)
// -> .msh (v2.2) at opt.out_msh. The caller applies its modelling choices to the magnetic first
// (OMFEM: retype_residual_gaps). Returns the written path. Throws on failure.
std::string mesh2d_from_mas(nlohmann::json magnetic, const MeshOptions& opt);

}  // namespace mvb::mesh
