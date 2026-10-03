#pragma once
// Moved verbatim from OMFEM src/meshing/MasMesher.cpp (ABT #1588, step 2): namespace omfem -> mvb::mesh,
// MasMeshOptions -> MeshOptions, file-local helpers shared with OMFEM's 2D mesher made external. No logic changed.

#include <string>

namespace mvb::mesh {

struct MeshOptions {
    std::string out_msh;              // output path (gmsh v2.2)
    bool   axisymmetric   = true;     // round-post -> axisymmetric (r,z) +X half
    // E-type cores -> planar (Cartesian) model: cut the FULL XY cross-section (both
    // lateral legs, centre leg solid through x=0) and solve it as a 2D planar problem
    // extruded by the core depth. No symmetry reduction, no hand-applied factor -- the
    // geometry is the complete MVB++ section. Mutually exclusive with axisymmetric.
    bool   planar         = false;
    int    polygon_segments = 16;     // MVB++ curve faceting
    double air_margin     = 1.5;      // air box = component bbox * this
    double air_target     = 1.0e-3;   // max element size in air [m]
    double core_target    = 3.0e-4;   // min element size [m]
    double conductor_target = 0.0;    // >0: refine near conductors to this size; 0 (default):
                                      // AUTO -- derived from the enriched layout's real pairwise
                                      // turn clearances (chord bound sqrt(4 r g)); <0: disable
                                      // conductor refinement entirely
                                      // [m] (set ~skin_depth/3 for AC eddy accuracy)
    // Electrostatic (stray-capacitance) mesh: split each turn face into an inner
    // copper electrode (at the wire conducting diameter) and the surrounding enamel
    // coating annulus (eps_r), tagged per turn as copper_<winding>_<index> /
    // coating_<winding>_<index>. Requires the MVB++ turns to be drawn at the wire
    // OUTER (insulated) diameter, which they are.
    bool   electrostatic  = false;
    double coating_target = 0.0;      // if >0, max element size inside the coating annulus [m]
    // Magnetostatic/harmonic conducting-diameter conductor: split each turn so the
    // CONDUCTING (copper) region carries the current and the enamel annulus is air.
    // Fixes the eddy/DC loss being computed over the outer (insulated) diameter.
    // The inner copper is tagged turn_<w>_<t> (the solver's normal conductor name);
    // the annulus becomes air. Mutually exclusive with electrostatic.
    bool   conductor_conducting_diameter = false;

    // ---- 3D meshing (tetrahedral) -------------------------------------------
    // When set, mesh3d_from_mas builds the FULL 3D solids (real turns) and tets
    // them, producing volume physical groups core / winding / air plus the
    // far-field boundary surface 'outer' that the 3D solvers (omfem_3d*) expect.
    // WIRE DIMENSIONS -- the one setting deciding what size the conductors are drawn at,
    // everywhere, with no per-path exceptions (Alf, 2026-08-06):
    //   false (DEFAULT) = CONDUCTING: the bare copper cross-section -- diameter for round
    //                     wire, conducting width x height for rectangular/planar;
    //   true            = OUTER: the insulated envelope (coating included), same shape rule.
    // Conducting is what winding-loss FEM needs (the eddy currents live in the copper);
    // outer is for electrostatics/clearance visualisation. Tools may set it from the
    // OMFEM_WIRE_DIM env var ("conducting" | "outer"); any other value must throw.
    bool        paint_coating = false;
    // MVB++ symmetry reduction: "none"/"half"/"quarter"/"eighth" ("auto" -> none).
    // analyze_symmetry tests the REAL solids (turns included) and keeps only the
    // valid coordinate planes. NOTE: a reduced domain needs the solver to impose
    // matching symmetry BCs on the cut planes; until that exists keep "none".
    std::string symmetry = "none";
};

}  // namespace mvb::mesh
