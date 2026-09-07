#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"
#include "constructive_models/Wire.h"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BOPAlgo_GlueEnum.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GC_MakeCircle.hxx>
#include <gp_Circ.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Ax1.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <Standard_Failure.hxx>
#include <StdFail_NotDone.hxx>
#include <iostream>
#include <cmath>
#include <numbers>
#include <map>
#include <vector>
#include <string>
#include <cstdlib>

namespace mvb {

// --- GLOBAL turn-part fusing toggle (see TurnBuilder.h) --------------------------------------
// Default OFF (fast compound for STL/STEP visualisation). The FEM meshing path turns it ON so each
// turn is a single clean solid (overlapping compound sub-solids otherwise crash tet meshing).
static bool s_fuseTurnParts = [](){ const char* e = std::getenv("MVB_FUSE_TURNS");
                                    return e && *e && std::string(e) != "0"; }();
void TurnBuilder::setFuseTurnParts(bool on) { s_fuseTurnParts = on; }
bool TurnBuilder::fuseTurnParts() { return s_fuseTurnParts; }

// Fuse all SOLIDs in a shape (a turn's compound of segments+corners) into one. The sub-solids
// genuinely OVERLAP at the corners (a straight cylinder runs into the corner torus), so this is a
// real boolean union (no glue). Returns the input unchanged if there is <=1 solid or the fuse
// fails (a warning is emitted; an unfused compound is still valid geometry for callers that don't
// mesh). One boolean per adjacent pair -> O(nparts) per turn.
static TopoDS_Shape fuse_turn_solids(const TopoDS_Shape& shape) {
    std::vector<TopoDS_Shape> solids;
    for (TopExp_Explorer ex(shape, TopAbs_SOLID); ex.More(); ex.Next()) solids.push_back(ex.Current());
    if (solids.size() <= 1) return shape;
    TopoDS_Shape acc = solids[0];
    for (std::size_t i = 1; i < solids.size(); ++i) {
        try {
            BRepAlgoAPI_Fuse f(acc, solids[i]);
            if (!f.IsDone()) { std::cerr << "WARN TurnBuilder: turn-part fuse not done; returning unfused compound\n"; return shape; }
            acc = f.Shape();
        } catch (const Standard_Failure& e) {
            std::cerr << "WARN TurnBuilder: turn-part fuse threw (" << e.GetMessageString() << "); returning unfused compound\n";
            return shape;
        }
    }
    return acc;
}
static inline TopoDS_Shape maybe_fuse_turn(const TopoDS_Shape& s) {
    return s_fuseTurnParts ? fuse_turn_solids(s) : s;
}

// Collapse a dimensionWithTolerance to a scalar with MKF's preference order
// (nominal -> (min+max)/2 -> max -> min) via OpenMagnetics::resolve_dimensional_values,
// but THROW when the tolerance object is empty — the MKF resolver silently returns 0.
static double resolve_dim_or_throw(const MAS::DimensionWithTolerance& dwt, const std::string& what) {
    if (!dwt.get_nominal() && !dwt.get_maximum() && !dwt.get_minimum()) {
        throw std::runtime_error(what + " has an empty dimensionWithTolerance "
                                 "(no nominal/minimum/maximum)");
    }
    return OpenMagnetics::resolve_dimensional_values(dwt);
}

static double get_wire_diameter(const MAS::Wire& wire) {
    const std::string wireName = wire.get_name().value_or("<unnamed>");
    if (auto cd = wire.get_conducting_diameter()) {
        return resolve_dim_or_throw(*cd, "wire '" + wireName + "' conductingDiameter");
    }
    if (auto od = wire.get_outer_diameter()) {
        return resolve_dim_or_throw(*od, "wire '" + wireName + "' outerDiameter");
    }
    throw std::runtime_error("get_wire_diameter: wire '" + wireName
                             + "' has neither conductingDiameter nor outerDiameter");
}

// Copper (conducting) cross-section of the wire — the geometry FEM winding-loss
// extraction must mesh, because the eddy/proximity currents and the skin-depth
// boundary live in the copper, not the insulation/serving. This intentionally
// ignores turn.get_dimensions() (those carry the OUTER/insulation footprint).
//
// For LITZ the meshable "copper" is the bare bundle treated as a single solid
// conductor: its diameter is computed from the strand geometry by MKF
// (Wire::get_outer_diameter_bare_litz), the single source of truth for that
// formula. Throws (no silent fallback) when the data needed is missing.
static std::pair<double, double> get_conducting_dimensions(const MAS::Wire& wire) {
    const std::string wireName = wire.get_name().value_or("<unnamed>");
    switch (wire.get_type()) {
        case MAS::WireType::ROUND: {
            auto cd = wire.get_conducting_diameter();
            if (!cd) {
                throw std::runtime_error(
                    "get_conducting_dimensions: ROUND wire '" + wireName
                    + "' has no conductingDiameter (required to draw the copper cross-section)");
            }
            double d = resolve_dim_or_throw(*cd, "ROUND wire '" + wireName + "' conductingDiameter");
            return {d, d};
        }
        case MAS::WireType::RECTANGULAR:
        case MAS::WireType::PLANAR: {
            auto cw = wire.get_conducting_width();
            auto ch = wire.get_conducting_height();
            if (!cw || !ch) {
                throw std::runtime_error(
                    "get_conducting_dimensions: RECTANGULAR/PLANAR wire '" + wireName
                    + "' has no conductingWidth / conductingHeight");
            }
            return {resolve_dim_or_throw(*cw, "wire '" + wireName + "' conductingWidth"),
                    resolve_dim_or_throw(*ch, "wire '" + wireName + "' conductingHeight")};
        }
        case MAS::WireType::LITZ: {
            OpenMagnetics::Wire omWire(wire);
            if (!omWire.get_number_conductors()) {
                throw std::runtime_error(
                    "get_conducting_dimensions: LITZ wire '" + wireName
                    + "' has no numberConductors (required for the bare-bundle diameter)");
            }
            int numberConductors = static_cast<int>(omWire.get_number_conductors().value());
            auto strand = omWire.resolve_strand();
            double strandConductingDiameter = resolve_dim_or_throw(
                strand.get_conducting_diameter(),
                "LITZ wire '" + wireName + "' strand conductingDiameter");
            auto strandCoating = OpenMagnetics::Wire::resolve_coating(strand);
            if (!strandCoating || !strandCoating->get_grade()) {
                throw std::runtime_error(
                    "get_conducting_dimensions: LITZ wire '" + wireName
                    + "' strand has no coating grade (required for the bare-bundle diameter)");
            }
            int grade = strandCoating->get_grade().value();
            auto standard = omWire.get_standard().value_or(MAS::WireStandard::IEC_60317);
            // Bare bundle outer diameter, treated as one solid copper conductor.
            double bundleDiameter = OpenMagnetics::Wire::get_outer_diameter_bare_litz(
                strandConductingDiameter, numberConductors, grade, standard);
            // PHYSICAL INVARIANT: the conducting envelope can never exceed the wire's OUTER envelope
            // — the bundle physically fits inside the catalogued wire. The bare-litz packing formula
            // is an ESTIMATE and can come out a few um larger than the wire's own outerDiameter
            // (e.g. Litz 20x0.1 Unserved: estimate 644um vs catalog outer 639um); MKF's turn blocking
            // spaces turns by the catalog outer, so an oversized estimate makes even perfectly placed
            // turns register copper overlaps. The catalog outer is authoritative for the envelope.
            if (auto od = wire.get_outer_diameter()) {
                double outer = resolve_dim_or_throw(*od, "LITZ wire '" + wireName + "' outerDiameter");
                if (bundleDiameter > outer) bundleDiameter = outer;
            }
            return {bundleDiameter, bundleDiameter};
        }
        case MAS::WireType::FOIL: {
            // ABT #208 (2026-08-23): a foil turn IS a swept thin rectangle -- conductingWidth is
            // the foil thickness (radial), conductingHeight the sheet height (axial), exactly the
            // frame the RECTANGULAR path already sweeps (dimensions[0]=width=RADIAL,
            // dimensions[1]=height=AXIAL). The 2KW-IH user MAS carries 0.4 x 33 mm, an 82:1
            // section, and the rect profile/revolve handles any aspect ratio. The former throw
            // here declared FOIL "not a swept section", which was true of the never-written
            // dedicated sheet builder but not of the geometry: per-turn drawing needs nothing a
            // rectangle does not have. REAL WINDING still refuses foil explicitly
            // (ConductorBuilder: a continuous foil conductor is a genuinely different
            // construction); this covers the ideal per-turn path the 3D viewer uses.
            auto cw = wire.get_conducting_width();
            auto ch = wire.get_conducting_height();
            if (!cw || !ch) {
                throw std::runtime_error(
                    "get_conducting_dimensions: FOIL wire '" + wireName
                    + "' has no conductingWidth / conductingHeight");
            }
            return {resolve_dim_or_throw(*cw, "FOIL wire '" + wireName + "' conductingWidth"),
                    resolve_dim_or_throw(*ch, "FOIL wire '" + wireName + "' conductingHeight")};
        }
        default:
            throw std::runtime_error(
                "get_conducting_dimensions: unsupported wire type for wire '" + wireName + "'");
    }
}

std::pair<double, double> conducting_dimensions(const MAS::Wire& wire) {
    return get_conducting_dimensions(wire);
}

// Resolve the wire cross-section used to build a turn solid.
//   paintCoating == true  → OUTER (insulation) footprint  [default, unchanged]
//   paintCoating == false → CONDUCTING (copper) footprint  [FEM winding-loss]
static std::pair<double, double> get_wire_dimensions(const MAS::Wire& wire, const MAS::Turn& turn,
                                                     bool paintCoating) {
    if (!paintCoating) {
        return get_conducting_dimensions(wire);
    }
    auto td = turn.get_dimensions();
    if (td && td->size() >= 2) {
        return {(*td)[0], (*td)[1]};
    }
    const std::string wireName = wire.get_name().value_or("<unnamed>");
    auto ow = wire.get_outer_width();
    auto oh = wire.get_outer_height();
    if (ow && oh) {
        return {resolve_dim_or_throw(*ow, "wire '" + wireName + "' outerWidth"),
                resolve_dim_or_throw(*oh, "wire '" + wireName + "' outerHeight")};
    }
    auto cw = wire.get_conducting_width();
    auto ch = wire.get_conducting_height();
    if (cw && ch) {
        return {resolve_dim_or_throw(*cw, "wire '" + wireName + "' conductingWidth"),
                resolve_dim_or_throw(*ch, "wire '" + wireName + "' conductingHeight")};
    }
    double d = get_wire_diameter(wire);
    return {d, d};
}

static bool is_rectangular_wire(const MAS::Wire& wire) {
    // FOIL dispatches the rectangular profile too (ABT #208): thickness x height IS its section.
    return wire.get_type() == MAS::WireType::RECTANGULAR ||
           wire.get_type() == MAS::WireType::PLANAR ||
           wire.get_type() == MAS::WireType::FOIL;
}

static TopoDS_Face build_circle_profile(const gp_Ax2& plane, double radius, int segments = 16) {
    gp_Pnt origin = plane.Location();
    gp_Dir dx = plane.XDirection();
    gp_Dir dy = plane.YDirection();

    if (segments <= 0) {
        gp_Circ circ(plane, radius);
        BRepBuilderAPI_MakeEdge edge(circ);
        BRepBuilderAPI_MakeWire wire(edge);
        BRepBuilderAPI_MakeFace face(wire.Wire());
        return face.Face();
    }

    BRepBuilderAPI_MakePolygon poly;
    double offset = std::numbers::pi / segments;
    for (int i = 0; i < segments; ++i) {
        double angle = 2.0 * std::numbers::pi * i / segments + offset;
        gp_Vec vx(dx);
        gp_Vec vy(dy);
        vx.Scale(radius * std::cos(angle));
        vy.Scale(radius * std::sin(angle));
        gp_Pnt p = origin.Translated(vx).Translated(vy);
        poly.Add(p);
    }
    poly.Close();
    BRepBuilderAPI_MakeFace face(BRepBuilderAPI_MakeWire(poly.Wire()).Wire());
    return face.Face();
}

static TopoDS_Face build_rect_profile(const gp_Ax2& plane, double width, double height) {
    gp_Pnt origin = plane.Location();
    gp_Dir dx = plane.XDirection();
    gp_Dir dy = plane.YDirection();
    double hw = width / 2.0;
    double hh = height / 2.0;

    gp_Pnt p1 = origin.Translated(gp_Vec(dx).Scaled(-hw)).Translated(gp_Vec(dy).Scaled(-hh));
    gp_Pnt p2 = origin.Translated(gp_Vec(dx).Scaled(-hw)).Translated(gp_Vec(dy).Scaled( hh));
    gp_Pnt p3 = origin.Translated(gp_Vec(dx).Scaled( hw)).Translated(gp_Vec(dy).Scaled( hh));
    gp_Pnt p4 = origin.Translated(gp_Vec(dx).Scaled( hw)).Translated(gp_Vec(dy).Scaled(-hh));

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(p1);
    poly.Add(p2);
    poly.Add(p3);
    poly.Add(p4);
    poly.Close();
    BRepBuilderAPI_MakeFace face(poly.Wire());
    return face.Face();
}

static TopoDS_Shape build_concentric_round_column_turn(double turn_radius, double wire_radius,
                                                         double height_pos, bool rect_wire,
                                                         double wire_width, double wire_height,
                                                         int segments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                                         int revolution_segments = DEFAULT_WIRE_REVOLUTION_SEGMENTS) {
    if (rect_wire) {
        // Rectangular wire on round column: revolve a rectangle face 360° around Y axis.
        // This matches MVB.js behavior and avoids BRepOffsetAPI_MakePipe degeneracies.
        double hw = wire_width / 2.0;
        double hh = wire_height / 2.0;
        gp_Pnt p1(turn_radius - hw, height_pos - hh, 0.0);
        gp_Pnt p2(turn_radius + hw, height_pos - hh, 0.0);
        gp_Pnt p3(turn_radius + hw, height_pos + hh, 0.0);
        gp_Pnt p4(turn_radius - hw, height_pos + hh, 0.0);

        BRepBuilderAPI_MakePolygon poly;
        poly.Add(p1);
        poly.Add(p2);
        poly.Add(p3);
        poly.Add(p4);
        poly.Close();

        BRepBuilderAPI_MakeFace faceMaker(poly.Wire());
        if (!faceMaker.IsDone()) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeFace failed for rectangle profile\n";
            return TopoDS_Shape();
        }
        TopoDS_Face profile = faceMaker.Face();

        gp_Ax1 rev_axis(gp_Pnt(0.0, height_pos, 0.0), gp_Dir(0, 1, 0));
        BRepPrimAPI_MakeRevol revol(profile, rev_axis, 2.0 * std::numbers::pi - 1e-6);
        if (!revol.IsDone()) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol failed for rect wire turn\n";
            return TopoDS_Shape();
        }
        try {
            return revol.Shape();
        } catch (const Standard_Failure& e) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol threw Standard_Failure: " << e.GetMessageString() << "\n";
            return TopoDS_Shape();
        } catch (...) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol threw unknown exception\n";
            return TopoDS_Shape();
        }
    }

    // Round wire on round column.
    // segments = 0         → exact torus (analytic surfaces, legacy)
    // segments > 0, rev=0  → revolve a P-gon cross-section (analytic faces of
    //                         revolution — also CYLINDRICAL_SURFACE in STEP)
    // segments > 0, rev > 0 → N-sector ThruSections loft (fully planar,
    //                         needed for node-conformal Gmsh meshes)
    if (segments <= 0) {
        gp_Pnt torus_center(0.0, height_pos, 0.0);
        gp_Ax2 torus_axis(torus_center, gp_Dir(0, 1, 0), gp_Dir(1, 0, 0));
        try {
            return BRepPrimAPI_MakeTorus(torus_axis, turn_radius, wire_radius).Shape();
        } catch (const Standard_Failure& e) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeTorus threw Standard_Failure: " << e.GetMessageString() << "\n";
            return TopoDS_Shape();
        }
    }

    if (revolution_segments > 0) {
        TopoDS_Shape ring = build_polygon_ring(turn_radius, wire_radius,
                                                height_pos, segments,
                                                revolution_segments);
        if (!ring.IsNull()) return ring;
        std::cerr << "WARN build_concentric_round_column_turn: polygon ring "
                     "failed, falling back to MakeRevol\n";
    }

    gp_Ax2 profile_plane(gp_Pnt(turn_radius, height_pos, 0.0),
                          gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    TopoDS_Face profile = build_circle_profile(profile_plane, wire_radius, segments);
    gp_Ax1 rev_axis(gp_Pnt(0.0, height_pos, 0.0), gp_Dir(0, 1, 0));
    BRepPrimAPI_MakeRevol revol(profile, rev_axis, 2.0 * std::numbers::pi - 1e-6);
    if (!revol.IsDone()) {
        std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol failed for round wire turn\n";
        return TopoDS_Shape();
    }
    try {
        return revol.Shape();
    } catch (const Standard_Failure& e) {
        std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol threw Standard_Failure: " << e.GetMessageString() << "\n";
        return TopoDS_Shape();
    }
}

static TopoDS_Shape build_concentric_rect_column_turn(double radial_pos, double wire_radius,
                                                       double height_pos,
                                                       double half_col_width, double half_col_depth,
                                                       bool rect_wire,
                                                       double wire_width, double wire_height) {
    // The straight runs sit at exactly the MAS radial position: coordinates[0] is the
    // cross-section centre on the X faces, and the same clearance past the column is used
    // on the Z faces. This offset is MKF data and must never be adjusted here.
    double clearance = radial_pos - half_col_width;
    if (clearance <= 0.0) {
        throw std::runtime_error(
            "build_concentric_rect_column_turn: turn radial position " + std::to_string(radial_pos)
            + " m lies inside the column (half-width " + std::to_string(half_col_width)
            + " m) — inconsistent MAS turn/bobbin data");
    }

    double wire_x_pos = half_col_width + clearance;   // == radial_pos, exactly
    double wire_z_pos = half_col_depth + clearance;
    double y = height_pos;

    // Corner bend radius. Normally equal to the clearance (arc axis at the column corner).
    // When the wire sits closer to the column than half a wire the quarter-torus would
    // self-intersect (major < minor), so grow ONLY the corner radius and pull the arc axis
    // inward along the corner diagonal — the straights stay tangent at their exact MAS
    // positions and the arc still meets them. Geometrically infeasible data throws instead
    // of silently repositioning the turn.
    // Minimum feasible corner radius: the swept corner self-intersects when the arc radius
    // is smaller than the profile's RADIAL half-extent — wire_width/2 for rectangular/foil
    // wire (the width lies radially in the racetrack; a foil bends about its thin axis), the
    // wire radius for round wire. max(w,h)/2 would wrongly reject foil windings whose HEIGHT
    // is large but which bend trivially.
    double turn_turn_radius = clearance;
    double corner_inset = 0.0;
    double min_bend = (rect_wire ? wire_width / 2.0 : wire_radius) * 1.02;
    if (turn_turn_radius < min_bend) {
        corner_inset = min_bend - clearance;
        turn_turn_radius = min_bend;
        if (corner_inset > std::min(half_col_width, half_col_depth)) {
            throw std::runtime_error(
                "build_concentric_rect_column_turn: minimum bend radius " + std::to_string(min_bend)
                + " m cannot be accommodated (clearance " + std::to_string(clearance)
                + " m, column half-dims " + std::to_string(half_col_width) + "/"
                + std::to_string(half_col_depth) + " m) — adjacent corner arcs would cross");
        }
    }
    double seg_x = half_col_width - corner_inset;   // straight-run half-extent along X
    double seg_z = half_col_depth - corner_inset;   // straight-run half-extent along Z

    gp_Pnt pts[8] = {
        gp_Pnt(+seg_x,      y, +wire_z_pos),
        gp_Pnt(-seg_x,      y, +wire_z_pos),
        gp_Pnt(-wire_x_pos, y, +seg_z),
        gp_Pnt(-wire_x_pos, y, -seg_z),
        gp_Pnt(-seg_x,      y, -wire_z_pos),
        gp_Pnt(+seg_x,      y, -wire_z_pos),
        gp_Pnt(+wire_x_pos, y, -seg_z),
        gp_Pnt(+wire_x_pos, y, +seg_z),
    };

    std::pair<gp_Pnt, gp_Dir> corners[4] = {
        {gp_Pnt(-seg_x, y, +seg_z), gp_Dir(0, 0, 1)},
        {gp_Pnt(-seg_x, y, -seg_z), gp_Dir(-1, 0, 0)},
        {gp_Pnt(+seg_x, y, -seg_z), gp_Dir(0, 0, -1)},
        {gp_Pnt(+seg_x, y, +seg_z), gp_Dir(1, 0, 0)},
    };

    // For round wire, build the turn manually from straight cylinders and corner quarter-tori.
    // This avoids BRepOffsetAPI_MakePipe failures on closed paths in OCCT 7.9.
    if (!rect_wire) {
        BRep_Builder builder;
        TopoDS_Compound compound;
        builder.MakeCompound(compound);

        // 4 straight segments: extrude a circle along each straight edge
        int straight_pairs[4][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}};
        for (int i = 0; i < 4; ++i) {
            gp_Pnt p1 = pts[straight_pairs[i][0]];
            gp_Pnt p2 = pts[straight_pairs[i][1]];
            gp_Vec vec(p1, p2);
            double len = vec.Magnitude();
            if (len < 1e-12) continue;

            // Cylinder along the segment, centered at p1, axis = vec
            gp_Ax2 cylAx2(p1, gp_Dir(vec));
            TopoDS_Shape cyl = BRepPrimAPI_MakeCylinder(cylAx2, wire_radius, len).Shape();
            builder.Add(compound, cyl);
        }

        // 4 corner quarter-tori: revolve a circle around the Y axis at each corner center
        int corner_incoming[4] = {0, 2, 4, 6}; // start index of straight segment before each corner
        for (int i = 0; i < 4; ++i) {
            auto [c_center, c_xref] = corners[i];
            gp_Pnt circ_center = c_center.Translated(gp_Vec(c_xref).Scaled(turn_turn_radius));
            gp_Dir incoming_dir = gp_Dir(gp_Vec(pts[corner_incoming[i]], pts[corner_incoming[i] + 1]));
            gp_Ax2 circ_ax2(circ_center, incoming_dir, gp_Dir(0, 1, 0));
            gp_Circ circ(circ_ax2, wire_radius);
            TopoDS_Edge circ_edge = BRepBuilderAPI_MakeEdge(circ).Edge();
            TopoDS_Wire circ_wire = BRepBuilderAPI_MakeWire(circ_edge).Wire();
            TopoDS_Face circ_face = BRepBuilderAPI_MakeFace(circ_wire).Face();

            gp_Ax1 rev_axis(c_center, gp_Dir(0, 1, 0));
            // All corners sweep -90° around +Y (clockwise when viewed from +Y)
            BRepPrimAPI_MakeRevol revol(circ_face, rev_axis, -std::numbers::pi / 2.0);
            if (revol.IsDone()) {
                builder.Add(compound, revol.Shape());
            }
        }

        // Return compound directly — fusing is unnecessary for STL visualization
        // and BRepAlgoAPI_Fuse per-turn is prohibitively slow (7 fuses × N turns).
        return compound;
    }

    // Rectangular wire: build the turn from 4 axis-aligned boxes + 4 corner
    // quarter-swept rectangles. Mirrors the MVB.js approach in
    // replicadBuilder.js::_makeQuarterSweptRectangle and avoids BRepOffsetAPI_MakePipe
    // entirely (it segfaults on closed planar paths in OCCT 7.9).
    double hw_wire = wire_width / 2.0;
    double hh_wire = wire_height / 2.0;

    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    // Straight segments (skipped when the corner inset consumed the whole face,
    // i.e. adjacent corner arcs meet exactly at the face midpoint).
    if (seg_x > 1e-12) {
        // Segment 0: along -X at z=+wire_z_pos
        gp_Pnt corner0(-seg_x, y - hh_wire, wire_z_pos - hw_wire);
        builder.Add(compound, BRepPrimAPI_MakeBox(corner0, 2.0 * seg_x, wire_height, wire_width).Shape());
        // Segment 2: along +X at z=-wire_z_pos
        gp_Pnt corner2(-seg_x, y - hh_wire, -wire_z_pos - hw_wire);
        builder.Add(compound, BRepPrimAPI_MakeBox(corner2, 2.0 * seg_x, wire_height, wire_width).Shape());
    }
    if (seg_z > 1e-12) {
        // Segment 1: along -Z at x=-wire_x_pos
        gp_Pnt corner1(-wire_x_pos - hw_wire, y - hh_wire, -seg_z);
        builder.Add(compound, BRepPrimAPI_MakeBox(corner1, wire_width, wire_height, 2.0 * seg_z).Shape());
        // Segment 3: along +Z at x=+wire_x_pos
        gp_Pnt corner3(wire_x_pos - hw_wire, y - hh_wire, -seg_z);
        builder.Add(compound, BRepPrimAPI_MakeBox(corner3, wire_width, wire_height, 2.0 * seg_z).Shape());
    }

    // 4 corner quarter-swept rectangles: revolve a rect face -90° around Y
    // centered at each corner. The profile sits at circ_center = c_center + c_xref * bend_radius,
    // oriented with normal = incoming_dir, local X = c_xref (radial), local Y = +Y (column axis).
    int corner_incoming[4] = {0, 2, 4, 6};
    for (int i = 0; i < 4; ++i) {
        auto [c_center, c_xref] = corners[i];
        gp_Pnt circ_center = c_center.Translated(gp_Vec(c_xref).Scaled(turn_turn_radius));
        gp_Dir incoming_dir = gp_Dir(gp_Vec(pts[corner_incoming[i]], pts[corner_incoming[i] + 1]));

        gp_Ax2 profPlane(circ_center, incoming_dir, c_xref);
        TopoDS_Face rect_face = build_rect_profile(profPlane, wire_width, wire_height);

        gp_Ax1 rev_axis(c_center, gp_Dir(0, 1, 0));
        BRepPrimAPI_MakeRevol revol(rect_face, rev_axis, -std::numbers::pi / 2.0);
        if (revol.IsDone()) {
            builder.Add(compound, revol.Shape());
        } else {
            std::cerr << "ERROR build_concentric_rect_column_turn: quarter-swept rect "
                      << i << " MakeRevol failed\n";
        }
    }

    // Return compound directly — fusing is unnecessary for STL visualization.
    return compound;
}

// Helpers for toroidal turn construction (mirrors MVB.js _createToroidalTurn).
static TopoDS_Shape rotate_about_axis(const TopoDS_Shape& s, const gp_Pnt& pt,
                                       const gp_Dir& dir, double rad) {
    if (s.IsNull() || std::abs(rad) < 1e-12) return s;
    gp_Trsf t;
    t.SetRotation(gp_Ax1(pt, dir), rad);
    return BRepBuilderAPI_Transform(s, t).Shape();
}

static TopoDS_Shape make_quarter_torus_at(double major_R, double minor_r,
                                            const gp_Pnt& center, double startAngleDeg) {
    double a = startAngleDeg * std::numbers::pi / 180.0;
    gp_Ax2 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1),
              gp_Dir(std::cos(a), std::sin(a), 0.0));
    try {
        // Match MVB.js _makeQuarterSweptRectangle pattern but with a circular
        // wire cross-section: disk profile in the XZ plane at (major_R, 0, 0)
        // (so the revolution axis Z lies IN the face's plane), revolve around
        // Z by 90°, then rotate by startAngle, then translate to center.
        gp_Ax2 circ_ax(gp_Pnt(major_R, 0.0, 0.0), gp_Dir(0, 1, 0));
        gp_Circ circ(circ_ax, minor_r);
        TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circ).Edge();
        TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge).Wire();
        TopoDS_Face face = BRepBuilderAPI_MakeFace(wire, Standard_True).Face();
        gp_Ax1 rev_axis(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
        BRepPrimAPI_MakeRevol mk(face, rev_axis, std::numbers::pi / 2.0);
        if (!mk.IsDone()) return TopoDS_Shape();
        TopoDS_Shape shape = mk.Shape();
        if (std::abs(a) > 1e-12) {
            gp_Trsf rot;
            rot.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), a);
            shape = BRepBuilderAPI_Transform(shape, rot).Shape();
        }
        return translate_shape(shape, center.X(), center.Y(), center.Z());
    } catch (const Standard_Failure& e) {
        std::cerr << "ERROR make_quarter_torus_at: " << e.GetMessageString() << "\n";
        return TopoDS_Shape();
    }
}

// Cache key for canonical toroidal-turn geometry (everything except turnRotationRad).
// Turns sharing the same key are identical shapes, differing only in Y-axis rotation.
// wire_w/wire_h/is_rect distinguish rectangular wires (which have separate width and
// height) from round ones; without them, two different rect wires with the same
// min-dimension would collide and render with each other's geometry.
struct ToroidalTurnKey {
    double wire_radius, wire_w, wire_h;
    double innerRadial, outerRadial, angleDiffRad, halfDepth;
    bool is_rect;
    bool operator<(const ToroidalTurnKey& o) const {
        auto snap = [](double v) { return std::round(v * 1e6) / 1e6; };
        if (is_rect != o.is_rect) return is_rect < o.is_rect;
        if (snap(wire_radius)  != snap(o.wire_radius))  return snap(wire_radius)  < snap(o.wire_radius);
        if (snap(wire_w)       != snap(o.wire_w))       return snap(wire_w)       < snap(o.wire_w);
        if (snap(wire_h)       != snap(o.wire_h))       return snap(wire_h)       < snap(o.wire_h);
        if (snap(innerRadial)  != snap(o.innerRadial))  return snap(innerRadial)  < snap(o.innerRadial);
        if (snap(outerRadial)  != snap(o.outerRadial))  return snap(outerRadial)  < snap(o.outerRadial);
        if (snap(angleDiffRad) != snap(o.angleDiffRad)) return snap(angleDiffRad) < snap(o.angleDiffRad);
        return snap(halfDepth) < snap(o.halfDepth);
    }
};

// Mirror of MVB.js _makeToroidalQuarterSweptRectangle: build a rectangle face
// (size wireWidth × wireHeight) at distance bendRadius from the Z axis in the
// XZ plane, rotate it by startAngleRad around Z, sweep ±π/2 around Z, then
// translate to `center`. The sign of ySign chooses CCW (+) or CW (−) sweep so
// the corner connects the right tube/radial segment of each Y-half.
static TopoDS_Shape make_toroidal_quarter_swept_rectangle(
    double bendRadius, double wireWidth, double wireHeight,
    const gp_Pnt& center, double startAngleRad, double ySign) {

    double hw = wireWidth / 2.0;
    double hh = wireHeight / 2.0;

    try {
        // 4 explicit edges → wire → face, matching MVB.js construction step-for-step.
        gp_Pnt p1(bendRadius - hw, 0, -hh);
        gp_Pnt p2(bendRadius + hw, 0, -hh);
        gp_Pnt p3(bendRadius + hw, 0,  hh);
        gp_Pnt p4(bendRadius - hw, 0,  hh);

        TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(p1, p2).Edge();
        TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(p2, p3).Edge();
        TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(p3, p4).Edge();
        TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(p4, p1).Edge();

        BRepBuilderAPI_MakeWire wireMaker;
        wireMaker.Add(e1);
        wireMaker.Add(e2);
        wireMaker.Add(e3);
        wireMaker.Add(e4);
        if (!wireMaker.IsDone()) {
            std::cerr << "ERROR make_toroidal_quarter_swept_rectangle: MakeWire failed\n";
            return TopoDS_Shape();
        }
        TopoDS_Wire rectWire = wireMaker.Wire();

        BRepBuilderAPI_MakeFace faceMaker(rectWire, Standard_True);
        if (!faceMaker.IsDone()) {
            std::cerr << "ERROR make_toroidal_quarter_swept_rectangle: MakeFace failed (status="
                      << static_cast<int>(faceMaker.Error()) << ")\n";
            return TopoDS_Shape();
        }
        TopoDS_Shape profile = faceMaker.Face();

        if (std::abs(startAngleRad) > 1e-9) {
            profile = rotate_about_axis(profile, gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), startAngleRad);
            if (profile.IsNull()) {
                std::cerr << "ERROR make_toroidal_quarter_swept_rectangle: rotate produced null\n";
                return TopoDS_Shape();
            }
        }

        gp_Ax1 revAxis(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
        BRepPrimAPI_MakeRevol revol(profile, revAxis, (std::numbers::pi / 2.0) * ySign, Standard_True);
        if (!revol.IsDone()) {
            std::cerr << "ERROR make_toroidal_quarter_swept_rectangle: MakeRevol not done\n";
            return TopoDS_Shape();
        }
        TopoDS_Shape swept = revol.Shape();
        if (swept.IsNull()) {
            std::cerr << "ERROR make_toroidal_quarter_swept_rectangle: revol.Shape() null\n";
            return TopoDS_Shape();
        }

        return translate_shape(swept, center.X(), center.Y(), center.Z());
    } catch (const Standard_Failure& e) {
        std::cerr << "ERROR make_toroidal_quarter_swept_rectangle: Standard_Failure: "
                  << e.GetMessageString()
                  << " (bendRadius=" << bendRadius << " w=" << wireWidth << " h=" << wireHeight
                  << " startAngle=" << startAngleRad << " ySign=" << ySign << ")\n";
        return TopoDS_Shape();
    } catch (const std::exception& e) {
        std::cerr << "ERROR make_toroidal_quarter_swept_rectangle: std::exception: " << e.what() << "\n";
        return TopoDS_Shape();
    } catch (...) {
        std::cerr << "ERROR make_toroidal_quarter_swept_rectangle: unknown exception\n";
        return TopoDS_Shape();
    }
}
static std::map<ToroidalTurnKey, TopoDS_Shape> s_toroidalTurnCache;

static TopoDS_Shape build_toroidal_turn(const MAS::Turn& turn, const MAS::Wire& wire,
                                         const MAS::CoreBobbinProcessedDescription& bobbin,
                                         bool paintCoating) {
    bool rect_wire = is_rectangular_wire(wire);
    auto [wire_w, wire_h] = get_wire_dimensions(wire, turn, paintCoating);
    double wire_radius = std::min(wire_w, wire_h) / 2.0;

    const std::string turnName = turn.get_name().empty() ? "<unnamed>" : turn.get_name();
    const auto& coords = turn.get_coordinates();
    if (coords.size() < 2) {
        throw std::runtime_error("build_toroidal_turn: turn '" + turnName
                                 + "' has fewer than 2 coordinates");
    }
    double cx = coords[0], cy = coords[1];
    double innerRadial = std::sqrt(cx * cx + cy * cy);
    double innerAngleDeg = std::atan2(cy, cx) * 180.0 / std::numbers::pi;

    auto addCoords = turn.get_additional_coordinates();
    const auto& wws = bobbin.get_winding_windows();
    if (wws.empty() || !wws[0].get_radial_height()) {
        throw std::runtime_error("build_toroidal_turn: toroidal bobbin has no winding-window "
                                 "radialHeight (required for the turn's face-run clearance)");
    }
    double wwRadialHeight = wws[0].get_radial_height().value();

    // The outer XY-plane crossing is MAS data (additionalCoordinates[0], written by MKF's
    // wind_toroidal_additional_turns with its own collision avoidance) — never invented here.
    if (!(addCoords && !addCoords->empty() && (*addCoords)[0].size() >= 2)) {
        throw std::runtime_error(
            "build_toroidal_turn: toroidal turn '" + turnName
            + "' has no additionalCoordinates (outer XY-plane crossing). MKF emits it when "
              "settings coil_include_additional_coordinates is enabled (the default); refusing "
              "to invent the outer crossing");
    }
    double ax = (*addCoords)[0][0], ay = (*addCoords)[0][1];
    double outerRadial = std::sqrt(ax * ax + ay * ay);
    double outerAngleDeg = std::atan2(ay, ax) * 180.0 / std::numbers::pi;

    // Special case: when the inner cartesian is at (or very near) the toroid
    // origin, atan2(0, 0) returns 0 — but the racetrack is meant to extend
    // toward the outer leg, not along +X. Inherit the outer leg's angle so
    // the angleDiff is zero and both halves are built coherently along the
    // single line connecting origin to the outer leg.
    if (innerRadial < 1e-9) {
        innerAngleDeg = outerAngleDeg;
    }

    // Azimuthal skew between the two crossings, normalized to (-pi, pi] so the outer
    // leg always takes the short way around, and the turn's placement azimuth — from
    // the MAS coordinates, NEVER from turn.rotation: MKF's angular compaction shifts
    // coordinates without updating rotation (stale; MKF ABT #186), and per the MAS
    // schema rotation is the cross-section orientation, not the placement.
    //
    // Frames: the canonical racetrack is built with the inner leg on the -X axis and the
    // toroid hole axis along +Y; buildAllNamed later counter-rotates the whole assembly by
    // -pi/2 about X, mapping build (x, y, z) -> final (x, z, -y). For the inner crossing to
    // land exactly at the MAS point (cx, cy, 0) the whole-turn rotation about +Y must be
    // beta = 180 deg - innerAngle (NOT innerAngle - 180: that mirrors the y sign — the
    // chirality bug this replaces). Same derivation makes the outer pieces rotate about the
    // GLOBAL axis by -angleDiff (see below), landing the outer leg exactly on (ax, ay, 0).
    double angleDiffRad = std::remainder(
        (outerAngleDeg - innerAngleDeg) * std::numbers::pi / 180.0, 2.0 * std::numbers::pi);
    double turnRotationRad = (180.0 - innerAngleDeg) * std::numbers::pi / 180.0;

    double halfDepth = bobbin.get_column_depth();
    double bendRadius = rect_wire ? std::max(wire_w, wire_h) / 2.0 : wire_radius;
    // The Y-clearance between a turn and the core's top/bottom face must
    // match the X-Z clearance between that turn's outer surface and the
    // core's inner hole surface. So:
    //   - layer-0 turn (against inner wall, innerRadial = B/2 - wire_radius):
    //       radial gap = 0  ⇒  Y gap = 0  ⇒  vertical tube length = C/2.
    //   - layer-1 turn (one wire OD further inward):
    //       radial gap = wire_OD  ⇒  Y gap = wire_OD  ⇒  vertical tube longer
    //       by wire_OD so the top/bottom horizontal segments clear layer 0.
    //   - a wire centered in the hole (innerRadial = 0):
    //       radial gap = B/2 - wire_radius (the whole inner hole minus the
    //       wire half) ⇒ tube length stretched to wrap the full ring.
    //
    // layerOffset = radial distance from this turn's outer surface to the
    //               core's inner hole boundary. The wire's "radial" half-
    //               dimension is wire_w/2 for rectangular wire (the face
    //               that points toward the toroidal axis) and wire_radius
    //               for round wire. Using wire_radius for rect wires would
    //               subtract the wrong amount whenever wire_w ≠ wire_h.
    // Both inner AND outer vertical tubes share the same length, so the
    // top/bottom horizontal segments stay parallel to the core's faces.
    double radialHalf = rect_wire ? wire_w / 2.0 : wire_radius;
    double layerOffset = std::max(0.0, (wwRadialHeight - innerRadial) - radialHalf);
    double radialHeight = halfDepth + layerOffset + bendRadius;

    double tubeLength = std::max(1e-7, radialHeight - bendRadius);
    double radialDistance = outerRadial - innerRadial;
    double innerX = -innerRadial;
    double outerX = innerX - radialDistance;

    // --- Exact outer-leg placement ---------------------------------------------------
    // The outer corner + tube rotate about the GLOBAL toroid axis by -angleDiff, landing
    // the outer vertical leg exactly on additionalCoordinates[0]. The run between the two
    // corners is then the straight CHORD between their radial-side faces (not a pure
    // radial segment), and each corner yaws about its own tube axis so its exit meets the
    // chord. The previous code rotated the whole outer half about the INNER leg, which
    // misplaces the outer crossing whenever the two crossing angles differ.
    // All horizontal (x, z) in the canonical frame; rotY = right-handed rotation about +Y.
    auto rotY = [](double x, double z, double ang) -> std::pair<double, double> {
        return {x * std::cos(ang) + z * std::sin(ang), -x * std::sin(ang) + z * std::cos(ang)};
    };
    const double deltaOut = -angleDiffRad;
    // Corner exit faces before yawing: inner at (innerX - b, 0); outer rotated by deltaOut.
    const double exInX = innerX - bendRadius;
    const auto [eOutX, eOutZ] = rotY(outerX + bendRadius, 0.0, deltaOut);
    const double vX = eOutX - exInX, vZ = eOutZ - 0.0;
    const double vLen = std::hypot(vX, vZ);
    double psiIn = 0.0, psiOut = 0.0;
    if (vLen > 1e-12) {
        // Yaw of the inner corner so its exit direction (-X canonically) points along the
        // chord: R_Y(psi)(-1,0) = (-cos psi, sin psi) == v/|v|.
        psiIn = std::atan2(vZ, -vX);
        // Yaw of the outer corner so its radial entry direction u points back along -v.
        const auto [uX, uZ] = rotY(1.0, 0.0, deltaOut);
        const double bX = -vX / vLen, bZ = -vZ / vLen;
        psiOut = std::atan2(uZ * bX - uX * bZ, uX * bX + uZ * bZ);
    }
    // Chord endpoints after the corners yaw (horizontal parts; y = ±radialHeight per half).
    const auto [rInX, rInZ] = rotY(-bendRadius, 0.0, psiIn);
    const double chordInX = innerX + rInX, chordInZ = rInZ;
    const auto [tOutX, tOutZ] = rotY(outerX, 0.0, deltaOut);
    const auto [rOutX, rOutZ] = rotY(eOutX - tOutX, eOutZ - tOutZ, psiOut);
    const double chordOutX = tOutX + rOutX, chordOutZ = tOutZ + rOutZ;
    const double chordLen = std::hypot(chordOutX - chordInX, chordOutZ - chordInZ);

    // Check cache: turns sharing the geometry parameters have identical
    // canonical shapes — only turnRotationRad differs.
    ToroidalTurnKey key{wire_radius, wire_w, wire_h,
                        innerRadial, outerRadial, angleDiffRad, halfDepth, rect_wire};
    auto it = s_toroidalTurnCache.find(key);

    TopoDS_Shape canonical;
    if (it != s_toroidalTurnCache.end()) {
        canonical = it->second;
    } else {
        BRep_Builder builder;
        TopoDS_Compound compound;
        builder.MakeCompound(compound);

        if (rect_wire) {
            // Rectangular toroidal turn — port of MVB.js _createToroidalTurn rect branch.
            // Each Y-half of the racetrack is 5 axis-aligned pieces:
            //   1) inner Y-tube (box wireWidth × tubeLength × wireHeight)
            //   2) inner corner (quarter swept rectangle, sweep around Z at corner center)
            //   3) radial X-tube (box radialLength × wireWidth × wireHeight) at top/bottom
            //   4) outer corner (quarter swept rectangle)
            //   5) outer Y-tube
            // Wire flat face lies tangent to the core: wireWidth tangent to the path,
            // wireHeight perpendicular (Z).
            auto build_rect_half = [&](double ySign) {
                // 1. Inner tube
                {
                    gp_Pnt corner(-wire_w / 2.0, -tubeLength / 2.0, -wire_h / 2.0);
                    TopoDS_Shape box = BRepPrimAPI_MakeBox(corner, wire_w, tubeLength, wire_h).Shape();
                    box = translate_shape(box, innerX, (tubeLength / 2.0) * ySign, 0.0);
                    if (!box.IsNull()) builder.Add(compound, box);
                }
                // 2. Inner corner: sweep starts at +X (tube end), turns toward the chord
                {
                    gp_Pnt center(innerX - bendRadius, tubeLength * ySign, 0.0);
                    TopoDS_Shape s = make_toroidal_quarter_swept_rectangle(
                        bendRadius, wire_w, wire_h, center, 0.0, ySign);
                    s = rotate_about_axis(s, gp_Pnt(innerX, 0, 0), gp_Dir(0,1,0), psiIn);
                    if (!s.IsNull()) builder.Add(compound, s);
                }
                // 3. Chord run between the two corner faces (top view), wire_w vertical
                if (chordLen > 1e-9) {
                    gp_Pnt p0(chordInX, radialHeight * ySign, chordInZ);
                    gp_Dir chordDir((chordOutX - chordInX) / chordLen, 0.0,
                                    (chordOutZ - chordInZ) / chordLen);
                    gp_Ax2 profPlane(p0, chordDir, gp_Dir(0, 1, 0));
                    TopoDS_Face face = build_rect_profile(profPlane, wire_w, wire_h);
                    TopoDS_Shape box = BRepPrimAPI_MakePrism(
                        face, gp_Vec(chordDir.XYZ() * chordLen)).Shape();
                    if (!box.IsNull()) builder.Add(compound, box);
                }
                // 4. Outer corner: sweep starts at radial connection (±90°), turns toward outer tube;
                //    placed by the GLOBAL-axis rotation, then yawed about its own tube axis
                {
                    double startAngleRad = (ySign > 0) ? (std::numbers::pi / 2.0) : (3.0 * std::numbers::pi / 2.0);
                    gp_Pnt center(outerX + bendRadius, tubeLength * ySign, 0.0);
                    TopoDS_Shape s = make_toroidal_quarter_swept_rectangle(
                        bendRadius, wire_w, wire_h, center, startAngleRad, ySign);
                    s = rotate_about_axis(s, gp_Pnt(0, 0, 0), gp_Dir(0,1,0), deltaOut);
                    s = rotate_about_axis(s, gp_Pnt(tOutX, 0, tOutZ), gp_Dir(0,1,0), psiOut);
                    if (!s.IsNull()) builder.Add(compound, s);
                }
                // 5. Outer tube: exactly on additionalCoordinates[0] via the global-axis rotation
                {
                    gp_Pnt corner(-wire_w / 2.0, -tubeLength / 2.0, -wire_h / 2.0);
                    TopoDS_Shape box = BRepPrimAPI_MakeBox(corner, wire_w, tubeLength, wire_h).Shape();
                    box = translate_shape(box, outerX, (tubeLength / 2.0) * ySign, 0.0);
                    box = rotate_about_axis(box, gp_Pnt(0, 0, 0), gp_Dir(0,1,0), deltaOut);
                    if (!box.IsNull()) builder.Add(compound, box);
                }
            };

            build_rect_half(+1.0);
            build_rect_half(-1.0);
        } else {
            auto build_half = [&](double ySign) {
                // 1. Inner tube along Y from origin
                {
                    TopoDS_Shape c = BRepPrimAPI_MakeCylinder(wire_radius, tubeLength).Shape();
                    c = rotate_about_axis(c, gp_Pnt(0,0,0), gp_Dir(1,0,0), -std::numbers::pi/2.0 * ySign);
                    c = translate_shape(c, innerX, 0.0, 0.0);
                    if (!c.IsNull()) builder.Add(compound, c);
                }
                // 2. Inner corner: turns from the tube toward the chord
                {
                    double sa = (ySign > 0) ? 0.0 : 270.0;
                    gp_Pnt center(innerX - bendRadius, tubeLength * ySign, 0.0);
                    TopoDS_Shape s = make_quarter_torus_at(bendRadius, wire_radius, center, sa);
                    s = rotate_about_axis(s, gp_Pnt(innerX, 0, 0), gp_Dir(0,1,0), psiIn);
                    if (!s.IsNull()) builder.Add(compound, s);
                }
                // 3. Chord run between the two corner faces (top view)
                if (chordLen > 1e-9) {
                    gp_Pnt p0(chordInX, radialHeight * ySign, chordInZ);
                    gp_Dir chordDir((chordOutX - chordInX) / chordLen, 0.0,
                                    (chordOutZ - chordInZ) / chordLen);
                    TopoDS_Shape c = BRepPrimAPI_MakeCylinder(
                        gp_Ax2(p0, chordDir), wire_radius, chordLen).Shape();
                    if (!c.IsNull()) builder.Add(compound, c);
                }
                // 4. Outer corner: placed by the GLOBAL-axis rotation, then yawed to meet the chord
                {
                    double sa = (ySign > 0) ? 90.0 : 180.0;
                    gp_Pnt center(outerX + bendRadius, tubeLength * ySign, 0.0);
                    TopoDS_Shape s = make_quarter_torus_at(bendRadius, wire_radius, center, sa);
                    s = rotate_about_axis(s, gp_Pnt(0, 0, 0), gp_Dir(0,1,0), deltaOut);
                    s = rotate_about_axis(s, gp_Pnt(tOutX, 0, tOutZ), gp_Dir(0,1,0), psiOut);
                    if (!s.IsNull()) builder.Add(compound, s);
                }
                // 5. Outer tube: exactly on additionalCoordinates[0] via the global-axis rotation
                {
                    TopoDS_Shape c = BRepPrimAPI_MakeCylinder(wire_radius, tubeLength).Shape();
                    c = rotate_about_axis(c, gp_Pnt(0,0,0), gp_Dir(1,0,0), -std::numbers::pi/2.0 * ySign);
                    c = translate_shape(c, outerX, 0.0, 0.0);
                    c = rotate_about_axis(c, gp_Pnt(0, 0, 0), gp_Dir(0,1,0), deltaOut);
                    if (!c.IsNull()) builder.Add(compound, c);
                }
            };

            build_half(+1.0);
            build_half(-1.0);
        }

        canonical = compound;
        s_toroidalTurnCache[key] = canonical;
    }

    // Apply the per-turn Y-axis rotation (cheap BRepBuilderAPI_Transform on cached compound)
    return rotate_about_axis(canonical, gp_Pnt(0,0,0), gp_Dir(0,1,0), turnRotationRad);
}

static TopoDS_Shape build_concentric_oblong_turn(double radial_pos, double wire_radius,
                                                  double height_pos,
                                                  double half_col_width, double half_col_depth,
                                                  bool rect_wire,
                                                  double wire_width, double wire_height,
                                                  int segments = DEFAULT_WIRE_POLYGON_SEGMENTS) {
    double straight_half = half_col_depth - half_col_width;
    if (straight_half <= 0.0) {
        // Effectively round column
        return build_concentric_round_column_turn(radial_pos, wire_radius, height_pos,
                                                   rect_wire, wire_width, wire_height, segments);
    }

    double tube_z_length = 2.0 * straight_half;
    double prof_r = rect_wire ? std::min(wire_width, wire_height) / 2.0 : wire_radius;
    // Rectangular/planar wire cross-section: wire_width is the radial extent (X),
    // wire_height the axial/stacking extent (Y). Round wire uses prof_r.
    double hw = wire_width / 2.0;   // radial half-width (X)
    double hh = wire_height / 2.0;  // vertical half-height (Y)

    // Straight segments run along Z at x = ±radial_pos, centred on z = 0 so their
    // ends meet the bends at z = ±straight_half. (Previously based at z = 0 and
    // extended +Z for the whole length, leaving the straight part shifted by
    // +straight_half — i.e. not centred in Z.) Rectangular/planar wire → a box of
    // the true foil cross-section; round wire → a cylinder.
    auto make_straight = [&](double cx) -> TopoDS_Shape {
        if (rect_wire) {
            gp_Pnt corner(cx - hw, height_pos - hh, -straight_half);
            return BRepPrimAPI_MakeBox(corner, wire_width, wire_height, tube_z_length).Shape();
        }
        gp_Ax2 ax(gp_Pnt(cx, height_pos, -straight_half), gp_Dir(0, 0, 1));
        return BRepPrimAPI_MakeCylinder(ax, prof_r, tube_z_length).Shape();
    };

    TopoDS_Shape tube_px = make_straight(radial_pos);
    TopoDS_Shape tube_nx = make_straight(-radial_pos);

    // Bends at the ±Z ends. The turn lies flat in the X-Z plane at height y, so the
    // bend must sweep in X-Z: revolve the wire cross-section 180° about the Y axis
    // through the end-centre. (Previously this revolved about the Z axis, which
    // swept the bend up in the X-Y plane — the "fountain" artefact.)
    // The cross-section matches the straight segment it continues: a rectangle
    // (radial wire_width × vertical wire_height) in the X-Y plane for rect/planar
    // wire — NOT a round tube — or a circle for round wire. The profile lies in the
    // X-Y plane (z = cz), coplanar with the Y revolve axis.
    // y_axis_sign selects the outward bulge: +Z end bulges +Z, -Z end bulges -Z.
    auto make_bend = [&](double cz, double y_axis_sign) -> TopoDS_Shape {
        TopoDS_Face profile;
        if (rect_wire) {
            BRepBuilderAPI_MakePolygon poly;
            poly.Add(gp_Pnt(radial_pos - hw, height_pos - hh, cz));
            poly.Add(gp_Pnt(radial_pos + hw, height_pos - hh, cz));
            poly.Add(gp_Pnt(radial_pos + hw, height_pos + hh, cz));
            poly.Add(gp_Pnt(radial_pos - hw, height_pos + hh, cz));
            poly.Close();
            BRepBuilderAPI_MakeFace fm(poly.Wire());
            if (!fm.IsDone()) {
                std::cerr << "ERROR build_concentric_oblong_turn: MakeFace failed for rect bend\n";
                return TopoDS_Shape();
            }
            profile = fm.Face();
        } else {
            gp_Ax2 circ_axis(gp_Pnt(radial_pos, height_pos, cz), gp_Dir(0, 0, 1));
            Handle(Geom_Circle) circle = GC_MakeCircle(circ_axis, prof_r).Value();
            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circle).Edge();
            TopoDS_Wire circ_wire = BRepBuilderAPI_MakeWire(edge).Wire();
            profile = BRepBuilderAPI_MakeFace(circ_wire).Face();
        }

        gp_Ax1 rev_axis(gp_Pnt(0, height_pos, cz), gp_Dir(0, y_axis_sign, 0));
        BRepPrimAPI_MakeRevol revol(profile, rev_axis, std::numbers::pi);
        try {
            return revol.Shape();
        } catch (const Standard_Failure& e) {
            std::cerr << "ERROR build_concentric_oblong_turn: MakeRevol threw Standard_Failure: " << e.GetMessageString() << "\n";
            return TopoDS_Shape();
        } catch (...) {
            std::cerr << "ERROR build_concentric_oblong_turn: MakeRevol threw unknown exception\n";
            return TopoDS_Shape();
        }
    };

    TopoDS_Shape torus_pz = make_bend(straight_half, -1.0);
    TopoDS_Shape torus_nz = make_bend(-straight_half, 1.0);

    // Return compound directly — fusing unnecessary for STL visualization.
    BRep_Builder cb;
    TopoDS_Compound comp;
    cb.MakeCompound(comp);
    if (!tube_px.IsNull()) cb.Add(comp, tube_px);
    if (!tube_nx.IsNull()) cb.Add(comp, tube_nx);
    if (!torus_pz.IsNull()) cb.Add(comp, torus_pz);
    if (!torus_nz.IsNull()) cb.Add(comp, torus_nz);
    return comp;
}

std::pair<double, double> TurnBuilder::wireDimensions(const MAS::Wire& wire, bool paintCoating) {
    MAS::Turn contextlessTurn;   // no per-turn dimensions: the wire's own data decides
    return get_wire_dimensions(wire, contextlessTurn, paintCoating);
}

std::pair<double, double> TurnBuilder::wireDimensions(const MAS::Wire& wire,
                                                      const MAS::Turn& turn,
                                                      bool paintCoating) {
    return get_wire_dimensions(wire, turn, paintCoating);
}

TopoDS_Shape TurnBuilder::fuseSolids(const TopoDS_Shape& compound) {
    return fuse_turn_solids(compound);
}

TopoDS_Shape TurnBuilder::buildTurn(const MAS::Turn& turn,
                                    const MAS::Wire& wire,
                                    const MAS::CoreBobbinProcessedDescription& bobbin,
                                    bool isToroidal,
                                    int wirePolygonSegments,
                                    int wireRevolutionSegments,
                                    bool paintCoating,
                                    const std::optional<WoundColumnSpec>& woundColumn) {
    const auto& coords = turn.get_coordinates();
    if (coords.size() < 2) {
        throw std::runtime_error("TurnBuilder::buildTurn: turn '"
                                 + (turn.get_name().empty() ? std::string("<unnamed>") : turn.get_name())
                                 + "' has fewer than 2 coordinates");
    }
    double radial_pos = coords[0];
    double height_pos = coords[1];

    bool rect_wire = is_rectangular_wire(wire);
    auto [wire_w, wire_h] = get_wire_dimensions(wire, turn, paintCoating);
    double wire_radius = std::min(wire_w, wire_h) / 2.0;

    // Dispatch by topology/column shape; maybe_fuse_turn() collapses the per-turn sub-solid
    // compound into ONE solid when fusing is enabled (FEM meshing), else returns it unchanged
    // (fast compound for visualisation). Applied here so every column shape + toroidal is covered.
    if (isToroidal) {
        if (woundColumn) {
            throw std::runtime_error(
                "TurnBuilder::buildTurn: turn '" + turn.get_name()
                + "' carries a wound-column spec on a toroidal core — toroids have a single "
                  "column, so multi-column placement data here means inconsistent MAS");
        }
        return maybe_fuse_turn(build_toroidal_turn(turn, wire, bobbin, paintCoating));
    }

    // Multi-column placement: the turn wraps a NON-MAIN column. Build the same turn shapes
    // relative to that column's axis (the radial offset |x - axisX| preserves the clearance
    // MKF computed — turn coordinates are absolute MAS data and must not be adjusted here),
    // then translate the loop to the leg position. additionalCoordinates (the second window
    // crossing MKF emits for 2D consumers) are intentionally unused: the closed 3D loop
    // already passes through both crossings.
    if (woundColumn) {
        const double rel_radial = std::abs(radial_pos - woundColumn->axisX);
        if (rel_radial <= woundColumn->halfWidth &&
            woundColumn->shape == MAS::ColumnShape::ROUND) {
            // The rect/oblong builders validate this themselves; ROUND takes no column dims.
            throw std::runtime_error(
                "TurnBuilder::buildTurn: turn '" + turn.get_name() + "' radial offset "
                + std::to_string(rel_radial) + " m lies inside its wound column (half-width "
                + std::to_string(woundColumn->halfWidth) + " m at x = "
                + std::to_string(woundColumn->axisX) + " m) — inconsistent MAS turn/column data");
        }
        TopoDS_Shape loop;
        switch (woundColumn->shape) {
            case MAS::ColumnShape::ROUND:
                loop = build_concentric_round_column_turn(rel_radial, wire_radius, height_pos,
                                                          rect_wire, wire_w, wire_h,
                                                          wirePolygonSegments, wireRevolutionSegments);
                break;
            case MAS::ColumnShape::OBLONG:
                loop = build_concentric_oblong_turn(rel_radial, wire_radius, height_pos,
                                                    woundColumn->halfWidth, woundColumn->halfDepth,
                                                    rect_wire, wire_w, wire_h, wirePolygonSegments);
                break;
            default:  // RECTANGULAR / IRREGULAR — same default the main-column dispatch uses
                loop = build_concentric_rect_column_turn(rel_radial, wire_radius, height_pos,
                                                         woundColumn->halfWidth, woundColumn->halfDepth,
                                                         rect_wire, wire_w, wire_h);
                break;
        }
        if (!loop.IsNull() && std::abs(woundColumn->axisX) > 1e-12) {
            gp_Trsf shift;
            shift.SetTranslation(gp_Vec(woundColumn->axisX, 0.0, 0.0));
            loop = BRepBuilderAPI_Transform(loop, shift).Shape();
        }
        return maybe_fuse_turn(loop);
    }

    double half_col_width = bobbin.get_column_width().value_or(0.0);
    double half_col_depth = bobbin.get_column_depth();

    if (bobbin.get_column_shape() == MAS::ColumnShape::ROUND) {
        return maybe_fuse_turn(build_concentric_round_column_turn(radial_pos, wire_radius, height_pos,
                                                   rect_wire, wire_w, wire_h,
                                                   wirePolygonSegments,
                                                   wireRevolutionSegments));
    }

    if (bobbin.get_column_shape() == MAS::ColumnShape::OBLONG) {
        return maybe_fuse_turn(build_concentric_oblong_turn(radial_pos, wire_radius, height_pos,
                                             half_col_width, half_col_depth,
                                             rect_wire, wire_w, wire_h, wirePolygonSegments));
    }

    // Default: rectangular column
    return maybe_fuse_turn(build_concentric_rect_column_turn(radial_pos, wire_radius, height_pos,
                                              half_col_width, half_col_depth,
                                              rect_wire, wire_w, wire_h));
}

void TurnBuilder::clearCache() {
    s_toroidalTurnCache.clear();
}

TopoDS_Shape TurnBuilder::buildFromTurnAlone(const MAS::Turn& turn,
                                              int wirePolygonSegments,
                                              bool paintCoating) {
    if (!paintCoating) {
        // A bare Turn carries only its OUTER footprint (turn.dimensions); the
        // copper / strand geometry needed to draw the conducting cross-section
        // is not present. Surface this loudly instead of silently emitting the
        // outer solid mislabelled as copper.
        throw std::runtime_error(
            "TurnBuilder::buildFromTurnAlone: conducting-diameter geometry "
            "(paintCoating=false) requires wire/strand data that a standalone "
            "Turn does not carry. Use drawMagnetic with a full Magnetic JSON "
            "(coil + core + wire) for copper-cross-section turns.");
    }
    // Validate required Turn fields explicitly — no silent fallbacks.
    const auto& coords = turn.get_coordinates();
    if (coords.size() < 2) {
        throw std::runtime_error(
            "TurnBuilder::buildFromTurnAlone: Turn '" + turn.get_name()
            + "' must have at least 2 coordinates");
    }
    auto dimsOpt = turn.get_dimensions();
    if (!dimsOpt || dimsOpt->size() < 2) {
        throw std::runtime_error(
            "TurnBuilder::buildFromTurnAlone: Turn '" + turn.get_name()
            + "' must carry dimensions [width, height] in standalone mode");
    }
    auto crossOpt = turn.get_cross_sectional_shape();
    if (!crossOpt) {
        throw std::runtime_error(
            "TurnBuilder::buildFromTurnAlone: Turn '" + turn.get_name()
            + "' must carry crossSectionalShape (ROUND/RECTANGULAR/OVAL) in standalone mode");
    }

    // Synthesize a minimal Wire matching the cross-section, so we can reuse
    // the well-tested buildTurn() pipeline. Dimensions on the Turn already
    // take precedence inside get_wire_dimensions().
    MAS::Wire wire;
    if (*crossOpt == MAS::TurnCrossSectionalShape::RECTANGULAR ||
        *crossOpt == MAS::TurnCrossSectionalShape::OVAL) {
        wire.set_type(MAS::WireType::RECTANGULAR);
        MAS::DimensionWithTolerance w; w.set_nominal((*dimsOpt)[0]);
        MAS::DimensionWithTolerance h; h.set_nominal((*dimsOpt)[1]);
        wire.set_outer_width(std::optional<MAS::DimensionWithTolerance>(w));
        wire.set_outer_height(std::optional<MAS::DimensionWithTolerance>(h));
        wire.set_conducting_width(std::optional<MAS::DimensionWithTolerance>(w));
        wire.set_conducting_height(std::optional<MAS::DimensionWithTolerance>(h));
    } else {
        wire.set_type(MAS::WireType::ROUND);
        MAS::DimensionWithTolerance d; d.set_nominal((*dimsOpt)[0]);
        wire.set_outer_diameter(std::optional<MAS::DimensionWithTolerance>(d));
        wire.set_conducting_diameter(std::optional<MAS::DimensionWithTolerance>(d));
    }

    bool hasSecondCrossing = turn.get_additional_coordinates().has_value()
                              && !turn.get_additional_coordinates()->empty();

    if (hasSecondCrossing) {
        // additionalCoordinates mark either a toroidal turn's outer-ring crossing or a
        // concentric turn's second window crossing (multi-column placement). Neither can be
        // resolved without the surrounding bobbin/core context, so refuse loudly instead of
        // guessing a topology.
        throw std::runtime_error(
            "TurnBuilder::buildFromTurnAlone: turn '" + turn.get_name()
            + "' carries additional_coordinates (toroidal outer crossing or multi-column "
              "second window crossing) which needs full bobbin/core context. Use drawWinding "
              "or drawMagnetic with a full Magnetic JSON instead of standalone drawTurns.");
    }

    // Concentric round-column turn: bobbin's column_shape is the only field
    // buildTurn consults for ROUND, and turn.coordinates[0] becomes the
    // turn radius — exactly what we want.
    MAS::CoreBobbinProcessedDescription bobbin;
    bobbin.set_column_shape(MAS::ColumnShape::ROUND);
    bobbin.set_column_depth(0.0);
    bobbin.set_column_thickness(0.0);
    bobbin.set_wall_thickness(0.0);
    bobbin.set_winding_windows(std::vector<MAS::WindingWindowElement>{});

    return buildTurn(turn, wire, bobbin, /*isToroidal=*/false, wirePolygonSegments);
}

} // namespace mvb
