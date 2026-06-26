// 3D regression test for concentric turns on an OBLONG (racetrack) column,
// as used by planar/rectangular foil windings (e.g. planar EL cores).
//
// Guards three bugs fixed in build_concentric_oblong_turn:
//   1. Bends were revolved about the Z axis, sweeping the arc up in the X-Y
//      plane (a vertical "fountain") instead of in the X-Z plane (the flat
//      plane of the turn). A correct turn is THIN in Y (≈ wire height), not
//      tall by ~the turn radius.
//   2. Straight segments were based at z = 0 and extended +Z for the whole
//      length, leaving them shifted by +straight_half (not centred in Z).
//   3. Rectangular/planar wire was approximated as a ROUND tube of radius
//      min(width,height)/2 instead of a true rectangular foil cross-section.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mvb/MagneticBuilder.h"
#include "MAS.hpp"

#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Shape.hxx>

#include <vector>
#include <cmath>

using Catch::Matchers::WithinAbs;

namespace {

// Oblong column: half-width (X) and half-depth (Z); the straight legs run
// along Z for 2*(half_depth - half_width). These mirror a planar EL bobbin.
constexpr double HALF_COL_WIDTH = 0.00196;    // m
constexpr double HALF_COL_DEPTH = 0.004505;   // m
constexpr double STRAIGHT_HALF  = HALF_COL_DEPTH - HALF_COL_WIDTH;  // 0.002545 m

constexpr double RADIAL_POS = 0.00258;   // turn radial position (m)
constexpr double HEIGHT_POS = 0.0004;    // turn height in Y (m)

MAS::DimensionWithTolerance nominal_dim(double v) {
    MAS::DimensionWithTolerance d;
    d.set_nominal(v);
    return d;
}

MAS::MagneticCore make_oblong_core() {
    MAS::CoreShape shape;
    shape.set_family(MAS::CoreShapeFamily::PLANAR_EL);
    shape.set_type(MAS::FunctionalDescriptionType::STANDARD);
    MAS::CoreFunctionalDescription fd;
    fd.set_material(MAS::CoreMaterialDataOrNameUnion{std::string("3F4")});
    fd.set_type(MAS::CoreType::TWO_PIECE_SET);   // not toroidal
    fd.set_number_stacks(int64_t(1));
    fd.set_shape(shape);
    MAS::MagneticCore core;
    core.set_functional_description(fd);
    return core;
}

MAS::Wire make_planar_wire(double width, double height, const std::string& name) {
    MAS::Wire wire;
    wire.set_type(MAS::WireType::PLANAR);
    wire.set_conducting_width(nominal_dim(width));
    wire.set_conducting_height(nominal_dim(height));
    wire.set_outer_width(nominal_dim(width));
    wire.set_outer_height(nominal_dim(height));
    wire.set_name(std::optional<std::string>(name));
    return wire;
}

MAS::Wire make_round_wire(double diameter, const std::string& name) {
    MAS::Wire wire;
    wire.set_type(MAS::WireType::ROUND);
    wire.set_conducting_diameter(nominal_dim(diameter));
    wire.set_outer_diameter(nominal_dim(diameter));
    wire.set_name(std::optional<std::string>(name));
    return wire;
}

MAS::Magnetic make_magnetic(const MAS::Wire& wire) {
    MAS::CoreBobbinProcessedDescription bobbinPd;
    bobbinPd.set_column_width(HALF_COL_WIDTH);
    bobbinPd.set_column_depth(HALF_COL_DEPTH);
    bobbinPd.set_column_thickness(0.0);
    bobbinPd.set_wall_thickness(0.0);
    bobbinPd.set_column_shape(MAS::ColumnShape::OBLONG);

    MAS::Bobbin bobbin;
    bobbin.set_processed_description(
        std::optional<MAS::CoreBobbinProcessedDescription>(bobbinPd));

    MAS::Turn turn;
    turn.set_name("Primary parallel 0 turn 0");
    turn.set_winding("Primary");
    turn.set_coordinates({RADIAL_POS, HEIGHT_POS});
    turn.set_rotation(std::optional<double>(0.0));

    MAS::CoilFunctionalDescription coilFunc;
    coilFunc.set_name("Primary");
    coilFunc.set_number_turns(1);
    coilFunc.set_number_parallels(1);
    coilFunc.set_isolation_side(MAS::IsolationSide::PRIMARY);
    coilFunc.set_wire(MAS::WireDataOrNameUnion(wire));

    MAS::Coil coil;
    coil.set_bobbin(MAS::BobbinDataOrNameUnion(bobbin));
    coil.set_turns_description(std::optional<std::vector<MAS::Turn>>(std::vector<MAS::Turn>{turn}));
    coil.set_functional_description(std::vector<MAS::CoilFunctionalDescription>{coilFunc});

    MAS::Magnetic magnetic;
    magnetic.set_core(make_oblong_core());
    magnetic.set_coil(coil);
    return magnetic;
}

}  // namespace

TEST_CASE("Oblong planar turn: flat racetrack in X-Z, not a vertical fountain",
          "[3d][oblong][planar][turn]") {
    const double wireWidth  = 0.000562;
    const double wireHeight = 0.00010439;
    auto magnetic = make_magnetic(make_planar_wire(wireWidth, wireHeight, "Planar 104um"));

    mvb::MagneticBuilder builder;
    auto turns = builder.buildTurnsNamed(magnetic.get_coil().value(), magnetic.get_core().value(),
                                         /*wirePolygonSegments=*/0);
    REQUIRE(turns.size() == 1);
    REQUIRE(!turns[0].shape.IsNull());

    GProp_GProps props;
    BRepGProp::VolumeProperties(turns[0].shape, props);
    REQUIRE(props.Mass() > 0.0);

    Bnd_Box bb;
    BRepBndLib::Add(turns[0].shape, bb);
    double xmin, ymin, zmin, xmax, ymax, zmax;
    bb.Get(xmin, ymin, zmin, xmax, ymax, zmax);

    const double dX = xmax - xmin;
    const double dY = ymax - ymin;
    const double dZ = zmax - zmin;
    INFO("turn bbox (mm): dX=" << dX*1000 << " dY=" << dY*1000 << " dZ=" << dZ*1000);

    // (1) FLAT: vertical extent is just the foil thickness, NOT ~the turn
    // diameter. The fountain bug made dY ≈ 2*RADIAL_POS (~5 mm).
    REQUIRE_THAT(dY, WithinAbs(wireHeight, 5e-6));
    REQUIRE(dY < RADIAL_POS);            // emphatically not a fountain

    // The loop lies in X-Z and is longer along Z (the oblong long axis) than X.
    REQUIRE(dZ > dX);
    REQUIRE(dX > 2.0 * RADIAL_POS);      // spans both ±X legs

    // (2) CENTRED IN Z: the racetrack is symmetric about z = 0 (the straights
    // span [-straight_half, +straight_half], bends bulge equally each end).
    REQUIRE_THAT(zmin + zmax, WithinAbs(0.0, 1e-4));

    // Vertically centred on the turn's height coordinate.
    REQUIRE_THAT(0.5 * (ymin + ymax), WithinAbs(HEIGHT_POS, 1e-5));
}

TEST_CASE("Oblong planar turn: straight legs are rectangular foil, centred in Z",
          "[3d][oblong][planar][turn]") {
    const double wireWidth  = 0.000562;
    const double wireHeight = 0.00010439;
    auto magnetic = make_magnetic(make_planar_wire(wireWidth, wireHeight, "Planar 104um"));

    mvb::MagneticBuilder builder;
    auto turns = builder.buildTurnsNamed(magnetic.get_coil().value(), magnetic.get_core().value(),
                                         /*wirePolygonSegments=*/0);
    REQUIRE(turns.size() == 1);

    // The straight legs are the sub-solids that cross z = 0 (the bends sit at
    // one Z end each). Each straight must be a true rectangular foil:
    //   radial width  = wire_width   (NOT min(w,h) → would be the round tube)
    //   stack height  = wire_height
    //   z-span        = [-straight_half, +straight_half]   (centred)
    int straightCount = 0;
    for (TopExp_Explorer exp(turns[0].shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        Bnd_Box bb;
        BRepBndLib::Add(exp.Current(), bb);
        double x0, y0, z0, x1, y1, z1;
        bb.Get(x0, y0, z0, x1, y1, z1);
        if (z0 < 0.0 && z1 > 0.0) {       // crosses z = 0 → a straight leg
            ++straightCount;
            INFO("straight leg bbox (mm): x[" << x0*1000 << "," << x1*1000
                 << "] y[" << y0*1000 << "," << y1*1000
                 << "] z[" << z0*1000 << "," << z1*1000 << "]");
            REQUIRE_THAT(x1 - x0, WithinAbs(wireWidth, 2e-5));   // rectangular foil width
            REQUIRE_THAT(y1 - y0, WithinAbs(wireHeight, 2e-5));  // foil thickness
            REQUIRE_THAT(z0, WithinAbs(-STRAIGHT_HALF, 1e-5));   // centred in Z
            REQUIRE_THAT(z1, WithinAbs(+STRAIGHT_HALF, 1e-5));
        }
    }
    REQUIRE(straightCount == 2);   // +X and -X legs
}

TEST_CASE("Oblong round wire turn: round cross-section preserved",
          "[3d][oblong][round][turn]") {
    const double diameter = 0.0004;
    auto magnetic = make_magnetic(make_round_wire(diameter, "Round 0.40"));

    mvb::MagneticBuilder builder;
    auto turns = builder.buildTurnsNamed(magnetic.get_coil().value(), magnetic.get_core().value(),
                                         /*wirePolygonSegments=*/0);
    REQUIRE(turns.size() == 1);
    REQUIRE(!turns[0].shape.IsNull());

    // Round wire still flows through the cylinder/torus path: a straight leg's
    // cross-section is a circle (radial width ≈ stack height ≈ diameter).
    for (TopExp_Explorer exp(turns[0].shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        Bnd_Box bb;
        BRepBndLib::Add(exp.Current(), bb);
        double x0, y0, z0, x1, y1, z1;
        bb.Get(x0, y0, z0, x1, y1, z1);
        if (z0 < 0.0 && z1 > 0.0) {       // a straight leg
            REQUIRE_THAT(x1 - x0, WithinAbs(diameter, 5e-5));
            REQUIRE_THAT(y1 - y0, WithinAbs(diameter, 5e-5));
            REQUIRE_THAT(x1 - x0, WithinAbs(y1 - y0, 1e-5));  // round → square bbox
        }
    }
}
