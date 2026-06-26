#include <catch2/catch_test_macros.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/BobbinBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <numbers>
#include <vector>
#include <map>
#include <filesystem>

static MAS::Dimension make_dim(double v) { return MAS::Dimension(v); }

TEST_CASE("E core with bobbin and turns builds complete assembly", "[integration]") {
    MAS::Magnetic magnetic;

    // Core: E 19/8/5
    MAS::CoreShape shape;
    shape.set_family(MAS::CoreShapeFamily::E);
    shape.set_type(MAS::FunctionalDescriptionType::STANDARD);
    std::map<std::string, MAS::Dimension> dims;
    dims["A"] = make_dim(0.019);
    dims["B"] = make_dim(0.008);
    dims["C"] = make_dim(0.005);
    dims["D"] = make_dim(0.005);
    dims["E"] = make_dim(0.012);
    dims["F"] = make_dim(0.006);
    shape.set_dimensions(dims);

    MAS::CoreGeometricalDescriptionElement piece1;
    piece1.set_type(MAS::CoreGeometricalDescriptionElementType::HALF_SET);
    piece1.set_coordinates({0.0, 0.0, 0.0});
    piece1.set_rotation(std::optional<std::vector<double>>(std::vector<double>{0.0, 0.0, 0.0}));
    piece1.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    MAS::CoreGeometricalDescriptionElement piece2;
    piece2.set_type(MAS::CoreGeometricalDescriptionElementType::HALF_SET);
    piece2.set_coordinates({0.0, 0.0, 0.0});
    piece2.set_rotation(std::optional<std::vector<double>>(std::vector<double>{std::numbers::pi, 0.0, 0.0}));
    piece2.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    MAS::MagneticCore core;
    core.set_geometrical_description(std::optional<std::vector<MAS::CoreGeometricalDescriptionElement>>(std::vector<MAS::CoreGeometricalDescriptionElement>{piece1, piece2}));
    magnetic.set_core(core);

    // Bobbin processed description
    MAS::CoreBobbinProcessedDescription bobbinPd;
    bobbinPd.set_column_width(0.007);
    bobbinPd.set_column_depth(0.004);
    bobbinPd.set_column_thickness(0.006);
    bobbinPd.set_wall_thickness(0.0005);
    bobbinPd.set_column_shape(MAS::ColumnShape::RECTANGULAR);

    MAS::Bobbin bobbin;
    bobbin.set_processed_description(std::optional<MAS::CoreBobbinProcessedDescription>(bobbinPd));

    // Turns
    std::vector<MAS::Turn> turns;
    for (int i = 0; i < 3; ++i) {
        MAS::Turn turn;
        turn.set_name("Turn_" + std::to_string(i));
        turn.set_winding("Primary");
        turn.set_length(0.05);
        turn.set_parallel(1);
        // radial position, height position
        double radial = 0.005;
        double height = -0.002 + i * 0.0012;
        turn.set_coordinates({radial, height});
        turns.push_back(turn);
    }

    MAS::CoilFunctionalDescription coilFunc;
    coilFunc.set_name("Primary");
    coilFunc.set_number_turns(3);
    coilFunc.set_number_parallels(1);
    coilFunc.set_isolation_side(MAS::IsolationSide::PRIMARY);
    coilFunc.set_wire("Round 1.00 - Grade 1");

    MAS::Coil coil;
    coil.set_bobbin(MAS::BobbinDataOrNameUnion(bobbin));
    coil.set_turns_description(std::optional<std::vector<MAS::Turn>>(turns));
    coil.set_functional_description(std::vector<MAS::CoilFunctionalDescription>{coilFunc});
    magnetic.set_coil(coil);

    mvb::MagneticBuilder builder;
    auto coreNamed = builder.buildCoreNamed(magnetic.get_core().value());
    auto bobbinNamed = builder.buildBobbinNamed(magnetic.get_coil().value(), magnetic.get_core().value());
    auto turnsNamed = builder.buildTurnsNamed(magnetic.get_coil().value(), magnetic.get_core().value());

    REQUIRE(coreNamed.size() == 2);
    REQUIRE(!bobbinNamed.shape.IsNull());
    REQUIRE(turnsNamed.size() == 3);

    // Combine and sanity check
    std::vector<TopoDS_Shape> all;
    for (const auto& ns : coreNamed) all.push_back(ns.shape);
    all.push_back(bobbinNamed.shape);
    for (const auto& ns : turnsNamed) all.push_back(ns.shape);

    double totalVolume = 0.0;
    int solids = 0;
    for (const auto& s : all) {
        REQUIRE(!s.IsNull());
        GProp_GProps props;
        BRepGProp::VolumeProperties(s, props);
        totalVolume += props.Mass();
        for (TopExp_Explorer exp(s, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    }

    REQUIRE(totalVolume > 0.0);
    REQUIRE(solids >= 6); // 2 core + 3 bobbin (body+flanges as compound) + 3 turns
}

// Verifies the polygonSegments knob on BobbinBuilder controls bobbin
// tessellation. With round-column bobbins, segment count drives the number
// of planar lateral faces on the body cylinder; with `0` we fall back to
// NURBS (1 cylindrical face). This mirrors the per-core polygonSegments
// behaviour and is what drawMagnetic plumbs through `corePolygonSegments`.
TEST_CASE("BobbinBuilder polygonSegments controls round bobbin tessellation", "[bobbin][polygonSegments]") {
    MAS::CoreBobbinProcessedDescription pd;
    pd.set_column_width(0.005);
    pd.set_column_depth(0.005);
    pd.set_column_thickness(0.001);
    pd.set_wall_thickness(0.0005);
    pd.set_column_shape(MAS::ColumnShape::ROUND);
    std::vector<MAS::WindingWindowElement> wwList;
    MAS::WindingWindowElement ww;
    ww.set_width(0.004);
    ww.set_height(0.010);
    wwList.push_back(ww);
    pd.set_winding_windows(wwList);

    auto countFaces = [](const TopoDS_Shape& s) {
        int n = 0;
        for (TopExp_Explorer exp(s, TopAbs_FACE); exp.More(); exp.Next()) ++n;
        return n;
    };
    auto countSolids = [](const TopoDS_Shape& s) {
        int n = 0;
        for (TopExp_Explorer exp(s, TopAbs_SOLID); exp.More(); exp.Next()) ++n;
        return n;
    };

    TopoDS_Shape b8  = mvb::BobbinBuilder::buildBobbin(pd, 0.0005, /*axisIsY=*/false, /*polygonSegments=*/8);
    TopoDS_Shape b32 = mvb::BobbinBuilder::buildBobbin(pd, 0.0005, /*axisIsY=*/false, /*polygonSegments=*/32);
    TopoDS_Shape b0  = mvb::BobbinBuilder::buildBobbin(pd, 0.0005, /*axisIsY=*/false, /*polygonSegments=*/0);

    REQUIRE(!b8.IsNull());
    REQUIRE(!b32.IsNull());
    REQUIRE(!b0.IsNull());

    // Bobbin with flanges should now be a compound of 3 solids (body + 2 flanges).
    REQUIRE(countSolids(b8) == 3);
    REQUIRE(countSolids(b32) == 3);
    REQUIRE(countSolids(b0) == 3);

    int f8  = countFaces(b8);
    int f32 = countFaces(b32);
    int f0  = countFaces(b0);

    // More tessellation segments ⇒ strictly more faces.
    REQUIRE(f32 > f8);
    // NURBS path (segments=0) collapses each cylindrical wall to a single
    // surface, so it has the smallest face count.
    REQUIRE(f0 < f8);
}
