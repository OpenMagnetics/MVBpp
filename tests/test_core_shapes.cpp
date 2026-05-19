#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <vector>
#include <numbers>

using Catch::Matchers::WithinRel;

static MAS::Dimension make_dim(double v) {
    return MAS::Dimension(v);
}

static MAS::Magnetic make_simple_e_magnetic() {
    MAS::Magnetic magnetic;

    // Core shape: E 19/8/5
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

    // Geometrical description: two half sets
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

    MAS::Coil coil;
    magnetic.set_coil(coil);

    return magnetic;
}

TEST_CASE("E core builds manually", "[core][e]") {
    auto magnetic = make_simple_e_magnetic();

    mvb::MagneticBuilder builder;
    auto pieces = builder.buildCoreNamed(magnetic.get_core());

    REQUIRE(pieces.size() == 2);

    double totalVolume = 0.0;
    for (const auto& ns : pieces) {
        const auto& s = ns.shape;
        REQUIRE(!s.IsNull());

        Bnd_Box box;
        BRepBndLib::Add(s, box);
        double xmin, ymin, zmin, xmax, ymax, zmax;
        box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

        REQUIRE(xmax > xmin);
        REQUIRE(ymax > ymin);
        REQUIRE(zmax > zmin);

        GProp_GProps props;
        BRepGProp::VolumeProperties(s, props);
        totalVolume += props.Mass();
    }

    REQUIRE(totalVolume > 0.0);
}

static MAS::Magnetic make_simple_t_magnetic() {
    MAS::Magnetic magnetic;

    MAS::CoreShape shape;
    shape.set_family(MAS::CoreShapeFamily::T);
    shape.set_type(MAS::FunctionalDescriptionType::STANDARD);
    std::map<std::string, MAS::Dimension> dims;
    dims["A"] = make_dim(0.020); // outer diameter
    dims["B"] = make_dim(0.010); // inner diameter
    dims["C"] = make_dim(0.008); // height
    shape.set_dimensions(dims);

    MAS::CoreGeometricalDescriptionElement piece;
    piece.set_type(MAS::CoreGeometricalDescriptionElementType::TOROIDAL);
    piece.set_coordinates({0.0, 0.0, 0.0});
    piece.set_rotation(std::optional<std::vector<double>>(std::vector<double>{0.0, 0.0, 0.0}));
    piece.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    MAS::MagneticCore core;
    core.set_geometrical_description(std::optional<std::vector<MAS::CoreGeometricalDescriptionElement>>(std::vector<MAS::CoreGeometricalDescriptionElement>{piece}));
    magnetic.set_core(core);

    MAS::Coil coil;
    magnetic.set_coil(coil);

    return magnetic;
}

static MAS::Magnetic make_simple_pq_magnetic() {
    MAS::Magnetic magnetic;

    MAS::CoreShape shape;
    shape.set_family(MAS::CoreShapeFamily::PQ);
    shape.set_type(MAS::FunctionalDescriptionType::STANDARD);
    std::map<std::string, MAS::Dimension> dims;
    // PQ 32/12 — note e*cos_g > c, exercises the self-intersection fix
    dims["A"] = make_dim(0.033);
    dims["B"] = make_dim(0.00594);
    dims["C"] = make_dim(0.022);
    dims["D"] = make_dim(0.0034);
    dims["E"] = make_dim(0.027);
    dims["F"] = make_dim(0.0135);
    shape.set_dimensions(dims);

    MAS::CoreGeometricalDescriptionElement piece1, piece2;
    piece1.set_type(MAS::CoreGeometricalDescriptionElementType::HALF_SET);
    piece1.set_coordinates({0.0, 0.0, 0.0});
    piece1.set_rotation(std::optional<std::vector<double>>(std::vector<double>{std::numbers::pi, std::numbers::pi, 0.0}));
    piece1.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    piece2.set_type(MAS::CoreGeometricalDescriptionElementType::HALF_SET);
    piece2.set_coordinates({0.0, 0.0, 0.0});
    piece2.set_rotation(std::optional<std::vector<double>>(std::vector<double>{0.0, 0.0, 0.0}));
    piece2.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    MAS::MagneticCore core;
    core.set_geometrical_description(std::optional<std::vector<MAS::CoreGeometricalDescriptionElement>>(
        std::vector<MAS::CoreGeometricalDescriptionElement>{piece1, piece2}));
    magnetic.set_core(core);

    MAS::Coil coil;
    magnetic.set_coil(coil);

    return magnetic;
}

TEST_CASE("PQ core builds manually with non-zero volume", "[core][pq]") {
    auto magnetic = make_simple_pq_magnetic();

    mvb::MagneticBuilder builder;
    auto pieces = builder.buildCoreNamed(magnetic.get_core());

    REQUIRE(pieces.size() == 2);

    double totalVolume = 0.0;
    for (const auto& ns : pieces) {
        const auto& s = ns.shape;
        REQUIRE(!s.IsNull());

        Bnd_Box box;
        BRepBndLib::Add(s, box);
        REQUIRE(!box.IsVoid());

        GProp_GProps props;
        BRepGProp::VolumeProperties(s, props);
        totalVolume += props.Mass();
    }

    printf("PQ 32/12 total volume (2 pieces): %.8e m^3\n", totalVolume);
    REQUIRE(totalVolume > 0.0);
    REQUIRE(totalVolume > 5.32e-6);  // Python reference: 5.37779e-6 m^3, within 1%
    REQUIRE(totalVolume < 5.43e-6);
}

TEST_CASE("T core builds manually", "[core][t]") {
    auto magnetic = make_simple_t_magnetic();

    mvb::MagneticBuilder builder;
    auto pieces = builder.buildCoreNamed(magnetic.get_core());

    REQUIRE(pieces.size() == 1);

    const auto& s = pieces[0].shape;
    REQUIRE(!s.IsNull());

    Bnd_Box box;
    BRepBndLib::Add(s, box);
    double xmin, ymin, zmin, xmax, ymax, zmax;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

    REQUIRE(xmax > xmin);
    REQUIRE(ymax > ymin);
    REQUIRE(zmax > zmin);

    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props);
    double volume = props.Mass();
    REQUIRE(volume > 0.0);
}
