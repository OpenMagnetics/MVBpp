// The core's insulating coating (epoxy/parylene/nylon/glass) in the 3D drawing.
//
// MVB++ has always been able to build the conformal shell -- buildCoreCoatingShell offsets
// the core surface and subtracts the core -- but only the FEM/step path ever passed a
// non-zero thickness. drawMagnetic, the entry point every 3D viewer goes through, hardcoded
// 0.0, so a coated core drew identically to a bare one and a PM editing the coating saw
// nothing change. These tests pin both halves of the fix: a declared coating is drawn, and an
// undeclared one is NOT (MKF's Core::get_coating_thickness() invents a default jacket for any
// uncoated toroid, which is right for the dielectric path it exists for and wrong to draw).

#include <catch2/catch_test_macros.hpp>
#include "mvb/MagneticBuilder.h"
#include "constructive_models/Magnetic.h"
#include "MAS.hpp"

#include <filesystem>
#include <fstream>
#include <map>
#include <numbers>
#include <sstream>
#include <string>

namespace {

MAS::Dimension dim(double value) { return MAS::Dimension(value); }

MAS::MagneticCore make_toroid_core() {
    MAS::CoreShape shape;
    shape.set_family(MAS::CoreShapeFamily::T);
    shape.set_type(MAS::FunctionalDescriptionType::STANDARD);
    std::map<std::string, MAS::Dimension> dimensions;
    dimensions["A"] = dim(0.025);
    dimensions["B"] = dim(0.015);
    dimensions["C"] = dim(0.010);
    shape.set_dimensions(dimensions);

    MAS::CoreGeometricalDescriptionElement piece;
    piece.set_type(MAS::CoreGeometricalDescriptionElementType::TOROIDAL);
    piece.set_coordinates({0.0, 0.0, 0.0});
    piece.set_rotation(std::optional<std::vector<double>>(std::vector<double>{
        std::numbers::pi / 2.0, std::numbers::pi / 2.0, 0.0}));
    piece.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    MAS::CoreFunctionalDescription functional;
    functional.set_material(MAS::CoreMaterialDataOrNameUnion{std::string("N87")});
    functional.set_type(MAS::CoreType::TOROIDAL);
    functional.set_number_stacks(int64_t(1));
    functional.set_shape(shape);

    MAS::MagneticCore core;
    core.set_functional_description(functional);
    core.set_geometrical_description(
        std::optional<std::vector<MAS::CoreGeometricalDescriptionElement>>(
            std::vector<MAS::CoreGeometricalDescriptionElement>{piece}));
    return core;
}

MAS::Magnetic make_bare_toroid_magnetic() {
    MAS::CoreBobbinProcessedDescription bobbinPd;
    bobbinPd.set_column_width(0.0025);
    bobbinPd.set_column_depth(0.005);
    bobbinPd.set_column_thickness(0.0);
    bobbinPd.set_wall_thickness(0.0);
    bobbinPd.set_column_shape(MAS::ColumnShape::RECTANGULAR);
    MAS::WindingWindowElement window;
    window.set_radial_height(0.0075);
    window.set_angle(360.0);
    window.set_coordinates(std::vector<double>({0.0075, 0.0, 0.0}));
    bobbinPd.set_winding_windows(std::vector<MAS::WindingWindowElement>{window});

    MAS::Bobbin bobbin;
    bobbin.set_processed_description(
        std::optional<MAS::CoreBobbinProcessedDescription>(bobbinPd));

    MAS::Wire wire;
    wire.set_type(MAS::WireType::ROUND);
    {
        MAS::DimensionWithTolerance outer;
        outer.set_nominal(0.001);
        wire.set_outer_diameter(std::optional<MAS::DimensionWithTolerance>(outer));
        MAS::DimensionWithTolerance conducting;
        conducting.set_nominal(0.0009);
        wire.set_conducting_diameter(std::optional<MAS::DimensionWithTolerance>(conducting));
    }
    wire.set_name(std::optional<std::string>("Round 0.9 - Grade 1"));

    MAS::Turn turn;
    turn.set_name("Primary parallel 0 turn 0");
    turn.set_winding("Primary");
    turn.set_length(0.05);
    turn.set_parallel(0);
    turn.set_coordinates({0.0, 0.0});
    // The outer XY-plane crossing. MKF emits it for every toroidal turn and MVB++ refuses to
    // invent one, so a hand-built turn has to carry it: the wire threads the bore at the
    // origin and comes back round outside the 12.5 mm outer radius.
    turn.set_additional_coordinates(
        std::optional<std::vector<std::vector<double>>>(
            std::vector<std::vector<double>>{{0.0130, 0.0}}));

    MAS::CoilFunctionalDescription coilFunctional;
    coilFunctional.set_name("Primary");
    coilFunctional.set_number_turns(1);
    coilFunctional.set_number_parallels(1);
    coilFunctional.set_isolation_side(MAS::IsolationSide::PRIMARY);
    coilFunctional.set_wire(MAS::WireDataOrNameUnion(wire));

    MAS::Coil coil;
    coil.set_bobbin(MAS::BobbinDataOrNameUnion(bobbin));
    coil.set_turns_description(std::optional<std::vector<MAS::Turn>>(std::vector<MAS::Turn>{turn}));
    coil.set_functional_description(std::vector<MAS::CoilFunctionalDescription>{coilFunctional});

    MAS::Magnetic magnetic;
    magnetic.set_core(make_toroid_core());
    magnetic.set_coil(coil);
    return magnetic;
}

// STEP names every solid it writes, so the exported text is the honest record of what was
// drawn -- a shell that is built but never emitted would not appear here.
size_t count_coating_solids(const MAS::Magnetic& magnetic, const std::string& label) {
    mvb::MagneticBuilder builder;
    auto directory = std::filesystem::temp_directory_path() / ("mvb_core_coating_" + label);
    std::filesystem::create_directories(directory);
    const auto path = builder.drawMagnetic(magnetic, directory.string(), "step");
    REQUIRE(std::filesystem::exists(path));

    std::ifstream file(path);
    std::stringstream buffer;
    buffer << file.rdbuf();
    const std::string step = buffer.str();

    size_t found = 0;
    for (size_t at = step.find(" coating"); at != std::string::npos; at = step.find(" coating", at + 1)) {
        ++found;
    }
    std::filesystem::remove_all(directory);
    return found;
}

TEST_CASE("An undeclared core coating is not drawn", "[core-coating]") {
    // MKF would happily resolve a default jacket for this uncoated toroid. Drawing it would
    // put a shell on every toroid in every consumer, including designs whose author never
    // mentioned a coating.
    auto magnetic = make_bare_toroid_magnetic();
    REQUIRE(count_coating_solids(magnetic, "bare") == 0);
}

TEST_CASE("A declared core coating is drawn as its own solid", "[core-coating]") {
    auto magnetic = make_bare_toroid_magnetic();
    MAS::CoreCoating coating;
    coating.set_type(std::optional<MAS::CoatingType>(MAS::CoatingType::EPOXY));
    coating.set_thickness(0.0005);
    auto core = magnetic.get_core().value();
    auto functional = core.get_functional_description();
    functional.set_coating(std::optional<MAS::CoreCoatingDataOrNameUnion>(coating));
    core.set_functional_description(functional);
    magnetic.set_core(core);

    REQUIRE(count_coating_solids(magnetic, "epoxy") > 0);
}

TEST_CASE("A magnetic-epoxy shield cap is not drawn as a core coating", "[core-coating]") {
    // It is the powder-loaded cap moulded over the WINDING of a semishielded drum, not an
    // insulating jacket on the ferrite surface, and drawCoreShell already draws it
    // translucently. Same exclusion StrayCapacitance::resolve_core_jacket makes.
    auto magnetic = make_bare_toroid_magnetic();
    MAS::CoreCoating coating;
    coating.set_type(std::optional<MAS::CoatingType>(MAS::CoatingType::MAGNETIC_EPOXY));
    coating.set_thickness(0.0005);
    auto core = magnetic.get_core().value();
    auto functional = core.get_functional_description();
    functional.set_coating(std::optional<MAS::CoreCoatingDataOrNameUnion>(coating));
    core.set_functional_description(functional);
    magnetic.set_core(core);

    REQUIRE(count_coating_solids(magnetic, "magnetic_epoxy") == 0);
}

}  // namespace
