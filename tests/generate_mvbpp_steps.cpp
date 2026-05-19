#include "mvb/MagneticBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <numbers>
#include <vector>
#include <map>
#include <filesystem>
#include <iostream>

static MAS::Dimension make_dim(double v) { return MAS::Dimension(v); }

static MAS::Magnetic make_simple_e_magnetic() {
    MAS::Magnetic magnetic;
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
    magnetic.set_coil(MAS::Coil());
    return magnetic;
}

static MAS::Magnetic make_simple_t_magnetic() {
    MAS::Magnetic magnetic;
    MAS::CoreShape shape;
    shape.set_family(MAS::CoreShapeFamily::T);
    shape.set_type(MAS::FunctionalDescriptionType::STANDARD);
    std::map<std::string, MAS::Dimension> dims;
    dims["A"] = make_dim(0.020);
    dims["B"] = make_dim(0.010);
    dims["C"] = make_dim(0.008);
    shape.set_dimensions(dims);

    MAS::CoreGeometricalDescriptionElement piece;
    piece.set_type(MAS::CoreGeometricalDescriptionElementType::TOROIDAL);
    piece.set_coordinates({0.0, 0.0, 0.0});
    piece.set_rotation(std::optional<std::vector<double>>(std::vector<double>{0.0, 0.0, 0.0}));
    piece.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    MAS::MagneticCore core;
    core.set_geometrical_description(std::optional<std::vector<MAS::CoreGeometricalDescriptionElement>>(std::vector<MAS::CoreGeometricalDescriptionElement>{piece}));
    magnetic.set_core(core);
    magnetic.set_coil(MAS::Coil());
    return magnetic;
}

static void exportCore(const std::string& outPath, const MAS::Magnetic& magnetic) {
    mvb::MagneticBuilder builder;
    auto named = builder.buildCoreNamed(magnetic.get_core());

    gp_Trsf trsf;
    trsf.SetScale(gp_Pnt(0, 0, 0), 1000.0);
    for (auto& ns : named) {
        if (!ns.shape.IsNull()) {
            ns.shape = BRepBuilderAPI_Transform(ns.shape, trsf).Shape();
        }
    }

    mvb::exportSTEP(named, outPath);
}

int main(int argc, char* argv[]) {
    std::string outDir = (argc > 1) ? argv[1] : ".";
    std::filesystem::create_directories(outDir);

    std::string ePath = std::filesystem::path(outDir) / "mvbpp_e_core.step";
    exportCore(ePath, make_simple_e_magnetic());
    std::cout << "MVB++ E core: " << ePath << "\n";

    std::string tPath = std::filesystem::path(outDir) / "mvbpp_t_core.step";
    exportCore(tPath, make_simple_t_magnetic());
    std::cout << "MVB++ T core: " << tPath << "\n";

    return 0;
}
