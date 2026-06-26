#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "mvb/StepExporter.h"
#include <STEPControl_Reader.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <numbers>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;
using json = nlohmann::json;

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

    MAS::CoreFunctionalDescription coreFunc;
    coreFunc.set_material(MAS::CoreMaterialDataOrNameUnion{std::string("N87")});
    coreFunc.set_type(MAS::CoreType::TWO_PIECE_SET);
    coreFunc.set_number_stacks(int64_t(1));
    coreFunc.set_shape(shape);

    MAS::MagneticCore core;
    core.set_functional_description(coreFunc);
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

    MAS::CoreFunctionalDescription coreFunc;
    coreFunc.set_material(MAS::CoreMaterialDataOrNameUnion{std::string("N87")});
    coreFunc.set_type(MAS::CoreType::TOROIDAL);
    coreFunc.set_number_stacks(int64_t(1));
    coreFunc.set_shape(shape);

    MAS::MagneticCore core;
    core.set_functional_description(coreFunc);
    core.set_geometrical_description(std::optional<std::vector<MAS::CoreGeometricalDescriptionElement>>(std::vector<MAS::CoreGeometricalDescriptionElement>{piece}));
    magnetic.set_core(core);
    magnetic.set_coil(MAS::Coil());
    return magnetic;
}

struct StepStats {
    double volume = 0.0;
    double xmin = 0.0, ymin = 0.0, zmin = 0.0;
    double xmax = 0.0, ymax = 0.0, zmax = 0.0;
    int solidCount = 0;
};

static TopoDS_Shape loadSTEP(const std::string& path) {
    STEPControl_Reader reader;
    IFSelect_ReturnStatus status = reader.ReadFile(path.c_str());
    REQUIRE(status == IFSelect_RetDone);
    reader.TransferRoots();
    return reader.OneShape();
}

static StepStats analyzeShape(const TopoDS_Shape& shape) {
    StepStats stats;
    REQUIRE(!shape.IsNull());

    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    stats.volume = props.Mass();

    Bnd_Box box;
    BRepBndLib::Add(shape, box);
    box.Get(stats.xmin, stats.ymin, stats.zmin, stats.xmax, stats.ymax, stats.zmax);

    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        ++stats.solidCount;
    }
    return stats;
}

static void exportAndCompare(const std::string& refPath,
                              const std::vector<TopoDS_Shape>& shapes,
                              const std::vector<std::string>& names,
                              const std::string& testName,
                              double volumeTol = 0.05,
                              double dimTol = 0.05) {
    REQUIRE(std::filesystem::exists(refPath));
    REQUIRE(!shapes.empty());

    StepStats ref = analyzeShape(loadSTEP(refPath));

    std::string outDir = std::filesystem::temp_directory_path().string();
    std::string generated = outDir + "/" + testName + "_generated.step";
    REQUIRE(mvb::exportSTEP(shapes, names, generated));

    REQUIRE(std::filesystem::exists(generated));
    StepStats gen = analyzeShape(loadSTEP(generated));

    INFO("Reference volume: " << ref.volume << ", Generated volume: " << gen.volume);
    INFO("Reference solids: " << ref.solidCount << ", Generated solids: " << gen.solidCount);
    INFO("Reference bbox: [" << ref.xmin << ", " << ref.ymin << ", " << ref.zmin << "] - ["
         << ref.xmax << ", " << ref.ymax << ", " << ref.zmax << "]");
    INFO("Generated bbox: [" << gen.xmin << ", " << gen.ymin << ", " << gen.zmin << "] - ["
         << gen.xmax << ", " << gen.ymax << ", " << gen.zmax << "]");

    // Volume tolerance
    REQUIRE_THAT(gen.volume, WithinRel(ref.volume, volumeTol));

    // Bounding box size comparison (sort dimensions to be orientation-agnostic)
    double refDims[3] = {ref.xmax - ref.xmin, ref.ymax - ref.ymin, ref.zmax - ref.zmin};
    double genDims[3] = {gen.xmax - gen.xmin, gen.ymax - gen.ymin, gen.zmax - gen.zmin};
    std::sort(std::begin(refDims), std::end(refDims));
    std::sort(std::begin(genDims), std::end(genDims));

    REQUIRE_THAT(genDims[0], WithinRel(refDims[0], dimTol));
    REQUIRE_THAT(genDims[1], WithinRel(refDims[1], dimTol));
    REQUIRE_THAT(genDims[2], WithinRel(refDims[2], dimTol));

    std::filesystem::remove(generated);
}

static void compareSTEPAgainstReference(const std::string& refPath,
                                         const MAS::Magnetic& magnetic,
                                         const std::string& testName) {
    mvb::MagneticBuilder builder;
    auto coreNamed = builder.buildCoreNamed(magnetic.get_core().value());
    REQUIRE(!coreNamed.empty());

    std::vector<TopoDS_Shape> coreShapes;
    coreShapes.reserve(coreNamed.size());
    for (auto& ns : coreNamed) coreShapes.push_back(ns.shape);

    // Scale to mm to match Python MVB export convention
    {
        gp_Trsf trsf;
        trsf.SetScale(gp_Pnt(0, 0, 0), 1000.0);
        for (auto& s : coreShapes) {
            s = BRepBuilderAPI_Transform(s, trsf).Shape();
        }
    }

    std::vector<std::string> names;
    for (size_t i = 0; i < coreShapes.size(); ++i) {
        names.push_back("Core_" + std::to_string(i));
    }

    exportAndCompare(refPath, coreShapes, names, testName);
}

static bool isShapeUsable(const TopoDS_Shape& shape) {
    if (shape.IsNull()) return false;
    try {
        Bnd_Box bb;
        BRepBndLib::Add(shape, bb);
        return !bb.IsVoid();
    } catch (...) {
        return false;
    }
}

static TopoDS_Shape cutBobbin(const TopoDS_Shape& bobbin, const std::vector<TopoDS_Shape>& cutters) {
    TopoDS_Shape result = bobbin;
    if (result.IsNull()) return result;
    for (const auto& tool : cutters) {
        if (tool.IsNull()) continue;
        BRepAlgoAPI_Cut cutter(result, tool);
        if (cutter.IsDone() && !cutter.Shape().IsNull()) {
            TopoDS_Shape candidate = cutter.Shape();
            if (isShapeUsable(candidate)) {
                result = candidate;
            }
        }
    }
    return result;
}

static void compareFullAssemblyAgainstReference(const std::string& refPath,
                                                 const MAS::Magnetic& magnetic,
                                                 const std::string& testName) {
    mvb::MagneticBuilder builder;
    std::vector<mvb::NamedShape> coreNamed;
    mvb::NamedShape bobbinNamed;
    std::vector<mvb::NamedShape> turnsNamed;
    try {
        coreNamed = builder.buildCoreNamed(magnetic.get_core().value());
    } catch (const std::exception& e) {
        FAIL("buildCoreNamed exception: " << e.what());
    } catch (...) {
        FAIL("buildCoreNamed unknown exception");
    }
    try {
        bobbinNamed = builder.buildBobbinNamed(magnetic.get_coil().value(), magnetic.get_core().value());
    } catch (const std::exception& e) {
        FAIL("buildBobbinNamed exception: " << e.what());
    } catch (...) {
        FAIL("buildBobbinNamed unknown exception");
    }
    try {
        turnsNamed = builder.buildTurnsNamed(magnetic.get_coil().value(), magnetic.get_core().value());
    } catch (const std::exception& e) {
        FAIL("buildTurnsNamed exception: " << e.what());
    } catch (...) {
        FAIL("buildTurnsNamed unknown exception");
    }

    REQUIRE(!coreNamed.empty());

    std::vector<TopoDS_Shape> coreShapes;
    coreShapes.reserve(coreNamed.size());
    for (auto& ns : coreNamed) coreShapes.push_back(ns.shape);
    TopoDS_Shape bobbinShape = bobbinNamed.shape;
    std::vector<TopoDS_Shape> turnShapes;
    turnShapes.reserve(turnsNamed.size());
    for (auto& ns : turnsNamed) turnShapes.push_back(ns.shape);

    // Cut bobbin with cores and turns to match Python MVB behavior
    std::vector<TopoDS_Shape> cutters = coreShapes;
    cutters.insert(cutters.end(), turnShapes.begin(), turnShapes.end());
    if (!bobbinShape.IsNull()) {
        bobbinShape = cutBobbin(bobbinShape, cutters);
    }

    std::vector<TopoDS_Shape> allShapes;
    std::vector<std::string> allNames;

    // Scale to mm to match Python MVB export convention
    try {
        gp_Trsf trsf;
        trsf.SetScale(gp_Pnt(0, 0, 0), 1000.0);
        for (auto& s : coreShapes) {
            allShapes.push_back(BRepBuilderAPI_Transform(s, trsf).Shape());
            allNames.push_back("Core_" + std::to_string(allNames.size()));
        }
        if (!bobbinShape.IsNull()) {
            allShapes.push_back(BRepBuilderAPI_Transform(bobbinShape, trsf).Shape());
            allNames.push_back("Bobbin");
        }
        for (auto& s : turnShapes) {
            allShapes.push_back(BRepBuilderAPI_Transform(s, trsf).Shape());
            allNames.push_back("Turn_" + std::to_string(allNames.size()));
        }
    } catch (const std::exception& e) {
        FAIL("Scale/transform exception: " << e.what());
    } catch (...) {
        FAIL("Scale/transform unknown exception");
    }

    try {
        exportAndCompare(refPath, allShapes, allNames, testName);
    } catch (const std::exception& e) {
        FAIL("exportAndCompare exception: " << e.what());
    } catch (...) {
        FAIL("exportAndCompare unknown exception");
    }
}

static MAS::Magnetic loadMagneticFromJson(const std::string& jsonPath) {
    std::ifstream f(jsonPath);
    REQUIRE(f.is_open());
    json j;
    try {
        f >> j;
    } catch (const std::exception& e) {
        FAIL("JSON parse error: " << e.what());
    }
    try {
        mvb::patch_dimension_nominals(j);
    } catch (const std::exception& e) {
        FAIL("patch_dimension_nominals error: " << e.what());
    }
    try {
        return j.at("magnetic").get<MAS::Magnetic>();
    } catch (const std::exception& e) {
        FAIL("MAS magnetic parse error: " << e.what());
    }
    return MAS::Magnetic(); // unreachable
}

static OpenMagnetics::Magnetic loadEnrichedMagneticFromJson(const std::string& jsonPath) {
    std::ifstream f(jsonPath);
    REQUIRE(f.is_open());
    json j;
    try {
        f >> j;
    } catch (const std::exception& e) {
        FAIL("JSON parse error: " << e.what());
    }
    try {
        mvb::patch_dimension_nominals(j);
    } catch (const std::exception& e) {
        FAIL("patch_dimension_nominals error: " << e.what());
    }

    MAS::Magnetic magnetic;
    try {
        magnetic = j.at("magnetic").get<MAS::Magnetic>();
    } catch (const std::exception& e) {
        FAIL("MAS magnetic parse error: " << e.what());
    }

    OpenMagnetics::Magnetic enriched;
    try {
        enriched = mvb::magnetic_autocomplete_safe(magnetic);
    } catch (const std::exception& e) {
        FAIL("magnetic_autocomplete_safe error: " << e.what());
    }

    // After enrichment, ensure bobbin functional description has flange dimensions
    // matching the Python StandardBobbin defaults (0.001, 0.002). This keeps the
    // MVB++ output aligned with the Python reference used in comparison tests.
    auto& coil = enriched.get_mutable_coil();
    auto& bobbinVar = coil.get_mutable_bobbin();
    if (std::holds_alternative<OpenMagnetics::Bobbin>(bobbinVar)) {
        auto& bobbin = std::get<OpenMagnetics::Bobbin>(bobbinVar);
        auto funcDesc = bobbin.get_functional_description();
        if (!funcDesc) {
            funcDesc = MAS::BobbinFunctionalDescription();
        }
        auto dims = funcDesc->get_dimensions();
        if (dims.find("flangeThickness") == dims.end()) {
            dims["flangeThickness"] = 0.001;
        }
        if (dims.find("flangeExtension") == dims.end()) {
            dims["flangeExtension"] = 0.002;
        }
        funcDesc->set_dimensions(dims);
        bobbin.set_functional_description(funcDesc);
    }

    return enriched;
}

static void compareFullAssemblyAgainstReference(const std::string& refPath,
                                                 const OpenMagnetics::Magnetic& magnetic,
                                                 const std::string& testName) {
    mvb::MagneticBuilder builder;
    std::vector<mvb::NamedShape> coreNamed;
    mvb::NamedShape bobbinNamed;
    std::vector<mvb::NamedShape> turnsNamed;
    try {
        coreNamed = builder.buildCoreNamed(magnetic.get_core());
    } catch (const std::exception& e) {
        FAIL("buildCoreNamed exception: " << e.what());
    } catch (...) {
        FAIL("buildCoreNamed unknown exception");
    }
    try {
        bobbinNamed = builder.buildBobbinNamed(magnetic.get_coil(), magnetic.get_core());
    } catch (const std::exception& e) {
        FAIL("buildBobbinNamed exception: " << e.what());
    } catch (...) {
        FAIL("buildBobbinNamed unknown exception");
    }
    try {
        turnsNamed = builder.buildTurnsNamed(magnetic.get_coil(), magnetic.get_core());
    } catch (const std::exception& e) {
        FAIL("buildTurnsNamed exception: " << e.what());
    } catch (...) {
        FAIL("buildTurnsNamed unknown exception");
    }

    REQUIRE(!coreNamed.empty());

    std::vector<TopoDS_Shape> coreShapes;
    coreShapes.reserve(coreNamed.size());
    for (auto& ns : coreNamed) coreShapes.push_back(ns.shape);
    TopoDS_Shape bobbinShape = bobbinNamed.shape;
    std::vector<TopoDS_Shape> turnShapes;
    turnShapes.reserve(turnsNamed.size());
    for (auto& ns : turnsNamed) turnShapes.push_back(ns.shape);

    std::vector<TopoDS_Shape> cutters = coreShapes;
    cutters.insert(cutters.end(), turnShapes.begin(), turnShapes.end());
    if (!bobbinShape.IsNull()) {
        bobbinShape = cutBobbin(bobbinShape, cutters);
    }

    std::vector<TopoDS_Shape> allShapes;
    std::vector<std::string> allNames;

    try {
        gp_Trsf trsf;
        trsf.SetScale(gp_Pnt(0, 0, 0), 1000.0);
        for (auto& s : coreShapes) {
            allShapes.push_back(BRepBuilderAPI_Transform(s, trsf).Shape());
            allNames.push_back("Core_" + std::to_string(allNames.size()));
        }
        if (!bobbinShape.IsNull()) {
            allShapes.push_back(BRepBuilderAPI_Transform(bobbinShape, trsf).Shape());
            allNames.push_back("Bobbin");
        }
        for (auto& s : turnShapes) {
            allShapes.push_back(BRepBuilderAPI_Transform(s, trsf).Shape());
            allNames.push_back("Turn_" + std::to_string(allNames.size()));
        }
    } catch (const std::exception& e) {
        FAIL("Scale/transform exception: " << e.what());
    } catch (...) {
        FAIL("Scale/transform unknown exception");
    }

    try {
        exportAndCompare(refPath, allShapes, allNames, testName);
    } catch (const std::exception& e) {
        FAIL("exportAndCompare exception: " << e.what());
    } catch (...) {
        FAIL("exportAndCompare unknown exception");
    }
}

TEST_CASE("E core STEP matches Python MVB reference", "[step][e]") {
    auto magnetic = make_simple_e_magnetic();
    std::string refPath = std::string(REFERENCE_STEPS_DIR) + "/reference_e_core_core.step";
    compareSTEPAgainstReference(refPath, magnetic, "e_core");
}

TEST_CASE("T core STEP matches Python MVB reference", "[step][t]") {
    auto magnetic = make_simple_t_magnetic();
    std::string refPath = std::string(REFERENCE_STEPS_DIR) + "/reference_t_core_core.step";
    compareSTEPAgainstReference(refPath, magnetic, "t_core");
}

TEST_CASE("Rectangular one-turn assembly STEP matches Python MVB reference", "[step][assembly][json]") {
    auto magnetic = loadEnrichedMagneticFromJson("testData/concentric_rectangular_column_one_turn.json");
    std::string refPath = std::string(REFERENCE_STEPS_DIR) + "/reference_rect_one_turn.step";
    compareFullAssemblyAgainstReference(refPath, magnetic, "rect_one_turn");
}

TEST_CASE("ETD49 5-turn assembly STEP matches Python MVB reference", "[step][assembly][json]") {
    auto magnetic = loadEnrichedMagneticFromJson("testData/ETD49_N87_10uH_5T.json");
    std::string refPath = std::string(REFERENCE_STEPS_DIR) + "/reference_etd49_5t.step";
    compareFullAssemblyAgainstReference(refPath, magnetic, "etd49_5t");
}

// ── STL export tests ────────────────────────────────────────────────────────

TEST_CASE("E core STL export produces non-empty binary data", "[stl][e]") {
    auto magnetic = make_simple_e_magnetic();
    mvb::MagneticBuilder builder;
    auto tmpDir = std::filesystem::temp_directory_path();
    std::string outDir = (tmpDir / "mvb_stl_test_e").string();
    std::filesystem::create_directories(outDir);
    std::string path = builder.drawMagnetic(magnetic, outDir, "stl");
    REQUIRE(std::filesystem::exists(path));
    auto sz = std::filesystem::file_size(path);
    CHECK(sz > 80);  // STL binary header is 80 bytes
    std::filesystem::remove_all(outDir);
}

TEST_CASE("Rectangular one-turn STL export produces non-empty data", "[stl][assembly][json]") {
    auto magnetic = loadEnrichedMagneticFromJson("testData/concentric_rectangular_column_one_turn.json");
    mvb::MagneticBuilder builder;
    auto tmpDir = std::filesystem::temp_directory_path();
    std::string outDir = (tmpDir / "mvb_stl_test_rect").string();
    std::filesystem::create_directories(outDir);
    std::string path = builder.drawMagnetic(magnetic, outDir, "stl");
    REQUIRE(std::filesystem::exists(path));
    auto sz = std::filesystem::file_size(path);
    CHECK(sz > 80);
    std::filesystem::remove_all(outDir);
}

TEST_CASE("ETD49 5-turn STL export with mm scale", "[stl][assembly][json]") {
    auto magnetic = loadEnrichedMagneticFromJson("testData/ETD49_N87_10uH_5T.json");
    mvb::MagneticBuilder builder;
    auto tmpDir = std::filesystem::temp_directory_path();
    std::string outDir = (tmpDir / "mvb_stl_test_etd49").string();
    std::filesystem::create_directories(outDir);
    std::string path = builder.drawMagnetic(magnetic, outDir, "stl",
                                             /*includeBobbin=*/true, /*scale=*/1000.0);
    REQUIRE(std::filesystem::exists(path));
    auto sz = std::filesystem::file_size(path);
    CHECK(sz > 80);
    std::filesystem::remove_all(outDir);
}
