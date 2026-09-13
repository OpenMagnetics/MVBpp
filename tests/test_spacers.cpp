// ABT #1170 (WP1 "Spacers end to end"): an additive-gapped core draws its plastic shims.
//
// MKF emits one SPACER element per column of a two-piece set (Core.cpp); MVB++ has had a
// SpacerBuilder since the WASM `_drawSpacer` was added, but nothing called it from
// buildAllNamed, so every STEP export, every 2D section and every FEM mesh of a spacer-gapped
// part was missing the very thing that HOLDS the halves apart.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "mvb/MagneticBuilder.h"
#include "mvb/SpacerBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Core.h"

#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRep_Tool.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <limits>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

double volume_of(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}

// The additive-gap fixture MKF's ABT #1170 test generates: E 42/21/15, 0.5 mm additive gap,
// three columns, autocompleted.
MAS::Magnetic load_additive_fixture() {
    auto path = std::filesystem::path(MAS_COMPLETE_DIR) / "additive_gap_e_core.json";
    REQUIRE(std::filesystem::exists(path));
    std::ifstream file(path);
    json j;
    file >> j;
    return j.at("magnetic").get<MAS::Magnetic>();
}

}  // namespace

TEST_CASE("Additive-gapped E core draws one shim per column", "[spacer][abt1170]") {
    auto magnetic = load_additive_fixture();
    const auto core = magnetic.get_core().value();

    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(magnetic, /*includeBobbin=*/false);

    std::vector<mvb::NamedShape> spacerSolids;
    for (const auto& ns : all) {
        if (ns.name.rfind("Spacer_", 0) == 0) spacerSolids.push_back(ns);
    }
    INFO("solids: " << all.size() << ", spacers: " << spacerSolids.size());
    REQUIRE(spacerSolids.size() == 3u);

    // Every shim must carry its MAS material, and the role WP0's OMFEM classify() keys on.
    for (const auto& ns : spacerSolids) {
        CHECK(ns.role == mvb::Role::Spacer);
        CHECK_FALSE(ns.materialName.empty());
    }

    // ---- volume: width x thickness x depth of the emitting MAS element, exactly.
    // The fixture is functional-only (like its siblings): buildAllNamed has MKF autocomplete it,
    // so the reference elements come from the same MKF processing, not from the file.
    json coreJson = core;
    OpenMagnetics::Core mkfCore(coreJson);
    REQUIRE(mkfCore.get_geometrical_description().has_value());
    auto geometricalDescription = mkfCore.get_geometrical_description().value();
    std::vector<std::vector<double>> expectedBoxes;
    for (const auto& element : geometricalDescription) {
        if (element.get_type() != MAS::CoreGeometricalDescriptionElementType::SPACER) continue;
        expectedBoxes.push_back(element.get_dimensions().value());
    }
    REQUIRE(expectedBoxes.size() == spacerSolids.size());
    for (size_t i = 0; i < spacerSolids.size(); ++i) {
        const double expected = expectedBoxes[i][0] * expectedBoxes[i][1] * expectedBoxes[i][2];
        UNSCOPED_INFO("spacer " << i << " expected volume " << expected << " m^3");
        CHECK_THAT(volume_of(spacerSolids[i].shape),
                   Catch::Matchers::WithinAbs(expected, 1e-9));
    }

    // ---- the shim touches BOTH half-sets: its faces sit exactly on the parting planes.
    // MKF puts the halves at y = +/- thickness/2, so a shim of that thickness centred on
    // y = 0 has zero clearance to either.
    // Extremes are read off the VERTICES, not a Bnd_Box: Bnd_Box pads every bound by the
    // shape tolerance (1e-7 m here), which is 100x the tolerance this check asserts.
    std::vector<mvb::NamedShape> corePieces;
    for (const auto& ns : all) if (ns.role == mvb::Role::Core) corePieces.push_back(ns);
    REQUIRE(corePieces.size() == 2u);   // one two-piece set: the two halves the shim separates
    for (size_t i = 0; i < spacerSolids.size(); ++i) {
        double ylo = std::numeric_limits<double>::max();
        double yhi = std::numeric_limits<double>::lowest();
        for (TopExp_Explorer ex(spacerSolids[i].shape, TopAbs_VERTEX); ex.More(); ex.Next()) {
            const double y = BRep_Tool::Pnt(TopoDS::Vertex(ex.Current())).Y();
            ylo = std::min(ylo, y);
            yhi = std::max(yhi, y);
        }
        // and it TOUCHES each half: minimum distance zero to both core pieces.
        for (const auto& piece : corePieces) {
            BRepExtrema_DistShapeShape dist(spacerSolids[i].shape, piece.shape);
            REQUIRE(dist.IsDone());
            UNSCOPED_INFO("spacer " << i << " to " << piece.name << ": " << dist.Value() << " m");
            CHECK(dist.Value() <= 1e-9);
        }
        const double thickness = expectedBoxes[i][1];
        UNSCOPED_INFO("spacer " << i << " spans y [" << ylo << ", " << yhi << "]");
        CHECK_THAT(ylo, Catch::Matchers::WithinAbs(-thickness / 2, 1e-9));
        CHECK_THAT(yhi, Catch::Matchers::WithinAbs(+thickness / 2, 1e-9));
    }

    // ---- no shim may eat into the winding.
    std::size_t turnSolids = 0;
    for (const auto& ns : all) if (ns.role == mvb::Role::Turn) ++turnSolids;
    REQUIRE(turnSolids > 0);   // otherwise the check below is vacuous
    for (const auto& spacer : spacerSolids) {
        for (const auto& ns : all) {
            if (ns.name.rfind("Spacer_", 0) == 0) continue;
            if (ns.role != mvb::Role::Turn) continue;
            BRepAlgoAPI_Common common(spacer.shape, ns.shape);
            REQUIRE(common.IsDone());   // an unevaluated intersection is not "no intersection"
            UNSCOPED_INFO(spacer.name << " vs " << ns.name);
            CHECK(volume_of(common.Shape()) < 1e-12);
        }
    }
}

TEST_CASE("A core with no additive gap has no shims", "[spacer][abt1170]") {
    auto path = std::filesystem::path(MAS_COMPLETE_DIR) / "buck_inductor_complete.json";
    std::ifstream file(path);
    json j;
    file >> j;
    auto magnetic = j.at("magnetic").get<MAS::Magnetic>();

    mvb::MagneticBuilder builder;
    auto all = builder.buildAllNamed(magnetic, /*includeBobbin=*/false);
    for (const auto& ns : all) {
        CHECK(ns.name.rfind("Spacer_", 0) != 0);
    }
}
