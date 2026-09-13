// ABT #1176 (WP7, MAS-RFC 0015): magnetic shunts drawn as solids of their own.
//
// Fixture tests/realwinding_fixtures/shunt_li2018_e32_magnetic.json (a bare MAS Magnetic, like its
// neighbours; mas_complete_fixtures holds full MAS files only) is the magnetic of MKF
// tests/TestMagneticShunt.cpp, make_planar_shunt_transformer(li_2018_case()) at MKF 00118eea, written
// out with OpenMagnetics::to_json: M. Li, Z. Ouyang, M. A. E. Andersen, "High Frequency LLC Resonant
// Converter with Magnetic Shunt Integrated Planar Transformer", IEEE Trans. Power Electron. 34(3),
// 2019, Table V -- E 32/6/20 3F46, additive gap 2la + tsh = 0.28 mm on every column, 4 primary
// layers of 2 turns and 4 secondary layers of 1 turn (0.07 mm copper, 0.25 mm insulation), and one
// betweenSections sheet 0.1 mm thick spanning the core's width and depth at the mating plane. The
// sheet material is TDK FPC C350 (the grade MKF validates sheets on, owner decision 2026-09-13)
// instead of the paper's IFL04, which is not in the material database; MVB++ only carries the name.
// The MKF fixture's assumptions (track width filling the window, sheet between the halves) apply.
//
// The other geometries (two sheets with gaps to the columns, a segmented sheet, a sheet through a
// winding) are derived from that fixture here, with the column faces read from MKF's processed core
// -- the same derivation MKF's own "sheets with gaps to the columns" test makes.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "mvb/MagneticBuilder.h"
#include "mvb/SectionDrawing.h"
#include "mvb/ShuntBuilder.h"
#include "MAS.hpp"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

using json = nlohmann::json;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

double volume_of(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}

// Exact extents from the vertices: a Bnd_Box is padded by the shape tolerance (1e-7 m), 100x the
// tolerances asserted below.
struct Extents {
    double lo[3] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                    std::numeric_limits<double>::max()};
    double hi[3] = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                    std::numeric_limits<double>::lowest()};
};
Extents extents_of(const TopoDS_Shape& shape) {
    Extents e;
    for (TopExp_Explorer ex(shape, TopAbs_VERTEX); ex.More(); ex.Next()) {
        const gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(ex.Current()));
        const double c[3] = {p.X(), p.Y(), p.Z()};
        for (int i = 0; i < 3; ++i) {
            e.lo[i] = std::min(e.lo[i], c[i]);
            e.hi[i] = std::max(e.hi[i], c[i]);
        }
    }
    return e;
}

json load_li_fixture() {
    auto path = std::filesystem::path(MAS_COMPLETE_DIR).parent_path() / "realwinding_fixtures" /
                "shunt_li2018_e32_magnetic.json";
    REQUIRE(std::filesystem::exists(path));
    std::ifstream file(path);
    json j;
    file >> j;
    return j;
}

std::vector<mvb::NamedShape> with_role(const std::vector<mvb::NamedShape>& all, mvb::Role role) {
    std::vector<mvb::NamedShape> out;
    for (const auto& ns : all)
        if (ns.role == role) out.push_back(ns);
    return out;
}

// Window faces of the +x window, from MKF's processed columns (never re-derived from the shape).
struct WindowFaces {
    double inner = 0.0;   // +x face of the central column
    double outer = 0.0;   // -x face of the +x lateral column
    double depth = 0.0;   // core depth
};
WindowFaces window_faces(const json& magneticJson) {
    OpenMagnetics::Core core(magneticJson.at("core"));
    auto columns = core.get_columns();
    REQUIRE(columns.size() == 3u);
    WindowFaces w;
    w.inner = columns[core.get_main_column_index()].get_width() / 2.0;
    bool found = false;
    for (const auto& column : columns) {
        if (column.get_coordinates()[0] > 0) {
            w.outer = column.get_coordinates()[0] - column.get_width() / 2.0;
            found = true;
        }
    }
    REQUIRE(found);
    w.depth = core.get_depth();
    return w;
}

// A shared volume with any turn or core solid (checked independently of the builder's own gate).
void require_no_overlap_with_turns_and_core(const mvb::NamedShape& shunt,
                                            const std::vector<mvb::NamedShape>& all) {
    std::size_t turns = 0, cores = 0;
    for (const auto& ns : all) {
        if (ns.role != mvb::Role::Turn && ns.role != mvb::Role::Core) continue;
        (ns.role == mvb::Role::Turn ? turns : cores)++;
        BRepAlgoAPI_Common common(shunt.shape, ns.shape);
        REQUIRE(common.IsDone());
        UNSCOPED_INFO(shunt.name << " vs " << ns.name);
        CHECK(volume_of(common.Shape()) < 1e-18);
    }
    REQUIRE(turns > 0);   // otherwise the check is vacuous
    REQUIRE(cores == 2u);
}

}  // namespace

TEST_CASE("Li 2018 E32/6/20: the betweenSections sheet is one Shunt solid of the MAS box",
          "[shunt][abt1176]") {
    const json magneticJson = load_li_fixture();
    const auto& shuntJson = magneticJson.at("shunts").at(0);
    const auto coords = shuntJson.at("coordinates").get<std::vector<double>>();
    const auto dims = shuntJson.at("dimensions").get<std::vector<double>>();

    mvb::MagneticBuilder builder;
    // Both overloads: the MAS one (pre-enriched geometricalDescription) and the MKF one.
    const auto viaMas = builder.buildAllNamed(magneticJson.get<MAS::Magnetic>(), /*includeBobbin=*/false);
    const auto viaMkf = builder.buildAllNamed(OpenMagnetics::Magnetic(magneticJson.get<MAS::Magnetic>()),
                                              /*includeBobbin=*/false);

    for (const auto* all : {&viaMas, &viaMkf}) {
        const auto shunts = with_role(*all, mvb::Role::Shunt);
        REQUIRE(shunts.size() == 1u);
        const auto& shunt = shunts.front();
        CHECK(shunt.name == "Shunt_0");
        CHECK(shunt.materialName == "C350");

        std::size_t solids = 0;
        for (TopExp_Explorer ex(shunt.shape, TopAbs_SOLID); ex.More(); ex.Next()) ++solids;
        CHECK(solids == 1u);

        // No gapToColumns: the sheet runs into every column through its gap, so nothing is cut back.
        CHECK_THAT(volume_of(shunt.shape), WithinRel(dims[0] * dims[1] * dims[2], 1e-12));
        const auto e = extents_of(shunt.shape);
        for (int i = 0; i < 3; ++i) {
            CHECK_THAT(e.lo[i], WithinAbs(coords[i] - dims[i] / 2.0, 1e-9));
            CHECK_THAT(e.hi[i], WithinAbs(coords[i] + dims[i] / 2.0, 1e-9));
        }
        require_no_overlap_with_turns_and_core(shunt, *all);

        // The sheet takes its volume out of every shim it passes through (MKF: g - t of the gap
        // left, t of sheet): each spacer keeps exactly its MAS box minus the box intersection with
        // the sheet. MKF's shims stand proud of the core in depth, so the sheet does not cover them.
        const auto spacers = with_role(*all, mvb::Role::Spacer);
        REQUIRE(spacers.size() == 3u);
        OpenMagnetics::Core core(magneticJson.at("core"));
        const auto geometricalDescription = core.get_geometrical_description();
        REQUIRE(geometricalDescription.has_value());
        std::vector<std::pair<std::vector<double>, std::vector<double>>> spacerBoxes;   // centre, size
        for (const auto& element : *geometricalDescription)
            if (element.get_type() == MAS::CoreGeometricalDescriptionElementType::SPACER)
                spacerBoxes.push_back({element.get_coordinates(), element.get_dimensions().value()});
        REQUIRE(spacerBoxes.size() == 3u);
        for (std::size_t i = 0; i < spacers.size(); ++i) {
            const auto& [c, b] = spacerBoxes[i];
            double overlap = 1.0;
            for (int k = 0; k < 3; ++k) {
                const double lo = std::max(c[k] - b[k] / 2.0, coords[k] - dims[k] / 2.0);
                const double hi = std::min(c[k] + b[k] / 2.0, coords[k] + dims[k] / 2.0);
                overlap *= std::max(0.0, hi - lo);
            }
            REQUIRE(overlap > 0.0);   // the sheet does pass through every shim of this fixture
            const double expected = b[0] * b[1] * b[2] - overlap;
            UNSCOPED_INFO(spacers[i].name << " expected " << expected << " m^3");
            CHECK_THAT(volume_of(spacers[i].shape), WithinRel(expected, 1e-9));
            BRepAlgoAPI_Common common(shunt.shape, spacers[i].shape);
            REQUIRE(common.IsDone());
            CHECK(volume_of(common.Shape()) < 1e-18);
        }
    }
}

TEST_CASE("Two inWindow sheets stop at gapToColumns from the column faces", "[shunt][abt1176]") {
    json magneticJson = load_li_fixture();
    const auto w = window_faces(magneticJson);
    const double thickness = magneticJson.at("shunts").at(0).at("dimensions").at(1).get<double>();
    const double gap = 0.1e-3;
    json sheets = json::array();
    for (int side : {1, -1}) {
        sheets.push_back({{"placement", "inWindow"},
                          {"coordinates", {side * (w.inner + w.outer) / 2.0, 0.0, 0.0}},
                          {"dimensions", {w.outer - w.inner - 2.0 * gap, thickness, w.depth}},
                          {"gapToColumns", {{"inner", gap}, {"outer", gap}}},
                          {"material", "C350"}});
    }
    magneticJson["shunts"] = sheets;

    mvb::MagneticBuilder builder;
    const auto all = builder.buildAllNamed(magneticJson.get<MAS::Magnetic>(), /*includeBobbin=*/true);
    const auto shunts = with_role(all, mvb::Role::Shunt);
    REQUIRE(shunts.size() == 2u);
    CHECK(shunts[0].name == "Shunt_0");
    CHECK(shunts[1].name == "Shunt_1");
    for (std::size_t i = 0; i < 2; ++i) {
        const double side = i == 0 ? 1.0 : -1.0;
        const auto e = extents_of(shunts[i].shape);
        // side coordinate u = side * x runs from the central column outwards
        const double u0 = std::min(side * e.lo[0], side * e.hi[0]);
        const double u1 = std::max(side * e.lo[0], side * e.hi[0]);
        CHECK_THAT(u0 - w.inner, WithinAbs(gap, 1e-9));
        CHECK_THAT(w.outer - u1, WithinAbs(gap, 1e-9));
        CHECK_THAT(volume_of(shunts[i].shape),
                   WithinRel((w.outer - w.inner - 2.0 * gap) * thickness * w.depth, 1e-12));
        require_no_overlap_with_turns_and_core(shunts[i], all);
        // A sheet clear of the columns leaves the shims whole.
        for (const auto& spacer : with_role(all, mvb::Role::Spacer)) {
            BRepAlgoAPI_Common common(shunts[i].shape, spacer.shape);
            REQUIRE(common.IsDone());
            CHECK(volume_of(common.Shape()) < 1e-18);
        }
    }
}

TEST_CASE("A segmented sheet becomes one solid per piece with the stated gaps", "[shunt][abt1176]") {
    json magneticJson = load_li_fixture();
    const auto w = window_faces(magneticJson);
    const double thickness = magneticJson.at("shunts").at(0).at("dimensions").at(1).get<double>();
    const double sideGap = 0.1e-3;
    const double width = w.outer - w.inner - 2.0 * sideGap;
    const std::array<double, 3> gaps = {0.2e-3, 0.3e-3, 0.0};
    const double length = (width - gaps[0] - gaps[1]) / 3.0;
    json segments = json::array();
    for (double g : gaps) segments.push_back({{"length", length}, {"gap", g}});
    magneticJson["shunts"] = json::array({{{"name", "segmented sheet"},
                                           {"placement", "inWindow"},
                                           {"coordinates", {(w.inner + w.outer) / 2.0, 0.0, 0.0}},
                                           {"dimensions", {width, thickness, w.depth}},
                                           {"gapToColumns", {{"inner", sideGap}, {"outer", sideGap}}},
                                           {"segments", segments},
                                           {"material", "C350"}}});

    mvb::MagneticBuilder builder;
    const auto all = builder.buildAllNamed(magneticJson.get<MAS::Magnetic>(), /*includeBobbin=*/false);
    const auto pieces = with_role(all, mvb::Role::Shunt);
    REQUIRE(pieces.size() == 3u);
    double previousEnd = 0.0;
    for (std::size_t k = 0; k < pieces.size(); ++k) {
        CHECK(pieces[k].name == "Shunt_0_" + std::to_string(k));   // the MAS name is not used
        CHECK(pieces[k].materialName == "C350");
        const auto e = extents_of(pieces[k].shape);
        CHECK_THAT(e.hi[0] - e.lo[0], WithinAbs(length, 1e-9));
        CHECK_THAT(volume_of(pieces[k].shape), WithinRel(length * thickness * w.depth, 1e-12));
        if (k == 0) CHECK_THAT(e.lo[0] - w.inner, WithinAbs(sideGap, 1e-9));
        else        CHECK_THAT(e.lo[0] - previousEnd, WithinAbs(gaps[k - 1], 1e-9));
        previousEnd = e.hi[0];
        require_no_overlap_with_turns_and_core(pieces[k], all);
    }
    CHECK_THAT(w.outer - previousEnd, WithinAbs(sideGap, 1e-9));
}

TEST_CASE("A sheet laid through a winding throws, naming the shunt and the turn", "[shunt][abt1176]") {
    json magneticJson = load_li_fixture();
    const auto w = window_faces(magneticJson);
    const double thickness = magneticJson.at("shunts").at(0).at("dimensions").at(1).get<double>();
    // The centre height of the first primary turn in the +x window, from MKF's turns description.
    const auto& turns = magneticJson.at("coil").at("turnsDescription");
    double turnY = std::numeric_limits<double>::quiet_NaN();
    std::string turnName;
    for (const auto& t : turns) {
        if (t.at("winding") == "Primary" && t.at("coordinates").at(0).get<double>() > 0) {
            turnY = t.at("coordinates").at(1).get<double>();
            turnName = t.at("name").get<std::string>();
            break;
        }
    }
    REQUIRE(std::isfinite(turnY));
    const double gap = 0.1e-3;
    magneticJson["shunts"] = json::array({{{"placement", "inWindow"},
                                           {"coordinates", {(w.inner + w.outer) / 2.0, turnY, 0.0}},
                                           {"dimensions", {w.outer - w.inner - 2.0 * gap, thickness, w.depth}},
                                           {"gapToColumns", {{"inner", gap}, {"outer", gap}}},
                                           {"material", "C350"}}});

    mvb::MagneticBuilder builder;
    INFO("turn " << turnName << " at y = " << turnY);
    REQUIRE_THROWS_WITH(builder.buildAllNamed(magneticJson.get<MAS::Magnetic>(), /*includeBobbin=*/false),
                        ContainsSubstring("Shunt collision gate: Shunt_0 intersects") &&
                            ContainsSubstring("(turn)"));
}

TEST_CASE("A shunt MKF refuses is not drawn: its exception propagates", "[shunt][abt1176]") {
    json magneticJson = load_li_fixture();
    mvb::MagneticBuilder builder;
    SECTION("outsideWindow") {
        magneticJson["shunts"][0]["placement"] = "outsideWindow";
        CHECK_THROWS_WITH(builder.buildAllNamed(magneticJson.get<MAS::Magnetic>(), false),
                          ContainsSubstring("placed outsideWindow"));
    }
    SECTION("a sheet short of the columns without gapToColumns") {
        const auto w = window_faces(magneticJson);
        magneticJson["shunts"][0]["coordinates"] = {(w.inner + w.outer) / 2.0, 0.0, 0.0};
        magneticJson["shunts"][0]["dimensions"][0] = (w.outer - w.inner) / 2.0;
        CHECK_THROWS_WITH(builder.buildAllNamed(magneticJson.get<MAS::Magnetic>(), false),
                          ContainsSubstring("gapToColumns"));
    }
}

TEST_CASE("Shunts are additional solids: without them the assembly is unchanged", "[shunt][abt1176]") {
    json withShunt = load_li_fixture();
    json without = withShunt;
    without.erase("shunts");

    mvb::MagneticBuilder builder;
    const auto a = builder.buildAllNamed(withShunt.get<MAS::Magnetic>(), /*includeBobbin=*/true);
    const auto b = builder.buildAllNamed(without.get<MAS::Magnetic>(), /*includeBobbin=*/true);
    CHECK(with_role(b, mvb::Role::Shunt).empty());

    // Same solids in the same order, the shunt aside; identical volumes except the shims the sheet
    // passes through.
    std::vector<const mvb::NamedShape*> rest;
    for (const auto& ns : a)
        if (ns.role != mvb::Role::Shunt) rest.push_back(&ns);
    REQUIRE(rest.size() == b.size());
    for (std::size_t i = 0; i < b.size(); ++i) {
        CHECK(rest[i]->name == b[i].name);
        CHECK(rest[i]->role == b[i].role);
        if (b[i].role == mvb::Role::Spacer) continue;
        if (b[i].role == mvb::Role::FR4) {
            // The planar board is one slab over the whole printed group, so the sheet at the
            // mating plane takes exactly its shared volume out of it.
            BRepAlgoAPI_Common common(b[i].shape, with_role(a, mvb::Role::Shunt).at(0).shape);
            REQUIRE(common.IsDone());
            const double shared = volume_of(common.Shape());
            CHECK(shared > 0.0);
            CHECK_THAT(volume_of(rest[i]->shape), WithinRel(volume_of(b[i].shape) - shared, 1e-9));
            continue;
        }
        CHECK(volume_of(rest[i]->shape) == volume_of(b[i].shape));
    }

    // The 2D gapping drawing gains the sheet's outline and loses nothing.
    OpenMagnetics::Magnetic mkfWith(withShunt.get<MAS::Magnetic>());
    OpenMagnetics::Magnetic mkfWithout(without.get<MAS::Magnetic>());
    auto count_paths = [](const std::string& svg) {
        std::size_t n = 0;
        for (auto pos = svg.find("<path"); pos != std::string::npos; pos = svg.find("<path", pos + 1)) ++n;
        return n;
    };
    const auto svgWith = mvb::SectionDrawing::drawCoreGappingTechnicalDrawing(mkfWith);
    const auto svgWithout = mvb::SectionDrawing::drawCoreGappingTechnicalDrawing(mkfWithout);
    CHECK(count_paths(svgWith) > count_paths(svgWithout));
}
