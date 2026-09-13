// ABT #1173 ([toroidbase], WP4): a toroid seated on a catalogue base (MAS bobbin family t with
// functionalDescription.base) draws the base and its pins. The terminal drops are NOT re-aimed at the
// pins here: toroid terminal orientation is being reworked in its own ticket (Alf, 2026-09-13), so the
// drops keep today's construction and the base sits under their common terminal plane.
//
// The bases are the REAL catalogue records (MAS b7cf597, Lodestone Pacific HTM460-6 / HTM600-6 /
// HTM850-6: flat plates, height = standoff = 2.032 mm, 3+3 rectangular tabs). The seat is chosen with
// MKF's find_toroid_bases_for_core, taking the smallest base that holds the coated ring.
//   - buck: buck_inductor_complete as is (T 10/6/4, 8 turns x 3 parallels of Round 0.63): HTM460-6.
//   - CMC: common_mode_choke_complete's T 25.3/14.8/10 coated exceeds every HTM base's 25.4 mm limit,
//     so the CMC case is that fixture on a T 16/9.6/6.3 N30 with two windings of 10 turns of Round 0.80
//     (the design is derived for the test and says so): HTM600-6.
// Geometry frame: the builder's toroid frame (ring axis Y, ring in XZ, drops along -Y). buildAllNamed
// rotates a toroidal assembly by -pi/2 about X at the end; the solids are rotated back here.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "mvb/BaseBuilder.h"
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <numbers>
#include <set>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

json load_complete(const std::string& file) {
    auto path = std::filesystem::path(MAS_COMPLETE_DIR) / file;
    REQUIRE(std::filesystem::exists(path));
    std::ifstream in(path);
    json j;
    in >> j;
    return j.at("magnetic");
}

double volume_of(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}

struct Box { double xlo, ylo, zlo, xhi, yhi, zhi; };
Box box_of(const TopoDS_Shape& shape) {
    Bnd_Box box;
    BRepBndLib::AddOptimal(shape, box, /*useTriangulation=*/false, /*useShapeTolerance=*/false);
    Box b{};
    box.Get(b.xlo, b.ylo, b.zlo, b.xhi, b.yhi, b.zhi);
    return b;
}

// Seat the magnetic's ring on the smallest catalogue base that holds it (MKF's selection), and
// return the magnetic json with that seated bobbin.
json seat_on_smallest_base(json magneticJson, std::string& reference) {
    OpenMagnetics::Core core(magneticJson.at("core"));
    auto bases = OpenMagnetics::find_toroid_bases_for_core(core);
    REQUIRE(!bases.empty());
    std::stable_sort(bases.begin(), bases.end(), [](const OpenMagnetics::Bobbin& a, const OpenMagnetics::Bobbin& b) {
        return a.get_base().get_maximum_core_outer_diameter().value() < b.get_base().get_maximum_core_outer_diameter().value();
    });
    reference = bases.front().get_manufacturer_info()->get_reference().value();
    auto seated = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, bases.front());
    json bobbinJson;
    to_json(bobbinJson, seated);
    magneticJson["coil"]["bobbin"] = bobbinJson;
    return magneticJson;
}

json cmc_on_t16() {
    json magneticJson = load_complete("common_mode_choke_complete.json");
    magneticJson["core"].erase("name");
    magneticJson["core"]["functionalDescription"]["shape"] = "T 16/9.6/6.3";
    for (auto& winding : magneticJson["coil"]["functionalDescription"]) {
        winding["numberTurns"] = 10;
        winding["wire"] = "Round 0.80 - Grade 1";
    }
    return magneticJson;
}

// The whole contract on one seated design.
void require_seated(const json& magneticJson, size_t expectedPins, size_t expectedTips) {
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    const auto& bobbinVariant = enriched.get_coil().get_bobbin();
    const auto* bobbin = std::get_if<OpenMagnetics::Bobbin>(&bobbinVariant);
    REQUIRE(bobbin != nullptr);
    REQUIRE(bobbin->has_base());
    const std::string bobbinName = bobbin->get_name().value();
    const auto base = bobbin->get_base();
    const double height = OpenMagnetics::resolve_dimensional_values(base.get_height());
    const double length = OpenMagnetics::resolve_dimensional_values(base.get_length());
    const double width = OpenMagnetics::resolve_dimensional_values(base.get_width());
    const auto ring = OpenMagnetics::flatten_dimensions(bobbin->get_functional_description()->get_dimensions());
    std::map<std::string, MAS::Pin> pins;
    // By value: the getters hand the optionals back by value, so a range-for over .value() of the
    // temporary would iterate a destroyed vector.
    const auto processedPins = bobbin->get_processed_description()->get_pins().value();
    for (const auto& pin : processedPins) pins.emplace(pin.get_name().value(), pin);
    REQUIRE(pins.size() == expectedPins);

    mvb::MagneticBuilder builder;
    std::vector<mvb::ConductorBuilder::PathPolyline> paths;
    REQUIRE_NOTHROW(paths = builder.buildRealWindingPaths(enriched));
    std::vector<mvb::NamedShape> all;
    REQUIRE_NOTHROW(all = builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                                mvb::DEFAULT_WIRE_POLYGON_SEGMENTS, mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                                /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                                /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                                /*useRealWindingGeometry=*/true, /*femReady=*/true));
    for (auto& ns : all) ns.shape = mvb::rotate_shape(ns.shape, std::numbers::pi / 2.0, 0.0, 0.0);

    const mvb::NamedShape* baseSolid = nullptr;
    std::map<std::string, const mvb::NamedShape*> pinSolids;
    for (const auto& ns : all) {
        if (ns.role == mvb::Role::Base) {
            REQUIRE(baseSolid == nullptr);
            baseSolid = &ns;
        }
        if (ns.role == mvb::Role::Pin) {
            const std::string prefix = bobbinName + " pin ";
            REQUIRE(ns.name.rfind(prefix, 0) == 0);
            pinSolids[ns.name.substr(prefix.size())] = &ns;
        }
    }
    REQUIRE(baseSolid != nullptr);
    CHECK(baseSolid->name == bobbinName + " base");
    REQUIRE(pinSolids.size() == expectedPins);

    // The terminal plane: every tip shares it.
    REQUIRE(paths.size() * 2 == expectedTips);
    const double plane = paths.front().end0[1];
    const Box b = box_of(baseSolid->shape);
    INFO("base " << bobbinName << ": x [" << b.xlo * 1e3 << ", " << b.xhi * 1e3 << "] y [" << b.ylo * 1e3 << ", "
                 << b.yhi * 1e3 << "] z [" << b.zlo * 1e3 << ", " << b.zhi * 1e3 << "] mm; terminal plane "
                 << plane * 1e3 << " mm; ring underside " << -ring.at("C") / 2 * 1e3 << " mm");
    // Below the ring, top face on the terminal (seating) plane, the record's size, centred on the axis.
    CHECK(std::abs(b.yhi - plane) <= 1e-9);
    CHECK(std::abs((b.yhi - b.ylo) - height) <= 1e-9);
    CHECK(b.yhi < -ring.at("C") / 2);
    CHECK(std::abs((b.xhi - b.xlo) - width) <= 1e-9);
    CHECK(std::abs((b.zhi - b.zlo) - length) <= 1e-9);
    CHECK(std::abs(b.xhi + b.xlo) <= 1e-9);
    CHECK(std::abs(b.zhi + b.zlo) <= 1e-9);
    CHECK(std::abs(volume_of(baseSolid->shape) - length * width * height) <= 1e-15);

    // Pins: MKF's x, z and size; top on the base's bottom face.
    std::set<std::pair<long, long>> centres;
    for (const auto& [name, solid] : pinSolids) {
        const auto& mkf = pins.at(name);
        const auto c = mkf.get_coordinates().value();
        const auto& d = mkf.get_dimensions();
        const Box p = box_of(solid->shape);
        UNSCOPED_INFO("pin " << name << " x [" << p.xlo * 1e3 << ", " << p.xhi * 1e3 << "] y [" << p.ylo * 1e3 << ", "
                             << p.yhi * 1e3 << "] z [" << p.zlo * 1e3 << ", " << p.zhi * 1e3 << "] mm");
        CHECK(std::abs(0.5 * (p.xlo + p.xhi) - c[0]) <= 1e-9);
        CHECK(std::abs(0.5 * (p.zlo + p.zhi) - c[2]) <= 1e-9);
        CHECK(std::abs(p.yhi - b.ylo) <= 1e-9);
        CHECK(std::abs((p.yhi - p.ylo) - d[2]) <= 1e-9);
        centres.insert({std::lround(c[0] * 1e7), std::lround(c[2] * 1e7)});
    }
    CHECK(centres.size() == expectedPins);

    // Every drop still ends on the terminal plane, i.e. on the base's top face, pointing down.
    for (const auto& path : paths) {
        CHECK(std::abs(path.end0[1] - b.yhi) <= 1e-9);
        CHECK(std::abs(path.end1[1] - b.yhi) <= 1e-9);
        CHECK(path.dir0[1] < -(1.0 - 1e-9));
        CHECK(path.dir1[1] < -(1.0 - 1e-9));
    }

    // No shared volume between the base (or a pin) and the core or the copper. buildAllNamed's own
    // base gate already ran; this is the independent check.
    size_t cores = 0, turns = 0;
    for (const auto& ns : all) {
        if (ns.role != mvb::Role::Core && ns.role != mvb::Role::CoreCoating && ns.role != mvb::Role::Turn) continue;
        cores += ns.role == mvb::Role::Core;
        turns += ns.role == mvb::Role::Turn;
        std::vector<const mvb::NamedShape*> accessories = {baseSolid};
        for (const auto& [name, solid] : pinSolids) accessories.push_back(solid);
        for (const auto* accessory : accessories) {
            BRepAlgoAPI_Common common(accessory->shape, ns.shape);
            REQUIRE(common.IsDone());
            UNSCOPED_INFO(accessory->name << " vs " << ns.name);
            CHECK(volume_of(common.Shape()) < 1e-18);
        }
    }
    CHECK(cores > 0);
    CHECK(turns > 0);
}

} // namespace

TEST_CASE("buck_inductor_complete seated on HTM460-6: base under the ring, six tabs, nothing overlaps",
          "[toroidbase][abt1173]") {
    // RED until ABT #1253: on the base MKF winds this design in the coated ring window and the real-winding
    // build refuses its face chords (corner certifier). The test states the requirement; it is not relaxed.
    std::string reference;
    const json seated = seat_on_smallest_base(load_complete("buck_inductor_complete.json"), reference);
    CHECK(reference == "HTM460-6");
    require_seated(seated, 6, 6);
}

TEST_CASE("A CMC on T 16/9.6/6.3 seated on HTM600-6: base under the ring, six tabs, nothing overlaps",
          "[toroidbase][abt1173]") {
    std::string reference;
    const json seated = seat_on_smallest_base(cmc_on_t16(), reference);
    CHECK(reference == "HTM600-6");
    require_seated(seated, 6, 4);
}

TEST_CASE("A vertically mounted toroid base is refused, not drawn", "[toroidbase][abt1173]") {
    using Catch::Matchers::ContainsSubstring;
    // TEST FIXTURE (not catalogue data): the HTM600-6 record turned on edge, under the T 16/9.6/6.3 CMC.
    // No vertical base with a boatWidth exists in the catalogue, and the builder must not invent the boat.
    json magneticJson = cmc_on_t16();
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Bobbin record;
    for (const auto& candidate : OpenMagnetics::find_toroid_bases_for_core(core))
        if (candidate.get_manufacturer_info()->get_reference().value() == "HTM600-6") record = candidate;
    REQUIRE(record.get_name());
    json recordJson;
    to_json(recordJson, record);
    recordJson["name"] = "fixture: HTM600-6 on edge";
    recordJson["functionalDescription"]["base"]["mounting"] = "vertical";
    auto seated = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, OpenMagnetics::Bobbin(recordJson, false));
    json bobbinJson;
    to_json(bobbinJson, seated);
    magneticJson["coil"]["bobbin"] = bobbinJson;

    mvb::MagneticBuilder builder;
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    CHECK_THROWS_WITH(builder.buildAllNamed(enriched, /*includeBobbin=*/true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                            mvb::DEFAULT_CORE_POLYGON_SEGMENTS, false, false, false, 0.0,
                                            /*useRealWindingGeometry=*/true, /*femReady=*/true),
                      ContainsSubstring("not implemented for vertical bases without data"));
}
