// Slab cores DS / HS / RS (MKF CorePieceSlab, ABT #263, MKF f31bca46).
//
// MKF learned the slab families and the adviser started choosing them (FSBB's default design
// advises a gapped DS 14/08), but MVB++'s shape factory had no builder for them. buildCoreShapes
// skipped the unknown family, every piece of the core was dropped, and the WASM drawCore
// reported the empty list as "side='' filtered out all geometry".
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "mvb/shapes/ShapeBuilder.h"
#include "MAS.hpp"
#include "constructive_models/Core.h"
#include <nlohmann/json.hpp>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <cmath>
#include <fstream>
#include <numbers>
#include <string>

using json = nlohmann::json;

namespace {

json slab_shape(const std::string& name) {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json j = json::parse(line, nullptr, false);
        if (!j.is_discarded() && j.value("name", "") == name) {
            mvb::patch_dimension_nominals(j);
            return j;
        }
    }
    FAIL("shape " << name << " is not in core_shapes.ndjson");
    return {};
}

double nominal(const json& shape, const char* letter) {
    const auto& d = shape.at("dimensions").at(letter);
    return d.is_number() ? d.get<double>() : d.at("nominal").get<double>();
}

// Area of the disc of radius r that survives the two flats at |y| = c.
double clipped_disc(double r, double c) {
    if (r <= c) return std::numbers::pi * r * r;
    return std::numbers::pi * r * r - 2.0 * (r * r * std::acos(c / r) - c * std::sqrt(r * r - c * c));
}

OpenMagnetics::Core gapped_core(const json& shape) {
    json coreJson;
    coreJson["functionalDescription"] = {
        {"name", "slab"}, {"type", "two-piece set"}, {"material", "N97"}, {"shape", shape},
        {"gapping", json::array({ json{{"length", 0.0001}, {"type", "subtractive"}} })},
        {"numberStacks", 1},
    };
    return OpenMagnetics::Core(coreJson);
}

} // namespace

TEST_CASE("DS and HS pieces are pot cores cut to the flats", "[shapes][slab]") {
    for (const std::string name : {"DS 14/08", "HS 14/08", "DS 26/16", "HS 26/16"}) {
        INFO(name);
        json j = slab_shape(name);
        auto shape = j.get<MAS::CoreShape>();
        // Exact circles, so the closed form below is the whole check.
        auto builder = mvb::shapes::createShapeBuilder(shape.get_family(), "", 0);
        REQUIRE(builder != nullptr);
        auto piece = builder->buildPiece(shape);
        REQUIRE_FALSE(piece.IsNull());

        const double a = nominal(j, "A"), b = nominal(j, "B"), c = nominal(j, "C");
        const double d = nominal(j, "D"), e = nominal(j, "E"), f = nominal(j, "F");
        const double h = j["dimensions"].contains("H") ? nominal(j, "H") : 0.0;

        // Envelope: width A, height B, depth C (MKF CorePieceSlab: set_width(A), set_depth(C)).
        Bnd_Box box;
        BRepBndLib::AddOptimal(piece, box, false, false);
        double x0, y0, z0, x1, y1, z1;
        box.Get(x0, y0, z0, x1, y1, z1);
        CHECK((x1 - x0) == Catch::Approx(a).epsilon(0.01));
        CHECK((y1 - y0) == Catch::Approx(b).epsilon(0.01));
        CHECK((z1 - z0) == Catch::Approx(c).epsilon(0.01));

        // Volume against the closed form.
        const double hc = c / 2.0;
        const double expected = b * clipped_disc(a / 2.0, hc)
                              - d * (clipped_disc(e / 2.0, hc) - std::numbers::pi * f * f / 4.0)
                              - b * std::numbers::pi * h * h / 4.0;
        GProp_GProps props;
        BRepGProp::VolumeProperties(piece, props);
        CHECK(props.Mass() == Catch::Approx(expected).epsilon(0.002));
    }
}

TEST_CASE("A gapped DS core as MKF lays it out draws both halves", "[shapes][slab]") {
    auto core = gapped_core(slab_shape("DS 14/08"));
    mvb::MagneticBuilder builder;
    auto pieces = builder.buildCoreNamed(core);
    REQUIRE(pieces.size() == 2);
    for (const auto& p : pieces) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(p.shape, props);
        CHECK(props.Mass() > 0.0);
    }
}

TEST_CASE("A core family with no 3D builder is refused by name, not drawn empty", "[shapes][slab]") {
    // RS pairs a slab half with a plain pot round, and MKF's geometricalDescription carries two
    // identical RS halfSets: nothing says which half is round, so there is no builder for it.
    auto core = gapped_core(slab_shape("RS 14/08"));
    mvb::MagneticBuilder builder;
    REQUIRE_THROWS_WITH(builder.buildCoreNamed(core),
                        Catch::Matchers::ContainsSubstring("no 3D geometry builder")
                            && Catch::Matchers::ContainsSubstring("rs"));
}
