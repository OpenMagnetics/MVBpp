// Mirrors MVB Python's test_builder.py::test_all_shapes_generated.
// Iterates over MAS/data/core_shapes.ndjson and verifies every shape
// (excluding ui/pqi/ut, same exclusions as MVB Python) builds a non-empty
// solid via the C++ shape factory.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/catch_approx.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "mvb/shapes/ShapeBuilder.h"
#include "MAS.hpp"
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <gp_Pnt.hxx>
#include <fstream>
#include <map>
#include <numbers>
#include <optional>
#include <sstream>
#include <set>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

static const std::set<std::string> EXCLUDED_FAMILIES = {"ui", "pqi", "ut"};

#ifndef MAS_DATA_DIR
#define MAS_DATA_DIR "."
#endif

TEST_CASE("All MAS core shapes build a non-empty solid",
          "[shapes][all_shapes]") {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());

    int total = 0, skipped = 0, failed = 0;
    std::vector<std::string> failures;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json j = json::parse(line, nullptr, false);
        if (j.is_discarded()) continue;

        std::string family = j.value("family", "");
        std::string name = j.value("name", "?");
        if (EXCLUDED_FAMILIES.count(family)) { ++skipped; continue; }

        // Patch null dim subfields (same as the step generator does)
        mvb::patch_dimension_nominals(j);

        MAS::CoreShape shape;
        try {
            shape = j.get<MAS::CoreShape>();
        } catch (const std::exception& e) {
            failures.push_back(name + " (parse: " + e.what() + ")");
            ++failed; ++total;
            continue;
        }

        auto subtype = shape.get_family_subtype().value_or("");
        auto builder = mvb::shapes::createShapeBuilder(shape.get_family(), subtype);
        if (!builder) {
            failures.push_back(name + " (no builder)");
            ++failed; ++total;
            continue;
        }

        TopoDS_Shape piece;
        try {
            piece = builder->buildPiece(shape);
        } catch (const std::exception& e) {
            failures.push_back(name + " (build: " + e.what() + ")");
            ++failed; ++total;
            continue;
        }

        if (piece.IsNull()) {
            failures.push_back(name + " (null)");
            ++failed; ++total;
            continue;
        }

        GProp_GProps props;
        BRepGProp::VolumeProperties(piece, props);
        if (props.Mass() <= 0.0) {
            failures.push_back(name + " (zero volume)");
            ++failed;
        }
        ++total;
    }

    std::cerr << "[all_shapes] total=" << total << " skipped=" << skipped
              << " failed=" << failed << "\n";
    for (const auto& f : failures) std::cerr << "  FAIL: " << f << "\n";

    REQUIRE(total > 800);          // sanity: should hit ~882-952 shapes
    REQUIRE(failures.empty());
}

// ABT #1126. "Non-empty solid" above is a weak test: RM 7LP passed it for years while
// rendering as a mangled polyhedron, because its MAS record had no C and ShapeRM defaulted
// c to 0, extruding the profile with zero half-depth. Two things pin the fix.
//
// First, the data. TDK's RM 7 LP (B65819P) and Ferroxcube's RM7/ILP publish the same core
// to every digit -- l/A 0.52 mm^-1, le 23.5 mm, Ae 45.3 mm^2, Amin 39.6 mm^2, Ve 1060 mm^3 --
// and differ only in the height tolerance. So the two pieces must come out the same size.
// With C absent, RM 7LP's volume was a small fraction of its twin's.
//
// Second, the code. A shape whose profile needs a dimension it has not got must say so
// rather than build something from zero.
TEST_CASE("RM 7LP is the same core as its twin, and a missing dimension is refused",
          "[shapes][abt1126]") {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());

    auto volumeOf = [](json j) {
        mvb::patch_dimension_nominals(j);
        auto shape = j.get<MAS::CoreShape>();
        auto builder = mvb::shapes::createShapeBuilder(shape.get_family(),
                                                      shape.get_family_subtype().value_or(""));
        REQUIRE(builder != nullptr);
        auto piece = builder->buildPiece(shape);
        REQUIRE_FALSE(piece.IsNull());
        GProp_GProps props;
        BRepGProp::VolumeProperties(piece, props);
        return props.Mass();
    };

    json lowProfile, twin;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json j = json::parse(line, nullptr, false);
        if (j.is_discarded()) continue;
        if (j.value("name", "") == "RM 7LP") lowProfile = j;
        if (j.value("name", "") == "RM 7/ILP") twin = j;
    }
    REQUIRE_FALSE(lowProfile.is_null());
    REQUIRE_FALSE(twin.is_null());

    double lowProfileVolume = volumeOf(lowProfile);
    double twinVolume = volumeOf(twin);
    INFO("RM 7LP = " << lowProfileVolume * 1e9 << " mm3, RM 7/ILP = " << twinVolume * 1e9 << " mm3");
    // 5%: the two records differ only in D's upper tolerance, so the pieces are the same
    // size to well within this. The degenerate build was out by far more.
    CHECK(lowProfileVolume == Catch::Approx(twinVolume).epsilon(0.05));

    // Take C away again and the builder must refuse, naming what is missing.
    json crippled = lowProfile;
    crippled["dimensions"].erase("C");
    mvb::patch_dimension_nominals(crippled);
    auto shape = crippled.get<MAS::CoreShape>();
    auto builder = mvb::shapes::createShapeBuilder(shape.get_family(),
                                                  shape.get_family_subtype().value_or(""));
    REQUIRE(builder != nullptr);
    REQUIRE_THROWS_WITH(builder->buildPiece(shape),
                        Catch::Matchers::ContainsSubstring("RM 7LP")
                            && Catch::Matchers::ContainsSubstring("dimension C"));
}

// Mirrors MVB Python's test_get_families: every family in core_shapes.ndjson
// (except ui/pqi) must be routable via the factory.
TEST_CASE("All MAS families are supported by the factory",
          "[shapes][get_families]") {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());

    std::set<std::string> families;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json j = json::parse(line, nullptr, false);
        if (j.is_discarded()) continue;
        families.insert(j.value("family", ""));
    }

    std::vector<std::string> unsupported;
    for (const auto& fam : families) {
        if (fam == "ui" || fam == "pqi") continue;  // MVB Python excludes these
        // Look up via a minimal CoreShape just to test factory routing
        json stub = {{"family", fam}, {"name", "stub"},
                     {"type", "standard"}, {"dimensions", json::object()}};
        MAS::CoreShape s;
        try { s = stub.get<MAS::CoreShape>(); }
        catch (...) { unsupported.push_back(fam + " (parse failed)"); continue; }
        auto b = mvb::shapes::createShapeBuilder(s.get_family(), "");
        if (!b) unsupported.push_back(fam);
    }

    std::cerr << "[get_families] families=" << families.size()
              << " unsupported=" << unsupported.size() << "\n";
    for (const auto& u : unsupported) std::cerr << "  UNSUPPORTED: " << u << "\n";

    REQUIRE(unsupported.empty());
}

// ABT #1143. The RM fix above was one builder; every other one still read an absent dimension
// as zero and built something from it: a PQ whose lateral legs took an angle made up from E and
// F when G was missing, an E core with no window when D was, a shielded drum drawn as a bare drum
// when its ring letters were. A piece built from a missing letter is a different part, drawn and
// meshed as if it were the real one. This table is the contract: per family, the letters its
// piece cannot be built without. Take any one of them away from any catalogue shape of the family
// and the builder must refuse, naming the shape and the letter.
//
// Letters that are NOT here are the features only some records have, read as "no such feature"
// when absent, the way MKF reads the same records: F2 (the centre column's own depth; E, EL),
// G (the slot of an ER / EQ / EP / EPX), K of an EP (MKF: F / 2), J and L of a PQ (MKF: F / 2 and
// F + (C - F) / 3, both or neither), H (a centre hole or bore: RM, PM, drum), A2 (a drum's second
// flange), and the UR letters its subtype does not use.
TEST_CASE("No shape builder reads an absent dimension as zero", "[shapes][abt1143]") {
    const std::map<std::string, std::vector<std::string>> requiredOf = {
        {"e", {"A", "B", "C", "D", "E", "F"}},
        {"ei", {"A", "B", "C", "D", "E", "F"}},
        {"planarE", {"A", "B", "C", "D", "E", "F"}},
        {"planarEL", {"A", "B", "C", "D", "E", "F"}},
        {"efd", {"A", "B", "C", "D", "E", "F", "F2", "K", "q"}},
        {"er", {"A", "B", "C", "D", "E", "F"}},
        {"etd", {"A", "B", "C", "D", "E", "F"}},
        {"eq", {"A", "B", "C", "D", "E", "F"}},
        {"eer", {"A", "B", "C", "D", "E", "F"}},
        {"planarER", {"A", "B", "C", "D", "E", "F"}},
        {"ec", {"A", "B", "C", "D", "E", "F", "T", "s"}},
        {"ep", {"A", "B", "C", "D", "E", "F"}},
        {"epx", {"A", "B", "C", "D", "E", "F", "K"}},
        {"epc", {"A", "B", "C", "D", "E", "F", "F2"}},
        {"lp", {"A", "B", "C", "D", "E", "F"}},
        {"p", {"A", "B", "D", "E", "F"}},
        {"pm", {"A", "B", "C", "D", "E", "F", "alpha"}},
        {"pq", {"A", "B", "C", "D", "E", "F", "G"}},
        {"pqi", {"A", "B", "C", "D", "E", "F", "G"}},
        {"rm", {"A", "B", "C", "D", "E", "F", "G", "J"}},
        {"t", {"A", "B", "C"}},
        {"u", {"A", "B", "C", "D", "E"}},
        {"c", {"A", "B", "C", "D", "E"}},
        {"ur", {"A", "B", "C", "D"}},
        {"ut", {"A", "B", "C", "D", "E", "F"}},
        {"drum", {"A", "B", "C", "D", "E", "F"}},
        {"drumRing", {"A", "B", "C", "D", "E", "F", "J", "K", "L"}},
        {"drumSemishielded", {"A", "B", "C", "D", "E", "F"}},
    };

    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());

    size_t shapes = 0, refusals = 0;
    std::set<std::string> familiesSeen;
    std::vector<std::string> failures;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json record = json::parse(line, nullptr, false);
        if (record.is_discarded()) continue;
        const std::string family = record.value("family", "");
        const std::string name = record.value("name", "?");
        const auto required = requiredOf.find(family);
        if (required == requiredOf.end()) continue;
        familiesSeen.insert(family);
        ++shapes;
        for (const auto& letter : required->second) {
            if (!record["dimensions"].contains(letter)) {
                failures.push_back(name + ": the catalogue record itself has no " + letter);
                continue;
            }
            json crippled = record;
            crippled["dimensions"].erase(letter);
            mvb::patch_dimension_nominals(crippled);
            const auto shape = crippled.get<MAS::CoreShape>();
            auto builder = mvb::shapes::createShapeBuilder(shape.get_family(), shape.get_family_subtype().value_or(""));
            if (!builder) {
                failures.push_back(name + ": no builder");
                break;
            }
            try {
                builder->buildPiece(shape);
                failures.push_back(name + " without " + letter + ": built anyway");
            } catch (const std::exception& e) {
                const std::string what = e.what();
                if (what.find(name) == std::string::npos || what.find("dimension " + letter) == std::string::npos) {
                    failures.push_back(name + " without " + letter + ": refused, but the message names neither: " + what);
                }
                ++refusals;
            }
        }
    }
    std::cerr << "[abt1143] " << shapes << " shapes of " << familiesSeen.size() << " families, " << refusals
              << " refusals, " << failures.size() << " failures\n";
    for (size_t i = 0; i < failures.size() && i < 40; ++i) std::cerr << "  FAIL: " << failures[i] << "\n";
    CHECK(failures.empty());
    // Every family of the table is in the catalogue, and the walk was not empty.
    CHECK(familiesSeen.size() == requiredOf.size());
    CHECK(shapes > 2000);
    CHECK(refusals > 8000);
}

namespace {

json catalogue_shape(const std::string& wanted) {
    std::ifstream f(std::string(MAS_DATA_DIR) + "/core_shapes.ndjson");
    REQUIRE(f.is_open());
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json record = json::parse(line, nullptr, false);
        if (!record.is_discarded() && record.value("name", "") == wanted) return record;
    }
    FAIL("no catalogue shape named " << wanted);
    return json();
}

TopoDS_Shape piece_of(json record, int corePolygonSegments = 0) {
    mvb::patch_dimension_nominals(record);
    const auto shape = record.get<MAS::CoreShape>();
    auto builder = mvb::shapes::createShapeBuilder(shape.get_family(), shape.get_family_subtype().value_or(""), corePolygonSegments);
    REQUIRE(builder != nullptr);
    return builder->buildPiece(shape);
}

// The record's dimensions as the builders read them.
std::map<std::string, double> dimensions_of(json record) {
    mvb::patch_dimension_nominals(record);
    const auto dimensions = record.get<MAS::CoreShape>().get_dimensions();
    REQUIRE(dimensions.has_value());
    return mvb::flatten_dimensions(*dimensions);
}

double volume_of(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}

}  // namespace

// The features that only some records carry are features: taking one away removes that feature
// and nothing else, and where MKF fills an absent letter by rule the piece is drawn with MKF's value.
TEST_CASE("A feature dimension is a feature, and MKF's rule fills what MKF fills", "[shapes][abt1143]") {
    SECTION("an RM's centre hole H") {
        const json withHole = catalogue_shape("RM 4");
        REQUIRE(withHole["dimensions"].contains("H"));
        json solid = withHole;
        solid["dimensions"].erase("H");
        const double b = dimensions_of(withHole).at("B");
        const double h = dimensions_of(withHole).at("H");
        const double hole = std::numbers::pi * h * h / 4.0 * b;
        CHECK(volume_of(piece_of(solid)) - volume_of(piece_of(withHole)) == Catch::Approx(hole).epsilon(1e-6));
    }
    SECTION("a PQ's lateral flank J and L") {
        // PQ 32/20 carries neither. MKF then takes J = F / 2 and L = F + (C - F) / 3; the same
        // record with those two written out is the same piece.
        const json plain = catalogue_shape("PQ 32/20");
        REQUIRE_FALSE(plain["dimensions"].contains("J"));
        REQUIRE_FALSE(plain["dimensions"].contains("L"));
        const auto dims = dimensions_of(plain);
        json written = plain;
        written["dimensions"]["J"] = {{"nominal", dims.at("F") / 2.0}};
        written["dimensions"]["L"] = {{"nominal", dims.at("F") + (dims.at("C") - dims.at("F")) / 3.0}};
        CHECK(volume_of(piece_of(written)) == Catch::Approx(volume_of(piece_of(plain))).epsilon(1e-9));
        // One of the two without the other is no flank at all.
        json half = plain;
        half["dimensions"]["J"] = {{"nominal", dims.at("F") / 2.0}};
        CHECK_THROWS_WITH(piece_of(half), Catch::Matchers::ContainsSubstring("PQ 32/20") && Catch::Matchers::ContainsSubstring("J and L"));
    }
}

// ABT #2586. EPC went through the EP pot builder, whose K (the post's place in the depth of the
// pot) no EPC record carries; read as zero it put the whole A x C profile on one side of the
// column axis and cut the EP's round window on its edge: EPC 25 came out with a half-disc pole of
// radius F / 2 at z = -5.75 .. 0 mm and a core spanning z = -8 .. 0. MKF models, and winds around,
// a stadium pole F x F2 centred on the axis of an E-like core (CorePieceEpc), so the bobbin and
// half of every turn lay outside the drawn core.
TEST_CASE("An EPC core is an E with a stadium pole centred on its axis", "[shapes][abt2586]") {
    size_t checked = 0;
    for (const std::string name : {"EPC 13", "EPC 17", "EPC 25", "EPC 30"}) {
        INFO(name);
        const json record = catalogue_shape(name);
        const auto dims = dimensions_of(record);
        const double a = dims.at("A"), b = dims.at("B"), c = dims.at("C"), d = dims.at("D"), e = dims.at("E"), f = dims.at("F"),
                     f2 = dims.at("F2");
        const TopoDS_Shape piece = piece_of(record);   // exact curves

        // Where it stands: A wide about the axis, B high below the mating plane y = 0, C deep about the axis.
        Bnd_Box box;
        BRepBndLib::AddOptimal(piece, box, /*useTriangulation=*/false, /*useShapeTolerance=*/false);
        double x0, y0, z0, x1, y1, z1;
        box.Get(x0, y0, z0, x1, y1, z1);
        const double micron = 1e-6;
        CHECK(x0 == Catch::Approx(-a / 2).margin(micron));
        CHECK(x1 == Catch::Approx(a / 2).margin(micron));
        CHECK(y0 == Catch::Approx(-b).margin(micron));
        CHECK(y1 == Catch::Approx(0.0).margin(micron));
        CHECK(z0 == Catch::Approx(-c / 2).margin(micron));
        CHECK(z1 == Catch::Approx(c / 2).margin(micron));

        // What a cut through the window shows: the two outer legs, (A - E) / 2 wide and C deep, and
        // the pole, a stadium F x F2 (MKF's column area: F F2 - (4 - pi) / 4 F2^2).
        const double slab = d / 10;
        const TopoDS_Shape cutter = BRepPrimAPI_MakeBox(gp_Pnt(-a, -d / 2 - slab / 2, -c), 2 * a, slab, 2 * c).Shape();
        BRepAlgoAPI_Common section(piece, cutter);
        REQUIRE(section.IsDone());
        const double stadium = f * f2 - (4 - std::numbers::pi) / 4 * f2 * f2;
        INFO("section " << volume_of(section.Shape()) / slab * 1e6 << " mm2, legs " << (a - e) * c * 1e6 << " mm2, stadium "
                        << stadium * 1e6 << " mm2");
        CHECK(volume_of(section.Shape()) / slab == Catch::Approx((a - e) * c + stadium).epsilon(1e-6));

        // The pole alone, about the axis: nothing of it beyond F / 2 along the width or F2 / 2 along the depth.
        const TopoDS_Shape poleRoom = BRepPrimAPI_MakeBox(gp_Pnt(-f / 2, -d / 2 - slab / 2, -f2 / 2), f, slab, f2).Shape();
        BRepAlgoAPI_Common pole(piece, poleRoom);
        REQUIRE(pole.IsDone());
        CHECK(volume_of(pole.Shape()) / slab == Catch::Approx(stadium).epsilon(1e-6));
        ++checked;
    }
    CHECK(checked == 4);
}
