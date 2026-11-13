// Piece-and-plate cores draw their plate (ABT #2372, web report #193).
//
// An EI is an E piece closed by an I bar: MKF emits it as a HALF_SET piece plus a PLATE element.
// buildCoreShapes skipped every PLATE, so the 3D showed a bare E, and the front view labelled only
// B and D: the I bar's thickness B2, one of the seven dimensions the user edits, was on no drawing.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/SectionDrawing.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include <nlohmann/json.hpp>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <string>

using json = nlohmann::json;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// MAS catalogue EI 12.7/3.4/1.25/12.7/2.77/7.11/3.56: A 12.7, B 3.4, B2 1.25, C 12.7, D 2.77 mm.
constexpr double kA = 12.7e-3, kB = 3.4e-3, kB2 = 1.25e-3, kC = 12.7e-3;

OpenMagnetics::Magnetic ei_magnetic() {
    json magnetic = json::parse(R"({
        "core": {"name": "EI plate test", "functionalDescription": {
            "type": "pieceAndPlate", "material": "3C95", "numberStacks": 1,
            "shape": "EI 12.7/3.4/1.25/12.7/2.77/7.11/3.56",
            "gapping": [{"type": "subtractive", "length": 0.00005}]}},
        "coil": {"bobbin": "Basic", "functionalDescription": [
            {"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary",
             "wire": "Round 0.2 - Grade 1"}]}
    })");
    return mvb::magnetic_autocomplete_safe(magnetic);
}

struct Box { double x0, y0, z0, x1, y1, z1; };

Box box_of(const TopoDS_Shape& shape) {
    Bnd_Box b;
    BRepBndLib::AddOptimal(shape, b, false, false);
    Box r{};
    b.Get(r.x0, r.y0, r.z0, r.x1, r.y1, r.z1);
    return r;
}

} // namespace

TEST_CASE("An EI core draws its I bar as an A x B2 x C plate on the legs of the E", "[shapes][piece-and-plate][abt-2372]") {
    auto magnetic = ei_magnetic();
    REQUIRE(magnetic.get_core().get_functional_description().get_type() == MAS::CoreType::PIECE_AND_PLATE);

    mvb::MagneticBuilder builder;
    auto pieces = builder.buildCoreNamed(magnetic.get_core());
    REQUIRE(pieces.size() == 2);

    // The plate is the piece above the mating plane; the E hangs below it.
    const Box first = box_of(pieces[0].shape), second = box_of(pieces[1].shape);
    const bool firstIsPlate = first.y0 > second.y0;
    const auto& plateShape = firstIsPlate ? pieces[0].shape : pieces[1].shape;
    const Box plate = firstIsPlate ? first : second;
    const Box piece = firstIsPlate ? second : first;

    GProp_GProps props;
    BRepGProp::VolumeProperties(plateShape, props);
    CHECK_THAT(props.Mass(), WithinRel(kA * kB2 * kC, 1e-9));

    const double tol = 1e-7;
    CHECK_THAT(plate.x0, WithinAbs(-kA / 2, tol));
    CHECK_THAT(plate.x1, WithinAbs(kA / 2, tol));
    CHECK_THAT(plate.z0, WithinAbs(-kC / 2, tol));
    CHECK_THAT(plate.z1, WithinAbs(kC / 2, tol));
    CHECK_THAT(plate.y1 - plate.y0, WithinAbs(kB2, tol));
    // It closes the E: it starts where the E's legs end, and the E is B tall below it.
    CHECK_THAT(plate.y0, WithinAbs(piece.y1, tol));
    CHECK_THAT(piece.y1 - piece.y0, WithinAbs(kB, tol));
}

TEST_CASE("The front view of an EI dimensions the I bar's thickness B2", "[2d][frontview][piece-and-plate][abt-2372]") {
    auto magnetic = ei_magnetic();
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedFrontView(magnetic, 800, 12));
    CHECK(svg.find("B2 = 1.25") != std::string::npos);
    CHECK(svg.find("B = 3.40") != std::string::npos);
    CHECK(svg.find("D = 2.77") != std::string::npos);
}
