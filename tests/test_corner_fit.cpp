// A leg MKF made exactly one tangent long holds its bend (ABT #2527).
//
// MKF plans a layer change (dragback) with radial legs exactly one bend long, on coordinates it ships rounded to 1 nm.
// The corner test compared the bend's tangent length with the leg's length to 8 ulps of the larger (3.3e-19 m on a
// 0.1835 mm leg), which is less than the round-off of a length taken between two such coordinates: on
// realwinding_interleaved_full_section_e16 four of the Secondary's corners missed it by 3.5e-19 and 7.0e-19 m and were
// drawn sharp (bisection mitre) where MKF had planned them round.
#include <catch2/catch_test_macros.hpp>
#include "mvb/TerminalFillet.h"
#include "mvb/WireAssembler.h"
#include <gp_Pnt.hxx>
#include <string>
#include <vector>

namespace {

mvb::Primitive dragbackStraight(const gp_Pnt& a, const gp_Pnt& b, int index) {
    mvb::Primitive p;
    p.kind = mvb::Primitive::SEG;
    p.seg = {a, b};
    p.isConnection = true;
    p.label = "wrap 'turn 17_ending' -> 'turn 18' (dragback) seg " + std::to_string(index);
    return p;
}

// An L: a leg `first` long along +z (the axis the fixture's coordinates are largest on), a right-angle corner at
// `corner`, a leg `second` long along +x.
std::vector<mvb::Primitive> dragbackL(const gp_Pnt& corner, double first, double second) {
    const gp_Pnt start(corner.X(), corner.Y(), corner.Z() - first);
    const gp_Pnt end(corner.X() + second, corner.Y(), corner.Z());
    return {dragbackStraight(start, corner, 0), dragbackStraight(corner, end, 1)};
}

// A Z: a long leg along +x, a step `step` long along +z, a long leg along +x. Both corners need the step.
std::vector<mvb::Primitive> dragbackZ(const gp_Pnt& corner, double step, double run) {
    const gp_Pnt start(corner.X() - run, corner.Y(), corner.Z());
    const gp_Pnt stepEnd(corner.X(), corner.Y(), corner.Z() + step);
    const gp_Pnt end(corner.X() + run, corner.Y(), corner.Z() + step);
    return {dragbackStraight(start, corner, 0), dragbackStraight(corner, stepEnd, 1), dragbackStraight(stepEnd, end, 2)};
}

size_t arcs(const std::vector<mvb::Primitive>& prims) {
    size_t n = 0;
    for (const auto& p : prims) n += p.kind == mvb::Primitive::ARC3;
    return n;
}

}  // namespace

TEST_CASE("A dragback leg exactly one tangent long holds MKF's bend", "[realwinding][cornerfit][abt2527]") {
    const double bend = 0.0001835, wireRadius = 0.00017;   // the fixture's Secondary: MKF's bend on a 0.34 mm copper
    const double run = 0.00474745203;
    // The corner walks the nanometre grid around where the fixture has it. A leg's length is taken between two
    // coordinates, so it is one bend (or two) long only to the round-off of that difference.
    size_t oneBendLegs = 0, twoBendSteps = 0;
    for (int i = 0; i < 200; ++i) {
        const gp_Pnt corner(0.000367, -0.00238404353 + i * 1e-9, -0.0135595 + i * 3e-9);
        INFO("corner at z = " << corner.Z() << " m");
        {
            auto prims = dragbackL(corner, bend, run);
            CHECK(mvb::filletDragbackInteriorCorners(prims, bend, wireRadius, "Secondary parallel 0") == 1);
            oneBendLegs += arcs(prims) == 1;
        }
        {
            auto prims = dragbackZ(corner, 2 * bend, run);
            CHECK(mvb::filletDragbackInteriorCorners(prims, bend, wireRadius, "Secondary parallel 0") == 2);
            twoBendSteps += arcs(prims) == 2;
        }
    }
    CHECK(oneBendLegs == 200);
    CHECK(twoBendSteps == 200);

    SECTION("a leg a nanometre short of the tangent does not") {
        auto prims = dragbackL(gp_Pnt(0.000367, -0.00238404353, -0.0135595), bend - 1e-9, run);
        CHECK(mvb::filletDragbackInteriorCorners(prims, bend, wireRadius, "Secondary parallel 0") == 0);
        CHECK(arcs(prims) == 0);
        CHECK(prims.size() == 2);
    }
}
