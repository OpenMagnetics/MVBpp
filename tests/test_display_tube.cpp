// The display tube (ABT #2489): a conductor's display drawing triangulated straight from its path, no CAD solids.
// What is asserted is the drawn surface itself: that it is closed, where it lies, and how much it encloses.
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "mvb/DisplayTube.h"
#include "mvb/MagneticBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"

#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <nlohmann/json.hpp>
#include <numbers>

using json = nlohmann::json;

namespace {

// The triangles a shape carries, read as a viewer reads them.
struct Drawn {
    size_t faces = 0, triangles = 0, degenerate = 0;
    size_t openEdges = 0;        // edges only one triangle uses: a hole in the skin
    size_t crowdedEdges = 0;     // edges more than two triangles use
    double volume = 0.0;         // enclosed, signed by the triangles' own orientation [m^3]
    std::array<double, 3> lo{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    std::array<double, 3> hi{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
};

Drawn drawnOf(const TopoDS_Shape& shape) {
    Drawn d;
    std::map<std::array<long long, 3>, int> vertexOf;   // welded on a nanometre grid
    std::map<std::pair<int, int>, int> uses;
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
        ++d.faces;
        TopLoc_Location location;
        const Handle(Poly_Triangulation) mesh = BRep_Tool::Triangulation(TopoDS::Face(ex.Current()), location);
        REQUIRE(!mesh.IsNull());
        for (int t = 1; t <= mesh->NbTriangles(); ++t) {
            int n[3];
            mesh->Triangle(t).Get(n[0], n[1], n[2]);
            gp_XYZ p[3];
            int id[3];
            for (int k = 0; k < 3; ++k) {
                p[k] = mesh->Node(n[k]).Transformed(location.Transformation()).XYZ();
                const std::array<long long, 3> key{std::llround(p[k].X() * 1e9), std::llround(p[k].Y() * 1e9), std::llround(p[k].Z() * 1e9)};
                id[k] = vertexOf.emplace(key, static_cast<int>(vertexOf.size())).first->second;
                for (int a = 0; a < 3; ++a) {
                    d.lo[a] = std::min(d.lo[a], p[k].Coord(a + 1));
                    d.hi[a] = std::max(d.hi[a], p[k].Coord(a + 1));
                }
            }
            ++d.triangles;
            d.volume += p[0].Dot(p[1].Crossed(p[2])) / 6.0;
            if (id[0] == id[1] || id[1] == id[2] || id[0] == id[2]) {
                ++d.degenerate;
                continue;
            }
            for (int k = 0; k < 3; ++k) ++uses[std::minmax(id[k], id[(k + 1) % 3])];
        }
    }
    for (const auto& [edge, count] : uses) {
        if (count == 1) ++d.openEdges;
        if (count > 2) ++d.crowdedEdges;
    }
    return d;
}

mvb::TubePiece straight(const gp_Pnt& a, const gp_Pnt& b) {
    mvb::TubePiece piece;
    piece.prim.kind = mvb::Primitive::SEG;
    piece.prim.seg.a = a;
    piece.prim.seg.b = b;
    return piece;
}

double areaOf(const std::vector<gp_XY>& outline) {
    double twice = 0.0;
    for (size_t k = 0; k < outline.size(); ++k) twice += outline[k].Crossed(outline[(k + 1) % outline.size()]);
    return 0.5 * twice;
}

} // namespace

TEST_CASE("A display tube along one straight is the closed prism of the wire's outline", "[displaytube][abt2489]") {
    const double radius = 0.5e-3, length = 10e-3;
    const auto outline = mvb::roundTubeOutline(radius, mvb::kDrawingAngularStep);
    REQUIRE(outline.size() == 32);   // 2 pi / 0.2, rounded up
    const Drawn d = drawnOf(mvb::displayTube({straight(gp_Pnt(0, 0, 0), gp_Pnt(length, 0, 0))}, outline, mvb::kDrawingAngularStep));
    CHECK(d.faces == 1);
    CHECK(d.triangles == 2 * 32 + 2 * 32);   // the skin of one chord, and the two caps
    CHECK(d.openEdges == 0);
    CHECK(d.crowdedEdges == 0);
    CHECK(std::abs(d.volume - areaOf(outline) * length) < 1e-9 * areaOf(outline) * length);   // outward, and all of it
    CHECK(std::abs(d.lo[0]) < 1e-12);
    CHECK(std::abs(d.hi[0] - length) < 1e-12);
    for (int a = 1; a < 3; ++a) {
        CHECK(std::abs(d.hi[a] - radius) < 1e-12);
        CHECK(d.lo[a] > -radius - 1e-12);
    }
}

TEST_CASE("Where its path turns a corner a display tube's two pieces share one section on the bisecting plane",
          "[displaytube][abt2489]") {
    const double radius = 0.4e-3, leg = 6e-3;
    const auto outline = mvb::roundTubeOutline(radius, mvb::kDrawingAngularStep);
    const Drawn d = drawnOf(mvb::displayTube({straight(gp_Pnt(0, 0, 0), gp_Pnt(leg, 0, 0)), straight(gp_Pnt(leg, 0, 0), gp_Pnt(leg, 0, leg))},
                                             outline, mvb::kDrawingAngularStep));
    CHECK(d.faces == 2);
    CHECK(d.openEdges == 0);      // the two pieces meet on the very same points: one skin
    CHECK(d.crowdedEdges == 0);
    // A prism cut by a plane through its centreline holds its section times its centreline length: the mitre neither
    // adds copper nor leaves a wedge open.
    CHECK(std::abs(d.volume - areaOf(outline) * 2.0 * leg) < 1e-9 * areaOf(outline) * 2.0 * leg);
    // The outer corner of the mitre reaches one radius past each leg's end.
    CHECK(std::abs(d.hi[0] - (leg + radius)) < 1e-9);
    CHECK(std::abs(d.lo[2] + radius) < 1e-9);
}

TEST_CASE("A display tube follows a bend as the torus the wire makes of it", "[displaytube][abt2489]") {
    const double radius = 0.5e-3, bend = 5e-3;
    mvb::TubePiece quarter;
    quarter.prim.kind = mvb::Primitive::ARC3;
    quarter.prim.arc.c = gp_Pnt(0, 0, 0);
    quarter.prim.arc.axis = gp_XYZ(0, 1, 0);
    quarter.prim.arc.v0 = gp_XYZ(bend, 0, 0);
    quarter.prim.arc.sweep = 0.5 * std::numbers::pi;
    const auto outline = mvb::roundTubeOutline(radius, mvb::kDrawingAngularStep);
    const Drawn d = drawnOf(mvb::displayTube({quarter}, outline, mvb::kDrawingAngularStep));
    CHECK(d.openEdges == 0);
    CHECK(d.crowdedEdges == 0);
    CHECK(d.triangles == 2 * 32 * 8 + 2 * 32);   // a quarter turn in chords of 0.2 rad is 8 chords
    const double torus = std::numbers::pi * radius * radius * bend * 0.5 * std::numbers::pi;
    INFO("drawn " << d.volume * 1e9 << " mm^3 of the bend's " << torus * 1e9);
    CHECK(d.volume < torus);              // chords lie inside the curve they follow
    CHECK(d.volume > 0.985 * torus);      // 32 sides and 0.2 rad chords are within 1.5 % of it
    CHECK(std::abs(d.hi[0] - (bend + radius)) < 1e-9);   // it starts on +X and turns about +Y towards -Z
    CHECK(std::abs(d.lo[2] + (bend + radius)) < 1e-9);
}

TEST_CASE("A display tube lays a rectangular section with its height on the axial direction", "[displaytube][abt2489]") {
    const double width = 3e-3, height = 0.5e-3, length = 8e-3;
    mvb::TubePiece piece = straight(gp_Pnt(0, 0, 0), gp_Pnt(length, 0, 0));
    piece.oriented = true;   // axial: the column axis Y, the default
    SECTION("sharp corners") {
        const auto outline = mvb::rectTubeOutline(width, height, 0.0, mvb::kDrawingAngularStep);
        REQUIRE(outline.size() == 4);
        const Drawn d = drawnOf(mvb::displayTube({piece}, outline, mvb::kDrawingAngularStep));
        CHECK(d.openEdges == 0);
        CHECK(std::abs(d.volume - width * height * length) < 1e-9 * width * height * length);
        CHECK(std::abs(d.hi[1] - 0.5 * height) < 1e-12);   // height along Y
        CHECK(std::abs(d.lo[1] + 0.5 * height) < 1e-12);
        CHECK(std::abs(d.hi[2] - 0.5 * width) < 1e-12);    // width across it, as rectProfileWire lays it
        CHECK(std::abs(d.lo[2] + 0.5 * width) < 1e-12);
    }
    SECTION("rounded corners") {
        const double corner = 0.2e-3;
        const auto outline = mvb::rectTubeOutline(width, height, corner, mvb::kDrawingAngularStep);
        const Drawn d = drawnOf(mvb::displayTube({piece}, outline, mvb::kDrawingAngularStep));
        CHECK(d.openEdges == 0);
        const double exact = (width * height - (4.0 - std::numbers::pi) * corner * corner) * length;
        INFO("drawn " << d.volume * 1e9 << " mm^3 of " << exact * 1e9);
        CHECK(d.volume < exact);
        CHECK(d.volume > 0.999 * exact);
        CHECK(std::abs(d.hi[1] - 0.5 * height) < 1e-12);
        CHECK(std::abs(d.hi[2] - 0.5 * width) < 1e-12);
    }
    SECTION("a corner radius that does not fit the section is refused") {
        CHECK_THROWS(mvb::rectTubeOutline(width, height, 0.3e-3, mvb::kDrawingAngularStep));
    }
}

// The whole path: MagneticBuilder::setDisplayTriangles, as the web's 3D view asks for it. Each conductor is then one
// closed skin of triangles, no solid, that lies where the solids of the same drawing lie and encloses their copper.
TEST_CASE("A display drawing asked for as triangles draws each conductor where its solids are", "[displaytube][realwinding][abt2489]") {
    // `within`: how far the triangles' bounds may lie from the solids' [m]. Five microns where both draw the same
    // thing. The flat wire swept whole is the exception, and there it is the SOLIDS that are off the wire: its leads
    // are built 10 % taller than the wire (5 % proud each side, 30 um here, so that the lead meets the wrap
    // transversally for the mesher) and overlap the wrap instead of meeting it; the tube draws the wire's own section
    // and mitres that corner, whose outer tip follows the sloping stub 50 um further.
    const auto [file, rectangular, within] = GENERATE(table<std::string, bool, double>({
        {std::string(MAS_COMPLETE_DIR) + "/single_switch_forward_transformer_complete.json", false, 5e-6},   // round wire, round column
        {std::string(MAS_EXAMPLES_DIR) + "/03_buck_inductor_pq3230_n95.json", true, 60e-6},                   // flat wire, swept whole on a round column
        {std::string(MAS_EXAMPLES_DIR) + "/18_stacked_inductor_e7033_n27.json", true, 5e-6},                  // flat wire, piece by piece on a rectangular column
    }));
    INFO(file);
    std::ifstream f(file);
    REQUIRE(f.is_open());
    json j;
    f >> j;
    const auto enriched = mvb::magnetic_autocomplete_safe(j.at("magnetic"), /*useRealWindingGeometry=*/true);
    const auto draw = [&](bool triangles) {
        mvb::MagneticBuilder builder;
        builder.setDisplayTriangles(triangles);
        return builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0, /*wirePolygonSegments=*/0,
                                     mvb::DEFAULT_CORE_POLYGON_SEGMENTS, /*paintCoating=*/true, /*emitCoatingShells=*/false,
                                     /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                     /*useRealWindingGeometry=*/true, /*femReady=*/false);
    };
    const auto solids = draw(false), triangles = draw(true);
    REQUIRE(solids.size() == triangles.size());
    int conductors = 0;
    for (size_t k = 0; k < solids.size(); ++k) {
        REQUIRE(solids[k].name == triangles[k].name);
        if (solids[k].role != mvb::Role::Turn) continue;
        ++conductors;
        INFO(solids[k].name);
        CHECK_FALSE(TopExp_Explorer(triangles[k].shape, TopAbs_SOLID).More());
        const Drawn d = drawnOf(triangles[k].shape);
        CHECK(d.triangles > 100);
        CHECK(d.openEdges == 0);
        CHECK(d.crowdedEdges == 0);
        GProp_GProps exact;
        BRepGProp::VolumeProperties(solids[k].shape, exact);
        INFO("drawn " << d.volume * 1e9 << " mm^3, the solids hold " << exact.Mass() * 1e9);
        CHECK(std::abs(d.volume - exact.Mass()) < 0.025 * exact.Mass());
        Bnd_Box box;
        BRepBndLib::AddOptimal(solids[k].shape, box, /*useTriangulation=*/false);
        double lo[3], hi[3];
        box.Get(lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
        for (int a = 0; a < 3; ++a) {
            INFO("axis " << a << ": triangles " << d.lo[a] * 1e3 << " .. " << d.hi[a] * 1e3 << " mm, solids " << lo[a] * 1e3 << " .. " << hi[a] * 1e3);
            CHECK(std::abs(d.lo[a] - lo[a]) < within);
            CHECK(std::abs(d.hi[a] - hi[a]) < within);
        }
    }
    CHECK(conductors >= 1);
    CHECK((rectangular || conductors > 1));
}
