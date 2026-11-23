// Multi-column winding placement (MAS windingWindow / winding window `column` edge):
// turns whose section resolves to a LATERAL core column must be built AROUND that
// column (loop centred on the leg axis, sized to the leg + bobbin wall), while
// main-column turns keep the legacy origin-centred loops byte-identical.
//
// Fixture: an E42/21/20 transformer wound by MKF with per-column winding windows —
// 24 primary turns on the central column, 12 secondary turns on the lateral column
// at x = -18.0625 mm (see MKF tests/TestMultiColumnWindingWindows.cpp for the
// producing side). Turn coordinates are ABSOLUTE MAS data; the secondary crossings
// sit at x = -13.6 mm (window side) and -22.5 mm (outer side) of the leg.
#include <catch2/catch_test_macros.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

#ifndef MAS_COMPLETE_DIR
#define MAS_COMPLETE_DIR "."
#endif

namespace {
struct BBox { double xmin, ymin, zmin, xmax, ymax, zmax; };
BBox bbox_of(const TopoDS_Shape& s) {
    Bnd_Box b;
    BRepBndLib::Add(s, b);
    BBox r{};
    b.Get(r.xmin, r.ymin, r.zmin, r.xmax, r.ymax, r.zmax);
    return r;
}
}  // namespace

TEST_CASE("Lateral-column turns wrap their wound column, main-column turns stay centred",
          "[turns][multicolumn]") {
    std::ifstream f(std::string(MAS_COMPLETE_DIR) + "/multicolumn_e42_transformer.json");
    REQUIRE(f.is_open());
    json j;
    f >> j;
    MAS::Magnetic magnetic = mvb::magnetic_from_json(j);
    REQUIRE(magnetic.get_coil().has_value());
    REQUIRE(magnetic.get_core().has_value());
    const MAS::Coil coil = *magnetic.get_coil();            // by-value optional getters:
    const MAS::MagneticCore core = *magnetic.get_core();    // copy before use, never bind refs

    mvb::MagneticBuilder builder;
    std::vector<mvb::NamedShape> turns = builder.buildTurnsNamed(coil, core);

    // 24 primary + 12 secondary turns, every one a real shape.
    REQUIRE(turns.size() == 36);

    const double lateralAxisX = -0.0180625;  // core columns[2] axis (E42, MAS metres)
    int nPrimary = 0, nSecondary = 0;
    for (const auto& ns : turns) {
        REQUIRE(!ns.shape.IsNull());
        const BBox b = bbox_of(ns.shape);
        const double cx = 0.5 * (b.xmin + b.xmax);
        const double cz = 0.5 * (b.zmin + b.zmax);
        INFO(ns.name << "  x=[" << b.xmin << "," << b.xmax << "] z=[" << b.zmin << "," << b.zmax << "]");
        if (ns.name.rfind("Primary", 0) == 0) {
            ++nPrimary;
            // Legacy loop around the origin: x/z centred on the main column axis.
            CHECK(std::abs(cx) < 1e-4);
            CHECK(std::abs(cz) < 1e-4);
        } else {
            ++nSecondary;
            // Loop around the lateral leg: centred on its axis, both crossings inside
            // the physical span (window side ~-13.2 mm, outer side ~-22.9 mm).
            CHECK(std::abs(cx - lateralAxisX) < 1e-4);
            CHECK(std::abs(cz) < 1e-4);
            CHECK(b.xmin > -0.0235);
            CHECK(b.xmax < -0.0125);
        }
    }
    CHECK(nPrimary == 24);
    CHECK(nSecondary == 12);
}

// ABT #2489: d3fc7fe gave a design wound on a lateral leg that leg's own former, and refused to build it from a
// MAS::Coil. A processed MAS document (one carrying a geometricalDescription) is drawn from exactly that type, so
// the idealised export of every such design threw instead of drawing.
TEST_CASE("A processed MAS magnetic wound on a lateral leg exports with that leg's former",
          "[bobbin][multicolumn][abt2489]") {
    std::ifstream f(std::string(MAS_COMPLETE_DIR) + "/multicolumn_e42_transformer.json");
    REQUIRE(f.is_open());
    json j;
    f >> j;
    MAS::Magnetic magnetic = mvb::magnetic_from_json(j);    // a full MAS document, read as the bindings do
    const MAS::MagneticCore core = magnetic.get_core().value();
    REQUIRE(core.get_geometrical_description().has_value());    // the path that skips MKF's autocomplete

    mvb::MagneticBuilder builder;
    const std::vector<mvb::NamedShape> all = builder.buildAllNamed(magnetic);

    const double lateralAxisX = -0.0180625;  // core columns[2] axis (E42, MAS metres)
    int formers = 0, onLateralLeg = 0, onMainColumn = 0;
    for (const auto& ns : all) {
        if (ns.role != mvb::Role::Bobbin) continue;
        for (TopExp_Explorer it(ns.shape, TopAbs_SOLID); it.More(); it.Next()) {
            const BBox b = bbox_of(it.Current());
            const double cx = 0.5 * (b.xmin + b.xmax);
            INFO(ns.name << "  x=[" << b.xmin << "," << b.xmax << "]");
            ++formers;
            if (std::abs(cx - lateralAxisX) < 1e-3) ++onLateralLeg;
            if (std::abs(cx) < 1e-3) ++onMainColumn;
        }
    }
    INFO(formers << " former solid(s)");
    CHECK(onLateralLeg > 0);    // the Secondary's leg
    CHECK(onMainColumn > 0);    // the Primary's column
}

// ABT #2489. Where two wound columns share a winding window each former's flanges used to fill it: 607 mm^3 of one
// former inside the other. Alf, 2026-11-23, on the builder STEP of flanges cut at the midplane: "the flanges have to
// have all the width (in X and Z axis) of the minimum of all the sizes (the interior part, the ones that almost
// collide with the other). More than that doesn't make sense". A flange stands out from its tube by the same width on
// every side, x and z, and that width is the room towards the neighbouring former: half the gap between the tubes.
TEST_CASE("Formers of two wound columns in one window have flanges as wide, on every side, as half the gap between their tubes",
          "[bobbin][multicolumn][abt2489]") {
    std::ifstream f(std::string(MAS_COMPLETE_DIR) + "/multicolumn_e42_transformer.json");
    REQUIRE(f.is_open());
    json j;
    f >> j;
    MAS::Magnetic magnetic = mvb::magnetic_from_json(j);
    mvb::MagneticBuilder builder;
    const std::vector<mvb::NamedShape> all = builder.buildAllNamed(magnetic);

    std::vector<TopoDS_Shape> solids;
    for (const auto& ns : all) {
        if (ns.role != mvb::Role::Bobbin) continue;
        for (TopExp_Explorer it(ns.shape, TopAbs_SOLID); it.More(); it.Next()) solids.push_back(it.Current());
    }
    REQUIRE(solids.size() == 6);   // two formers: a tube and two flanges each

    double overlap = 0.0;
    for (std::size_t i = 0; i < solids.size(); ++i) {
        for (std::size_t k = i + 1; k < solids.size(); ++k) {
            BRepAlgoAPI_Common common(solids[i], solids[k]);
            REQUIRE(common.IsDone());
            GProp_GProps props;
            BRepGProp::VolumeProperties(common.Shape(), props);
            overlap += std::abs(props.Mass());
        }
    }
    INFO("former in former: " << overlap * 1e9 << " mm^3");
    CHECK(overlap < 1e-13);    // 1e-4 mm^3

    // The tubes are the two tall solids; every other solid is a flange of the tube it is centred on. The boxes are
    // the exact ones: Bnd_Box's default pads each face by its tolerance (0.15 um here), which is what is measured.
    const auto exactBox = [](const TopoDS_Shape& s) {
        Bnd_Box b;
        BRepBndLib::AddOptimal(s, b, /*useTriangulation=*/false, /*useShapeTolerance=*/false);
        BBox r{};
        b.Get(r.xmin, r.ymin, r.zmin, r.xmax, r.ymax, r.zmax);
        return r;
    };
    std::vector<BBox> tubes, flanges;
    for (const auto& solid : solids) {
        const BBox b = exactBox(solid);
        ((b.ymax - b.ymin) > 0.01 ? tubes : flanges).push_back(b);
    }
    REQUIRE(tubes.size() == 2);
    REQUIRE(flanges.size() == 4);
    std::sort(tubes.begin(), tubes.end(), [](const BBox& a, const BBox& b) { return a.xmin < b.xmin; });
    // Leg tube wall at -14.05 mm, main tube wall at -6.975 mm: 7.075 mm between them, 3.5375 mm for each former.
    const double room = 0.5 * (tubes[1].xmin - tubes[0].xmax);
    INFO("half the gap between the tubes: " << room * 1e3 << " mm");
    CHECK(std::abs(room - 0.0035375) < 1e-9);
    for (const BBox& flange : flanges) {
        const BBox& tube = 0.5 * (flange.xmin + flange.xmax) < 0.5 * (tubes[0].xmax + tubes[1].xmin) ? tubes[0] : tubes[1];
        INFO("flange x " << flange.xmin * 1e3 << " .. " << flange.xmax * 1e3 << " mm, z " << flange.zmin * 1e3 << " .. "
                         << flange.zmax * 1e3 << " mm round the tube x " << tube.xmin * 1e3 << " .. " << tube.xmax * 1e3
                         << " mm, z " << tube.zmin * 1e3 << " .. " << tube.zmax * 1e3 << " mm");
        CHECK(std::abs((tube.xmin - flange.xmin) - room) < 1e-7);
        CHECK(std::abs((flange.xmax - tube.xmax) - room) < 1e-7);
        CHECK(std::abs((tube.zmin - flange.zmin) - room) < 1e-7);
        CHECK(std::abs((flange.zmax - tube.zmax) - room) < 1e-7);
    }
}

// ABT #2492 (ABT #2489, regression 2). Primary on the main column, Secondary on the left leg. MKF sent every terminal
// lead to winding window 0's border, so the Secondary's lead ran from its leg across the main column and the Primary;
// the drawing was refused ("the run would fold back"). Each lead now leaves by its own window. And every conductor's
// leads finish on ONE tip plane: the Secondary's stopped 1 mm short of the Primary's, inside the clearance the
// Primary's port needs ("0.317 mm inside its port plane").
TEST_CASE("Real winding on the main column and a lateral leg draws certified, every lead ending on one tip plane",
          "[realwinding][multicolumn][abt2492]") {
    std::ifstream f(std::string(MAS_COMPLETE_DIR) + "/multicolumn_e42_transformer.json");
    REQUIRE(f.is_open());
    json j;
    f >> j;
    const auto enriched = mvb::magnetic_autocomplete_safe(j.at("magnetic"), /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    std::vector<mvb::NamedShape> named;
    // The enamel gate and the port-plane rule run inside the build: a clean build is the certificate.
    REQUIRE_NOTHROW(named = builder.buildAllNamed(enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                                  mvb::DEFAULT_WIRE_POLYGON_SEGMENTS, mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                                  /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                                  /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                                  /*useRealWindingGeometry=*/true, /*femReady=*/false));
    std::vector<double> tips;
    double mainColumnHalfWidth = 0.0;
    for (const auto& ns : named) {
        if (ns.name.find("parallel") == std::string::npos || ns.name.find("terminal") != std::string::npos) continue;
        Bnd_Box box;
        BRepBndLib::AddOptimal(ns.shape, box, /*useTriangulation=*/false, /*useShapeTolerance=*/false);
        double x0, y0, z0, x1, y1, z1;
        box.Get(x0, y0, z0, x1, y1, z1);
        INFO(ns.name << ": x " << x0 * 1e3 << " .. " << x1 * 1e3 << " mm, lead tips at z = " << z0 * 1e3 << " mm");
        tips.push_back(z0);
        // The Secondary wraps the left leg (axis -18.06 mm): none of its copper reaches the main column's side of x = 0.
        if (ns.name.find("Secondary") != std::string::npos) CHECK(x1 < 0.0);
        if (ns.name.find("Primary") != std::string::npos) mainColumnHalfWidth = std::max(mainColumnHalfWidth, x1);
    }
    REQUIRE(tips.size() == 2);
    CHECK(mainColumnHalfWidth > 0.0);
    CHECK(std::abs(tips[0] - tips[1]) < 1e-9);
}

