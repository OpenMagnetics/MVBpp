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
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <TopExp_Explorer.hxx>
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
