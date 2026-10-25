// A rect-column former and the copper wound on it share ONE corner: MKF's.
//
// MKF charges every turn on a rectangular (or irregular, drawn as its bounding rectangle) column
// for a corner of radius  frame.cornerRadius + standoff,  where frame.cornerRadius is the former's
// outside corner (Bobbin::get_column_corner_radius), and ConductorBuilder draws exactly that wrap
// (ABT #1295). The plastic under it must be drawn on the same corner, concentric with every turn.
// Until 2026-10-25 BobbinBuilder drew an IRREGULAR former SHARP and a RECTANGULAR one at 1/4 of
// its smaller side, so on 02_flyback_efd25_3c95 (EFD 25/13/9) and isolated_buck (EFD 30/15/9) the
// square former's corner sat 0.47 mm inside the copper of every innermost turn (gate run
// 2026-10-25_1659: 160 and 64 Bobbin-vs-Primary overlaps).
//
// Both tests run on the REAL corpus design, MAS examples/02_flyback_efd25_3c95.json, unaltered,
// through the same real-winding autocomplete the corpus gate uses. They measure the bobbin as
// BobbinBuilder draws it, BEFORE buildAllNamed cuts it with the copper: the cut hides a former
// that disagrees with the copper (when it succeeds) and leaves slivers in the copper (when it
// does not), so the uncut former is the only honest witness of the corner.
#include <catch2/catch_test_macros.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Bobbin.h"
#include "support/Utils.h"
#include "json.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Bnd_Box.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Lin.hxx>

#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

struct CornerFixture {
    OpenMagnetics::Magnetic enriched;
    TopoDS_Shape formerBody;                 // the bobbin's column tube, uncut
    double mkfFormerCorner = 0.0;            // Bobbin::get_column_corner_radius
    std::vector<gp_Cylinder> formerCorners;  // the body's outer corner cylinders
    // Every bare-copper solid of every conductor, with the turn it belongs to.
    struct Part { std::string conductor, part, turn; TopoDS_Shape solid; };
    std::vector<Part> copper;
};

const CornerFixture& flyback02() {
    static const CornerFixture fx = [] {
        CornerFixture f;
        std::ifstream in(std::string(MAS_EXAMPLES_DIR) + "/02_flyback_efd25_3c95.json");
        REQUIRE(in.good());
        json raw = json::parse(in);
        f.enriched = mvb::magnetic_autocomplete_safe(raw.contains("magnetic") ? raw.at("magnetic") : raw,
                                                     /*useRealWindingGeometry=*/true);

        auto bobbinVar = f.enriched.get_coil().get_bobbin();
        auto* bobbin = std::get_if<OpenMagnetics::Bobbin>(&bobbinVar);
        REQUIRE(bobbin != nullptr);
        REQUIRE(bobbin->get_processed_description());
        // The premise: the corpus design's former is the cornered kind this guards.
        REQUIRE(bobbin->get_processed_description()->get_column_shape() == MAS::ColumnShape::IRREGULAR);
        OpenMagnetics::Bobbin forCorner = *bobbin;
        f.mkfFormerCorner = forCorner.get_column_corner_radius();
        REQUIRE(f.mkfFormerCorner > 0.0);

        mvb::MagneticBuilder builder;
        auto former = builder.buildBobbinNamed(f.enriched.get_coil(), f.enriched.get_core());
        REQUIRE_FALSE(former.shape.IsNull());
        // The column tube is the solid with the smallest footprint across the column axis (Y).
        double smallest = std::numeric_limits<double>::max();
        for (TopExp_Explorer e(former.shape, TopAbs_SOLID); e.More(); e.Next()) {
            Bnd_Box bb;
            BRepBndLib::AddOptimal(e.Current(), bb, false, false);
            double x0, y0, z0, x1, y1, z1;
            bb.Get(x0, y0, z0, x1, y1, z1);
            if ((x1 - x0) * (z1 - z0) < smallest) {
                smallest = (x1 - x0) * (z1 - z0);
                f.formerBody = e.Current();
            }
        }
        REQUIRE_FALSE(f.formerBody.IsNull());
        for (TopExp_Explorer e(f.formerBody, TopAbs_FACE); e.More(); e.Next()) {
            BRepAdaptor_Surface s(TopoDS::Face(e.Current()));
            if (s.GetType() != GeomAbs_Cylinder) continue;
            const gp_Cylinder c = s.Cylinder();
            if (c.Axis().Direction().IsParallel(gp_Dir(0, 1, 0), 1e-9)) f.formerCorners.push_back(c);
        }

        auto named = builder.buildAllNamed(f.enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                           mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                           mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                           /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                           /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                           /*useRealWindingGeometry=*/true, /*femReady=*/false);
        for (const auto& ns : named) {
            if (ns.name.find(" parallel ") == std::string::npos) continue;
            size_t k = 0;
            for (TopExp_Explorer e(ns.shape, TopAbs_SOLID); e.More(); e.Next(), ++k) {
                const std::string part = k < ns.partNames.size() ? ns.partNames[k] : std::string();
                // "<conductor> / turn N -> ..." : the corner is drawn at turn N's station.
                std::string turn;
                const auto slash = part.find(" / turn ");
                const auto arrow = part.find(" -> ");
                if (slash != std::string::npos && arrow != std::string::npos && arrow > slash)
                    turn = ns.name + part.substr(slash + 2, arrow - slash - 2);
                f.copper.push_back({ns.name, part, turn, e.Current()});
            }
        }
        REQUIRE_FALSE(f.copper.empty());
        return f;
    }();
    return fx;
}

double conductingRadiusOf(const CornerFixture& f, const std::string& conductor) {
    OpenMagnetics::Coil coil = f.enriched.get_coil();
    const auto windings = coil.get_functional_description();
    for (size_t i = 0; i < windings.size(); ++i) {
        if (conductor.rfind(windings[i].get_name() + " parallel ", 0) != 0) continue;
        const auto wire = coil.resolve_wire(i);
        REQUIRE(wire.get_conducting_diameter());
        return 0.5 * OpenMagnetics::resolve_dimensional_values(wire.get_conducting_diameter().value());
    }
    FAIL("no winding owns conductor " << conductor);
    return 0.0;
}

double outerRadiusOf(const CornerFixture& f, const std::string& conductor) {
    OpenMagnetics::Coil coil = f.enriched.get_coil();
    const auto windings = coil.get_functional_description();
    for (size_t i = 0; i < windings.size(); ++i) {
        if (conductor.rfind(windings[i].get_name() + " parallel ", 0) != 0) continue;
        const auto wire = coil.resolve_wire(i);
        REQUIRE(wire.get_outer_diameter());
        return 0.5 * OpenMagnetics::resolve_dimensional_values(wire.get_outer_diameter().value());
    }
    FAIL("no winding owns conductor " << conductor);
    return 0.0;
}

}  // namespace

TEST_CASE("02_flyback: the former is drawn with MKF's corner and no copper lies inside it",
          "[realwinding][bobbincorner]") {
    const auto& f = flyback02();

    // The former's four outer corners are rounded, at MKF's radius.
    REQUIRE(f.formerCorners.size() == 4);
    for (const auto& c : f.formerCorners) CHECK(std::abs(c.Radius() - f.mkfFormerCorner) < 1e-12);

    // Every bare-copper solid near the former stays at least its own enamel away from it: the
    // innermost turn sits on the former at its coated radius, so the bare copper is one coating
    // thickness off the plastic on the faces AND round the corners. A former corner inside the
    // copper measures 0 here.
    Bnd_Box formerBox;
    BRepBndLib::AddOptimal(f.formerBody, formerBox, false, false);
    size_t checked = 0, corners = 0;
    double worst = std::numeric_limits<double>::max();
    std::string worstPart;
    for (const auto& p : f.copper) {
        Bnd_Box bb;
        BRepBndLib::Add(p.solid, bb);
        if (bb.Distance(formerBox) > 1e-4) continue;
        const double enamel = outerRadiusOf(f, p.conductor) - conductingRadiusOf(f, p.conductor);
        BRepExtrema_DistShapeShape d(f.formerBody, p.solid);
        REQUIRE(d.IsDone());
        ++checked;
        if (p.part.find(" corner ") != std::string::npos) ++corners;
        if (d.Value() - enamel < worst) { worst = d.Value() - enamel; worstPart = p.part; }
        INFO(p.part << ": " << d.Value() * 1e6 << " um from the former, enamel " << enamel * 1e6 << " um");
        CHECK(d.Value() >= enamel - 1e-7);
    }
    INFO("worst: " << worstPart << " at " << worst * 1e6 << " um beyond its enamel");
    // Non-vacuous: the innermost layer's corners are among the solids judged.
    REQUIRE(corners >= 4);
    REQUIRE(checked > corners);
}

TEST_CASE("02_flyback: a wire corner on the former is concentric with the former's corner, at the "
          "radius MKF charged it",
          "[realwinding][bobbincorner]") {
    const auto& f = flyback02();
    REQUIRE(f.formerCorners.size() == 4);
    OpenMagnetics::Coil coil = f.enriched.get_coil();

    size_t onFormer = 0;
    for (const auto& p : f.copper) {
        if (p.part.find(" corner ") == std::string::npos || p.turn.empty()) continue;
        const auto charged = coil.get_turn_bend_radius(p.turn);
        INFO(p.part << " (turn '" << p.turn << "')");
        REQUIRE(charged);
        // The former corner this piece wraps: the nearest corner axis.
        Bnd_Box bb;
        BRepBndLib::Add(p.solid, bb);
        double x0, y0, z0, x1, y1, z1;
        bb.Get(x0, y0, z0, x1, y1, z1);
        const gp_Pnt mid(0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1));
        const gp_Cylinder* axis = nullptr;
        double best = std::numeric_limits<double>::max();
        for (const auto& c : f.formerCorners) {
            const double dd = gp_Lin(c.Axis()).Distance(mid);
            if (dd < best) { best = dd; axis = &c; }
        }
        REQUIRE(axis != nullptr);
        // A concentric quarter arc of centreline radius R and copper radius r comes exactly R - r
        // from its own centre axis; centred anywhere else, it comes farther.
        const gp_Pnt o = axis->Axis().Location();
        const gp_Dir dir = axis->Axis().Direction();
        const double span = 1.0;   // 1 m each way: far beyond any bobbin
        TopoDS_Edge axisEdge = BRepBuilderAPI_MakeEdge(o.Translated(gp_Vec(dir) * -span),
                                                       o.Translated(gp_Vec(dir) * span)).Edge();
        BRepExtrema_DistShapeShape d(axisEdge, p.solid);
        REQUIRE(d.IsDone());
        const double r = conductingRadiusOf(f, p.conductor);
        const double concentric = charged.value() - r;
        INFO("axis distance " << d.Value() * 1e3 << " mm, MKF bend " << charged.value() * 1e3
                              << " mm, copper radius " << r * 1e3 << " mm, former corner "
                              << axis->Radius() * 1e3 << " mm");
        // MKF's bend for a turn laid ON the former is the former's corner plus the turn's coated
        // radius. That turn wraps the former itself, so its arc is centred on the former's corner.
        const bool laidOnFormer =
            std::abs(charged.value() - axis->Radius() - outerRadiusOf(f, p.conductor)) < 1e-12;
        if (laidOnFormer) {
            CHECK(std::abs(d.Value() - concentric) < 1e-6);
            ++onFormer;
        }
        else {
            // A turn farther out keeps MKF's radius, but a face that rides over a dragback return
            // (ConductorBuilder rectRideFor) carries its corner centres outward with it, so the
            // arc is never CLOSER to the former's corner than the concentric one.
            CHECK(d.Value() >= concentric - 1e-6);
        }
    }
    // Non-vacuous: every innermost-layer turn of 02 contributes its four corners.
    REQUIRE(onFormer >= 4);
}
