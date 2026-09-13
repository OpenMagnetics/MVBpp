// ABT #1248 ([toroidmount]): toroid mounting VERTICAL (MKF Settings toroid_mounting default) /
// HORIZONTAL. In the EXPORTED frame (buildAllNamed / the STEP / buildRealWindingPaths) every toroid
// terminal exits toward -Y and ends in a planar cap on ONE XZ plane, which sits far enough below
// every other solid for OMFEM's port snap (face inset 0.5 cap radii, refused within 0.9 cap radii
// of other copper). VERTICAL stands the ring on its rim (hole axis Z) turned so the terminals sit
// at the bottom; HORIZONTAL lays it flat (hole axis Y), with the entrance lead entering straight
// through the middle of the hole. Only the leads and the rigid placement differ between the two:
// the turns are the same copper.
#include <catch2/catch_test_macros.hpp>

#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"
#include "support/Settings.h"
#include "json.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pln.hxx>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <string>
#include <vector>

#ifndef MAS_EXAMPLES_DIR
#define MAS_EXAMPLES_DIR "."
#endif

using json = nlohmann::json;
using Mounting = mvb::ConductorBuilder::ToroidMounting;
using P3 = std::array<double, 3>;

namespace {

json loadMountFixture(const std::string& name) {
    const auto root = std::filesystem::path{__FILE__}.parent_path();
    std::ifstream f(root / "mas_complete_fixtures" / name);
    if (!f.good()) f = std::ifstream(root / "realwinding_fixtures" / name);
    if (!f.good()) f = std::ifstream(std::filesystem::path(MAS_EXAMPLES_DIR) / name);
    REQUIRE(f.good());
    json j = json::parse(f);
    mvb::patch_dimension_nominals(j);
    return j.contains("magnetic") ? j.at("magnetic") : j;
}

struct MountGuard {
    OpenMagnetics::SettingsGuard<MAS::OrientationEnum> guard;
    explicit MountGuard(Mounting m)
        : guard(OpenMagnetics::Settings::GetInstance(),
                &OpenMagnetics::Settings::get_toroid_mounting,
                &OpenMagnetics::Settings::set_toroid_mounting,
                m == Mounting::Vertical ? MAS::OrientationEnum::VERTICAL
                                        : MAS::OrientationEnum::HORIZONTAL) {}
};

const char* nameOf(Mounting m) { return m == Mounting::Vertical ? "VERTICAL" : "HORIZONTAL"; }

double dist(const P3& a, const P3& b) {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                     (a[2] - b[2]) * (a[2] - b[2]));
}

// Everything one mounting produces for one design, in the exported frame.
struct Mounted {
    OpenMagnetics::Magnetic enriched;
    mvb::ConductorBuilder::ToroidMountingFrame frame;
    std::vector<mvb::ConductorBuilder::PathPolyline> paths;
    std::map<std::string, mvb::ConductorBuilder::TerminalLeadLength> leads;
    Bnd_Box coreBox;   // exported frame
};

Mounted mount(const std::string& fixture, Mounting m) {
    MountGuard g(m);
    Mounted r;
    r.enriched = mvb::magnetic_autocomplete_safe(loadMountFixture(fixture),
                                                 /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    REQUIRE(mvb::MagneticBuilder::toroidMountingOf(r.enriched) == m);
    r.frame = mvb::MagneticBuilder::toroidMountingFrameOf(r.enriched);
    REQUIRE(r.frame.mounting == m);
    // Includes the collision gate: every lead is checked against every other conductor.
    REQUIRE_NOTHROW(r.paths = builder.buildRealWindingPaths(r.enriched));
    REQUIRE_NOTHROW(r.leads = builder.measureTerminalLeadLengths(r.enriched));
    for (const auto& ns : builder.buildCoreNamed(r.enriched.get_core())) {
        const TopoDS_Shape moved = BRepBuilderAPI_Transform(ns.shape, r.frame.toExported, true).Shape();
        BRepBndLib::AddOptimal(moved, r.coreBox, false, false);
    }
    REQUIRE_FALSE(r.coreBox.IsVoid());
    return r;
}

struct Plane { double y; size_t tips; double capRadiusMax; };

// Every terminal cap: normal (0,-1,0), all on one XZ plane; the plane at least
// 0.9 + 0.5 cap radii below every non-drop copper surface and below the core.
Plane requireCapsOnOneXZPlaneBelowEverything(const Mounted& r, const std::string& what) {
    std::vector<std::pair<P3, double>> tips;   // point, wire radius
    double planeLo = std::numeric_limits<double>::max(), planeHi = std::numeric_limits<double>::lowest();
    for (const auto& p : r.paths) {
        for (int k = 0; k < 2; ++k) {
            const P3& e = k == 0 ? p.end0 : p.end1;
            const P3& d = k == 0 ? p.dir0 : p.dir1;
            INFO(what << ": " << p.name << " end" << k << " at (" << e[0] * 1e3 << ", " << e[1] * 1e3
                      << ", " << e[2] * 1e3 << ") mm, outward (" << d[0] << ", " << d[1] << ", "
                      << d[2] << ")");
            CHECK(std::fabs(d[0]) < 1e-9);
            CHECK(std::fabs(d[2]) < 1e-9);
            CHECK(d[1] < -(1.0 - 1e-9));
            planeLo = std::min(planeLo, e[1]);
            planeHi = std::max(planeHi, e[1]);
            tips.push_back({e, p.wireRadius});
        }
    }
    INFO(what << ": tip plane y in [" << planeLo * 1e3 << ", " << planeHi * 1e3 << "] mm");
    CHECK(planeHi - planeLo < 1e-9);
    double capR = 0.0;
    for (const auto& t : tips) capR = std::max(capR, t.second);
    const double guard = (0.9 + 0.5) * capR;   // OMFEM: face inset 0.5 cap radii, 0.9 cap radii clear

    // Non-drop copper: a sampled point whose (x, z) is not a tip's (a drop runs straight in -Y).
    double lowestCopper = std::numeric_limits<double>::max();
    for (const auto& p : r.paths)
        for (const auto& prim : p.prims)
            for (const auto& q : prim) {
                bool onDrop = false;
                for (const auto& t : tips)
                    if (std::hypot(q[0] - t.first[0], q[2] - t.first[2]) < 1e-6) { onDrop = true; break; }
                if (!onDrop) lowestCopper = std::min(lowestCopper, q[1] - p.wireRadius);
            }
    double xmin, ymin, zmin, xmax, ymax, zmax;
    r.coreBox.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    INFO(what << ": plane " << planeLo * 1e3 << " mm; lowest non-drop copper " << lowestCopper * 1e3
              << " mm; core bottom " << ymin * 1e3 << " mm; OMFEM guard 1.4 cap radii = "
              << guard * 1e3 << " mm");
    CHECK(planeLo <= lowestCopper - guard);
    CHECK(planeLo <= ymin - guard);
    return {planeLo, tips.size(), capR};
}

// Hole axis in the exported frame from the core's tight bbox: the thin extent is the axis.
void requireHoleAxis(const Mounted& r, const std::string& what) {
    double xmin, ymin, zmin, xmax, ymax, zmax;
    r.coreBox.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    const double dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
    INFO(what << ": core bbox extents x " << dx * 1e3 << " y " << dy * 1e3 << " z " << dz * 1e3 << " mm");
    if (r.frame.mounting == Mounting::Vertical) {
        CHECK(dz < 0.9 * dx);
        CHECK(dz < 0.9 * dy);
        CHECK(std::fabs(dx - dy) < 1e-6);   // the ring's round outline stands in XY
    }
    else {
        CHECK(dy < 0.9 * dx);
        CHECK(dy < 0.9 * dz);
        CHECK(std::fabs(dx - dz) < 1e-6);
    }
}

// Is a sampled path primitive one of the terminal-lead pieces (same endpoints)?
bool isLeadPrim(const std::vector<P3>& prim, const Mounted& r) {
    for (const auto& [w, t] : r.leads)
        for (const auto& e : t.per_end)
            for (const auto& pc : e.pieces)
                if (dist(prim.front(), pc.start) < 1e-9 && dist(prim.back(), pc.end) < 1e-9)
                    return true;
    return false;
}

void printLeadNumbers(const Mounted& r, const std::string& what) {
    std::cout << std::setprecision(9);
    for (const auto& [w, t] : r.leads) {
        std::cout << "[toroidmount] " << what << " '" << w << "': terminal lead total "
                  << t.total_m * 1e3 << " mm over " << t.parallels << " parallel(s)\n";
        for (const auto& e : t.per_end) {
            std::cout << "[toroidmount]   parallel " << e.parallel << " " << e.end << ": "
                      << e.length_m * 1e3 << " mm =";
            for (const auto& pc : e.pieces) std::cout << " " << pc.kind << "(" << pc.length_m * 1e3 << ")";
            std::cout << "\n";
        }
    }
}

void requireMountingContract(const std::string& fixture, Mounting m, size_t expectedTips) {
    const std::string what = fixture + " " + nameOf(m);
    const Mounted r = mount(fixture, m);
    requireHoleAxis(r, what);
    const Plane pl = requireCapsOnOneXZPlaneBelowEverything(r, what);
    CHECK(pl.tips == expectedTips);
    printLeadNumbers(r, what);
}

} // namespace

TEST_CASE("Toroid mounting: MKF's default is VERTICAL", "[toroidmount]") {
    OpenMagnetics::Settings::GetInstance().reset();
    CHECK(OpenMagnetics::Settings::GetInstance().get_toroid_mounting() == MAS::OrientationEnum::VERTICAL);
}

TEST_CASE("Toroid mounting: buck caps on one XZ plane, -Y, both mountings", "[toroidmount]") {
    requireMountingContract("buck_inductor_complete.json", Mounting::Vertical, 6);
    requireMountingContract("buck_inductor_complete.json", Mounting::Horizontal, 6);
}

TEST_CASE("Toroid mounting: CMC caps on one XZ plane, -Y, both mountings", "[toroidmount]") {
    requireMountingContract("common_mode_choke_complete.json", Mounting::Vertical, 4);
    requireMountingContract("common_mode_choke_complete.json", Mounting::Horizontal, 4);
}

TEST_CASE("Toroid mounting: current transformer caps on one XZ plane, -Y, both mountings",
          "[toroidmount]") {
    requireMountingContract("current_transformer_complete.json", Mounting::Vertical, 4);
    requireMountingContract("current_transformer_complete.json", Mounting::Horizontal, 4);
}

TEST_CASE("Toroid mounting: 05_pfc and 12_boost caps on one XZ plane, both mountings",
          "[toroidmount]") {
    for (const char* f : {"05_pfc_inductor_t4020_hf60.json", "12_boost_inductor_t5026_26.json"}) {
        requireMountingContract(f, Mounting::Vertical, 2);
        requireMountingContract(f, Mounting::Horizontal, 2);
    }
}

TEST_CASE("Toroid mounting: 3-winding CMC in both mountings", "[toroidmount]") {
    requireMountingContract("realwinding_cmc_3w_1layer.json", Mounting::Vertical, 6);
    requireMountingContract("realwinding_cmc_3w_1layer.json", Mounting::Horizontal, 6);
}

TEST_CASE("Toroid mounting: vertical terminals sit at the bottom of the standing ring",
          "[toroidmount]") {
    // Buck (one winding, terminal gap) and current transformer: every terminal crossing is in the
    // lower half of the ring. A lead's crossing is where its piece chain meets the turn.
    for (const char* f : {"buck_inductor_complete.json", "current_transformer_complete.json"}) {
        const Mounted r = mount(f, Mounting::Vertical);
        for (const auto& [w, t] : r.leads)
            for (const auto& e : t.per_end) {
                const P3& junction = e.end == "entrance" ? e.pieces.back().end : e.pieces.front().start;
                INFO(f << " '" << w << "' parallel " << e.parallel << " " << e.end << " crossing at ("
                       << junction[0] * 1e3 << ", " << junction[1] * 1e3 << ", " << junction[2] * 1e3
                       << ") mm");
                CHECK(junction[1] < 0.0);
            }
    }
}

TEST_CASE("Toroid mounting: vertical 2-winding CMC has one winding on each side",
          "[toroidmount]") {
    const Mounted r = mount("common_mode_choke_complete.json", Mounting::Vertical);
    std::map<std::string, std::pair<double, size_t>> xSum;   // winding -> (sum x, count) of turn copper
    for (const auto& p : r.paths) {
        const std::string w = p.name.substr(0, p.name.find(" parallel "));
        for (const auto& prim : p.prims) {
            if (isLeadPrim(prim, r)) continue;
            for (const auto& q : prim) { xSum[w].first += q[0]; xSum[w].second += 1; }
        }
    }
    REQUIRE(xSum.size() == 2);
    const double x1 = xSum.at("Winding 1").first / xSum.at("Winding 1").second;
    const double x2 = xSum.at("Winding 2").first / xSum.at("Winding 2").second;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    r.coreBox.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    const double outerR = 0.5 * (xmax - xmin);
    INFO("winding copper centroid x: Winding 1 " << x1 * 1e3 << " mm, Winding 2 " << x2 * 1e3
                                                 << " mm, ring outer radius " << outerR * 1e3 << " mm");
    CHECK(x1 < -0.25 * outerR);   // winding 1 on the left
    CHECK(x2 > 0.25 * outerR);    // winding 2 on the right
    // The terminals sit in the gaps between the windings, on the vertical centreline.
    for (const auto& [w, t] : r.leads)
        for (const auto& e : t.per_end) {
            const P3& junction = e.end == "entrance" ? e.pieces.back().end : e.pieces.front().start;
            INFO("'" << w << "' " << e.end << " crossing x " << junction[0] * 1e3 << " mm");
            CHECK(std::fabs(junction[0]) < 0.25 * outerR);
        }
}

TEST_CASE("Toroid mounting: the turns are the same copper in both mountings (rigid motion)",
          "[toroidmount]") {
    for (const char* f : {"buck_inductor_complete.json", "common_mode_choke_complete.json",
                          "current_transformer_complete.json"}) {
        const Mounted v = mount(f, Mounting::Vertical);
        const Mounted h = mount(f, Mounting::Horizontal);
        REQUIRE(v.paths.size() == h.paths.size());
        const gp_Trsf backV = v.frame.toExported.Inverted();
        const gp_Trsf backH = h.frame.toExported.Inverted();
        size_t compared = 0;
        double worst = 0.0;
        for (size_t i = 0; i < v.paths.size(); ++i) {
            REQUIRE(v.paths[i].name == h.paths[i].name);
            std::vector<std::vector<P3>> tv, th;
            for (const auto& prim : v.paths[i].prims) if (!isLeadPrim(prim, v)) tv.push_back(prim);
            for (const auto& prim : h.paths[i].prims) if (!isLeadPrim(prim, h)) th.push_back(prim);
            INFO(f << " '" << v.paths[i].name << "': " << tv.size() << " vertical vs " << th.size()
                   << " horizontal turn primitives");
            REQUIRE(tv.size() == th.size());
            // A one-turn bore-through conductor (the CT primary) is all lead: nothing to compare.
            for (size_t k = 0; k < tv.size(); ++k) {
                REQUIRE(tv[k].size() == th[k].size());
                for (size_t q = 0; q < tv[k].size(); ++q) {
                    gp_Pnt a(tv[k][q][0], tv[k][q][1], tv[k][q][2]);
                    gp_Pnt b(th[k][q][0], th[k][q][1], th[k][q][2]);
                    a.Transform(backV);
                    b.Transform(backH);
                    worst = std::max(worst, a.Distance(b));
                    ++compared;
                }
            }
        }
        INFO(f << ": " << compared << " turn centreline samples, worst build-frame mismatch "
               << worst * 1e9 << " nm");
        CHECK(compared > 0);
        CHECK(worst < 1e-9);
    }
}

TEST_CASE("Toroid mounting: horizontal entrance enters straight through the middle of the hole",
          "[toroidmount]") {
    for (const char* f : {"current_transformer_complete.json", "buck_inductor_complete.json"}) {
        const Mounted h = mount(f, Mounting::Horizontal);
        const auto windows = h.enriched.get_core().get_winding_windows();
        REQUIRE(!windows.empty());
        REQUIRE(windows[0].get_radial_height());
        const double holeR = *windows[0].get_radial_height();
        for (const auto& [w, t] : h.leads)
            for (const auto& e : t.per_end) {
                if (e.end != "entrance") continue;
                INFO(f << " '" << w << "' parallel " << e.parallel << " entrance: " << e.pieces.size()
                       << " piece(s)");
                REQUIRE(e.pieces.size() == 1);
                const auto& pc = e.pieces.front();
                CHECK(pc.kind == "SEG");
                const double len = dist(pc.start, pc.end);
                REQUIRE(len > 0.0);
                // Along the hole axis (Y), from the plane below up to the first turn's crossing.
                CHECK(std::fabs(pc.end[0] - pc.start[0]) < 1e-9 * len);
                CHECK(std::fabs(pc.end[2] - pc.start[2]) < 1e-9 * len);
                CHECK(pc.end[1] > pc.start[1]);
                // Never outside the hole.
                for (const P3& q : {pc.start, pc.end})
                    CHECK(std::hypot(q[0], q[2]) <= holeR);
                // Collinear with the first turn's inner axial leg: the turn primitive it hands over
                // to starts where the lead ends and leaves along +Y.
                const std::string name = w + " parallel " + std::to_string(e.parallel);
                bool found = false;
                for (const auto& p : h.paths) {
                    if (p.name != name) continue;
                    for (const auto& prim : p.prims) {
                        // The piece it hands over to: the first turn's inner leg, or -- on a
                        // one-turn bore-through conductor, which has no wrap -- the exit's axial leg.
                        if (dist(prim.front(), pc.end) > 1e-9 || dist(prim.back(), pc.end) < 1e-9) continue;
                        const P3 u = {prim[1][0] - prim[0][0], prim[1][1] - prim[0][1], prim[1][2] - prim[0][2]};
                        const double ul = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
                        INFO("first turn primitive leaves along (" << u[0] / ul << ", " << u[1] / ul
                                                                   << ", " << u[2] / ul << ")");
                        CHECK(std::fabs(u[0]) < 1e-9 * ul);
                        CHECK(std::fabs(u[2]) < 1e-9 * ul);
                        CHECK(u[1] > 0.0);
                        found = true;
                        break;
                    }
                }
                CHECK(found);
            }
        // Vertical is unaffected: its entrance is an axial leg plus a drop, not a straight segment.
        const Mounted v = mount(f, Mounting::Vertical);
        for (const auto& [w, t] : v.leads)
            for (const auto& e : t.per_end)
                if (e.end == "entrance") CHECK(e.pieces.size() >= 2);
    }
}

TEST_CASE("Toroid mounting: STEP assembly (FEM) caps are planar, -Y, one XZ plane", "[toroidmount]") {
    for (const char* f : {"buck_inductor_complete.json", "common_mode_choke_complete.json"}) {
        for (Mounting m : {Mounting::Vertical, Mounting::Horizontal}) {
            MountGuard g(m);
            const std::string what = std::string(f) + " " + nameOf(m);
            auto enriched = mvb::magnetic_autocomplete_safe(loadMountFixture(f), true);
            mvb::MagneticBuilder builder;
            std::vector<mvb::NamedShape> named;
            REQUIRE_NOTHROW(named = builder.buildAllNamed(
                                enriched, /*includeBobbin=*/false, /*symmetryPlanes=*/0,
                                /*wirePolygonSegments=*/12, mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                /*paintCoating=*/false, /*emitCoatingShells=*/false,
                                /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                /*useRealWindingGeometry=*/true, /*femReady=*/true));
            Bnd_Box core;
            size_t caps = 0;
            double capY = std::numeric_limits<double>::quiet_NaN(), capR = 0.0;
            for (const auto& ns : named) {
                if (ns.role == mvb::Role::Core) BRepBndLib::AddOptimal(ns.shape, core, false, false);
                // buildAllNamed re-roles appended conductor shapes by name, so select caps by name.
                if (ns.name.find(" terminal ") == std::string::npos) continue;
                const TopoDS_Face face = TopoDS::Face(ns.shape);
                BRepAdaptor_Surface ad(face);
                REQUIRE(ad.GetType() == GeomAbs_Plane);
                gp_Dir n = ad.Plane().Axis().Direction();
                if (face.Orientation() == TopAbs_REVERSED) n.Reverse();
                GProp_GProps props;
                BRepGProp::SurfaceProperties(face, props);
                const gp_Pnt c = props.CentreOfMass();
                INFO(what << ": " << ns.name << " centre (" << c.X() * 1e3 << ", " << c.Y() * 1e3 << ", "
                          << c.Z() * 1e3 << ") mm, outward normal (" << n.X() << ", " << n.Y() << ", "
                          << n.Z() << ")");
                CHECK(n.Y() < -(1.0 - 1e-9));
                if (std::isnan(capY)) capY = c.Y();
                CHECK(std::fabs(c.Y() - capY) < 1e-9);
                capR = std::max(capR, std::sqrt(props.Mass() / std::numbers::pi));
                ++caps;
            }
            CHECK(caps == (std::string(f).find("buck") != std::string::npos ? 6u : 4u));
            REQUIRE_FALSE(core.IsVoid());
            double xmin, ymin, zmin, xmax, ymax, zmax;
            core.Get(xmin, ymin, zmin, xmax, ymax, zmax);
            INFO(what << ": caps at y " << capY * 1e3 << " mm, core bottom " << ymin * 1e3
                      << " mm, cap radius " << capR * 1e3 << " mm; core extents x " << (xmax - xmin) * 1e3
                      << " y " << (ymax - ymin) * 1e3 << " z " << (zmax - zmin) * 1e3 << " mm");
            CHECK(capY <= ymin - 1.4 * capR);
            if (m == Mounting::Vertical) CHECK((zmax - zmin) < 0.9 * (ymax - ymin));
            else CHECK((ymax - ymin) < 0.9 * (zmax - zmin));
        }
    }
}
