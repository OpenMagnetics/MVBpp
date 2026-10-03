#include "mvb/mesh/Mesh2d.h"
// Moved verbatim from OMFEM src/meshing/MasMesher.cpp at 59cc050 (ABT #1588, step 6): namespace omfem ->
// mvb::mesh, MasMeshOptions -> MeshOptions, mesh_from_mas -> mesh2d_from_mas taking the parsed magnetic
// (the caller applies its modelling choices, retype_residual_gaps, first; the OMFEM-binary freshness
// check stays in OMFEM, as for step 2). No logic changed.


#include "mvb/mesh/MeshSupport.h"
#include "mvb/mesh/SizeField.h"
#include <cctype>
#include <cerrno>
#include <cstring>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <unordered_map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <dlfcn.h>
#include <sys/stat.h>
#include <gmsh.h>
#include <mmg/mmg3d/libmmg3d.h>
#include "mvb/MagneticBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/SectionBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Symmetry.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/Resistivity.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/StrayCapacitance.h"   // get_wire_insulation_relative_permittivity
#include "physical_models/Temperature.h"          // lumped thermal network (temperature from losses)
#include "processors/MagneticSimulator.h"         // MKF full loss simulation (for the MKF reference)
#include "support/Utils.h"                         // inputs_autocomplete
#include "support/Settings.h"                      // accessories OUT: pin routing off (WP3)
#include "Definitions.h"                          // resolve_dimensional_values
#include "processors/Inputs.h"
#include "support/Utils.h"
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Bnd_Box.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Tool.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <cfloat>

namespace mvb::mesh {

using nlohmann::json;

// Elements per unit distance from a gap mouth in the graded fringing field (ABT #1305), overridable
// with OMFEM_GAP_GRADE_N. MEASURED 2026-09-21 on two geometries (Five_Turns PQ 20/16, Seven_Turns
// PQ 27/17, 5 um residual gaps, mkf-57's Ampere loops), N = 2/4/8/16:
//   per-turn P_cu, max change N -> 2N:  five 0.05 / 0.21 / 0.25 %   seven 0.09 / 0.37 %
//   gap-adjacent Ampere loops (expect 0 / 1 / 3 A):
//       mid-height strip beside the outer leg  +0.107  -0.023  +0.004  -0.005
//       turn 2 alone                           +1.090  +0.985  +0.993  +1.002
//       turns 1, 2, 3                          +3.143  +2.969  +3.002  +2.989
//   elements (five): 331,643 / 391,022 / 503,464 / 760,190; Seven_Turns at N=16 ran out of a
//   16 GB memory cap.
// 8 is the smallest N with per-turn loss converged AND every gap-adjacent loop within ~0.01 A. The
// full-height strip still reads +0.14 A at N=8 (+0.07 at 16), but it does not depend on the distance
// to the plates and the mid-height strip is clean, so it is the outer-leg WALL mm away from the gap,
// resolved at the ordinary core/air size -- a separate matter from the gap grading.
constexpr double kGapGradeN = 8.0;
// Along-gap element size on the gap faces in CORNER mode, in gap lengths (ABT #1292; see the gap block).
constexpr double kGapFaceCapLengths = 2.0;

// The core's columns in the section plane, from the enriched processedDescription: x of the column
// centre and its width along x. Mirrors GapRefinement's column table (MasMesher.hpp) so the 2D and
// 3D paths put a gap in the SAME column -- the nearest by x, which is where MKF places a lateral gap.
struct CoreColumnX { double x; double width; bool central; };
static std::vector<CoreColumnX> core_columns_x(const json& enr) {
    std::vector<CoreColumnX> cols;
    if (!(enr.contains("core") && enr["core"].contains("processedDescription") &&
          enr["core"]["processedDescription"].contains("columns")))
        throw std::runtime_error("core_columns_x: the enriched core has no processedDescription.columns");
    for (const auto& col : enr["core"]["processedDescription"]["columns"]) {
        if (!(col.contains("coordinates") && col["coordinates"].is_array() && !col["coordinates"].empty()))
            throw std::runtime_error("core_columns_x: a column has no coordinates");
        cols.push_back({col["coordinates"].at(0).get<double>(), col.at("width").get<double>(),
                        col.value("type", std::string()) == "central"});
    }
    if (cols.empty()) throw std::runtime_error("core_columns_x: the core has no columns");
    return cols;
}


std::string mesh2d_from_mas(json magnetic, const MeshOptions& opt) {

    // --- MVB++: enrich + build 3D solids + 2D section (axisymmetric +X half). ---
    auto step = [](const char* s){ std::fprintf(stderr, "[mesh] %s\n", s); };
    // Stacked cores (numberStacks N>1) are N identical cores stacked in DEPTH (z). The
    // geometry is generated from numberStacks during enrichment (the per-piece
    // geometricalDescription, which buildAllNamed consumes), so reduce to a SINGLE stack
    // BEFORE enriching. MKF's processedDescription for the full N-stack core already reports
    // the TOTAL depth and effective area (N x single); problem_json_from_mas enriches the
    // UNMODIFIED magnetic separately and reads that total into component_depth. So one stack's
    // cross-section extruded by the total depth = the correct effective area. Building all N
    // stacks instead puts the mid-depth section plane (z=0) exactly between stacks for even N,
    // and the section misses the core entirely (B=0, L collapses ~14x).
    int64_t n_stacks = 1;
    if (magnetic.contains("core") && magnetic["core"].contains("functionalDescription"))
        n_stacks = magnetic["core"]["functionalDescription"].value("numberStacks", int64_t{1});
    if (n_stacks > 1) {
        magnetic["core"]["functionalDescription"]["numberStacks"] = 1;
        magnetic["core"].erase("geometricalDescription");   // force regeneration for 1 stack
        std::fprintf(stderr, "[mesh] stacked core: building 1 of %ld stacks "
                     "(component_depth stays total via the unmodified magnetic)\n",
                     static_cast<long>(n_stacks));
    }
    step("enrich");
    OpenMagnetics::Magnetic enriched = mvb::magnetic_autocomplete_safe(magnetic);
    // Enriched json (gap coordinates + processed columns filled) for adaptive gap sizing. Cheap:
    // serialises the already-enriched object, no second autocomplete.
    json enr_j; OpenMagnetics::to_json(enr_j, enriched);
    const std::vector<CoreGap> core_gaps = extract_core_gaps(enr_j);
    const double post_half = central_column_half_width(enr_j);
    step("buildAllNamed");
    mvb::MagneticBuilder builder;
    // Geometry comes from MVB++. For the magnetics mesh at the conducting diameter we ask
    // MVB++ for the COPPER footprint directly (paintCoating=false) — it resolves the
    // conducting cross-section for round, rectangular, foil AND litz wire. OMFEM therefore
    // no longer carves a copper disk itself (which only worked for centred round wire and
    // gave up on foil/rect). The electrostatic mesh still needs the OUTER footprint so the
    // enamel annulus exists as a dielectric, so it keeps painting the coating.
    const bool magnetics_copper = opt.conductor_conducting_diameter && !opt.electrostatic;
    const bool paint_coating = !magnetics_copper;
    // A toroid core is a smooth annulus — build it with TRUE circles (0 polygon segments)
    // so the ring and bore are exact, not faceted (the section then yields an exact annular
    // area instead of a polygon approximation). Other cores keep the requested faceting.
    const bool is_toroid_core = enriched.get_core().get_type() == MAS::CoreType::TOROIDAL;
    const int core_segs = is_toroid_core ? 0 : opt.polygon_segments;
    // Mesh the bobbin too (magnetically inert -> the field solver treats any non-core/non-turn
    // region as air) so the THERMAL FEA can conduct winding->bobbin->core. Opt out with
    // OMFEM_NO_BOBBIN if a bobbin geometry ever breaks meshing.
    const bool include_bobbin = !std::getenv("OMFEM_NO_BOBBIN");
    std::vector<mvb::NamedShape> named =
        builder.buildAllNamed(enriched, include_bobbin, /*symmetryPlanes=*/0,
                              opt.polygon_segments, core_segs, paint_coating);
    std::fprintf(stderr, "[mesh] built %zu solids\n", named.size());

    // Core angular coverage around the winding, for the 2.5D angular-weighting thermal model:
    // the winding is a rotationally-symmetric loop, but the core is NOT -- at some angles a leg/
    // yoke sits behind the turn (covered -> conducts into the core) and at others it is the open
    // window (exposed -> convects to air). f_covered = that covered fraction. Found by marching a
    // ray outward from just past the winding OD at many angles around the post axis and testing
    // whether core material lies behind the winding there. Written to <mesh>.coverage.
    const NameContext nctx = name_context(enriched);   // ABT #1169: positive core/bobbin names
    {
        std::vector<mvb::NamedShape> coreS, turnS;
        for (auto& ns : named) {
            std::string w; int i = 0;
            const bool isb = ns.name.rfind("Bobbin",0)==0 || ns.name.rfind("FR4",0)==0;
            const std::string r = isb ? "bobbin" : classify(ns.name, w, i, nctx);
            if (r == "core") coreS.push_back(ns); else if (r == "turn") turnS.push_back(ns);
        }
        if (!coreS.empty() && !turnS.empty()) {
            const auto tb = mvb::aggregate_bbox({turnS[0]});         // single turn -> loop normal = axis
            const double td[3] = {tb.xmax-tb.xmin, tb.ymax-tb.ymin, tb.zmax-tb.zmin};
            const int axis = (td[0]<=td[1] && td[0]<=td[2]) ? 0 : (td[1]<=td[2] ? 1 : 2);
            const auto wb = mvb::aggregate_bbox(turnS);
            const double ext[3] = {wb.xmax-wb.xmin, wb.ymax-wb.ymin, wb.zmax-wb.zmin};
            const double mid[3] = {0.5*(wb.xmin+wb.xmax), 0.5*(wb.ymin+wb.ymax), 0.5*(wb.zmin+wb.zmax)};
            const int d1 = (axis+1)%3, d2 = (axis+2)%3;
            const double r_in = 0.5*std::max(ext[d1], ext[d2]);     // winding OD/2
            const auto ab = mvb::aggregate_bbox(named);
            const double aext[3] = {ab.xmax-ab.xmin, ab.ymax-ab.ymin, ab.zmax-ab.zmin};
            const double r_max = 0.5*std::max(aext[d1], aext[d2]);  // out to the core OD
            std::vector<BRepClass3d_SolidClassifier*> cls;   // one per actual solid (shapes may be compounds)
            for (auto& ns : coreS)
                for (TopExp_Explorer ex(ns.shape, TopAbs_SOLID); ex.More(); ex.Next())
                    cls.push_back(new BRepClass3d_SolidClassifier(TopoDS::Solid(ex.Current())));
            if (std::getenv("OMFEM_DEBUG_FACES")) {
                std::fprintf(stderr, "[coverage] extracted %zu core solids from %zu shapes\n", cls.size(), coreS.size());
                auto test = [&](double X,double Y,double Z){ gp_Pnt pt(X,Y,Z); int st=-1;
                    for(auto*c:cls){c->Perform(pt,1e-7); if(c->State()==TopAbs_IN){st=1;break;} if(c->State()==TopAbs_ON)st=0;}
                    std::fprintf(stderr,"[cov-test] (%.1f,%.1f,%.1f)mm -> %s\n",X*1e3,Y*1e3,Z*1e3, st==1?"IN":st==0?"ON":"OUT"); };
                test(0,0,0); test(0.014,0,0); test(0,0,0.0); test(0.015,0,0); test(0,0,0.012);
            }
            const int N = 180, M = 12; int ncov = 0;
            // Sample off the mid-plane: the two core halves meet at axis=mid, so points there read
            // ON the interface, not IN. A quarter of the winding axial extent up sits inside one half.
            const double axis_off = mid[axis] + 0.25 * ext[axis];
            for (int k = 0; k < N; k++) {
                const double phi = 2.0*M_PI*k/N; bool covered = false;
                for (int m = 0; m <= M && !covered; m++) {
                    const double rr = r_in + (r_max - r_in)*m/M;
                    double p[3]; p[axis] = axis_off;
                    p[d1] = mid[d1] + rr*std::cos(phi); p[d2] = mid[d2] + rr*std::sin(phi);
                    const gp_Pnt pt(p[0], p[1], p[2]);
                    for (auto* c : cls) { c->Perform(pt, 1e-7);
                        if (c->State() == TopAbs_IN || c->State() == TopAbs_ON) { covered = true; break; } }
                }
                if (covered) ++ncov;
            }
            for (auto* c : cls) delete c;
            const double f_cov = double(ncov)/N;
            if (std::getenv("OMFEM_DEBUG_FACES")) {
                const auto cb = mvb::aggregate_bbox(coreS);
                std::fprintf(stderr, "[coverage] %zu core solids, core bbox x[%.1f,%.1f] y[%.1f,%.1f] z[%.1f,%.1f] mm\n",
                    coreS.size(), cb.xmin*1e3,cb.xmax*1e3,cb.ymin*1e3,cb.ymax*1e3,cb.zmin*1e3,cb.zmax*1e3);
                std::fprintf(stderr, "[coverage] axis=%d d1=%d d2=%d r_in=%.2f r_max=%.2f mid=(%.2f,%.2f,%.2f)mm\n",
                    axis,d1,d2, r_in*1e3, r_max*1e3, mid[0]*1e3,mid[1]*1e3,mid[2]*1e3);
            }
            std::fprintf(stderr, "[coverage] core angular coverage f_covered=%.3f (axis=%d, winding OD/2=%.2f mm)\n",
                         f_cov, axis, r_in*1e3);
            std::ofstream cf(opt.out_msh + ".coverage"); if (cf) cf << f_cov << "\n";
        }
    }

    if (std::getenv("OMFEM_DEBUG_FACES")) {
        auto bb = [](const std::vector<mvb::NamedShape>& v){ auto b=mvb::aggregate_bbox(v);
            std::fprintf(stderr, "  bbox x[%.2f,%.2f] y[%.2f,%.2f] z[%.2f,%.2f] mm\n",
                b.xmin*1e3,b.xmax*1e3,b.ymin*1e3,b.ymax*1e3,b.zmin*1e3,b.zmax*1e3); };
        std::vector<mvb::NamedShape> cores, turns;
        for (auto& ns : named) { std::string w; int i; auto r=classify(ns.name,w,i,nctx);
            if (r=="core") cores.push_back(ns); else if (r=="turn"){ if(turns.size()<2) turns.push_back(ns); } }
        std::fprintf(stderr, "[bbox] all:"); bb(named);
        if(!cores.empty()){ std::fprintf(stderr,"[bbox] core:"); bb(cores); }
        if(!turns.empty()){ std::fprintf(stderr,"[bbox] turn0:"); bb({turns[0]}); }
    }
    if (const char* sd = std::getenv("OMFEM_DUMP_SECTIONS")) {
        for (const char* pl : {"XY", "XZ", "YZ"}) {
            try {
                auto faces = mvb::SectionBuilder::cut2DFaces(named, mvb::parseSectionPlane(pl), 0.0);
                int ncore = 0, nturn = 0;
                for (auto& ns : faces) { std::string w; int i; auto r = classify(ns.name, w, i, nctx);
                    if (r == "core") ncore++; else if (r == "turn") nturn++; }
                const std::string path = std::string(sd) + "_" + pl + ".step";
                std::ostringstream sink; auto* old = std::cout.rdbuf(sink.rdbuf());
                mvb::exportSTEP(faces, path); std::cout.rdbuf(old);
                std::fprintf(stderr, "[sections] %s: %zu faces (core=%d turn=%d) -> %s\n",
                             pl, faces.size(), ncore, nturn, path.c_str());
            } catch (const std::exception& e) {
                std::fprintf(stderr, "[sections] %s: FAILED (%s)\n", pl, e.what());
            }
        }
    }
    if (const char* p3d = std::getenv("OMFEM_DUMP_3D_STEP")) {
        // Export the full 3D named solids (pre-symmetry, pre-section) for review. This is
        // the geometry buildAllNamed produces — a route to the 3D STEP when the standalone
        // mvbpp_step_generator (drawMagnetic path) crashes on the same core.
        std::ostringstream sink; std::streambuf* old = std::cout.rdbuf(sink.rdbuf());
        const bool ok = mvb::exportSTEP(named, p3d);
        std::cout.rdbuf(old);
        std::fprintf(stderr, "[mesh] dumped 3D STEP -> %s (%s)\n", p3d, ok ? "ok" : "FAILED");
    }
    // Section/cut plane (default XY). OMFEM_SECTION_PLANE=YZ|XZ is a diagnostic to see how the 2D
    // temperature depends on which plane is sliced. Rather than section a different plane (which
    // collapses a coordinate and degenerates the 2D mesh), ROTATE the 3D geometry 90 deg so the
    // requested plane lands in XY, then run the normal (validated) XY pipeline -- boundary tagging,
    // half-symmetry and all. YZ (normal +X) -> rotate +90 about Y; XZ (normal +Y) -> +90 about X.
    if (const char* sp = std::getenv("OMFEM_SECTION_PLANE")) {
        const std::string secpl = sp;
        if (secpl == "YZ" || secpl == "XZ") {
            gp_Trsf t;
            t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0),
                                 secpl == "YZ" ? gp_Dir(0, 1, 0) : gp_Dir(1, 0, 0)), M_PI / 2.0);
            for (auto& ns : named) ns.shape = BRepBuilderAPI_Transform(ns.shape, t, true).Shape();
            std::fprintf(stderr, "[mesh] rotated geometry so %s cross-section lands in XY\n", secpl.c_str());
        }
    }
    if (opt.planar) {
        // Planar (E-type) model: section at the mid-depth XY plane, solved as 2D Cartesian
        // with the core depth out-of-plane. HALF-SYMMETRY: if the solids are mirror-
        // symmetric about the centre-leg plane (YZ, x=0) — as a centred E-core is — keep
        // only the +X half (the centre leg is cut at x=0). The x=0 cut edge is tagged
        // "axis" downstream; both 2D solvers (magnetostatic + harmonic) pin A=0 there and
        // scale L/loss by the multiplicity, so the result is unchanged while the mesh is
        // halved. Falls back to the full section when not YZ-symmetric (e.g. an asymmetric
        // winding layout). OMFEM_NO_SYM forces the full section (A/B validation).
        // A TOROID must never be half-cut: although the ring is YZ-symmetric, the flux runs
        // azimuthally around the FULL ring and the winding wraps the whole circumference —
        // cutting at x=0 and pinning A=0 there severs the closed flux path and discards half
        // the turns (CMC example 07 collapsed to -99.8%). Keep the full ring.
        bool half = !std::getenv("OMFEM_NO_SYM") && !is_toroid_core;
        if (half) {
            half = false;
            for (auto pl : mvb::analyze_symmetry(named).valid_planes)
                if (pl == mvb::SymmetryPlane::YZ) { half = true; break; }
        }
        if (half) {
            step("half-symmetry: cut YZ (+X half)");
            named = mvb::cut_to_region(named, {{mvb::SymmetryPlane::YZ, mvb::SymmetryHalf::Positive}},
                                       mvb::aggregate_bbox(named));
        }
        const double cutoff = std::getenv("OMFEM_CUT_OFFSET") ? std::atof(std::getenv("OMFEM_CUT_OFFSET")) : 0.0;
        step(half ? "cut2DFaces XY (+X half)" : "cut2DFaces XY (full planar section)");
        named = mvb::SectionBuilder::cut2DFaces(named, mvb::parseSectionPlane("XY"), cutoff);
    } else if (opt.axisymmetric) {
        // Axisymmetric meridian (r, axial) half-plane. The rotation axis is the in-plane
        // Y axis (x = r), recovered by sectioning the XY plane and keeping +X. Reduce ONLY
        // across YZ (x=0, the radial half) — NEVER across XZ (y=0), which cuts the AXIAL
        // direction and discards half the core (so a gap between two pieces vanishes).
        //
        // The old apply_symmetry("quarter") cut the first 2 valid_planes. For a core whose
        // DEPTH mid-plane (XY) is a symmetry plane (PQ/RM/pot: planes=[XY,YZ,XZ]) it picked
        // XY+YZ and kept full axial — fine. But for a core where XY is NOT symmetric (EP, or
        // any asymmetric two-piece: planes=[YZ,XZ]) it picked YZ+XZ, cutting the axial extent
        // in half -> the gap was lost -> L collapsed (~14x low). Cutting only YZ is correct
        // for every axisymmetric core; the axial extent stays full so the gap sits inside.
        auto planes = mvb::analyze_symmetry(named).valid_planes;
        if (std::getenv("OMFEM_DEBUG_FACES")) {
            std::string s; for (auto p : planes) s += (p==mvb::SymmetryPlane::XY?"XY ":p==mvb::SymmetryPlane::YZ?"YZ ":p==mvb::SymmetryPlane::XZ?"XZ ":"? ");
            std::fprintf(stderr, "[sym] axisym valid_planes = [ %s]\n", s.c_str());
        }
        if (std::find(planes.begin(), planes.end(), mvb::SymmetryPlane::YZ) != planes.end()) {
            step("cut YZ (radial half, +X)");
            named = mvb::cut_to_region(named, {{mvb::SymmetryPlane::YZ, mvb::SymmetryHalf::Positive}},
                                       mvb::aggregate_bbox(named));
        }
        step("cut2DFaces XY");
        named = mvb::SectionBuilder::cut2DFaces(named, mvb::parseSectionPlane("XY"), 0.0);
        step("filter_by_side +X");
        named = mvb::filter_by_side(named, mvb::parse_side_spec("+X"));
    } else {
        named = mvb::SectionBuilder::cut2DFaces(named, mvb::parseSectionPlane("XY"), 0.0);
        named = mvb::filter_by_side(named, mvb::parse_side_spec("+X+Y+Z"));
    }
    std::fprintf(stderr, "[mesh] %zu 2D faces\n", named.size());
    if (named.empty()) throw std::runtime_error("mesh_from_mas: MVB++ produced no 2D faces");

    // For the electrostatic (capacitance) mesh: per-winding wire conducting RADIUS,
    // so each turn face can be split into an inner copper electrode + enamel coating
    // annulus. The conducting diameter comes from the resolved (enriched) wire.
    std::map<std::string, double> conducting_radius_by_winding;
    // Carve a copper electrode + coating annulus only for the ELECTROSTATIC mesh, which
    // needs both regions. The magnetics conducting-diameter mesh already has copper faces
    // straight from MVB++ (paint_coating=false above), so it needs no carve and no
    // hand-resolved conducting radius — this also removes the old non-round-wire fallback.
    const bool split_conductor = opt.electrostatic;
    if (split_conductor) {
        OpenMagnetics::Coil coil = enriched.get_coil();   // mutable copy: resolve_wire is non-const
        const auto& fd = coil.get_functional_description();
        for (size_t i = 0; i < fd.size(); ++i) {
            const std::string wname = fd[i].get_name();
            const auto cd = coil.resolve_wire(i).get_conducting_diameter();
            if (!cd)
                throw std::runtime_error("mesh_from_mas(electrostatic): winding '" + wname +
                                         "' has no resolved wire conductingDiameter");
            // MKF's resolver (nominal -> (min+max)/2 -> max -> min; throws if none), not a
            // hand-read .nominal.
            const double dcond = OpenMagnetics::resolve_dimensional_values(cd.value());
            if (dcond <= 0.0)
                throw std::runtime_error("mesh_from_mas(electrostatic): winding '" + wname +
                                         "' has a non-positive conducting diameter");
            conducting_radius_by_winding[wname] = 0.5 * dcond;
        }
    }

    // Record each named face's (centroid, region) BEFORE the STEP round-trip, so we
    // can classify the imported gmsh surfaces by GEOMETRY (nearest centroid) rather
    // than STEP labels. This works even where gmsh can't read STEP labels (e.g. the
    // wasm build without OCC_CAF), and is robust to label renaming.
    struct NamedCentroid { double cx, cy; std::string region; };
    std::vector<NamedCentroid> name_centroids;
    // The bobbin solid carries the MAS bobbin's own name (e.g. "Bobbin_0"), which the role
    // classifier would otherwise mistake for a core piece. Resolve that name so we can tag it
    // "bobbin": meshed (for the thermal FEA's winding->bobbin->core conduction) but magnetically
    // inert (the field solver treats any non-core/non-turn region as air).
    const std::string& bobbin_name = nctx.bobbin_name;
    int bobbin_idx = 0;
    for (const auto& ns : named) {
        std::string winding; int index = 0;
        std::string role;
        if (ns.name.rfind("Bobbin", 0) == 0 || ns.name.rfind("FR4Board", 0) == 0 ||
            (!bobbin_name.empty() && ns.name.rfind(bobbin_name, 0) == 0)) role = "bobbin";
        else role = classify(ns.name, winding, index, nctx);
        if (role.empty()) continue;
        // ABT #1169: accessory solids (spacer, shunt, divider, sleeve, pin) get their own 2D
        // region so a per-region material can be attached to them; a pin joins the bobbin's
        // thermal bucket. The turn/bobbin/core spellings are unchanged.
        if (role == "bobbin" || role == "pin") index = bobbin_idx++;
        // The 2D section numbers its conductors per turn; a continuous real-winding conductor
        // ("winding") keeps that spelling here so the 2D tagging is unchanged by ABT #1169
        // (the 2D mesher never asks MVB++ for real-winding geometry, but the mapping stays
        // where it was rather than being quietly redefined).
        const std::string turnRegion = "turn_" + winding + "_" + std::to_string(index);
        const std::string region =
            region_for_role(role == "pin" ? "bobbin" : (role == "winding" ? "turn" : role),
                            winding, index, turnRegion);
        for (TopExp_Explorer ex(ns.shape, TopAbs_FACE); ex.More(); ex.Next()) {
            GProp_GProps gp; BRepGProp::SurfaceProperties(TopoDS::Face(ex.Current()), gp);
            const gp_Pnt c = gp.CentreOfMass();
            if (std::getenv("OMFEM_DEBUG_FACES") && role == "core")
                std::fprintf(stderr, "[coreface] %s centroid=(%.3f,%.3f)mm area=%.4g mm^2\n",
                             region.c_str(), c.X()*1e3, c.Y()*1e3, gp.Mass()*1e6);
            name_centroids.push_back({c.X(), c.Y(), region});
        }
    }
    // The section MUST contain a core face. If it doesn't, the core did not intersect the
    // cut plane and the bare-winding inductance (B_core = 0, L ~ 1/400 of the true value)
    // would be silently, confidently wrong — surface it. Stacked cores (the historical cause:
    // an even numberStacks put the mid-plane between stacks) are now reduced to a single stack
    // before the section, so this should only fire on a genuinely degenerate geometry.
    if (std::none_of(name_centroids.begin(), name_centroids.end(),
                     [](const NamedCentroid& nc){ return nc.region.rfind("core_", 0) == 0; }))
        throw std::runtime_error(
            "mesh_from_mas: the 2D section contains no core region — the core did not intersect "
            "the cut plane. The geometry may be degenerate; OMFEM_CUT_OFFSET can move the section.");
    // Diagnostic: radial/axial extent of each winding's turn centroids (leakage is
    // governed by where MVB++ places primary vs secondary in the window).
    {
        std::map<std::string, std::array<double,5>> ext;  // [rmin,rmax,ymin,ymax,count]
        for (const auto& nc : name_centroids) {
            if (nc.region.rfind("turn_", 0) != 0) continue;
            std::string w = nc.region.substr(5); w = w.substr(0, w.find_last_of('_'));
            auto& e = ext.try_emplace(w, std::array<double,5>{DBL_MAX,-DBL_MAX,DBL_MAX,-DBL_MAX,0}).first->second;
            e[0]=std::min(e[0],nc.cx); e[1]=std::max(e[1],nc.cx);
            e[2]=std::min(e[2],nc.cy); e[3]=std::max(e[3],nc.cy); e[4]+=1;
        }
        for (auto& kv : ext)
            std::fprintf(stderr, "[mesh] winding %-12s turns=%.0f  r=[%.2f,%.2f]mm  y=[%.2f,%.2f]mm\n",
                         kv.first.c_str(), kv.second[4], kv.second[0]*1e3, kv.second[1]*1e3,
                         kv.second[2]*1e3, kv.second[3]*1e3);
    }

    // Export the 2D faces to a temporary STEP — the STEP round-trip heals the
    // OCCT geometry so gmsh's kernel imports it cleanly (importing the raw faces
    // directly trips GeomAdaptor_Surface::UContinuity).
    const std::string step_path = opt.out_msh + ".step";
    {
        std::ostringstream sink;                       // silence OCCT STEP-writer chatter
        std::streambuf* old = std::cout.rdbuf(sink.rdbuf());
        const bool ok = mvb::exportSTEP(named, step_path);
        std::cout.rdbuf(old);
        if (!ok) throw std::runtime_error("mesh_from_mas: exportSTEP failed");
    }

    // --- gmsh: import the STEP (with labels), classify each surface, then mesh. ---
    bool own_init = !gmsh::isInitialized();
    if (own_init) gmsh::initialize(); apply_gmsh_threads();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::option::setNumber("Geometry.OCCImportLabels", 1);
    gmsh::model::add("omfem_mas");
    gmsh::vectorpair imported_dimtags;
    // MVB++ exportSTEP owns the metres->millimetres conversion (ABT #317), so the file on disk
    // declares MILLIMETRES while everything in OMFEM -- mesh size targets, sigma, mu0, the
    // solvers -- is SI metres. Without this the round-trip silently returns geometry 1000x too
    // big (measured: a core bbox of x[-6.325,6.325] mm came back as x[-6325,6325]), which also
    // wrecks the bbox-based region matching. Ask OCC to convert to metres on import so the unit
    // declared in the file is honoured whatever the exporter does.
    gmsh::option::setString("Geometry.OCCTargetUnit", "M");
    gmsh::model::occ::importShapes(step_path, imported_dimtags);
    gmsh::model::occ::synchronize();

    struct Imp { int tag; std::string region; double area; };
    std::vector<Imp> imps;
    {
        gmsh::vectorpair surfs; gmsh::model::occ::getEntities(surfs, 2);
        for (auto& dt : surfs) {
            double mass = 0.0; gmsh::model::occ::getMass(2, dt.second, mass);
            if (mass <= 1e-15) continue;
            // match this surface to a named face by nearest centroid (geometry, not labels)
            double gx, gy, gz; gmsh::model::occ::getCenterOfMass(2, dt.second, gx, gy, gz);
            double best = DBL_MAX; const std::string* region = nullptr; int group = -1;
            for (const auto& nc : name_centroids) {
                const double d = (gx-nc.cx)*(gx-nc.cx) + (gy-nc.cy)*(gy-nc.cy);
                if (d < best) { best = d; region = &nc.region; }
            }
            if (!region) continue;
            imps.push_back({dt.second, *region, mass});
        }
    }
    if (imps.empty()) { if (own_init) gmsh::finalize(); throw std::runtime_error("mesh_from_mas: no faces classified"); }

    // bounding box of all imported faces -> surrounding air rectangle.
    double xmin=1e30, ymin=1e30, xmax=-1e30, ymax=-1e30;
    { gmsh::vectorpair ents; gmsh::model::occ::getEntities(ents, 2);
      for (auto& dt : ents) { double a,b,c,d,e,f; gmsh::model::occ::getBoundingBox(2,dt.second,a,b,c,d,e,f);
        xmin=std::min(xmin,a); ymin=std::min(ymin,b); xmax=std::max(xmax,d); ymax=std::max(ymax,e);} }
    const double cx=0.5*(xmin+xmax), cy=0.5*(ymin+ymax);
    const double hx=0.5*(xmax-xmin)*opt.air_margin, hy=0.5*(ymax-ymin)*opt.air_margin;
    const double x_lo = opt.axisymmetric ? std::max(0.0, cx-hx) : (cx-hx);
    const int air_tag = gmsh::model::occ::addRectangle(x_lo, cy-hy, 0, (cx+hx)-x_lo, 2*hy);
    gmsh::model::occ::synchronize();

    // Inside each turn face (drawn at the wire OUTER diameter), add a concentric copper
    // disk at the conducting diameter. After the fragment the inner disk is the copper
    // (electrode for electrostatics; the current-carrying conductor for magnetics) and
    // the leftover ring is the enamel coating (dielectric / air).
    if (split_conductor) {
        constexpr double kPI = 3.14159265358979323846;
        const size_t n_imported = imps.size();
        for (size_t i = 0; i < n_imported; ++i) {
            if (imps[i].region.rfind("turn_", 0) != 0) continue;
            std::string w = imps[i].region.substr(5);     // "<winding>_<index>"
            w = w.substr(0, w.find_last_of('_'));          // "<winding>"
            auto rit = conducting_radius_by_winding.find(w);
            if (rit == conducting_radius_by_winding.end())
                throw std::runtime_error("mesh_from_mas(split): no conducting radius for winding '" + w + "'");
            const double rc = rit->second;
            double gx, gy, gz; gmsh::model::occ::getCenterOfMass(2, imps[i].tag, gx, gy, gz);
            const int cu = gmsh::model::occ::addDisk(gx, gy, 0, rc, rc);
            const std::string creg = "copper_" + imps[i].region.substr(5);
            imps.push_back({cu, creg, kPI * rc * rc});
        }
        gmsh::model::occ::synchronize();
    }

    // fragment air against the imported faces (makes everything conforming).
    gmsh::vectorpair obj{{2, air_tag}}, tools;
    for (auto& im : imps) tools.push_back({2, im.tag});
    gmsh::vectorpair out; std::vector<gmsh::vectorpair> outMap;
    gmsh::model::occ::fragment(obj, tools, out, outMap);
    gmsh::model::occ::synchronize();

    // assign region to each result surface; process imports in DESCENDING area so
    // small specific faces (turns) overwrite the large enclosing core.
    std::map<int,std::string> surf_region;
    std::vector<size_t> order(imps.size());
    for (size_t i=0;i<imps.size();++i) order[i]=i;
    std::sort(order.begin(), order.end(), [&](size_t a,size_t b){return imps[a].area>imps[b].area;});
    for (size_t k=0;k<order.size();++k) {
        size_t i = order[k];
        const auto& mp = outMap[i+1];                 // outMap[0]=air, [i+1]=tools[i]
        for (auto& dt : mp) if (dt.first==2) surf_region[dt.second]=imps[i].region;
    }

    // group surfaces; unassigned = air. Split turns by centroid x sign (+/-).
    std::map<std::string,std::vector<int>> groups;
    gmsh::vectorpair alls; gmsh::model::getEntities(alls, 2);
    // A toroid turn crosses the XY mid-plane at TWO radii — inner (hole-side) and outer — at
    // the SAME azimuth, so an x-sign split tags both the same. Split by RADIUS instead: the
    // wire carries +current through the hole and -current around the outside. Threshold =
    // midway between the closest and farthest turn crossing (cleanly between r_in and r_out).
    double r_mid_toroid = 0.0;
    if (is_toroid_core) {
        double rmn = 1e30, rmx = 0.0;
        for (auto& dt : alls) {
            auto it = surf_region.find(dt.second);
            if (it == surf_region.end() || it->second.rfind("turn_", 0) != 0) continue;
            double gx,gy,gz; gmsh::model::occ::getCenterOfMass(2, dt.second, gx,gy,gz);
            const double r = std::hypot(gx, gy); rmn = std::min(rmn, r); rmx = std::max(rmx, r);
        }
        if (rmx > rmn) r_mid_toroid = 0.5 * (rmn + rmx);
    }
    // Multi-column placement: a turn wound on a LATERAL column crosses the section plane
    // twice on the SAME side of x=0, so the legacy x-sign split would give both crossings
    // the same current direction — a broken loop. Tag by WHICH crossing this face is
    // instead: the one at the turn's own coordinates carries +I ("plus"), the
    // additionalCoordinates return crossing -I ("minus") — MKF's turn orientation. For a
    // main-column turn coordinates[0] is +x and the second crossing (when emitted) is its
    // mirror at -x, so this reproduces the x-sign rule exactly on legacy geometry; turns
    // without additionalCoordinates keep the x-sign rule outright.
    std::map<std::string, std::pair<double,double>> turn_crossings;  // region -> {x_own, x_return}
    {
        auto turnsOpt = enriched.get_coil().get_turns_description();
        if (turnsOpt) {
            const auto turns = *turnsOpt;   // local copy: by-value optional getter
            for (const auto& t : turns) {
                std::smatch m;
                const std::string tn = t.get_name();
                if (!std::regex_match(tn, m, kTurnRe)) continue;
                const auto addl = t.get_additional_coordinates();
                if (!addl || addl->empty() || (*addl)[0].empty()) continue;
                if (t.get_coordinates().empty()) continue;
                turn_crossings["turn_" + m[1].str() + "_" + m[3].str()] =
                    {t.get_coordinates()[0], (*addl)[0][0]};
            }
        }
    }
    auto leg_suffix = [&](const std::string& region, double gx, double gy) -> std::string {
        if (is_toroid_core && r_mid_toroid > 0.0)
            return (std::hypot(gx, gy) < r_mid_toroid) ? "_inner" : "_outer";
        auto it = turn_crossings.find(region);
        if (it != turn_crossings.end()) {
            const double d_own    = std::fabs(gx - it->second.first);
            const double d_return = std::fabs(gx - it->second.second);
            return (d_own <= d_return) ? "_plus" : "_minus";
        }
        return (gx >= 0.0) ? "_plus" : "_minus";
    };
    for (auto& dt : alls) {
        int s = dt.second;
        auto it = surf_region.find(s);
        if (it == surf_region.end()) { groups["air"].push_back(s); continue; }
        std::string reg = it->second;
        if (opt.electrostatic) {
            // turn_<w>_<i> -> coating_<w>_<i> (the ring left after the copper disk);
            // copper_<w>_<i> stays. No leg (+/-) split for the electrostatic solve.
            if (reg.rfind("turn_", 0) == 0) reg = "coating_" + reg.substr(5);
        } else if (opt.conductor_conducting_diameter && split_conductor) {
            // Magnetics: the inner copper disk IS the conductor -> rename to
            // turn_<w>_<i>_<leg> (the solver's normal name); the enamel annulus
            // (the original turn_ face) becomes air (non-conducting).
            if (reg.rfind("copper_", 0) == 0) {
                double gx,gy,gz; gmsh::model::occ::getCenterOfMass(2,s,gx,gy,gz);
                const std::string turnKey = "turn_" + reg.substr(7);
                reg = turnKey + leg_suffix(turnKey, gx, gy);
            } else if (reg.rfind("turn_", 0) == 0) {
                reg = "air";
            }
        } else if (reg.rfind("turn_",0)==0) {
            double gx,gy,gz; gmsh::model::occ::getCenterOfMass(2,s,gx,gy,gz);
            reg += leg_suffix(reg, gx, gy);
        }
        groups[reg].push_back(s);
    }

    int pg = 1;
    for (auto& [name, surfs] : groups) {
        if (surfs.empty()) continue;
        gmsh::model::addPhysicalGroup(2, surfs, pg);
        gmsh::model::setPhysicalName(2, pg, name);
        ++pg;
    }

    // boundary curves: far-field 'outer_air' (box edges, minus axis) and 'axis'.
    double ax_lo=1e30, ay_lo=1e30, ax_hi=-1e30, ay_hi=-1e30;
    for (auto& dt : alls) { double a,b,c,d,e,f; gmsh::model::getBoundingBox(2,dt.second,a,b,c,d,e,f);
        ax_lo=std::min(ax_lo,a); ay_lo=std::min(ay_lo,b); ax_hi=std::max(ax_hi,d); ay_hi=std::max(ay_hi,e);}
    const double eps = 1e-4*std::max({ax_hi-ax_lo, ay_hi-ay_lo, 1e-9});
    std::vector<int> outer, axis;
    gmsh::vectorpair curves; gmsh::model::getEntities(curves, 1);
    for (auto& dt : curves) {
        double a,b,c,d,e,f; gmsh::model::occ::getBoundingBox(1,dt.second,a,b,c,d,e,f);
        bool on_axis = std::fabs(a)<1e-6 && std::fabs(d)<1e-6;
        bool box = (std::fabs(a-ax_lo)<eps&&std::fabs(d-ax_lo)<eps) || (std::fabs(a-ax_hi)<eps&&std::fabs(d-ax_hi)<eps)
                || (std::fabs(b-ay_lo)<eps&&std::fabs(e-ay_lo)<eps) || (std::fabs(b-ay_hi)<eps&&std::fabs(e-ay_hi)<eps);
        if (on_axis) axis.push_back(dt.second);
        else if (box) outer.push_back(dt.second);
    }
    if (!outer.empty()) { gmsh::model::addPhysicalGroup(1, outer, pg); gmsh::model::setPhysicalName(1, pg, "outer_air"); ++pg; }
    if (!axis.empty())  { gmsh::model::addPhysicalGroup(1, axis,  pg); gmsh::model::setPhysicalName(1, pg, "axis"); ++pg; }

    // Adaptive mesh sizing via the shared SizeField module: conductor skin refinement + (optional)
    // coating refinement + air ceiling are composed into ONE Min background field. The old code
    // called setAsBackgroundMesh TWICE (conductor then coating) so only the LAST field survived --
    // when both were active the conductor refinement was silently lost. The Min composition fixes
    // that structurally. Curvature stays off here (2D wires are already refined by the conductor
    // field) to preserve the validated element distribution.
    mvb::mesh::SizingPolicy pol2d; pol2d.N_curv = 0;
    mvb::mesh::SizeFieldBuilder sizer(2, opt.air_target, pol2d);
    sizer.add_ceiling(opt.air_target);
    double cmin = opt.core_target;
    auto collect_curves = [&](auto pred) {
        std::vector<double> cv;
        for (auto& [name, surfs] : groups) {
            if (!pred(name)) continue;
            for (int s : surfs) {
                gmsh::vectorpair b; gmsh::model::getBoundary({{2, s}}, b, false, false, false);
                for (auto& dt : b) if (dt.first == 1) cv.push_back(std::abs(dt.second));
            }
        }
        return cv;
    };
    if (opt.conductor_target > 0.0) {
        auto cond = collect_curves([](const std::string& n){ return n.rfind("turn_", 0) == 0; });
        sizer.add_distance_refinement(cond, opt.conductor_target, opt.conductor_target,
                                      8.0 * opt.conductor_target, 200);
        if (!cond.empty()) cmin = std::min(cmin, opt.conductor_target);
    }
    if (opt.electrostatic && opt.coating_target > 0.0) {
        auto cu = collect_curves([](const std::string& n){
            return n.rfind("copper_", 0) == 0 || n.rfind("coating_", 0) == 0; });
        sizer.add_distance_refinement(cu, opt.coating_target, opt.coating_target,
                                      40.0 * opt.coating_target, 400);
        if (!cu.empty()) cmin = std::min(cmin, opt.coating_target);
    }
    // FERRITE SIZE (ABT #1292 study): without it the ferrite away from the gap and the conductors is meshed at the
    // air ceiling, and the iGSE loss (B^beta, beta ~ 2.9) is resolved only where some other feature's band happens
    // to reach. OMFEM_2D_FERRITE_SIZE [m] gives the ferrite a size of its own.
    if (const char* fs = std::getenv("OMFEM_2D_FERRITE_SIZE")) {
        const double h = std::atof(fs);
        if (!(h > 0.0)) throw std::runtime_error("OMFEM_2D_FERRITE_SIZE must be a positive size in metres");
        std::vector<double> ferrite;
        for (auto& [name, surfs] : groups) if (name.rfind("core", 0) == 0) for (int sf : surfs) ferrite.push_back(sf);
        if (ferrite.empty()) throw std::runtime_error("OMFEM_2D_FERRITE_SIZE: the section has no core face");
        sizer.add_region_size(ferrite, h);
        cmin = std::min(cmin, h);
        std::fprintf(stderr, "[mesh] ferrite: %zu face(s) at %.3g mm\n", ferrite.size(), h * 1e3);
    }
    // Adaptive GAP refinement: each functional gap resolved with >= N_gap elements across its
    // length, in a fringing band around the centre leg. The gap is explicit conforming geometry,
    // so this refines the FRINGING field (gross reluctance is already captured by the air slab);
    // gt only ever refines below core_target, never coarsens. Section is (r=x, axial=y).
    const double gapdiv = std::getenv("OMFEM_GAP_DIV") ? std::atof(std::getenv("OMFEM_GAP_DIV")) : 6.0;
    // ABT #1305: each gap is refined in ITS OWN column. The centre post keeps its box from the axis
    // (it carries the gap's return flux, unchanged); a lateral column gets its own footprint plus the
    // same 4-gap-length fringing margin. The section may hold only ONE of a mirrored pair of legs:
    // PQ 20/16's section is x in [-2.56, 12.81] mm -- a sliver past the axis, so "ax_lo >= 0" is not
    // the test -- and contains the +10.12 mm leg but not the -10.12 mm one. So the question is
    // geometric: if the column's own footprint misses the section and its MIRROR lands in it, the
    // section's leg IS that mirror. A column that misses the section both ways is an error.
    const std::vector<CoreColumnX> gap_cols = core_gaps.empty() ? std::vector<CoreColumnX>{} : core_columns_x(enr_j);
    for (const auto& g : core_gaps) {
        const double gt = std::min(opt.core_target, g.len / gapdiv);
        if (!(gt > 0.0)) continue;
        const double yb = std::max(g.len, gt) + 4.0 * gt;            // axial half-band incl. fringing
        size_t best = 0; double bd = std::numeric_limits<double>::max();
        for (size_t k = 0; k < gap_cols.size(); ++k) {
            const double dd = std::fabs(gap_cols[k].x - g.xc);
            if (dd < bd) { bd = dd; best = k; }
        }
        const CoreColumnX& c = gap_cols[best];
        // cx: where this column sits IN THE SECTION -- its own x, or the mirror leg when only that
        // one is modelled. Decided once, used by both the box and the graded field below.
        double cx = c.x;
        if (!c.central) {
            const double reach = 0.5 * c.width + 4.0 * g.len;
            auto overlaps = [&](double x) { return x + reach > ax_lo && x - reach < ax_hi; };
            cx = overlaps(c.x) ? c.x : (overlaps(-c.x) ? -c.x : c.x);   // -c.x: the mirror leg
        }
        double x0, x1;
        if (c.central) {
            x0 = ax_lo;
            x1 = (post_half > 0.0) ? std::min(ax_hi, post_half + 4.0 * g.len) : ax_hi;
        } else {
            const double reach = 0.5 * c.width + 4.0 * g.len;
            x0 = std::max(ax_lo, cx - reach);
            x1 = std::min(ax_hi, cx + reach);
        }
        if (!(x1 > x0))
            throw std::runtime_error("mesh_from_mas: the refinement box of the gap at x=" + std::to_string(g.xc * 1e3) +
                                     " mm, y=" + std::to_string(g.yc * 1e3) + " mm lies outside the section [" +
                                     std::to_string(ax_lo * 1e3) + ", " + std::to_string(ax_hi * 1e3) + "] mm");
        std::fprintf(stderr, "[mesh] gap y=%.3fmm len=%.4gmm in %s column at x=%.3fmm: box x[%.3f, %.3f]mm, size %.3gum\n",
                     g.yc * 1e3, g.len * 1e3, c.central ? "the CENTRAL" : "a LATERAL", c.x * 1e3, x0 * 1e3, x1 * 1e3, gt * 1e6);
        const double gradeN = std::getenv("OMFEM_GAP_GRADE_N") ? std::atof(std::getenv("OMFEM_GAP_GRADE_N")) : kGapGradeN;
        if (!(gradeN > 0.0)) throw std::runtime_error("OMFEM_GAP_GRADE_N must be positive");
        const double span_lo = c.central ? -post_half : cx - 0.5 * c.width;
        const double span_hi = c.central ?  post_half : cx + 0.5 * c.width;
        const double lo = std::min(span_lo, span_hi), hi = std::max(span_lo, span_hi);
        // The gap's entities, found by position. The tolerance scales with the gap (a quarter of its
        // length), never a fixed epsilon -- OCC bounding boxes inflate by ~0.1 um, well inside it.
        const double tol = 0.25 * g.len;
        // the gap's faces: curves lying within the gap's axial extent, over its column
        std::vector<double> faces;
        {
            gmsh::vectorpair ents; gmsh::model::getEntities(ents, 1);
            for (const auto& e : ents) {
                double bx0, by0, bz0, bx1, by1, bz1;
                gmsh::model::getBoundingBox(e.first, e.second, bx0, by0, bz0, bx1, by1, bz1);
                if (by0 >= g.yc - 0.5 * g.len - tol && by1 <= g.yc + 0.5 * g.len + tol &&
                    bx1 > lo - tol && bx0 < hi + tol)
                    faces.push_back(static_cast<double>(e.second));
            }
        }
        if (faces.empty())
            throw std::runtime_error("mesh_from_mas: no curve of the gap at x=" + std::to_string(g.xc * 1e3) +
                                     " mm, y=" + std::to_string(g.yc * 1e3) + " mm found in the section -- the gap "
                                     "MKF placed is not in the geometry being meshed");
        // CORNER is the default (ABT #1292, measured 2026-09-21 on Five_Turns PQ 20/16, identical input and solver):
        //                 elements   wall    L [uH]   P_cu [mW] (middle turn)   P_core [W]   gap-adjacent loops
        //   face  (#1305)  503,464   6 min   48.570   2.61280 (1.75607)         1.2490       <= 0.0071 A
        //   corner, cap 2   48,276   19 s    48.570   2.61280 (1.75606)         1.2489       <= 0.0088 A
        // (P_core with the volume-quantile iGSE cap; with the old count cap the two read 3.4 % apart.) The
        // multicolumn E42 section went from 4.08M unknowns (bad_alloc / >90 min) to a 13-min, 2.8 GB pass.
        // FACE mode is kept only to reproduce the #1305 reference.
        const char* modeEnv = std::getenv("OMFEM_GAP_REFINE");
        const std::string mode = modeEnv ? modeEnv : "corner";
        if (mode == "face") {
            // FACE MODE (ABT #1305, the measured reference): a box over the whole column plus a size graded
            // from the gap's whole FACE curves. Correct, but it carries the 0.83 um size along the full
            // column width -- 4.08M unknowns on multicolumn E42 (ABT #1292).
            sizer.add_box(gt, x0, g.yc - yb, -1.0, x1, g.yc + yb, 1.0, 6.0 * gt);
            // GRADED FRINGING FIELD (ABT #1305, measured 2026-09-21). The box above is only
            // max(len, gt) + 4 gt tall -- +-8 um for a 5 um gap -- and gmsh evaluates size at VERTICES,
            // so a triangle larger than the band can straddle it with no vertex inside and never see it.
            // On PQ 20/16 ONE 203 um element covered everything from 376 um to 76 um off the outer-gap
            // mouth (h/r = 2.7 at the failing Ampere loop, which read -0.62 A where it must read 0).
            // Near a slot mouth the field goes like MMF/(pi r) (mkf-57), so P1 needs h well below r out
            // to wherever copper and flux paths sit: a size that GROWS WITH DISTANCE from the gap's own
            // faces, h(d) ~ gt + d/N, defined everywhere. gmsh's Threshold ramps linearly from gt at the
            // faces to the air ceiling at DistMax, so DistMax = N * air_target gives slope 1/N. N is the
            // number of elements per distance-to-mouth, OMFEM_GAP_GRADE_N; its default is the smallest
            // value that passes the Ampere loops (see the commit that set it), not a chosen number.
            sizer.add_distance_refinement(faces, gt, gt, gradeN * opt.air_target, 400);
            std::fprintf(stderr, "[mesh]   FACE mode: graded from %zu gap-face curve(s): h = %.3gum + d/%.3g up to %.3gmm\n",
                         faces.size(), gt * 1e6, gradeN, opt.air_target * 1e3);
        } else if (mode == "corner") {
            // CORNER MODE (Alf 2026-09-21: "a small halo around small parts instead of meshing the whole
            // design with a fine mesh"). The field in a thin gap is uniform except within a few gap lengths
            // of its MOUTHS; the singularities are the mouth corners. So the fine size is a POINT halo
            // around each corner, graded at the same measured slope 1/N, and the gap faces carry only a
            // along-gap cap of OMFEM_GAP_FACE_CAP x the gap length (the strip between them is then a few
            // layers of elongated elements, which represent its uniform field). No box over the column.
            // The mouth corners, from the GEOMETRY: a column's MAS width is not where its leg meets the section
            // (a PQ leg is not a rectangle -- PQ 20/16's outer leg starts at 9.176 mm, not the 9.000 mm of
            // x - width/2). A mouth corner is an endpoint of the gap's curves where the boundary TURNS: the
            // end of a face that opens into the window (one gap curve there), or the meeting of a face and a
            // short mouth line where the section closes the gap as its own face (two non-collinear curves).
            // A point that merely splits a straight face is collinear and is not a corner. A point on the
            // section's outer boundary is not a mouth either (a symmetry cut is not a singularity).
            std::map<int, std::vector<std::array<double, 2>>> ends;   // point -> unit directions of its gap curves
            for (double fc : faces) {
                gmsh::vectorpair pts; gmsh::model::getBoundary({{1, static_cast<int>(fc)}}, pts, false, false, false);
                if (pts.size() != 2) continue;                        // a closed curve has no corner of its own
                std::vector<double> pa, pb;
                gmsh::model::getValue(0, std::abs(pts[0].second), {}, pa);
                gmsh::model::getValue(0, std::abs(pts[1].second), {}, pb);
                const double dx = pb[0] - pa[0], dy = pb[1] - pa[1], dn = std::hypot(dx, dy);
                if (!(dn > 0.0)) continue;
                ends[std::abs(pts[0].second)].push_back({dx / dn, dy / dn});
                ends[std::abs(pts[1].second)].push_back({-dx / dn, -dy / dn});
            }
            double sx0, sy0, sz0, sx1, sy1, sz1;
            gmsh::model::getBoundingBox(-1, -1, sx0, sy0, sz0, sx1, sy1, sz1);
            std::vector<double> corners;
            for (const auto& [pt, dirs] : ends) {
                bool turns = dirs.size() == 1;
                for (size_t i = 0; i < dirs.size() && !turns; ++i)
                    for (size_t k = i + 1; k < dirs.size() && !turns; ++k)
                        turns = dirs[i][0] * dirs[k][0] + dirs[i][1] * dirs[k][1] > -0.999;   // not opposite = a turn
                if (!turns) continue;
                std::vector<double> xyz; gmsh::model::getValue(0, pt, {}, xyz);
                if (std::fabs(xyz[0] - sx0) <= tol || std::fabs(xyz[0] - sx1) <= tol ||
                    std::fabs(xyz[1] - sy0) <= tol || std::fabs(xyz[1] - sy1) <= tol) continue;
                corners.push_back(static_cast<double>(pt));
                std::fprintf(stderr, "[mesh]     mouth corner at (%.4f, %.4f) mm\n", xyz[0] * 1e3, xyz[1] * 1e3);
            }
            if (corners.empty())
                throw std::runtime_error("mesh_from_mas: the gap at x=" + std::to_string(g.xc * 1e3) + " mm, y=" +
                                         std::to_string(g.yc * 1e3) + " mm has no mouth corner among the ends of its " +
                                         std::to_string(faces.size()) + " face curve(s)");
            sizer.add_point_refinement(corners, gt, gt, gradeN * opt.air_target);
            // Along-gap size in gap lengths. 2 is measured, not chosen: caps of 10 / 5 / 2 / 1 / 0.5 / 0.17 move
            // L by < 0.001 % and P_cu by < 0.002 % (Five_Turns); 2 is where the gap-adjacent loops settle.
            const char* capEnv = std::getenv("OMFEM_GAP_FACE_CAP");
            const double cap = (capEnv ? std::atof(capEnv) : kGapFaceCapLengths) * g.len;
            if (!(cap >= gt)) throw std::runtime_error("OMFEM_GAP_FACE_CAP must give a size >= the corner size");
            sizer.add_distance_refinement(faces, cap, cap, gradeN * opt.air_target, 400);
            std::fprintf(stderr, "[mesh]   CORNER mode: %zu mouth corner(s) h = %.3gum + d/%.3g; %zu gap face(s) capped at %.3gum\n",
                         corners.size(), gt * 1e6, gradeN, faces.size(), cap * 1e6);
        } else {
            throw std::runtime_error("OMFEM_GAP_REFINE must be 'face' or 'corner', got '" + mode + "'");
        }
        cmin = std::min(cmin, gt);
    }
    sizer.finalize(cmin);   // MeshSizeMin=cmin, MeshSizeMax=air_target
    mvb::mesh::generate_robust(2, mvb::mesh::algorithm_chain_2d(false));
    mvb::mesh::optimize_mesh(2);
    mvb::mesh::validate_regions_nonempty(2);   // throw if any face region is empty
    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
    write_mesh_atomically(opt.out_msh);
    gmsh::model::remove();
    if (own_init) gmsh::finalize();
    if (!std::getenv("OMFEM_KEEP_STEP"))
        std::remove(step_path.c_str());                // drop the temporary STEP (keep it for review)
    return opt.out_msh;
}

}  // namespace mvb::mesh
