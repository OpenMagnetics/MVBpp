#include "mvb/mesh/SizeField.h"
#include <limits>

#include <gmsh.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>

#include <sys/wait.h>
#include <unistd.h>

namespace mvb::mesh {

namespace gf = gmsh::model::mesh::field;
namespace go = gmsh::option;

double skin_depth(double f_hz, double sigma, double mur) {
    if (!(f_hz > 0.0))  throw std::invalid_argument("skin_depth: frequency must be > 0");
    if (!(sigma > 0.0)) throw std::invalid_argument("skin_depth: conductivity must be > 0");
    if (!(mur > 0.0))   throw std::invalid_argument("skin_depth: mu_r must be > 0");
    const double mu0 = 4.0e-7 * std::numbers::pi;
    return 1.0 / std::sqrt(std::numbers::pi * f_hz * mu0 * mur * sigma);
}

SizeFieldBuilder::SizeFieldBuilder(int dim, double air_target, const SizingPolicy& policy)
    : dim_(dim), air_target_(air_target), policy_(policy) {
    if (dim_ != 2 && dim_ != 3) throw std::invalid_argument("SizeFieldBuilder: dim must be 2 or 3");
    if (!(air_target_ > 0.0))   throw std::invalid_argument("SizeFieldBuilder: air_target must be > 0");
}

void SizeFieldBuilder::add_distance_refinement(const std::vector<double>& entities, double size_min,
                                               double dist_min, double dist_max, int sampling) {
    if (entities.empty() || !(size_min > 0.0)) return;   // nothing to refine toward
    if (dist_min <= 0.0) dist_min = size_min;
    if (dist_max <= 0.0) dist_max = policy_.dist_band * size_min;
    const int df = gf::add("Distance");
    gf::setNumbers(df, dim_ == 2 ? "CurvesList" : "SurfacesList", entities);
    gf::setNumber(df, "Sampling", static_cast<double>(sampling));
    const int tf = gf::add("Threshold");
    gf::setNumber(tf, "InField", df);
    gf::setNumber(tf, "SizeMin", size_min);
    gf::setNumber(tf, "SizeMax", air_target_);
    gf::setNumber(tf, "DistMin", dist_min);
    gf::setNumber(tf, "DistMax", dist_max);
    gf::setNumber(tf, "StopAtDistMax", 1);   // far field defers to the ceiling, not this threshold
    fields_.push_back(tf);
    finest_ = std::min(finest_, size_min);
}

void SizeFieldBuilder::add_point_refinement(const std::vector<double>& points, double size_min,
                                            double dist_min, double dist_max) {
    if (points.empty() || !(size_min > 0.0)) return;
    if (!(dist_max > dist_min)) throw std::invalid_argument("add_point_refinement: dist_max must exceed dist_min");
    const int df = gf::add("Distance");
    gf::setNumbers(df, "PointsList", points);
    const int tf = gf::add("Threshold");
    gf::setNumber(tf, "InField", df);
    gf::setNumber(tf, "SizeMin", size_min);
    gf::setNumber(tf, "SizeMax", air_target_);
    gf::setNumber(tf, "DistMin", dist_min);
    gf::setNumber(tf, "DistMax", dist_max);
    gf::setNumber(tf, "StopAtDistMax", 1);
    fields_.push_back(tf);
    finest_ = std::min(finest_, size_min);
}

void SizeFieldBuilder::add_region_size(const std::vector<double>& regions, double size_in) {
    if (regions.empty() || !(size_in > 0.0)) return;
    const int fc = gf::add("Constant");
    gf::setNumbers(fc, dim_ == 2 ? "SurfacesList" : "VolumesList", regions);
    gf::setNumber(fc, "VIn", size_in);
    gf::setNumber(fc, "VOut", air_target_);
    gf::setNumber(fc, "IncludeBoundary", 1);
    fields_.push_back(fc);
    finest_ = std::min(finest_, size_in);
}

void SizeFieldBuilder::add_box(double size_in, double x0, double y0, double z0,
                               double x1, double y1, double z1, double thickness) {
    if (!(size_in > 0.0)) return;
    const int fb = gf::add("Box");
    gf::setNumber(fb, "VIn", size_in);
    gf::setNumber(fb, "VOut", air_target_);
    gf::setNumber(fb, "XMin", x0); gf::setNumber(fb, "XMax", x1);
    gf::setNumber(fb, "YMin", y0); gf::setNumber(fb, "YMax", y1);
    gf::setNumber(fb, "ZMin", z0); gf::setNumber(fb, "ZMax", z1);
    gf::setNumber(fb, "Thickness", thickness);
    fields_.push_back(fb);
    finest_ = std::min(finest_, size_in);
}

void SizeFieldBuilder::add_cylinder(double size_in, double radius, double thickness,
                                    double cx, double cy, double cz, double ax, double ay, double az) {
    if (!(size_in > 0.0) || !(radius > 0.0)) return;
    const int fc = gf::add("Cylinder");
    gf::setNumber(fc, "VIn", size_in);
    gf::setNumber(fc, "VOut", air_target_);
    gf::setNumber(fc, "Radius", radius);
    gf::setNumber(fc, "XCenter", cx); gf::setNumber(fc, "YCenter", cy); gf::setNumber(fc, "ZCenter", cz);
    gf::setNumber(fc, "XAxis", ax); gf::setNumber(fc, "YAxis", ay); gf::setNumber(fc, "ZAxis", az);
    (void)thickness;   // gmsh Cylinder has a hard VIn/VOut step; the surrounding Min grades it
    fields_.push_back(fc);
    finest_ = std::min(finest_, size_in);
}

void SizeFieldBuilder::add_ring_bands(char axis, const std::vector<RingBand>& bands) {
    if (bands.empty()) return;
    // Cylindrical coordinates about the revolution axis, written in gmsh's MathEval variables.
    std::string rad, along;
    switch (axis) {
        case 'X': case 'x': rad = "sqrt(y^2+z^2)"; along = "x"; break;
        case 'Y': case 'y': rad = "sqrt(x^2+z^2)"; along = "y"; break;
        case 'Z': case 'z': rad = "sqrt(x^2+y^2)"; along = "z"; break;
        default: throw std::invalid_argument(
            "SizeFieldBuilder::add_ring_bands: axis must be X, Y or Z");
    }
    auto F = [](double v) { char s[64]; std::snprintf(s, sizeof s, "%.10g", v); return std::string(s); };
    std::string expr;
    std::size_t used = 0;
    for (const auto& b : bands) {
        if (!(b.size_in > 0.0) || !(b.tube > 0.0) || !(b.thickness > 0.0) || !(b.radius > 0.0))
            throw std::invalid_argument("SizeFieldBuilder::add_ring_bands: size_in, radius, tube "
                                        "and thickness must all be > 0");
        if (b.size_in >= air_target_) continue;   // coarser than the ceiling: refines nothing
        // Distance from the band circle, in the (radius, along-axis) half-plane.
        const std::string d = "sqrt((" + rad + "-(" + F(b.radius) + "))^2+((" + along + ")-(" +
                              F(b.axis_pos) + "))^2)";
        const std::string ramp = "min(1,max(0,(" + d + "-(" + F(b.tube) + "))/(" +
                                 F(b.thickness) + ")))";
        const std::string term = "(" + F(b.size_in) + "+(" + F(air_target_ - b.size_in) + ")*" +
                                 ramp + ")";
        expr = expr.empty() ? term : ("min(" + expr + "," + term + ")");
        finest_ = std::min(finest_, b.size_in);
        ++used;
    }
    if (!used) return;
    const int fe = gf::add("MathEval");
    gf::setString(fe, "F", expr);
    fields_.push_back(fe);
}

void SizeFieldBuilder::add_ceiling(double size) {
    if (!(size > 0.0)) return;
    const int fc = gf::add("MathEval");
    char buf[64]; std::snprintf(buf, sizeof(buf), "%.10g", size);
    gf::setString(fc, "F", buf);
    fields_.push_back(fc);
}

int SizeFieldBuilder::finalize(double min_floor) {
    if (fields_.empty())
        throw std::runtime_error("SizeFieldBuilder::finalize: no size sub-fields were added");
    std::vector<double> tags(fields_.begin(), fields_.end());
    const int fmin = gf::add("Min");
    gf::setNumbers(fmin, "FieldsList", tags);
    gf::setAsBackgroundMesh(fmin);
    // The Min field is now authoritative: silence every competing size source so it isn't
    // overridden by point sizes / boundary extension.
    go::setNumber("Mesh.MeshSizeFromPoints", 0);
    go::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    // Curvature: cheap and robust in 2D (clamped below); in 3D it inflates count and fights HXT.
    go::setNumber("Mesh.MeshSizeFromCurvature", dim_ == 2 ? static_cast<double>(policy_.N_curv) : 0.0);
    // Derived floor: half the finest feature size, so gmsh never coarsens below what a field asks.
    if (min_floor <= 0.0) min_floor = (finest_ < 1e29) ? 0.5 * finest_ : 0.0;
    if (min_floor > 0.0) go::setNumber("Mesh.MeshSizeMin", min_floor);
    go::setNumber("Mesh.MeshSizeMax", air_target_);
    return fmin;
}

std::vector<int> algorithm_chain_3d() {
    if (const char* e = std::getenv("OMFEM_ALGO3D")) return { std::atoi(e) };  // debug override
    // OMFEM_ALGO3D_CHAIN="1,10" sets the whole chain and its ORDER, keeping a fallback. The default
    // below puts HXT first for its failure SHAPE, which is the right call in general -- but on the
    // real-winding corpus HXT wins rarely and crashes often (2026-09-20: of six meshes, FIVE were
    // produced by Delaunay and one by HXT, at the cost of a segfault and a fallback on most designs),
    // so a run over that corpus is better off asking Delaunay first and keeping HXT as the escape.
    if (const char* c = std::getenv("OMFEM_ALGO3D_CHAIN")) {
        std::vector<int> chain; std::string tok;
        for (std::istringstream is(c); std::getline(is, tok, ','); )
            if (!tok.empty()) chain.push_back(std::atoi(tok.c_str()));
        if (chain.empty())
            throw std::runtime_error("OMFEM_ALGO3D_CHAIN is set but empty -- give a comma-separated "
                                     "list of gmsh 3D algorithm ids, e.g. \"1,10\"");
        return chain;
    }
    // HXT first, Delaunay as fallback. The order matters because of their FAILURE SHAPES, not
    // their success rates: when HXT can't mesh, it fails in seconds ("HXT 3D mesh failed"), so
    // the fallback still gets its turn; when tetgen-Delaunay can't cope it can STALL for hours
    // in boundary recovery, starving every algorithm after it (measured on e138 at the 88 um
    // auto target: algo 1 burned a full 3600 s timeout flat at 0.5 GB; HXT then meshed the same
    // model in 458 s / 202k tets). Never 7 (MMG3D) for first pass.
    return { 10, 1 };
}

std::vector<int> algorithm_chain_2d(bool boundary_layer) {
    if (const char* e = std::getenv("OMFEM_ALGO2D")) return { std::atoi(e) };
    if (boundary_layer) return { 7, 6 };   // BAMG honours the anisotropic BL metric
    return { 6, 5 };                        // Frontal-Delaunay -> Delaunay
}

static std::size_t element_count(int dim) {
    std::vector<int> et; std::vector<std::vector<std::size_t>> en, ev;
    gmsh::model::mesh::getElements(et, en, ev, dim);
    std::size_t n = 0;
    for (std::size_t i = 0; i < et.size(); ++i) n += en[i].size();
    return n;
}

// Algorithms whose FAILURE MODE can kill the process instead of throwing. HXT (10) corrupts
// its heap on some surface triangulations -- observed as SIGSEGV, "free(): invalid",
// "double free", "malloc():" and "Fatal glibc error" on the SAME design across runs -- so the
// fallback below, which only catches C++ exceptions, never got its turn: the process was gone.
// The corpus compensated by restarting the whole mesher at a COARSER conductor target, i.e. it
// answered a crash with a mesh nobody asked for.
// A 3D meshing algorithm can kill the process instead of throwing: gmsh's HXT (10) corrupts its
// heap on some surface triangulations -- SIGSEGV, "free(): invalid", "double free", "malloc():",
// "Fatal glibc error" on the SAME design across runs. generate_robust's fallback only catches C++
// exceptions, so the process was gone before the fallback could run and the corpus compensated by
// re-meshing at a different target.
//
// So the PARENT never meshes. It forks a child per algorithm and supervises: the child returns
// normally and goes on to finish the whole run (it is a full copy of this process), so on success
// the parent has nothing left to do and just propagates the child's exit status. Only a death BY
// SIGNAL brings the parent back for the next algorithm.
//
// Why the parent must not mesh the fallback itself, which is how this was first written: after an
// isolated HXT crash the parent DID mesh successfully with Delaunay and then died silently inside
// the optimiser (buck_inductor_complete, twice, 2026-09-19) -- while the same Delaunay mesh at the
// same target in a FRESH process completed and wrote its file. Whatever the forked crash leaves
// behind, a parent that has never meshed does not carry it. Keeping the parent a pure supervisor
// costs one fork per attempt and no CAD rework, because fork copies the built model.
// Returns: 0 = this attempt crashed, try the next. 1 = this IS the child, carry on with the run.
static int fork_attempt(int dim, int algo, const std::string& algo_opt, bool escalate) {
    std::fflush(nullptr);                 // no buffered output duplicated into the child
    const pid_t pid = fork();
    if (pid < 0)
        throw std::runtime_error(std::string("generate_robust: fork failed: ") + std::strerror(errno));
    if (pid == 0) {
        // TEST HOOK: OMFEM_MESH_FAULT_INJECT=<algo> makes this child die by SIGSEGV on that
        // algorithm. The crash is nondeterministic in the wild (HXT meshed the same input
        // cleanly on one run and segfaulted on the next), so without this the RECOVERY path
        // cannot be exercised on demand -- and an untested recovery path is the one that fails
        // when it is finally needed.
        if (const char* fi = std::getenv("OMFEM_MESH_FAULT_INJECT")) {
            if (std::atoi(fi) == algo) {
                std::fprintf(stderr, "[mesh] FAULT INJECTION: killing this child on algo %d\n", algo);
                std::fflush(nullptr);
                std::raise(SIGSEGV);
            }
        }
        try {
            go::setNumber(algo_opt, algo);
            if (escalate) go::setNumber("Mesh.ToleranceInitialDelaunay", 1e-9);
            gmsh::model::mesh::generate(dim);      // may die here; the parent survives
            if (element_count(dim) == 0) _exit(70);          // empty mesh: distinct from a crash
        } catch (const std::exception& e) {
            std::fprintf(stderr, "[mesh] algo %d threw: %s\n", algo, e.what());
            _exit(71);                                       // a refusal, not a crash
        }
        return 1;                               // child carries on and finishes the run
    }
    int status = 0;
    while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
    if (WIFSIGNALED(status)) {
        std::fprintf(stderr, "[mesh] algo %d died on signal %d (%s) [recovered: next algorithm, "
                             "same target]\n", algo, WTERMSIG(status), strsignal(WTERMSIG(status)));
        return 0;
    }
    if (WIFEXITED(status) && (WEXITSTATUS(status) == 70 || WEXITSTATUS(status) == 71)) {
        std::fprintf(stderr, "[mesh] algo %d did not mesh (exit %d) [recovered: next algorithm]\n",
                     algo, WEXITSTATUS(status));
        return 0;
    }
    // The child WAS the run: it meshed, wrote the output and exited. Do not redo any of it.
    _exit(WIFEXITED(status) ? WEXITSTATUS(status) : 1);
}

int generate_robust(int dim, const std::vector<int>& chain) {
    const std::string algo_opt = (dim == 3) ? "Mesh.Algorithm3D" : "Mesh.Algorithm";
    // ISOLATION IS OFF BY DEFAULT -- it CAUSED more crashes than it caught (2026-09-20).
    //
    // fork() here happens AFTER gmsh has gone multithreaded (OMFEM_GMSH_THREADS=12, plus OCC's own
    // pool), so the child inherits malloc arenas and thread-pool state left by threads that do not
    // exist in it, and then runs the whole remaining job in that state. That is undefined
    // behaviour, and it showed: with isolation on, BOTH algorithms died by signal on the same
    // design at the same target -- algo 1 SIGSEGV and algo 10 "corrupted double-linked list" --
    // 6 signal deaths each across the pass. Two entirely separate code paths corrupting the heap on
    // identical input is not a property of either algorithm.
    // The control settles it: 05_pfc_inductor_t4020_hf60 at the identical target with
    // OMFEM_MESH_ISOLATE=0 meshed cleanly, algo 1, 2474451 tets in 631 s, no crash, while the
    // corpus with isolation on crashed on that same design repeatedly.
    // OMFEM_MESH_ISOLATE=1 re-enables it. Making it safe means forking BEFORE any thread is
    // created (a supervisor at process start that re-runs the whole job), not here.
    const bool isolate = dim == 3 && std::getenv("OMFEM_MESH_ISOLATE")
                      && std::atoi(std::getenv("OMFEM_MESH_ISOLATE")) != 0;
    std::string failures;
    for (std::size_t i = 0; i < chain.size(); ++i) {
        const int algo = chain[i];
        try {
            if (isolate) {
                if (fork_attempt(dim, algo, algo_opt, i > 0) == 0) {
                    failures += "algo " + std::to_string(algo) + ": crashed or refused (isolated); ";
                    continue;
                }
                std::fprintf(stderr, "[mesh] %dD meshed by algo %d%s\n", dim, algo,
                             i ? " (after an earlier algorithm failed)" : "");
                return algo;                    // this IS the child, carrying on with the run
            }
            go::setNumber(algo_opt, algo);
            if (i > 0) go::setNumber("Mesh.ToleranceInitialDelaunay", 1e-9);  // escalate on retry
            gmsh::model::mesh::generate(dim);
            if (element_count(dim) > 0) {
                std::fprintf(stderr, "[mesh] %dD meshed by algo %d%s\n", dim, algo,
                             i ? " (after an earlier algorithm failed)" : "");
                return algo;
            }
            failures += "algo " + std::to_string(algo) + ": empty mesh; ";
        } catch (const std::exception& e) {
            failures += "algo " + std::to_string(algo) + ": " + e.what() + "; ";
            try { gmsh::model::mesh::clear(); } catch (...) {}
        }
    }
    throw std::runtime_error("generate_robust(" + std::to_string(dim) + "D): all algorithms failed: " + failures);
}

void validate_regions_nonempty(int topdim) {
    gmsh::vectorpair groups;
    gmsh::model::getPhysicalGroups(groups, topdim);
    std::string empty;
    for (const auto& [dim, tag] : groups) {
        std::string name; gmsh::model::getPhysicalName(dim, tag, name);
        std::vector<int> ents; gmsh::model::getEntitiesForPhysicalGroup(dim, tag, ents);
        std::size_t n = 0;
        for (int e : ents) {
            std::vector<int> et; std::vector<std::vector<std::size_t>> en, ev;
            gmsh::model::mesh::getElements(et, en, ev, dim, e);
            for (std::size_t i = 0; i < et.size(); ++i) n += en[i].size();
        }
        if (n == 0) empty += (empty.empty() ? "" : ", ") + (name.empty() ? ("tag " + std::to_string(tag)) : name);
    }
    if (!empty.empty())
        throw std::runtime_error("mesh validation: region(s) have ZERO elements (swallowed by the "
                                 "fragment or mis-classified): " + empty);
}

void optimize_mesh(int topdim) {
    try {
        gmsh::option::setNumber("Mesh.Optimize", 1);
        if (topdim == 3) {
            gmsh::option::setNumber("Mesh.OptimizeNetgen", 1);
            gmsh::model::mesh::optimize("Netgen");
        }
        gmsh::model::mesh::optimize("");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[mesh] optimize warning: %s\n", e.what());
    }
}


ThinSolidRefinement add_thin_solid_refinement(SizeFieldBuilder& builder, double thin_cut,
                                              int elements_across, double min_aspect) {
    if (!(thin_cut > 0.0))
        throw std::invalid_argument("add_thin_solid_refinement: thin_cut must be > 0");
    if (elements_across < 1)
        throw std::invalid_argument("add_thin_solid_refinement: elements_across must be >= 1");
    ThinSolidRefinement out;
    gmsh::vectorpair vols;
    gmsh::model::occ::getEntities(vols, 3);
    if (vols.empty()) gmsh::model::getEntities(vols, 3);
    std::vector<double> faces;
    double thinnest = std::numeric_limits<double>::max();
    for (const auto& dt : vols) {
        double x0, y0, z0, x1, y1, z1;
        try {
            gmsh::model::getBoundingBox(dt.first, dt.second, x0, y0, z0, x1, y1, z1);
        } catch (const std::exception&) {
            continue;
        }
        const double bboxMin = std::min(std::min(x1 - x0, y1 - y0), z1 - z0);
        const double dmax = std::max(std::max(x1 - x0, y1 - y0), z1 - z0);
        gmsh::vectorpair bnd;
        gmsh::model::getBoundary({dt}, bnd, /*combined=*/false, /*oriented=*/false);
        // A BOUNDING BOX CANNOT SEE A RING'S WALL. The bbox is the right thickness for a wire
        // piece, and blind for any solid that wraps around something: a foil band 25 um thick on
        // a 10 mm radius has a bbox 21 mm on its smallest side, so the detector reported "no
        // solid thinner than 61.578 mm; nothing to refine" and left every band one element thick
        // (ABT #970 spike, 2026-09-04). The same blindness applies to a toroid's conductor rings.
        //
        // 4V/S is thickness measured on the solid itself, not on a box around it: for a slab of
        // thickness t the surface is ~2A and 4V/S = 2t; for a cylinder of radius r it is exactly
        // the diameter 2r. Taking the MINIMUM of the two keeps every existing case bit-identical
        // (a wire's bbox already gives its diameter, which is what elements_across means) and
        // makes a wrapped solid visible at twice its wall, which is finite where the bbox was
        // infinite. Exactness for the wrapped case is the caller's elements_across to spend.
        double area = 0.0;
        for (const auto& f : bnd) {
            if (f.first != 2) continue;
            double a = 0.0;
            try {
                gmsh::model::occ::getMass(2, f.second, a);
            } catch (const std::exception&) {
                a = 0.0;
            }
            area += std::abs(a);
        }
        double volume = 0.0;
        try {
            gmsh::model::occ::getMass(3, dt.second, volume);
        } catch (const std::exception&) {
            volume = 0.0;
        }
        volume = std::abs(volume);
        const double hydraulic = area > 0.0 ? 4.0 * volume / area
                                            : std::numeric_limits<double>::max();
        const double dmin = std::min(bboxMin, hydraulic);
        if (!(dmin > 0.0) || dmin > thin_cut) continue;
        if (!(dmax >= min_aspect * dmin)) continue;   // chunky: a core half or a bobbin
        ++out.solids;
        thinnest = std::min(thinnest, dmin);
        for (const auto& f : bnd)
            if (f.first == 2) faces.push_back(static_cast<double>(std::abs(f.second)));
    }
    if (faces.empty()) return out;
    std::sort(faces.begin(), faces.end());
    faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
    out.faces = faces.size();
    out.thinnest = thinnest;
    out.size_min = thinnest / elements_across;
    builder.add_distance_refinement(faces, out.size_min);
    return out;
}
}  // namespace mvb::mesh
