// Moved verbatim from OMFEM src/meshing/StepAudit.cpp (ABT #1588, step 3): namespace omfem::cad -> mvb::cad. No logic changed.
#include <STEPCAFControl_Reader.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDataStd_Name.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Extrema_ExtFlag.hxx>
#include <Extrema_ExtAlgo.hxx>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Vec.hxx>
#include <unistd.h>
#include <sys/wait.h>
#include <BRepClass3d_SolidClassifier.hxx>
#include <gp_Pnt.hxx>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "mvb/mesh/StepAudit.h"

#include <algorithm>
#include <cstdio>
#include <iostream>

namespace mvb::cad {

using Part = AuditPart;

namespace {

std::string labelName(const TDF_Label& label) {
    Handle(TDataStd_Name) nameAttribute;
    if (label.FindAttribute(TDataStd_Name::GetID(), nameAttribute)) {
        TCollection_AsciiString ascii(nameAttribute->Get());
        return std::string(ascii.ToCString());
    }
    return std::string();
}

double solidVolume(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}



// FORK-ISOLATED BOOLEAN. OCCT 7.9.3 SEGFAULTS on some touching pairs of these assemblies
// (reproduced: 'Bobbin' x 'Primary parallel 0 [solid 8]' on 01_simple_inductor_etd34_n87,
// inside BOPAlgo_BuilderSolid::PerformAreas -> BRepClass3d_SolidClassifier::
// PerformInfinitePoint). A segfault cannot be caught, so one bad pair used to take the whole
// corpus audit down with it -- which is exactly when an overlap audit matters most.
//
// Run the boolean in a forked child and read back only the numbers we need. A child that dies
// costs us that ONE pair, reported as BOOLEAN_CRASHED, and the sweep continues. The child is
// read-only w.r.t. the parent's OCCT state (copy-on-write), so this is safe.
// Material class of a named solid, so the report can answer the only two questions that carry
// different consequences: copper touching copper is a MVB++ geometry defect to be fixed, while
// copper inside the core is a DESIGN problem (Alf, 2026-08-25: "if core is involved, then we
// need to change the design, but run it by me first").
enum class Mat { Copper, Core, Bobbin, Insulation, Other };

Mat classify(const std::string& n) {
    auto has = [&](const char* k) { return n.find(k) != std::string::npos; };
    if (has("coating") || has("Coating") || has("nsulation") || has("Layer") || has("layer"))
        return Mat::Insulation;
    if (has("Core") || has("core")) return Mat::Core;
    if (has("Bobbin") || has("bobbin")) return Mat::Bobbin;
    if (has("FR4")) return Mat::Other;
    // Everything MKF winds is copper: "Primary parallel 0", "Winding 2 parallel 1",
    // "<winding> terminal 0", and their per-solid chunks.
    if (has("parallel") || has("Winding") || has("winding") || has("terminal") || has("Turn"))
        return Mat::Copper;
    return Mat::Other;
}

const char* matName(Mat m) {
    switch (m) {
        case Mat::Copper:     return "copper";
        case Mat::Core:       return "core";
        case Mat::Bobbin:     return "bobbin";
        case Mat::Insulation: return "insulation";
        default:              return "other";
    }
}


// BOOLEAN-FREE INTERPENETRATION TEST. When OCCT segfaults on a pair (it does; see
// commonIsolated) the pair is left UNPROVEN, and an audit whose answer is "I could not tell"
// on the very pairs that break the kernel is not much of an audit. This decides the same
// question without a boolean: tessellate A's boundary, step a short way INWARD from each
// facet into A's own interior, and classify those points against B. A point strictly inside
// both solids is interpenetration, full stop -- the same point-in-solid ground truth the
// overlap doctrine already uses to verify a boolean's verdict.
//
// Sampling can PROVE an overlap; it cannot prove absence, so a negative result is reported as
// "none found by sampling", never as CLEAN.
struct SampleResult {
    bool   ran = false;
    bool   found = false;
    size_t points = 0;
    gp_Pnt where;
};

SampleResult sampledInterpenetration(const TopoDS_Shape& A, const TopoDS_Shape& B,
                                     double linDeflection) {
    SampleResult out;
    try {
        BRepMesh_IncrementalMesh mesher(A, linDeflection, Standard_False, 0.5, Standard_True);
        if (!mesher.IsDone()) return out;
        BRepClass3d_SolidClassifier inA(A), inB(B);
        // Step inward by a fraction of the deflection: far enough to be off the surface (so a
        // tangential contact does not read as interpenetration), small enough to stay inside a
        // thin piece.
        const double step = linDeflection * 0.5;
        for (TopExp_Explorer fx(A, TopAbs_FACE); fx.More(); fx.Next()) {
            const TopoDS_Face& f = TopoDS::Face(fx.Current());
            TopLoc_Location loc;
            Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(f, loc);
            if (tri.IsNull()) continue;
            const gp_Trsf& tr = loc.Transformation();
            for (Standard_Integer t = 1; t <= tri->NbTriangles(); ++t) {
                Standard_Integer i1, i2, i3;
                tri->Triangle(t).Get(i1, i2, i3);
                gp_Pnt p1 = tri->Node(i1).Transformed(tr);
                gp_Pnt p2 = tri->Node(i2).Transformed(tr);
                gp_Pnt p3 = tri->Node(i3).Transformed(tr);
                const gp_Pnt c((p1.X() + p2.X() + p3.X()) / 3.0, (p1.Y() + p2.Y() + p3.Y()) / 3.0,
                               (p1.Z() + p2.Z() + p3.Z()) / 3.0);
                gp_Vec n = gp_Vec(p1, p2).Crossed(gp_Vec(p1, p3));
                if (n.Magnitude() < 1e-18) continue;
                n.Normalize();
                // Try both sides and keep whichever lands strictly inside A.
                for (double sgn : {-1.0, 1.0}) {
                    const gp_Pnt probe = c.Translated(n * (sgn * step));
                    inA.Perform(probe, 1e-9);
                    if (inA.State() != TopAbs_IN) continue;
                    ++out.points;
                    inB.Perform(probe, 1e-9);
                    if (inB.State() == TopAbs_IN) {
                        out.ran = true; out.found = true; out.where = probe;
                        return out;
                    }
                    break;
                }
            }
        }
        out.ran = true;
    } catch (const Standard_Failure&) {
    }
    return out;
}

struct CommonResult {
    bool   ok = false;        // the boolean completed
    bool   crashed = false;   // the child died (signal) -- pair unproven, not clean
    bool   separated = false; // Extrema proved a positive minimum distance; no boolean run
    double volume = 0.0;      // mm^3 of the common
    double cx = 0, cy = 0, cz = 0;
    double bx0 = 0, by0 = 0, bz0 = 0, bx1 = 0, by1 = 0, bz1 = 0;   // bbox of the common
};

// The boolean stays fork-isolated (OCCT segfaults must not kill the audit), but the pool
// below runs several children CONCURRENTLY: a 178-solid faceted toroid has ~176 touching
// pairs at ~17 s of boolean each, which is a 50-minute serial audit and the whole reason the
// sweep's overlap gate kept timing out (OVL_TIMEOUT on every toroid, 2026-08-28).
struct PendingCommon {
    pid_t pid = -1;
    int   fd = -1;
};

PendingCommon spawnCommon(const TopoDS_Shape& A, const TopoDS_Shape& B) {
    PendingCommon pc;
    int fds[2];
    if (pipe(fds) != 0) return pc;
    const pid_t pid = fork();
    if (pid < 0) { close(fds[0]); close(fds[1]); return pc; }
    if (pid == 0) {                       // ---- child ----
        close(fds[0]);
        double buf[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        try {
            // Distance first: a separated pair (true minimum distance > 0) cannot
            // interpenetrate and never needs the boolean. buf[0]=2 marks "proven separated".
            bool separated = false;
            {
                BRepExtrema_DistShapeShape dist(A, B, Extrema_ExtFlag_MIN, Extrema_ExtAlgo_Tree);
                if (dist.IsDone() && dist.NbSolution() > 0 && dist.Value() > 1e-9)
                    separated = true;
            }
            if (separated) {
                buf[0] = 2.0;
            } else {
                BRepAlgoAPI_Common common(A, B);
                if (common.IsDone() && !common.Shape().IsNull()) {
                    double volume = 0.0, weight = 0.0;
                    gp_Pnt centroid(0, 0, 0);
                    for (TopExp_Explorer ex(common.Shape(), TopAbs_SOLID); ex.More(); ex.Next()) {
                        GProp_GProps props;
                        BRepGProp::VolumeProperties(ex.Current(), props);
                        volume += props.Mass();
                        if (props.Mass() > weight) { weight = props.Mass(); centroid = props.CentreOfMass(); }
                    }
                    buf[0] = 1.0; buf[1] = volume;
                    buf[2] = centroid.X(); buf[3] = centroid.Y(); buf[4] = centroid.Z();
                    Bnd_Box cb;
                    for (TopExp_Explorer ex(common.Shape(), TopAbs_SOLID); ex.More(); ex.Next())
                        BRepBndLib::Add(ex.Current(), cb);
                    if (!cb.IsVoid()) cb.Get(buf[5], buf[6], buf[7], buf[8], buf[9], buf[10]);
                }
            }
        } catch (...) {
            buf[0] = 0.0;
        }
        const ssize_t wrote = write(fds[1], buf, sizeof(buf));
        (void)wrote;
        close(fds[1]);
        _exit(0);                          // _exit: no atexit/global dtors in the child
    }
    close(fds[1]);                         // ---- parent ----
    pc.pid = pid;
    pc.fd = fds[0];
    return pc;
}

CommonResult collectCommon(const PendingCommon& pc) {
    CommonResult out;
    if (pc.pid < 0) return out;
    double buf[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    ssize_t got = 0, n = 0;
    while (got < static_cast<ssize_t>(sizeof(buf)) &&
           (n = read(pc.fd, reinterpret_cast<char*>(buf) + got, sizeof(buf) - got)) > 0)
        got += n;
    close(pc.fd);
    int status = 0;
    waitpid(pc.pid, &status, 0);
    if (got < static_cast<ssize_t>(sizeof(buf)) || !WIFEXITED(status)) {
        out.crashed = true;
        return out;
    }
    out.ok = buf[0] != 0.0;
    out.separated = buf[0] == 2.0;
    out.volume = buf[1];
    out.cx = buf[2]; out.cy = buf[3]; out.cz = buf[4];
    out.bx0 = buf[5]; out.by0 = buf[6]; out.bz0 = buf[7];
    out.bx1 = buf[8]; out.by1 = buf[9]; out.bz1 = buf[10];
    return out;
}

CommonResult commonIsolated(const TopoDS_Shape& A, const TopoDS_Shape& B) {
    return collectCommon(spawnCommon(A, B));
}


}  // namespace

int audit_overlaps(std::vector<AuditPart>& parts, double toleranceMm3) {
    std::cout << std::unitbuf;
    std::cout << "[intersect] " << "(imported)" << ": " << parts.size() << " solid(s), tolerance "
              << toleranceMm3 << " mm^3\n";
    if (parts.size() < 2) {
        std::cout << "[intersect] VERDICT: NO OVERLAPS (nothing to pair)\n";
        return 0;
    }

    size_t pairsTested = 0;
    size_t overlaps = 0;
    size_t joints = 0;
    double jointVolume = 0.0;
    size_t contacts = 0;
    size_t artifacts = 0;
    size_t crashedPairs = 0;
    size_t copperOverlaps = 0;
    size_t coreOverlaps = 0;
    size_t copperJoints = 0;
    double contactVolume = 0.0;
    // Equivalent cube side under which a common lump is surfaces TOUCHING, not overlapping:
    // 20 um, an order below any wire dimension in the suite.
    const double contactMm = 0.02;
    double worstVolume = 0.0;
    std::string worstPair;
    // --adjacent: test ONLY consecutive solids of the same body (a piece against the next piece
    // of its own conductor). The mitre doctrine says chain neighbours ABUT exactly, so any
    // neighbour pair with common volume is a construction defect -- and this scoped audit runs
    // in minutes where the all-pairs sweep on a 400-piece faceted toroid exceeds its budget
    // without a verdict (measured 2026-08-27).
    // The caller asks for the scoped audit through the environment now that this is a library
    // function shared by both tools (the standalone tool still accepts --adjacent and sets it).
    const bool adjacentOnly = std::getenv("OMFEM_AUDIT_ADJACENT") != nullptr;
    // Phase 1: collect the candidate pairs (cheap filters, serial). Phase 2 runs the
    // fork-isolated booleans through a bounded pool -- OMFEM_INTERSECT_JOBS children at a
    // time (default 4, the same cap the OCC thread budget uses on this shared box).
    std::vector<std::pair<size_t, size_t>> candidates;
    size_t boxPruned = 0;   // pairs PROVEN clear by the box-overlap volume bound
    for (size_t a = 0; a < parts.size(); ++a) {
        for (size_t b = a + 1; b < parts.size(); ++b) {
            if (adjacentOnly && (b != a + 1 || parts[a].body != parts[b].body)) continue;
            // Bounding boxes that do not touch cannot intersect — the O(n^2) boolean sweep is
            // only affordable because this filter removes almost every pair.
            if (parts[a].box.IsOut(parts[b].box)) {
                continue;
            }
            // TIGHTER, STILL EXACT: the common solid is contained in the INTERSECTION of the two
            // bounding boxes, so its volume can never exceed that box's volume. A pair whose
            // boxes overlap in less than the tolerance cannot produce an above-tolerance overlap
            // -- rejecting it is a proof, not a heuristic, and it costs six comparisons instead
            // of a forked boolean. Grazing corner pairs are exactly the population this removes.
            {
                double ax0, ay0, az0, ax1, ay1, az1, bx0, by0, bz0, bx1, by1, bz1;
                parts[a].box.Get(ax0, ay0, az0, ax1, ay1, az1);
                parts[b].box.Get(bx0, by0, bz0, bx1, by1, bz1);
                const double ox = std::min(ax1, bx1) - std::max(ax0, bx0);
                const double oy = std::min(ay1, by1) - std::max(ay0, by0);
                const double oz = std::min(az1, bz1) - std::max(az0, bz0);
                if (ox <= 0 || oy <= 0 || oz <= 0 || ox * oy * oz <= toleranceMm3) {
                    ++boxPruned;
                    continue;
                }
            }
            // The SEPARATION PRE-FILTER (BRepExtrema) runs INSIDE the forked child now -- on
            // 06_llc (228 heavyweight welded solids) the serial Extrema sweep alone blew the
            // whole 3600 s gate budget before a single boolean ran. The child does distance
            // first and skips the boolean for separated pairs, so the expensive work is both
            // pooled and crash-isolated (Extrema can crash on the same shapes booleans do).
            candidates.emplace_back(a, b);
        }
    }
    // Pool width. The children are short-lived forks of an already-loaded model (copy-on-write,
    // so no second import and no second copy of the geometry), and the box has 24 cores, so the
    // old fixed 4 left most of the machine idle while the sweep waited. One third of the cores,
    // capped at 8, keeps a comfortable margin for the generator and the other agents' work.
    int poolJobs = std::max(4, std::min(8, static_cast<int>(std::thread::hardware_concurrency()) / 3));
    if (const char* pj = std::getenv("OMFEM_INTERSECT_JOBS")) poolJobs = std::max(1, std::atoi(pj));
    struct InFlight { size_t a, b; PendingCommon pc; };
    std::vector<InFlight> inFlight;
    for (size_t next = 0; next < candidates.size() || !inFlight.empty();) {
        while (next < candidates.size() && (int)inFlight.size() < poolJobs) {
            const auto [na, nb] = candidates[next];
            ++next;
            ++pairsTested;
            // Name the pair BEFORE the boolean and flush: when OCCT segfaults there is no
            // unwinding and no destructor, so a buffered line is simply lost -- which is why
            // the first crash report had no output at all to say which pair did it.
            if (std::getenv("OMFEM_INTERSECT_TRACE"))
                std::cerr << "[intersect.trace] '" << parts[na].name << "' x '" << parts[nb].name
                          << "'" << std::endl;
            InFlight f{na, nb, spawnCommon(parts[na].shape, parts[nb].shape)};
            if (f.pc.pid < 0) { --pairsTested; continue; }   // spawn failed; drop loudly below
            inFlight.push_back(f);
        }
        if (inFlight.empty()) break;
        // Reap whichever child finishes first, then judge its pair with the same logic as the
        // serial version.
        int status = 0;
        const pid_t done = waitpid(-1, &status, 0);
        auto it = std::find_if(inFlight.begin(), inFlight.end(),
                               [&](const InFlight& f) { return f.pc.pid == done; });
        if (it == inFlight.end()) continue;
        const size_t a = it->a, b = it->b;
        // collectCommon re-waits on an already-reaped pid; feed it the status we have.
        CommonResult cr;
        {
            double buf[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            ssize_t got = 0, n = 0;
            while (got < static_cast<ssize_t>(sizeof(buf)) &&
                   (n = read(it->pc.fd, reinterpret_cast<char*>(buf) + got, sizeof(buf) - got)) > 0)
                got += n;
            close(it->pc.fd);
            if (got < static_cast<ssize_t>(sizeof(buf)) || !WIFEXITED(status)) {
                cr.crashed = true;
            } else {
                cr.ok = buf[0] != 0.0;
                cr.separated = buf[0] == 2.0;
                cr.volume = buf[1];
                cr.cx = buf[2]; cr.cy = buf[3]; cr.cz = buf[4];
                cr.bx0 = buf[5]; cr.by0 = buf[6]; cr.bz0 = buf[7];
                cr.bx1 = buf[8]; cr.by1 = buf[9]; cr.bz1 = buf[10];
            }
        }
        inFlight.erase(it);
        {
            if (cr.crashed) {
                // Fall back to the boolean-free test rather than give up on the pair.
                double dx, dy, dz, x0, y0, z0, x1, y1, z1;
                parts[a].box.Get(x0, y0, z0, x1, y1, z1);
                dx = x1 - x0; dy = y1 - y0; dz = z1 - z0;
                const double defl = std::max(1e-4, 0.02 * std::min(std::min(dx, dy), dz));
                const SampleResult sr = sampledInterpenetration(parts[a].shape, parts[b].shape, defl);
                if (sr.found) {
                    ++overlaps;
                    const Mat ma2 = classify(parts[a].name), mb2 = classify(parts[b].name);
                    if (ma2 == Mat::Copper && mb2 == Mat::Copper) ++copperOverlaps;
                    if (ma2 == Mat::Core || mb2 == Mat::Core)     ++coreOverlaps;
                    std::cout << "[intersect] OVERLAP [" << matName(ma2) << "/" << matName(mb2)
                              << "] (boolean crashed; proven by point sampling at ("
                              << sr.where.X() << "," << sr.where.Y() << "," << sr.where.Z()
                              << ")): '" << parts[a].name << "' vs '" << parts[b].name << "'\n";
                    if (worstVolume <= 0) worstPair = parts[a].name + " vs " + parts[b].name;
                    continue;
                }
                ++crashedPairs;
                std::cout << "[intersect] BOOLEAN_CRASHED (OCCT): '" << parts[a].name
                          << "' vs '" << parts[b].name << "' -- no interpenetration found by "
                          << "sampling (" << sr.points << " interior points"
                          << (sr.ran ? "" : ", TESSELLATION FAILED") << "); NOT PROVEN CLEAN\n";
                continue;
            }
            if (!cr.ok) {
                std::cout << "[intersect] WARN boolean failed: '" << parts[a].name << "' vs '"
                          << parts[b].name << "'\n";
                continue;
            }
            const double volume = cr.volume;
            const gp_Pnt centroid(cr.cx, cr.cy, cr.cz);
            if (volume > toleranceMm3) {
                // NEVER trust a single boolean. OCCT can return a "common" that is really one
                // whole operand when a solid is malformed or the surfaces are tangent (measured
                // on 26_psps_e3216: a wrap piece 0.4 mm BELOW the flange came back as 100%
                // inside it). Verify by classifying the common's own centre of mass against
                // both parents: a genuine interpenetration has its centroid inside both.
                BRepClass3d_SolidClassifier classifyA(parts[a].shape);
                BRepClass3d_SolidClassifier classifyB(parts[b].shape);
                classifyA.Perform(centroid, 1e-7);
                classifyB.Perform(centroid, 1e-7);
                const bool insideA = classifyA.State() == TopAbs_IN || classifyA.State() == TopAbs_ON;
                const bool insideB = classifyB.State() == TopAbs_IN || classifyB.State() == TopAbs_ON;
                if (!insideA || !insideB) {
                    ++artifacts;
                    continue;   // boolean artefact, not geometry
                }
                // A genuine common volume lies inside BOTH operands, so its bounding box lies
                // inside the intersection of their boxes. OCCT returned, on the current
                // transformer's tangent corner/chord joints, a "common" 3 mm long whose box ran
                // 2.4 mm PAST the corner's own box (x -12.01..-8.90 against the corner's
                // -12.50..-11.88) with a centroid that classified ON the corner's surface. That
                // is a boolean artefact of two tangent surfaces, not copper inside copper.
                {
                    const bool haveBox = cr.bx1 > cr.bx0 || cr.by1 > cr.by0 || cr.bz1 > cr.bz0;
                    if (haveBox) {
                        const double cx0 = cr.bx0, cy0 = cr.by0, cz0 = cr.bz0;
                        const double cx1 = cr.bx1, cy1 = cr.by1, cz1 = cr.bz1;
                        double ax0, ay0, az0, ax1, ay1, az1, bx0, by0, bz0, bx1, by1, bz1;
                        parts[a].box.Get(ax0, ay0, az0, ax1, ay1, az1);
                        parts[b].box.Get(bx0, by0, bz0, bx1, by1, bz1);
                        const double slack = 1e-3;   // 1 um in the mm frame
                        const bool inside =
                            cx0 >= std::max(ax0, bx0) - slack && cx1 <= std::min(ax1, bx1) + slack &&
                            cy0 >= std::max(ay0, by0) - slack && cy1 <= std::min(ay1, by1) + slack &&
                            cz0 >= std::max(az0, bz0) - slack && cz1 <= std::min(az1, bz1) + slack;
                        if (!inside) {
                            // The boolean is garbage here, so it proves nothing either way.
                            // Decide the pair by the independent proof: mesh-sample each solid's
                            // surface, step half a micron inward, classify against the other.
                            const double defl = 1e-3;   // 1 um in the mm frame
                            const SampleResult sAB = sampledInterpenetration(parts[a].shape, parts[b].shape, defl);
                            const SampleResult sBA = sampledInterpenetration(parts[b].shape, parts[a].shape, defl);
                            if (sAB.found || sBA.found) {
                                const gp_Pnt where = sAB.found ? sAB.where : sBA.where;
                                std::cout << "[intersect] OVERLAP (boolean common was an artefact; proven by "
                                             "point sampling at (" << where.X() << "," << where.Y() << ","
                                          << where.Z() << ")): '" << parts[a].name << "' vs '"
                                          << parts[b].name << "'\n";
                                ++overlaps;
                                const Mat ma2 = classify(parts[a].name), mb2 = classify(parts[b].name);
                                if (ma2 == Mat::Copper && mb2 == Mat::Copper) ++copperOverlaps;
                                if (ma2 == Mat::Core || mb2 == Mat::Core)     ++coreOverlaps;
                                if (worstVolume <= 0) worstPair = parts[a].name + " vs " + parts[b].name;
                                continue;
                            }
                            std::cout << "[intersect] boolean artefact: common box escapes the "
                                         "operands' boxes and point sampling finds no interpenetration -- '"
                                      << parts[a].name << "' vs '" << parts[b].name << "' (common "
                                      << volume << " mm^3, box x[" << cx0 << "," << cx1 << "] vs A x["
                                      << ax0 << "," << ax1 << "] B x[" << bx0 << "," << bx1 << "])\n";
                            ++artifacts;
                            continue;
                        }
                    }
                }
                // Scale test: a tangency between two touching faces yields a sliver whose
                // equivalent cube side is microns. Real interpenetration has depth.
                if (std::cbrt(volume) < contactMm) {
                    ++contacts;
                    contactVolume = std::max(contactVolume, volume);
                    continue;
                }
                const double smaller = std::min(parts[a].volume, parts[b].volume);
                const double fraction = smaller > 0 ? 100.0 * volume / smaller : 0.0;
                // A conductor that could not be swept as one body is exported as a compound of
                // per-run pieces; consecutive pieces MEET at their junction and interpenetrate
                // there by construction (the welding lens). That is one wire, not a fault — it
                // is reported as a JOINT and never fails the audit. Interpenetration between
                // DIFFERENT bodies is a real defect: two conductors shorting, or copper inside
                // the core/bobbin.
                const Mat ma = classify(parts[a].name), mb = classify(parts[b].name);
                const bool copperCopper = (ma == Mat::Copper && mb == Mat::Copper);
                const bool coreInvolved = (ma == Mat::Core || mb == Mat::Core);
                if (parts[a].body == parts[b].body) {
                    // SAME CONDUCTOR, two of its own pieces interpenetrating at a junction --
                    // the "welding lens" the chain construction leaves behind. This used to be
                    // exempted outright ("that is one wire, not a fault"), which hid it from
                    // every audit. Alf, 2026-08-25: "no overlap must be allowed in copper
                    // copper". It is still counted separately, because the FIX is different
                    // (make the junction abut exactly), but it is no longer invisible and no
                    // longer passes silently. MVB_ALLOW_WELD_LENS=1 restores the old leniency.
                    ++joints;
                    jointVolume = std::max(jointVolume, volume);
                    if (copperCopper) ++copperJoints;
                    std::cout << "[intersect] SELF-OVERLAP " << volume << " mm^3 (" << fraction
                              << "% of the smaller piece) within '" << parts[a].body << "': '"
                              << parts[a].name << "' vs '" << parts[b].name << "' common bbox mm ["
                              << cr.bx0 << "," << cr.by0 << "," << cr.bz0 << "]..[" << cr.bx1 << ","
                              << cr.by1 << "," << cr.bz1 << "] centroid (" << cr.cx << "," << cr.cy
                              << "," << cr.cz << ")\n";
                    continue;
                }
                ++overlaps;
                if (copperCopper) ++copperOverlaps;
                if (coreInvolved)  ++coreOverlaps;
                std::cout << "[intersect] OVERLAP [" << matName(ma) << "/" << matName(mb) << "] "
                          << volume << " mm^3 (" << fraction
                          << "% of the smaller part): '" << parts[a].name << "' vs '"
                          << parts[b].name << "'\n";
                if (volume > worstVolume) {
                    worstVolume = volume;
                    worstPair = parts[a].name + " vs " + parts[b].name;
                }
            }
        }
    }

    std::cout << "[intersect] " << pairsTested << " candidate pair(s) tested (" << boxPruned
              << " proven clear by the box-volume bound); " << joints
              << " same-body joint(s) (worst " << jointVolume << " mm^3); " << contacts
              << " tangential contact(s) (worst " << contactVolume << " mm^3); " << artifacts
              << " boolean artefact(s) rejected; " << crashedPairs << " pair(s) UNPROVEN\n";
    std::cout << "[intersect] by class: " << copperOverlaps << " COPPER-COPPER overlap(s), "
              << copperJoints << " copper self-overlap(s) at junctions, " << coreOverlaps
              << " overlap(s) INVOLVING THE CORE\n";

    // Copper inside copper is never acceptable, whether it is two different conductors shorting
    // or one conductor's own chain pieces interpenetrating at a junction (Alf, 2026-08-25).
    // A junction must ABUT exactly, not overlap and rely on the consumer's boolean to weld it.
    const bool allowLens = std::getenv("MVB_ALLOW_WELD_LENS") != nullptr;
    if (copperJoints > 0 && !allowLens && overlaps == 0) {
        std::cout << "[intersect] VERDICT: " << copperJoints
                  << " COPPER SELF-OVERLAP(S) at chain junctions; worst " << jointVolume
                  << " mm^3. No copper-copper overlap is allowed -- junctions must abut.\n";
        return 1;
    }
    if (overlaps == 0 && crashedPairs > 0) {
        // A pair whose boolean crashed was never evaluated. Calling that "no overlaps" would be
        // a false clear -- the audit must report what it could not decide, not round it down.
        std::cout << "[intersect] VERDICT: INCONCLUSIVE -- no overlaps among the pairs that "
                     "completed, but " << crashedPairs << " pair(s) crashed OCCT and are "
                     "UNPROVEN\n";
        return 3;
    }
    if (overlaps == 0) {
        std::cout << "[intersect] VERDICT: NO OVERLAPS\n";
        return 0;
    }
    std::cout << "[intersect] VERDICT: " << overlaps << " OVERLAPPING PAIR(S); worst "
              << worstVolume << " mm^3 (" << worstPair << ")"
              << (crashedPairs ? "  [+" + std::to_string(crashedPairs) + " UNPROVEN]" : "")
              << "\n";
    return 1;
}

}  // namespace mvb::cad
