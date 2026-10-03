// mvbpp_solidstats -- per-solid geometry statistics, indexed by the SAME volume tags that
// Moved verbatim from OMFEM tools/omfem_solidstats.cpp at d92be53 (ABT #1588, step 8); renamed. No logic changed.
// omfem_stepcheck --isolate uses (both go through gmsh::model::occ::importShapes, so tag N
// here is tag N there).
//
// Written 2026-08-31 to answer a specific question: 21_interleaved_flyback's LATE solids cost
// 6-13x more to fragment than its early ones (early 1..26 = 216 s, late 51..76 = 586 s at equal
// count; marginal solids 45-46 ~63 s each vs 4-9 s early; baseline 02 = 4.9 s/solid). Face
// count and surface type are already excluded at the WHOLE-FILE level (21 has FEWER faces than
// the green 02: 4237 vs 10630, both zero periodic/seam/micro-edge/sliver) -- but that aggregate
// hides per-solid structure, which is what this dumps.
//
// OCC's fragment cost is driven by the number of face-face intersection CANDIDATES, i.e. by how
// many other solids each solid's faces must be tested against, so the neighbour count (how many
// other solids' bounding boxes a solid overlaps) is reported alongside the face counts.
//
// Exit 0 on success, 2 on usage/read error.
#include <gmsh.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace {

struct SolidStat {
    int tag = 0;
    int faces = 0;
    int edges = 0;
    double vol = 0.0;          // mm^3, from the mass-properties gmsh exposes
    double bbox[6] = {0, 0, 0, 0, 0, 0};
    double diag = 0.0;         // mm
    int neighbours = 0;        // solids whose bbox overlaps this one's
    std::map<std::string, int> surf;   // surface-type histogram
};

bool bboxOverlap(const double a[6], const double b[6], double tol) {
    for (int i = 0; i < 3; ++i) {
        if (a[i] > b[i + 3] + tol) return false;
        if (b[i] > a[i + 3] + tol) return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <in.step> [--range A-B] [--csv]\n", argv[0]);
        return 2;
    }
    const std::string path = argv[1];
    int lo = 0, hi = 0;
    bool csv = false;
    for (int i = 2; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--csv")) csv = true;
        else if (!std::strcmp(argv[i], "--range") && i + 1 < argc) {
            std::sscanf(argv[i + 1], "%d-%d", &lo, &hi);
        }
    }

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);

    gmsh::vectorpair imported;
    try {
        gmsh::model::occ::importShapes(path, imported);
        gmsh::model::occ::synchronize();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "read failed: %s\n", e.what());
        gmsh::finalize();
        return 2;
    }

    gmsh::vectorpair vols;
    gmsh::model::getEntities(vols, 3);

    std::vector<SolidStat> stats;
    stats.reserve(vols.size());
    for (auto& v : vols) {
        SolidStat s;
        s.tag = v.second;
        try {
            gmsh::model::getBoundingBox(3, s.tag, s.bbox[0], s.bbox[1], s.bbox[2],
                                        s.bbox[3], s.bbox[4], s.bbox[5]);
            const double dx = s.bbox[3] - s.bbox[0];
            const double dy = s.bbox[4] - s.bbox[1];
            const double dz = s.bbox[5] - s.bbox[2];
            s.diag = std::sqrt(dx * dx + dy * dy + dz * dz);
        } catch (...) {}

        gmsh::vectorpair bnd;
        try {
            gmsh::model::getBoundary({{3, s.tag}}, bnd, false, false, false);
        } catch (...) {}
        s.faces = static_cast<int>(bnd.size());
        for (auto& f : bnd) {
            std::string ty = "?";
            try { gmsh::model::getType(2, std::abs(f.second), ty); } catch (...) {}
            s.surf[ty]++;
            gmsh::vectorpair ed;
            try { gmsh::model::getBoundary({{2, std::abs(f.second)}}, ed, false, false, false); } catch (...) {}
            s.edges += static_cast<int>(ed.size());
        }
        try {
            double cx, cy, cz;
            gmsh::model::occ::getCenterOfMass(3, s.tag, cx, cy, cz);
            gmsh::model::occ::getMass(3, s.tag, s.vol);
        } catch (...) {}
        stats.push_back(s);
    }

    // Neighbour count: how many other solids each one's bbox touches. This is the population
    // OCC must test faces against, so it -- not the raw face count -- is what scales the boolean.
    // Counted WITHIN --range, because that is the set the fragment actually sees: stepcheck
    // --isolate removes everything else before fragmenting, so neighbours outside the range
    // cost nothing. (Counting against the whole file instead hid the early/late difference
    // entirely -- 11.5 vs 11.0 -- and made the two blocks look identical.)
    auto inRange = [&](int tag) { return !lo || (tag >= lo && tag <= hi); };
    for (size_t i = 0; i < stats.size(); ++i) {
        if (!inRange(stats[i].tag)) continue;
        for (size_t j = 0; j < stats.size(); ++j) {
            if (i == j || !inRange(stats[j].tag)) continue;
            if (bboxOverlap(stats[i].bbox, stats[j].bbox, 1e-9)) stats[i].neighbours++;
        }
    }

    // --pairs: for every bbox-overlapping pair WITHIN the range, the true minimum distance
    // (OCC BRepExtrema via gmsh). Bbox overlap says two solids MIGHT interact; this says HOW
    // they meet, which is what actually costs the boolean. Tangential contact (two curved
    // faces meeting at ~0 separation without crossing) is OCC's worst case: the intersector
    // marches numerically instead of solving, and separation grows only quadratically away
    // from the touch point. Needed because pair COUNT does not explain the cost -- early
    // 6..31 has MORE within-subset pairs (~69) than late 51..76 (~57) yet costs a third as much.
    if (std::getenv("OMFEM_PAIRS") || [&] {
            for (int i = 2; i < argc; ++i) if (!std::strcmp(argv[i], "--pairs")) return true;
            return false; }()) {
        int touch = 0, near_ = 0, apart = 0, pairs = 0;
        double sumTouch = 0.0;
        for (size_t i = 0; i < stats.size(); ++i) {
            if (!inRange(stats[i].tag)) continue;
            for (size_t j = i + 1; j < stats.size(); ++j) {
                if (!inRange(stats[j].tag)) continue;
                if (!bboxOverlap(stats[i].bbox, stats[j].bbox, 1e-9)) continue;
                double d = -1, ax, ay, az, bx, by, bz;
                try {
                    gmsh::model::occ::getDistance(3, stats[i].tag, 3, stats[j].tag,
                                                  d, ax, ay, az, bx, by, bz);
                } catch (...) { continue; }
                ++pairs;
                if (d <= 1e-9) { ++touch; sumTouch += d; }
                else if (d <= 1e-3) ++near_;
                else ++apart;
            }
        }
        std::printf("PAIRS range=%d-%d bboxPairs=%d  touching(<=1e-9)=%d  near(<=1e-3mm)=%d  "
                    "apart=%d\n", lo, hi, pairs, touch, near_, apart);
    }

    if (csv) std::printf("tag,faces,edges,vol_mm3,diag_mm,neighbours,types\n");
    long long fSum = 0, nSum = 0;
    int counted = 0;
    for (auto& s : stats) {
        if (lo && (s.tag < lo || s.tag > hi)) continue;
        std::string types;
        for (auto& kv : s.surf) {
            if (!types.empty()) types += " ";
            types += kv.first + ":" + std::to_string(kv.second);
        }
        if (csv) {
            std::printf("%d,%d,%d,%.6f,%.4f,%d,%s\n", s.tag, s.faces, s.edges,
                        s.vol * 1e9, s.diag, s.neighbours, types.c_str());
        } else {
            std::printf("solid %3d  faces=%4d edges=%5d vol=%12.6f mm3 diag=%8.4f mm "
                        "nbrs=%3d  [%s]\n",
                        s.tag, s.faces, s.edges, s.vol * 1e9, s.diag, s.neighbours,
                        types.c_str());
        }
        fSum += s.faces;
        nSum += s.neighbours;
        ++counted;
    }
    if (counted && !csv) {
        std::printf("SUMMARY n=%d meanFaces=%.1f meanNbrs=%.1f\n",
                    counted, double(fSum) / counted, double(nSum) / counted);
    }

    gmsh::finalize();
    return 0;
}
