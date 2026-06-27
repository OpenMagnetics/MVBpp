#include "mvb/MeshableModel.h"

#include "mvb/MagneticBuilder.h"

#include <algorithm>
#include <regex>
#include <utility>

namespace mvb {

namespace {
// Turn label: "<winding> parallel <par> turn <idx>" (MAS Turn::get_name). Core pieces
// are "<base>" or "<base>_<i>". Bobbin is "Bobbin" (planar coils also "FR4Board").
const std::regex kTurnRe(R"(^(.+) parallel (\d+) turn (\d+)$)");
const std::regex kSuffixRe(R"(^(.+)_(\d+)$)");
} // namespace

RegionTag classify_region(const std::string& raw) {
    std::string name = raw;
    // STEP labels can be "Shapes/<name>" or "<name>/<dup>/COMPOUND" — reduce to the
    // canonical short name.
    if (name.rfind("Shapes/", 0) == 0) name = name.substr(7);
    if (auto slash = name.find('/'); slash != std::string::npos) name = name.substr(0, slash);

    RegionTag tag;
    tag.name = name;

    std::smatch m;
    auto match_turn = [&](const std::string& s) -> bool {
        if (std::regex_match(s, m, kTurnRe)) {
            tag.role     = Role::Turn;
            tag.winding  = m[1].str();
            tag.parallel = std::stoi(m[2].str());
            tag.turn     = std::stoi(m[3].str());
            return true;
        }
        return false;
    };

    if (match_turn(name)) return tag;
    // A symmetry/section cut can split a solid and append "_<int>"; strip one and retry
    // the turn pattern.
    std::smatch ms;
    if (std::regex_match(name, ms, kSuffixRe) && match_turn(ms[1].str())) return tag;

    if (name == "Bobbin" || name == "FR4Board") { tag.role = Role::Bobbin; return tag; }

    // Anything else from buildAllNamed (includeBobbin=false) is a core piece. The piece
    // index is a trailing "_<int>" if present, else 0.
    tag.role = Role::Core;
    if (std::regex_match(name, ms, kSuffixRe)) {
        try { tag.corePiece = std::stoi(ms[2].str()); return tag; }
        catch (...) {}
    }
    tag.corePiece = 0;
    return tag;
}

MeshableModel buildMeshable(const OpenMagnetics::Magnetic& magnetic,
                            const MeshableOptions& opt) {
    MagneticBuilder builder;
    std::vector<NamedShape> named = builder.buildAllNamed(
        magnetic, opt.includeBobbin, /*symmetryPlanes=*/0,
        opt.wirePolygonSegments, opt.corePolygonSegments, opt.paintCoating);

    MeshableModel model;

    // Symmetry reduction with reporting (mirrors apply_symmetry(), but records which
    // planes/halves it used so the consumer can apply the matching BCs + multiplicity).
    if (opt.symmetryPlanes > 0 && !named.empty()) {
        const int n = std::min(opt.symmetryPlanes, 3);
        const SymmetryResult res = analyze_symmetry(named);
        std::vector<std::pair<SymmetryPlane, SymmetryHalf>> cuts;
        for (int i = 0; i < n && i < static_cast<int>(res.valid_planes.size()); ++i) {
            cuts.emplace_back(res.valid_planes[i], SymmetryHalf::Positive);
            model.symmetry.push_back(
                {res.valid_planes[i], SymmetryHalf::Positive, WallBC::FluxTangential});
        }
        if (!cuts.empty()) {
            const ShapeBBox bb = aggregate_bbox(named);
            named = cut_to_region(named, cuts, bb);
            model.multiplicity = 1 << static_cast<int>(cuts.size());
        }
    }

    model.regions.reserve(named.size());
    for (auto& ns : named) {
        RegionTag tag = classify_region(ns.name);
        model.regions.push_back({std::move(ns), std::move(tag)});
    }
    return model;
}

} // namespace mvb
