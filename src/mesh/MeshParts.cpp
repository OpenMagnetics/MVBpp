// The meshed magnetic without its air (see partsOnly in Mesher.h).
#include "mvb/mesh/Mesher.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace mvb::mesh {

namespace {

int nodesOf(ElementType t) {
    switch (t) {
        case ElementType::Point: return 1;
        case ElementType::Line2: return 2;
        case ElementType::Tri3: return 3;
        case ElementType::Quad4: return 4;
        case ElementType::Tet4: return 4;
        case ElementType::Pyramid5: return 5;
        case ElementType::Prism6: return 6;
        case ElementType::Hex8: return 8;
    }
    throw std::runtime_error("partsOnly: unknown element type " + std::to_string(int(t)));
}

bool isVolume(ElementType t) {
    return t == ElementType::Tet4 || t == ElementType::Pyramid5 || t == ElementType::Prism6 || t == ElementType::Hex8;
}

// The faces of a volume element, as local node indices in gmsh's node order.
const std::vector<std::vector<int>>& facesOf(ElementType t) {
    static const std::vector<std::vector<int>> tet{{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    static const std::vector<std::vector<int>> pyramid{{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};
    static const std::vector<std::vector<int>> prism{{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
    static const std::vector<std::vector<int>> hex{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                                   {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    switch (t) {
        case ElementType::Tet4: return tet;
        case ElementType::Pyramid5: return pyramid;
        case ElementType::Prism6: return prism;
        case ElementType::Hex8: return hex;
        default: throw std::runtime_error("partsOnly: element type " + std::to_string(int(t)) + " has no faces");
    }
}

using FaceKey = std::vector<std::int64_t>;   // the face's node ids, sorted

FaceKey keyOf(std::vector<std::int64_t> nodes) {
    std::sort(nodes.begin(), nodes.end());
    return nodes;
}

}  // namespace

Mesh partsOnly(const Mesh& mesh) {
    std::set<int> air, parts;
    for (const auto& r : mesh.regions) {
        if (r.dimension != 3) continue;
        (r.name == "air" ? air : parts).insert(r.tag);
    }
    if (air.empty())
        throw std::runtime_error("partsOnly: the mesh has no volume region named 'air' to leave out");
    if (parts.empty())
        throw std::runtime_error("partsOnly: the mesh has no volume region besides the air");

    std::unordered_map<std::int64_t, std::size_t> at;
    for (std::size_t i = 0; i < mesh.nodeIds.size(); ++i) at[mesh.nodeIds[i]] = i;
    auto point = [&](std::int64_t id) {
        auto it = at.find(id);
        if (it == at.end()) throw std::runtime_error("partsOnly: an element names node " + std::to_string(id) + ", which the mesh lacks");
        return std::array<double, 3>{mesh.xyz[3 * it->second], mesh.xyz[3 * it->second + 1], mesh.xyz[3 * it->second + 2]};
    };

    // Every face of every part element: which faces the part elements have (a face set is kept
    // only on them), and per part region how often each face occurs and, the first time, its nodes
    // ordered so the normal points out of the element. A face a region's elements share twice is
    // inside it; once, it is on the region's skin.
    struct SkinFace {
        std::vector<std::int64_t> nodes;
        int count = 0;
    };
    std::set<FaceKey> partFaces;
    std::set<std::int64_t> partNodes;
    std::map<int, std::map<FaceKey, std::size_t>> faceIndex;   // region -> face -> index into skins[region]
    std::map<int, std::vector<SkinFace>> skins;
    for (const auto& b : mesh.blocks) {
        if (!parts.count(b.region)) continue;
        if (!isVolume(b.type))
            throw std::runtime_error("partsOnly: a volume region holds a non-volume element block");
        const int n = nodesOf(b.type);
        for (std::size_t e = 0; e < b.ids.size(); ++e) {
            std::vector<std::int64_t> el(b.nodes.begin() + std::ptrdiff_t(n * e), b.nodes.begin() + std::ptrdiff_t(n * (e + 1)));
            std::array<double, 3> centre{0, 0, 0};
            for (auto id : el) {
                partNodes.insert(id);
                const auto p = point(id);
                for (int c = 0; c < 3; ++c) centre[c] += p[c] / n;
            }
            for (const auto& local : facesOf(b.type)) {
                std::vector<std::int64_t> face;
                for (int k : local) face.push_back(el[std::size_t(k)]);
                // Newell's normal against the element centre: reversed if it points inwards.
                std::array<double, 3> normal{0, 0, 0}, mid{0, 0, 0};
                for (std::size_t k = 0; k < face.size(); ++k) {
                    const auto p = point(face[k]), q = point(face[(k + 1) % face.size()]);
                    normal[0] += (p[1] - q[1]) * (p[2] + q[2]);
                    normal[1] += (p[2] - q[2]) * (p[0] + q[0]);
                    normal[2] += (p[0] - q[0]) * (p[1] + q[1]);
                    for (int c = 0; c < 3; ++c) mid[c] += p[c] / double(face.size());
                }
                double outward = 0;
                for (int c = 0; c < 3; ++c) outward += normal[c] * (mid[c] - centre[c]);
                if (outward < 0) std::reverse(face.begin() + 1, face.end());
                auto key = keyOf(face);
                partFaces.insert(key);
                auto& index = faceIndex[b.region];
                auto [it, inserted] = index.emplace(std::move(key), skins[b.region].size());
                if (inserted) skins[b.region].push_back({face, 0});
                ++skins[b.region][it->second].count;
            }
        }
    }

    Mesh out;
    out.lengthUnit = mesh.lengthUnit;
    std::set<int> keptRegions(parts.begin(), parts.end());
    std::int64_t maxElement = 0;
    for (const auto& b : mesh.blocks)
        for (auto id : b.ids) maxElement = std::max(maxElement, id);
    for (const auto& b : mesh.blocks) {
        if (air.count(b.region)) continue;
        if (parts.count(b.region)) {
            out.blocks.push_back(b);
            continue;
        }
        // A face, edge or point set: kept where it lies on the parts (a terminal is a conductor's
        // face; the air box's 'outer' is not), element by element.
        const int n = nodesOf(b.type);
        ElementBlock kept{b.type, b.region, b.entity, {}, {}};
        for (std::size_t e = 0; e < b.ids.size(); ++e) {
            std::vector<std::int64_t> el(b.nodes.begin() + std::ptrdiff_t(n * e), b.nodes.begin() + std::ptrdiff_t(n * (e + 1)));
            bool on;
            if (b.type == ElementType::Tri3 || b.type == ElementType::Quad4) on = partFaces.count(keyOf(el)) > 0;
            else if (isVolume(b.type)) throw std::runtime_error("partsOnly: a face region holds a volume element block");
            else on = std::all_of(el.begin(), el.end(), [&](std::int64_t id) { return partNodes.count(id) > 0; });
            if (!on) continue;
            kept.ids.push_back(b.ids[e]);
            kept.nodes.insert(kept.nodes.end(), el.begin(), el.end());
        }
        if (kept.ids.empty()) continue;
        keptRegions.insert(b.region);
        out.blocks.push_back(std::move(kept));
    }

    int maxTag = 0;
    std::set<std::string> names;
    for (const auto& r : mesh.regions) {
        maxTag = std::max(maxTag, r.tag);
        names.insert(r.name);
        if (keptRegions.count(r.tag)) out.regions.push_back(r);
    }
    // Each part's skin, a face set of its own named after the part. Elementary tag 0: the skin is
    // not one of the geometry's surfaces but the union of the part's.
    for (const auto& r : mesh.regions) {
        if (!parts.count(r.tag)) continue;
        const std::string name = r.name + " skin";
        if (names.count(name))
            throw std::runtime_error("partsOnly: the mesh already has a region named '" + name + "'");
        Region skin{++maxTag, 2, name, {MaterialKind::None, {}}};
        ElementBlock tris{ElementType::Tri3, skin.tag, 0, {}, {}}, quads{ElementType::Quad4, skin.tag, 0, {}, {}};
        for (const auto& f : skins[r.tag]) {
            if (f.count > 2)
                throw std::runtime_error("partsOnly: a face of region '" + r.name + "' is shared by " +
                                         std::to_string(f.count) + " of its elements (a non-manifold mesh)");
            if (f.count != 1) continue;
            auto& block = f.nodes.size() == 3 ? tris : quads;
            block.ids.push_back(++maxElement);
            block.nodes.insert(block.nodes.end(), f.nodes.begin(), f.nodes.end());
        }
        if (tris.ids.empty() && quads.ids.empty())
            throw std::runtime_error("partsOnly: region '" + r.name + "' has no skin (no elements?)");
        out.regions.push_back(skin);
        if (!tris.ids.empty()) out.blocks.push_back(std::move(tris));
        if (!quads.ids.empty()) out.blocks.push_back(std::move(quads));
    }

    // The nodes the kept elements use, in the mesh's own order and numbering.
    std::set<std::int64_t> used;
    for (const auto& b : out.blocks) used.insert(b.nodes.begin(), b.nodes.end());
    for (std::size_t i = 0; i < mesh.nodeIds.size(); ++i) {
        if (!used.count(mesh.nodeIds[i])) continue;
        out.nodeIds.push_back(mesh.nodeIds[i]);
        out.xyz.insert(out.xyz.end(), mesh.xyz.begin() + std::ptrdiff_t(3 * i), mesh.xyz.begin() + std::ptrdiff_t(3 * i + 3));
    }
    return out;
}

}  // namespace mvb::mesh
