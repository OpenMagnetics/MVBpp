// Mesh exporters, the importer the round-trip gate reads them back with, and the identity diff
// (ABT #1588). Every writer reads only the neutral mvb::mesh::Mesh, so every format sees the same
// mesh, and every refusal throws: a file that cannot say what the mesh is must not be written.
#include "mvb/mesh/Mesher.h"

#include "Definitions.h"
#include "MAS.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace mvb::mesh {

namespace {

struct TypeInfo {
    ElementType type;
    int nodes;
    int dimension;
    int gmsh;                // msh2 element type code
    const char* nastran;     // bulk-data card, nullptr where Nastran has none we write
    const char* abaqus;      // INP element type, as gmsh's own INP writer names it
    int vtk;                 // VTK cell type
    std::array<int, 8> vtkOrder;  // VTK node k = gmsh node vtkOrder[k] (gmsh's getVertexVTK)
};

constexpr TypeInfo kTypes[] = {
    {ElementType::Point,    1, 0, 15, nullptr,  nullptr, 1,  {0}},
    {ElementType::Line2,    2, 1, 1,  nullptr,  nullptr, 3,  {0, 1}},
    {ElementType::Tri3,     3, 2, 2,  "CTRIA3", "CPS3",  5,  {0, 1, 2}},
    {ElementType::Quad4,    4, 2, 3,  "CQUAD4", "CPS4",  9,  {0, 1, 2, 3}},
    {ElementType::Tet4,     4, 3, 4,  "CTETRA", "C3D4",  10, {0, 1, 2, 3}},
    {ElementType::Pyramid5, 5, 3, 7,  "CPYRAM", "C3D5",  14, {0, 1, 2, 3, 4}},
    {ElementType::Prism6,   6, 3, 6,  "CPENTA", "C3D6",  13, {0, 2, 1, 3, 5, 4}},
    {ElementType::Hex8,     8, 3, 5,  "CHEXA",  "C3D8",  12, {0, 1, 2, 3, 4, 5, 6, 7}},
};

const TypeInfo& info(ElementType t) {
    for (const auto& i : kTypes)
        if (i.type == t) return i;
    throw std::runtime_error("mvb::mesh: unknown element type " + std::to_string(int(t)));
}

const TypeInfo& infoFromGmsh(int code) {
    for (const auto& i : kTypes)
        if (i.gmsh == code) return i;
    throw std::runtime_error("importMesh(msh2): element type " + std::to_string(code) +
                             " is not one this mesh model carries");
}

const TypeInfo& infoFromNastran(const std::string& card) {
    for (const auto& i : kTypes)
        if (i.nastran && card == i.nastran) return i;
    throw std::runtime_error("importMesh(bdf): unknown element card '" + card + "'");
}

// Structural checks every writer runs first, so no format writes a mesh the model cannot hold.
std::map<int, const Region*> checkMesh(const Mesh& mesh, const char* who) {
    const std::string w = who;
    if (mesh.lengthUnit != "m")
        throw std::runtime_error(w + ": mesh coordinates must be in metres, the mesh says '" +
                                 mesh.lengthUnit + "'");
    if (mesh.xyz.size() != 3 * mesh.nodeIds.size())
        throw std::runtime_error(w + ": " + std::to_string(mesh.nodeIds.size()) + " node ids but " +
                                 std::to_string(mesh.xyz.size()) + " coordinates");
    std::set<std::int64_t> nodeSet(mesh.nodeIds.begin(), mesh.nodeIds.end());
    if (nodeSet.size() != mesh.nodeIds.size())
        throw std::runtime_error(w + ": duplicate node id");
    std::map<int, const Region*> regions;
    for (const auto& r : mesh.regions) {
        if (!regions.emplace(r.tag, &r).second)
            throw std::runtime_error(w + ": region tag " + std::to_string(r.tag) +
                                     " used twice (tags must be unique across dimensions)");
        if (r.name.empty())
            throw std::runtime_error(w + ": region " + std::to_string(r.tag) + " has no name");
    }
    std::set<std::int64_t> elementIds;
    for (const auto& b : mesh.blocks) {
        const auto& ti = info(b.type);
        auto it = regions.find(b.region);
        if (it == regions.end())
            throw std::runtime_error(w + ": a block refers to region " + std::to_string(b.region) +
                                     ", which the mesh does not declare");
        if (it->second->dimension != ti.dimension)
            throw std::runtime_error(w + ": region '" + it->second->name + "' is " +
                                     std::to_string(it->second->dimension) +
                                     "D but holds " + std::to_string(ti.dimension) + "D elements");
        if (b.nodes.size() != b.ids.size() * std::size_t(ti.nodes))
            throw std::runtime_error(w + ": block of region '" + it->second->name + "' has " +
                                     std::to_string(b.ids.size()) + " elements but " +
                                     std::to_string(b.nodes.size()) + " node references");
        for (auto id : b.ids)
            if (!elementIds.insert(id).second)
                throw std::runtime_error(w + ": element id " + std::to_string(id) + " used twice");
        for (auto n : b.nodes)
            if (!nodeSet.count(n))
                throw std::runtime_error(w + ": an element of region '" + it->second->name +
                                         "' uses node " + std::to_string(n) + ", which does not exist");
    }
    return regions;
}

double unitScale(const std::string& unit, const char* who) {
    if (unit == "m") return 1.0;
    if (unit == "mm") return 1000.0;
    throw std::runtime_error(std::string(who) + ": unit must be \"m\" or \"mm\", got \"" + unit + "\"");
}

std::string g17(double v) {
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

// ---- msh2 ----------------------------------------------------------------------------------

std::string writeMsh2(const Mesh& mesh, const std::string& unit) {
    // msh2 has no unit field: a millimetre file would be read as metres by every reader.
    if (unit != "m")
        throw std::runtime_error("exportMesh(msh2): the format has no unit field and is written in "
                                 "metres only; asked for \"" + unit + "\"");
    checkMesh(mesh, "exportMesh(msh2)");
    std::string out = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
    out += "$PhysicalNames\n" + std::to_string(mesh.regions.size()) + "\n";
    for (const auto& r : mesh.regions) {
        if (r.name.find('"') != std::string::npos)
            throw std::runtime_error("exportMesh(msh2): region name '" + r.name + "' contains a quote");
        out += std::to_string(r.dimension) + " " + std::to_string(r.tag) + " \"" + r.name + "\"\n";
    }
    out += "$EndPhysicalNames\n$Nodes\n" + std::to_string(mesh.nodeIds.size()) + "\n";
    for (std::size_t i = 0; i < mesh.nodeIds.size(); ++i)
        out += std::to_string(mesh.nodeIds[i]) + " " + g17(mesh.xyz[3 * i]) + " " +
               g17(mesh.xyz[3 * i + 1]) + " " + g17(mesh.xyz[3 * i + 2]) + "\n";
    std::size_t count = 0;
    for (const auto& b : mesh.blocks) count += b.ids.size();
    out += "$EndNodes\n$Elements\n" + std::to_string(count) + "\n";
    for (const auto& b : mesh.blocks) {
        const auto& ti = info(b.type);
        for (std::size_t e = 0; e < b.ids.size(); ++e) {
            out += std::to_string(b.ids[e]) + " " + std::to_string(ti.gmsh) + " 2 " +
                   std::to_string(b.region) + " " + std::to_string(b.entity);
            for (int k = 0; k < ti.nodes; ++k) out += " " + std::to_string(b.nodes[e * ti.nodes + k]);
            out += "\n";
        }
    }
    out += "$EndElements\n";
    return out;
}

void expectLine(std::istream& in, const std::string& want, const char* who) {
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line != want)
            throw std::runtime_error(std::string(who) + ": expected '" + want + "', found '" + line + "'");
        return;
    }
    throw std::runtime_error(std::string(who) + ": file ends before '" + want + "'");
}

Mesh readMsh2(const std::string& text) {
    const char* who = "importMesh(msh2)";
    std::istringstream in(text);
    Mesh mesh;
    expectLine(in, "$MeshFormat", who);
    std::string version;
    int fileType = -1, dataSize = 0;
    in >> version >> fileType >> dataSize;
    if (version.rfind("2.", 0) != 0 || fileType != 0)
        throw std::runtime_error(std::string(who) + ": not an ASCII msh 2 file (version '" + version +
                                 "', file type " + std::to_string(fileType) + ")");
    expectLine(in, "$EndMeshFormat", who);
    std::string section;
    bool sawNodes = false, sawElements = false;
    while (in >> section) {
        if (section == "$PhysicalNames") {
            std::size_t n;
            in >> n;
            for (std::size_t i = 0; i < n; ++i) {
                Region r;
                in >> r.dimension >> r.tag;
                std::string rest;
                std::getline(in, rest);
                auto a = rest.find('"'), z = rest.rfind('"');
                if (a == std::string::npos || z == a)
                    throw std::runtime_error(std::string(who) + ": physical name without quotes: " + rest);
                r.name = rest.substr(a + 1, z - a - 1);
                mesh.regions.push_back(r);
            }
            expectLine(in, "$EndPhysicalNames", who);
        } else if (section == "$Nodes") {
            std::size_t n;
            in >> n;
            mesh.nodeIds.resize(n);
            mesh.xyz.resize(3 * n);
            for (std::size_t i = 0; i < n; ++i)
                in >> mesh.nodeIds[i] >> mesh.xyz[3 * i] >> mesh.xyz[3 * i + 1] >> mesh.xyz[3 * i + 2];
            if (!in) throw std::runtime_error(std::string(who) + ": truncated $Nodes");
            expectLine(in, "$EndNodes", who);
            sawNodes = true;
        } else if (section == "$Elements") {
            std::size_t n;
            in >> n;
            for (std::size_t i = 0; i < n; ++i) {
                std::int64_t id;
                int code, ntags;
                in >> id >> code >> ntags;
                if (!in || ntags < 2)
                    throw std::runtime_error(std::string(who) + ": element " + std::to_string(id) +
                                             " needs physical and elementary tags");
                int phys, entity, extra;
                in >> phys >> entity;
                for (int t = 2; t < ntags; ++t) in >> extra;  // partition tags: not in the model
                const auto& ti = infoFromGmsh(code);
                if (mesh.blocks.empty() || mesh.blocks.back().type != ti.type ||
                    mesh.blocks.back().region != phys || mesh.blocks.back().entity != entity)
                    mesh.blocks.push_back(ElementBlock{ti.type, phys, entity, {}, {}});
                auto& b = mesh.blocks.back();
                b.ids.push_back(id);
                for (int k = 0; k < ti.nodes; ++k) {
                    std::int64_t v;
                    in >> v;
                    b.nodes.push_back(v);
                }
            }
            if (!in) throw std::runtime_error(std::string(who) + ": truncated $Elements");
            expectLine(in, "$EndElements", who);
            sawElements = true;
        } else if (section.size() > 1 && section[0] == '$') {
            // A section this model does not carry ($NodeData, $Periodic...): refuse rather than
            // drop it, so a comparison never passes on half a file.
            throw std::runtime_error(std::string(who) + ": unsupported section " + section);
        } else {
            throw std::runtime_error(std::string(who) + ": unexpected token '" + section + "'");
        }
    }
    if (!sawNodes || !sawElements)
        throw std::runtime_error(std::string(who) + ": missing $Nodes or $Elements");
    return mesh;
}

// ---- Nastran bulk data (large field) -----------------------------------------------------------

// A real in a 16-character large field, with the most significant digits that fit, in
// Nastran's exponent form without the E ("1.234567890123-5"). Reads back with readReal.
std::string real16(double v) {
    if (!std::isfinite(v)) throw std::runtime_error("exportMesh(bdf): non-finite value");
    if (v == 0.0) return "0.";
    for (int digits = 16; digits >= 1; --digits) {
        char buf[48];
        std::snprintf(buf, sizeof buf, "%.*E", digits - 1, v);
        std::string s = buf;
        auto e = s.find('E');
        std::string mant = s.substr(0, e);
        if (mant.find('.') == std::string::npos) mant += ".";
        while (mant.back() == '0') mant.pop_back();   // 7.500 -> 7.5: same value, shorter field
        int exp = std::stoi(s.substr(e + 1));
        std::string out = mant + (exp < 0 ? "-" : "+") + std::to_string(std::abs(exp));
        if (out.size() <= 16) return out;
    }
    throw std::runtime_error("exportMesh(bdf): cannot fit a real into 16 characters");
}

double readReal(std::string f, const char* who) {
    while (!f.empty() && f.back() == ' ') f.pop_back();
    while (!f.empty() && f.front() == ' ') f.erase(0, 1);
    if (f.empty()) throw std::runtime_error(std::string(who) + ": blank real field");
    // Nastran exponent without E: a sign after the first character that is not after E/D.
    for (std::size_t i = 1; i < f.size(); ++i)
        if ((f[i] == '-' || f[i] == '+') && f[i - 1] != 'E' && f[i - 1] != 'e' && f[i - 1] != 'D') {
            f.insert(i, "E");
            break;
        }
    std::size_t used = 0;
    double v = std::stod(f, &used);
    if (used != f.size()) throw std::runtime_error(std::string(who) + ": bad real '" + f + "'");
    return v;
}

std::string field16(const std::string& s) {
    if (s.size() > 16) throw std::runtime_error("exportMesh(bdf): field '" + s + "' is over 16 characters");
    return s + std::string(16 - s.size(), ' ');
}

// One large-field card: a name and any number of 16-char data fields, four per line.
std::string card(const std::string& name, const std::vector<std::string>& fields) {
    std::string head = name + "*";
    std::string out = head + std::string(8 - head.size(), ' ');
    for (std::size_t i = 0; i < fields.size(); ++i) {
        if (i && i % 4 == 0) out += "\n*       ";
        out += field16(fields[i]);
    }
    // Trailing blanks are not data; a line ends at its last field.
    while (!out.empty() && out.back() == ' ') out.pop_back();
    return out + "\n";
}

std::string id16(std::int64_t id, const char* what) {
    if (id < 1 || id > 99999999)
        throw std::runtime_error(std::string("exportMesh(bdf): ") + what + " id " + std::to_string(id) +
                                 " is outside Nastran's 1..99999999");
    return std::to_string(id);
}

// MAT4 values from the region's MAS record: K, CP, RHO in SI, or a throw that names every value
// MAS lacks (the BDF writes a complete thermal card or none; a blank field is Nastran's default,
// K = 0, RHO = 1, never "unknown").
struct Mat4 { double k, cp, rho; };

Mat4 mat4From(const Region& r) {
    const auto& m = r.material.record;
    const std::string name = m.value("name", std::string("(unnamed)"));
    std::vector<std::string> missing;
    double k = 0, cp = 0, rho = 0;
    auto dim = [&](const char* key, double& out, const char* what) {
        if (!m.contains(key) || m[key].is_null()) { missing.push_back(what); return; }
        out = OpenMagnetics::resolve_dimensional_values(m[key].get<MAS::Dimension>());
    };
    auto num = [&](const char* key, double& out, const char* what) {
        if (!m.contains(key) || m[key].is_null()) { missing.push_back(what); return; }
        if (!m[key].is_number())
            throw std::runtime_error("exportMesh(bdf): material '" + name + "' " + key + " is not a number");
        out = m[key].get<double>();
    };
    switch (r.material.kind) {
    case MaterialKind::Core:
        dim("heatConductivity", k, "heatConductivity");
        dim("heatCapacity", cp, "heatCapacity");
        num("density", rho, "density");
        break;
    case MaterialKind::Insulation:
        num("thermalConductivity", k, "thermalConductivity");
        num("specificHeat", cp, "specificHeat");
        num("density", rho, "density");
        break;
    case MaterialKind::Wire:
        // MAS gives a conductor's conductivity as a temperature table; MAT4's K is one number and
        // choosing the temperature is not this writer's call.
        if (m.contains("thermalConductivity"))
            missing.push_back("thermalConductivity at one temperature (MAS gives a table; MATT4 export "
                              "needs a reference temperature)");
        else
            missing.push_back("thermalConductivity");
        missing.push_back("specific heat (the MAS wire material schema has no field for it)");
        num("density", rho, "density");
        break;
    default:
        throw std::runtime_error("exportMesh(bdf): mat4From on a region without a material");
    }
    if (!missing.empty()) {
        std::string list;
        for (const auto& s : missing) list += (list.empty() ? "" : ", ") + s;
        throw std::runtime_error("exportMesh(bdf): region '" + r.name + "' (MAS material '" + name +
                                 "') cannot get a complete MAT4: MAS lacks " + list);
    }
    return {k, cp, rho};
}

const char* kindName(MaterialKind k) {
    switch (k) {
    case MaterialKind::None: return "none";
    case MaterialKind::Core: return "core";
    case MaterialKind::Wire: return "wire";
    case MaterialKind::Insulation: return "insulation";
    default: return "unset";
    }
}

nlohmann::json emData(const Region& r) {
    nlohmann::json em = nlohmann::json::object();
    for (const char* key : {"permeability", "resistivity", "permittivity", "relativePermittivity"})
        if (r.material.record.contains(key)) em[key] = r.material.record[key];
    return em;
}

void requireMaterial(const Region& r, const char* who) {
    if (r.dimension == 3 && r.material.kind == MaterialKind::Unset)
        throw std::runtime_error(std::string(who) + ": volume region '" + r.name +
                                 "' has no material set (MaterialKind::None for one with no MAS material)");
    if (r.material.kind != MaterialKind::Unset && r.material.kind != MaterialKind::None &&
        !r.material.record.contains("name"))
        throw std::runtime_error(std::string(who) + ": region '" + r.name + "' has a " +
                                 kindName(r.material.kind) + " material with no MAS record");
}

std::string writeBdf(const Mesh& mesh, const std::string& unit) {
    const auto regions = checkMesh(mesh, "exportMesh(bdf)");
    const double s = unitScale(unit, "exportMesh(bdf)");
    // Material values in the unit system that goes with the length: m-kg-s-K, or mm-t-s-K
    // (K unchanged, CP x 1e6, RHO x 1e-12).
    const double cpScale = unit == "mm" ? 1e6 : 1.0, rhoScale = unit == "mm" ? 1e-12 : 1.0;
    std::string out;
    out += "$ Written by MVB++ mvbpp_mesh (Nastran bulk data, large field).\n";
    out += "$ LENGTH UNIT: " + unit + "\n";
    out += std::string("$ UNIT SYSTEM: ") + (unit == "mm" ? "mm-t-s-K" : "m-kg-s-K") + "\n";
    out += "$ Materials are MAS's. MAT4 carries K, CP, RHO; the electromagnetic data Nastran has no\n";
    out += "$ card for is in the $ EM lines and in the .materials.json companion. A region whose\n";
    out += "$ material is 'none' (air) has a property and no MAT card, deliberately.\n";
    out += "BEGIN BULK\n";
    for (const auto& r : mesh.regions) {
        requireMaterial(r, "exportMesh(bdf)");
        if (r.dimension < 2)
            throw std::runtime_error("exportMesh(bdf): region '" + r.name + "' is " +
                                     std::to_string(r.dimension) + "D; no Nastran card is written for it");
        const std::string pid = id16(r.tag, "property");
        out += "$ REGION " + std::to_string(r.tag) + " " + std::to_string(r.dimension) + " " + r.name + "\n";
        out += std::string("$ MATERIAL ") + kindName(r.material.kind) + " " +
               (r.material.kind == MaterialKind::None || r.material.kind == MaterialKind::Unset
                    ? std::string("-") : r.material.record["name"].get<std::string>()) + "\n";
        if (r.dimension == 3) {
            out += card("PSOLID", {pid, pid});
            if (r.material.kind != MaterialKind::None) {
                const auto em = emData(r).dump();
                for (std::size_t i = 0; i < em.size(); i += 72) out += "$ EM " + em.substr(i, 72) + "\n";
                const auto m = mat4From(r);
                out += card("MAT4", {pid, real16(m.k), real16(m.cp * cpScale), real16(m.rho * rhoScale)});
            }
        } else {
            // A boundary or port face set: named faces, not a structural shell (no MID, no T).
            out += card("PSHELL", {pid});
        }
    }
    for (std::size_t i = 0; i < mesh.nodeIds.size(); ++i)
        out += card("GRID", {id16(mesh.nodeIds[i], "node"), "", real16(mesh.xyz[3 * i] * s),
                             real16(mesh.xyz[3 * i + 1] * s), real16(mesh.xyz[3 * i + 2] * s)});
    for (const auto& b : mesh.blocks) {
        const auto& ti = info(b.type);
        if (!ti.nastran)
            throw std::runtime_error("exportMesh(bdf): no Nastran card for element type " +
                                     std::to_string(int(b.type)));
        for (std::size_t e = 0; e < b.ids.size(); ++e) {
            std::vector<std::string> f{id16(b.ids[e], "element"), std::to_string(b.region)};
            for (int k = 0; k < ti.nodes; ++k) f.push_back(std::to_string(b.nodes[e * ti.nodes + k]));
            out += card(ti.nastran, f);
        }
    }
    out += "ENDDATA\n";
    return out;
}

Mesh readBdf(const std::string& text) {
    const char* who = "importMesh(bdf)";
    std::istringstream in(text);
    std::string line, unit;
    std::vector<std::vector<std::string>> cards;   // name, fields...
    std::map<int, Region> regions;
    std::vector<int> regionOrder;
    bool inBulk = false;
    auto fieldsOf = [](const std::string& l, std::vector<std::string>& f) {
        for (std::size_t at = 8; at < l.size(); at += 16) f.push_back(l.substr(at, 16));
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("$ LENGTH UNIT: ", 0) == 0) { unit = line.substr(15); continue; }
        if (line.rfind("$ REGION ", 0) == 0) {
            std::istringstream r(line.substr(9));
            Region reg;
            r >> reg.tag >> reg.dimension;
            std::getline(r >> std::ws, reg.name);
            regions[reg.tag] = reg;
            regionOrder.push_back(reg.tag);
            continue;
        }
        if (line.empty() || line[0] == '$') continue;
        if (line == "BEGIN BULK") { inBulk = true; continue; }
        if (line == "ENDDATA") break;
        if (!inBulk) throw std::runtime_error(std::string(who) + ": data before BEGIN BULK");
        if (line[0] == '*') {
            if (cards.empty()) throw std::runtime_error(std::string(who) + ": continuation with no card");
            fieldsOf(line, cards.back());
            continue;
        }
        auto name = line.substr(0, std::min<std::size_t>(8, line.size()));
        while (!name.empty() && name.back() == ' ') name.pop_back();
        if (name.empty() || name.back() != '*')
            throw std::runtime_error(std::string(who) + ": only large-field cards are read, got '" + name + "'");
        name.pop_back();
        cards.push_back({name});
        fieldsOf(line, cards.back());
    }
    if (unit.empty()) throw std::runtime_error(std::string(who) + ": no '$ LENGTH UNIT:' line");
    const double s = unitScale(unit, who);
    auto integer = [&](const std::string& f) {
        std::size_t used = 0;
        std::string t = f;
        while (!t.empty() && t.back() == ' ') t.pop_back();
        auto v = std::stoll(t, &used);
        if (used != t.size()) throw std::runtime_error(std::string(who) + ": bad integer '" + f + "'");
        return std::int64_t(v);
    };
    Mesh mesh;
    for (int tag : regionOrder) mesh.regions.push_back(regions[tag]);
    for (const auto& c : cards) {
        const auto& name = c[0];
        if (name == "GRID") {
            if (c.size() < 6) throw std::runtime_error(std::string(who) + ": short GRID");
            mesh.nodeIds.push_back(integer(c[1]));
            for (int k = 0; k < 3; ++k) mesh.xyz.push_back(readReal(c[3 + k], who) / s);
        } else if (name == "PSOLID" || name == "PSHELL" || name == "MAT4") {
            continue;
        } else {
            const auto& ti = infoFromNastran(name);
            if (c.size() < std::size_t(3 + ti.nodes))
                throw std::runtime_error(std::string(who) + ": short " + name);
            const int pid = int(integer(c[2]));
            if (mesh.blocks.empty() || mesh.blocks.back().type != ti.type || mesh.blocks.back().region != pid)
                mesh.blocks.push_back(ElementBlock{ti.type, pid, 0, {}, {}});
            auto& b = mesh.blocks.back();
            b.ids.push_back(integer(c[1]));
            for (int k = 0; k < ti.nodes; ++k) b.nodes.push_back(integer(c[3 + k]));
        }
    }
    return mesh;
}


// ---- Abaqus INP ------------------------------------------------------------------------------
// The mesh and its named element sets (one *ELEMENT block per mesh block, ELSET = the region's
// name), coordinates %.17g so they read back bit-exact. INP has no unit field: the unit is
// stated in a ** comment, as are each region's tag/dimension and each block's gmsh entity, so the
// reader rebuilds the mesh exactly. No *MATERIAL or *SECTION: the materials go in the
// .materials.json companion (exportMaterials), named per region in ** MATERIAL lines.

std::string abaqusName(const std::string& n) {
    for (char c : n)
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-'))
            return "\"" + n + "\"";
    return n;
}

std::string writeInp(const Mesh& mesh, const std::string& unit) {
    const auto regions = checkMesh(mesh, "exportMesh(inp)");
    const double s = unitScale(unit, "exportMesh(inp)");
    std::string out = "*HEADING\nMVB++ mvbpp_mesh mesh\n";
    out += "** LENGTH UNIT: " + unit + "\n";
    for (const auto& r : mesh.regions) {
        if (r.name.find('"') != std::string::npos)
            throw std::runtime_error("exportMesh(inp): region name '" + r.name + "' contains a quote");
        out += "** REGION " + std::to_string(r.tag) + " " + std::to_string(r.dimension) + " " + r.name + "\n";
        if (r.material.kind != MaterialKind::Unset && r.material.kind != MaterialKind::None)
            out += std::string("** MATERIAL ") + r.name + " = " + r.material.record["name"].get<std::string>() + "\n";
    }
    out += "*NODE\n";
    for (std::size_t i = 0; i < mesh.nodeIds.size(); ++i)
        out += std::to_string(mesh.nodeIds[i]) + ", " + g17(mesh.xyz[3 * i] * s) + ", " +
               g17(mesh.xyz[3 * i + 1] * s) + ", " + g17(mesh.xyz[3 * i + 2] * s) + "\n";
    for (const auto& b : mesh.blocks) {
        const auto& ti = info(b.type);
        if (!ti.abaqus)
            throw std::runtime_error("exportMesh(inp): no Abaqus element type written for element type " +
                                     std::to_string(int(b.type)));
        out += "** ENTITY " + std::to_string(b.entity) + "\n*ELEMENT, TYPE=" + ti.abaqus +
               ", ELSET=" + abaqusName(regions.at(b.region)->name) + "\n";
        for (std::size_t e = 0; e < b.ids.size(); ++e) {
            out += std::to_string(b.ids[e]);
            for (int k = 0; k < ti.nodes; ++k) out += ", " + std::to_string(b.nodes[e * ti.nodes + k]);
            out += "\n";
        }
    }
    return out;
}

std::vector<std::string> splitCommas(const std::string& line) {
    std::vector<std::string> f;
    std::string cur;
    bool quoted = false;
    for (char c : line) {
        if (c == '"') { quoted = !quoted; continue; }
        if (c == ',' && !quoted) { f.push_back(cur); cur.clear(); continue; }
        cur += c;
    }
    f.push_back(cur);
    for (auto& x : f) {
        while (!x.empty() && std::isspace(static_cast<unsigned char>(x.back()))) x.pop_back();
        while (!x.empty() && std::isspace(static_cast<unsigned char>(x.front()))) x.erase(0, 1);
    }
    return f;
}

Mesh readInp(const std::string& text) {
    const char* who = "importMesh(inp)";
    std::istringstream in(text);
    std::string line, unit;
    Mesh mesh;
    std::map<std::string, int> tagByName;
    enum { None, Nodes, Elements } section = None;
    const TypeInfo* ti = nullptr;
    int entity = 0;
    bool haveEntity = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("** LENGTH UNIT: ", 0) == 0) { unit = line.substr(16); continue; }
        if (line.rfind("** REGION ", 0) == 0) {
            std::istringstream r(line.substr(10));
            Region reg;
            r >> reg.tag >> reg.dimension;
            std::getline(r >> std::ws, reg.name);
            tagByName[reg.name] = reg.tag;
            mesh.regions.push_back(reg);
            continue;
        }
        if (line.rfind("** ENTITY ", 0) == 0) { entity = std::stoi(line.substr(10)); haveEntity = true; continue; }
        if (line.rfind("**", 0) == 0 || line.empty()) continue;
        if (line[0] == '*') {
            auto f = splitCommas(line);
            std::string kw = f[0];
            for (auto& c : kw) c = char(std::toupper(static_cast<unsigned char>(c)));
            if (kw == "*NODE") { section = Nodes; continue; }
            if (kw == "*ELEMENT") {
                std::string type, set;
                for (std::size_t i = 1; i < f.size(); ++i) {
                    auto eq = f[i].find('=');
                    std::string k = f[i].substr(0, eq), v = eq == std::string::npos ? "" : f[i].substr(eq + 1);
                    for (auto& c : k) c = char(std::toupper(static_cast<unsigned char>(c)));
                    if (k == "TYPE") type = v;
                    if (k == "ELSET") set = v;
                }
                ti = nullptr;
                for (const auto& t : kTypes)
                    if (t.abaqus && type == t.abaqus) ti = &t;
                if (!ti) throw std::runtime_error(std::string(who) + ": element type '" + type + "' is not one this model reads");
                auto it = tagByName.find(set);
                if (it == tagByName.end())
                    throw std::runtime_error(std::string(who) + ": ELSET '" + set + "' has no ** REGION line");
                if (!haveEntity) throw std::runtime_error(std::string(who) + ": *ELEMENT without its ** ENTITY line");
                mesh.blocks.push_back(ElementBlock{ti->type, it->second, entity, {}, {}});
                haveEntity = false;
                section = Elements;
                continue;
            }
            if (kw == "*HEADING") { std::getline(in, line); section = None; continue; }
            throw std::runtime_error(std::string(who) + ": keyword " + kw + " is not one this model reads");
        }
        auto f = splitCommas(line);
        if (section == Nodes) {
            if (f.size() != 4) throw std::runtime_error(std::string(who) + ": node line '" + line + "'");
            mesh.nodeIds.push_back(std::stoll(f[0]));
            for (int k = 1; k <= 3; ++k) mesh.xyz.push_back(std::stod(f[k]));
        } else if (section == Elements) {
            if (f.size() != std::size_t(1 + ti->nodes))
                throw std::runtime_error(std::string(who) + ": element line '" + line + "'");
            mesh.blocks.back().ids.push_back(std::stoll(f[0]));
            for (int k = 1; k <= ti->nodes; ++k) mesh.blocks.back().nodes.push_back(std::stoll(f[k]));
        } else {
            throw std::runtime_error(std::string(who) + ": data outside *NODE/*ELEMENT: '" + line + "'");
        }
    }
    if (unit.empty()) throw std::runtime_error(std::string(who) + ": no '** LENGTH UNIT:' line");
    const double s = unitScale(unit, who);
    if (s != 1.0) for (auto& x : mesh.xyz) x /= s;
    return mesh;
}

// ---- VTK legacy (ASCII unstructured grid) ----------------------------------------------------
// ParaView-readable. Coordinates %.17g; gmsh's node and element ids, each cell's region tag and
// gmsh entity travel as POINT_DATA/CELL_DATA arrays, and each region as a field array
// "region_<tag>" holding [dimension, character codes of the name] (legacy VTK has no string
// arrays), so the reader rebuilds the mesh exactly. The unit is in the title line.

std::string writeVtk(const Mesh& mesh, const std::string& unit) {
    checkMesh(mesh, "exportMesh(vtk)");
    const double s = unitScale(unit, "exportMesh(vtk)");
    std::map<std::int64_t, std::size_t> index;
    for (std::size_t i = 0; i < mesh.nodeIds.size(); ++i) index[mesh.nodeIds[i]] = i;
    std::string out = "# vtk DataFile Version 3.0\nMVB++ mvbpp_mesh mesh, LENGTH UNIT: " + unit +
                      "\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out += "FIELD FieldData " + std::to_string(mesh.regions.size()) + "\n";
    for (const auto& r : mesh.regions) {
        out += "region_" + std::to_string(r.tag) + " 1 " + std::to_string(r.name.size() + 1) + " int\n" +
               std::to_string(r.dimension);
        for (unsigned char c : r.name) out += " " + std::to_string(int(c));
        out += "\n";
    }
    out += "POINTS " + std::to_string(mesh.nodeIds.size()) + " double\n";
    for (std::size_t i = 0; i < mesh.nodeIds.size(); ++i)
        out += g17(mesh.xyz[3 * i] * s) + " " + g17(mesh.xyz[3 * i + 1] * s) + " " + g17(mesh.xyz[3 * i + 2] * s) + "\n";
    std::size_t cells = 0, size = 0;
    for (const auto& b : mesh.blocks) { cells += b.ids.size(); size += b.ids.size() * (1 + info(b.type).nodes); }
    out += "CELLS " + std::to_string(cells) + " " + std::to_string(size) + "\n";
    for (const auto& b : mesh.blocks) {
        const auto& ti = info(b.type);
        for (std::size_t e = 0; e < b.ids.size(); ++e) {
            out += std::to_string(ti.nodes);
            for (int k = 0; k < ti.nodes; ++k)
                out += " " + std::to_string(index.at(b.nodes[e * ti.nodes + ti.vtkOrder[k]]));
            out += "\n";
        }
    }
    out += "CELL_TYPES " + std::to_string(cells) + "\n";
    for (const auto& b : mesh.blocks)
        for (std::size_t e = 0; e < b.ids.size(); ++e) out += std::to_string(info(b.type).vtk) + "\n";
    out += "CELL_DATA " + std::to_string(cells) + "\n";
    auto cellArray = [&](const char* name, auto value) {
        out += std::string("SCALARS ") + name + " long 1\nLOOKUP_TABLE default\n";
        for (const auto& b : mesh.blocks)
            for (std::size_t e = 0; e < b.ids.size(); ++e) out += std::to_string(value(b, e)) + "\n";
    };
    cellArray("region", [](const ElementBlock& b, std::size_t) { return std::int64_t(b.region); });
    cellArray("entity", [](const ElementBlock& b, std::size_t) { return std::int64_t(b.entity); });
    cellArray("elementId", [](const ElementBlock& b, std::size_t e) { return b.ids[e]; });
    out += "POINT_DATA " + std::to_string(mesh.nodeIds.size()) + "\nSCALARS nodeId long 1\nLOOKUP_TABLE default\n";
    for (auto id : mesh.nodeIds) out += std::to_string(id) + "\n";
    return out;
}

Mesh readVtk(const std::string& text) {
    const char* who = "importMesh(vtk)";
    std::istringstream in(text);
    std::string line;
    std::getline(in, line);
    if (line.rfind("# vtk DataFile", 0) != 0) throw std::runtime_error(std::string(who) + ": not a legacy VTK file");
    std::getline(in, line);
    const auto at = line.find("LENGTH UNIT: ");
    if (at == std::string::npos) throw std::runtime_error(std::string(who) + ": the title states no LENGTH UNIT");
    const double s = unitScale(line.substr(at + 13), who);
    std::string word;
    in >> word;
    if (word != "ASCII") throw std::runtime_error(std::string(who) + ": only ASCII files are read");
    Mesh mesh;
    std::vector<std::vector<std::int64_t>> conn;
    std::vector<int> types;
    std::map<std::string, std::vector<std::int64_t>> cellData;
    std::vector<std::int64_t> nodeIds;
    std::string mode;
    while (in >> word) {
        if (word == "DATASET") { in >> word; if (word != "UNSTRUCTURED_GRID") throw std::runtime_error(std::string(who) + ": not an unstructured grid"); }
        else if (word == "FIELD") {
            std::string name; std::size_t n; in >> name >> n;
            for (std::size_t i = 0; i < n; ++i) {
                std::string arr, type; std::size_t comps, tuples; in >> arr >> comps >> tuples >> type;
                if (arr.rfind("region_", 0) != 0 || comps != 1 || type != "int" || tuples < 1)
                    throw std::runtime_error(std::string(who) + ": unexpected field array " + arr);
                Region r; r.tag = std::stoi(arr.substr(7)); in >> r.dimension;
                for (std::size_t k = 1; k < tuples; ++k) { int c; in >> c; r.name += char(c); }
                mesh.regions.push_back(r);
            }
        } else if (word == "POINTS") {
            std::size_t n; std::string type; in >> n >> type;
            mesh.xyz.resize(3 * n);
            for (auto& x : mesh.xyz) { in >> x; x /= s; }
        } else if (word == "CELLS") {
            std::size_t n, size; in >> n >> size;
            conn.resize(n);
            for (auto& c : conn) { std::size_t k; in >> k; c.resize(k); for (auto& v : c) in >> v; }
        } else if (word == "CELL_TYPES") {
            std::size_t n; in >> n; types.resize(n);
            for (auto& t : types) in >> t;
        } else if (word == "CELL_DATA" || word == "POINT_DATA") {
            mode = word; std::size_t n; in >> n;
        } else if (word == "SCALARS") {
            std::string name, type; int comps; in >> name >> type >> comps >> word >> word;  // LOOKUP_TABLE default
            std::vector<std::int64_t> v(mode == "CELL_DATA" ? conn.size() : mesh.xyz.size() / 3);
            for (auto& x : v) in >> x;
            if (mode == "POINT_DATA" && name == "nodeId") nodeIds = v;
            else if (mode == "CELL_DATA") cellData[name] = v;
            else throw std::runtime_error(std::string(who) + ": unexpected point array " + name);
        } else {
            throw std::runtime_error(std::string(who) + ": unexpected keyword " + word);
        }
        if (!in) throw std::runtime_error(std::string(who) + ": truncated at " + word);
    }
    for (const char* a : {"region", "entity", "elementId"})
        if (!cellData.count(a)) throw std::runtime_error(std::string(who) + ": no cell array " + a);
    if (nodeIds.size() != mesh.xyz.size() / 3) throw std::runtime_error(std::string(who) + ": no nodeId array");
    mesh.nodeIds = nodeIds;
    for (std::size_t c = 0; c < conn.size(); ++c) {
        const TypeInfo* ti = nullptr;
        for (const auto& t : kTypes) if (t.vtk == types[c]) ti = &t;
        if (!ti || conn[c].size() != std::size_t(ti->nodes))
            throw std::runtime_error(std::string(who) + ": cell " + std::to_string(c) + " has VTK type " + std::to_string(types[c]));
        const int region = int(cellData["region"][c]), entity = int(cellData["entity"][c]);
        if (mesh.blocks.empty() || mesh.blocks.back().type != ti->type || mesh.blocks.back().region != region ||
            mesh.blocks.back().entity != entity)
            mesh.blocks.push_back(ElementBlock{ti->type, region, entity, {}, {}});
        auto& b = mesh.blocks.back();
        b.ids.push_back(cellData["elementId"][c]);
        std::array<std::int64_t, 8> nodes{};
        for (int k = 0; k < ti->nodes; ++k) nodes[ti->vtkOrder[k]] = nodeIds.at(std::size_t(conn[c][k]));
        for (int k = 0; k < ti->nodes; ++k) b.nodes.push_back(nodes[k]);
    }
    return mesh;
}

}  // namespace

std::string exportMesh(const Mesh& mesh, const std::string& format, const std::string& unit) {
    if (format == "msh2") return writeMsh2(mesh, unit);
    if (format == "bdf") return writeBdf(mesh, unit);
    if (format == "inp") return writeInp(mesh, unit);
    if (format == "vtk") return writeVtk(mesh, unit);
    throw std::runtime_error("exportMesh: unknown format \"" + format + "\" (msh2, bdf, inp, vtk)");
}

nlohmann::json exportMaterials(const Mesh& mesh) {
    checkMesh(mesh, "exportMaterials");
    nlohmann::json out = nlohmann::json::array();
    for (const auto& r : mesh.regions) {
        requireMaterial(r, "exportMaterials");
        nlohmann::json e = {{"pid", r.tag}, {"name", r.name}, {"dimension", r.dimension},
                            {"kind", kindName(r.material.kind)}};
        if (r.material.kind != MaterialKind::None && r.material.kind != MaterialKind::Unset) {
            e["material"] = r.material.record["name"];
            e["electromagnetic"] = emData(r);
            nlohmann::json th = nlohmann::json::object();
            for (const char* key : {"heatConductivity", "heatCapacity", "thermalConductivity",
                                    "specificHeat", "density"})
                if (r.material.record.contains(key)) th[key] = r.material.record[key];
            e["thermal"] = th;
            if (r.material.record.contains("provenance")) e["provenance"] = r.material.record["provenance"];
        }
        out.push_back(e);
    }
    return out;
}

Mesh importMesh(const std::string& text, const std::string& format) {
    if (format == "msh2") return readMsh2(text);
    if (format == "bdf") return readBdf(text);
    if (format == "inp") return readInp(text);
    if (format == "vtk") return readVtk(text);
    throw std::runtime_error("importMesh: unknown format \"" + format + "\" (msh2, bdf, inp, vtk)");
}

std::string diffMeshes(const Mesh& a, const Mesh& b) {
    if (a.lengthUnit != b.lengthUnit) return "length unit " + a.lengthUnit + " vs " + b.lengthUnit;
    if (a.nodeIds.size() != b.nodeIds.size())
        return "node count " + std::to_string(a.nodeIds.size()) + " vs " + std::to_string(b.nodeIds.size());
    for (std::size_t i = 0; i < a.nodeIds.size(); ++i) {
        if (a.nodeIds[i] != b.nodeIds[i])
            return "node #" + std::to_string(i) + " id " + std::to_string(a.nodeIds[i]) + " vs " +
                   std::to_string(b.nodeIds[i]);
        for (int k = 0; k < 3; ++k)
            if (a.xyz[3 * i + k] != b.xyz[3 * i + k])
                return "node " + std::to_string(a.nodeIds[i]) + " coordinate " + std::to_string(k) + " " +
                       g17(a.xyz[3 * i + k]) + " vs " + g17(b.xyz[3 * i + k]);
    }
    if (a.regions.size() != b.regions.size())
        return "region count " + std::to_string(a.regions.size()) + " vs " + std::to_string(b.regions.size());
    for (std::size_t i = 0; i < a.regions.size(); ++i) {
        const auto &x = a.regions[i], &y = b.regions[i];
        if (x.tag != y.tag || x.dimension != y.dimension || x.name != y.name)
            return "region #" + std::to_string(i) + " " + std::to_string(x.dimension) + "/" +
                   std::to_string(x.tag) + " '" + x.name + "' vs " + std::to_string(y.dimension) + "/" +
                   std::to_string(y.tag) + " '" + y.name + "'";
    }
    if (a.blocks.size() != b.blocks.size())
        return "block count " + std::to_string(a.blocks.size()) + " vs " + std::to_string(b.blocks.size());
    for (std::size_t i = 0; i < a.blocks.size(); ++i) {
        const auto &x = a.blocks[i], &y = b.blocks[i];
        const std::string at = "block #" + std::to_string(i) + ": ";
        if (x.type != y.type) return at + "type";
        if (x.region != y.region) return at + "region " + std::to_string(x.region) + " vs " + std::to_string(y.region);
        if (x.entity != y.entity) return at + "entity " + std::to_string(x.entity) + " vs " + std::to_string(y.entity);
        if (x.ids != y.ids) return at + "element ids";
        if (x.nodes != y.nodes) return at + "element nodes";
    }
    return "";
}

}  // namespace mvb::mesh
