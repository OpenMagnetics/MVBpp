#include "mvb/Symmetry.h"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace mvb {

const char* to_string(SymmetryPlane p) {
    switch (p) {
        case SymmetryPlane::None: return "None";
        case SymmetryPlane::XY:   return "XY";
        case SymmetryPlane::YZ:   return "YZ";
        case SymmetryPlane::XZ:   return "XZ";
    }
    return "?";
}

void plane_axes(SymmetryPlane p, int& i_normal, int& i_in_a, int& i_in_b) {
    switch (p) {
        case SymmetryPlane::XY: i_normal = 2; i_in_a = 0; i_in_b = 1; return;
        case SymmetryPlane::YZ: i_normal = 0; i_in_a = 1; i_in_b = 2; return;
        case SymmetryPlane::XZ: i_normal = 1; i_in_a = 0; i_in_b = 2; return;
        case SymmetryPlane::None:
        default:                i_normal = -1; i_in_a = -1; i_in_b = -1;
    }
}

namespace {

double shape_volume(const TopoDS_Shape& s) {
    if (s.IsNull()) return 0.0;
    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props);
    return std::max(props.Mass(), 0.0);
}

// Return the axis index normal to a symmetry plane (the symmetry planes all
// pass through the origin so a single coordinate is enough to classify a
// point or AABB against the plane). XY -> z, YZ -> x, XZ -> y.
int plane_normal_axis(SymmetryPlane p) {
    switch (p) {
        case SymmetryPlane::XY: return 2;
        case SymmetryPlane::YZ: return 0;
        case SymmetryPlane::XZ: return 1;
        case SymmetryPlane::None:
        default:                return -1;
    }
}

// Cheap AABB extents along a single axis, or {NaN, NaN} for a null shape.
std::pair<double, double> shape_axis_extent(const TopoDS_Shape& s, int axis) {
    if (s.IsNull() || axis < 0) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }
    Bnd_Box bb;
    BRepBndLib::Add(s, bb);
    if (bb.IsVoid()) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }
    double xmin, ymin, zmin, xmax, ymax, zmax;
    bb.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    const double mins[3] = {xmin, ymin, zmin};
    const double maxs[3] = {xmax, ymax, zmax};
    return {mins[axis], maxs[axis]};
}

// Side classification of a shape's AABB against a plane through the origin
// on the given axis. +1: fully on positive side. -1: fully on negative side.
// 0: straddles the plane (needs a boolean Common to measure the half-volume).
// The eps absorbs OCCT bbox padding noise: a shape whose bbox sits in [eps,
// large] is treated as fully positive (the bbox is overestimated already).
int classify_axis(double bb_min, double bb_max, double eps = 1e-9) {
    if (bb_min >= -eps) return +1;
    if (bb_max <=  eps) return -1;
    return 0;
}

// Build a solid axis-aligned box occupying one half of the bounding box.
TopoDS_Shape make_halfspace(SymmetryPlane plane, SymmetryHalf half,
                            const ShapeBBox& b) {
    const double dx = b.xmax - b.xmin;
    const double dy = b.ymax - b.ymin;
    const double dz = b.zmax - b.zmin;
    const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
    const double pad = std::max(diag * 2.0, 1e-3);
    const double cx = 0.5 * (b.xmin + b.xmax);
    const double cy = 0.5 * (b.ymin + b.ymax);
    const double cz = 0.5 * (b.zmin + b.zmax);

    double lx_min = cx - pad, lx_max = cx + pad;
    double ly_min = cy - pad, ly_max = cy + pad;
    double lz_min = cz - pad, lz_max = cz + pad;

    switch (plane) {
        case SymmetryPlane::YZ:
            if (half == SymmetryHalf::Positive) lx_min = 0.0;
            else                                 lx_max = 0.0;
            break;
        case SymmetryPlane::XZ:
            if (half == SymmetryHalf::Positive) ly_min = 0.0;
            else                                 ly_max = 0.0;
            break;
        case SymmetryPlane::XY:
            if (half == SymmetryHalf::Positive) lz_min = 0.0;
            else                                 lz_max = 0.0;
            break;
        case SymmetryPlane::None:
            break;
    }

    BRepPrimAPI_MakeBox mk(gp_Pnt(lx_min, ly_min, lz_min),
                           gp_Pnt(lx_max, ly_max, lz_max));
    return mk.Solid();
}

double volume_on_side(const std::vector<NamedShape>& shapes,
                      const TopoDS_Shape& halfspace) {
    // Kept for binary compatibility / potential future callers. The analyze
    // path now does AABB-prefiltered per-shape ops directly.
    double v = 0.0;
    for (const auto& ns : shapes) {
        if (ns.shape.IsNull()) continue;
        BRepAlgoAPI_Common op(ns.shape, halfspace);
        op.Build();
        if (!op.IsDone()) continue;
        v += shape_volume(op.Shape());
    }
    return v;
}

// Enumerate all solid sub-shapes of a (possibly compound) result. Returns
// each as a separate TopoDS_Shape so consumers can name them individually.
std::vector<TopoDS_Shape> flatten_solids(const TopoDS_Shape& s) {
    std::vector<TopoDS_Shape> out;
    if (s.IsNull()) return out;
    if (s.ShapeType() == TopAbs_SOLID) {
        out.push_back(s);
        return out;
    }
    for (TopExp_Explorer ex(s, TopAbs_SOLID); ex.More(); ex.Next()) {
        out.push_back(ex.Current());
    }
    return out;
}

} // namespace

ShapeBBox aggregate_bbox(const std::vector<NamedShape>& shapes) {
    Bnd_Box box;
    for (const auto& ns : shapes) {
        if (!ns.shape.IsNull()) BRepBndLib::Add(ns.shape, box);
    }
    ShapeBBox out{0, 0, 0, 0, 0, 0};
    if (!box.IsVoid()) {
        box.Get(out.xmin, out.ymin, out.zmin, out.xmax, out.ymax, out.zmax);
    }
    return out;
}

SymmetryResult analyze_symmetry(const std::vector<NamedShape>& shapes,
                                double volume_tolerance) {
    SymmetryResult out;

    // Precompute per-shape volume and per-axis AABB extents once. The
    // original implementation ran 6 boolean Common ops per plane (positive
    // and negative halfspace per shape) and recomputed shape volumes from
    // scratch. Most ETD/EP/PQ core halves and bobbin flanges sit fully on
    // one side of every symmetry plane, so an AABB pre-filter eliminates
    // the booleans entirely for those shapes.
    struct PerShape {
        const TopoDS_Shape* shape = nullptr;
        double V = 0.0;
        double bbmin[3] = {0, 0, 0};
        double bbmax[3] = {0, 0, 0};
        bool   bb_valid = false;
    };
    std::vector<PerShape> per; per.reserve(shapes.size());
    double V_total = 0.0;
    for (const auto& ns : shapes) {
        PerShape ps;
        ps.shape = &ns.shape;
        ps.V     = shape_volume(ns.shape);
        if (!ns.shape.IsNull()) {
            Bnd_Box bb;
            BRepBndLib::Add(ns.shape, bb);
            if (!bb.IsVoid()) {
                bb.Get(ps.bbmin[0], ps.bbmin[1], ps.bbmin[2],
                       ps.bbmax[0], ps.bbmax[1], ps.bbmax[2]);
                ps.bb_valid = true;
            }
        }
        V_total += ps.V;
        per.push_back(ps);
    }
    if (V_total <= 0.0) return out;

    const ShapeBBox bb = aggregate_bbox(shapes);

    const std::array<SymmetryPlane, 3> planes = {
        SymmetryPlane::XY, SymmetryPlane::YZ, SymmetryPlane::XZ};

    for (SymmetryPlane p : planes) {
        SymmetryCandidate c;
        c.plane = p;
        const int axis = plane_normal_axis(p);
        TopoDS_Shape pos_halfspace;  // built lazily on first straddling shape
        bool pos_halfspace_built = false;

        double V_pos = 0.0;
        for (const auto& ps : per) {
            if (ps.V <= 0.0 || !ps.bb_valid) continue;
            int side = classify_axis(ps.bbmin[axis], ps.bbmax[axis]);
            if (side > 0) { V_pos += ps.V; continue; }
            if (side < 0) { continue; }
            if (!pos_halfspace_built) {
                pos_halfspace = make_halfspace(p, SymmetryHalf::Positive, bb);
                pos_halfspace_built = true;
            }
            BRepAlgoAPI_Common op(*ps.shape, pos_halfspace);
            op.Build();
            if (!op.IsDone()) continue;
            V_pos += shape_volume(op.Shape());
        }
        c.volume_plus  = V_pos;
        c.volume_minus = V_total - V_pos;  // each shape's vol splits into the two halves
        c.volume_asymmetry =
            std::abs(c.volume_plus - c.volume_minus) /
            std::max(V_total, 1e-30);
        c.valid = c.volume_asymmetry < volume_tolerance;
        out.candidates.push_back(c);
        if (c.valid) out.valid_planes.push_back(p);
    }
    return out;
}

std::vector<NamedShape> cut_half_space(const NamedShape& in,
                                       SymmetryPlane plane,
                                       SymmetryHalf  half,
                                       const ShapeBBox& bbox) {
    std::vector<NamedShape> out;
    if (plane == SymmetryPlane::None || in.shape.IsNull()) {
        out.push_back(in);
        return out;
    }
    // AABB short-circuit: if the shape is entirely on the keep side, copy it
    // through; if entirely on the discard side, drop it. Same logic as the
    // analyze_symmetry pre-filter; saves boolean Common on every untouched
    // core half / bobbin flange / outer turn.
    const int axis = plane_normal_axis(plane);
    auto extent = shape_axis_extent(in.shape, axis);
    if (!std::isnan(extent.first)) {
        int side = classify_axis(extent.first, extent.second);
        const int keepSign = (half == SymmetryHalf::Positive) ? +1 : -1;
        if (side == keepSign)        { out.push_back(in); return out; }
        if (side == -keepSign)       { return out; }
        // side == 0: shape straddles, fall through to boolean cut
    }
    TopoDS_Shape halfspace = make_halfspace(plane, half, bbox);
    BRepAlgoAPI_Common op(in.shape, halfspace);
    op.Build();
    if (!op.IsDone()) return out;
    const TopoDS_Shape result = op.Shape();
    const auto solids = flatten_solids(result);
    if (solids.empty()) return out;  // whole shape was on the other side

    if (solids.size() == 1) {
        out.emplace_back(solids.front(), in.name);
    } else {
        for (std::size_t i = 0; i < solids.size(); ++i) {
            out.emplace_back(solids[i],
                             in.name + "_" + std::to_string(i));
        }
    }
    return out;
}

std::vector<NamedShape> cut_half_space(const std::vector<NamedShape>& in,
                                       SymmetryPlane plane,
                                       SymmetryHalf  half,
                                       const ShapeBBox& bbox) {
    std::vector<NamedShape> out;
    out.reserve(in.size());
    for (const auto& ns : in) {
        auto r = cut_half_space(ns, plane, half, bbox);
        out.insert(out.end(),
                   std::make_move_iterator(r.begin()),
                   std::make_move_iterator(r.end()));
    }
    return out;
}

std::vector<NamedShape> cut_to_region(
    const std::vector<NamedShape>& in,
    const std::vector<std::pair<SymmetryPlane, SymmetryHalf>>& cuts,
    const ShapeBBox& bbox) {
    std::vector<NamedShape> current = in;
    for (const auto& [plane, half] : cuts) {
        current = cut_half_space(current, plane, half, bbox);
        if (current.empty()) break;
    }
    return current;
}

// ---- High-level helpers --------------------------------------------------

namespace {

// Centroid of a shape via mass properties — works for solids, faces, and
// edges (BRepGProp picks the right integral based on dimensionality).
gp_Pnt shape_centroid(const TopoDS_Shape& s) {
    GProp_GProps props;
    if (TopExp_Explorer(s, TopAbs_SOLID).More()) {
        BRepGProp::VolumeProperties(s, props);
    } else if (TopExp_Explorer(s, TopAbs_FACE).More()) {
        BRepGProp::SurfaceProperties(s, props);
    } else {
        BRepGProp::LinearProperties(s, props);
    }
    return props.CentreOfMass();
}

} // namespace

std::vector<NamedShape> filter_by_side(const std::vector<NamedShape>& in,
                                       const std::array<int, 3>& axisSign,
                                       double eps) {
    if (axisSign[0] == 0 && axisSign[1] == 0 && axisSign[2] == 0) return in;
    std::vector<NamedShape> out;
    out.reserve(in.size());
    for (const auto& ns : in) {
        if (ns.shape.IsNull()) continue;
        gp_Pnt c = shape_centroid(ns.shape);
        const double coord[3] = { c.X(), c.Y(), c.Z() };
        bool keep = true;
        for (int i = 0; i < 3; ++i) {
            if (axisSign[i] == 0) continue;
            // |coord| within eps of plane → on the plane, keep it.
            if (std::abs(coord[i]) <= eps) continue;
            const double signed_v = coord[i] * axisSign[i];
            if (signed_v < 0) { keep = false; break; }
        }
        if (keep) out.push_back(ns);
    }
    return out;
}

std::array<int, 3> parse_side_spec(const std::string& spec) {
    std::array<int, 3> out{0, 0, 0};
    std::string s = spec;
    // strip whitespace + lowercase
    s.erase(std::remove_if(s.begin(), s.end(),
                            [](unsigned char c) { return std::isspace(c); }),
            s.end());
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (lower.empty() || lower == "auto" || lower == "none") return out;

    int sign = 0;
    for (char c : s) {
        if (c == '+') { sign = +1; continue; }
        if (c == '-') { sign = -1; continue; }
        int axis;
        switch (c) {
            case 'X': case 'x': axis = 0; break;
            case 'Y': case 'y': axis = 1; break;
            case 'Z': case 'z': axis = 2; break;
            default:
                throw std::invalid_argument(
                    "parse_side_spec: bad token '" + std::string(1, c) +
                    "' in '" + spec + "' (expected +X/-X/+Y/-Y/+Z/-Z)");
        }
        if (sign == 0) {
            throw std::invalid_argument(
                "parse_side_spec: missing +/- before axis in '" + spec + "'");
        }
        out[axis] = sign;
        sign = 0;
    }
    return out;
}

std::vector<NamedShape> apply_symmetry(const std::vector<NamedShape>& in,
                                       int nPlanes) {
    if (nPlanes <= 0 || in.empty()) return in;
    if (nPlanes > 3) nPlanes = 3;
    SymmetryResult res = analyze_symmetry(in);
    if (res.valid_planes.empty()) return in;
    std::vector<std::pair<SymmetryPlane, SymmetryHalf>> cuts;
    for (int i = 0; i < nPlanes && i < static_cast<int>(res.valid_planes.size()); ++i) {
        cuts.emplace_back(res.valid_planes[i], SymmetryHalf::Positive);
    }
    if (cuts.empty()) return in;
    ShapeBBox bb = aggregate_bbox(in);
    return cut_to_region(in, cuts, bb);
}

int parse_symmetry_spec(const std::string& spec) {
    std::string s = spec;
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (s.empty() || s == "auto" || s == "none" || s == "0") return 0;
    if (s == "half"    || s == "1") return 1;
    if (s == "quarter" || s == "2") return 2;
    if (s == "eighth"  || s == "3") return 3;
    // Try numeric
    try {
        int n = std::stoi(spec);
        if (n < 0 || n > 3) {
            throw std::invalid_argument(
                "parse_symmetry_spec: out of range " + spec +
                " (expected 0..3 or auto/half/quarter/eighth)");
        }
        return n;
    } catch (const std::exception&) {
        throw std::invalid_argument(
            "parse_symmetry_spec: bad value '" + spec +
            "' (expected auto/half/quarter/eighth or 0..3)");
    }
}

} // namespace mvb
