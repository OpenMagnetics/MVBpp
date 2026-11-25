// A conductor's display drawing triangulated straight from its path (ABT #2489).
// Building every piece of the wire as an exact CAD solid (sweep, mitre cut, validity check) and then tessellating those
// solids is what a real-winding drawing spent its time on: 345 s over the 38 corpus designs against 121 s drawn this way,
// with a quarter of the triangles and a rounder wire. This file draws the tube that follows the path as triangles, one
// triangulated face per piece.
#include "mvb/DisplayTube.h"

#include <BRep_Builder.hxx>
#include <Poly_Triangulation.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Vec.hxx>

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

namespace mvb {
namespace {

struct Station {
    gp_XYZ point;
    gp_XYZ tangent;   // unit
};

// The point and the unit tangent of a path piece at its parameter t in [0, 1]. The same closed forms the assembler
// sweeps and the certifier measures (cert::evalPrim); the spiral's profile comes from its one definition.
Station stationOf(const Primitive& pr, double t) {
    Station s;
    gp_XYZ rate;
    switch (pr.kind) {
        case Primitive::SEG: {
            rate = pr.seg.b.XYZ() - pr.seg.a.XYZ();
            s.point = pr.seg.a.XYZ() + rate * t;
            break;
        }
        case Primitive::ARC3: {
            const gp_XYZ v = rotateXYZ(pr.arc.v0, pr.arc.axis, pr.arc.sweep * t);
            s.point = pr.arc.c.XYZ() + v;
            rate = pr.arc.axis.Crossed(v) * pr.arc.sweep;
            break;
        }
        case Primitive::SPIRAL: {
            const Spiral& sp = pr.spiral;
            const double sweep = sp.az1 - sp.az0;
            const double az = sp.az0 + sweep * t;
            const double r = sp.r0 + (sp.r1 - sp.r0) * spiralRadiusFrac(sp, t);
            const double y = sp.y0 + (sp.y1 - sp.y0) * spiralHeightFrac(sp, t);
            const double rRate = (sp.r1 - sp.r0) * spiralRadiusFracRate(sp, t);
            const double yRate = (sp.y1 - sp.y0) * spiralHeightFracRate(sp, t);
            s.point = gp_XYZ(sp.cx + r * std::cos(az), y, sp.cz - r * std::sin(az));
            rate = gp_XYZ(rRate * std::cos(az) - r * std::sin(az) * sweep, yRate,
                          -rRate * std::sin(az) - r * std::cos(az) * sweep);
            break;
        }
        case Primitive::BLEND: {
            const gp_XYZ d = pr.blendc.b.XYZ() - pr.blendc.a.XYZ();
            const double along = d.Dot(pr.blendc.u);
            const gp_XYZ across = d - pr.blendc.u * along;
            s.point = pr.blendc.a.XYZ() + pr.blendc.u * (along * t) + across * (0.5 * (1.0 - std::cos(std::numbers::pi * t)));
            rate = pr.blendc.u * along + across * (0.5 * std::numbers::pi * std::sin(std::numbers::pi * t));
            break;
        }
        default:
            throw std::runtime_error("displayTube: unsupported path piece kind");
    }
    const double speed = rate.Modulus();
    if (!(speed > 0.0))
        throw std::runtime_error("displayTube: the path piece '" + pr.label + "' has no direction at t = " + std::to_string(t));
    s.tangent = rate / speed;
    return s;
}

// How many chords a piece is drawn with: its turning angle over the angular step, as a tessellator would choose.
int chordsOf(const Primitive& pr, double angularStep) {
    switch (pr.kind) {
        case Primitive::SEG:
            return 1;
        case Primitive::ARC3:
            return std::max(1, static_cast<int>(std::ceil(std::abs(pr.arc.sweep) / angularStep)));
        case Primitive::SPIRAL: {
            const int turning = static_cast<int>(std::ceil(std::abs(pr.spiral.az1 - pr.spiral.az0) / angularStep));
            const bool shaped = pr.spiral.blend || pr.spiral.levelOut || pr.spiral.easeIn > 0.0;
            return std::max(shaped ? 8 : 1, turning);
        }
        case Primitive::BLEND:
            return 12;
        default:
            throw std::runtime_error("displayTube: unsupported path piece kind");
    }
}

// v turned by the least rotation that takes the unit direction `from` to the unit direction `to`.
gp_XYZ carried(const gp_XYZ& v, const gp_XYZ& from, const gp_XYZ& to) {
    const gp_XYZ axis = from.Crossed(to);
    const double sine = axis.Modulus();
    const double cosine = from.Dot(to);
    if (sine < 1e-14) {
        if (cosine < 0.0) throw std::runtime_error("displayTube: the path reverses on itself");
        return v;
    }
    return rotateXYZ(v, axis / sine, std::atan2(sine, cosine));
}

} // namespace

std::vector<gp_XY> roundTubeOutline(double radius, double angularStep) {
    if (!(radius > 0.0) || !(angularStep > 0.0)) throw std::runtime_error("roundTubeOutline: radius and angular step must be positive");
    const int sides = std::max(6, static_cast<int>(std::ceil(2.0 * std::numbers::pi / angularStep)));
    std::vector<gp_XY> outline(sides);
    for (int k = 0; k < sides; ++k) {
        const double angle = 2.0 * std::numbers::pi * k / sides;
        outline[k] = gp_XY(radius * std::cos(angle), radius * std::sin(angle));
    }
    return outline;
}

std::vector<gp_XY> rectTubeOutline(double width, double height, double cornerRadius, double angularStep) {
    const double halfWidth = 0.5 * width, halfHeight = 0.5 * height;
    if (!(width > 0.0) || !(height > 0.0) || !(angularStep > 0.0))
        throw std::runtime_error("rectTubeOutline: width, height and angular step must be positive");
    if (!(cornerRadius >= 0.0) || cornerRadius > std::min(halfWidth, halfHeight) + 1e-12)
        throw std::runtime_error("rectTubeOutline: corner radius " + std::to_string(cornerRadius * 1e3) + " mm does not fit a " +
                                 std::to_string(width * 1e3) + " x " + std::to_string(height * 1e3) + " mm section");
    if (cornerRadius == 0.0)
        return {gp_XY(halfWidth, -halfHeight), gp_XY(halfWidth, halfHeight), gp_XY(-halfWidth, halfHeight), gp_XY(-halfWidth, -halfHeight)};
    // Four quarter circles, counter-clockwise from the one at (+width, -height); a side the corners use up is not drawn.
    const int sides = std::max(2, static_cast<int>(std::ceil(0.5 * std::numbers::pi / angularStep)));
    const double innerWidth = std::max(halfWidth - cornerRadius, 0.0), innerHeight = std::max(halfHeight - cornerRadius, 0.0);
    const double centres[4][2] = {{innerWidth, -innerHeight}, {innerWidth, innerHeight}, {-innerWidth, innerHeight}, {-innerWidth, -innerHeight}};
    std::vector<gp_XY> outline;
    for (int corner = 0; corner < 4; ++corner)
        for (int k = 0; k <= sides; ++k) {
            const double angle = (corner - 1 + static_cast<double>(k) / sides) * 0.5 * std::numbers::pi;
            const gp_XY point(centres[corner][0] + cornerRadius * std::cos(angle), centres[corner][1] + cornerRadius * std::sin(angle));
            if (outline.empty() || (point - outline.back()).Modulus() > 1e-12) outline.push_back(point);
        }
    if (outline.size() > 1 && (outline.front() - outline.back()).Modulus() <= 1e-12) outline.pop_back();
    return outline;
}

TopoDS_Shape displayTube(const std::vector<TubePiece>& path, const std::vector<gp_XY>& outline, double angularStep) {
    if (path.empty()) return TopoDS_Shape();
    if (outline.size() < 3 || !(angularStep > 0.0)) throw std::runtime_error("displayTube: an outline of at least three points and a positive angular step are needed");
    const int around = static_cast<int>(outline.size());
    double reach = 0.0;   // the outline's farthest point from the centreline
    for (const auto& q : outline) reach = std::max(reach, q.Modulus());

    // The frame (normal, tangent x normal) is carried along the whole conductor, so that consecutive pieces share
    // their junction ring. An oriented piece sets it at every station; the carried one settles its sign.
    Station first = stationOf(path.front().prim, 0.0);
    gp_XYZ normal = std::abs(first.tangent.Y()) < 0.9 ? gp_XYZ(0, 1, 0) : gp_XYZ(1, 0, 0);
    normal = normal - first.tangent * normal.Dot(first.tangent);
    normal.Normalize();
    gp_XYZ previousTangent = first.tangent;
    bool settled = false;   // an oriented frame has been taken: its sign is kept from then on

    // The section's width axis where the piece says how its section lies: tangent x axial, as rectProfileWire has it.
    // Where the path runs along the axial direction itself there is no such axis and the carried frame stands.
    gp_XYZ height;   // the direction the outline's second coordinate runs along at the current station
    const auto orient = [&](const TubePiece& piece, double t, const Station& s) {
        height = s.tangent.Crossed(normal);
        if (!piece.oriented) return;
        const gp_XYZ along = piece.axialStart * (1.0 - t) + piece.axialEnd * t;
        const gp_XYZ axial = along - s.tangent * along.Dot(s.tangent);
        if (axial.Modulus() < 1e-6 || along.Modulus() < 1e-12) return;
        gp_XYZ width = s.tangent.Crossed(axial);
        width.Normalize();
        if (settled && width.Dot(normal) < 0.0) width = width * -1.0;   // the section is symmetric: no half turn
        normal = width;
        settled = true;
        height = s.tangent.Crossed(normal);
        // Upright: the height stays on the axial direction itself, on the side the square one points to.
        if (piece.upright) height = along / along.Modulus() * (along.Dot(height) < 0.0 ? -1.0 : 1.0);
    };
    const auto ringAt = [&](const Station& s) {
        std::vector<gp_XYZ> ring(around);
        for (int k = 0; k < around; ++k) ring[k] = s.point + normal * outline[k].X() + height * outline[k].Y();
        return ring;
    };

    BRep_Builder builder;
    TopoDS_Compound tube;
    builder.MakeCompound(tube);
    std::vector<gp_XYZ> junctionRing;   // the ring the previous piece ended on
    for (size_t i = 0; i < path.size(); ++i) {
        const Primitive& pr = path[i].prim;
        const int chords = chordsOf(pr, angularStep);
        std::vector<std::vector<gp_XYZ>> rings;
        rings.reserve(chords + 1);
        for (int c = 0; c <= chords; ++c) {
            const double t = static_cast<double>(c) / chords;
            const Station s = stationOf(pr, t);
            normal = carried(normal, previousTangent, s.tangent);
            normal = normal - s.tangent * normal.Dot(s.tangent);
            normal.Normalize();
            previousTangent = s.tangent;
            orient(path[i], t, s);
            // The piece starts on the ring its predecessor ended on: one mitre section, shared.
            if (c == 0 && i > 0) rings.push_back(junctionRing);
            else rings.push_back(ringAt(s));
        }
        // Where the next piece leaves in another direction, this one ends on the plane that bisects the corner: every
        // generator of the tube is run on to that plane, and the next piece starts from the same points.
        if (i + 1 < path.size()) {
            const Station end = stationOf(pr, 1.0);
            const Station next = stationOf(path[i + 1].prim, 0.0);
            gp_XYZ bisector = end.tangent + next.tangent;
            const double size = bisector.Modulus();
            if (size < 0.4)
                throw std::runtime_error("displayTube: the corner between '" + pr.label + "' and '" + path[i + 1].prim.label +
                                         "' turns back by more than 150 degrees; no mitre section exists");
            bisector /= size;
            const double slope = bisector.Dot(end.tangent);
            for (auto& q : rings.back()) q = q + end.tangent * (bisector.Dot(end.point - q) / slope);
            junctionRing = rings.back();
        }
        const bool capStart = i == 0, capEnd = i + 1 == path.size();
        const int nodes = static_cast<int>(rings.size()) * around + (capStart ? 1 : 0) + (capEnd ? 1 : 0);
        const int triangles = 2 * around * chords + (capStart ? around : 0) + (capEnd ? around : 0);
        Handle(Poly_Triangulation) mesh = new Poly_Triangulation(nodes, triangles, Standard_False);
        int node = 0;
        for (const auto& ring : rings)
            for (const auto& q : ring) mesh->SetNode(++node, gp_Pnt(q));
        const auto at = [&](int ringIndex, int k) { return ringIndex * around + (k % around) + 1; };
        int triangle = 0;
        for (int c = 0; c < chords; ++c)
            for (int k = 0; k < around; ++k) {
                mesh->SetTriangle(++triangle, Poly_Triangle(at(c, k), at(c, k + 1), at(c + 1, k + 1)));
                mesh->SetTriangle(++triangle, Poly_Triangle(at(c, k), at(c + 1, k + 1), at(c + 1, k)));
            }
        if (capStart) {
            mesh->SetNode(++node, gp_Pnt(stationOf(pr, 0.0).point));
            for (int k = 0; k < around; ++k) mesh->SetTriangle(++triangle, Poly_Triangle(node, at(0, k + 1), at(0, k)));
        }
        if (capEnd) {
            mesh->SetNode(++node, gp_Pnt(stationOf(pr, 1.0).point));
            for (int k = 0; k < around; ++k) mesh->SetTriangle(++triangle, Poly_Triangle(node, at(chords, k), at(chords, k + 1)));
        }
        // The chord sag of a tube of this size at this step: what a tessellator would call its deflection.
        mesh->Deflection(reach * (1.0 - std::cos(0.5 * angularStep)));
        TopoDS_Face face;
        builder.MakeFace(face, mesh);
        builder.Add(tube, face);
    }
    return tube;
}

} // namespace mvb
