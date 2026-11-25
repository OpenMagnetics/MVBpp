#pragma once
// A conductor's display drawing triangulated straight from its path, with no CAD solids (ABT #2489).
// A display drawing is looked at, never meshed or exported as a part. The path already says everything a picture
// needs: a centreline of SEG / ARC3 / SPIRAL / BLEND pieces and the wire's section.
#include "mvb/WireAssembler.h"

#include <TopoDS_Shape.hxx>
#include <gp_XY.hxx>
#include <gp_XYZ.hxx>

#include <vector>

namespace mvb {

// One piece of a conductor's path as the display tube draws it.
struct TubePiece {
    Primitive prim;
    // A section that is not a circle has ONE correct orientation: its height follows `axialStart` at the start of the
    // piece and `axialEnd` at its end (the direction is made square to the path there), as rectProfileWire places it.
    // Not oriented: the section is carried along the path without twist, which is all a circle needs.
    bool oriented = false;
    gp_XYZ axialStart{0, 1, 0}, axialEnd{0, 1, 0};
    // How an oriented section stands on a path that climbs along its axial direction. Square (the default): in the
    // plane square to the path, as a piece built on its own is (rectProfileWire). Upright: its height stays ON the
    // axial direction and its width square to both, as a conductor swept whole with a fixed binormal has it; the
    // section is then sheared by the climb, and cut through the axis it shows the rectangle the layout placed.
    bool upright = false;
};

// The outline of a round wire of `radius` [m]: a regular polygon whose sides span `angularStep` [rad].
std::vector<gp_XY> roundTubeOutline(double radius, double angularStep);
// The outline of a `width` x `height` [m] wire with corners of `cornerRadius` [m] (0: sharp), width first: the same
// section rectProfileWire draws. Each corner arc is drawn with sides that span `angularStep`.
std::vector<gp_XY> rectTubeOutline(double width, double height, double cornerRadius, double angularStep);

// The tube of section `outline` that follows `path`, as a compound of triangulated faces, one per path piece, in order.
// The faces carry triangles only (no surface): a picture, not a part. Consecutive pieces share their junction ring;
// where the path turns a corner that ring lies on the corner's bisecting plane. The two free ends are capped.
// `outline` is a convex polygon given counter-clockwise, in metres from the centreline. `angularStep` [rad] is the
// turning angle one chord of the path spans.
TopoDS_Shape displayTube(const std::vector<TubePiece>& path, const std::vector<gp_XY>& outline, double angularStep);

} // namespace mvb
