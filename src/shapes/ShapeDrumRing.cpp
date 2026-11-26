#include "mvb/shapes/ShapeDrumRing.h"
#include "mvb/Utils.h"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>

namespace mvb {
namespace shapes {

// Build frame: identical to ShapeDrum (extruded along Z, then centred by the base class), so
// the ring is added in the already-centred frame — z from -L/2 to +L/2 keeps it concentric
// and vertically centred on the drum, which is how the matched pairs are assembled.
TopoDS_Shape ShapeDrumRing::applyExtras(const std::map<std::string, double>& dims,
                                        const TopoDS_Shape& piece) const {
    // Drum first: asymmetric-flange trim, bore, and the centring translation.
    TopoDS_Shape result = ShapeDrum::applyExtras(dims, piece);

    // The ring is the shield: a record without it is not this part, and a bare drum drawn (and
    // meshed) in its place would be an unshielded inductor passed off as a shielded one.
    const double j = required(dims, "J");
    const double k = required(dims, "K");
    const double l = required(dims, "L");
    if (!(k < j)) {
        throw std::runtime_error("shape '" + m_shapeName + "' has a ring bore K = " + std::to_string(k) +
                                 " that is not inside its outer diameter J = " + std::to_string(j));
    }

    TopoDS_Shape ringOuter = build_polygon_cylinder(l, j / 2.0, m_corePolygonSegments);
    ringOuter = translate_shape(ringOuter, 0.0, 0.0, -l / 2.0);
    // The ring bore faces the drum winding: circumscribe so faceting never bites in.
    TopoDS_Shape ringBore =
        build_polygon_cylinder(l, k / 2.0, m_corePolygonSegments, /*circumscribed=*/true);
    ringBore = translate_shape(ringBore, 0.0, 0.0, -l / 2.0);

    TopoDS_Shape ring = done(BRepAlgoAPI_Cut(ringOuter, ringBore), "boring the shield ring");

    // Disjoint solids by construction (the clearance gap): the fuse yields the two-solid
    // assembly, which is what a shielded drum physically is.
    return done(BRepAlgoAPI_Fuse(result, ring), "assembling the drum and its shield ring");
}

} // namespace shapes
} // namespace mvb
