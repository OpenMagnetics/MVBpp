#include "mvb/shapes/ShapeE.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace mvb {
namespace shapes {

bool ShapeE::hasEfdPole(const std::map<std::string, double>& dims) const {
    if (featureDimension(dims, "q")) return true;
    const auto f2 = featureDimension(dims, "F2");
    return f2 && *f2 + 1e-9 < required(dims, "F");
}

TopoDS_Face ShapeE::buildProfile(const std::map<std::string, double>& dims) const {
    const double a = required(dims, "A") / 2.0;
    const double cFull = required(dims, "C");

    // EFD variant: a record whose pole is shallower than it is wide (F2 < F), or that carries the
    // chamfer q, describes the EFD pole, and then needs all of it: the pole's width F, its depth
    // F2, its chamfer q and its offset K (signed: three EFD sizes have it below the base). The
    // profile is asymmetric with a triangular dent notch at the top-center (and a mini-indent at
    // the bottom-center when K > 0). Matches MVB Python Efd.get_shape_base.
    if (hasEfdPole(dims)) {
        const double q = required(dims, "q");
        const double f = required(dims, "F");
        const double f2 = required(dims, "F2");
        const double k = requiredSigned(dims, "K");
        double topC = cFull - k - f2 / 2.0;
        double bottomC = k + f2 / 2.0;
        double dentHeight = cFull * 2.0 / 5.0;
        double dentTopW = f / 2.0;
        double dentBottomW = f / 2.0 - q;

        BRepBuilderAPI_MakePolygon poly;
        poly.Add(gp_Pnt(-a, topC, 0.0));
        poly.Add(gp_Pnt(-dentTopW, topC, 0.0));
        poly.Add(gp_Pnt(-dentBottomW, topC - dentHeight, 0.0));
        poly.Add(gp_Pnt(dentBottomW, topC - dentHeight, 0.0));
        poly.Add(gp_Pnt(dentTopW, topC, 0.0));
        poly.Add(gp_Pnt(a, topC, 0.0));
        poly.Add(gp_Pnt(a, -bottomC, 0.0));
        if (k > 0.0) {
            double miniW = f / 2.0 - q;
            if (miniW < 0.0) miniW = 0.0;
            poly.Add(gp_Pnt(miniW, -bottomC, 0.0));
            poly.Add(gp_Pnt(miniW, -bottomC + k, 0.0));
            poly.Add(gp_Pnt(-miniW, -bottomC + k, 0.0));
            poly.Add(gp_Pnt(-miniW, -bottomC, 0.0));
        }
        poly.Add(gp_Pnt(-a, -bottomC, 0.0));
        poly.Close();
        BRepBuilderAPI_MakeWire wire(poly.Wire());
        BRepBuilderAPI_MakeFace face(wire.Wire());
        return face.Face();
    }

    double c = cFull / 2.0;
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a, c, 0.0));
    poly.Add(gp_Pnt(a, c, 0.0));
    poly.Add(gp_Pnt(a, -c, 0.0));
    poly.Add(gp_Pnt(-a, -c, 0.0));
    poly.Close();

    BRepBuilderAPI_MakeWire wire(poly.Wire());
    BRepBuilderAPI_MakeFace face(wire.Wire());
    return face.Face();
}

TopoDS_Shape ShapeE::buildWindingWindow(const std::map<std::string, double>& dims) const {
    const double b = required(dims, "B");
    const double d = required(dims, "D");
    const double e = required(dims, "E");
    const double f = required(dims, "F");
    const double c = required(dims, "C");

    // Winding window box: E × 2C × D, Y-centered. The 2C height (vs C)
    // ensures the cut fully covers asymmetric profiles (EFD with K>0 has
    // top_c ≠ bottom_c). For symmetric profiles the extra Y outside the
    // profile is clipped away by the profile intersection, so it's safe.
    gp_Pnt wwCorner(-e / 2.0, -c, b - d);
    TopoDS_Shape wwBox = BRepPrimAPI_MakeBox(wwCorner, e, 2.0 * c, d).Shape();

    // F2 is the centre column's depth where it is not the core's own (EL: longer than F, a
    // stadium; EFD: shorter). A record without it has the plain F x C column.
    const double f2 = featureDimension(dims, "F2").value_or(f);
    bool oblong = std::abs(f2 - f) > 0.0001;

    TopoDS_Shape column;
    if (oblong && f2 > f) {
        // Stadium-shaped central column for EL cores
        double half_width = f / 2.0;      // X direction, semicircle radius
        double half_depth = f2 / 2.0;     // Y direction
        double rect_half_length = half_depth - half_width;
        double zCenter = b - d / 2.0;

        if (rect_half_length <= 0.0) {
            // Actually round
            column = BRepPrimAPI_MakeCylinder(half_width, d).Shape();
            // MakeCylinder is base-at-zero; center it at zCenter
            column = translate_shape(column, 0.0, 0.0, zCenter - d / 2.0);
        } else {
            gp_Pnt rectCorner(-half_width, -rect_half_length, b - d);
            TopoDS_Shape centerRect = BRepPrimAPI_MakeBox(rectCorner, f, rect_half_length * 2.0, d).Shape();

            TopoDS_Shape topCyl = BRepPrimAPI_MakeCylinder(half_width, d).Shape();
            topCyl = translate_shape(topCyl, 0.0, rect_half_length, zCenter - d / 2.0);

            TopoDS_Shape bottomCyl = BRepPrimAPI_MakeCylinder(half_width, d).Shape();
            bottomCyl = translate_shape(bottomCyl, 0.0, -rect_half_length, zCenter - d / 2.0);

            TopoDS_Shape fused = done(BRepAlgoAPI_Fuse(centerRect, topCyl), "joining the stadium column's first end");
            column = done(BRepAlgoAPI_Fuse(fused, bottomCyl), "joining the stadium column's second end");
        }
    } else if (oblong && f2 < f) {
        // EFD: leave the WW cut as a plain E × 2C × D box (no column carved
        // out). The column is added back separately in applyExtras so it
        // survives the profile dent intact at full F × F2 × B size.
        return wwBox;
    } else {
        // Rectangular central column: F x C x D, centered at same Z
        gp_Pnt ccCorner(-f / 2.0, -c / 2.0, b - d);
        column = BRepPrimAPI_MakeBox(ccCorner, f, c, d).Shape();
    }

    return done(BRepAlgoAPI_Cut(wwBox, column), "leaving the centre column in the winding window");
}

TopoDS_Shape ShapeE::applyExtras(const std::map<std::string, double>& dims,
                                 const TopoDS_Shape& piece) const {
    const double b = required(dims, "B");
    const double f = required(dims, "F");

    TopoDS_Shape result = piece;
    // EFD: add central column with chamfered (octagonal)
    // corners — F × F2 with q-sized 45° chamfers at each corner. Matches
    // MVB.js EFDShape.getShapeExtras.
    if (hasEfdPole(dims)) {
        const double f2 = required(dims, "F2");
        const double q = required(dims, "q");
        double hf = f / 2.0;
        double hf2 = f2 / 2.0;
        TopoDS_Shape column;
        {
            BRepBuilderAPI_MakePolygon poly;
            poly.Add(gp_Pnt( hf - q,  hf2,    0.0));
            poly.Add(gp_Pnt(-hf + q,  hf2,    0.0));
            poly.Add(gp_Pnt(-hf,      hf2 - q, 0.0));
            poly.Add(gp_Pnt(-hf,     -hf2 + q, 0.0));
            poly.Add(gp_Pnt(-hf + q, -hf2,    0.0));
            poly.Add(gp_Pnt( hf - q, -hf2,    0.0));
            poly.Add(gp_Pnt( hf,     -hf2 + q, 0.0));
            poly.Add(gp_Pnt( hf,      hf2 - q, 0.0));
            poly.Close();
            BRepBuilderAPI_MakeFace face(poly.Wire());
            column = BRepPrimAPI_MakePrism(face.Face(), gp_Vec(0, 0, b)).Shape();
        }
        result = done(BRepAlgoAPI_Fuse(result, column), "adding the EFD centre column");
    }
    return translate_shape(result, 0.0, 0.0, -b);
}

} // namespace shapes
} // namespace mvb
