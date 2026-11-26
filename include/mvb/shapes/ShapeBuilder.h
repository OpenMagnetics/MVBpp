#pragma once

#include "MAS.hpp"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <BRepAlgoAPI_BooleanOperation.hxx>
#include <stdexcept>
#include <map>
#include <optional>
#include <string>
#include <memory>

namespace mvb {
namespace shapes {

class ShapeBuilder {
public:
    virtual ~ShapeBuilder() = default;

    // Polygon segment count used by derived builders when approximating
    // circles/cylinders. Set by the factory so subclasses can consult it
    // via member access without changing their virtual signatures.
    void setCorePolygonSegments(int segments) { m_corePolygonSegments = segments; }
    int  corePolygonSegments() const { return m_corePolygonSegments; }

    // Build a single core piece (includes intrinsic family rotation)
    virtual TopoDS_Shape buildPiece(const MAS::CoreShape& shapeData) const;

    // Apply machining (gap) to an already-rotated piece
    virtual TopoDS_Shape applyMachining(const TopoDS_Shape& piece,
                                        const MAS::Machining& machining,
                                        const std::map<std::string, double>& dims) const;

protected:
    // Subclasses implement the 2D profile in the XY plane
    virtual TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const = 0;

    // Subclasses implement the winding window cutout (default: none)
    virtual TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const;

    // Optional post-processing before intrinsic rotation (default: identity)
    virtual TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                                     const TopoDS_Shape& piece) const;

    // True for toroidal families (no -90° X rotation)
    virtual bool isToroidal() const { return false; }

    // Helper: extrude a face along Z by height
    static TopoDS_Shape extrude(const TopoDS_Face& face, double height);

    // Helper: create a rectangular box centered at origin
    static TopoDS_Shape makeBox(double x, double y, double z);

    // Helper: create a cylinder along Z axis, centered at origin
    static TopoDS_Shape makeCylinder(double height, double radius, int segments);

    // No dimension is ever read as zero because the record lacks it (ABT #1126, #1143): a piece
    // built from a missing letter is a different part, drawn and meshed as if it were the real one.
    //
    // The shape's dimensions as plain numbers; throws when the record carries none. Also names the
    // shape for the messages of the readers below (one builder builds one shape at a time).
    std::map<std::string, double> dimensionsOf(const MAS::CoreShape& shapeData) const;
    // A dimension the piece cannot be built without: throws when it is absent or not positive.
    double required(const std::map<std::string, double>& dims, const std::string& key) const;
    // A dimension that must be given but may be zero or negative (an offset).
    double requiredSigned(const std::map<std::string, double>& dims, const std::string& key) const;
    // A dimension that only the records with that feature carry (a centre hole, a slot). Absent or
    // zero is "no such feature", which is how MKF reads the same records; a negative one throws.
    std::optional<double> featureDimension(const std::map<std::string, double>& dims, const std::string& key) const;
    // The result of a boolean, or a throw naming the step: a cut or a union that did not complete
    // is never skipped, because what remains is not the shape.
    TopoDS_Shape done(BRepAlgoAPI_BooleanOperation& operation, const char* step) const;
    TopoDS_Shape done(BRepAlgoAPI_BooleanOperation&& operation, const char* step) const { return done(operation, step); }

    int m_corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS;
    mutable std::string m_shapeName = "<unnamed shape>";
};

// Factory — `corePolygonSegments` is applied to the returned builder.
std::unique_ptr<ShapeBuilder> createShapeBuilder(MAS::CoreShapeFamily family,
                                                   const std::string& subtype = "",
                                                   int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS);

// Return all shape family strings that have a registered builder.
std::vector<std::string> getSupportedShapeFamilies();

} // namespace shapes
} // namespace mvb
