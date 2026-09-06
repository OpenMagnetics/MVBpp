#pragma once

#include "mvb/NamedShape.h"

#include <TopoDS_Shape.hxx>
#include <string>
#include <vector>

namespace mvb {

// Export a vector of named shapes as a STEP file using OCCT's XCAF layer
// so each shape's `name` becomes a TDataStd_Name label readable by any
// XCAF-aware reader (FreeCAD, OCCT's own STEPCAFControl_Reader, gmsh's
// OCAF path). Replaces the legacy single-"Assembly"-label behaviour.
bool exportSTEP(const std::vector<NamedShape>& shapes,
                const std::string& filepath);

struct StepExportOptions {
    // FEM product only (ABT #490 class): re-express every solid that carries a PERIODIC surface
    // (cylinder, cone, torus, surface of revolution) as B-splines, so gmsh does not route its
    // faces to the periodic mesher. Done HERE, on the millimetre geometry the file is written in:
    // a B-spline converted in metres and scaled x1000 keeps metre-scale knots under mm-scale
    // poles, and BOPAlgo then finds sporadic self-intersections on faceted revolves (ABT #1111).
    // The conversion is exact (a conic is a rational B-spline); it moves no point.
    bool nurbsPeriodicSolids = false;
};

bool exportSTEP(const std::vector<NamedShape>& shapes,
                const std::string& filepath,
                const StepExportOptions& options);

// Legacy overload — pairs shapes and names positionally, delegates to the
// NamedShape overload. Kept so existing callers compile; new code should
// use NamedShape directly.
bool exportSTEP(const std::vector<TopoDS_Shape>& shapes,
                const std::vector<std::string>& names,
                const std::string& filepath);

// Export a single compound shape as an STL file.
bool exportSTL(const TopoDS_Shape& compound,
               const std::string& filepath);

// Export a list of shapes as binary STL bytes (all shapes fused into one
// compound before tessellation). Returns an empty string if the input is
// empty or all shapes are null.
//
// toleranceMm       : linear deflection for meshing (mm, in the shape's
//                     current units — pass e.g. 0.1 for mm-scaled geometry)
// angularTolerance  : angular deflection for meshing (radians)
// binary            : if true, produces binary STL; otherwise ASCII.
std::string exportSTLToBytes(const std::vector<TopoDS_Shape>& shapes,
                             double toleranceMm      = 0.1,
                             double angularTolerance = 0.2,
                             bool binary             = true);

// Read a STEP file written with XCAF labels back into a vector of named
// shapes. Unnamed top-level shapes get an empty name string. Returns an
// empty vector on read failure.
std::vector<NamedShape> importSTEP(const std::string& filepath);

} // namespace mvb
