#pragma once

#include "MAS.hpp"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <numbers>

// Forward declaration for MKF Magnetic
namespace OpenMagnetics { class Magnetic; }

namespace mvb {

constexpr int DEFAULT_CORE_POLYGON_SEGMENTS = 16;
// WIRE CROSS-SECTIONS ARE EXACT BY DEFAULT (Alf, 2026-08-23). 0 = the true circle, which sweeps
// to ANALYTIC surfaces (cylinder/torus); N > 0 approximates it by an N-gon prism, which sweeps to
// N BSpline faces each carrying a pole per spine sample.
//
// The 16-gon was costing two orders of magnitude of memory for accuracy it did not buy. Measured
// 2026-08-23 on the real-winding path, same design, only this value changed:
//     01_simple_inductor_etd34_n87   23 turns    2,809 MB  ->    200 MB   (14x)
//     14_dab_xfmr_pm8770_n97         96 turns   27,091 MB  ->     92 MB  (300x)
// -- 120 to 280 MB of RAM PER TURN of 0.8 mm round wire, which is what took this box to 0 GB
// available and, at a 16 GB ulimit, turns an allocation failure into a SEGFAULT inside OCC (it
// certifies the design, then dies). The two discretisations were also unbalanced: the spine is
// sampled to a sag of 2% of the wire radius (113-165 samples per revolution here) while the
// 16-gon section is itself 1.9% off the true circle, so the fine sampling was holding an accuracy
// the profile had already spent. See ABT #860.
//
// TRADE-OFF, deliberately taken: an exact circle revolves into a PERIODIC surface (the seam of a
// cylinder/torus), and the faceted paths exist partly because gmsh refuses those ("Impossible to
// mesh periodic surface" -- see MVB_SPLIT_PROFILE and DEFAULT_WIRE_REVOLUTION_SEGMENTS below).
// Any consumer that needs node-conformal meshes should pass its own segment count rather than
// rely on this default; the web already passes one explicitly.
constexpr int DEFAULT_WIRE_POLYGON_SEGMENTS = 0;
// Azimuthal segmentation of revolved wire turns. Even with a polygonal
// cross-section, MakeRevol produces analytic surfaces of revolution that
// STEP serialises as CYLINDRICAL/TOROIDAL_SURFACE. Setting this > 0 replaces
// the revolve with a ThruSections loft over N angular slices → fully planar.
// Set to 0 for the legacy MakeRevol behaviour.
constexpr int DEFAULT_WIRE_REVOLUTION_SEGMENTS = 12;

// Extract nominal double from MAS Dimension variant
double flatten_dimension(const MAS::Dimension& dim);

// Extract all nominal dimensions from a shape dimension map
std::map<std::string, double> flatten_dimensions(const std::map<std::string, MAS::Dimension>& dims);

// The value of a dimension the profile cannot be built without.
//
// A dimension that silently defaults to zero does not fail: it extrudes a
// degenerate solid that reads as a modelling bug rather than the data gap it
// is. RM 7LP shipped for years with no C and rendered as a mangled polyhedron
// (ABT #1126). Throws naming the shape and the dimension.
double require_dimension(const std::map<std::string, double>& dims,
                         const std::string& key,
                         const std::string& shapeName);

// Build a polygon-approximated circle wire in the XY plane centered at origin
// segments = 0 yields a perfect BRep circle edge.
//
// Faceting direction rule: a faceted solid must be a SUBSET of its nominal solid so
// faceting can never create a collision the true geometry does not have. Boundaries
// that ADD material (outer walls, posts, columns) use the default INSCRIBED polygon
// (vertices on the circle, chords inside). Boundaries that CARVE a void another part
// lives in (a toroid bore, a winding-window tool, a gap-cutting tool) must pass
// circumscribed = true: the polygon's apothem equals the radius (vertices at
// r / cos(pi/N)), so the carved void contains the nominal void.
TopoDS_Wire build_polygon_circle(double radius, int segments, bool circumscribed = false);

// Build a polygon-approximated cylinder solid along Z axis
// segments = 0 yields exact revolved circle; circumscribed as build_polygon_circle
TopoDS_Shape build_polygon_cylinder(double height, double radius, int segments,
                                    bool circumscribed = false);

// Build a ring (torus approximation) as a solid by lofting a polygonal
// circular cross-section (cross_segments sides) at `revolution_segments`
// azimuthal stations around the Y axis. Produces only PLANAR faces — STEP
// export contains no CYLINDRICAL / TOROIDAL / SURFACE_OF_REVOLUTION.
// turn_radius = major radius, wire_radius = minor, y = ring plane height.
// revolution_segments <= 0 falls back to a classical MakeRevol (analytic).
TopoDS_Shape build_polygon_ring(double turn_radius, double wire_radius,
                                 double y, int cross_segments,
                                 int revolution_segments);

// Apply a 3D rotation (radians) to a shape about X, Y, Z axes in order
TopoDS_Shape rotate_shape(const TopoDS_Shape& shape, double rx, double ry, double rz);

// Translate a shape
gp_Trsf translation_trsf(double x, double y, double z);
TopoDS_Shape translate_shape(const TopoDS_Shape& shape, double x, double y, double z);

// Get family string from enum
std::string core_shape_family_to_string(MAS::CoreShapeFamily family);

// Enumerate every core shape family that has a builder registered in the
// factory (mirrors MVB.js getSupportedFamilies()).
std::vector<std::string> get_supported_families();

// Preprocess JSON to add missing "nominal" fields to dimension objects
// (MVB test data often has only min/max)
void patch_dimension_nominals(nlohmann::json& j);

// Safely enrich a raw MAS::Magnetic using MKF's magnetic_autocomplete.
// This avoids object-slicing and Coil::wind() crashes for raw MAS files
// (e.g. with "Basic" bobbins) by constructing OpenMagnetics::Coil with
// windInConstructor=false before calling MKF enrichment.
//
// useRealWindingGeometry flips MKF's Settings::set_coil_use_real_winding_geometry for the
// duration of the enrichment (exception-safe RAII guard), so Coil::wind() places turns
// with connection-lead turn blocking — the positions ConductorBuilder then honours
// verbatim. Default false = byte-identical legacy behaviour.
OpenMagnetics::Magnetic magnetic_autocomplete_safe(const MAS::Magnetic& magnetic,
                                                   bool useRealWindingGeometry = false);
OpenMagnetics::Magnetic magnetic_autocomplete_safe(const nlohmann::json& magneticJson,
                                                   bool useRealWindingGeometry = false);

// Helpers for cutting bobbin with cores/turns to match Python MVB behavior
bool is_shape_usable(const TopoDS_Shape& shape);
TopoDS_Shape cut_bobbin(const TopoDS_Shape& bobbin, const std::vector<TopoDS_Shape>& cutters);

} // namespace mvb
