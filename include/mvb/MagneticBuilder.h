#pragma once

#include "MAS.hpp"
#include "Utils.h"
#include "mvb/NamedShape.h"
#include "mvb/Symmetry.h"
#include <TopoDS_Shape.hxx>
#include <vector>
#include <string>
#include <map>

// Forward declarations for MKF-enriched magnetic overloads
namespace OpenMagnetics { class Magnetic; class Coil; }

namespace mvb {

// Unified configuration for drawMagnetic / drawMagneticToBytes.
// All bindings (Python, WASM) accept this struct so the API surface
// is identical across languages.
struct DrawConfig {
    std::string format                = "step";  // "step" or "stl"
    bool        includeBobbin         = true;
    double      scale                 = 1.0;     // 1000 for mm export
    int         symmetryPlanes        = 0;       // 0=full, 1=half, 2=quarter
    int         wirePolygonSegments    = DEFAULT_WIRE_POLYGON_SEGMENTS;
    int         corePolygonSegments    = DEFAULT_CORE_POLYGON_SEGMENTS;
    // true  → turns drawn at the OUTER (insulation) diameter (visualisation).
    // false → turns drawn at the CONDUCTING (copper) diameter — required for
    //         FEM winding-loss meshing (LITZ → bare bundle as a solid).
    bool        paintCoating           = true;
};

class MagneticBuilder {
public:
    // Build geometry and export to STEP ("step") or STL ("stl").
    //
    // Parameters
    // ----------
    //   outputPath            : directory where magnetic.step / magnetic.stl is written
    //   format                : "step" or "stl"
    //   includeBobbin         : include bobbin geometry
    //   scale                 : uniform scale factor (use 1000 for mm export)
    //   symmetryPlanes        : 0=full, 1=half, 2=quarter domain
    //   wirePolygonSegments   : <=0 = exact torus, >0 = faceted polygon (wire cross-section)
    //   corePolygonSegments   : polygon segments for core cylinders/circles
    std::string drawMagnetic(const MAS::Magnetic& magnetic,
                             const std::string& outputPath,
                             const std::string& format = "step",
                             bool includeBobbin = true,
                             double scale = 1.0,
                             int symmetryPlanes = 0,
                             int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                             int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Config-based overload — preferred for bindings.
    std::string drawMagnetic(const MAS::Magnetic& magnetic,
                             const std::string& outputPath,
                             const DrawConfig& cfg) const;

    // Overload that accepts an already-enriched OpenMagnetics::Magnetic
    // to avoid object-slicing issues with MAS::Magnetic
    std::string drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                             const std::string& outputPath,
                             const std::string& format = "step",
                             bool includeBobbin = true,
                             double scale = 1.0,
                             int symmetryPlanes = 0,
                             int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                             int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Config-based overload — preferred for bindings.
    std::string drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                             const std::string& outputPath,
                             const DrawConfig& cfg) const;

    // Named-shape overloads. Each returned element carries the logical
    // name (core name / "Turn_<i>" or Turn::get_name / bobbin name) so the
    // identity survives downstream operations (symmetry cut, STEP export
    // with XCAF labels, mesh tagging). These are the only public builders;
    // the prior unnamed variants (buildCore / buildBobbin / buildTurns) have
    // been removed in favour of this Named API.
    std::vector<NamedShape> buildCoreNamed(const MAS::MagneticCore& core,
                                            int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;
    std::vector<NamedShape> buildTurnsNamed(const MAS::Coil& coil,
                                            const MAS::MagneticCore& core,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                            bool paintCoating = true) const;
    std::vector<NamedShape> buildTurnsNamed(const OpenMagnetics::Coil& coil,
                                            const MAS::MagneticCore& core,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                            bool paintCoating = true) const;
    NamedShape buildBobbinNamed(const MAS::Coil& coil,
                                const MAS::MagneticCore& core,
                                int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;
    NamedShape buildBobbinNamed(const OpenMagnetics::Coil& coil,
                                const MAS::MagneticCore& core,
                                int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Assemble all geometry as named shapes (no export). Optionally applies
    // symmetry cuts according to symmetryPlanes (0=full, 1=half, 2=quarter).
    // emitCoatingShells: emit BOTH the bare-copper turn ("<turn>") and the outer insulated
    // footprint ("<turn> coating") per turn, so a thermal mesh can resolve the low-k wire enamel
    // (the mesher fragments the overlap into copper core + coating annulus). Implies copper turns.
    // includeInsulation: also emit INSULATION-layer solids ("insulation_layer_<i>") for the
    // inter-layer/inter-section tape, when they carry real thickness (zero-thickness placeholders
    // are skipped). For thermal FEA -- a low-k conduction barrier between windings.
    std::vector<NamedShape> buildAllNamed(const MAS::Magnetic& magnetic,
                                          bool includeBobbin = true,
                                          int symmetryPlanes = 0,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                          int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS,
                                          bool paintCoating = true,
                                          bool emitCoatingShells = false,
                                          bool includeInsulation = false) const;
    std::vector<NamedShape> buildAllNamed(const OpenMagnetics::Magnetic& magnetic,
                                          bool includeBobbin = true,
                                          int symmetryPlanes = 0,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                          int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS,
                                          bool paintCoating = true,
                                          bool emitCoatingShells = false,
                                          bool includeInsulation = false) const;

    // ---- Standalone builders for the unified bindings API -----------------
    //
    // Build a single core piece from a MAS::CoreShape. Validates the shape
    // through OpenMagnetics::CorePiece::factory (which also fills in derived
    // parameters), then dispatches to mvb's ShapeBuilder for geometry.
    // Returns a single NamedShape with name = shape.name (or family code).
    NamedShape buildCorePieceNamed(const MAS::CoreShape& shape,
                                    int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Build a bobbin from a fully-populated MAS::Bobbin. Throws if
    // `processed_description` is not present. `axisIsY=true` orients the
    // bobbin tube along Y (matches concentric-core convention); pass false
    // for toroidal bobbins (kept along Z).
    NamedShape buildBobbinNamedFromBobbin(const MAS::Bobbin& bobbin,
                                          bool axisIsY = true,
                                          int polygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Build a list of turns where every Turn carries its own dimensions and
    // cross_sectional_shape (no Wire/bobbin lookup). Throws if any turn is
    // missing required fields. Toroidal layout is auto-detected from the
    // presence of `additional_coordinates`; otherwise concentric round
    // column is assumed.
    // paintCoating must stay true: standalone turns carry only their outer
    // footprint, so the conducting cross-section cannot be recovered here and
    // paintCoating=false throws (use a full Magnetic JSON via drawMagnetic).
    std::vector<NamedShape> buildTurnsNamedFromTurns(
        const std::vector<MAS::Turn>& turns,
        int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
        bool paintCoating = true) const;
};

} // namespace mvb
