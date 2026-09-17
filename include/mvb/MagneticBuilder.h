#pragma once
#include <limits>

#include "MAS.hpp"
#include "Utils.h"
#include "mvb/ConductorBuilder.h"
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
    double      scale                 = 1.0;     // extra USER scale; exporters emit mm natively (ABT #317) — keep 1.0
    int         symmetryPlanes        = 0;       // 0=full, 1=half, 2=quarter
    int         wirePolygonSegments    = DEFAULT_WIRE_POLYGON_SEGMENTS;
    int         corePolygonSegments    = DEFAULT_CORE_POLYGON_SEGMENTS;
    // true  → turns drawn at the OUTER (insulation) diameter (the web viewer).
    // false → turns drawn at the CONDUCTING (copper) diameter — DEFAULT (ABT #1261),
    //         and what FEM winding-loss meshing requires (LITZ → bare bundle as a solid).
    bool        paintCoating           = false;
    // true → real winding: ONE continuous conductor per (winding, parallel) instead of
    //        independent per-turn loops; MKF enriches with real-winding turn blocking on.
    bool        useRealWindingGeometry = false;
    // Real-winding only. false (OM / drawing) → fast per-run compound conductors for the 3D viewer.
    // true (FEM export) → the slow one-piece/conformal machinery: single body per parallel where a
    // sweep closes, a mitre-jointed conformal compound for dense toroids. Meshable but far slower.
    bool        femReady               = false;
    // Minimum radius for ANY drawn corner/fillet, metres; 0 = each site's policy radius.
    // Floors, never tightens (see ConductorBuilder::Options::minBendRadius).
    double      minBendRadius          = 0.0;
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
    //   scale                 : extra uniform USER scale on top of the exporters'
    //                           native millimetre output (ABT #317). Keep 1.0 for a
    //                           correct mm file; 1000 now double-scales (1e6x too big).
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
                                            bool paintCoating = false) const;
    std::vector<NamedShape> buildTurnsNamed(const OpenMagnetics::Coil& coil,
                                            const MAS::MagneticCore& core,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                            bool paintCoating = false) const;
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
    // useRealWindingGeometry: replace the per-turn closed loops with ONE continuous copper
    // body per (winding, parallel) ("<winding> parallel <k>", ConductorBuilder), enriching
    // through MKF with its real-winding turn blocking on. A magnetic that already carries a
    // geometricalDescription THROWS with the flag on — MKF must (re-)wind; caller-provided
    // turn positions are never silently re-interpreted.
    // Real-winding conductor CENTRELINES (sampled, collision-checked, seam-aimed), for
    // implicit/level-set winding meshing. No solids are built. See
    // ConductorBuilder::PathPolyline.
    // ABT #1248: a toroid's centrelines are returned in the EXPORTED frame (the frame
    // buildAllNamed's assembly and its STEP are in), so their terminal tips and directions are the
    // caps OMFEM sees. Other cores are unaffected.
    std::vector<ConductorBuilder::PathPolyline> buildRealWindingPaths(
        const OpenMagnetics::Magnetic& magnetic) const;

    // ---- ABT #1248: toroid mounting ------------------------------------------------------------
    // The mounting a toroid is built with: its bobbin's base.mounting when the bobbin is a toroid
    // base, otherwise MKF Settings toroid_mounting (default VERTICAL). Throws for a non-toroid.
    static ConductorBuilder::ToroidMounting toroidMountingOf(const OpenMagnetics::Magnetic& magnetic);
    static ConductorBuilder::ToroidMounting toroidMountingOf(const MAS::Magnetic& magnetic);
    // The full frame (mounting, build-frame lead direction, build -> exported rigid motion).
    // The magnetic must carry its wound turns (see ConductorBuilder::resolveToroidMountingFrame).
    static ConductorBuilder::ToroidMountingFrame toroidMountingFrameOf(
        const OpenMagnetics::Magnetic& magnetic);
    static ConductorBuilder::ToroidMountingFrame toroidMountingFrameOf(const MAS::Magnetic& magnetic);

    // ---- ABT #1215: terminal-lead copper length per winding ------------------------------
    // ABT #1248: a toroid's lead piece endpoints are in the exported frame, like the paths.
    // The lead copper the real-winding conductors carry beyond their turns, measured on the
    // finished centreline (see ConductorBuilder::measureTerminalLeadLengths for exactly which
    // primitives count and how they are split into entrance/exit). Pass the SAME settings the
    // conductors are built with -- buildAllNamed(..., paintCoating, ..., coreCoatingThickness,
    // useRealWindingGeometry=true, femReady) -- because the bend radii follow the drawn wire
    // radius and femReady decides which corners are filleted. With emitCoatingShells the
    // conductor that carries current is the bare one: pass paintCoating=false.
    // The magnetic must already be enriched through MKF's real-winding autocomplete.
    // Runs the conductor path planning once more (no solids), so it costs a fraction of a build.
    std::map<std::string, ConductorBuilder::TerminalLeadLength> measureTerminalLeadLengths(
        const OpenMagnetics::Magnetic& magnetic,
        bool paintCoating = false,
        bool femReady = true,
        int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
        int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS,
        double coreCoatingThickness = 0.0) const;
    // {"winding_<name>": {"terminal_lead_length_m": x, "parallels": n,
    //                     "ends": [{"parallel": k, "end": "entrance"|"exit", "length_m": l,
    //                               "pieces": [{"label", "kind", "length_m", "start_m", "end_m",
    //                                           "radius_m", "sweep_rad"}, ...]}, ...]}}
    static nlohmann::json terminalLeadLengthsToJson(
        const std::map<std::string, ConductorBuilder::TerminalLeadLength>& leads);
    // The sidecar sits next to the STEP with the extension swapped: out/design.step ->
    // out/design.leads.json. Returns the path written; throws when it cannot be written.
    static std::string terminalLeadSidecarPath(const std::string& stepPath);
    static std::string writeTerminalLeadSidecar(
        const std::map<std::string, ConductorBuilder::TerminalLeadLength>& leads,
        const std::string& stepPath);

    // The core-coating thickness to DRAW for this core, in metres, or 0 for "no coating layer".
    //
    // Only a coating the design DECLARES is drawn. MKF's Core::get_coating_thickness() also
    // invents a default jacket for any uncoated toroid -- right for the winding-to-core
    // dielectric path it exists for, wrong to put on screen for a design whose author never
    // mentioned a coating. Presence is decided from functionalDescription.coating; the
    // thickness itself is still resolved by MKF, which is what turns the name-only form
    // ("epoxy") into a datasheet number.
    //
    // MAGNETIC_EPOXY returns 0: that is the powder-loaded shield cap moulded over the WINDING
    // of a semishielded drum, drawn separately (and translucently) by drawCoreShell, not an
    // insulating jacket on the ferrite surface.
    //
    // Public because every entry point that assembles a magnetic has to make the same call --
    // drawMagnetic here, and the WASM/Python bindings, which reach buildAllNamed directly and
    // are what every 3D viewer actually goes through.
    static double declaredCoreCoatingThickness(const MAS::MagneticCore& core);

    // The coating layer itself: the core surface offset outward by `thickness`, minus the
    // core, i.e. a uniform-thickness conformal wrap (outer, inner, top and bottom). Returns a
    // null shape when the offset or the cut fails -- a sharp-cornered core OCCT cannot offset
    // -- and the caller then draws nothing rather than fabricating a layer.
    //
    // Public for the same reason: the 3D viewers do not call drawMagnetic, they compose the
    // assembly from drawCore / drawCoreShell / drawBobbin / drawTurns, so the coating has to
    // be reachable as its own product.
    static TopoDS_Shape buildCoreCoatingShell(const TopoDS_Shape& core, double thickness);

    std::vector<NamedShape> buildAllNamed(const MAS::Magnetic& magnetic,
                                          bool includeBobbin = true,
                                          int symmetryPlanes = 0,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                          int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS,
                                          bool paintCoating = false,
                                          bool emitCoatingShells = false,
                                          bool includeInsulation = false,
                                          double coreCoatingThickness = 0.0,
                                          bool useRealWindingGeometry = false,
                                          bool femReady = false) const;
    std::vector<NamedShape> buildAllNamed(const OpenMagnetics::Magnetic& magnetic,
                                          bool includeBobbin = true,
                                          int symmetryPlanes = 0,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                          int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS,
                                          bool paintCoating = false,
                                          bool emitCoatingShells = false,
                                          bool includeInsulation = false,
                                          // >0: build the core's insulating coating as a conformal
                                          // shell solid ("<core> coating") of this thickness [m].
                                          double coreCoatingThickness = 0.0,
                                          bool useRealWindingGeometry = false,
                                          bool femReady = false) const;

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

    // Turns ONLY, in real-winding form: one continuous copper body per (winding, parallel),
    // with the real leads, pitch and dragbacks. The turns-only counterpart of
    // buildAllNamed(..., useRealWindingGeometry=true), for a consumer that draws core,
    // bobbin and turns as SEPARATE meshes (the web 3D viewer, which colours and toggles
    // them independently) and so cannot use the single-call assembly.
    //
    // The magnetic must already be enriched through MKF's real-winding autocomplete —
    // same precondition as buildAllNamed with the flag on. The core is built internally
    // (not returned) because the conductor builder aims the terminal leads at the real
    // window opening; see buildRealWindingConductorsNamed.
    //
    // diagnosticSkipCollisionCheck builds the conductors even when two of them overlap, so
    // the overlap can be LOOKED AT instead of only read about in an exception. The result is
    // interpenetrating copper: a picture, never a part. Tools and tests investigating a
    // specific refusal only — nothing in the library sets it.
    std::vector<NamedShape> buildRealWindingTurnsNamed(
        const OpenMagnetics::Magnetic& magnetic,
        int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
        int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS,
        bool paintCoating = false,
        bool femReady = false,
        bool diagnosticSkipCollisionCheck = false) const;

    // ---- ABT #1169 (WP0): the accessory-solid attachment point ------------------------
    //
    // Everything a real magnetic carries that is neither a core piece nor a conductor:
    // core spacers (WP1), bobbin dividers and pins (WP2/WP4), magnetic shunts (WP6), lead
    // sleeves (WP7), a toroid base (WP4). Those work packages add their builders HERE and
    // nowhere else, so the whole assembly, the STEP export, the 2D sections and the OMFEM
    // meshers all pick them up from one place.
    //
    // Called once from each buildAllNamed overload, immediately after the bobbin block.
    // Names must follow D2 of the manufacturing-fidelity proposal ("Spacer_<i>",
    // "<bobbin> pin <name>", "<bobbin> divider <i>", "Shunt_<i>",
    // "<winding> parallel <p> <entrance|exit> sleeve", "<bobbin> base") because names are
    // the only channel that survives the STEP round-trip into OMFEM's classify(); each
    // solid must also carry its Role. A solid that overlaps the former (a divider, a pin)
    // has to be added to the bobbin cutter list in buildAllNamed as well — the bobbin is
    // already cut by the time this runs.
    //
    // Empty in WP0 on purpose: the hook lands with the strict classifier so that the first
    // accessory to be drawn cannot be silently meshed as core.
    struct AccessoryOptions {
        bool includeBobbin          = false;
        int  wirePolygonSegments    = DEFAULT_WIRE_POLYGON_SEGMENTS;
        int  corePolygonSegments    = DEFAULT_CORE_POLYGON_SEGMENTS;
        bool paintCoating           = false;
        bool useRealWindingGeometry = false;
        bool femReady               = false;
        // ABT #1173 (WP4): the real-winding toroid terminal plane of the conductors this assembly
        // drew (ConductorBuilder::Options::toroidTerminalPlaneOut); NaN when there is none. A
        // toroid base is drawn with its top face on it (BaseBuilder.h).
        double toroidTerminalPlaneY = std::numeric_limits<double>::quiet_NaN();
    };
    void appendAccessorySolids(std::vector<NamedShape>& all,
                               const MAS::Magnetic& magnetic,
                               const AccessoryOptions& opts) const;
    void appendAccessorySolids(std::vector<NamedShape>& all,
                               const OpenMagnetics::Magnetic& magnetic,
                               const AccessoryOptions& opts) const;

  private:
    // Shared by buildRealWindingConductorsNamed and measureTerminalLeadLengths.
    ConductorBuilder::Options realWindingConductorOptions(
        const OpenMagnetics::Magnetic& magnetic, const std::vector<NamedShape>& coreShapes,
        int wirePolygonSegments, bool femReady, MAS::CoreBobbinProcessedDescription& bobbinPd,
        bool& toroidalCore) const;
    // Single implementation of real-winding conductor emission, shared by buildAllNamed
    // and buildRealWindingTurnsNamed so the assembly and the viewer cannot drift apart.
    // coreShapes are handed to ConductorBuilder as lead-aiming obstacles (ignored for
    // toroids); they are not part of the result.
    std::vector<NamedShape> buildRealWindingConductorsNamed(
        const OpenMagnetics::Magnetic& magnetic,
        const std::vector<NamedShape>& coreShapes,
        int wirePolygonSegments,
        bool paintCoating,
        bool emitCoatingShells,
        bool femReady,
        bool diagnosticSkipCollisionCheck = false,
        // ABT #1173: receives the LOWEST toroid terminal plane over the conductor builds emitted
        // (bare and coating shells lay their leads at their own radii); NaN when none.
        double* toroidTerminalPlaneOut = nullptr) const;
};

} // namespace mvb
