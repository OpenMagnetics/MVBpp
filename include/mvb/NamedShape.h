#pragma once

#include <TopoDS_Shape.hxx>
#include <string>
#include <utility>
#include <vector>

namespace mvb {

// What a drawn solid IS, independent of how it happens to be named (ABT #1169, WP0/D2).
//
// Names remain the only channel that survives a STEP round-trip, so OMFEM still classifies
// imported geometry by name. This enum is the SECOND, in-process channel: a consumer that
// holds the NamedShape vector (the Python/WASM bindings, the mesher when it calls MVB++
// directly, the geometry gates) reads the role instead of re-parsing a string. It is never
// written to STEP.
//
// Every producer sets it explicitly. The default is Core only because a default-constructed
// NamedShape has no shape either; a producer that forgets is a bug, and the role tests over
// the complete MAS fixtures are what catch it.
enum class Role {
    Core,         // a core piece
    CoreCoating,  // conformal coating shell around a core piece
    Turn,         // conductor: per-turn loop, or one continuous real-winding parallel
    TurnCoating,  // the insulated outer footprint of a conductor (enamel/serving/tape)
    Insulation,   // an insulation layer between sections/layers
    Bobbin,       // the former (tube + flanges), or a toroid base
    Pin,          // a bobbin/base terminal pin
    Divider,      // a chamber wall inside a multi-chamber bobbin
    Base,         // reserved: a standoff base that is not the former itself
    Spacer,       // a core spacer (additive gap)
    Shunt,        // a magnetic shunt
    Sleeve,       // insulating sleeve over a lead run
    Terminal,     // a planar port/cap face (FEM current injection), not a solid
    Solder,       // solder body joining a foil terminal
    FR4,          // planar-coil substrate board
};

// Stable, lower-case spelling of a Role — diagnostics and test failure messages only.
// Not a serialisation format and not what OMFEM classifies on.
inline const char* role_name(Role r) {
    switch (r) {
        case Role::Core:        return "core";
        case Role::CoreCoating: return "core_coating";
        case Role::Turn:        return "turn";
        case Role::TurnCoating: return "turn_coating";
        case Role::Insulation:  return "insulation";
        case Role::Bobbin:      return "bobbin";
        case Role::Pin:         return "pin";
        case Role::Divider:     return "divider";
        case Role::Base:        return "base";
        case Role::Spacer:      return "spacer";
        case Role::Shunt:       return "shunt";
        case Role::Sleeve:      return "sleeve";
        case Role::Terminal:    return "terminal";
        case Role::Solder:      return "solder";
        case Role::FR4:         return "fr4";
    }
    return "unknown";
}

// A geometric shape paired with a human-readable identifier. Names are
// produced at build time from MAS metadata (`MagneticCore::get_name`,
// `Turn::get_name`, bobbin name, etc.) with fallbacks, then carried
// through every downstream operation — boolean cuts, STEP export, mesh
// tagging — so the original logical identity survives the pipeline.
struct NamedShape {
    TopoDS_Shape shape;
    std::string  name;
    // ABT #685: names for the INDIVIDUAL SOLIDS of a multi-solid shape, in the order
    // TopExp_Explorer visits them. A round-wire conductor is a compound of per-primitive solids
    // (the conformal mitre assembly), and with only the compound named a viewer falls back to
    // numbering them itself — "Primary parallel 001", "…012" — a convention that lives in the
    // viewer, not in the file, so a solid cannot be named unambiguously in conversation. Filling
    // this makes exportSTEP write the compound as an assembly with a real name on every part.
    // Empty (the default) keeps the previous single-product behaviour exactly.
    std::vector<std::string> partNames;
    // ABT #1169: what this solid is. Set by every producer; carried through cuts, symmetry and
    // sectioning. NOT exported to STEP (names are the only channel that survives that).
    Role role = Role::Core;

    NamedShape() = default;
    NamedShape(TopoDS_Shape s, std::string n)
        : shape(std::move(s)), name(std::move(n)) {}
    NamedShape(TopoDS_Shape s, std::string n, Role r)
        : shape(std::move(s)), name(std::move(n)), role(r) {}
    NamedShape(TopoDS_Shape s, std::string n, std::vector<std::string> parts)
        : shape(std::move(s)), name(std::move(n)), partNames(std::move(parts)) {}
    NamedShape(TopoDS_Shape s, std::string n, std::vector<std::string> parts, Role r)
        : shape(std::move(s)), name(std::move(n)), partNames(std::move(parts)), role(r) {}
};

} // namespace mvb
