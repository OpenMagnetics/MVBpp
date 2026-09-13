#pragma once

#include "MAS.hpp"
#include "mvb/NamedShape.h"

#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace mvb {

// ABT #1173 / WP4: the base a toroid sits on (MAS bobbin family t with functionalDescription.base,
// seated by MKF's Bobbin::create_toroid_bobbin_on_base) and the pins it carries.
//
// WHERE THE BASE GOES. MKF places the base's pins for the bare coated ring: seating plane
// y_s = -(C/2 + standoff), pin centres half a pin below it (Bobbin::get_toroid_base_pin_rail_distance).
// A WOUND ring does not rest on its core: the copper under the ring and the entrance leads run
// below -C/2, and a plate whose top face sat at -C/2 would cut through them. The base's top face is
// therefore the real-winding terminal plane, where every terminal drop ends
// (ConductorBuilder: max(2 OD of the thinnest wire, 1 OD of the thickest) below the lowest copper,
// decision 7 of the manufacturing-features proposal). So the base is drawn with its top face on that
// plane and its bottom face `height` below it, and every pin keeps MKF's x and z and is moved along
// Y by the one offset that puts its top on the base's bottom face (the board). Nothing else of MKF's
// footprint is re-derived.
//
// FOOTPRINT. MKF separates a horizontal base's two pin rows along Z (rows at z = -+rowDistance/2,
// pins along X). The base's `length` therefore runs along Z and its `width` along X, centred on the
// ring axis; a pin outside that footprint throws (the base would not carry it).
//
// POCKET. Cut only when the record states pocketInnerDiameter AND pocketDepth (a cylinder down
// from the top face, axis Y). Without a pocket the base is a flat plate and its height must equal its
// standoff; with one, height - pocketDepth must equal the standoff. Anything else throws.
//
// VERTICAL MOUNTING is not drawn: with boatWidth absent the boat is not described at all, and with
// it present the board normal of an on-edge ring (-Z) belongs to the toroid terminal-orientation
// rework. Both throw.
//
// NOT HERE: the terminal drops are not re-aimed at the pins (moved to the terminal-orientation
// ticket); they keep landing on the terminal plane at their own rim azimuths.
class BaseBuilder {
public:
    // The base of a coil's bobbin, when it has one (family t with functionalDescription.base).
    // BobbinT is the bobbin alternative of the coil's bobbin variant (MAS::Bobbin or OpenMagnetics::Bobbin);
    // a bobbin given by name, or any other alternative, has no base.
    template <typename BobbinT, typename BobbinVariantT>
    static std::optional<MAS::BobbinBase> baseOf(const BobbinVariantT& bobbinVariant) {
        const BobbinT* bobbin = std::get_if<BobbinT>(&bobbinVariant);
        if (!bobbin) return std::nullopt;
        const auto functional = bobbin->get_functional_description();
        if (!functional || functional->get_family() != MAS::BobbinFamily::T) return std::nullopt;
        return functional->get_base();
    }

    // Throws for a base this builder cannot draw: vertical mounting.
    static void requireDrawableMounting(const MAS::BobbinBase& base, const std::string& bobbinName);

    // The base solid ("<bobbinName> base", Role::Base) and its pins ("<bobbinName> pin <n>",
    // Role::Pin). `functional` must carry the seated ring's A and C and the base; `processed` the
    // pins MKF placed. `terminalPlaneY` is the conductor builder's toroid terminal plane.
    static std::vector<NamedShape> buildBaseNamed(const MAS::BobbinFunctionalDescription& functional,
                                                  const MAS::CoreBobbinProcessedDescription& processed,
                                                  double terminalPlaneY, const std::string& bobbinName);
};

// Collision gate: a Base solid sharing more than 1 um^3 with any other solid (core, coating, copper,
// another accessory) throws naming the pair. The base's own pins only touch its bottom face.
void checkBaseCollisions(const std::vector<NamedShape>& all);

} // namespace mvb
