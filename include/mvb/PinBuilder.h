#pragma once

#include "MAS.hpp"
#include "mvb/NamedShape.h"

#include <TopoDS_Shape.hxx>
#include <string>
#include <vector>

namespace mvb {

// ABT #1171 / WP2: the solder pins of a bobbin.
//
// Every number comes from MKF: `Bobbin::expand_pinout` turns a catalogue `pinout`
// (pin count, rows, pitches) into `processedDescription.pins[]`, each pin carrying its
// coordinates in the bobbin frame, its shape, its dimensions and, for a horizontal
// former, its rotation. This builder draws those pins and nothing else - it does not
// re-derive a pitch, a row distance or a protrusion, and it does not invent a pin for a
// bobbin whose record only states how many pins it has (MKF leaves `pins` absent there,
// and absent means nothing to draw).
class PinBuilder {
public:
    // One solid per pin, named "<bobbinName> pin <pin name>".
    //
    // `bobbinName` is the name the bobbin solid itself carries (MagneticBuilder's
    // getBobbinNameT result), so a pin and its former share a prefix in the STEP tree.
    //
    // Absent `pins` means nothing to draw. But once `pins` is present every entry must be
    // drawable, so this THROWS (never skips) for a pin with no name, no [x, y, z]
    // coordinates, or dimensions that cannot describe a solid (fewer than three, or a
    // non-positive diameter or length): a pin MKF placed and MVB++ silently dropped is a
    // missing terminal in the real part.
    static std::vector<NamedShape> buildPinsNamed(
        const MAS::CoreBobbinProcessedDescription& bobbinProcessedDescription,
        const std::string& bobbinName);

    // The same pins as bare shapes, for the bobbin cutter list.
    static std::vector<TopoDS_Shape> buildPins(
        const MAS::CoreBobbinProcessedDescription& bobbinProcessedDescription);

    // One pin. Exposed for tests, which check a single pin's placement against MKF's
    // coordinate without building a whole magnetic.
    static TopoDS_Shape buildPin(const MAS::Pin& pin);
};

} // namespace mvb
