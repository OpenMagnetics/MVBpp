#include "mvb/PinBuilder.h"

#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace mvb {

namespace {

// Degrees about X, Y then Z, applied about the pin's own centre. MKF emits (90, 0, 0)
// for a horizontal former, which takes the default -Y pin axis to -Z.
gp_Trsf rotation_of(const MAS::Pin& pin) {
    gp_Trsf combined;
    const auto rotation = pin.get_rotation();
    if (!rotation) return combined;
    const gp_Pnt origin(0.0, 0.0, 0.0);
    const gp_Dir axes[3] = {gp_Dir(1.0, 0.0, 0.0), gp_Dir(0.0, 1.0, 0.0), gp_Dir(0.0, 0.0, 1.0)};
    for (std::size_t axis = 0; axis < 3 && axis < rotation->size(); ++axis) {
        const double degrees = (*rotation)[axis];
        if (degrees == 0.0) continue;
        gp_Trsf step;
        step.SetRotation(gp_Ax1(origin, axes[axis]), degrees * std::numbers::pi / 180.0);
        combined.PreMultiply(step);
    }
    return combined;
}

} // namespace

TopoDS_Shape PinBuilder::buildPin(const MAS::Pin& pin) {
    const auto& dimensions = pin.get_dimensions();
    const std::string named = pin.get_name() ? (" '" + pin.get_name().value() + "'") : std::string();
    if (dimensions.size() < 3)
        throw std::runtime_error(
            "PinBuilder: pin" + named + " carries " + std::to_string(dimensions.size()) +
            " dimensions; a pin solid needs three ([width/diameter, depth, length]).");

    const double width  = dimensions[0];
    const double depth  = dimensions[1];
    const double length = dimensions[2];
    if (!(width > 0.0) || !(length > 0.0))
        throw std::runtime_error(
            "PinBuilder: pin" + named + " has width " + std::to_string(width) + " m and length " +
            std::to_string(length) + " m; both have to be positive to be a solid.");

    // Built centred on the origin and pointing along -Y, so a rotation about the centre
    // turns the pin without moving it, and the final translation is exactly MKF's
    // coordinate - the pin's CENTRE, as bobbin.json $defs/pin.coordinates defines it.
    TopoDS_Shape shape;
    if (pin.get_shape() == MAS::PinShape::RECTANGULAR) {
        if (!(depth > 0.0))
            throw std::runtime_error(
                "PinBuilder: rectangular pin" + named + " has depth " + std::to_string(depth) +
                " m; a box needs a positive depth.");
        shape = BRepPrimAPI_MakeBox(gp_Pnt(-width / 2.0, -length / 2.0, -depth / 2.0),
                                    width, length, depth).Shape();
    }
    else {
        // ROUND, and IRREGULAR: a shape MAS does not describe further is drawn as the
        // circular pin it almost always is, with its stated diameter, rather than
        // guessed at.
        shape = BRepPrimAPI_MakeCylinder(
                    gp_Ax2(gp_Pnt(0.0, length / 2.0, 0.0), gp_Dir(0.0, -1.0, 0.0)),
                    width / 2.0, length).Shape();
    }
    if (shape.IsNull())
        throw std::runtime_error("PinBuilder: OCCT produced no solid for pin" + named + ".");

    // A pin MKF could not place has no position to be drawn at: the origin would put a copper
    // rod through the core. That is a data fault to hear about, never a pin to drop silently.
    const auto coordinates = pin.get_coordinates();
    if (!coordinates || coordinates->size() < 3)
        throw std::runtime_error(
            "PinBuilder: pin" + named + " has " +
            (coordinates ? std::to_string(coordinates->size()) + " coordinates" : std::string("no coordinates")) +
            "; a pin in processedDescription.pins[] needs its [x, y, z] centre to be drawn.");
    gp_Trsf placement = rotation_of(pin);
    gp_Trsf translation;
    translation.SetTranslation(gp_Vec((*coordinates)[0], (*coordinates)[1], (*coordinates)[2]));
    placement.PreMultiply(translation);
    return BRepBuilderAPI_Transform(shape, placement, true).Shape();
}

std::vector<NamedShape> PinBuilder::buildPinsNamed(
    const MAS::CoreBobbinProcessedDescription& bobbinProcessedDescription,
    const std::string& bobbinName)
{
    std::vector<NamedShape> out;
    const auto pins = bobbinProcessedDescription.get_pins();
    if (!pins) return out;
    out.reserve(pins->size());
    for (std::size_t index = 0; index < pins->size(); ++index) {
        const auto& pin = (*pins)[index];
        // A name is the pin's identity in the STEP tree and in OMFEM's classifier; numbering an
        // unnamed pin here would invent one that MKF's pin assignment can never refer to.
        if (!pin.get_name() || pin.get_name()->empty())
            throw std::runtime_error(
                "PinBuilder: pin " + std::to_string(index) + " of bobbin '" + bobbinName +
                "' has no name; MKF names every pin it places.");
        out.emplace_back(buildPin(pin), bobbinName + " pin " + pin.get_name().value(), Role::Pin);
    }
    return out;
}

std::vector<TopoDS_Shape> PinBuilder::buildPins(
    const MAS::CoreBobbinProcessedDescription& bobbinProcessedDescription)
{
    std::vector<TopoDS_Shape> out;
    for (auto& named : buildPinsNamed(bobbinProcessedDescription, std::string("Bobbin")))
        out.push_back(named.shape);
    return out;
}

} // namespace mvb
