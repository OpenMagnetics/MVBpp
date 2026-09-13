#include "mvb/BaseBuilder.h"
#include "mvb/PinBuilder.h"

#include "Definitions.h"
#include "support/Utils.h"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace mvb {

namespace {

// Same meaning as the shunt gate's: faces that merely touch share boolean noise, anything above a
// 1 um cube is material occupying the base.
constexpr double kSharedVolumeTolerance = 1e-18;
// Two lengths of one record that must agree (height vs standoff + pocket): the catalogue writes
// them to the micrometre, so agreement means equal to rounding, not "close".
constexpr double kRecordAgreement = 1e-9;

double volume_of(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}

std::string mm(double metres) {
    std::ostringstream s;
    s << metres * 1e3 << " mm";
    return s.str();
}

} // namespace

void BaseBuilder::requireDrawableMounting(const MAS::BobbinBase& base, const std::string& bobbinName) {
    if (base.get_mounting() == MAS::OrientationEnum::HORIZONTAL) return;
    if (!base.get_boat_width())
        throw std::runtime_error(
            "BaseBuilder: toroid base '" + bobbinName + "' is mounted vertically and states no boatWidth; "
            "drawing is not implemented for vertical bases without data (the slot that holds the ring on "
            "edge is not described, and inventing one would put a body through the ring).");
    throw std::runtime_error(
        "BaseBuilder: toroid base '" + bobbinName + "' is mounted vertically; drawing a vertical base (the "
        "boat that holds the ring on edge, board normal -Y under the standing ring since ABT #1248) is not "
        "implemented yet.");
}

std::vector<NamedShape> BaseBuilder::buildBaseNamed(const MAS::BobbinFunctionalDescription& functional,
                                                    const MAS::CoreBobbinProcessedDescription& processed,
                                                    double terminalPlaneY, const std::string& bobbinName) {
    if (functional.get_family() != MAS::BobbinFamily::T || !functional.get_base())
        throw std::runtime_error("BaseBuilder: bobbin '" + bobbinName + "' is not a toroid base (family t with "
                                 "functionalDescription.base).");
    const MAS::BobbinBase base = functional.get_base().value();
    requireDrawableMounting(base, bobbinName);

    const auto ring = OpenMagnetics::flatten_dimensions(functional.get_dimensions());
    if (!ring.count("C") || !(ring.at("C") > 0))
        throw std::runtime_error("BaseBuilder: toroid base '" + bobbinName + "' has no seated ring height 'C'; "
                                 "seat the core with MKF's Bobbin::create_toroid_bobbin_on_base first.");
    const double ringHalfHeight = ring.at("C") / 2.0;
    const double length = OpenMagnetics::resolve_dimensional_values(base.get_length());
    const double width = OpenMagnetics::resolve_dimensional_values(base.get_width());
    const double height = OpenMagnetics::resolve_dimensional_values(base.get_height());
    const double standoff = OpenMagnetics::resolve_dimensional_values(base.get_standoff());
    if (!(length > 0) || !(width > 0) || !(height > 0))
        throw std::runtime_error("BaseBuilder: toroid base '" + bobbinName + "' states length " + mm(length) +
                                 ", width " + mm(width) + ", height " + mm(height) + "; all must be positive.");

    const bool pocketDiameterStated = base.get_pocket_inner_diameter().has_value();
    const bool pocketDepthStated = base.get_pocket_depth().has_value();
    if (pocketDiameterStated != pocketDepthStated)
        throw std::runtime_error("BaseBuilder: toroid base '" + bobbinName + "' states a pocket " +
                                 (pocketDiameterStated ? "diameter but no depth" : "depth but no diameter") +
                                 "; a pocket needs both to be cut.");
    const double pocketDepth = pocketDepthStated ? OpenMagnetics::resolve_dimensional_values(base.get_pocket_depth().value()) : 0.0;
    if (std::abs(height - pocketDepth - standoff) > kRecordAgreement)
        throw std::runtime_error(
            "BaseBuilder: toroid base '" + bobbinName + "' has height " + mm(height) +
            (pocketDepthStated ? " minus pocket depth " + mm(pocketDepth) : std::string(" and no pocket")) +
            ", which does not equal its standoff " + mm(standoff) +
            ". The base's top face (or pocket floor) is where the wound part rests, so these must agree; "
            "posts, feet or a seat the record does not describe cannot be drawn.");

    // The top face is the terminal plane, which lies below the ring's underside.
    if (!(terminalPlaneY < -ringHalfHeight))
        throw std::runtime_error("BaseBuilder: the toroid terminal plane y = " + mm(terminalPlaneY) +
                                 " is not below the ring's underside y = " + mm(-ringHalfHeight) +
                                 "; the base '" + bobbinName + "' would cut the ring.");
    const double top = terminalPlaneY;
    const double bottom = top - height;

    // Pins: MKF's, moved along Y so their tops sit on the bottom face.
    const auto pins = processed.get_pins();
    if (!pins || pins->empty())
        throw std::runtime_error("BaseBuilder: toroid base '" + bobbinName + "' carries no pins in "
                                 "processedDescription.pins[]; MKF places them when the ring is seated.");
    const double mkfSeatingPlane = -(ringHalfHeight + standoff);
    std::vector<NamedShape> out;
    std::vector<NamedShape> pinSolids;
    for (const auto& pin : pins.value()) {
        const std::string pinName = pin.get_name().value_or("");
        if (pinName.empty())
            throw std::runtime_error("BaseBuilder: a pin of toroid base '" + bobbinName + "' has no name.");
        if (pin.get_rotation())
            throw std::runtime_error("BaseBuilder: pin '" + pinName + "' of the horizontal toroid base '" + bobbinName +
                                     "' carries a rotation; a horizontal base's pins hang straight down -Y.");
        const auto coordinates = pin.get_coordinates();
        const auto& dimensions = pin.get_dimensions();
        if (!coordinates || coordinates->size() < 3 || dimensions.size() < 3)
            throw std::runtime_error("BaseBuilder: pin '" + pinName + "' of toroid base '" + bobbinName +
                                     "' has no [x, y, z] centre or fewer than three dimensions.");
        const double pinTop = (*coordinates)[1] + dimensions[2] / 2.0;
        if (std::abs(pinTop - mkfSeatingPlane) > kRecordAgreement)
            throw std::runtime_error("BaseBuilder: pin '" + pinName + "' of toroid base '" + bobbinName + "' starts at y = " +
                                     mm(pinTop) + ", not on MKF's seating plane -(C/2 + standoff) = " + mm(mkfSeatingPlane) +
                                     "; the pins and the base record disagree.");
        // Inside the footprint: width along X, length along Z (see the header).
        const double halfAlongX = dimensions[0] / 2.0;   // PinBuilder: box width along X, or the diameter
        const double halfAlongZ = pin.get_shape() == MAS::PinShape::RECTANGULAR ? dimensions[1] / 2.0 : dimensions[0] / 2.0;
        if (std::abs((*coordinates)[0]) + halfAlongX > width / 2.0 + kRecordAgreement ||
            std::abs((*coordinates)[2]) + halfAlongZ > length / 2.0 + kRecordAgreement)
            throw std::runtime_error("BaseBuilder: pin '" + pinName + "' of toroid base '" + bobbinName + "' at (x, z) = (" +
                                     mm((*coordinates)[0]) + ", " + mm((*coordinates)[2]) +
                                     ") lies outside the base's footprint (width " + mm(width) + " along X, length " +
                                     mm(length) + " along Z).");
        gp_Trsf shift;
        shift.SetTranslation(gp_Vec(0.0, bottom - pinTop, 0.0));
        TopoDS_Shape solid = BRepBuilderAPI_Transform(PinBuilder::buildPin(pin), shift, true).Shape();
        pinSolids.emplace_back(solid, bobbinName + " pin " + pinName, Role::Pin);
    }

    TopoDS_Shape body = BRepPrimAPI_MakeBox(gp_Pnt(-width / 2.0, bottom, -length / 2.0), width, height, length).Shape();
    if (pocketDiameterStated) {
        const double pocketDiameter = OpenMagnetics::resolve_dimensional_values(base.get_pocket_inner_diameter().value());
        if (!(pocketDiameter > 0) || !(pocketDepth > 0) || !(pocketDepth < height) || pocketDiameter >= std::min(width, length))
            throw std::runtime_error("BaseBuilder: toroid base '" + bobbinName + "' states a pocket of diameter " +
                                     mm(pocketDiameter) + " and depth " + mm(pocketDepth) + " that does not fit its " +
                                     mm(width) + " x " + mm(length) + " x " + mm(height) + " body.");
        TopoDS_Shape pocket = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0.0, top - pocketDepth, 0.0), gp_Dir(0.0, 1.0, 0.0)),
                                                       pocketDiameter / 2.0, pocketDepth).Shape();
        BRepAlgoAPI_Cut cut(body, pocket);
        if (!cut.IsDone())
            throw std::runtime_error("BaseBuilder: cutting the pocket of toroid base '" + bobbinName + "' failed (OCCT boolean).");
        body = cut.Shape();
    }
    out.emplace_back(body, bobbinName + " base", Role::Base);
    for (auto& p : pinSolids) out.push_back(std::move(p));
    return out;
}

void checkBaseCollisions(const std::vector<NamedShape>& all) {
    for (const auto& base : all) {
        if (base.role != Role::Base || base.shape.IsNull()) continue;
        Bnd_Box baseBox;
        BRepBndLib::Add(base.shape, baseBox);
        for (const auto& other : all) {
            if (&other == &base || other.shape.IsNull() || other.role == Role::Terminal) continue;
            for (TopExp_Explorer ex(other.shape, TopAbs_SOLID); ex.More(); ex.Next()) {
                Bnd_Box box;
                BRepBndLib::Add(ex.Current(), box);
                if (box.IsOut(baseBox)) continue;
                BRepAlgoAPI_Common common(base.shape, ex.Current());
                if (!common.IsDone())
                    throw std::runtime_error("Base collision gate: could not intersect " + base.name + " with " +
                                             other.name + " (OCCT boolean failed).");
                const double shared = volume_of(common.Shape());
                if (shared > kSharedVolumeTolerance) {
                    std::ostringstream msg;
                    msg << "Base collision gate: " << base.name << " intersects " << other.name << " ("
                        << role_name(other.role) << "), sharing " << shared * 1e9 << " mm^3.";
                    throw std::runtime_error(msg.str());
                }
            }
        }
    }
}

} // namespace mvb
