#include "mvb/ShuntBuilder.h"

#include "constructive_models/Magnetic.h"
#include "physical_models/MagneticShunt.h"
#include "Defaults.h"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>

namespace mvb {

namespace {

// A shared volume below this is boolean noise on faces that merely touch (a sheet laid against a
// cut bobbin face, a shim face on a sheet face); anything above it is copper, ferrite or plastic
// occupying the sheet. 1e-18 m^3 is a cube of 1 um side.
constexpr double kSharedVolumeTolerance = 1e-18;

double volume_of(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}

std::string shunt_label(const MAS::MagneticShunt& shunt, std::size_t index) {
    std::string label = "shunt " + std::to_string(index);
    if (shunt.get_name()) label += " ('" + shunt.get_name().value() + "')";
    return label;
}

std::string material_name(const MAS::MagneticShunt& shunt, std::size_t index) {
    const auto& material = shunt.get_material();
    if (const auto* name = std::get_if<std::string>(&material)) {
        if (name->empty())
            throw std::runtime_error("ShuntBuilder: " + shunt_label(shunt, index) +
                                     " names an empty material.");
        return *name;
    }
    const std::string& name = std::get<MAS::CoreMaterial>(material).get_name();
    if (name.empty())
        throw std::runtime_error("ShuntBuilder: " + shunt_label(shunt, index) +
                                 " carries an inline core-material record with no name.");
    return name;
}

// MKF's validation, exactly the calls its own models make. Throws whatever MKF throws.
void validate_with_mkf(const OpenMagnetics::Magnetic& magnetic) {
    OpenMagnetics::MagneticShuntModel::check_supported_placements(magnetic);
    // The generated getter returns the optional BY VALUE: bind a copy, never a reference.
    const std::vector<MAS::MagneticShunt> shunts = magnetic.get_shunts().value();
    // Only the core-loop reluctance of the network depends on these, and it is discarded here:
    // the call is made for its geometric checks.
    const double temperature = OpenMagnetics::Defaults().ambientTemperature;
    const double frequency = OpenMagnetics::Defaults().measurementFrequency;

    std::vector<MAS::MagneticShunt> onColumn;
    OpenMagnetics::Magnetic copy = magnetic;
    for (std::size_t i = 0; i < shunts.size(); ++i) {
        if (OpenMagnetics::MagneticShuntModel::is_leakage_shunt(shunts[i])) {
            OpenMagnetics::MagneticShuntModel::extract_network_inputs(copy, i, temperature, frequency);
        }
        else {
            onColumn.push_back(shunts[i]);   // outsideWindow already threw above
        }
    }
    if (!onColumn.empty()) {
        // apply_shunts_to_gapping walks every shunt it is given; hand it the onColumn ones only,
        // so a leakage sheet crossing a column gap is not asked for a gap-filler permeability.
        OpenMagnetics::Magnetic columnOnly = magnetic;
        columnOnly.set_shunts(onColumn);
        OpenMagnetics::MagneticShuntModel::apply_shunts_to_gapping(columnOnly, temperature, frequency);
    }
}

TopoDS_Shape make_box(double x0, double y0, double z0, double dx, double dy, double dz,
                      const std::string& what) {
    TopoDS_Shape box = BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), dx, dy, dz).Shape();
    if (box.IsNull()) throw std::runtime_error("ShuntBuilder: OCCT returned a null box for " + what);
    return box;
}

}  // namespace

std::vector<NamedShape> ShuntBuilder::buildShuntsNamed(const OpenMagnetics::Magnetic& magnetic) {
    std::vector<NamedShape> out;
    if (!magnetic.get_shunts() || magnetic.get_shunts()->empty()) return out;
    validate_with_mkf(magnetic);

    // The generated getter returns the optional BY VALUE: bind a copy, never a reference.
    const std::vector<MAS::MagneticShunt> shunts = magnetic.get_shunts().value();
    for (std::size_t i = 0; i < shunts.size(); ++i) {
        const auto& shunt = shunts[i];
        // MKF validated these (size and sign) for in-window sheets; an onColumn sheet went through
        // apply_shunts_to_gapping, which reads the same box, so the same guarantees hold.
        const auto& c = shunt.get_coordinates();
        const auto& d = shunt.get_dimensions();
        if (c.size() < 2 || d.size() < 3)
            throw std::runtime_error("ShuntBuilder: " + shunt_label(shunt, i) +
                                     " needs x and y coordinates and three dimensions.");
        const double width = d[0], height = d[1], depth = d[2];
        if (!(width > 0.0) || !(height > 0.0) || !(depth > 0.0))
            throw std::runtime_error("ShuntBuilder: " + shunt_label(shunt, i) +
                                     " has a non-positive dimension.");
        const double z = c.size() > 2 ? c[2] : 0.0;
        const double x0 = c[0] - width / 2.0;
        const double y0 = c[1] - height / 2.0;
        const double z0 = z - depth / 2.0;
        const std::string material = material_name(shunt, i);
        const std::string base = "Shunt_" + std::to_string(i);

        if (!shunt.get_segments() || shunt.get_segments()->empty()) {
            NamedShape ns{make_box(x0, y0, z0, width, height, depth, base), base, Role::Shunt};
            ns.materialName = material;
            out.push_back(std::move(ns));
            continue;
        }

        double x = x0;
        const std::vector<MAS::MagneticShuntSegment> segments = shunt.get_segments().value();
        for (std::size_t k = 0; k < segments.size(); ++k) {
            // A missing gap is 0, as in MKF's sum (MKF throws for a negative one). A missing length
            // has no reading at all.
            if (!segments[k].get_length())
                throw std::runtime_error("ShuntBuilder: " + shunt_label(shunt, i) + " segment " +
                                         std::to_string(k) + " has no length.");
            const double length = segments[k].get_length().value();
            const double gap = segments[k].get_gap().value_or(0.0);
            if (!(length > 0.0) || gap < 0.0)
                throw std::runtime_error("ShuntBuilder: " + shunt_label(shunt, i) + " segment " +
                                         std::to_string(k) + " has a non-positive length or a negative gap.");
            const std::string name = base + "_" + std::to_string(k);
            NamedShape ns{make_box(x, y0, z0, length, height, depth, name), name, Role::Shunt};
            ns.materialName = material;
            out.push_back(std::move(ns));
            x += length + gap;
        }
        // extract_network_inputs made this check for an in-window sheet; apply_shunts_to_gapping
        // (onColumn) reads no segments, so repeat it with MKF's own tolerance rather than draw
        // pieces past the box.
        if (std::fabs((x - x0) - width) > std::max(1e-6, 0.01 * width)) {
            std::ostringstream msg;
            msg << "ShuntBuilder: " << shunt_label(shunt, i) << " segments add up to " << (x - x0)
                << " m but the sheet is " << width << " m wide.";
            throw std::runtime_error(msg.str());
        }
    }
    return out;
}

std::vector<NamedShape> ShuntBuilder::buildShuntsNamed(const MAS::Magnetic& magnetic) {
    if (!magnetic.get_shunts() || magnetic.get_shunts()->empty()) return {};
    // The MKF model needs the processed core (columns, gap coordinates); the MAS -> MKF conversion
    // keeps the shunts (ABT #1176).
    return buildShuntsNamed(OpenMagnetics::Magnetic(magnetic));
}

void cutByShunts(NamedShape& host, const std::vector<NamedShape>& all) {
    if (host.shape.IsNull()) return;
    Bnd_Box hostBox;
    BRepBndLib::Add(host.shape, hostBox);
    for (const auto& shunt : all) {
        if (shunt.role != Role::Shunt || &shunt == &host) continue;
        Bnd_Box shuntBox;
        BRepBndLib::Add(shunt.shape, shuntBox);
        if (hostBox.IsOut(shuntBox)) continue;
        BRepAlgoAPI_Common common(host.shape, shunt.shape);
        if (!common.IsDone())
            throw std::runtime_error("ShuntBuilder: could not intersect " + shunt.name + " with " +
                                     host.name + " (OCCT boolean failed).");
        const double shared = volume_of(common.Shape());
        if (shared <= kSharedVolumeTolerance) continue;
        const double before = volume_of(host.shape);
        BRepAlgoAPI_Cut cut(host.shape, shunt.shape);
        if (!cut.IsDone() || cut.Shape().IsNull())
            throw std::runtime_error("ShuntBuilder: cutting " + shunt.name + " out of " + host.name +
                                     " failed (OCCT boolean).");
        const double after = volume_of(cut.Shape());
        if (std::fabs(before - shared - after) > 1e-9 * before) {
            std::ostringstream msg;
            msg << "ShuntBuilder: cutting " << shunt.name << " out of " << host.name << " removed "
                << (before - after) << " m^3 but the two share " << shared << " m^3.";
            throw std::runtime_error(msg.str());
        }
        host.shape = cut.Shape();
    }
}

namespace {

void append_and_cut_spacers(std::vector<NamedShape>& all, std::vector<NamedShape> shunts) {
    if (shunts.empty()) return;
    for (auto& s : shunts) all.push_back(std::move(s));
    for (auto& ns : all)
        if (ns.role == Role::Spacer) cutByShunts(ns, all);
}

}  // namespace

void appendShuntSolids(std::vector<NamedShape>& all, const OpenMagnetics::Magnetic& magnetic) {
    append_and_cut_spacers(all, ShuntBuilder::buildShuntsNamed(magnetic));
}

void appendShuntSolids(std::vector<NamedShape>& all, const MAS::Magnetic& magnetic) {
    append_and_cut_spacers(all, ShuntBuilder::buildShuntsNamed(magnetic));
}

void checkShuntCollisions(const std::vector<NamedShape>& all) {
    std::vector<const NamedShape*> shunts;
    for (const auto& ns : all)
        if (ns.role == Role::Shunt && !ns.shape.IsNull()) shunts.push_back(&ns);
    if (shunts.empty()) return;

    for (std::size_t s = 0; s < shunts.size(); ++s) {
        const NamedShape& shunt = *shunts[s];
        Bnd_Box shuntBox;
        BRepBndLib::Add(shunt.shape, shuntBox);
        for (const auto& other : all) {
            if (&other == &shunt || other.shape.IsNull() || other.role == Role::Terminal) continue;
            // Each shunt pair once.
            if (other.role == Role::Shunt) {
                bool earlier = false;
                for (std::size_t t = 0; t <= s; ++t) earlier = earlier || shunts[t] == &other;
                if (earlier) continue;
            }
            // Per solid: a real-winding conductor is a compound of many, and only the few near
            // the sheet are worth a boolean.
            for (TopExp_Explorer ex(other.shape, TopAbs_SOLID); ex.More(); ex.Next()) {
                Bnd_Box box;
                BRepBndLib::Add(ex.Current(), box);
                if (box.IsOut(shuntBox)) continue;
                BRepAlgoAPI_Common common(shunt.shape, ex.Current());
                if (!common.IsDone())
                    throw std::runtime_error("Shunt collision gate: could not intersect " + shunt.name +
                                             " with " + other.name + " (OCCT boolean failed).");
                const double shared = volume_of(common.Shape());
                if (shared > kSharedVolumeTolerance) {
                    std::ostringstream msg;
                    msg << "Shunt collision gate: " << shunt.name << " intersects " << other.name
                        << " (" << role_name(other.role) << "), sharing " << shared * 1e9
                        << " mm^3. A magnetic shunt must not occupy copper, core, insulation or "
                           "another shunt; fix the shunt coordinates/dimensions or the winding.";
                    throw std::runtime_error(msg.str());
                }
            }
        }
    }
}

}  // namespace mvb
