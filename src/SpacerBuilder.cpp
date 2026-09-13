#include "mvb/SpacerBuilder.h"
#include "mvb/Utils.h"

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pnt.hxx>

#include <stdexcept>
#include <string>

namespace mvb {

namespace {

std::string insulation_material_name(const MAS::CoreGeometricalDescriptionElement& part,
                                     std::size_t index) {
    auto materialOpt = part.get_insulation_material();
    if (!materialOpt) {
        throw std::runtime_error(
            "SpacerBuilder: spacer " + std::to_string(index) +
            " has no insulationMaterial. MAS core/spacer.json REQUIRES it; a shim with no "
            "declared material cannot be given a dielectric constant or a thermal "
            "conductivity downstream.");
    }
    if (const auto* name = std::get_if<std::string>(&*materialOpt)) {
        return *name;
    }
    const auto& record = std::get<MAS::InsulationMaterial>(*materialOpt);
    // MAS models InsulationMaterial::name as REQUIRED, so this is a plain string, not an
    // optional; an empty one is still a record that cannot be looked up downstream.
    const std::string& name = record.get_name();
    if (name.empty()) {
        throw std::runtime_error(
            "SpacerBuilder: spacer " + std::to_string(index) +
            " carries an inline insulationMaterial record with no name.");
    }
    return name;
}

}  // namespace

std::vector<Spacer> SpacerBuilder::buildSpacers(
    const std::vector<MAS::CoreGeometricalDescriptionElement>& geometricalDescription)
{
    std::vector<Spacer> out;
    std::size_t index = 0;
    for (const auto& part : geometricalDescription) {
        if (part.get_type() != MAS::CoreGeometricalDescriptionElementType::SPACER) continue;

        // MAS core/spacer.json: dimensions = [X, Y, Z] of the cube, Y being the column axis
        // (the shim thickness = the additive gap); coordinates = the CENTRE of that cube.
        auto dimsOpt = part.get_dimensions();
        const auto& coords = part.get_coordinates();
        if (!dimsOpt || dimsOpt->size() < 3) {
            throw std::runtime_error(
                "SpacerBuilder: spacer " + std::to_string(index) +
                " has no three-element `dimensions`; MAS core/spacer.json requires them and "
                "there is no sane box to invent.");
        }
        if (coords.size() < 3) {
            throw std::runtime_error(
                "SpacerBuilder: spacer " + std::to_string(index) +
                " has fewer than three `coordinates`.");
        }
        const auto& dims = *dimsOpt;

        const double width     = dims[0];  // X
        const double thickness = dims[1];  // Y — the column axis, i.e. the additive gap
        const double depth     = dims[2];  // Z

        if (width <= 0.0 || thickness <= 0.0 || depth <= 0.0) {
            throw std::runtime_error(
                "SpacerBuilder: spacer " + std::to_string(index) +
                " has a non-positive dimension [" + std::to_string(width) + ", " +
                std::to_string(thickness) + ", " + std::to_string(depth) + "].");
        }

        const std::string materialName = insulation_material_name(part, index);

        gp_Pnt corner(coords[0] - width / 2.0, coords[1] - thickness / 2.0, coords[2] - depth / 2.0);
        TopoDS_Shape box = BRepPrimAPI_MakeBox(corner, width, thickness, depth).Shape();
        if (box.IsNull()) {
            throw std::runtime_error(
                "SpacerBuilder: OCCT returned a null solid for spacer " + std::to_string(index));
        }
        out.push_back(Spacer{box, materialName});
        ++index;
    }
    return out;
}

TopoDS_Shape SpacerBuilder::buildSpacersCompound(
    const std::vector<MAS::CoreGeometricalDescriptionElement>& geometricalDescription)
{
    auto spacers = buildSpacers(geometricalDescription);
    if (spacers.empty()) return TopoDS_Shape();
    if (spacers.size() == 1) return spacers.front().shape;

    BRep_Builder b;
    TopoDS_Compound compound;
    b.MakeCompound(compound);
    for (const auto& s : spacers) b.Add(compound, s.shape);
    return compound;
}

void appendSpacerSolids(std::vector<NamedShape>& all, const MAS::MagneticCore& core) {
    auto geometricalDescription = core.get_geometrical_description();
    if (!geometricalDescription) return;   // un-enriched core: nothing to draw, not an error
    auto spacers = SpacerBuilder::buildSpacers(*geometricalDescription);
    for (std::size_t i = 0; i < spacers.size(); ++i) {
        NamedShape ns{spacers[i].shape, "Spacer_" + std::to_string(i)};
        // WP0 (ABT #1169) added Role to NamedShape and made OMFEM's classify() throw on an
        // unrecognised name; `Spacer_` is its spacer branch.
        ns.role = Role::Spacer;   // mvb::Role, added by WP0
        ns.materialName = spacers[i].insulationMaterialName;
        all.push_back(std::move(ns));
    }
}

} // namespace mvb
