#include "mvb/shapes/ShapeBuilder.h"
#include "mvb/shapes/ShapeC.h"
#include "mvb/shapes/ShapeDrum.h"
#include "mvb/shapes/ShapeDrumRing.h"
#include "mvb/shapes/ShapeMolded.h"
#include "mvb/shapes/ShapeDrumSemishielded.h"
#include "mvb/shapes/ShapeE.h"
#include "mvb/shapes/ShapeEp.h"
#include "mvb/shapes/ShapeLp.h"
#include "mvb/shapes/ShapeEpx.h"
#include "mvb/shapes/ShapeEr.h"
#include "mvb/shapes/ShapeEc.h"
#include "mvb/shapes/ShapeP.h"
#include "mvb/shapes/ShapePm.h"
#include "mvb/shapes/ShapePQ.h"
#include "mvb/shapes/ShapeRM.h"
#include "mvb/shapes/ShapeSlab.h"
#include "mvb/shapes/ShapeT.h"
#include "mvb/shapes/ShapeU.h"
#include "mvb/shapes/ShapeUr.h"
#include "mvb/shapes/ShapeUt.h"
#include <memory>

namespace mvb {
namespace shapes {

std::unique_ptr<ShapeBuilder> createShapeBuilder(MAS::CoreShapeFamily family,
                                                  const std::string& subtype,
                                                  int corePolygonSegments) {
    std::unique_ptr<ShapeBuilder> builder;
    switch (family) {
        case MAS::CoreShapeFamily::E:
        case MAS::CoreShapeFamily::EFD:
        case MAS::CoreShapeFamily::EI:
        case MAS::CoreShapeFamily::EL:
        case MAS::CoreShapeFamily::ELP:
        case MAS::CoreShapeFamily::EF:  // ABT #266: EF is an E-core variant
        case MAS::CoreShapeFamily::PLANAR_E:
        case MAS::CoreShapeFamily::PLANAR_EL:
            builder = std::make_unique<ShapeE>(); break;
        case MAS::CoreShapeFamily::ER:
        case MAS::CoreShapeFamily::ETD:
        case MAS::CoreShapeFamily::EQ:
        case MAS::CoreShapeFamily::EER:  // ABT #265: ER-type with round centre leg
        case MAS::CoreShapeFamily::PLANAR_ER:
            builder = std::make_unique<ShapeEr>(); break;
        case MAS::CoreShapeFamily::DRUM_SEMISHIELDED:  // ABT #362: drum; the epoxy shell is a coating
            builder = std::make_unique<ShapeDrumSemishielded>();
            break;
        case MAS::CoreShapeFamily::MOLDED:  // ABT #357: SMC block moulded around the coil
            builder = std::make_unique<ShapeMolded>();
            break;
        case MAS::CoreShapeFamily::T:
            builder = std::make_unique<ShapeT>(); break;
        case MAS::CoreShapeFamily::DRUM:
            builder = std::make_unique<ShapeDrum>();
            break;
        case MAS::CoreShapeFamily::DRUM_RING:  // ABT #366/#367: drum + concentric shield ring
            builder = std::make_unique<ShapeDrumRing>();
            break;
        case MAS::CoreShapeFamily::U:
            builder = std::make_unique<ShapeU>(); break;
        case MAS::CoreShapeFamily::UR:
            builder = std::make_unique<ShapeUr>(); break;
        case MAS::CoreShapeFamily::UT:
            builder = std::make_unique<ShapeUt>(); break;
        case MAS::CoreShapeFamily::C:
            builder = std::make_unique<ShapeC>(); break;
        case MAS::CoreShapeFamily::P:
            builder = std::make_unique<ShapeP>(); break;
        case MAS::CoreShapeFamily::PM:
            builder = std::make_unique<ShapePm>(subtype); break;
        case MAS::CoreShapeFamily::EP:
        case MAS::CoreShapeFamily::EPC:  // ABT #270: EP derivatives share EP geometry
        case MAS::CoreShapeFamily::EPQ:  // ABT #268
        case MAS::CoreShapeFamily::EPW:  // ABT #269
        case MAS::CoreShapeFamily::EPT:  // ABT #271
        case MAS::CoreShapeFamily::LEP:  // ABT #272
            builder = std::make_unique<ShapeEp>(); break;
        case MAS::CoreShapeFamily::EPX:
            builder = std::make_unique<ShapeEpx>(); break;
        case MAS::CoreShapeFamily::LP:
            builder = std::make_unique<ShapeLp>(); break;
        case MAS::CoreShapeFamily::EC:
            builder = std::make_unique<ShapeEc>(); break;
        case MAS::CoreShapeFamily::PQ:
        case MAS::CoreShapeFamily::PQI:
            builder = std::make_unique<ShapePQ>(); break;
        case MAS::CoreShapeFamily::RM:
            builder = std::make_unique<ShapeRM>(); break;
        // Slab cores (MKF CorePieceSlab, ABT #263): pot cores with two flats. RS is NOT here:
        // its set pairs a slab half with a plain pot round, and the geometricalDescription
        // MKF emits carries two identical RS halfSets, so nothing says which half is round.
        case MAS::CoreShapeFamily::DS:
        case MAS::CoreShapeFamily::HS:
            builder = std::make_unique<ShapeSlab>(); break;
        default:
            return nullptr;
    }
    if (builder) builder->setCorePolygonSegments(corePolygonSegments);
    return builder;
}

} // namespace shapes

std::vector<std::string> get_supported_families() {
    // Derive from the factory: try every CoreShapeFamily and keep the ones
    // that produce a non-null builder.
    std::vector<std::string> result;
    for (int i = 0; i <= static_cast<int>(MAS::CoreShapeFamily::UT); ++i) {
        auto family = static_cast<MAS::CoreShapeFamily>(i);
        auto builder = shapes::createShapeBuilder(family);
        if (builder) {
            result.push_back(core_shape_family_to_string(family));
        }
    }
    return result;
}

} // namespace mvb
