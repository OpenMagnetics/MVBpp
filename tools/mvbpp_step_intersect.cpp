// The overlap audit's front end: read the STEP, name the solids, hand them to the shared
// Moved verbatim from OMFEM tools/omfem_step_intersect.cpp at 00740b3 (ABT #1588, step 8); renamed; StepAudit.hpp -> mvb; omfem::cad:: -> mvb. No logic changed.
// implementation (mvb::cad::audit_overlaps). The logic itself moved into the library so the
// CAD check can run it on the solids it has ALREADY imported -- reading a 10-20 MB STEP twice
// was the most expensive thing the two checks did (Alf, 2026-09-03).
#include "mvb/mesh/StepAudit.h"

#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <TCollection_AsciiString.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDataStd_Name.hxx>
#include <TDocStd_Document.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Shape.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: mvbpp_step_intersect <file.step> [volumeTolerance_mm3] [--adjacent]\n";
        return 2;
    }
    // Default tolerance: 1e-6 mm^3 -- a cube 0.01 mm on a side. Below that a "common" solid is
    // numerical noise from surfaces that merely touch, which is contact, not interpenetration.
    const double toleranceMm3 = (argc > 2 && argv[2][0] != '-') ? std::atof(argv[2]) : 1e-6;
    for (int i = 2; i < argc; ++i)
        if (!std::strcmp(argv[i], "--adjacent")) setenv("OMFEM_AUDIT_ADJACENT", "1", 1);

    Handle(TDocStd_Document) doc;
    XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
    STEPCAFControl_Reader reader;
    reader.SetNameMode(true);
    if (reader.ReadFile(argv[1]) != IFSelect_RetDone || !reader.Transfer(doc)) {
        std::cerr << "ERROR: cannot read " << argv[1] << "\n";
        return 2;
    }
    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
    TDF_LabelSequence labels;
    shapeTool->GetFreeShapes(labels);

    std::vector<mvb::cad::AuditPart> parts;
    for (Standard_Integer i = 1; i <= labels.Length(); ++i) {
        std::string name = "shape" + std::to_string(i);
        Handle(TDataStd_Name) nm;
        if (labels.Value(i).FindAttribute(TDataStd_Name::GetID(), nm))
            name = std::string(TCollection_AsciiString(nm->Get()).ToCString());
        TopoDS_Shape shape = shapeTool->GetShape(labels.Value(i));
        if (shape.IsNull()) continue;
        int localIndex = 0;
        for (TopExp_Explorer explorer(shape, TopAbs_SOLID); explorer.More(); explorer.Next()) {
            mvb::cad::AuditPart part;
            part.shape = explorer.Current();
            part.body = name;
            part.name = name + (localIndex ? " [solid " + std::to_string(localIndex) + "]" : "");
            BRepBndLib::Add(part.shape, part.box);
            GProp_GProps props;
            BRepGProp::VolumeProperties(part.shape, props);
            part.volume = props.Mass();
            parts.push_back(std::move(part));
            ++localIndex;
        }
    }
    std::cout << "[intersect] " << argv[1] << ":\n";
    return mvb::cad::audit_overlaps(parts, toleranceMm3);
}
