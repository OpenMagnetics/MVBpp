#pragma once
// Moved verbatim from OMFEM include/omfem/StepAudit.hpp (ABT #1588, step 3): namespace omfem::cad -> mvb::cad. No logic changed.
// THE COPPER OVERLAP AUDIT, as a function of the solids -- so the tool that has already imported
// a STEP can run it without importing it again (Alf, 2026-09-03: "we are losing a lot of time
// waiting for checks"). The audit and the CAD check each used to read the same 10-20 MB file
// through OCC's STEP reader, which is the single most expensive thing either of them does.
//
// The logic is unchanged and lives in StepAudit.cpp: never trust one boolean (classify the
// common's centroid against both parents, and reject a common whose bounding box escapes the
// operands'), decide a rejected common by surface point sampling, and separate a same-conductor
// junction from two different bodies interpenetrating.
#include <TopoDS_Shape.hxx>
#include <Bnd_Box.hxx>

#include <string>
#include <vector>

namespace mvb::cad {

struct AuditPart {
    TopoDS_Shape shape;
    std::string name;     // as named in the STEP, e.g. "Primary parallel 0 [solid 3]"
    std::string body;     // the conductor/part the solid belongs to
    Bnd_Box box;
    double volume = 0.0;  // mm^3
};

// Prints the same report the standalone tool prints. Returns 0 = no overlaps, 1 = overlaps,
// 3 = at least one pair could not be proven either way.
int audit_overlaps(std::vector<AuditPart>& parts, double toleranceMm3 = 1e-6);

}  // namespace mvb::cad
