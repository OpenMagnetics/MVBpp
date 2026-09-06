#include "mvb/StepExporter.h"

#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepCheck.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepTools.hxx>
#include <BRepTools_Modification.hxx>
#include <BRepTools_Modifier.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2dConvert.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomLib.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_DataMap.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopAbs.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_GTrsf2d.hxx>
#include <gp_Mat2d.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Bnd_Box.hxx>
#include <gp_Trsf.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <StlAPI_Writer.hxx>
#include <TDF_Label.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDataStd_Name.hxx>
#include <TDocStd_Application.hxx>
#include <TDocStd_Document.hxx>
#include <TCollection_ExtendedString.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <random>
#include <sstream>
#include <stdexcept>

namespace mvb {

namespace {

std::string label_name(const TDF_Label& label) {
    Handle(TDataStd_Name) attr;
    if (label.FindAttribute(TDataStd_Name::GetID(), attr)) {
        const TCollection_ExtendedString& n = attr->Get();
        std::string out;
        out.reserve(n.Length());
        for (int i = 1; i <= n.Length(); ++i) {
            out.push_back(static_cast<char>(n.Value(i)));
        }
        return out;
    }
    return {};
}


// FEM product only: re-express every solid carrying a PERIODIC surface as B-splines (ABT #490
// class, 2026-08-25; moved here from buildAllNamed 2026-09-06, ABT #1111).
//
// gmsh decides whether to send a face to its markedly more fragile periodic mesher by reading
// the UNDERLYING SURFACE, not the face -- OCCFace.cpp:142, `_periodic[0] = surface.IsUPeriodic()`
// -- so a 60-degree panel of a cylinder, or a 20-degree slice of a torus, is treated exactly like
// a full closed one and routed to meshGeneratorPeriodic. That path produced 16 of the 23 corpus
// mesh failures measured on 2026-08-24. A conic IS exactly representable as a rational B-spline,
// so this conversion MOVES NO POINT: it only drops the periodicity flag on faces that do not
// actually wrap.
//
// WHY IT RUNS ON THE MILLIMETRE SHAPE, INSIDE THE EXPORTER. A B-spline's knots are not scaled by
// a transform; its poles are. Converting in metres (as buildAllNamed did) and then scaling x1000
// for the STEP left every length-parametrised direction (a facet strip's generatrix, a planar
// cap) with a knot span 1000x smaller than the geometry it spans, i.e. |dS/du| = 1000 mm per
// parameter unit. BOPAlgo's face/face work converts 3D tolerances to parameter space through
// exactly that ratio, and on faceted revolves it then found sporadic self-intersections between a
// strip and its cap: 11 of the 38 corpus designs at --segments 12, one solid each, every one
// clean again once rescaled to metres (12_boost solid 51, single_switch solid 511). Angle-
// parametrised surfaces (a torus in both directions) are immune, which is why exact revolves
// never showed it. Converting after the scale keeps knots and geometry in the same unit.
//
// ONLY SOLIDS THAT ACTUALLY CARRY A PERIODIC FACE -- per SOLID, not per part. NurbsConvert
// rewrites EVERY surface, so run unconditionally it re-expressed 9579 PLANAR faces of the
// faceted 14_dab as B-splines; per-PART gating was still too coarse (one periodic lead pipe
// condemned every all-planar prism of the conductor, and two abutting prisms' independently
// converted caps are what the fragment then imprints as nm slivers). A planar solid stays
// analytic: its junction caps remain the same exact planes as its neighbour's.
TopoDS_Shape nurbsConvertPeriodicSolidsMm(const TopoDS_Shape& shapeMm, const std::string& name) {
    auto hasPeriodicFace = [](const TopoDS_Shape& s) {
        for (TopExp_Explorer fx(s, TopAbs_FACE); fx.More(); fx.Next()) {
            Handle(Geom_Surface) srf = BRep_Tool::Surface(TopoDS::Face(fx.Current()));
            if (!srf.IsNull() && (srf->IsUPeriodic() || srf->IsVPeriodic())) return true;
        }
        return false;
    };
    if (!hasPeriodicFace(shapeMm)) return shapeMm;
    const bool wasValid = BRepCheck_Analyzer(shapeMm).IsValid();
    try {
        TopoDS_Shape out;
        if (shapeMm.ShapeType() == TopAbs_COMPOUND) {
            TopoDS_Compound rebuilt;
            BRep_Builder rb;
            rb.MakeCompound(rebuilt);
            for (TopoDS_Iterator it(shapeMm); it.More(); it.Next()) {
                const TopoDS_Shape& child = it.Value();
                if (hasPeriodicFace(child)) {
                    BRepBuilderAPI_NurbsConvert cc(child, /*Copy=*/Standard_True);
                    rb.Add(rebuilt, cc.Shape().IsNull() ? child : cc.Shape());
                } else {
                    rb.Add(rebuilt, child);
                }
            }
            out = rebuilt;
        } else {
            BRepBuilderAPI_NurbsConvert conv(shapeMm, /*Copy=*/Standard_True);
            out = conv.Shape();
        }
        if (out.IsNull()) return shapeMm;
        // Conversion is exact, so both volume and validity must survive it. Anything else is OCC
        // damage, not a representation change, and it is refused loudly (below). ADAPTIVE
        // integration on BOTH sides: the default
        // VolumeProperties integrates an analytic face exactly and a B-spline face by fixed-order
        // quadrature, so comparing the two measures the INTEGRATION SCHEME, not the geometry
        // (every part came back ~0.17% "changed" and was reverted, which silently disabled this
        // whole pass). The Eps overload refines until the requested relative accuracy is met and
        // reports the accuracy it achieved; the difference is judged against that, plus what the
        // B-Rep can even express (every face is located only to within the shape tolerance, so
        // the volume is defined only to about tolerance x surface area).
        GProp_GProps g0, g1;
        const double eps = 1e-7;
        const double err0 = std::abs(BRepGProp::VolumeProperties(shapeMm, g0, eps));
        const double err1 = std::abs(BRepGProp::VolumeProperties(out, g1, eps));
        const double v0 = std::abs(g0.Mass());
        GProp_GProps sp;
        BRepGProp::SurfaceProperties(shapeMm, sp);
        // The model's own resolution is its FACES' tolerance, not Precision::Confusion(): in
        // this millimetre frame the scale carried the metre-frame 1e-7 to 1e-4 (100 nm, the
        // same physical tolerance), while Precision::Confusion() would read 0.1 nm and reject
        // an exact conversion whose quadrature differs in the 7th digit (12_boost, 2389.44 ->
        // 2389.44 mm3 "changed").
        double faceTol = Precision::Confusion();
        for (TopExp_Explorer fx(shapeMm, TopAbs_FACE); fx.More(); fx.Next())
            faceTol = std::max(faceTol, BRep_Tool::Tolerance(TopoDS::Face(fx.Current())));
        const double brepTol = faceTol * std::abs(sp.Mass());
        const double noise = std::max((err0 + err1) * v0, brepTol);
        // A failed conversion is not something to ship quietly: the FEM product would carry
        // periodic faces that the mesher then fails on, far from here. Say what happened.
        if (v0 > 0 && std::abs(g1.Mass() - v0) > noise) {
            std::ostringstream m;
            m.precision(10);
            m << "exportSTEP: the B-spline conversion of '" << name << "' changed its volume "
              << g0.Mass() << " -> " << g1.Mass() << " mm3, beyond the model's own " << noise
              << " mm3 (face tolerance " << faceTol << " mm, integrator errors " << err0 << " / "
              << err1 << "); an exact conversion cannot do that";
            throw std::runtime_error(m.str());
        }
        if (wasValid && !BRepCheck_Analyzer(out).IsValid())
            throw std::runtime_error("exportSTEP: the B-spline conversion of '" + name +
                                     "' produced an invalid solid");
        return out;
    } catch (const Standard_Failure& e) {
        throw std::runtime_error("exportSTEP: the B-spline conversion of '" + name + "' threw: " +
                                 std::string(e.GetMessageString()));
    }
}


// UNIT PARAMETRISATION OF THE WRITTEN B-SPLINES (ABT #1111, 2026-09-06).
//
// The model is built in metres and scaled x1000 here for the millimetre STEP. A transform moves
// a B-spline's poles and leaves its knots alone (OCCT: Geom_BSplineSurface::Transform, and
// ParametricTransformation is the identity for B-splines, so the pcurves stay too), so after the
// scale every B-spline surface that was parametrised by length in metres -- swept pipes, the
// converted facet strips and caps above -- spans ~1000 mm of geometry per parameter unit.
// Nothing in OCCT is truly scale-free in parameter space: BOPAlgo maps 3D tolerances to UV
// through exactly that ratio and reported sporadic face/face self-intersections on such solids.
// This pass rescales the knots of every B-spline surface so that one parameter unit spans about
// one millimetre in its steepest direction (an affine reparametrisation: EXACT, no point moves)
// and maps each pcurve on it by the same affinity, so every edge still lands on the same points.
// 3D curve parametrisation is left alone; the edge ranges do not change.
class UnitParametrisation : public BRepTools_Modification {
public:
    DEFINE_STANDARD_RTTI_INLINE(UnitParametrisation, BRepTools_Modification)

    Standard_Boolean NewSurface(const TopoDS_Face& F, Handle(Geom_Surface)& S, TopLoc_Location& L,
                                Standard_Real& Tol, Standard_Boolean& RevWires,
                                Standard_Boolean& RevFace) override {
        Handle(Geom_Surface) s0 = BRep_Tool::Surface(F, L);
        Handle(Geom_BSplineSurface) bs = Handle(Geom_BSplineSurface)::DownCast(s0);
        if (bs.IsNull()) return Standard_False;
        GeomAdaptor_Surface ad(bs);
        // Resolution(1 mm) is the parameter step that moves the surface by one millimetre at
        // its steepest, i.e. 1 / (mm per unit). Knots x that ratio makes the ratio one.
        const double ru = ad.UResolution(1.0), rv = ad.VResolution(1.0);
        const double ku = (ru > 0.0 && std::isfinite(ru)) ? 1.0 / ru : 1.0;
        const double kv = (rv > 0.0 && std::isfinite(rv)) ? 1.0 / rv : 1.0;
        if (std::abs(ku - 1.0) < 0.5 && std::abs(kv - 1.0) < 0.5) return Standard_False;
        // Only pcurves whose parametrisation survives an affinity exactly: B-spline/Bezier
        // (poles move, knots stay) and lines (converted first to a degree-1 B-spline, which
        // keeps the trim parameters as its knots). A conic pcurve would have to be
        // re-approximated, so a face carrying one is left alone and said so.
        for (TopExp_Explorer ex(F, TopAbs_EDGE); ex.More(); ex.Next()) {
            Standard_Real f, l;
            Handle(Geom2d_Curve) pc =
                BRep_Tool::CurveOnSurface(TopoDS::Edge(ex.Current()), F, f, l);
            if (pc.IsNull()) continue;
            Handle(Geom2d_Curve) basis = pc;
            while (Handle(Geom2d_TrimmedCurve) tc = Handle(Geom2d_TrimmedCurve)::DownCast(basis))
                basis = tc->BasisCurve();
            const bool ok = basis->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve)) ||
                            basis->IsKind(STANDARD_TYPE(Geom2d_BezierCurve)) ||
                            basis->IsKind(STANDARD_TYPE(Geom2d_Line));
            if (!ok) {
                std::cerr << "[reparam] face left at its metre parametrisation (" << ku << " x "
                          << kv << " mm/unit): pcurve of type "
                          << basis->DynamicType()->Name() << "\n";
                return Standard_False;
            }
        }
        Handle(Geom_BSplineSurface) nb = Handle(Geom_BSplineSurface)::DownCast(bs->Copy());
        TColStd_Array1OfReal uk(1, nb->NbUKnots()), vk(1, nb->NbVKnots());
        nb->UKnots(uk);
        nb->VKnots(vk);
        for (Standard_Integer i = uk.Lower(); i <= uk.Upper(); ++i) uk(i) *= ku;
        for (Standard_Integer i = vk.Lower(); i <= vk.Upper(); ++i) vk(i) *= kv;
        nb->SetUKnots(uk);
        nb->SetVKnots(vk);
        myFactors.Bind(F, std::make_pair(ku, kv));
        S = nb;
        Tol = BRep_Tool::Tolerance(F);
        RevWires = Standard_False;
        RevFace = Standard_False;
        return Standard_True;
    }
    Standard_Boolean NewCurve2d(const TopoDS_Edge& E, const TopoDS_Face& F, const TopoDS_Edge&,
                                const TopoDS_Face&, Handle(Geom2d_Curve)& C,
                                Standard_Real& Tol) override {
        const std::pair<double, double>* k = myFactors.Seek(F);
        if (!k) return Standard_False;
        Standard_Real f, l;
        Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(E, F, f, l);
        if (pc.IsNull()) return Standard_False;
        gp_GTrsf2d g;
        g.SetVectorialPart(gp_Mat2d(gp_XY(k->first, 0.0), gp_XY(0.0, k->second)));
        Handle(Geom2d_Curve) trimmed = new Geom2d_TrimmedCurve(pc, f, l);
        Handle(Geom2d_Curve) basis = pc;
        while (Handle(Geom2d_TrimmedCurve) tc = Handle(Geom2d_TrimmedCurve)::DownCast(basis))
            basis = tc->BasisCurve();
        if (basis->IsKind(STANDARD_TYPE(Geom2d_Line))) {
            // GeomLib::GTransform re-trims a line by projecting the transformed end points onto
            // the (unit-speed) transformed line, which changes the parameter range. A degree-1
            // B-spline with knots {f, l} is the same segment with the edge's own parameters.
            Handle(Geom2d_BSplineCurve) bs = Geom2dConvert::CurveToBSplineCurve(trimmed);
            if (bs.IsNull() || std::abs(bs->FirstParameter() - f) > 1e-9 * std::max(1.0, std::abs(f)) ||
                std::abs(bs->LastParameter() - l) > 1e-9 * std::max(1.0, std::abs(l)))
                throw std::runtime_error("exportSTEP: line pcurve did not convert to a B-spline"
                                         " on the edge's own parameters");
            trimmed = bs;
        }
        C = GeomLib::GTransform(trimmed, g);
        if (C.IsNull() || C->FirstParameter() > f + 1e-9 * std::max(1.0, std::abs(f)) ||
            C->LastParameter() < l - 1e-9 * std::max(1.0, std::abs(l)))
            throw std::runtime_error("exportSTEP: pcurve affinity changed the parameter range"
                                     " while reparametrising a B-spline face");
        Tol = BRep_Tool::Tolerance(E);
        return Standard_True;
    }
    Standard_Boolean NewCurve(const TopoDS_Edge&, Handle(Geom_Curve)&, TopLoc_Location&,
                              Standard_Real&) override { return Standard_False; }
    Standard_Boolean NewPoint(const TopoDS_Vertex&, gp_Pnt&, Standard_Real&) override {
        return Standard_False;
    }
    Standard_Boolean NewParameter(const TopoDS_Vertex&, const TopoDS_Edge&, Standard_Real&,
                                  Standard_Real&) override { return Standard_False; }
    GeomAbs_Shape Continuity(const TopoDS_Edge& E, const TopoDS_Face& F1, const TopoDS_Face& F2,
                             const TopoDS_Edge&, const TopoDS_Face&, const TopoDS_Face&) override {
        return BRep_Tool::Continuity(E, F1, F2);
    }

private:
    NCollection_DataMap<TopoDS_Shape, std::pair<double, double>, TopTools_ShapeMapHasher> myFactors;
};

TopoDS_Shape unitParametrised(const TopoDS_Shape& shapeMm, const std::string& name) {
    bool anyBSpline = false;
    for (TopExp_Explorer fx(shapeMm, TopAbs_FACE); fx.More() && !anyBSpline; fx.Next())
        anyBSpline = !Handle(Geom_BSplineSurface)::DownCast(
                          BRep_Tool::Surface(TopoDS::Face(fx.Current()))).IsNull();
    if (!anyBSpline) return shapeMm;
    const bool wasValid = BRepCheck_Analyzer(shapeMm).IsValid();
    Handle(UnitParametrisation) mod = new UnitParametrisation;
    BRepTools_Modifier modifier(shapeMm);
    modifier.Perform(mod);
    if (!modifier.IsDone())
        throw std::runtime_error("exportSTEP: B-spline reparametrisation of '" + name + "' failed");
    const TopoDS_Shape out = modifier.ModifiedShape(shapeMm);
    BRepCheck_Analyzer ana(out);
    if (wasValid && !ana.IsValid()) {
        // Name the statuses: a bare "invalid" cannot be acted on.
        std::ostringstream why;
        int listed = 0;
        const TopAbs_ShapeEnum kinds[] = {TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE, TopAbs_WIRE,
                                          TopAbs_EDGE,  TopAbs_VERTEX};
        for (TopAbs_ShapeEnum kind : kinds) {
            int idx = 0;
            for (TopExp_Explorer e(out, kind); e.More() && listed < 12; e.Next(), ++idx) {
                const Handle(BRepCheck_Result)& res = ana.Result(e.Current());
                if (res.IsNull()) continue;
                bool bad = false;
                std::ostringstream st;
                for (BRepCheck_ListIteratorOfListOfStatus it(res->Status()); it.More(); it.Next()) {
                    if (it.Value() == BRepCheck_NoError) continue;
                    bad = true;
                    BRepCheck::Print(it.Value(), st);
                }
                if (!bad) continue;
                std::string sts = st.str();
                for (char& c : sts) if (c == '\n') c = ' ';
                why << " [" << TopAbs::ShapeTypeToString(kind) << " " << idx << ": " << sts << "]";
                ++listed;
            }
        }
        if (const char* dump = std::getenv("MVB_REPARAM_DUMP")) {
            BRepTools::Write(shapeMm, (std::string(dump) + ".before.brep").c_str());
            BRepTools::Write(out, (std::string(dump) + ".after.brep").c_str());
        }
        throw std::runtime_error("exportSTEP: B-spline reparametrisation left '" + name +
                                 "' invalid (an exact affine change):" + why.str());
    }
    return out;
}

} // namespace

bool exportSTEP(const std::vector<NamedShape>& shapes, const std::string& filepath) {
    return exportSTEP(shapes, filepath, StepExportOptions{});
}

bool exportSTEP(const std::vector<NamedShape>& shapes,
                const std::string& filepath,
                const StepExportOptions& options) {
    if (shapes.empty()) return false;

    Handle(TDocStd_Document) doc =
        new TDocStd_Document(TCollection_ExtendedString("BinXCAF"));
    Handle(XCAFDoc_ShapeTool) shapeTool =
        XCAFDoc_DocumentTool::ShapeTool(doc->Main());

    // Model coordinates are in SI (metres) — the whole MVB++/MKF/MAS convention
    // (see exportSTLToBytes below). OCCT's length unit is millimetres and it
    // faithfully labels whatever numbers it is given: writing metre-valued
    // coordinates with STEPCAFControl_Writer's default emits
    // SI_UNIT(.MILLI.,.METRE.) over numbers that are actually metres, so every
    // STEP opens 1000x too small (ABT #317). Setting write.step.unit=M does NOT
    // fix it — OCCT then divides the coordinates by 1000 to match the metre
    // header, preserving the same wrong physical size (verified empirically).
    // The only correct fix is to convert the geometry metres->millimetres here,
    // once, so the emitted numbers match the default mm header. This is the
    // single source of truth for the conversion: callers pass metres and MUST
    // NOT pre-scale.
    gp_Trsf metresToMillimetres;
    metresToMillimetres.SetScale(gp_Pnt(0, 0, 0), 1000.0);

    // MVB_EXPORT_VALIDITY_DIAG: report per-solid B-Rep validity IMMEDIATELY BEFORE and
    // IMMEDIATELY AFTER the metres->millimetres scale. The assembler's own exit gate proves
    // every solid leaves construction valid, yet 06_llc reads back from STEP with one invalid
    // solid; the scale is the only geometric operation in between, and ABT #860 already
    // records that BRepBuilderAPI_Transform carries tolerances through the same factor.
    const bool validityDiag = std::getenv("MVB_EXPORT_VALIDITY_DIAG") != nullptr;

    for (const auto& ns : shapes) {
        if (ns.shape.IsNull()) continue;
        if (validityDiag) {
            int idx = 0, badBefore = 0;
            for (TopExp_Explorer e(ns.shape, TopAbs_SOLID); e.More(); e.Next(), ++idx) {
                if (!BRepCheck_Analyzer(e.Current()).IsValid()) {
                    ++badBefore;
                    std::fprintf(stderr,
                                 "[export-diag] PRE-SCALE invalid: '%s' solid %d\n",
                                 ns.name.c_str(), idx);
                }
            }
            if (badBefore == 0)
                std::fprintf(stderr, "[export-diag] PRE-SCALE all %d solid(s) valid: '%s'\n",
                             idx, ns.name.c_str());
        }
        TopoDS_Shape shapeMm =
            BRepBuilderAPI_Transform(ns.shape, metresToMillimetres).Shape();
        // Deliberately NOT applied to the drawing product: exact analytic surfaces are smaller
        // and render better, and a viewer never meets gmsh. MVB_NO_NURBS=1 disables it for
        // comparison.
        if (options.nurbsPeriodicSolids && !std::getenv("MVB_NO_NURBS"))
            shapeMm = nurbsConvertPeriodicSolidsMm(shapeMm, ns.name);
        // Every product: the file's B-splines are parametrised in the millimetres it is
        // written in (see UnitParametrisation). MVB_NO_REPARAM=1 disables it for comparison.
        if (!std::getenv("MVB_NO_REPARAM")) shapeMm = unitParametrised(shapeMm, ns.name);
        if (validityDiag) {
            int idx = 0, badAfter = 0;
            for (TopExp_Explorer e(shapeMm, TopAbs_SOLID); e.More(); e.Next(), ++idx) {
                if (!BRepCheck_Analyzer(e.Current()).IsValid()) {
                    ++badAfter;
                    std::fprintf(stderr,
                                 "[export-diag] POST-SCALE invalid: '%s' solid %d\n",
                                 ns.name.c_str(), idx);
                }
            }
            if (badAfter == 0)
                std::fprintf(stderr, "[export-diag] POST-SCALE all %d solid(s) valid: '%s'\n",
                             idx, ns.name.c_str());
        }
        // ABT #685: with per-solid names supplied, the shape goes in as an ASSEMBLY so every solid
        // becomes its own named product. Without them a multi-solid conductor arrives as one
        // unnamed-parts product and the viewer invents the numbering ("Primary parallel 001"),
        // which exists nowhere in the file — so a part cannot be named unambiguously when
        // discussing a build. Guarded on the count matching the solids actually present: a
        // mismatch means the caller's list is stale, and a wrong name is worse than none.
        std::size_t solidCount = 0;
        for (TopExp_Explorer e(shapeMm, TopAbs_SOLID); e.More(); e.Next()) ++solidCount;
        const bool nameParts = !ns.partNames.empty() && ns.partNames.size() == solidCount &&
                               solidCount > 1;
        // AddShape(isAssembly=false): register the shape as a free top-level
        // component so STEPCAFControl_Writer emits it as a discrete product
        // (with our name attached) rather than burying it in a compound.
        TDF_Label lab = shapeTool->AddShape(shapeMm, nameParts ? Standard_True : Standard_False);
        if (!ns.name.empty()) {
            TDataStd_Name::Set(lab, TCollection_ExtendedString(ns.name.c_str()));
        }
        if (std::getenv("MVB_STEP_NAME_DIAG")) {
            std::fprintf(stderr, "[step] '%s': solids=%zu partNames=%zu nameParts=%d\n",
                         ns.name.c_str(), solidCount, ns.partNames.size(), int(nameParts));
        }
        if (nameParts) {
            TDF_LabelSequence components;
            shapeTool->GetComponents(lab, components);
            if (std::getenv("MVB_STEP_NAME_DIAG")) {
                std::fprintf(stderr, "[step]   components=%d\n", components.Length());
            }
            for (Standard_Integer k = 1;
                 k <= components.Length() && std::size_t(k) <= ns.partNames.size(); ++k) {
                const std::string& partName = ns.partNames[k - 1];
                if (partName.empty()) continue;
                TDF_Label component = components.Value(k);
                TDataStd_Name::Set(component, TCollection_ExtendedString(partName.c_str()));
                // The name belongs on the referred PRODUCT too: STEP writes the component as an
                // occurrence of a product, and most viewers show the product's name.
                TDF_Label referred;
                if (shapeTool->GetReferredShape(component, referred)) {
                    TDataStd_Name::Set(referred, TCollection_ExtendedString(partName.c_str()));
                }
            }
        }
    }

    // MVB_STEP_WRITER_SWEEP: the geometry is provably valid in memory at this point (see the
    // validity diag above), so any invalid solid read back from the file is introduced by the
    // WRITE/READ round-trip. Writing and re-reading costs seconds against a ~25 min geometry
    // build, so try the candidate writer settings here, in one run, and report which of them
    // survives a round-trip. write.precision.mode: 0=average tolerance of the shape (OCC's
    // default, and the suspect -- an average UNDER-states the tolerance of the loosest face,
    // so on re-read that face violates its own written tolerance), 2=maximum tolerance,
    // 1=explicit value. write.surfacecurve.mode 0 drops pcurves so the reader recomputes them.
    if (std::getenv("MVB_STEP_WRITER_SWEEP")) {
        struct Cfg { const char* name; int precMode; double precVal; int curveMode; };
        const Cfg cfgs[] = {
            {"default(prec=0,pcurves=1)", 0, 0.0, 1},
            {"prec=2(max)", 2, 0.0, 1},
            {"prec=1(val=1e-5)", 1, 1e-5, 1},
            {"prec=2,nopcurves", 2, 0.0, 0},
        };
        for (const Cfg& c : cfgs) {
            Interface_Static::SetIVal("write.precision.mode", c.precMode);
            if (c.precMode == 1) Interface_Static::SetRVal("write.precision.val", c.precVal);
            Interface_Static::SetIVal("write.surfacecurve.mode", c.curveMode);
            const std::string tmp = filepath + ".sweep.step";
            STEPCAFControl_Writer w;
            if (!w.Transfer(doc) || w.Write(tmp.c_str()) != IFSelect_RetDone) {
                std::fprintf(stderr, "[writer-sweep] %-28s WRITE FAILED\n", c.name);
                continue;
            }
            int bad = 0, total = 0;
            for (const auto& rs : importSTEP(tmp)) {
                for (TopExp_Explorer e(rs.shape, TopAbs_SOLID); e.More(); e.Next()) {
                    ++total;
                    if (!BRepCheck_Analyzer(e.Current()).IsValid()) ++bad;
                }
            }
            std::fprintf(stderr, "[writer-sweep] %-28s invalid=%d of %d solids\n",
                         c.name, bad, total);
            std::error_code ec;
            std::filesystem::remove(tmp, ec);
        }
        // Leave the writer on OCC's defaults; the chosen setting is applied below.
        Interface_Static::SetIVal("write.precision.mode", 0);
        Interface_Static::SetIVal("write.surfacecurve.mode", 1);
    }

    STEPCAFControl_Writer writer;
    if (!writer.Transfer(doc)) {
        std::cerr << "mvb::exportSTEP: STEPCAFControl_Writer::Transfer failed\n";
        return false;
    }
    if (writer.Write(filepath.c_str()) != IFSelect_RetDone) {
        std::cerr << "mvb::exportSTEP: write failed to " << filepath << "\n";
        return false;
    }
    return true;
}

bool exportSTEP(const std::vector<TopoDS_Shape>& shapes,
                const std::vector<std::string>& names,
                const std::string& filepath) {
    std::vector<NamedShape> ns;
    ns.reserve(shapes.size());
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        const std::string n =
            (i < names.size()) ? names[i] : std::string{};
        ns.emplace_back(shapes[i], n);
    }
    return exportSTEP(ns, filepath);
}

std::vector<NamedShape> importSTEP(const std::string& filepath) {
    std::vector<NamedShape> out;

    Handle(TDocStd_Document) doc =
        new TDocStd_Document(TCollection_ExtendedString("BinXCAF"));

    STEPCAFControl_Reader reader;
    if (reader.ReadFile(filepath.c_str()) != IFSelect_RetDone) return out;
    if (!reader.Transfer(doc)) return out;

    Handle(XCAFDoc_ShapeTool) shapeTool =
        XCAFDoc_DocumentTool::ShapeTool(doc->Main());

    TDF_LabelSequence freeShapes;
    shapeTool->GetFreeShapes(freeShapes);
    for (Standard_Integer i = 1; i <= freeShapes.Length(); ++i) {
        const TDF_Label lab = freeShapes.Value(i);
        TopoDS_Shape s = shapeTool->GetShape(lab);
        if (s.IsNull()) continue;
        out.emplace_back(s, label_name(lab));
    }
    return out;
}

bool exportSTL(const TopoDS_Shape& compound, const std::string& filepath) {
    if (compound.IsNull()) return false;
    // Own the metres->millimetres conversion here, like exportSTEP and
    // exportSTLToBytes (ABT #317), so every exporter is uniform: callers hand over
    // native SI-metre geometry and MUST pass scale=1.0 (NOT the old 1000 — passing
    // 1000 now double-scales, 1e6x too big). The 0.001 deflection was always in the
    // caller's mm space, so it is unchanged.
    gp_Trsf metresToMillimetres;
    metresToMillimetres.SetScale(gp_Pnt(0, 0, 0), 1000.0);
    const TopoDS_Shape mm = BRepBuilderAPI_Transform(compound, metresToMillimetres).Shape();
    BRepMesh_IncrementalMesh mesh(mm, 0.001);
    StlAPI_Writer writer;
    writer.Write(mm, filepath.c_str());
    return true;
}

std::string exportSTLToBytes(const std::vector<TopoDS_Shape>& shapes,
                              double toleranceMm,
                              double angularTolerance,
                              bool binary)
{
    if (shapes.empty()) return {};

    // Model coordinates are in SI (metres); scale metres -> millimetres so the STL
    // opens at correct physical size in the CAD/slicer tools that consume it (they
    // assume mm — the frontend even passes its mesh tolerance as `tolMm`). Same
    // convention as exportSTEP (ABT #317): a raw-metre STL opened 1000x too small.
    gp_Trsf metresToMillimetres;
    metresToMillimetres.SetScale(gp_Pnt(0, 0, 0), 1000.0);

    // Combine into a single compound — StlAPI_Writer walks sub-solids.
    TopoDS_Compound compound;
    TopoDS_Builder b;
    b.MakeCompound(compound);
    bool any = false;
    for (const auto& s : shapes) {
        if (s.IsNull()) continue;
        b.Add(compound, BRepBuilderAPI_Transform(s, metresToMillimetres).Shape());
        any = true;
    }
    if (!any) return {};

    // Geometry is now in millimetres (scaled above). Compute a scale-appropriate
    // absolute linear deflection from the overall bounding-box diagonal. Using
    // relative deflection on tiny per-face boxes (e.g. individual wire turns) creates
    // sub-micron meshes and OOMs in WASM. Absolute deflection ~0.5 % of overall size,
    // floored at 1.0 mm (was 0.001 m before the mm scaling — physically identical),
    // keeps mesh counts bounded while preserving recognisable geometry.
    Bnd_Box bbox;
    BRepBndLib::Add(compound, bbox);
    Standard_Real xMin, yMin, zMin, xMax, yMax, zMax;
    bbox.Get(xMin, yMin, zMin, xMax, yMax, zMax);
    const Standard_Real dx = xMax - xMin;
    const Standard_Real dy = yMax - yMin;
    const Standard_Real dz = zMax - zMin;
    const Standard_Real diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
    const Standard_Real linDeflection = std::max(diagonal * 0.005, 1.0);
    BRepMesh_IncrementalMesh mesh(compound, linDeflection, Standard_False,
                                   angularTolerance, Standard_False);
    mesh.Perform();

    // StlAPI_Writer has no in-memory API; write to a temp file and slurp.
    std::filesystem::path tmpDir = std::filesystem::temp_directory_path();
    // Use random suffix to avoid race conditions with concurrent calls.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dist;
    std::string tmpName = "mvbpp_stl_" + std::to_string(dist(gen)) + ".stl";
    std::filesystem::path tmpPath = tmpDir / tmpName;
    StlAPI_Writer writer;
    writer.ASCIIMode() = !binary;
    if (!writer.Write(compound, tmpPath.string().c_str())) return {};

    std::ifstream f(tmpPath, std::ios::binary);
    if (!f) { std::filesystem::remove(tmpPath); return {}; }
    std::string data((std::istreambuf_iterator<char>(f)),
                      std::istreambuf_iterator<char>());
    std::filesystem::remove(tmpPath);
    return data;
}

} // namespace mvb
