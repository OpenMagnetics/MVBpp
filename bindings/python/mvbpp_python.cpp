// MVB++ Python bindings — unified 7-function API.
//
//     drawCore        - 3D/2D MAS::MagneticCore
//     drawCorePiece   - 3D/2D MAS::CoreShape (single piece via OpenMagnetics::CorePiece)
//     drawBobbin      - 3D/2D MAS::Bobbin (must carry processedDescription)
//     drawTurns       - 3D/2D list of self-sufficient MAS::Turn
//     drawWinding     - 3D/2D MAS::Coil filtered by winding name
//     drawMagnetic    - 3D/2D MAS::Magnetic (full assembly)
//     drawView        - SVG pictorial cross-section, optionally dimensioned
//
// All draw* (3D/2D) functions share the same call shape:
//     fn(json_str, output_path=None, *, mode="3D", plane="XY", offset=0.0,
//        format="step", scale=1.0, polygonSegments=32)
// outputPath=None  → return bytes
// outputPath=<str> → write file, return path
// mode="2D"        → cut at plane+offset, export planar Faces (STEP only;
//                     STL throws — meshing in 2D belongs in Gmsh/ElmerFEM).

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mvb/MagneticBuilder.h"
#include "mvb/SectionBuilder.h"
#include "mvb/SectionDrawing.h"
#include "mvb/StepExporter.h"
#include "mvb/Symmetry.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"

#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;
using nlohmann::json;

namespace {

// --------------------------------------------------------------------------
// Common helpers
// --------------------------------------------------------------------------

std::string upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

void apply_scale(std::vector<mvb::NamedShape>& named, double scale) {
    if (scale == 1.0) return;
    gp_Trsf trsf;
    trsf.SetScale(gp_Pnt(0, 0, 0), scale);
    for (auto& ns : named) {
        if (ns.shape.IsNull()) continue;
        ns.shape = BRepBuilderAPI_Transform(ns.shape, trsf).Shape();
    }
}

std::string slurp_file(const std::filesystem::path& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("mvbpp: failed to read " + p.string());
    return std::string{std::istreambuf_iterator<char>(f),
                       std::istreambuf_iterator<char>()};
}

void write_bytes(const std::filesystem::path& path, const std::string& data) {
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("mvbpp: cannot write to " + path.string());
    f.write(data.data(), static_cast<std::streamsize>(data.size()));
}

// Encode a vector of NamedShape to STEP/STL bytes.
std::string encode_to_bytes(const std::vector<mvb::NamedShape>& named,
                             const std::string& format) {
    std::string fmt = format;
    std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (fmt == "stl") {
        std::vector<TopoDS_Shape> shapes;
        shapes.reserve(named.size());
        for (const auto& ns : named) shapes.push_back(ns.shape);
        std::string data = mvb::exportSTLToBytes(shapes);
        if (data.empty()) throw std::runtime_error("mvbpp: STL export produced empty output");
        return data;
    }
    if (fmt == "step") {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<uint64_t> dist;
        auto tmp = std::filesystem::temp_directory_path() /
                   ("mvbpp_" + std::to_string(dist(gen)) + ".step");
        struct Cleanup {
            std::filesystem::path p;
            ~Cleanup() { std::error_code ec; std::filesystem::remove(p, ec); }
        } cleanup{tmp};
        if (!mvb::exportSTEP(named, tmp.string())) {
            throw std::runtime_error("mvbpp: exportSTEP failed");
        }
        return slurp_file(tmp);
    }
    throw std::runtime_error("mvbpp: unsupported format '" + format +
                             "' (expected 'step' or 'stl')");
}

// Deliver `named` as either bytes (output_path=None) or a written file path.
py::object deliver(std::vector<mvb::NamedShape> named,
                    const py::object& output_path,
                    const std::string& mode,
                    const std::string& plane,
                    double offset,
                    const std::string& format,
                    double scale,
                    const std::string& symmetry,
                    const std::string& side) {
    std::string m = upper(mode);
    if (m != "2D" && m != "3D") {
        throw std::runtime_error("mvbpp: mode must be '3D' or '2D' (got '" + mode + "')");
    }

    apply_scale(named, scale);

    // Symmetry cuts run on 3-D solids first (they're meaningless on 2-D
    // faces, and pre-cutting the solid produces cleaner 2-D sections).
    int nSym = mvb::parse_symmetry_spec(symmetry);
    if (nSym > 0) {
        named = mvb::apply_symmetry(named, nSym);
        if (named.empty()) {
            throw std::runtime_error(
                "mvbpp: symmetry='" + symmetry +
                "' removed all geometry (no valid planes detected)");
        }
    }

    std::string fmt = format;
    std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (m == "2D") {
        if (fmt == "stl") {
            throw std::runtime_error(
                "mvbpp: mode='2D' requires format='step' (STL is not valid for "
                "planar faces; use STEP for ElmerFEM/Gmsh 2D meshing).");
        }
        named = mvb::SectionBuilder::cut2DFaces(named,
                                                  mvb::parseSectionPlane(plane),
                                                  offset);
        if (named.empty()) {
            throw std::runtime_error(
                "mvbpp: 2D section produced no faces (plane=" + plane +
                ", offset=" + std::to_string(offset) + ")");
        }
    }

    // Side filter runs last so it works on 2-D faces and 3-D solids alike.
    auto axisSign = mvb::parse_side_spec(side);
    named = mvb::filter_by_side(named, axisSign);
    if (named.empty()) {
        throw std::runtime_error(
            "mvbpp: side='" + side + "' filtered out all geometry");
    }

    std::string data = encode_to_bytes(named, fmt);

    if (output_path.is_none()) {
        return py::bytes(data);
    }
    auto path_str = output_path.cast<std::string>();
    write_bytes(path_str, data);
    return py::cast(path_str);
}

// --------------------------------------------------------------------------
// Per-function builders — each returns a vector<NamedShape>
// --------------------------------------------------------------------------

std::vector<mvb::NamedShape> build_core(const std::string& json_str, int polygonSegments) {
    auto j = json::parse(json_str);
    auto core = j.get<MAS::MagneticCore>();
    if (!core.get_geometrical_description()) {
        throw std::runtime_error(
            "drawCore: MagneticCore.geometricalDescription is required. "
            "Run MKF magnetic_autocomplete on the parent Magnetic, or use "
            "drawMagnetic which performs autocomplete itself.");
    }
    mvb::MagneticBuilder b;
    return b.buildCoreNamed(core, polygonSegments);
}

std::vector<mvb::NamedShape> build_core_piece(const std::string& json_str, int polygonSegments) {
    auto j = json::parse(json_str);
    auto shape = j.get<MAS::CoreShape>();
    mvb::MagneticBuilder b;
    return {b.buildCorePieceNamed(shape, polygonSegments)};
}

std::vector<mvb::NamedShape> build_bobbin(const std::string& json_str, int polygonSegments) {
    auto j = json::parse(json_str);
    auto bobbin = j.get<MAS::Bobbin>();
    mvb::MagneticBuilder b;
    // axisIsY=true matches the concentric-core convention. Toroidal bobbins
    // are rare in standalone form; users should call drawMagnetic instead.
    return {b.buildBobbinNamedFromBobbin(bobbin, /*axisIsY=*/true, polygonSegments)};
}

std::vector<mvb::NamedShape> build_turns(const std::string& json_str, int polygonSegments,
                                         bool paintCoating) {
    auto j = json::parse(json_str);
    if (!j.is_array()) {
        throw std::runtime_error("drawTurns: input must be a JSON array of MAS::Turn objects");
    }
    auto turns = j.get<std::vector<MAS::Turn>>();
    mvb::MagneticBuilder b;
    // paintCoating=false throws here — a bare Turn carries no copper/strand
    // data; use drawMagnetic with a full Magnetic JSON for conductor geometry.
    return b.buildTurnsNamedFromTurns(turns, polygonSegments, paintCoating);
}

std::vector<mvb::NamedShape> build_winding(const std::string& json_str,
                                            const std::string& windingName,
                                            int polygonSegments,
                                            bool paintCoating) {
    auto j = json::parse(json_str);
    auto coil = j.get<MAS::Coil>();
    auto turnsOpt = coil.get_turns_description();
    if (!turnsOpt || turnsOpt->empty()) {
        throw std::runtime_error("drawWinding: Coil.turnsDescription is empty or missing");
    }
    std::string target = windingName;
    std::transform(target.begin(), target.end(), target.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    std::vector<MAS::Turn> filtered;
    for (const auto& t : *turnsOpt) {
        std::string w = t.get_winding();
        std::transform(w.begin(), w.end(), w.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (w == target) filtered.push_back(t);
    }
    if (filtered.empty()) {
        throw std::runtime_error("drawWinding: no turns match windingName='" + windingName + "'");
    }
    mvb::MagneticBuilder b;
    return b.buildTurnsNamedFromTurns(filtered, polygonSegments, paintCoating);
}

std::vector<mvb::NamedShape> build_magnetic(const std::string& json_str, int polygonSegments,
                                            bool paintCoating) {
    auto j = json::parse(json_str);
    auto magnetic = j.get<MAS::Magnetic>();
    mvb::MagneticBuilder b;
    // Use polygonSegments for both wire and core for simplicity; the user
    // can craft separate paths later if they need finer-grained control.
    // paintCoating: true → OUTER (insulation) diameter [default], false →
    // CONDUCTING (copper) diameter for FEM winding-loss meshing (LITZ → bare
    // bundle as a solid via MKF).
    return b.buildAllNamed(magnetic, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                            polygonSegments, polygonSegments, paintCoating);
}

// --------------------------------------------------------------------------
// drawView — SVG pictorial output (autodetect MagneticCore vs Magnetic)
// --------------------------------------------------------------------------

py::object draw_view_impl(const std::string& json_str,
                            const py::object& output_path,
                            bool dimensions,
                            const std::string& plane,
                            double offset,
                            double widthPx,
                            const std::string& format) {
    std::string fmt = format;
    std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (fmt != "svg") {
        throw std::runtime_error("drawView: only format='svg' is supported (got '" + format + "')");
    }

    auto j = json::parse(json_str);

    OpenMagnetics::Magnetic magnetic;
    bool is_magnetic = j.contains("coil") || j.contains("magnetizing_inductance");
    if (is_magnetic) {
        magnetic = mvb::magnetic_autocomplete_safe(j);
    } else if (j.contains("functionalDescription") || j.contains("functional_description")) {
        // Looks like a MagneticCore — wrap in a minimal Magnetic so we can
        // call SectionDrawing (it operates on Magnetic).
        json wrap;
        wrap["core"] = j;
        // Provide a minimal coil so MAS::Magnetic deserialisation does not
        // throw. We do not need turns/wires for a core-only view.
        wrap["coil"] = json{
            {"bobbin", "Dummy"},
            {"functionalDescription", json::array()}
        };
        magnetic = mvb::magnetic_autocomplete_safe(wrap);
    } else {
        throw std::runtime_error(
            "drawView: input JSON does not look like a MAS::Magnetic or "
            "MAS::MagneticCore (no functionalDescription / coil keys)");
    }

    std::string svg = mvb::SectionDrawing::drawView(magnetic, plane, offset,
                                                     dimensions, widthPx);

    if (output_path.is_none()) {
        return py::cast(svg);  // return as Python str
    }
    auto path_str = output_path.cast<std::string>();
    write_bytes(path_str, svg);
    return py::cast(path_str);
}

} // namespace

PYBIND11_MODULE(mvbpp, m) {
    m.doc() = R"doc(
MVB++ — magnetics 3-D geometry builder (Python bindings).

The public API consists of seven functions. Each draw* function returns
either bytes (when outputPath is None) or the file path that was written.
Inputs are MAS-1.0 JSON strings.

3D/2D draw functions
--------------------
    mvbpp.drawCore       (core_json,        outputPath=None, *, mode="3D", plane="XY", offset=0.0, format="step", scale=1.0, polygonSegments=32)
    mvbpp.drawCorePiece  (core_shape_json,  outputPath=None, *, mode="3D", plane="XY", offset=0.0, format="step", scale=1.0, polygonSegments=32)
    mvbpp.drawBobbin     (bobbin_json,      outputPath=None, *, mode="3D", plane="XY", offset=0.0, format="step", scale=1.0, polygonSegments=32)
    mvbpp.drawTurns      (turns_json,       outputPath=None, *, mode="3D", plane="XY", offset=0.0, format="step", scale=1.0, polygonSegments=32, paintCoating=True)
    mvbpp.drawWinding    (coil_json, windingName, outputPath=None, *, mode="3D", plane="XY", offset=0.0, format="step", scale=1.0, polygonSegments=32, paintCoating=True)
    mvbpp.drawMagnetic   (magnetic_json,    outputPath=None, *, mode="3D", plane="XY", offset=0.0, format="step", scale=1.0, polygonSegments=32, paintCoating=True)

Pictorial-view function
-----------------------
    mvbpp.drawView       (obj_json,         outputPath=None, *, dimensions=True, plane="XY", offset=0.0, widthPx=800, format="svg")

Notes
-----
* mode="2D" cuts at the requested plane (XY/XZ/YZ) shifted by `offset` (m)
  and exports a STEP compound of planar TopoDS_Faces (suitable for ElmerFEM
  surface meshing).  STL is not allowed in 2D mode.
* polygonSegments=0 attempts exact NURBS surfaces; some shapes (RM, PQ) are
  known to fail OCCT booleans in this mode.  When in doubt use the default.
* drawView with dimensions=True currently supports plane in {"XY","XZ"} and
  offset=0.0; other combinations throw.  dimensions=False renders a plain
  outline at any plane/offset.
* paintCoating (drawMagnetic/drawTurns/drawWinding): True → turns drawn at the
  OUTER (insulation) diameter [default, for visualisation]; False → turns drawn
  at the CONDUCTING (copper) diameter for FEM winding-loss meshing (LITZ → bare
  bundle treated as one solid conductor via MKF).  False requires full wire
  data: only drawMagnetic (full Magnetic JSON) supports it; drawTurns/drawWinding
  operate on bare turns and raise for paintCoating=False.
)doc";

    auto def_draw = [&](const char* name,
                         std::vector<mvb::NamedShape> (*build)(const std::string&, int)) {
        m.def(name,
              [build](const std::string& json_str,
                      py::object outputPath,
                      const std::string& mode,
                      const std::string& plane,
                      double offset,
                      const std::string& format,
                      double scale,
                      int polygonSegments,
                      const std::string& symmetry,
                      const std::string& side) {
                  auto named = build(json_str, polygonSegments);
                  return deliver(std::move(named), outputPath, mode, plane,
                                  offset, format, scale, symmetry, side);
              },
              py::arg("json_str"),
              py::arg("outputPath") = py::none(),
              py::kw_only(),
              py::arg("mode") = std::string("3D"),
              py::arg("plane") = std::string("XY"),
              py::arg("offset") = 0.0,
              py::arg("format") = std::string("step"),
              py::arg("scale") = 1.0,
              py::arg("polygonSegments") = mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
              py::arg("symmetry") = std::string("none"),
              py::arg("side") = std::string("+X+Y+Z"));
    };

    // Paint-aware variant for the turn-drawing functions: adds a trailing
    // keyword-only `paintCoating` (true → OUTER/insulation [default], false →
    // CONDUCTING/copper for FEM winding-loss meshing). pybind11 honours the
    // default, so existing callers are unaffected.
    auto def_draw_paint = [&](const char* name,
                              std::vector<mvb::NamedShape> (*build)(const std::string&, int, bool)) {
        m.def(name,
              [build](const std::string& json_str,
                      py::object outputPath,
                      const std::string& mode,
                      const std::string& plane,
                      double offset,
                      const std::string& format,
                      double scale,
                      int polygonSegments,
                      const std::string& symmetry,
                      const std::string& side,
                      bool paintCoating) {
                  auto named = build(json_str, polygonSegments, paintCoating);
                  return deliver(std::move(named), outputPath, mode, plane,
                                  offset, format, scale, symmetry, side);
              },
              py::arg("json_str"),
              py::arg("outputPath") = py::none(),
              py::kw_only(),
              py::arg("mode") = std::string("3D"),
              py::arg("plane") = std::string("XY"),
              py::arg("offset") = 0.0,
              py::arg("format") = std::string("step"),
              py::arg("scale") = 1.0,
              py::arg("polygonSegments") = mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
              py::arg("symmetry") = std::string("none"),
              py::arg("side") = std::string("+X+Y+Z"),
              py::arg("paintCoating") = true);
    };

    def_draw("drawCore",       &build_core);
    def_draw("drawCorePiece",  &build_core_piece);
    def_draw("drawBobbin",     &build_bobbin);
    def_draw_paint("drawTurns",    &build_turns);
    def_draw_paint("drawMagnetic", &build_magnetic);

    // drawWinding takes an extra positional `windingName`.
    m.def("drawWinding",
          [](const std::string& coil_json,
              const std::string& windingName,
              py::object outputPath,
              const std::string& mode,
              const std::string& plane,
              double offset,
              const std::string& format,
              double scale,
              int polygonSegments,
              const std::string& symmetry,
              const std::string& side,
              bool paintCoating) {
              auto named = build_winding(coil_json, windingName, polygonSegments, paintCoating);
              return deliver(std::move(named), outputPath, mode, plane,
                              offset, format, scale, symmetry, side);
          },
          py::arg("coil_json"),
          py::arg("windingName"),
          py::arg("outputPath") = py::none(),
          py::kw_only(),
          py::arg("mode") = std::string("3D"),
          py::arg("plane") = std::string("XY"),
          py::arg("offset") = 0.0,
          py::arg("format") = std::string("step"),
          py::arg("scale") = 1.0,
          py::arg("polygonSegments") = mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
          py::arg("symmetry") = std::string("none"),
          py::arg("side") = std::string("+X+Y+Z"),
          py::arg("paintCoating") = true);

    // Test-only helper: enrich a raw MAS::Magnetic JSON via MKF autocomplete
    // so test code can extract pre-enriched substructures for drawCore /
    // drawBobbin / drawTurns / drawWinding (which by contract reject raw
    // input).  Underscore-prefixed to mark it private / unstable.
    m.def("_enrichMagnetic",
          [](const std::string& json_str) {
              auto j = nlohmann::json::parse(json_str);
              auto enriched = mvb::magnetic_autocomplete_safe(j);
              nlohmann::json out;
              OpenMagnetics::to_json(out, enriched);
              return out.dump();
          },
          py::arg("magnetic_json"),
          "Test helper: run MKF magnetic_autocomplete and return enriched MAS JSON.");

    m.def("drawView", &draw_view_impl,
          py::arg("json_str"),
          py::arg("outputPath") = py::none(),
          py::kw_only(),
          py::arg("dimensions") = true,
          py::arg("plane") = std::string("XY"),
          py::arg("offset") = 0.0,
          py::arg("widthPx") = 800.0,
          py::arg("format") = std::string("svg"));
}
