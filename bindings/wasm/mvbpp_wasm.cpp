#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <emscripten/emscripten.h>

#include "mvb/MagneticBuilder.h"
#include "mvb/Symmetry.h"
#include "mvb/StepExporter.h"
#include "mvb/SectionDrawing.h"
#include "mvb/SectionBuilder.h"
#include "mvb/ProjectionDrawing.h"
#include "mvb/SpacerBuilder.h"
#include "mvb/FR4Builder.h"
#include "mvb/BobbinBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Coil.h"
#include <nlohmann/json.hpp>

#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace emscripten;
using json = nlohmann::json;

// ─────────────────────────────────────────────────────────────────────────────
// JS-throw helper + guard template.
// ─────────────────────────────────────────────────────────────────────────────
[[noreturn]] inline void throw_js_error(const std::string& msg) {
    EM_ASM({ throw new Error(UTF8ToString($0)); }, msg.c_str());
    __builtin_unreachable();
}

template <auto Fn> struct guard;

template <typename R, typename... Args, R(*Fn)(Args...)>
struct guard<Fn> {
    static R call(Args... args) {
        try {
            return Fn(args...);
        } catch (const std::exception& e) {
            throw_js_error(std::string("[mvbpp] ") + e.what());
        } catch (...) {
            throw_js_error("[mvbpp] unknown C++ exception");
        }
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Common helpers (shared logic with Python bindings)
// ─────────────────────────────────────────────────────────────────────────────
namespace {

val makeUint8Array(const std::string& data) {
    val memoryView = val(typed_memory_view(data.size(),
                                            reinterpret_cast<const uint8_t*>(data.data())));
    val out = val::global("Uint8Array").new_(data.size());
    out.call<void>("set", memoryView);
    return out;
}

// Remove NamedShapes whose geometry is null/empty. Used where an empty result
// is legitimate (e.g. a planar core has no bobbin) so we can return an empty
// buffer instead of throwing on the downstream empty-export guard.
void drop_null_shapes(std::vector<mvb::NamedShape>& named) {
    named.erase(std::remove_if(named.begin(), named.end(),
                    [](const mvb::NamedShape& ns) { return ns.shape.IsNull(); }),
                named.end());
}

// paintCoating arrives from JS as an optional trailing argument. To keep the
// long-standing default (turns drawn at the OUTER/insulation diameter), an
// omitted/undefined value MUST resolve to true — existing callers pass no such
// argument and must keep getting the insulation footprint. Pass false to draw
// the CONDUCTING (copper) cross-section for FEM winding-loss meshing.
bool resolvePaintCoating(const val& v) {
    if (v.isUndefined() || v.isNull()) return true;
    return v.as<bool>();
}

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

OpenMagnetics::Magnetic parseEnriched(const std::string& json_str) {
    nlohmann::json j;
    try {
        j = json::parse(json_str);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("parseEnriched: json::parse failed: ") + e.what());
    }
    bool alreadyEnriched =
        j.contains("core") &&
        j.at("core").contains("geometricalDescription") &&
        !j.at("core").at("geometricalDescription").is_null();
    if (alreadyEnriched) {
        json coreJson = j.at("core");
        json coilJson = j.contains("coil") ? j.at("coil") : json::object();
        try {
            OpenMagnetics::Core core(coreJson);
            if (core.get_geometrical_description().has_value() &&
                !core.get_geometrical_description()->empty()) {
                OpenMagnetics::Coil coil(coilJson, false);
                OpenMagnetics::Magnetic om;
                om.set_core(core);
                om.set_coil(coil);
                return om;
            }
        } catch (const std::exception& e) {
            // Fast path failed (e.g. Coil JSON has minor schema variance,
            // or Core ctor is stricter than what's in the fixture). Fall
            // through to the full MKF autocomplete pipeline below — it's
            // slower but tolerates more JSON shapes.
            std::fprintf(stderr,
                "[mvbpp] parseEnriched: fast path failed (%s); "
                "falling back to magnetic_autocomplete_safe\n", e.what());
        }
    }
    try {
        return mvb::magnetic_autocomplete_safe(j);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("parseEnriched[autocomplete]: ") + e.what());
    }
}

std::string encode_to_bytes(const std::vector<mvb::NamedShape>& named,
                            const std::string& format) {
    std::string fmt = format;
    std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    // Use .compare() and explicit char-by-char checks rather than
    // operator== — under some Emscripten/libc++ builds the implicit
    // const char*/std::string operator== mis-compiles to always-false
    // for short literals; .compare() is reliable.
    bool isStl  = (fmt.size() == 3 && fmt[0]=='s' && fmt[1]=='t' && fmt[2]=='l');
    bool isStep = (fmt.size() == 4 && fmt[0]=='s' && fmt[1]=='t' && fmt[2]=='e' && fmt[3]=='p');

    if (isStl) {
        std::vector<TopoDS_Shape> shapes;
        shapes.reserve(named.size());
        for (const auto& ns : named) shapes.push_back(ns.shape);
        std::string data = mvb::exportSTLToBytes(shapes);
        if (data.empty()) throw std::runtime_error("mvbpp: STL export produced empty output");
        return data;
    }
    if (isStep) {
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
    // Diagnostic: dump hex bytes of `fmt` so we can tell what's actually
    // being compared when this branch is unexpectedly reached.
    std::string hex;
    for (unsigned char c : fmt) {
        char buf[6]; std::snprintf(buf, sizeof(buf), "%02x ", c);
        hex += buf;
    }
    throw std::runtime_error("mvbpp: unsupported format '" + format +
                             "' (lc='" + fmt + "' len=" + std::to_string(fmt.size()) +
                             " hex=" + hex + ") (expected 'step' or 'stl')");
}

// Apply the full delivery pipeline: scale → symmetry → 2D cut → side filter → encode
std::string deliver(std::vector<mvb::NamedShape> named,
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

    auto axisSign = mvb::parse_side_spec(side);
    named = mvb::filter_by_side(named, axisSign);
    if (named.empty()) {
        throw std::runtime_error(
            "mvbpp: side='" + side + "' filtered out all geometry");
    }

    return encode_to_bytes(named, fmt);
}

// ── Per-function builders ───────────────────────────────────────────────────

std::vector<mvb::NamedShape> build_core(const std::string& json_str, int polygonSegments) {
    nlohmann::json j;
    try { j = json::parse(json_str); }
    catch (const std::exception& e) {
        throw std::runtime_error(std::string("drawCore: json::parse failed: ") + e.what());
    }

    // Accept either a Magnetic (has "core" key) or a MagneticCore directly.
    // In both cases enrich (if needed) so we always get a populated
    // geometricalDescription before calling buildCoreNamed.
    bool isMagnetic = j.contains("core");
    if (!isMagnetic) {
        // Wrap the core into a minimal magnetic so MKF autocomplete runs.
        nlohmann::json wrap;
        wrap["core"] = j;
        wrap["coil"] = json{
            {"bobbin", "Dummy"},
            {"functionalDescription", json::array()},
        };
        j = std::move(wrap);
    }

    auto enriched = parseEnriched(j.dump());
    const auto& core = enriched.get_core();
    mvb::MagneticBuilder b;
    return b.buildCoreNamed(core, polygonSegments);
}

std::vector<mvb::NamedShape> build_spacer(const std::string& json_str, int /*polygonSegments*/) {
    auto enriched = parseEnriched(json_str);
    const auto& core = enriched.get_core();
    auto gdOpt = core.get_geometrical_description();
    if (!gdOpt || gdOpt->empty()) {
        throw std::runtime_error("drawSpacer: core has no geometricalDescription");
    }
    auto shapes = mvb::SpacerBuilder::buildSpacers(*gdOpt);
    std::vector<mvb::NamedShape> out;
    out.reserve(shapes.size());
    int i = 0;
    for (auto& s : shapes) {
        if (s.IsNull()) continue;
        out.push_back(mvb::NamedShape{ std::move(s), "Spacer_" + std::to_string(i++) });
    }
    return out;
}

std::vector<mvb::NamedShape> build_board(const std::string& json_str, int /*polygonSegments*/) {
    auto enriched = parseEnriched(json_str);
    const auto& coilOm = enriched.get_coil();
    const auto& core   = enriched.get_core();

    auto groupsOpt = coilOm.get_groups_description();
    if (!groupsOpt || groupsOpt->empty()) return {};

    // Pull the bobbin's processed description and patch from the core's
    // centre column when MKF left the bobbin variant as a string (common
    // after enrichment of raw MAS planar fixtures). Mirrors the helper
    // used in MagneticBuilder::buildAllNamed; kept inline here so this
    // binding doesn't need a private MagneticBuilder header.
    MAS::CoreBobbinProcessedDescription bobbinPd;
    if (const auto* b = std::get_if<OpenMagnetics::Bobbin>(&coilOm.get_bobbin())) {
        if (auto pd = b->get_processed_description()) bobbinPd = *pd;
    }
    if (bobbinPd.get_column_width().value_or(0.0) <= 0.0) {
        auto corePd = core.get_processed_description();
        if (corePd && !corePd->get_columns().empty()) {
            const auto& centralCol = corePd->get_columns()[0];
            double wallThickness = bobbinPd.get_wall_thickness();
            if (std::isnan(wallThickness) || wallThickness < 0.0) wallThickness = 0.0;
            if (centralCol.get_width() > 0.0)
                bobbinPd.set_column_width(centralCol.get_width() / 2.0 + wallThickness);
            if (centralCol.get_depth() > 0.0)
                bobbinPd.set_column_depth(centralCol.get_depth() / 2.0 + wallThickness);
            bobbinPd.set_column_shape(centralCol.get_shape());
        }
    }

    auto shape = mvb::FR4Builder::buildFR4Board(*groupsOpt, bobbinPd);
    std::vector<mvb::NamedShape> out;
    if (!shape.IsNull()) out.push_back(mvb::NamedShape{ std::move(shape), "FR4Board" });
    return out;
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
    return {b.buildBobbinNamedFromBobbin(bobbin, /*axisIsY=*/true, polygonSegments)};
}

std::vector<mvb::NamedShape> build_turns(const std::string& json_str, int polygonSegments,
                                         bool paintCoating) {
    auto j = json::parse(json_str);
    mvb::MagneticBuilder b;
    if (j.is_array()) {
        // Standalone path: each Turn must carry its own dimensions and
        // cross_sectional_shape. Throws for toroidal turns (which need
        // bobbin context); callers should hand a full magnetic instead.
        // paintCoating=false also throws here — a bare Turn has no copper data.
        auto turns = j.get<std::vector<MAS::Turn>>();
        return b.buildTurnsNamedFromTurns(turns, polygonSegments, paintCoating);
    }
    if (j.is_object() && (j.contains("coil") || j.contains("core"))) {
        // Magnetic JSON: enrich + use bobbin-aware overload so toroidal
        // and concentric turns both work uniformly.
        auto enriched = parseEnriched(json_str);
        return b.buildTurnsNamed(enriched.get_coil(), enriched.get_core(),
                                  polygonSegments, paintCoating);
    }
    throw std::runtime_error(
        "drawTurns: input must be a JSON array of MAS::Turn objects, or a "
        "Magnetic JSON with coil + core (required for toroidal turns).");
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

// Adaptive wire-segment count: a coil with N turns produces N solids each
// with the requested polygonSegments around their wire cross-section, then
// every face goes through cut_bobbin (BOP) and BRepMesh. At wireSeg=32 with
// ~80 turns this is ~40 s and OOMs the WASM heap. Downscale aggressively so
// the geometry stays renderable while preserving recognisable wire shape.
// The core gets the caller's value unchanged (core is one big shape, cheap).
static int adaptive_wire_segments(int requested, std::size_t numTurns) {
    if (numTurns <= 20)  return requested;
    if (numTurns <= 40)  return std::min(requested, 12);
    if (numTurns <= 60)  return std::min(requested, 8);
    return std::min(requested, 6);
}

std::vector<mvb::NamedShape> build_magnetic(const std::string& json_str, int polygonSegments,
                                            bool paintCoating) {
    auto j = json::parse(json_str);
    auto magnetic = j.get<MAS::Magnetic>();

    std::size_t numTurns = 0;
    auto turnsOpt = magnetic.get_coil().get_turns_description();
    if (turnsOpt.has_value()) numTurns = turnsOpt->size();
    const int wireSeg = adaptive_wire_segments(polygonSegments, numTurns);

    mvb::MagneticBuilder b;
    return b.buildAllNamed(magnetic, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                            wireSeg, polygonSegments, paintCoating);
}

std::string draw_view_impl(const std::string& json_str,
                            bool dimensions,
                            const std::string& plane,
                            double offset,
                            double widthPx,
                            const std::string& format,
                            double labelFontPx,
                            const std::string& projColor,
                            const std::string& dimColor) {
    std::string fmt = format;
    std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (fmt != "svg") {
        throw std::runtime_error("drawView: only format='svg' is supported (got '" + format + "')");
    }

    auto j = json::parse(json_str);

    OpenMagnetics::Magnetic magnetic;
    bool is_magnetic = j.contains("coil") || j.contains("core") ||
                       j.contains("magnetizing_inductance");
    if (is_magnetic) {
        magnetic = parseEnriched(json_str);
    } else if (j.contains("functionalDescription") || j.contains("functional_description")) {
        json wrap;
        wrap["core"] = j;
        wrap["coil"] = json{
            {"bobbin", "Dummy"},
            {"functionalDescription", json::array()}
        };
        magnetic = mvb::magnetic_autocomplete_safe(wrap);
    } else {
        throw std::runtime_error(
            "drawView: input JSON does not look like a MAS::Magnetic or "
            "MAS::MagneticCore (no functionalDescription / coil / core keys)");
    }

    return mvb::SectionDrawing::drawView(magnetic, plane, offset,
                                          dimensions, widthPx,
                                          labelFontPx, projColor, dimColor);
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Embind entry points — free functions so guard<&fn> works
// ─────────────────────────────────────────────────────────────────────────────
namespace {

// drawCore
val _drawCore(const std::string& json_str,
              const std::string& mode,
              const std::string& plane,
              double offset,
              const std::string& format,
              double scale,
              int polygonSegments,
              const std::string& symmetry,
              const std::string& side) {
    auto named = build_core(json_str, polygonSegments);
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    return makeUint8Array(data);
}

std::string _drawCoreToPath(const std::string& json_str,
                            const std::string& path,
                            const std::string& mode,
                            const std::string& plane,
                            double offset,
                            const std::string& format,
                            double scale,
                            int polygonSegments,
                            const std::string& symmetry,
                            const std::string& side) {
    auto named = build_core(json_str, polygonSegments);
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    write_bytes(path, data);
    return path;
}

// drawCorePiece
val _drawCorePiece(const std::string& json_str,
                   const std::string& mode,
                   const std::string& plane,
                   double offset,
                   const std::string& format,
                   double scale,
                   int polygonSegments,
                   const std::string& symmetry,
                   const std::string& side) {
    auto named = build_core_piece(json_str, polygonSegments);
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    return makeUint8Array(data);
}

std::string _drawCorePieceToPath(const std::string& json_str,
                                 const std::string& path,
                                 const std::string& mode,
                                 const std::string& plane,
                                 double offset,
                                 const std::string& format,
                                 double scale,
                                 int polygonSegments,
                                 const std::string& symmetry,
                                 const std::string& side) {
    auto named = build_core_piece(json_str, polygonSegments);
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    write_bytes(path, data);
    return path;
}

// drawBobbin
val _drawBobbin(const std::string& json_str,
                const std::string& mode,
                const std::string& plane,
                double offset,
                const std::string& format,
                double scale,
                int polygonSegments,
                const std::string& symmetry,
                const std::string& side) {
    auto named = build_bobbin(json_str, polygonSegments);
    // A planar / PCB winding (and many MAS examples — ER25/15 planar, PQ3230,
    // E70/33 stacked) has no physical bobbin: buildBobbin returns an empty
    // shape. That is NOT an error, so return an empty buffer and let the caller
    // simply render no bobbin, rather than throwing "STL export produced empty
    // output" (which surfaced as "[mvbpp] unknown C++ exception").
    drop_null_shapes(named);
    if (named.empty()) return makeUint8Array(std::string());
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    return makeUint8Array(data);
}

std::string _drawBobbinToPath(const std::string& json_str,
                              const std::string& path,
                              const std::string& mode,
                              const std::string& plane,
                              double offset,
                              const std::string& format,
                              double scale,
                              int polygonSegments,
                              const std::string& symmetry,
                              const std::string& side) {
    auto named = build_bobbin(json_str, polygonSegments);
    // See _drawBobbin: an empty bobbin (planar / no-bobbin cores) is valid.
    drop_null_shapes(named);
    if (named.empty()) { write_bytes(path, std::string()); return path; }
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    write_bytes(path, data);
    return path;
}

// drawTurns
val _drawTurns(const std::string& json_str,
               const std::string& mode,
               const std::string& plane,
               double offset,
               const std::string& format,
               double scale,
               int polygonSegments,
               const std::string& symmetry,
               const std::string& side,
               val paintCoating) {
    auto named = build_turns(json_str, polygonSegments, resolvePaintCoating(paintCoating));
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    return makeUint8Array(data);
}

std::string _drawTurnsToPath(const std::string& json_str,
                             const std::string& path,
                             const std::string& mode,
                             const std::string& plane,
                             double offset,
                             const std::string& format,
                             double scale,
                             int polygonSegments,
                             const std::string& symmetry,
                             const std::string& side,
                             val paintCoating) {
    auto named = build_turns(json_str, polygonSegments, resolvePaintCoating(paintCoating));
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    write_bytes(path, data);
    return path;
}

// drawWinding
val _drawWinding(const std::string& json_str,
                 const std::string& windingName,
                 const std::string& mode,
                 const std::string& plane,
                 double offset,
                 const std::string& format,
                 double scale,
                 int polygonSegments,
                 const std::string& symmetry,
                 const std::string& side,
                 val paintCoating) {
    auto named = build_winding(json_str, windingName, polygonSegments,
                               resolvePaintCoating(paintCoating));
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    return makeUint8Array(data);
}

std::string _drawWindingToPath(const std::string& json_str,
                               const std::string& windingName,
                               const std::string& path,
                               const std::string& mode,
                               const std::string& plane,
                               double offset,
                               const std::string& format,
                               double scale,
                               int polygonSegments,
                               const std::string& symmetry,
                               const std::string& side,
                               val paintCoating) {
    auto named = build_winding(json_str, windingName, polygonSegments,
                               resolvePaintCoating(paintCoating));
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    write_bytes(path, data);
    return path;
}

// drawSpacer (takes a Magnetic JSON, builds spacer boxes from
// core.geometricalDescription entries of type SPACER)
val _drawSpacer(const std::string& json_str,
                const std::string& mode,
                const std::string& plane,
                double offset,
                const std::string& format,
                double scale,
                int polygonSegments,
                const std::string& symmetry,
                const std::string& side) {
    auto named = build_spacer(json_str, polygonSegments);
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    return makeUint8Array(data);
}

std::string _drawSpacerToPath(const std::string& json_str,
                              const std::string& path,
                              const std::string& mode,
                              const std::string& plane,
                              double offset,
                              const std::string& format,
                              double scale,
                              int polygonSegments,
                              const std::string& symmetry,
                              const std::string& side) {
    auto named = build_spacer(json_str, polygonSegments);
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    write_bytes(path, data);
    return path;
}

// drawBoard (takes a Magnetic JSON, builds the FR4 PCB plate for planar
// transformer coils — empty when the coil is not PRINTED.)
val _drawBoard(const std::string& json_str,
               const std::string& mode,
               const std::string& plane,
               double offset,
               const std::string& format,
               double scale,
               int polygonSegments,
               const std::string& symmetry,
               const std::string& side) {
    auto named = build_board(json_str, polygonSegments);
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    return makeUint8Array(data);
}

std::string _drawBoardToPath(const std::string& json_str,
                             const std::string& path,
                             const std::string& mode,
                             const std::string& plane,
                             double offset,
                             const std::string& format,
                             double scale,
                             int polygonSegments,
                             const std::string& symmetry,
                             const std::string& side) {
    auto named = build_board(json_str, polygonSegments);
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    write_bytes(path, data);
    return path;
}

// drawMagnetic
val _drawMagnetic(const std::string& json_str,
                  const std::string& mode,
                  const std::string& plane,
                  double offset,
                  const std::string& format,
                  double scale,
                  int polygonSegments,
                  const std::string& symmetry,
                  const std::string& side,
                  val paintCoating) {
    auto named = build_magnetic(json_str, polygonSegments, resolvePaintCoating(paintCoating));
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    return makeUint8Array(data);
}

std::string _drawMagneticToPath(const std::string& json_str,
                                const std::string& path,
                                const std::string& mode,
                                const std::string& plane,
                                double offset,
                                const std::string& format,
                                double scale,
                                int polygonSegments,
                                const std::string& symmetry,
                                const std::string& side,
                                val paintCoating) {
    auto named = build_magnetic(json_str, polygonSegments, resolvePaintCoating(paintCoating));
    auto data = deliver(std::move(named), mode, plane, offset,
                        format, scale, symmetry, side);
    write_bytes(path, data);
    return path;
}

// drawView
std::string _drawView(const std::string& json_str,
                      bool dimensions,
                      const std::string& plane,
                      double offset,
                      double widthPx,
                      const std::string& format) {
    return draw_view_impl(json_str, dimensions, plane, offset, widthPx, format,
                          /*labelFontPx=*/12.0, /*projColor=*/"#000000", /*dimColor=*/"#1976d2");
}

std::string _drawViewToPath(const std::string& json_str,
                            const std::string& path,
                            bool dimensions,
                            const std::string& plane,
                            double offset,
                            double widthPx,
                            const std::string& format) {
    auto svg = draw_view_impl(json_str, dimensions, plane, offset, widthPx, format,
                              12.0, "#000000", "#1976d2");
    write_bytes(path, svg);
    return path;
}

// Dimensioned 2D view with caller-supplied colors (used by Core2DVisualizer via
// the worker). Separate from drawView so the plain drawView keeps its arity and
// existing callers are unaffected.
std::string _drawDimensionedView(const std::string& json_str,
                                 const std::string& plane,
                                 double widthPx,
                                 double labelFontPx,
                                 const std::string& projColor,
                                 const std::string& dimColor) {
    return draw_view_impl(json_str, /*dimensions=*/true, plane, /*offset=*/0.0,
                          widthPx, "svg", labelFontPx, projColor, dimColor);
}

// Test helper
std::string _enrichMagnetic(const std::string& json_str) {
    auto j = nlohmann::json::parse(json_str);
    auto enriched = mvb::magnetic_autocomplete_safe(j);
    nlohmann::json out;
    OpenMagnetics::to_json(out, enriched);
    return out.dump();
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Embind registration
// ─────────────────────────────────────────────────────────────────────────────
EMSCRIPTEN_BINDINGS(mvbpp) {

    register_vector<std::string>("VectorString");

    // 7-function unified API
    function("drawCore",            &guard<&_drawCore>::call);
    function("drawCoreToPath",      &guard<&_drawCoreToPath>::call);
    function("drawCorePiece",       &guard<&_drawCorePiece>::call);
    function("drawCorePieceToPath", &guard<&_drawCorePieceToPath>::call);
    function("drawBobbin",          &guard<&_drawBobbin>::call);
    function("drawBobbinToPath",    &guard<&_drawBobbinToPath>::call);
    function("drawTurns",           &guard<&_drawTurns>::call);
    function("drawTurnsToPath",     &guard<&_drawTurnsToPath>::call);
    function("drawWinding",         &guard<&_drawWinding>::call);
    function("drawWindingToPath",   &guard<&_drawWindingToPath>::call);
    function("drawSpacer",          &guard<&_drawSpacer>::call);
    function("drawSpacerToPath",    &guard<&_drawSpacerToPath>::call);
    function("drawBoard",           &guard<&_drawBoard>::call);
    function("drawBoardToPath",     &guard<&_drawBoardToPath>::call);
    function("drawMagnetic",        &guard<&_drawMagnetic>::call);
    function("drawMagneticToPath",  &guard<&_drawMagneticToPath>::call);
    function("drawView",            &guard<&_drawView>::call);
    function("drawViewToPath",      &guard<&_drawViewToPath>::call);
    function("drawDimensionedView", &guard<&_drawDimensionedView>::call);

    // Test helper
    function("_enrichMagnetic",     &guard<&_enrichMagnetic>::call);

    // Metadata
    function("getSupportedFamilies", &guard<&mvb::get_supported_families>::call);
}
