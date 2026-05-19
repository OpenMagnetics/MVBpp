# MVB++ Project Plan

## Vision
**MVB++** is a C++23 equivalent of the Python/CadQuery-based `../MVB` library, using OpenCASCADE 7.9.3 directly to generate 3D geometry (cores, bobbins, turns, spacers, FR4 boards) and 2D drawings (sections, projections) from MAS JSON. Outputs match Python MVB references and are consumed by C++, Python (pybind11) and WebFrontend (WASM / Emscripten).

## High-level goals
1. **Self-contained build** — fetch and build OCCT 7.9.3 from source automatically.
2. **Direct MAS.hpp usage** — generated `MAS` C++ types, no wrappers.
3. **Meters internally** — scale to mm only at STEP/STL export time.
4. **MKF integration** — `OpenMagnetics::magnetic_autocomplete_safe` enriches raw MAS JSONs.
5. **Unified Named API** — every shape returned as `NamedShape{shape, name}` so names survive symmetry cuts, STEP export, and binding boundaries.
6. **DrawConfig** — single config struct shared by all bindings.

## Coordinate conventions
- **Concentric cores** (E, T, U, etc.): profile in XY, extruded along Z, rotated `-90°` around X → column axis aligns with **Y**.
- **Toroidal cores**: profile in XY, extruded along Z, origin at toroid center.

## Architecture
```
include/mvb/         — public headers (MagneticBuilder, BobbinBuilder, TurnBuilder,
                         SpacerBuilder, FR4Builder, SectionBuilder, Symmetry,
                         StepExporter, ProjectionDrawing, SectionDrawing, NamedShape, Utils)
src/                 — library implementation
src/shapes/          — one builder per core shape family (Factory + Shape*)
tests/               — Catch2 tests (steps, assemblies, all shapes, gapping, symmetry,
                         battery, toroidal top view) + Python reference generator
tools/               — CLI utilities (step_generator, section, bobbin_generator)
bindings/python/     — pybind11 module + pytest suite
bindings/wasm/       — Emscripten module + JS tests
```

## Status

### Done
- **Infrastructure** — CMake fetches OCCT 7.9.3, MKF (+ submodules `cci_coords`, `CAS`, `EAS`), MAS, nlohmann/json, Catch2, pybind11. Builds static `libmvb++.a`, shared `libmvb++.so`, Python wheel and WASM module.
- **Core shape builders** — `ShapeE`, `ShapeT`, `ShapeEr`, `ShapeP`, `ShapeEtd`, `ShapeU`, `ShapeToroidal`, plus families wired via `shapes/Factory.cpp` (drives `get_supported_families()` in `Utils.cpp`). 882/882 non-excluded MAS shapes build a non-empty solid.
- **Bobbin** — `BobbinBuilder` handles rectangular and round bobbins, hollow body + flanges + holes, matches Python `StandardBobbin` volumes/bboxes.
- **Turns** — `TurnBuilder` builds round and rectangular wire turns for concentric and toroidal cores. `buildFromTurnAlone` lets bindings render turns without a wire/bobbin lookup.
- **Symmetry** — plane detection, cutting, high-level binding helpers (`filter_by_side`, `apply_symmetry`, spec parsers).
- **Sections & projections** — `SectionBuilder::cut2DFaces`, `parseSectionPlane`, `drawView` dispatcher, `SectionDrawing`, `ProjectionDrawing`.
- **Spacers and FR4** — `SpacerBuilder` for gap spacers, `FR4Builder` for PCB boards.
- **STEP / STL export** — `StepExporter` with mm scaling; STL tests pass.
- **Bindings**
  - Python (pybind11 + scikit-build-core), wheel installable, pytest suite green.
  - WASM (Emscripten): exposes `buildMagneticSTEP`, `buildMagneticSTL`, drawCore/drawTurns/drawSpacer/drawBoard, `parseEnriched` fallback, tagged exception bridging, whole-archive cmrc for embedded data.
- **MAS 1.0 migration** — schema enum casing migrated; `magnetic_autocomplete_safe` logs failures to stderr.
- **No-fallback policy** — required MAS fields throw on missing data instead of using defaults.

### Test state
| Tag                | Result                                | Notes                                      |
|--------------------|---------------------------------------|--------------------------------------------|
| `[step]`           | pass                                  | E/T core + assembly STEPs match references |
| `[stl]`            | pass                                  |                                            |
| `[assembly]`       | pass (4/4)                            | rect_one_turn, etd49_5t                    |
| `[shapes]`         | pass (882 ran, 8 excluded, 0 failed)  | Excluded: 4 `ui` + 3 `pqi` + 1 `ut` (intentional, mirrors MVB Python) |
| `[get_families]`   | pass                                  |                                            |
| `[symmetry]`       | pass (8/8)                            | Can be slow on complex shapes              |
| `[battery]`        | slow (>15min, scans 25+ MAS examples through MKF) | Preexisting perf issue, not a correctness regression |
| `[topview]`        | pass (toroidal 2D)                    |                                            |
| `[json]`           | pass                                  |                                            |
| `[gapping][additive]`     | pass (448 ran, 0 failed)       |                                            |
| `[gapping][subtractive]`  | pass (448 ran, 0 failed)       | Fixed: U/UR/C zero-length gap + UR cylindrical column tool |
| `[gapping][distributed]`  | pass (448 ran, 0 failed)       | Fixed alongside subtractive                |

### Open items

1. **`[symmetry]` / `[battery]` performance**
   - `[battery]` test scans 25+ MAS examples through MKF and can exceed 15 minutes; blocks CI parallelism
   - `[symmetry]` boolean ops on PQ3230 + distributed gapping can exceed 2 minutes
   - Consider face-based gap rendering or cached cut shapes

2. **WASM JavaScript tests are stale**
   - `bindings/wasm/tests/test_mvbpp.js` calls `drawMagneticToBuffer` / `DrawConfig` which are not in the current embind registration
   - Actual exposed API: `buildMagneticSTEP`, `buildMagneticSTL`, `drawCore`, `drawTurns`, `drawSpacer`, `drawBoard`, `parseEnriched`
   - Tests need to be rewritten against the current API

3. **MKF bobbin dimension divergence**
   - `magnetic_autocomplete` produces flange dimensions that differ from Python `StandardBobbin` defaults when the input has no explicit bobbin
   - Currently MVB++ follows MKF values; document as intended behaviour or wire through a "Python-compatible defaults" switch

4. **README / docs**
   - README is light on the binding APIs, DrawConfig fields, and section/symmetry spec syntax
   - AGENTS.md is the up-to-date reference for build + gotchas; user docs lag behind

### Maintenance gotchas (see AGENTS.md for full list)
- `LD_LIBRARY_PATH=build/occt-install/lib` required for the test binary
- MKF `SHARED→STATIC` patch is a fragile string replace in `CMakeLists.txt`
- MKF submodules (`cci_coords`, `CAS`, `EAS`) auto-init in CMake configure; manual fix path documented
- MAS 1.0 enum casing changes — `python3 MAS/scripts/migrate-to-1.0.py` is the migration tool

## File map
| File / Dir | Purpose |
|------------|---------|
| `CMakeLists.txt` | Build, OCCT external project, MKF/MAS fetch + patches |
| `include/mvb/` + `src/` | Public headers and library implementation |
| `src/shapes/Factory.cpp` | Single source of truth for supported families |
| `src/MagneticBuilder.cpp` | Top-level assembly (core + bobbin + turns + spacers) |
| `src/Symmetry.cpp` | Symmetry plane detection + cutting + spec parser |
| `src/SectionBuilder.cpp` + `src/SectionDrawing.cpp` + `src/ProjectionDrawing.cpp` | 2D drawing pipeline |
| `tests/test_all_shapes.cpp` | 882-shape sanity battery |
| `tests/test_all_gapped.cpp` | Per-shape additive/subtractive/distributed gapping |
| `tests/test_symmetry.cpp` | Symmetry plane and cutting |
| `tests/test_mas_battery.cpp` | Wide MAS-fixture coverage |
| `tests/test_toroidal_topview.cpp` | Toroidal 2D top view |
| `tools/mvbpp_step_generator.cpp` | CLI single/batch STEP generator |
| `tools/mvbpp_section.cpp` | Section/projection CLI |
| `bindings/python/` | pybind11 module + pytest suite |
| `bindings/wasm/` | Emscripten module + (stale) JS tests |
| `AGENTS.md` | Build instructions, gotchas, conventions |
| `PLAN.md` | This file |
