# AGENTS.md — MVB++

## What this is

C++23 library that builds 3D magnetic component geometry (cores, bobbins, turns) from MAS JSON using OpenCASCADE 7.9.3. Produces STEP/STL files and SVG dimensioned drawings. Has Python (pybind11) and WASM (Emscripten) bindings.

## Build

```bash
# Configure (fetches OCCT, MKF, MAS, nlohmann/json, Catch2 automatically)
cmake -S . -B build -GNinja -DCMAKE_BUILD_TYPE=Release

# Build library
cmake --build build --target mvb++ -j$(nproc)

# Build and run tests
cmake --build build --target mvb_tests -j$(nproc)
cd build && LD_LIBRARY_PATH=$PWD/occt-install/lib:$LD_LIBRARY_PATH ./mvb_tests

# Run specific test tags
./mvb_tests "[step][e]"           # E core STEP comparison
./mvb_tests "[assembly][json]"    # Full assembly tests
./mvb_tests "[symmetry]"          # Symmetry analysis (can be slow)
./mvb_tests "[all_shapes]"        # All MAS core shapes
```

First build takes ~10min (OCCT builds from source). Subsequent builds are fast.

### Python bindings

```bash
pip install -e .                    # scikit-build-core + pybind11
pytest bindings/python/tests/ -v   # requires installed mvbpp wheel
```

### WASM bindings

```bash
cd bindings/wasm
npm run build   # requires emsdk (source ~/emsdk/emsdk_env.sh)
npm test
```

## Architecture

```
include/mvb/          Public headers (MagneticBuilder, Utils, Symmetry, StepExporter, etc.)
src/                  Library implementation
src/shapes/           One builder per core shape family (ShapeE, ShapeT, ShapeEr, etc.)
tests/                Catch2 tests + reference STEPs
tools/                CLI utilities (step_generator, section, bobbin_generator)
bindings/python/      pybind11 Python module
bindings/wasm/        Emscripten WASM module
testData/             JSON test fixtures (MAS format)
```

Key types:
- `mvb::MagneticBuilder` — main entry point, builds core/bobbin/turns geometry
- `mvb::NamedShape` — shape + string name (survives symmetry cuts, STEP export)
- `mvb::StepExporter` — STEP/STL export with mm scaling
- `mvb::Symmetry` — symmetry plane detection and cutting
- `MAS::Magnetic` — deserialized MAS JSON (generated from schema)
- `OpenMagnetics::Magnetic` — MKF-enriched magnetic (has geometricalDescription)

## Coordinate conventions

- **Concentric cores** (E, T, U, etc.): profile in XY, extruded along Z, rotated -90° around X → column axis aligns with Y.
- **Toroidal cores**: profile in XY, extruded along Z, origin at toroid center.
- **Meters internally** — scale to mm only at STEP export (Python MVB exports in mm).

## Dependencies

All fetched via CMake FetchContent (no manual install needed):
- **OCCT 7.9.3** — built from source as ExternalProject (static)
- **MKF** — `OpenMagnetics::magnetic_autocomplete` enriches raw MAS JSON
- **MAS** — auto-generated C++ types from JSON schema
- **nlohmann/json v3.11.3**, **Catch2 v3.8.0**, **pybind11 v2.13.6**

MKF has submodules (`cci_coords`, `CAS`, `EAS`) that are auto-initialized during CMake configure.

## Code conventions

- C++23 required (`std::numbers::pi`, `std::optional`, etc.)
- `get_supported_families()` in `Utils.cpp` is implemented in `shapes/Factory.cpp` — derives from the factory switch, don't maintain a parallel list
- OCCT boolean ops (`BRepAlgoAPI_Cut`, `BRepAlgoAPI_Fuse`) must check `IsDone()` before calling `.Shape()` — null shapes propagate silently
- `MagneticBuilder::buildAllNamed(MAS::Magnetic)` skips MKF autocomplete when `geometricalDescription` is already present — use this for pre-enriched data
- `DrawConfig` struct in `MagneticBuilder.h` is the unified config for all bindings — prefer it over positional parameters
- `magnetic_autocomplete_safe` logs failures to `stderr` (not silent `catch(...)`) — check stderr when debugging enrichment issues

## Gotchas

### LD_LIBRARY_PATH required for tests
OCCT builds both static and shared libs. The test binary links shared by default:
```bash
LD_LIBRARY_PATH=build/occt-install/lib:$LD_LIBRARY_PATH ./mvb_tests
```

### MAS schema enum casing (MAS 1.0)
MAS 1.0 changed all enum values to camelCase/lowercase. If you see `"Input JSON does not conform to schema!"`, check:
- `"two-piece set"` → `"twoPieceSet"` (CoreType)
- `"Wound"` → `"wound"` (WiringTechnology)
- `"inner or top"` → `"innerOrTop"` (CoilAlignment)
- `"outer or bottom"` → `"outerOrBottom"` (CoilAlignment)

Use `python3 MAS/scripts/migrate-to-1.0.py <file-or-dir>` to fix old JSON files.

### MKF SHARED→STATIC patching
CMakeLists.txt patches MKF's `add_library(MKF SHARED` → `STATIC` via string replace. If upstream MKF changes this string, the patch silently fails and you get link errors.

### MKF submodules
MKF needs `cci_coords`, `CAS`, `EAS` submodules. CMake auto-initializes them, but if you see `CAS.hpp: No such file or directory`, the submodule init failed. Run manually:
```bash
cd build/_deps/mkf-src && git submodule update --init cci_coords CAS EAS
```

### Symmetry/gapping tests can be slow
Boolean operations on complex shapes (PQ3230, distributed gapping) can take >2min and may timeout. The gapping tests iterate ~461 shapes through MKF. Run separately:
```bash
./mvb_tests "[symmetry]"    # ~3min
./mvb_tests "[gapping]"     # ~5min, iterates all MAS shapes
./mvb_tests "[step]"        # fast, ~10s
./mvb_tests "[stl]"         # fast, ~10s
./mvb_tests "[assembly]"    # fast, ~30s
```

### WASM build/test
```bash
cd bindings/wasm && npm run build   # requires emsdk: source ~/emsdk/emsdk_env.sh
node tests/test_mvbpp.js            # 27 tests, all pass
```
The embind surface is the `draw*` family: `drawMagnetic`, `drawCore`, `drawCorePiece`, `drawBobbin`, `drawTurns`, `drawWinding`, `drawSpacer`, `drawBoard`, `drawView`, plus `*ToPath` variants and the test helper `_enrichMagnetic`. All take the trailing 9-arg config pack `(mode, plane, offset, format, scale, polygonSegments, symmetry, side)` documented in `bindings/wasm/mvbpp_wasm.cpp`.

## Test data

- `testData/` — MAS JSON fixtures for MVB++ tests
- `tests/reference_steps/` — reference STEP/STL files for comparison tests
- `MAS/examples/` — MAS example files (fetched via FetchContent)
- `MAS/samples/` — MAS sample files (fetched via FetchContent)

Test data must conform to MAS 1.0 schema. Run `python3 MAS/scripts/migrate-to-1.0.py testData/` after schema changes.

## Current status

- Core shapes, bobbin, turns, sections/projections, symmetry, spacers, FR4 all implemented
- C++ tests: `[step]/[stl]/[assembly]/[shapes]/[json]/[topview]/[symmetry]/[gapping]` all pass; `[battery]` is slow (>15min) but functional
- Python (pybind11) and WASM (Emscripten) bindings: all tests pass
- See `PLAN.md` for detailed status and open performance/docs items
