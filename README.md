# MVB++

3-D magnetics geometry builder. Reads [MAS](https://github.com/OpenMagnetics/MAS)
JSON, produces solid models (STEP / STL) and 2-D projections (SVG) of cores,
bobbins and windings using OpenCASCADE 7.9.3.

C++23 library with Python (pybind11) and WebAssembly (Emscripten) bindings.
Mirrors the Python/CadQuery [MVB](https://github.com/OpenMagnetics/MVB)
library's API and conventions, but is ~10× faster on the hot paths and
deployable to the browser.

---

## Build

CMake fetches every dependency automatically — OpenCASCADE is built from
source on the first configure (~10 min). Subsequent builds are incremental.

```bash
cmake -S . -B build -GNinja -DCMAKE_BUILD_TYPE=Release
ninja -C build mvb++          # static library
ninja -C build mvb_tests      # Catch2 test runner
```

The test binary links OCCT's shared libs; the library path must be visible
at runtime:

```bash
LD_LIBRARY_PATH=build/occt-install/lib ./build/mvb_tests "[step]"
```

### Python bindings

```bash
pip install -e .                     # scikit-build-core + pybind11
pytest bindings/python/tests/ -v
```

The wheel exposes a `mvbpp` module.

### WebAssembly bindings

Requires the [emsdk](https://emscripten.org/docs/getting_started/downloads.html)
to be active in the shell (`source ~/emsdk/emsdk_env.sh`).

```bash
cd bindings/wasm
npm run build
node tests/test_mvbpp.js
```

The build emits `dist/mvbpp.{js,wasm,d.ts}` ready to be loaded by any
JavaScript host.

---

## Architecture

```
                       ┌──────────────────────┐
                       │   MAS JSON (input)   │
                       └──────────┬───────────┘
                                  │
                  ┌───────────────▼────────────────┐
                  │  MKF magnetic_autocomplete     │   (skipped if MAS
                  │  → enriched MAS::Magnetic      │    already has
                  └───────────────┬────────────────┘    geometricalDescription)
                                  │
        ┌─────────────────────────┼─────────────────────────┐
        │                         │                         │
        ▼                         ▼                         ▼
┌──────────────┐         ┌──────────────────┐      ┌────────────────┐
│ shapes/      │         │ BobbinBuilder    │      │ TurnBuilder    │
│ Factory.cpp  │         │ (rect / round,   │      │ (round / rect, │
│ → ShapeE,    │         │  flanges, holes) │      │  concentric /  │
│   ShapeP,    │         └──────┬───────────┘      │  toroidal)     │
│   ShapeEtd,  │                │                  └───────┬────────┘
│   ShapeT,    │                │                          │
│   ShapeU,    │                │     SpacerBuilder        │
│   ShapeEr,   │                │     FR4Builder           │
│   Toroidal…  │                │                          │
└──────┬───────┘                │                          │
       │                        │                          │
       └────────────┬───────────┴──────────────┬───────────┘
                    ▼                          ▼
            ┌───────────────────────────────────────┐
            │      MagneticBuilder::buildAllNamed   │
            │   → std::vector<NamedShape>           │
            │     (named TopoDS_Shape parts)        │
            └───────────────┬───────────────────────┘
                            │
       ┌────────────────────┼─────────────────────────┐
       ▼                    ▼                         ▼
┌──────────────┐   ┌────────────────────┐    ┌──────────────────┐
│  Symmetry    │   │  SectionBuilder    │    │  StepExporter    │
│  analyze /   │   │  cut2DFaces (XY/   │    │  STEP / STL      │
│  apply /     │   │  XZ/YZ, offset)    │    │  (m → mm)        │
│  filter_by_  │   │  → planar Faces    │    └────────┬─────────┘
│  side        │   └────────┬───────────┘             │
└──────┬───────┘            │                         │
       │                    ▼                         │
       │            ┌───────────────────┐             │
       │            │ Section/Projection│             │
       │            │ Drawing → SVG     │             │
       │            └────────┬──────────┘             │
       │                     │                        │
       └─────────────────────┴────────────┬───────────┘
                                          ▼
                       ┌──────────────────────────────────┐
                       │  Bindings: C++ / Python / WASM   │
                       │  (DrawConfig + draw* dispatch)   │
                       └──────────────────────────────────┘
```

Data flows top-to-bottom: MAS JSON in, enriched by MKF, decomposed into
named OCCT shapes by `MagneticBuilder`, then post-processed (symmetry,
2-D sectioning) and exported. All three binding surfaces (`mvb::`,
`mvbpp` Python module, `mvbpp` WASM module) consume the same
`DrawConfig` and return the same `NamedShape` set.

---

## C++ API

Everything lives in namespace `mvb`. The entry point is `MagneticBuilder`.

```cpp
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"

mvb::MagneticBuilder builder;

// Load and enrich a MAS JSON. The enrichment runs MKF's
// magnetic_autocomplete, which fills in derived geometric data
// (gap descriptions, processed bobbin, turn coordinates).
auto enriched = mvb::magnetic_autocomplete_safe(json_blob);

// 1) Write a STEP file (returns the output path).
mvb::DrawConfig cfg;
cfg.format                = "step";
cfg.includeBobbin         = true;
cfg.scale                 = 1000.0;       // m -> mm at export time
cfg.symmetryPlanes        = 0;            // 1 = half, 2 = quarter
cfg.wirePolygonSegments   = 8;
cfg.corePolygonSegments   = 32;
auto path = builder.drawMagnetic(enriched, "/tmp/out", cfg);

// 2) Build named TopoDS_Shape pieces (no I/O) for downstream use.
auto core   = builder.buildCoreNamed  (enriched.get_core());
auto bobbin = builder.buildBobbinNamed(enriched.get_coil(),
                                       enriched.get_core());
auto turns  = builder.buildTurnsNamed (enriched.get_coil(),
                                       enriched.get_core());
auto all    = builder.buildAllNamed   (enriched);
```

### `DrawConfig` fields

| Field                  | Default     | Meaning                                              |
|------------------------|-------------|------------------------------------------------------|
| `format`               | `"step"`    | `"step"` or `"stl"`                                  |
| `includeBobbin`        | `true`      | Skip bobbin geometry when `false`                    |
| `scale`                | `1.0`       | Uniform scale; pass `1000.0` to export in mm         |
| `symmetryPlanes`       | `0`         | `0` full, `1` half (1 plane), `2` quarter (2 planes) |
| `wirePolygonSegments`  | 8           | `≤0` = exact NURBS torus, `>0` = faceted wire        |
| `corePolygonSegments`  | 32          | Polygon facet count for round cores / bobbins        |

### Symmetry

`mvb::analyze_symmetry(shapes)` reports which of the XY/YZ/XZ planes the
assembly is mirror-symmetric across (volume-based, AABB-prefiltered).
`mvb::cut_to_region(shapes, cuts, bbox)` returns the assembly trimmed to the
selected half-spaces. `apply_symmetry(shapes, N)` is the shortcut used by
`buildAllNamed` when `symmetryPlanes > 0`.

```cpp
auto sym  = mvb::analyze_symmetry(all);                       // {planes, candidates}
auto half = mvb::apply_symmetry(all, /*nPlanes=*/1);          // = cut to half along sym.valid_planes[0]
```

---

## Python API

```python
import mvbpp, json
mas = open("examples/01_simple_inductor_etd34_n87.json").read()

# Write a STEP file.
mvbpp.drawMagnetic(mas, "/tmp/out",
                   mode="3D", format="step", scale=1000.0,
                   polygonSegments=32)

# Return bytes (omit outputPath).
step_bytes = mvbpp.drawMagnetic(mas, None, format="step")

# Render a 2-D section as SVG.
svg_bytes  = mvbpp.drawView(mas, None,
                            dimensions=True, plane="XZ", offset=0.0,
                            widthPx=800, format="svg")
```

All `draw*` functions share the same trailing keyword arguments:
`mode={"3D","2D"}`, `plane={"XY","XZ","YZ"}`, `offset`, `format`,
`scale`, `polygonSegments`, `symmetry`, `side`.

- `symmetry` (default `"none"`) — `"none"`/`"0"` full, `"half"`/`"1"` cuts
  along one symmetry plane, `"quarter"`/`"2"` cuts along two. Unknown
  values raise.
- `side` (default `"+X+Y+Z"`, FEM-friendly positive octant) — keeps only
  shapes whose centroid lies in the requested half-spaces. Combine axes
  like `"+X"`, `"+X-Y"`, `"+X+Y+Z"`. Pass `""`, `"none"`, or `"auto"` to
  disable filtering. Unknown tokens raise.

Available entry points:

| Function           | Input JSON                                      |
|--------------------|-------------------------------------------------|
| `drawCore`         | `MagneticCore`                                  |
| `drawCorePiece`    | `CoreShape`                                     |
| `drawBobbin`       | `Bobbin`                                        |
| `drawTurns`        | array of `Turn`                                 |
| `drawWinding`      | `Coil` (with `windingName=` arg)                |
| `drawMagnetic`     | `Magnetic` (full assembly)                      |
| `drawView`         | any of the above (2-D outline / dimensioned)    |

`mode="2D"` cuts at `plane` shifted by `offset` (m) and exports a STEP
compound of planar `TopoDS_Face`s, suitable for ElmerFEM surface meshing.
STL is not allowed in 2-D mode.

---

## WebAssembly API

The embind surface mirrors the Python API and adds two extras
(`drawSpacer`, `drawBoard`). Each `draw*` returns a `Uint8Array`;
each `*ToPath` variant writes to a virtual-FS path inside the WASM module.

```js
import mvbppFactory from "./dist/mvbpp.js";
const mvbpp = await mvbppFactory();

const stepBytes = mvbpp.drawMagnetic(
    masJson,
    /* mode */         "3D",
    /* plane */        "XY",
    /* offset */       0.0,
    /* format */       "step",
    /* scale */        1000.0,
    /* polygonSegments*/ 32,
    /* symmetry */     "none",
    /* side */         "+X+Y+Z"
);
```

The 9 trailing args (`mode, plane, offset, format, scale, polygonSegments,
symmetry, side`) are required positional arguments; they map 1-to-1 onto
the Python keyword arguments above (same accepted values, same defaults
when called from Python). See `bindings/wasm/mvbpp_wasm.cpp` for the full
signature and `bindings/wasm/tests/test_mvbpp.js` for end-to-end usage.

The test helper `_enrichMagnetic(masJson)` exposes `magnetic_autocomplete`
for callers that want the MKF-enriched MAS JSON without going through a
draw call.

---

## CLI tools

Built alongside the library when targeting `mvbpp_step_generator`,
`mvbpp_section`, `mvbpp_generate_bobbins`.

```bash
# Single STEP export.
./build/mvbpp_step_generator -o /tmp/out.step examples/01_etd34.json

# Front + top + XY + XZ SVG projections.
./build/mvbpp_section input.json /tmp/svgs/
```

---

## Tests

The test binary supports Catch2 tag filters. Common combinations:

```bash
./build/mvb_tests "[step]"               # STEP reference comparisons (~10 s)
./build/mvb_tests "[stl]"                # STL reference comparisons (~10 s)
./build/mvb_tests "[assembly]"           # Full magnetic assemblies (~30 s)
./build/mvb_tests "[shapes]"             # Iterates 882 MAS core shapes
./build/mvb_tests "[symmetry]"           # ~1:35 — boolean cuts on heavy ETD/PQ
./build/mvb_tests "[gapping]"            # ~5 min — additive/subtractive/distributed
./build/mvb_tests "[battery][simple]"    # 25 MAS examples, full pipeline
./build/mvb_tests "[topview]"            # Toroidal 2-D top view
./build/mvb_tests "[json]"               # JSON round-trip
```

Python and WASM tests live under `bindings/{python,wasm}/tests/` and are
invoked via `pytest` and `node` respectively.

### Battery env vars

`tests/test_mas_battery.cpp` honours a few overrides handy for debugging
perf:

| Var                              | Effect                                                |
|----------------------------------|-------------------------------------------------------|
| `MVB_BATTERY_FILTER=<substring>` | Run only fixtures whose filename contains substring   |
| `MVB_BATTERY_DEBUG_PAIRS=1`      | Dump every overlap-checked pair with its volume       |
| `MVB_BATTERY_TURN_CORE_STRIDE=N` | Override turn-vs-core sampling stride (default 16)    |

---

## Conventions and gotchas

* **Internal units are metres.** Scale to mm only at export (`DrawConfig.scale = 1000`).
* **Concentric cores** (E, T, U, …) are laid out with the column along Y; the
  generator builds the cross-section in XY, extrudes along Z, then rotates
  −90° around X. **Toroidal** cores stay in their natural XY/Z frame.
* **No silent fallbacks.** Missing MAS data raises an exception rather than
  substituting a default; see the project-wide policy in `AGENTS.md` and
  `CLAUDE.md`.
* OCCT boolean ops must always check `IsDone()` before `Shape()` — null
  shapes propagate silently otherwise. The library does this; callers
  reaching into the OCCT API themselves should too.
* `LD_LIBRARY_PATH=build/occt-install/lib` is required for the test binary
  and any consumer linking the shared OCCT libs.

Developer reference and build-system internals live in
[`AGENTS.md`](AGENTS.md); the current status, open performance items and
design decisions are tracked in [`PLAN.md`](PLAN.md).
