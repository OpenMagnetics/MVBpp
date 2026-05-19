"""
mvbpp Python binding tests — 7-function unified API.

Run with:  pytest bindings/python/tests/ -v
Build/install the wheel first:
    pip install dist/mvbpp-*.whl --force-reinstall
"""
import json
import os
import tempfile
from pathlib import Path

import pytest

import mvbpp

# Most tests below were written before the `side` kwarg landed and assume
# the full (un-filtered) geometry is returned. The new public default
# `side="+X+Y+Z"` keeps only shapes whose centroid lies in the +X+Y+Z
# octant (FEM-friendly). Pass `side=""` to disable the filter.
NO_FILTER = dict(side="")

TESTS_DIR = Path(__file__).parent
BASIC = TESTS_DIR / "concentric_basic.json"
ETD49 = TESTS_DIR / "ETD49_N87_10uH_5T.json"


# ── helpers ──────────────────────────────────────────────────────────────────

def _magnetic_str(path: Path) -> str:
    """Return the raw `magnetic` object from a MAS file as a JSON string."""
    return json.dumps(json.loads(path.read_text())["magnetic"])


@pytest.fixture(scope="session")
def basic_magnetic_raw() -> str:
    return _magnetic_str(BASIC)


@pytest.fixture(scope="session")
def etd49_magnetic_raw() -> str:
    return _magnetic_str(ETD49)


@pytest.fixture(scope="session")
def basic_enriched(basic_magnetic_raw) -> dict:
    """Enriched MAS::Magnetic dict (geometricalDescription, processed bobbin, …)."""
    return json.loads(mvbpp._enrichMagnetic(basic_magnetic_raw))


@pytest.fixture(scope="session")
def etd49_enriched(etd49_magnetic_raw) -> dict:
    return json.loads(mvbpp._enrichMagnetic(etd49_magnetic_raw))


def is_step(data: bytes) -> bool:
    return data[:12] == b"ISO-10303-21"


def is_stl_binary(data: bytes) -> bool:
    # 80-byte header + 4-byte triangle count
    if len(data) < 84:
        return False
    import struct
    n = struct.unpack("<I", data[80:84])[0]
    return len(data) == 84 + n * 50


def is_svg(text: str) -> bool:
    s = text.lstrip()
    return s.startswith("<?xml") or s.startswith("<svg")


# ── Public API surface ──────────────────────────────────────────────────────

class TestPublicSurface:
    """The module must expose exactly the 7 documented draw* functions."""

    EXPECTED = {
        "drawCore", "drawCorePiece", "drawBobbin", "drawTurns",
        "drawWinding", "drawMagnetic", "drawView",
    }

    def test_all_functions_present(self):
        for name in self.EXPECTED:
            assert hasattr(mvbpp, name), f"missing public function: {name}"

    def test_no_old_snake_case_api(self):
        for legacy in ("draw_magnetic", "draw_magnetic_to_bytes",
                       "draw_dimensioned_front_view",
                       "draw_dimensioned_top_view",
                       "write_dimensioned_front_view",
                       "get_symmetry_planes", "DrawConfig"):
            assert not hasattr(mvbpp, legacy), \
                f"legacy symbol still exposed: {legacy}"


# ── drawMagnetic ─────────────────────────────────────────────────────────────

class TestDrawMagnetic:

    def test_returns_step_bytes(self, basic_magnetic_raw):
        out = mvbpp.drawMagnetic(basic_magnetic_raw, **NO_FILTER)
        assert isinstance(out, bytes) and is_step(out)
        assert len(out) > 1000

    def test_contains_brep_geometry(self, basic_magnetic_raw):
        text = mvbpp.drawMagnetic(basic_magnetic_raw, **NO_FILTER).decode("ascii", "replace")
        assert "CLOSED_SHELL" in text or "ADVANCED_BREP" in text

    def test_etd49_larger_than_basic(self, basic_magnetic_raw, etd49_magnetic_raw):
        b = mvbpp.drawMagnetic(basic_magnetic_raw, **NO_FILTER)
        e = mvbpp.drawMagnetic(etd49_magnetic_raw, **NO_FILTER)
        assert len(e) > len(b)

    def test_writes_file_when_path_given(self, basic_magnetic_raw):
        with tempfile.NamedTemporaryFile(suffix=".step", delete=False) as f:
            path = f.name
        try:
            ret = mvbpp.drawMagnetic(basic_magnetic_raw, path, **NO_FILTER)
            assert ret == path
            data = Path(path).read_bytes()
            assert is_step(data) and len(data) > 1000
        finally:
            os.unlink(path)

    def test_stl_format(self, basic_magnetic_raw):
        out = mvbpp.drawMagnetic(basic_magnetic_raw, format="stl", **NO_FILTER)
        assert isinstance(out, bytes) and is_stl_binary(out)

    def test_scale_mm(self, basic_magnetic_raw):
        out = mvbpp.drawMagnetic(basic_magnetic_raw, scale=1000.0, **NO_FILTER)
        assert is_step(out)

    def test_polygon_segments_low(self, basic_magnetic_raw):
        out = mvbpp.drawMagnetic(basic_magnetic_raw, polygonSegments=8, **NO_FILTER)
        assert is_step(out)

    def test_2d_step_section(self, basic_magnetic_raw):
        out = mvbpp.drawMagnetic(basic_magnetic_raw, mode="2D", **NO_FILTER)
        assert is_step(out)

    def test_2d_stl_rejected(self, basic_magnetic_raw):
        with pytest.raises(Exception, match="2D"):
            mvbpp.drawMagnetic(basic_magnetic_raw, mode="2D", format="stl", **NO_FILTER)


# ── drawCore ─────────────────────────────────────────────────────────────────

class TestDrawCore:

    def test_basic(self, basic_enriched):
        core = json.dumps(basic_enriched["core"])
        out = mvbpp.drawCore(core, **NO_FILTER)
        assert is_step(out) and len(out) > 1000

    def test_raw_core_rejected(self, basic_magnetic_raw):
        raw_core = json.dumps(json.loads(basic_magnetic_raw)["core"])
        with pytest.raises(Exception, match="geometricalDescription"):
            mvbpp.drawCore(raw_core)

    def test_2d_section(self, basic_enriched):
        core = json.dumps(basic_enriched["core"])
        out = mvbpp.drawCore(core, mode="2D", **NO_FILTER)
        assert is_step(out)


# ── drawCorePiece ────────────────────────────────────────────────────────────

class TestDrawCorePiece:

    def test_from_e_shape(self, basic_enriched):
        # functionalDescription.shape is enriched into a full CoreShape dict.
        shape = basic_enriched["core"]["functionalDescription"]["shape"]
        out = mvbpp.drawCorePiece(json.dumps(shape), **NO_FILTER)
        assert is_step(out) and len(out) > 1000

    def test_2d_section(self, basic_enriched):
        shape = basic_enriched["core"]["functionalDescription"]["shape"]
        out = mvbpp.drawCorePiece(json.dumps(shape), mode="2D", **NO_FILTER)
        assert is_step(out)


# ── drawBobbin ───────────────────────────────────────────────────────────────

class TestDrawBobbin:

    def test_basic(self, basic_enriched):
        # After enrichment, coil.bobbin is a full Bobbin dict (not just a
        # name reference).  If it is still a string the test is skipped.
        bobbin = basic_enriched["coil"]["bobbin"]
        if isinstance(bobbin, str):
            pytest.skip("Enriched coil.bobbin is still a name reference")
        out = mvbpp.drawBobbin(json.dumps(bobbin), **NO_FILTER)
        assert is_step(out) and len(out) > 500


# ── drawTurns ────────────────────────────────────────────────────────────────

class TestDrawTurns:

    def test_basic(self, basic_enriched):
        turns = basic_enriched["coil"].get("turnsDescription")
        if not turns:
            pytest.skip("Enriched coil has no turnsDescription")
        out = mvbpp.drawTurns(json.dumps(turns), **NO_FILTER)
        assert is_step(out) and len(out) > 500


# ── drawWinding ──────────────────────────────────────────────────────────────

class TestDrawWinding:

    def test_primary(self, basic_enriched):
        coil_json = json.dumps(basic_enriched["coil"])
        if not basic_enriched["coil"].get("turnsDescription"):
            pytest.skip("Enriched coil has no turnsDescription")
        out = mvbpp.drawWinding(coil_json, "Primary", **NO_FILTER)
        assert is_step(out) and len(out) > 500

    def test_unknown_winding_rejected(self, basic_enriched):
        coil_json = json.dumps(basic_enriched["coil"])
        if not basic_enriched["coil"].get("turnsDescription"):
            pytest.skip("Enriched coil has no turnsDescription")
        with pytest.raises(Exception, match="windingName"):
            mvbpp.drawWinding(coil_json, "DoesNotExist", **NO_FILTER)


# ── drawView ─────────────────────────────────────────────────────────────────

class TestDrawView:

    def test_default_dimensioned_xy(self, basic_magnetic_raw):
        svg = mvbpp.drawView(basic_magnetic_raw)
        assert isinstance(svg, str) and is_svg(svg)

    def test_dimensioned_xz(self, basic_magnetic_raw):
        svg = mvbpp.drawView(basic_magnetic_raw, plane="XZ")
        assert is_svg(svg)

    def test_undimensioned_any_plane(self, basic_magnetic_raw):
        svg = mvbpp.drawView(basic_magnetic_raw, dimensions=False, plane="YZ")
        assert is_svg(svg)

    def test_writes_file(self, basic_magnetic_raw):
        with tempfile.NamedTemporaryFile(suffix=".svg", delete=False) as f:
            path = f.name
        try:
            ret = mvbpp.drawView(basic_magnetic_raw, path)
            assert ret == path
            assert "<svg" in Path(path).read_text()
        finally:
            os.unlink(path)

    def test_format_png_rejected(self, basic_magnetic_raw):
        with pytest.raises(Exception, match="svg"):
            mvbpp.drawView(basic_magnetic_raw, format="png")


# ── error handling ───────────────────────────────────────────────────────────

class TestErrors:

    def test_invalid_json(self):
        with pytest.raises(Exception):
            mvbpp.drawMagnetic("not json at all")

    def test_empty_object(self):
        with pytest.raises(Exception):
            mvbpp.drawMagnetic("{}")

    def test_unknown_format(self, basic_magnetic_raw):
        with pytest.raises(Exception, match="format"):
            mvbpp.drawMagnetic(basic_magnetic_raw, format="obj", **NO_FILTER)

    def test_unknown_mode(self, basic_magnetic_raw):
        with pytest.raises(Exception, match="mode"):
            mvbpp.drawMagnetic(basic_magnetic_raw, mode="4D", **NO_FILTER)


# ── symmetry & side ──────────────────────────────────────────────────────────

class TestSymmetryAndSide:
    """`symmetry` selects how many bisecting planes to apply; `side` filters
    surviving shapes by the sign of their centroid in each axis."""

    def test_default_side_filters_geometry(self, basic_magnetic_raw):
        full = mvbpp.drawMagnetic(basic_magnetic_raw, side="")
        filt = mvbpp.drawMagnetic(basic_magnetic_raw)  # default +X+Y+Z
        # Default must filter at least some shapes out.
        assert len(filt) < len(full)

    def test_symmetry_half_shrinks_3d(self, basic_magnetic_raw):
        full = mvbpp.drawMagnetic(basic_magnetic_raw, side="")
        half = mvbpp.drawMagnetic(basic_magnetic_raw, symmetry="half", side="")
        assert len(half) < len(full)

    def test_symmetry_quarter_shrinks_more_than_half(self, basic_magnetic_raw):
        half = mvbpp.drawMagnetic(basic_magnetic_raw, symmetry="half", side="")
        quart = mvbpp.drawMagnetic(basic_magnetic_raw, symmetry="quarter", side="")
        assert len(quart) < len(half)

    def test_symmetry_int_alias(self, basic_magnetic_raw):
        # "half" and "quarter" must equal their numeric forms when accepted.
        a = mvbpp.drawMagnetic(basic_magnetic_raw, symmetry="quarter", side="")
        # Numeric strings must also parse.
        b = mvbpp.drawMagnetic(basic_magnetic_raw, symmetry="2", side="")
        assert len(a) == len(b)

    def test_unknown_symmetry_rejected(self, basic_magnetic_raw):
        with pytest.raises(Exception, match="symmetry"):
            mvbpp.drawMagnetic(basic_magnetic_raw, symmetry="sextant", side="")

    def test_unknown_side_rejected(self, basic_magnetic_raw):
        with pytest.raises(Exception, match="side"):
            mvbpp.drawMagnetic(basic_magnetic_raw, side="+Q")

    def test_side_filters_in_2d(self, basic_magnetic_raw):
        full = mvbpp.drawMagnetic(basic_magnetic_raw, mode="2D", plane="XZ", side="")
        plus_x = mvbpp.drawMagnetic(basic_magnetic_raw, mode="2D", plane="XZ", side="+X")
        assert len(plus_x) < len(full) and is_step(plus_x)

    def test_fem_quadrant_combo(self, basic_magnetic_raw):
        # The FEM use-case: 2-D XZ section, quarter-symmetry, +X half plane.
        out = mvbpp.drawMagnetic(basic_magnetic_raw,
                                  mode="2D", plane="XZ",
                                  symmetry="quarter", side="+X")
        assert is_step(out)

    def test_side_empty_means_no_filter(self, basic_magnetic_raw):
        a = mvbpp.drawMagnetic(basic_magnetic_raw, side="")
        b = mvbpp.drawMagnetic(basic_magnetic_raw, side="none")
        c = mvbpp.drawMagnetic(basic_magnetic_raw, side="auto")
        assert len(a) == len(b) == len(c)
