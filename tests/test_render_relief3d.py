"""Tests for the 3D perspective relief rendering (render/relief3d.py) and the
chi-k_sn-over-3D workflow helper.

The offline tests build a small synthetic projected DEM and exercise the
vertical-exaggeration maths plus a real ``grdview`` render (needs GMT/PyGMT +
Ghostscript, like the other render smoke tests — run under the pixi env). The
chi workflow test is network- and lsdtt3-gated.
"""

import numpy as np
import pytest
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from pyproj import CRS

from lsdviztools3.driver import config
from lsdviztools3.render import relief3d
from lsdviztools3.render.relief3d import _zsize_from_ve, render_relief3d
from lsdviztools3.render.style import MapStyle


def _synthetic_dem(path, *, epsg=32630, n=40, res=30.0, height=500.0):
    """Write a small north-up projected DEM (a smooth Gaussian hill) to a tif."""
    x = 400000.0 + res * np.arange(n)
    y = 5650000.0 + res * np.arange(n)[::-1]  # descending => north-up
    xx, yy = np.meshgrid(x, y)
    cx, cy = x.mean(), y.mean()
    sigma = res * n / 5.0
    z = height * np.exp(-(((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sigma ** 2)))
    da = xr.DataArray(z.astype("float32"), dims=("y", "x"),
                      coords={"x": x, "y": y})
    da = da.rio.write_crs(f"EPSG:{epsg}")
    da.rio.to_raster(path)
    return path


# --- _zsize_from_ve: pure function ------------------------------------------

def test_zsize_from_ve_matches_formula():
    """zsize_cm = ve * zrange_m * width_cm / (xrange_km * 1000)."""
    crs = CRS.from_epsg(32630)
    # extent km: x-range 10 km; z-range 500 m; width 12 cm; ve 2 -> 1.2 cm.
    zsize = _zsize_from_ve(2.0, [0.0, 10.0, 0.0, 8.0], 500.0, "12c", crs)
    assert zsize.endswith("c")
    assert float(zsize[:-1]) == pytest.approx(1.2)


def test_zsize_from_ve_scales_with_ve():
    """Doubling the exaggeration doubles the plotted z height."""
    crs = CRS.from_epsg(32630)
    z1 = float(_zsize_from_ve(1.0, [0.0, 10.0, 0.0, 8.0], 500.0, "12c", crs)[:-1])
    z2 = float(_zsize_from_ve(2.0, [0.0, 10.0, 0.0, 8.0], 500.0, "12c", crs)[:-1])
    assert z2 == pytest.approx(2.0 * z1)


def test_zsize_from_ve_rejects_degenerate_ranges():
    crs = CRS.from_epsg(32630)
    with pytest.raises(ValueError):
        _zsize_from_ve(2.0, [0.0, 0.0, 0.0, 8.0], 500.0, "12c", crs)
    with pytest.raises(ValueError):
        _zsize_from_ve(2.0, [0.0, 10.0, 0.0, 8.0], 0.0, "12c", crs)


# --- render_relief3d: real grdview render -----------------------------------

def test_render_relief3d_writes_png(tmp_path):
    """A 3D relief surface with vertical exaggeration renders to a PNG."""
    dem = _synthetic_dem(tmp_path / "hill.tif")
    fig = render_relief3d(str(dem), coordinates="local",
                          vertical_exaggeration=2.0,
                          style=MapStyle(cmap="geo", transparent=False))
    out = tmp_path / "relief3d.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_render_relief3d_ve_overrides_zsize(tmp_path, monkeypatch):
    """When vertical_exaggeration is set, it drives zsize (zsize arg ignored)."""
    captured = {}
    real = relief3d._zsize_from_ve

    def spy(ve, extent, zrange_m, figure_size, crs):
        captured["ve"] = ve
        return real(ve, extent, zrange_m, figure_size, crs)

    monkeypatch.setattr(relief3d, "_zsize_from_ve", spy)
    dem = _synthetic_dem(tmp_path / "hill.tif")
    render_relief3d(str(dem), coordinates="local", vertical_exaggeration=3.0,
                    zsize="99c", style=MapStyle(transparent=False))
    assert captured.get("ve") == 3.0  # VE path taken, not the raw zsize


# --- render_chi_ksn_relief3d: network + lsdtt3 gated ------------------------

def _chi_available() -> bool:
    try:
        config.resolve_lsdtt_tool("lsdtt-chi-analysis")
        return True
    except FileNotFoundError:
        return False


@pytest.mark.skipif(not _chi_available(),
                    reason="lsdtt-chi-analysis not found.")
@pytest.mark.network
def test_render_chi_ksn_relief3d_eigg(small_dem, tmp_path):
    """Drape the log10 k_sn network over the 3D relief surface -> PNG."""
    from lsdviztools3.workflows import (
        lsdtt3_chi_analysis, render_chi_ksn_relief3d,
    )

    outputs = lsdtt3_chi_analysis(
        small_dem, tmp_path, stem="eigg",
        source_threshold=100, n_largest_basins=3, m_over_n=0.45,
    )
    dem_utm = tmp_path / "eigg_remove_seas.tif"
    assert dem_utm.exists() and outputs["segment_lines"].exists()

    png = render_chi_ksn_relief3d(
        dem_utm, outputs["segment_lines"], tmp_path / "eigg_ksn_3d.png",
        vertical_exaggeration=2.0,
    )
    assert png.exists() and png.stat().st_size > 0
