"""Tests for orthographic locator basemap rendering."""

import numpy as np
import rioxarray  # noqa: F401  (registers xarray .rio accessor)
import xarray as xr

from lsdviztools3.render.basemap import render_orthographic_basemap
from lsdviztools3.render.style import MapStyle


class _DummyFigure:
    def __init__(self):
        self.calls = []

    def coast(self, **kwargs):
        self.calls.append(("coast", kwargs))

    def plot(self, **kwargs):
        self.calls.append(("plot", kwargs))

    def basemap(self, **kwargs):
        self.calls.append(("basemap", kwargs))


def _tiny_geo_grid():
    da = xr.DataArray(
        np.array([[1.0, 2.0], [3.0, 4.0]]),
        dims=("y", "x"),
        coords={"x": [-6.2, -6.1], "y": [56.8, 56.9]},
    )
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")
    da = da.rio.write_crs("EPSG:4326")
    return da


def test_render_orthographic_basemap_plots_coast_and_bbox():
    fig = _DummyFigure()
    out = render_orthographic_basemap(
        _tiny_geo_grid(),
        style=MapStyle(figure_size="8c"),
        fig=fig,
        mode="aerial",
        show_center=True,
    )

    assert out is fig
    names = [name for name, _ in fig.calls]
    assert names.count("coast") == 1
    assert names.count("basemap") == 1
    assert names.count("plot") == 3  # bbox line + corners + center marker

    coast_kwargs = fig.calls[0][1]
    assert coast_kwargs["projection"].startswith("G")
    assert coast_kwargs["region"] == "g"
    assert coast_kwargs["frame"].startswith("g")
    basemap_kwargs = fig.calls[1][1]
    assert basemap_kwargs["rose"]


def test_render_orthographic_basemap_global_mode_uses_global_region():
    fig = _DummyFigure()
    render_orthographic_basemap(_tiny_geo_grid(), fig=fig, mode="global")
    coast_kwargs = fig.calls[0][1]
    assert coast_kwargs["region"] == "g"
    assert "+z" not in coast_kwargs["projection"]
    # No north arrow on a globe: GMT's rose there broke the PNG export.
    assert "basemap" not in [name for name, _ in fig.calls]


def test_render_orthographic_basemap_aerial_mode_honors_view_limit():
    fig = _DummyFigure()
    render_orthographic_basemap(
        _tiny_geo_grid(),
        fig=fig,
        mode="aerial",
        view_limit_deg=2.0,
        look_bearing=45.0,
        look_vertical_angle=20.0,
        elevation_km=50.0,
    )
    coast_kwargs = fig.calls[0][1]
    proj = coast_kwargs["projection"]
    assert coast_kwargs["region"] == "g"
    assert "+a45.0" in proj
    assert "+t20.0" in proj
    assert "+z50.0" in proj
    assert "+v" in proj
    assert coast_kwargs["frame"] == "g5.0"
    assert fig.calls[1][0] == "basemap"
    assert fig.calls[1][1]["rose"] is not None


def test_render_orthographic_basemap_aerial_no_north_arrow():
    fig = _DummyFigure()
    render_orthographic_basemap(_tiny_geo_grid(), fig=fig, mode="aerial", show_north_arrow=False)
    names = [name for name, _ in fig.calls]
    assert "basemap" not in names


def test_render_orthographic_basemap_local_alias_maps_to_aerial():
    fig = _DummyFigure()
    render_orthographic_basemap(_tiny_geo_grid(), fig=fig, mode="local", elevation_km=60)
    assert "+z60.0" in fig.calls[0][1]["projection"]


def test_render_orthographic_basemap_invalid_mode_raises():
    import pytest

    with pytest.raises(ValueError, match="mode must be one of"):
        render_orthographic_basemap(_tiny_geo_grid(), fig=_DummyFigure(), mode="bad")


def test_render_orthographic_basemap_missing_crs_raises():
    da = xr.DataArray(np.ones((2, 2)), dims=("y", "x"), coords={"x": [0, 1], "y": [0, 1]})
    import pytest

    with pytest.raises(ValueError, match="no CRS metadata"):
        render_orthographic_basemap(da, fig=_DummyFigure())


# --------------------------------------------------------------------------
# Regional mode, the automatic choice, insets, and the commands

import math  # noqa: E402

import pytest  # noqa: E402

from lsdviztools3.render.basemap import (  # noqa: E402
    _nice_length_km,
    add_locator_inset,
    regional_extent,
    render_regional_basemap,
)


class _RegionalDummy(_DummyFigure):
    def grdimage(self, **kwargs):
        self.calls.append(("grdimage", kwargs))


def _geo_grid(west, south, east, north):
    da = xr.DataArray(np.ones((2, 2)), dims=("y", "x"),
                      coords={"x": [west, east], "y": [south, north]})
    return da.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:4326")


def test_auto_mode_is_regional_for_small_and_global_for_large_rasters():
    fig = _RegionalDummy()
    render_orthographic_basemap(_tiny_geo_grid(), fig=fig)
    assert fig.calls[0][0] == "basemap" and fig.calls[0][1]["projection"].startswith("M")
    fig = _RegionalDummy()
    render_orthographic_basemap(_geo_grid(-20, 30, 20, 60), fig=fig)
    assert fig.calls[0][1]["region"] == "g"


def test_regional_extent_is_square_and_centred():
    region, extent = regional_extent(-115.74, 45.31, -115.40, 45.51)
    w, e, s, n = region
    assert extent == pytest.approx(10 * math.hypot(0.34 * 111.32 * math.cos(math.radians(45.41)),
                                                   0.20 * 110.57), rel=1e-3)
    width_km = (e - w) * 111.32 * math.cos(math.radians(45.41))
    height_km = (n - s) * 110.57
    assert width_km == pytest.approx(extent) and height_km == pytest.approx(extent)
    assert 0.5 * (w + e) == pytest.approx(-115.57) and 0.5 * (s + n) == pytest.approx(45.41)
    # Small footprints get at least 250 km; an extent can be given.
    assert regional_extent(0, 0, 0.01, 0.01)[1] == 250.0
    assert regional_extent(0, 0, 0.01, 0.01, 40)[1] == 40.0
    with pytest.raises(ValueError):
        regional_extent(0, 0, 1, 1, 0)


@pytest.mark.parametrize("target,expected", [(90, 50), (100, 100), (37, 20), (4.9, 2), (1200, 1000)])
def test_scale_bar_length_is_a_round_number(target, expected):
    assert _nice_length_km(target) == expected


def test_regional_map_layers_and_options():
    fig = _RegionalDummy()
    render_regional_basemap(_tiny_geo_grid(), fig=fig, style=MapStyle(figure_size="8c"))
    names = [n for n, _ in fig.calls]
    assert names == ["basemap", "coast", "plot", "basemap"]
    coast = fig.calls[1][1]
    assert coast["projection"] == "M8c" and coast["borders"] and coast["rivers"]
    assert coast["land"] == "gray90"
    last = fig.calls[3][1]
    assert last["map_scale"].startswith("jBL") and "+w" in last["map_scale"] and last["rose"]
    # Footprint outline: the raster's bounds.
    assert sorted(set(fig.calls[2][1]["x"])) == pytest.approx([-6.2, -6.1], abs=0.06)

    fig = _RegionalDummy()
    render_regional_basemap(_tiny_geo_grid(), fig=fig, borders=False, rivers=False,
                            scale_bar=False, show_north_arrow=False, frame=False)
    names = [n for n, _ in fig.calls]
    assert names == ["basemap", "coast", "plot"]
    assert "borders" not in fig.calls[1][1] and "rivers" not in fig.calls[1][1]
    assert fig.calls[0][1]["frame"] == ["wsne", "f0"]


# Real GMT drawing from here on.
pygmt = pytest.importorskip("pygmt")


def _utm_raster():
    da = xr.DataArray(np.random.default_rng(0).normal(1000, 50, (40, 40)), dims=("y", "x"),
                      coords={"y": 5040000.0 - np.arange(40) * 500.0,
                              "x": 600000.0 + np.arange(40) * 500.0})
    return da.rio.write_crs("EPSG:32611")


@pytest.mark.parametrize("mode", ["global", "regional", "aerial"])
def test_locator_modes_export_png(tmp_path, mode):
    """Each mode exports a PNG (global mode used to fail in Ghostscript)."""
    fig = render_orthographic_basemap(_utm_raster(), mode=mode,
                                      style=MapStyle(figure_size="8c", transparent=False))
    out = tmp_path / f"locator_{mode}.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


@pytest.mark.network
def test_regional_relief_exports_png(tmp_path):
    fig = render_regional_basemap(_utm_raster(), relief=True,
                                  style=MapStyle(figure_size="8c", transparent=False))
    out = tmp_path / "regional_relief.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


@pytest.mark.parametrize("mode,position", [("global", "TR"), ("regional", "BL")])
def test_locator_inset_on_a_map(tmp_path, mode, position):
    from lsdviztools3.render.hillshade import render_hillshade

    da = _utm_raster()
    fig = render_hillshade(da, coordinates="local", shade=True,
                           style=MapStyle(cmap="gray", transparent=False))
    add_locator_inset(fig, da, mode=mode, position=position, width="3c")
    out = tmp_path / f"inset_{mode}.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_locator_inset_refuses_bad_arguments():
    with pytest.raises(ValueError, match="global' or 'regional"):
        add_locator_inset(None, _utm_raster(), mode="aerial")
    with pytest.raises(ValueError, match="position"):
        add_locator_inset(None, _utm_raster(), position="middle")


def test_basemap_and_hillshade_commands(tmp_path):
    from lsdviztools3.cli import _common, plot_basemap, plot_hillshade

    dem = tmp_path / "dem.tif"
    _utm_raster().rio.to_raster(dem)
    out = tmp_path / "regional.png"
    assert plot_basemap.main(["--in", str(dem), "--out", str(out), "--mode", "regional",
                              "--extent-km", "120", "--no-rivers"]) == 0
    assert out.exists()
    out2 = tmp_path / "hs.png"
    assert plot_hillshade.main(["--in", str(dem), "--out", str(out2), "--locator", "global",
                                "--locator-position", "BR", "--locator-width", "3c"]) == 0
    assert out2.exists()
    args = plot_hillshade.build_parser().parse_args(["--in", "x", "--out", "y"])
    assert _common.locator_kwargs(args) is None
