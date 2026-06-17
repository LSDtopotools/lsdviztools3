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
