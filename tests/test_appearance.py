"""Unit tests for publication-appearance behaviour (no network / no GMT)."""

import numpy as np
import xarray as xr
from pyproj import CRS

from lsdviztools3.render import base
from lsdviztools3.render.style import MapStyle


def _pt(s):
    """Parse a GMT font size like '9p' to an int point size."""
    return int(s.rstrip("p"))


def test_style_font_and_resolution_defaults():
    s = MapStyle()
    assert _pt(s.font_annotation) == 9          # tick labels
    assert _pt(s.font_label) == 11              # axis labels
    assert _pt(s.font_annotation) >= 8          # 8 pt floor
    assert _pt(s.font_title) >= 8
    assert s.dpi == 300 and s.high_res_dpi == 600
    assert s.transparent is True
    assert MapStyle.SINGLE_COLUMN == "8.5c"
    assert MapStyle.DOUBLE_COLUMN == "17c"


def test_display_frame_projected_has_crs_title_and_km():
    frame = base.display_frame(CRS.from_epsg(32629), MapStyle())
    assert frame[0] == "WSne"
    assert any("Easting (km)" in part for part in frame)
    assert any("Northing (km)" in part for part in frame)


def test_display_frame_geographic_has_no_crs_title():
    frame = base.display_frame(CRS.from_epsg(4326), MapStyle())
    assert not any(part.startswith("WSne+t") for part in frame)
    assert any("Longitude" in part for part in frame)


def test_grid_to_display_scales_projected_to_km():
    da = xr.DataArray(
        np.zeros((2, 2)),
        dims=("y", "x"),
        coords={"x": [600000.0, 600030.0], "y": [6300000.0, 6300030.0]},
    )
    grid, extent = base.grid_to_display(da, CRS.from_epsg(32629))
    assert grid.x.values[0] == 600.0  # metres -> km
    assert extent == [600.0, 600.03, 6300.0, 6300.03]


def test_grid_to_display_leaves_geographic_in_degrees():
    da = xr.DataArray(
        np.zeros((2, 2)),
        dims=("y", "x"),
        coords={"x": [-6.2, -6.1], "y": [56.8, 56.9]},
    )
    grid, extent = base.grid_to_display(da, CRS.from_epsg(4326))
    assert grid.x.values[0] == -6.2  # unchanged
    assert extent == [-6.2, -6.1, 56.8, 56.9]
