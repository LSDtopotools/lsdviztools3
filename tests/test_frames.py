"""Maps of rasters in lsdtt3 oblique frames.

An lsdtt3 frame is a Hotine Oblique Mercator CRS whose central line runs
along a feature (+x) through the frame centre (0, 0). These tests check the
frame helpers against pyproj, the reference-grid graticule, the north arrow
at the convergence of its own position, the axis and CRS labels, that
``coordinates="local"`` keeps a frame instead of warping it to UTM, and an
end-to-end render of a synthetic frame DEM.

Expected convergences (true north, clockwise from grid up) are pyproj's
``meridian_convergence`` with the sign flipped: oblique frame at azimuth 30
gives 60 at its centre; UTM 33N at 45 N on the east zone edge -2.1223 (its
magnitude, 2.12, is the plan's figure); EPSG:3413 at 0 E, 80 N -45.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from pyproj import CRS, Proj, Transformer

from lsdviztools3 import frames
from lsdviztools3.render.base import crs_label_text, display_frame
from lsdviztools3.render.graticule import (
    add_reference_graticule,
    graticule_label,
    graticule_lines,
    nice_interval,
    placed_labels,
)
from lsdviztools3.render.style import MapStyle


def _frame(alpha=30.0, lat0=38.4, lonc=-110.0):
    return CRS.from_proj4(
        f"+proj=omerc +lat_0={lat0} +lonc={lonc} +alpha={alpha} +gamma=90 +k=1 "
        "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")


class _RecordingFig:
    def __init__(self):
        self.plots, self.texts = [], []

    def plot(self, **kwargs):
        self.plots.append(kwargs)

    def text(self, **kwargs):
        self.texts.append(kwargs)


# --- frame helpers ------------------------------------------------------------

def test_frame_detection_and_parameters():
    f = _frame()
    assert frames.is_oblique_frame(f)
    assert frames.frame_azimuth(f) == pytest.approx(30.0)
    assert frames.frame_centre(f) == pytest.approx((-110.0, 38.4))
    for other in (CRS.from_epsg(32612), CRS.from_epsg(4326), CRS.from_epsg(3413)):
        assert not frames.is_oblique_frame(other)
        assert frames.frame_azimuth(other) is None


def test_reference_crs_choice():
    assert frames.reference_crs_for(_frame()).to_epsg() == 32612
    assert frames.reference_crs_for(_frame(lat0=-33.9, lonc=18.4)).to_epsg() == 32734
    assert frames.reference_crs_for(_frame(lat0=80.0, lonc=-40.0)).to_epsg() == 3413
    assert frames.reference_crs_for(_frame(lat0=-78.0, lonc=160.0)).to_epsg() == 3976
    assert frames.reference_crs_for(_frame(), "EPSG:26912").to_epsg() == 26912
    assert frames.reference_crs_for(CRS.from_epsg(32612)) is None


@pytest.mark.parametrize("crs, xy, expected", [
    (_frame(), (0.0, 0.0), 60.0),
    (CRS.from_epsg(32633), (736446.026, 4987329.505), -2.1223),
    (CRS.from_epsg(32633), (500000.0, 4982950.400), 0.0),
    (CRS.from_epsg(3413), (767861.606, -767861.606), -45.0),
])
def test_grid_convergence_matches_pyproj(crs, xy, expected):
    assert frames.grid_convergence_degrees(crs, *xy) == pytest.approx(expected, abs=1e-3)
    lon, lat = Transformer.from_crs(crs, 4326, always_xy=True).transform(*xy)
    pyproj_mc = Proj(crs).get_factors(lon, lat).meridian_convergence
    assert frames.grid_convergence_degrees(crs, *xy) == pytest.approx(-pyproj_mc, abs=1e-3)


# --- graticule ------------------------------------------------------------------

def test_nice_interval():
    assert nice_interval(37_000.0) == 10_000.0
    assert nice_interval(9_000.0) == 5_000.0
    assert nice_interval(4_400.0) == 2_000.0


def test_graticule_lines_are_constant_reference_coordinates():
    f = _frame()
    extent = [-20_000.0, 20_000.0, -5_000.0, 5_000.0]
    lines = graticule_lines(extent, f)
    to_ref = Transformer.from_crs(f, CRS.from_epsg(32612), always_xy=True)
    assert {line["axis"] for line in lines} == {"E", "N"}
    for line in lines:
        e, n = to_ref.transform(line["x"], line["y"])
        held = e if line["axis"] == "E" else n
        assert np.max(np.abs(held - line["value"])) < 1e-3          # metres
        assert np.all((line["x"] >= extent[0] - 1e-6) & (line["x"] <= extent[1] + 1e-6))
        assert np.all((line["y"] >= extent[2] - 1e-6) & (line["y"] <= extent[3] + 1e-6))
        assert line["value"] % 5_000.0 == 0.0
    labelled = [line for line in lines if line["label"]]
    assert labelled
    for line in labelled:
        x, y, edge = line["label"]
        on = {"bottom": y == pytest.approx(extent[2]), "top": y == pytest.approx(extent[3]),
              "left": x == pytest.approx(extent[0]), "right": x == pytest.approx(extent[1])}
        assert on[edge]
    assert graticule_label("E", 512_000.0) == "512 km E"


def test_graticule_only_over_frames_and_toggle():
    extent_km = [-20.0, 20.0, -5.0, 5.0]
    fig = _RecordingFig()
    add_reference_graticule(fig, extent_km, _frame())
    assert fig.plots and fig.texts
    assert all(t["justify"] in ("CB", "CT", "LM", "RM") for t in fig.texts)
    utm = _RecordingFig()
    add_reference_graticule(utm, [600.0, 626.0, 4300.0, 4322.0], CRS.from_epsg(32617))
    assert not utm.plots and not utm.texts
    off = _RecordingFig()
    add_reference_graticule(off, extent_km, _frame(), MapStyle(graticule=False))
    assert not off.plots
    fixed = _RecordingFig()
    add_reference_graticule(fixed, extent_km, _frame(), MapStyle(graticule_interval_km=2.0))
    assert len(fixed.plots) > len(fig.plots)


def test_labels_on_one_edge_do_not_overlap():
    extent = [-62_000.0, 62_000.0, -20_000.0, 20_000.0]
    lines = [{"axis": "E", "value": 600e3, "label": (-21_000.0, -20_000.0, "bottom")},
             {"axis": "N", "value": 4750e3, "label": (-20_500.0, -20_000.0, "bottom")},
             {"axis": "N", "value": 4800e3, "label": (40_000.0, -20_000.0, "bottom")},
             {"axis": "N", "value": 4700e3, "label": (-62_000.0, 10_000.0, "left")}]
    kept = [line["value"] for line, _ in placed_labels(lines, extent)]
    assert kept == [600e3, 4800e3, 4700e3]


# --- labels ---------------------------------------------------------------------

def test_frame_axis_labels_and_crs_title():
    frame = display_frame(_frame(), MapStyle())
    assert "xaf+lAlong-feature (km)" in frame and "yaf+lAcross-feature (km)" in frame
    assert "xaf+lEasting (km)" in display_frame(CRS.from_epsg(32612), MapStyle())
    assert crs_label_text(_frame()) == "Oblique Mercator, azimuth 30.0@."
    assert crs_label_text(CRS.from_epsg(32612)) == "CRS is EPSG:32612"


# --- north arrow ----------------------------------------------------------------

def _arrow_bearing(fig):
    xs, ys = fig.plots[0]["x"], fig.plots[0]["y"]
    apex = (xs[1], ys[1])
    base = (0.5 * (xs[0] + xs[2]), 0.5 * (ys[0] + ys[2]))
    return math.degrees(math.atan2(apex[0] - base[0], apex[1] - base[1]))


def test_north_arrow_in_an_oblique_frame_points_60_degrees():
    from lsdviztools3.render.decorations import add_north_arrow

    fig = _RecordingFig()
    extent = [-2.0, 2.0, -2.0, 2.0]       # km, around the frame centre
    add_north_arrow(fig, extent, _frame())
    xs, ys = fig.plots[0]["x"], fig.plots[0]["y"]
    cx, cy = 1000.0 * sum(xs[:3]) / 3.0, 1000.0 * sum(ys[:3]) / 3.0
    expected = frames.grid_convergence_degrees(_frame(), cx, cy)
    assert expected == pytest.approx(60.0, abs=0.1)
    assert _arrow_bearing(fig) == pytest.approx(expected, abs=0.05)
    assert fig.texts[0]["text"] == "N"


def test_north_arrow_utm_tilts_slightly_and_polar_turns():
    from lsdviztools3.render.decorations import add_north_arrow

    edge = _RecordingFig()
    add_north_arrow(edge, [720.0, 750.0, 4975.0, 5000.0], CRS.from_epsg(32633))
    assert _arrow_bearing(edge) == pytest.approx(-2.1, abs=0.15)
    # Polar stereographic convergence changes fast: compare at the arrow itself
    # (top-right corner of a map around 0 E, 80 N, where it is -45).
    polar = _RecordingFig()
    add_north_arrow(polar, [740.0, 800.0, -800.0, -740.0], CRS.from_epsg(3413))
    xs, ys = polar.plots[0]["x"], polar.plots[0]["y"]
    at_arrow = frames.grid_convergence_degrees(
        CRS.from_epsg(3413), 1000.0 * sum(xs[:3]) / 3.0, 1000.0 * sum(ys[:3]) / 3.0)
    assert at_arrow == pytest.approx(-45.0, abs=3.0)
    assert _arrow_bearing(polar) == pytest.approx(at_arrow, abs=0.05)
    grid = _RecordingFig()
    add_north_arrow(grid, [720.0, 750.0, 4975.0, 5000.0], CRS.from_epsg(32633),
                    true_north=False)
    assert _arrow_bearing(grid) == pytest.approx(0.0, abs=1e-9)


# --- coordinates="local" keeps the frame -----------------------------------------

def _frame_dataarray(nx=60, ny=20, cell=100.0):
    import rioxarray  # noqa: F401
    import xarray as xr

    x = -nx * cell / 2 + (np.arange(nx) + 0.5) * cell
    y = ny * cell / 2 - (np.arange(ny) + 0.5) * cell
    X, Y = np.meshgrid(x, y)
    z = 500.0 + 0.05 * X + 20.0 * np.sin(X / 700.0) * np.cos(Y / 400.0)
    da = xr.DataArray(z.astype("float32"), dims=("y", "x"), coords={"x": x, "y": y})
    return da.rio.write_crs(_frame().to_wkt()).rio.write_nodata(np.nan)


def test_local_coordinates_keep_a_frame_raster_and_vector():
    import geopandas as gpd
    from shapely.geometry import LineString

    from lsdviztools3.io.raster import resolve_target_crs, to_coordinates
    from lsdviztools3.io.vector import to_coordinates as vector_to_coordinates

    da = _frame_dataarray()
    assert frames.is_oblique_frame(resolve_target_crs(da, "local"))
    same = to_coordinates(da, "local")
    assert same.shape == da.shape and np.array_equal(same.x.values, da.x.values)
    gdf = gpd.GeoDataFrame(geometry=[LineString([(-1000, 0), (1000, 0)])], crs=_frame().to_wkt())
    kept = vector_to_coordinates(gdf, "local")
    assert frames.is_oblique_frame(kept.crs)
    wgs = gdf.to_crs(4326)
    onto_frame = vector_to_coordinates(wgs, resolve_target_crs(da, "local"))
    assert np.allclose(np.asarray(onto_frame.geometry[0].coords),
                       np.asarray(gdf.geometry[0].coords), atol=1e-6)


def test_render_a_frame_dem(tmp_path):
    pytest.importorskip("pygmt")
    from lsdviztools3.render.hillshade import render_hillshade

    path = tmp_path / "frame_dem.tif"
    _frame_dataarray().rio.to_raster(path)
    out = tmp_path / "frame_map.png"
    fig = render_hillshade(str(path), coordinates="local", shade=True,
                           style=MapStyle(cmap="gray", transparent=False))
    fig.savefig(out)
    assert out.stat().st_size > 1000
