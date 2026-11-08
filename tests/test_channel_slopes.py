"""Smoothed channel profiles and slopes (lsdtt3 _channel_slope_points)."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from lsdviztools3.charts.channel_slopes import (
    load_channel_points,
    plot_channel_profiles,
    plot_slope_area,
    plot_slope_distance,
)
from lsdviztools3.cli import plot_channel_slopes


def _points(n=40):
    """One straight channel: 10 m cells, S = 0.02, area growing downstream.
    The first two cells have no valid fit (slope_valid 0, NaN)."""
    d = np.arange(n)[::-1] * 10.0                         # flow distance to the outlet
    z = 1.0 + 0.02 * d
    area = 1e5 * (n - np.arange(n))[::-1] + 1e5
    slope = np.full(n, 0.02)
    zs = z.copy()
    valid = np.ones(n, dtype=int)
    valid[:2] = 0
    slope[:2] = np.nan
    zs[:2] = np.nan
    return pd.DataFrame({"latitude": 56.0 + np.arange(n) * 1e-4, "longitude": -3.0 + np.arange(n) * 1e-4,
                         "flat": np.arange(n), "drainage_area_m2": area, "flow_distance_m": d,
                         "elevation_m": z, "smoothed_elevation_m": zs, "slope_m_per_m": slope,
                         "slope_valid": valid, "slope_window_length_m": 50.0})


def _field():
    return pd.DataFrame({"latitude": [56.001, 56.002], "longitude": [-2.999, -2.998], "site_id": ["a", "b"],
                         "d50_mm": [12.0, 80.0], "drainage_area_m2": [1e6, 2e6], "flow_distance_m": [100.0, 50.0],
                         "smoothed_elevation_m": [3.0, 2.0], "slope_m_per_m": [0.02, 0.03], "slope_valid": [1, 1]})


@pytest.fixture
def files(tmp_path):
    sp = tmp_path / "demo_channel_slope_points.csv"
    fp = tmp_path / "field.csv"
    _points().to_csv(sp, index=False)
    _field().to_csv(fp, index=False)
    return sp, fp


def test_load_csv_is_wgs84_points(files):
    g = load_channel_points(files[0])
    assert g.crs.to_epsg() == 4326
    assert len(g) == 40
    assert g.geometry.iloc[0].x == pytest.approx(-3.0)


def test_load_gpkg_round_trip(tmp_path, files):
    g = load_channel_points(files[0])
    p = tmp_path / "pts.gpkg"
    g.to_file(p, driver="GPKG")
    back = load_channel_points(p)
    assert len(back) == 40 and "slope_m_per_m" in back.columns


def test_load_refuses_missing_columns_and_files(tmp_path):
    p = tmp_path / "bad.csv"
    pd.DataFrame({"latitude": [1.0], "longitude": [2.0], "x": [3]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="drainage_area_m2"):
        load_channel_points(p)
    q = tmp_path / "nolatlon.csv"
    pd.DataFrame({"drainage_area_m2": [1.0], "slope_m_per_m": [0.1]}).to_csv(q, index=False)
    with pytest.raises(ValueError, match="latitude"):
        load_channel_points(q)
    with pytest.raises(FileNotFoundError):
        load_channel_points(tmp_path / "missing.csv")


def test_profile_draws_raw_and_only_valid_smoothed(files):
    pts = load_channel_points(files[0])
    fig, ax = plot_channel_profiles(pts)
    raw, smooth = ax.collections[:2]
    assert len(raw.get_offsets()) == 40
    assert len(smooth.get_offsets()) == 38                # the two invalid fits are left out
    # x in km
    assert raw.get_offsets()[:, 0].max() == pytest.approx(0.39)
    plt.close(fig)


def test_slope_plots_log_axes_and_valid_only(files):
    pts = load_channel_points(files[0])
    fig, ax = plot_slope_area(pts)
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"
    assert len(ax.collections[0].get_offsets()) == 38
    plt.close(fig)
    fig, ax = plot_slope_distance(pts)
    assert ax.get_yscale() == "log" and ax.get_xscale() == "linear"
    plt.close(fig)


def test_min_area_filter(files):
    pts = load_channel_points(files[0])
    cut = float(np.median(pts["drainage_area_m2"]))
    fig, ax = plot_slope_area(pts, min_area_m2=cut)
    xs = ax.collections[0].get_offsets()[:, 0]
    assert len(xs) > 0 and xs.min() >= cut
    plt.close(fig)


def test_field_points_single_colour_and_coloured(files):
    pts, field = load_channel_points(files[0]), load_channel_points(files[1])
    fig, ax = plot_slope_area(pts, field=field)
    assert len(ax.collections) == 2 and len(fig.axes) == 1   # no colourbar
    assert len(ax.collections[1].get_offsets()) == 2
    plt.close(fig)
    fig, ax = plot_channel_profiles(pts, field=field, field_color_by="d50_mm")
    assert len(fig.axes) == 2                                 # a colourbar
    assert np.allclose(sorted(ax.collections[-1].get_array()), [12.0, 80.0])
    plt.close(fig)
    with pytest.raises(ValueError, match="no column"):
        plot_slope_area(pts, field=field, field_color_by="nope")


def test_legend_handles_are_fixed_size(files):
    pts, field = load_channel_points(files[0]), load_channel_points(files[1])
    fig, ax = plot_slope_area(pts, field=field, field_color_by="d50_mm")
    sizes = {t.get_text(): h.get_markersize() for t, h in zip(ax.get_legend().get_texts(),
                                                             ax.get_legend().legend_handles)}
    assert sizes["field points"] == 7 and sizes["channel cells"] == 5
    plt.close(fig)


def test_cli_writes_three_pngs(tmp_path, files):
    out = tmp_path / "out"
    assert plot_channel_slopes.main(["--slope-points", str(files[0]), "--field-points", str(files[1]),
                                     "--field-color-by", "d50_mm", "--out-dir", str(out), "--dpi", "60"]) == 0
    names = sorted(p.name for p in out.glob("*.png"))
    assert names == ["demo_channel_profiles.png", "demo_slope_area.png", "demo_slope_distance.png"]


def test_cli_refusals(tmp_path, files):
    with pytest.raises(SystemExit, match="needs --field-points"):
        plot_channel_slopes.main(["--slope-points", str(files[0]), "--field-color-by", "d50_mm",
                                  "--out-dir", str(tmp_path)])
    with pytest.raises(SystemExit, match="no column"):
        plot_channel_slopes.main(["--slope-points", str(files[0]), "--field-points", str(files[1]),
                                  "--field-color-by", "nope", "--out-dir", str(tmp_path)])
    with pytest.raises(SystemExit):
        plot_channel_slopes.main(["--slope-points", str(tmp_path / "missing.csv"), "--out-dir", str(tmp_path)])


def test_map_reprojects_wgs84_points_onto_a_utm_dem(tmp_path):
    """CSV points are WGS84; the DEM is UTM 33N. The map must put them in the
    DEM's CRS before drawing (else nothing lands on the hillshade)."""
    import geopandas as gpd
    import rasterio
    from rasterio.transform import from_origin
    from lsdviztools3.charts.channel_slopes import render_channel_slope_map
    dem = tmp_path / "dem.tif"
    with rasterio.open(dem, "w", driver="GTiff", height=30, width=80, count=1, dtype="float32",
                       crs="EPSG:32633", transform=from_origin(499990, 4200100, 10, 10)) as d:
        d.write(np.add.outer(np.arange(30.0), np.arange(80.0)).astype("float32"), 1)
    utm = gpd.GeoDataFrame(_points(), geometry=gpd.points_from_xy(500000 + np.arange(40) * 15.0,
                                                                   np.full(40, 4199950.0)), crs="EPSG:32633")
    wgs = utm.to_crs("EPSG:4326")
    field = wgs.iloc[[10, 20]]
    out = render_channel_slope_map(dem, wgs, tmp_path / "map.png", field=field, slope_vmin=0.01, slope_vmax=0.1)
    assert out.exists() and out.stat().st_size > 1000
