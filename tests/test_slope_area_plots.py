"""lsdtt-slope-area outputs: readers, charts and the lsdviz3-slope-area CLI."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from lsdviztools3.charts.slope_area import (
    plot_binned_slope_area,
    plot_low_concavity_slope_area,
    plot_slope_area_method_comparison,
    plot_theta_distribution,
    stream_order_pen_width,
)
from lsdviztools3.io import slope_area as sa_io

BIN_COLUMNS = ["drainage_area", "d05_slope", "d25_slope", "d50_slope", "d75_slope", "d95_slope",
               "number_of_pixels"]


def _binned(scale=1.0, n=6):
    a = np.logspace(3, 7, n)
    d50 = 0.2 * scale * (a / 1e3) ** -0.1
    return pd.DataFrame({"drainage_area": a[::-1], "d05_slope": d50 * 0.3, "d25_slope": d50 * 0.7,
                         "d50_slope": d50, "d75_slope": d50 * 1.3, "d95_slope": d50 * 2.0,
                         "number_of_pixels": np.arange(n, 0, -1) * 10})[BIN_COLUMNS]


def _points():
    """10 channel pixels: 3 low concavity, 4 not, 3 no fit (one with nan slope)."""
    flag = np.array([1, 1, 1, 0, 0, 0, 0, -1, -1, -1])
    theta = np.array([0.05, -0.1, 0.12, 0.5, -0.6, 0.9, 30.0, np.nan, np.nan, np.nan])
    slope = np.array([0.1, 0.11, 0.12, 0.2, 0.05, 0.08, 0.3, 0.2, np.nan, 0.15])
    n = flag.size
    return pd.DataFrame({"flat": np.arange(n), "x": 500000.0 + 10 * np.arange(n),
                         "y": np.full(n, 4200050.0), "drainage_area_m2": np.logspace(4, 6, n),
                         "slope": slope, "theta": theta, "r2": np.full(n, 0.5),
                         "n_pixels": np.full(n, 10), "log10_area_range": np.full(n, 0.3),
                         "n_junctions_in_window": np.zeros(n, dtype=int), "low_concavity": flag})


def _run_dir(tmp_path, name="run", *, lambdas=("90", "7_5"), windows=(), low_concavity=False,
             points=True):
    d = tmp_path / name
    d.mkdir()
    for k, lam in enumerate(lambdas):
        _binned(1 + k).to_csv(d / f"sa_lambda_{lam}m.csv", index=False)
    for k, w in enumerate(windows):
        _binned(0.8 + 0.1 * k).to_csv(d / f"sa_d8window_{w}m.csv", index=False)
    # files that look similar but are not binned tables
    pd.DataFrame({"polyfit_lambda": [90], "drainage_area_with_max_d50_slope": [1e4]}).to_csv(
        d / "sa_maximum_d50_slope_area.csv", index=False)
    (d / "sa_ingested_params.param").write_text(
        "# lsdtt3 parameter file\nwrite_fname: sa\n\n# Default parameters used\n"
        "low_concavity_max_abs_concavity: 0.2\nbasin_mask_fname: \n", encoding="utf-8")
    if low_concavity:
        if points:
            _points().to_csv(d / "sa_low_concavity_points.csv", index=False, na_rep="nan")
        _segments().to_file(d / "sa_low_concavity_segments.fgb", driver="FlatGeobuf")
    return d


def _segments():
    import geopandas as gpd
    from shapely.geometry import LineString
    rows = {"segment_id": [0, 1, 2], "link_id": [0, 0, 1], "low_concavity": [1, 0, -1],
            "n_pixels": [3, 4, 3], "length_m": [30.0, 40.0, 30.0], "median_theta": [0.05, 0.6, np.nan],
            "median_slope": [0.1, 0.1, np.nan], "min_area_m2": [1e4, 2e4, 3e4],
            "max_area_m2": [2e4, 3e4, 4e4], "stream_order": [1, 2, 3], "upstream_flat": [0, 3, 7],
            "downstream_flat": [2, 6, 9]}
    geoms = [LineString([(500000, 4200050), (500200, 4200060)]),
             LineString([(500200, 4200060), (500400, 4200040)]),
             LineString([(500400, 4200040), (500700, 4200020)])]
    return gpd.GeoDataFrame(rows, geometry=geoms, crs="EPSG:32633")


def _dem(path):
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(path, "w", driver="GTiff", height=30, width=80, count=1, dtype="float32",
                       crs="EPSG:32633", transform=from_origin(499990, 4200100, 10, 10)) as d:
        d.write(np.add.outer(np.arange(30.0), np.arange(80.0)).astype("float32"), 1)
    return path


# -----------------------------------------------------------------------------
# Readers
# -----------------------------------------------------------------------------
def test_ridgeline_reexports_the_reader():
    from lsdviztools3.io import ridgeline
    assert ridgeline.load_slope_area is sa_io.load_slope_area
    assert ridgeline.SLOPE_AREA_COLUMNS == sa_io.SLOPE_AREA_COLUMNS


def test_load_slope_area_sorts_and_validates(tmp_path):
    p = tmp_path / "t.csv"
    _binned().to_csv(p, index=False)
    df = sa_io.load_slope_area(p)
    assert np.all(np.diff(df["drainage_area"]) > 0)
    _binned().drop(columns="d50_slope").to_csv(p, index=False)
    with pytest.raises(ValueError, match="d50_slope"):
        sa_io.load_slope_area(p)


def test_find_binned_slope_area(tmp_path):
    d = _run_dir(tmp_path, lambdas=("90", "7_5"), windows=("300", "150"))
    found = sa_io.find_binned_slope_area(d, "sa")
    assert [(t["method"], t["scale_m"]) for t in found] == [
        ("polyfit", 7.5), ("polyfit", 90.0), ("d8_flow_path", 150.0), ("d8_flow_path", 300.0)]
    assert found[0]["label"] == "polyfit, λ 7.5 m"
    assert found[2]["label"] == "D8 flow path, window 150 m"
    assert len(found[0]["table"]) == 6
    assert [t["scale_m"] for t in sa_io.find_binned_slope_area(d, "sa", method="d8_flow_path")] == [150, 300]
    assert sa_io.find_binned_slope_area(d, "other") == []


def test_low_concavity_readers(tmp_path):
    d = _run_dir(tmp_path, low_concavity=True)
    files = sa_io.low_concavity_files(d, "sa")
    pts = sa_io.load_low_concavity_points(files["points"])
    assert len(pts) == 10 and pts["theta"].isna().sum() == 3
    assert pts["low_concavity"].tolist().count(-1) == 3
    seg = sa_io.load_low_concavity_segments(files["segments"])
    assert sorted(seg["low_concavity"]) == [-1, 0, 1]
    params = sa_io.load_ingested_params(d / "sa_ingested_params.param")
    assert params["write_fname"] == "sa" and params["basin_mask_fname"] == ""
    assert sa_io.param_float(params, "low_concavity_max_abs_concavity") == 0.2
    assert sa_io.param_float(params, "missing", 0.15) == 0.15


def test_low_concavity_rasters(tmp_path):
    import rasterio
    from rasterio.transform import from_origin
    flag = np.array([[1, 0], [-1, -9999]], dtype="float32")
    theta = np.array([[0.1, 0.5], [-9999, -9999]], dtype="float32")
    for name, arr in (("f.tif", flag), ("t.tif", theta)):
        with rasterio.open(tmp_path / name, "w", driver="GTiff", height=2, width=2, count=1,
                           dtype="float32", nodata=-9999, crs="EPSG:32633",
                           transform=from_origin(0, 20, 10, 10)) as dst:
            dst.write(arr, 1)
    pix = sa_io.load_low_concavity_rasters(tmp_path / "t.tif", tmp_path / "f.tif")
    assert sorted(pix["low_concavity"]) == [-1, 0, 1]
    assert pix["theta"].isna().sum() == 1


# -----------------------------------------------------------------------------
# Charts
# -----------------------------------------------------------------------------
def _tables(labels=("a", "b")):
    return [{"label": lab, "table": _binned(1 + k)} for k, lab in enumerate(labels)]


def test_binned_chart():
    fig, ax = plot_binned_slope_area(_tables())
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"
    assert len(ax.get_lines()) == 2 and len(ax.collections) == 2     # d50 lines, d25-d75 bands
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["a", "b"]
    plt.close(fig)
    with pytest.raises(ValueError):
        plot_binned_slope_area([])


def test_method_comparison_chart():
    fig, ax = plot_slope_area_method_comparison(_tables(("p",)), _tables(("d1", "d2")),
                                                name="poly", compare_name="d8")
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"
    lines = ax.get_lines()
    assert len(lines) == 3 and len(ax.collections) == 0
    assert [ln.get_linestyle() for ln in lines] == ["-", "--", "--"]
    assert len({ln.get_color() for ln in lines}) == 3
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["poly: p", "d8: d1", "d8: d2"]
    plt.close(fig)
    fig, ax = plot_slope_area_method_comparison(_tables(("p",)), _tables(("d",)), bands=True)
    assert len(ax.collections) == 2
    plt.close(fig)


def test_theta_distribution_chart():
    pts = _points()
    fig, ax = plot_theta_distribution(pts["theta"], pts["low_concavity"], 0.15)
    xs = sorted(ln.get_xdata()[0] for ln in ax.get_lines())
    assert xs == pytest.approx([-0.15, 0.15])
    texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "low concavity (3)" in texts and "not low concavity (4)" in texts
    assert "no valid fit, not drawn (3)" in texts
    assert not any("outside" in t for t in texts)         # a small sample is shown whole
    assert "3 of 7" in ax.get_title()
    assert sum(p.get_height() for p in ax.patches) == pytest.approx(7)
    plt.close(fig)
    fig, ax = plot_theta_distribution(pts["theta"], pts["low_concavity"], 0.15, theta_range=(-1, 1))
    assert ax.get_xlim() == (-1, 1)
    assert "valid, outside the x range (1)" in [t.get_text() for t in ax.get_legend().get_texts()]
    assert sum(p.get_height() for p in ax.patches) == pytest.approx(6)    # theta = 30 not drawn
    plt.close(fig)


def test_low_concavity_slope_area_chart():
    fig, ax = plot_low_concavity_slope_area(_points())
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"
    sizes = [c.get_offsets().shape[0] for c in ax.collections]
    assert sizes == [2, 4, 3]                             # no fit (one nan slope), not, low
    assert "1 without a positive gradient" in ax.get_title()
    plt.close(fig)


def test_stream_order_pen_width():
    assert stream_order_pen_width(1) < stream_order_pen_width(2) < stream_order_pen_width(3)
    assert stream_order_pen_width(50) == 3.2


def test_low_concavity_map_writes_png(tmp_path):
    from lsdviztools3.charts.slope_area import render_low_concavity_map
    out = render_low_concavity_map(_dem(tmp_path / "dem.tif"), _segments(), tmp_path / "m.png")
    assert out.exists() and out.stat().st_size > 1000


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def test_cli_end_to_end(tmp_path, capsys):
    from lsdviztools3.cli import plot_slope_area
    d = _run_dir(tmp_path, "polyfit", lambdas=("90",), low_concavity=True)
    c = _run_dir(tmp_path, "d8", lambdas=(), windows=("150",))
    out = tmp_path / "out"
    assert plot_slope_area.main(["--data-dir", str(d), "--prefix", "sa", "--out-dir", str(out),
                                 "--out-prefix", "x", "--compare-data-dir", str(c),
                                 "--dem", str(_dem(tmp_path / "dem.tif")), "--dpi", "50"]) == 0
    names = sorted(p.name for p in out.iterdir())
    assert names == ["x_binned.png", "x_low_concavity_map.png", "x_low_concavity_slope_area.png",
                     "x_method_comparison.png", "x_theta_distribution.png"]
    printed = capsys.readouterr().out
    assert "low_concavity_max_abs_concavity: 0.2" in printed
    for n in names:
        assert str(out / n) in printed


def test_cli_without_points_or_dem(tmp_path, capsys):
    from lsdviztools3.cli import plot_slope_area
    d = _run_dir(tmp_path, low_concavity=True, points=False)
    out = tmp_path / "out"
    assert plot_slope_area.main(["--data-dir", str(d), "--prefix", "sa", "--out-dir", str(out),
                                 "--dpi", "50"]) == 0
    assert sorted(p.name for p in out.iterdir()) == ["sa_binned.png"]
    printed = capsys.readouterr().out
    assert "sa_low_concavity_points.csv not found" in printed
    assert "give --dem" in printed


def test_cli_no_outputs(tmp_path):
    from lsdviztools3.cli import plot_slope_area
    (tmp_path / "empty").mkdir()
    with pytest.raises(SystemExit):
        plot_slope_area.main(["--data-dir", str(tmp_path / "empty"), "--prefix", "sa",
                              "--out-dir", str(tmp_path / "o")])
