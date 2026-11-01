"""Perimeter valley sources map and profile (lsdtt3 run_perimeter_sources)."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import LineString, Point

from lsdviztools3.charts.perimeter_sources import (
    _ingested_value,
    plot_perimeter_profile,
    plot_perimeter_sources_map,
)
from lsdviztools3.cli import plot_perimeter_sources

CRS = "EPSG:32632"


def _write_run(tmp_path):
    """A 100 x 100 cell (10 m) square: a valley enters the north edge and
    leaves the south edge; a narrow notch on the west edge is rejected."""
    stem = tmp_path / "p"
    n = 100
    # Perimeter clockwise from the north-west corner, cell centres.
    xs, ys = [], []
    for c in range(n):
        xs.append(5.0 + 10 * c); ys.append(995.0)
    for r in range(1, n):
        xs.append(995.0); ys.append(995.0 - 10 * r)
    for c in range(n - 2, -1, -1):
        xs.append(5.0 + 10 * c); ys.append(5.0)
    for r in range(n - 2, 0, -1):
        xs.append(5.0); ys.append(995.0 - 10 * r)
    xs, ys = np.array(xs), np.array(ys)
    m = len(xs)
    elev = 100.0 + 0.1 * np.abs(xs - 500.0)
    run_id = np.full(m, -1)
    north = (ys == 995.0) & (np.abs(xs - 500.0) < 100)
    south = (ys == 5.0) & (np.abs(xs - 500.0) < 100)
    west = (xs == 5.0) & (np.abs(ys - 500.0) < 15)
    run_id[north], run_id[south], run_id[west] = 0, 1, 2
    elev[south] -= 20.0
    pd.DataFrame({
        "visit": np.arange(m), "s_m": 10.0 * np.arange(m),
        "row": 0, "col": 0, "x": xs, "y": ys, "elevation": elev,
        "gradient": np.where(run_id >= 0, 0.01, 0.2), "gentle": (run_id >= 0).astype(int),
        "run_id": run_id, "outlet_cell": south.astype(int),
        "outlet_run_id": np.where(south, 0, -1), "repeated_visit": 0,
    }).to_csv(f"{stem}_perimeter_profile.csv", index=False)

    def line(mask):
        return LineString(list(zip(xs[mask], ys[mask])))

    gpd.GeoDataFrame({
        "run_id": [0, 1, 2], "verdict": ["source", "outlet_crossing", "too_narrow"],
        "width_m": [190.0, 190.0, 20.0], "min_elevation": [100.0, 80.0, 149.5],
        "relief_m": [9.0, 9.0, 0.5], "relief_before_m": [9.0, 9.0, 0.5],
        "relief_after_m": [9.0, 9.0, 0.5], "inward_checked_m": [500.0, 0.0, 0.0],
        "inward_max_gradient": [0.03, 0.0, 0.0], "failure_x": [np.nan] * 3,
        "failure_y": [np.nan] * 3, "outlet_id": [0, -1, -1],
        "path_length_to_outlet_m": [990.0, 0.0, 0.0],
    }, geometry=[line(north), line(south), line(west)], crs=CRS).to_file(
        f"{stem}_perimeter_crossings.fgb", driver="FlatGeobuf")
    gpd.GeoDataFrame({
        "outlet_id": [0, 1, 2], "rank": [1, 2, 3],
        "summed_drainage_area_m2": [9e5, 1e4, 5e3], "n_cells": [19, 1, 1],
        "width_m": [190.0, 10.0, 10.0], "on_crossing": [1, 0, 0],
        "selected": [1, 0, 0],
    }, geometry=[Point(505, 5), Point(995, 305), Point(5, 805)], crs=CRS).to_file(
        f"{stem}_perimeter_outlets.fgb", driver="FlatGeobuf")
    gpd.GeoDataFrame({
        "source_id": [0], "outlet_id": [0], "row": [0], "col": [50],
        "elevation": [100.0], "crossing_width_m": [190.0], "crossing_relief_m": [9.0],
        "path_length_to_outlet_m": [990.0], "drop_to_outlet_m": [20.0],
    }, geometry=[Point(505, 995)], crs=CRS).to_file(
        f"{stem}_perimeter_sources.fgb", driver="FlatGeobuf")
    (tmp_path / "p_ingested_params.param").write_text(
        "run_perimeter_sources: true\nperimeter_max_gradient: 0.07\n")

    yy, xx = np.mgrid[0:n, 0:n]
    dem = (100.0 + 0.1 * np.abs(xx * 10.0 + 5 - 500.0) + 0.02 * yy * 10.0).astype("float32")
    with rasterio.open(tmp_path / "dem.tif", "w", driver="GTiff", height=n, width=n,
                       count=1, dtype="float32", crs=CRS,
                       transform=from_origin(0.0, 1000.0, 10.0, 10.0)) as dst:
        dst.write(dem, 1)
    return stem


def test_map_and_profile_write(tmp_path):
    stem = _write_run(tmp_path)
    m = plot_perimeter_sources_map(stem, tmp_path / "map.png", dem=tmp_path / "dem.tif")
    p = plot_perimeter_profile(stem, tmp_path / "prof.png")
    assert m.stat().st_size > 0 and p.stat().st_size > 0


def test_max_gradient_read_from_the_ingested_params(tmp_path):
    stem = _write_run(tmp_path)
    assert _ingested_value(stem, "perimeter_max_gradient", 0.05) == pytest.approx(0.07)
    assert _ingested_value(stem, "no_such_key", 0.05) == pytest.approx(0.05)
    assert _ingested_value(tmp_path / "missing", "perimeter_max_gradient", 0.05) == 0.05


def test_missing_outputs_name_the_parameter(tmp_path):
    with pytest.raises(FileNotFoundError, match="run_perimeter_sources"):
        plot_perimeter_profile(tmp_path / "none", tmp_path / "x.png")


def test_cli_writes_both(tmp_path):
    stem = _write_run(tmp_path)
    assert plot_perimeter_sources.main([
        "--stem", str(stem), "--dem", str(tmp_path / "dem.tif"),
        "--out-map", str(tmp_path / "m.png"), "--out-profile", str(tmp_path / "p.png"),
        "--dpi", "60"]) == 0
    assert (tmp_path / "m.png").exists() and (tmp_path / "p.png").exists()


def test_cli_map_needs_a_backdrop(tmp_path):
    with pytest.raises(SystemExit, match="--base or --dem"):
        plot_perimeter_sources.main(["--stem", "x", "--out-map", str(tmp_path / "m.png")])
