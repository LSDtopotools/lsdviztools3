"""lsdtt-slope-stability outputs: readers, charts, maps and the lsdviz3-slope-stability CLI."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from lsdviztools3.charts import slope_stability as ch
from lsdviztools3.io import slope_stability as ss_io

NR, NC = 30, 40
CATALOGUE = ["area_m2", "volume_m3", "mean_depth_m", "max_depth_m", "length_m", "width_m",
             "length_to_width"]


def _write(path, values):
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(path, "w", driver="GTiff", height=NR, width=NC, count=1, dtype="float32",
                       crs="EPSG:32633", transform=from_origin(499990, 4200100, 10, 10),
                       nodata=-9999.0) as d:
        d.write(np.where(np.isfinite(values), values, -9999.0).astype("float32"), 1)
    return path


def _catalogue(n, seed):
    rng = np.random.default_rng(seed)
    area = 100.0 * rng.integers(1, 30, n)
    depth = rng.uniform(1, 5, n)
    length = rng.uniform(10, 80, n)
    width = rng.uniform(10, 40, n)
    return pd.DataFrame({"area_m2": area, "volume_m3": area * depth, "mean_depth_m": depth,
                         "max_depth_m": depth * 1.5, "length_m": length, "width_m": width,
                         "length_to_width": length / width})


def _curves():
    a = np.logspace(1, 4, 20)
    rows = []
    for name, s in (("always_saturated_above", 0.001 * a), ("never_saturated_below", 0.002 * a),
                    ("fs_min_1.5", np.full_like(a, 0.3)), ("fs_min_1.0", np.full_like(a, 0.5)),
                    ("fs_max_1.0", np.full_like(a, 0.9))):
        rows += [{"region_id": 0, "curve": name, "slope": v, "specific_catchment_area_m": x}
                 for x, v in zip(a, s)]
    return pd.DataFrame(rows)


def _run_dir(tmp_path, *, sinmap=True, catalogues=True, slip=True):
    d = tmp_path / "run"
    d.mkdir()
    r, c = np.mgrid[0:NR, 0:NC]
    if sinmap:
        _write(d / "ss_sinmap_stability_class.tif", 1.0 + (c % 6))
        _write(d / "ss_dinf_slope.tif", 0.05 + 0.02 * r)
        _write(d / "ss_specific_catchment_area.tif", 10.0 * (1 + c))
        _curves().to_csv(d / "ss_sinmap_slope_area_curves.csv", index=False)
        pd.DataFrame({"region_id": 0, "stability_class": np.arange(1, 7), "n_cells": 200,
                      "area_m2": 20000.0, "fraction_of_region": 1 / 6}).to_csv(
            d / "ss_sinmap_class_areas.csv", index=False)
    depth = np.where(c < NC - 4, 0.5 + 0.1 * r, np.nan)
    _write(d / "ss_rupture_depth.tif", depth)
    _write(d / "ss_rupture_failure_probability.tif", np.where(np.isfinite(depth), (c % 5) / 4, np.nan))
    if catalogues:
        _catalogue(40, 1).to_csv(d / "ss_landslide_catalogue.csv", index=False)
        _catalogue(30, 2).to_csv(d / "ss_landslide_catalogue_amalgamated.csv", index=False)
    if slip:
        evaluated = np.where(np.isfinite(depth), (r % 3 != 0).astype(float), np.nan)
        F = np.where(evaluated == 1, 0.6 + 0.03 * c, np.nan)
        _write(d / "ss_slip_evaluated.tif", evaluated)
        _write(d / "ss_slip_safety_factor.tif", F)
        _write(d / "ss_slip_max_depth.tif", np.where(evaluated == 1, depth * 2, np.nan))
        _write(d / "ss_slip_body_depth.tif", np.where(F <= 1, depth * 2, np.nan))
        _catalogue(20, 3).to_csv(d / "ss_slip_body_catalogue.csv", index=False)
        _catalogue(15, 4).to_csv(d / "ss_slip_body_catalogue_amalgamated.csv", index=False)
    return d


def _dem(tmp_path):
    r, c = np.mgrid[0:NR, 0:NC]
    return _write(tmp_path / "dem.tif", 100.0 + 2.0 * r + np.sin(c / 3.0))


# -----------------------------------------------------------------------------
# Readers
# -----------------------------------------------------------------------------
def test_files_and_present(tmp_path):
    d = _run_dir(tmp_path, slip=False)
    f = ss_io.present(ss_io.slope_stability_files(d, "ss"))
    assert {"sinmap_class", "rupture_depth", "landslides"} <= set(f)
    assert "slip_safety_factor" not in f
    assert all(p.name.startswith("ss_") for p in f.values())


def test_read_band_nodata_and_catalogue_check(tmp_path):
    d = _run_dir(tmp_path)
    v = ss_io.read_band(d / "ss_rupture_depth.tif")
    assert np.isnan(v[:, -1]).all() and np.isfinite(v[:, 0]).all()
    bad = tmp_path / "bad.csv"
    pd.DataFrame({"area_m2": [1.0]}).to_csv(bad, index=False)
    with pytest.raises(KeyError):
        ss_io.load_catalogue(bad)
    assert set(CATALOGUE) <= set(ss_io.load_catalogue(d / "ss_landslide_catalogue.csv").columns)


# -----------------------------------------------------------------------------
# Charts
# -----------------------------------------------------------------------------
def test_sinmap_slope_area_draws_every_class_and_curve(tmp_path):
    d = _run_dir(tmp_path)
    fig, ax = ch.plot_sinmap_slope_area(ss_io.read_band(d / "ss_dinf_slope.tif"),
                                        ss_io.read_band(d / "ss_specific_catchment_area.tif"),
                                        ss_io.read_band(d / "ss_sinmap_stability_class.tif"),
                                        ss_io.load_sinmap_curves(d / "ss_sinmap_slope_area_curves.csv"))
    assert len(ax.collections) == 6          # one scatter per class
    assert len(ax.get_lines()) == 5          # the five class boundaries
    assert ax.get_xscale() == "log"
    plt.close(fig)


def test_sinmap_subsampling_is_deterministic():
    assert len(ch._subsample(10, 20)) == 10
    idx = ch._subsample(1000, 100)
    assert len(idx) <= 100 and idx[0] == 0 and np.all(np.diff(idx) == idx[1])


def test_catalogue_charts_one_series_per_catalogue():
    cats = [("a", _catalogue(40, 1)), ("b", _catalogue(30, 2))]
    for fn in (ch.plot_catalogue_area_volume, ch.plot_catalogue_area_depth):
        fig, ax = fn(cats)
        assert len(ax.collections) == 2 and ax.get_xscale() == "log" and ax.get_yscale() == "log"
        plt.close(fig)
    fig, ax = ch.plot_catalogue_area_distribution(cats)
    assert len(ax.patches) == 2              # one step line per catalogue
    plt.close(fig)
    fig, ax = ch.plot_catalogue_length_width(cats)
    assert len(ax.patches) == 2
    plt.close(fig)


def test_area_distribution_counts_every_landslide():
    df = _catalogue(50, 7)
    fig, ax = ch.plot_catalogue_area_distribution([("a", df)])
    heights = np.nan_to_num(ax.patches[0].get_data().values)
    assert heights.sum() == len(df)
    plt.close(fig)


def test_depth_comparison(tmp_path):
    d = _run_dir(tmp_path)
    fig, ax = ch.plot_depth_comparison(ss_io.read_band(d / "ss_rupture_depth.tif"),
                                       ss_io.read_band(d / "ss_slip_max_depth.tif"))
    assert len(ax.collections) == 1 and len(ax.get_lines()) == 1   # crowns and the 1:1 line
    plt.close(fig)
    with pytest.raises(ValueError):
        ch.plot_depth_comparison(np.full(4, np.nan), np.ones(4))


def test_screened_mask(tmp_path):
    d = _run_dir(tmp_path)
    m = ch.screened_mask(d / "ss_slip_evaluated.tif")
    ev = ss_io.read_band(d / "ss_slip_evaluated.tif")
    v = np.asarray(m.values, dtype=float)
    assert np.array_equal(np.isfinite(v), ev == 0)


def test_class_colours_are_ordered_and_distinct():
    assert list(ch.SINMAP_CLASS_COLOURS) == [1, 2, 3, 4, 5, 6]
    assert len(set(ch.SINMAP_CLASS_COLOURS.values())) == 6
    # The screened-out grey must not be pure grey (GMT paints NaN red then).
    r, g, b = (int(x) for x in ch.SCREENED_COLOUR.split("/"))
    assert not (r == g == b)


# -----------------------------------------------------------------------------
# Maps and CLI
# -----------------------------------------------------------------------------
def test_value_map_with_screened_overlay(tmp_path):
    d = _run_dir(tmp_path)
    out = ch.render_value_map(_dem(tmp_path), d / "ss_slip_safety_factor.tif", tmp_path / "f.png",
                              cmap=ch.SAFETY_FACTOR_CMAP, label="F", vmin=0.5, vmax=1.5, title="t",
                              overlay_mask=ch.screened_mask(d / "ss_slip_evaluated.tif"),
                              overlay_label="screened")
    assert out.exists() and out.stat().st_size > 1000


def test_cli_end_to_end(tmp_path, capsys):
    from lsdviztools3.cli import plot_slope_stability
    d = _run_dir(tmp_path)
    out = tmp_path / "out"
    assert plot_slope_stability.main(["--data-dir", str(d), "--prefix", "ss", "--out-dir", str(out),
                                      "--out-prefix", "x", "--dem", str(_dem(tmp_path)), "--dpi", "50"]) == 0
    names = sorted(p.name for p in out.iterdir())
    assert names == sorted([
        "x_sinmap_slope_area.png",
        "x_landslide_area_volume.png", "x_landslide_area_depth.png",
        "x_landslide_area_distribution.png", "x_landslide_length_width.png",
        "x_slip_body_area_volume.png", "x_slip_body_area_depth.png",
        "x_slip_body_area_distribution.png", "x_slip_body_length_width.png",
        "x_depth_comparison.png",
        "x_sinmap_class_map.png", "x_rupture_depth_map.png", "x_rupture_failure_probability_map.png",
        "x_slip_safety_factor_map.png", "x_slip_body_depth_map.png"])
    printed = capsys.readouterr().out
    for n in names:
        assert str(out / n) in printed


def test_cli_without_dem_or_slip(tmp_path, capsys):
    from lsdviztools3.cli import plot_slope_stability
    d = _run_dir(tmp_path, sinmap=False, slip=False)
    out = tmp_path / "out"
    assert plot_slope_stability.main(["--data-dir", str(d), "--prefix", "ss", "--out-dir", str(out),
                                      "--dpi", "50"]) == 0
    assert sorted(p.name for p in out.iterdir()) == sorted([
        "ss_landslide_area_volume.png", "ss_landslide_area_depth.png",
        "ss_landslide_area_distribution.png", "ss_landslide_length_width.png"])
    printed = capsys.readouterr().out
    assert "no SINMAP slope-area plot" in printed
    assert "no slip bodies catalogue charts" in printed
    assert "give --dem" in printed


def test_cli_no_outputs(tmp_path):
    from lsdviztools3.cli import plot_slope_stability
    (tmp_path / "empty").mkdir()
    with pytest.raises(SystemExit):
        plot_slope_stability.main(["--data-dir", str(tmp_path / "empty"), "--prefix", "ss",
                                   "--out-dir", str(tmp_path / "o")])
