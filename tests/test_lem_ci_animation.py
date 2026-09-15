"""Tests for the channelization-index LEM animation (render/lem_ci_animation)."""

import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin

from lsdviztools3.render import lem_ci_animation as lci


def _write_dem(path, z):
    ny, nx = z.shape
    with rasterio.open(
        path, "w", driver="GTiff", height=ny, width=nx, count=1, dtype="float32",
        crs="EPSG:32601", transform=from_origin(500000.0, ny * 3.0, 3.0, 3.0),
        nodata=-9999.0,
    ) as dst:
        dst.write(z.astype("float32"), 1)


@pytest.fixture
def ci_run(tmp_path):
    """A tiny run: initial surface, three frames, and a stepped CI series."""
    y, x = np.mgrid[0:30, 0:45]
    base = 50.0 * np.sin(np.pi * y / 29.0)
    _write_dem(tmp_path / "run_initial_topography.tif", base)
    for i, it in enumerate((20, 40, 60), start=1):
        _write_dem(tmp_path / f"run_model_iteration_{it}.tif",
                   base * (1 + 0.2 * i) + np.cos(x / 3.0) * i)
    t = np.arange(1, 61) * 500.0
    ci = np.where(t <= 10000, 177.0, np.where(t <= 20000, 250.0, 450.0))
    pd.DataFrame({"iteration": np.arange(1, 61), "time_yr": t,
                  "channelization_index": ci}).to_csv(
        tmp_path / "run_model_elevation_stats.csv", index=False)
    return tmp_path


def test_frames_are_time_ordered_and_start_at_zero(ci_run):
    frames = lci.find_ci_frames(ci_run, prefix="run")
    assert [f.iteration for f in frames] == [0, 20, 40, 60]
    assert frames[0].time_yr == 0.0
    assert [f.time_yr for f in frames[1:]] == [10000.0, 20000.0, 30000.0]


def test_elevation_limits_span_all_frames(ci_run):
    frames = lci.find_ci_frames(ci_run, prefix="run")
    lo, hi = lci.elevation_limits(frames)
    # The last frame is the tallest (base * 1.6 + cos term); the pinned range
    # must include it, not just the first frame's 50 m.
    assert hi > 79.0 and lo <= 0.0


def test_missing_ci_column_is_an_error(ci_run):
    stats = ci_run / "run_model_elevation_stats.csv"
    pd.read_csv(stats).drop(columns="channelization_index").to_csv(stats, index=False)
    with pytest.raises(ValueError, match="channelization_index"):
        lci.animate_ci_run(ci_run, ci_run / "x.gif", prefix="run")


@pytest.mark.parametrize("view", ["map", "oblique"])
def test_animation_writes_one_gif_frame_per_raster(ci_run, view):
    from PIL import Image

    out = lci.animate_ci_run(ci_run, ci_run / f"ci_{view}.gif", view=view,
                             prefix="run", fps=10, dpi=40)
    with Image.open(out) as gif:
        assert gif.n_frames == 4
        assert gif.info["duration"] == 100
