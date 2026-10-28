"""Method comparison maps and charts (lsdtt3 fluvial-surfaces, plan 7g)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin

from lsdviztools3.charts.surface_methods import plot_mapped_scores, plot_method_levels
from lsdviztools3.render.surface_methods import (
    agreement_combinations,
    method_stem,
    read_method_list,
)

METHODS = ["ladder", "clubb2017"]


def _write_run(tmp_path):
    """A two-method run's CSVs and rasters under stem tmp_path/r."""
    stem = tmp_path / "r"
    pd.DataFrame({"method": METHODS, "floodplain_cells": [10, 12]}).to_csv(
        f"{stem}_method_summary.csv", index=False)
    for m, med in zip(METHODS, ([2.0, 12.0, 12.5, 30.0], [3.0, 20.0])):
        n = len(med)
        pd.DataFrame({
            "label": range(1, n + 1), "band": range(n),
            "class": ["floodplain"] + ["terrace"] * (n - 1),
            "source_id": [0] * n, "n_pixels": [10] * n,
            "area_m2": [1000.0 * (k + 1) for k in range(n)],
            "median_field_m": med,
        }).to_csv(f"{stem}_{m}_surface_patches.csv", index=False)
    pd.DataFrame({
        "scope": ["all", "window"], "source_id": [-1, 0],
        "window_centre_m": [0, 500], "n_patches": [4, 4],
        "required_gap_m": [1.5, 1.5], "n_clusters": [3, 3],
        "empty_fraction": [0.9, 0.9], "null_mean_empty_fraction": [0.5, 0.5],
        "p_value": [0.02, 0.02], "cluster_centres_m": ["2 12.25 30", "2 12.25 30"],
        "cluster_n_patches": ["1 2 1", "1 2 1"],
    }).to_csv(f"{stem}_ladder_patch_level_clusters.csv", index=False)
    pd.DataFrame({
        "method": METHODS, "n_scored_cells": [100, 100], "cell_area_m2": [25, 25],
        "floodplain_precision": [0.9, 1.0], "floodplain_recall": [0.5, 0.2],
        "floodplain_iou": [0.47, 0.2], "terrace_precision": [0.8, np.nan],
        "terrace_recall": [0.4, 0.0], "terrace_iou": [0.36, 0.0],
    }).to_csv(f"{stem}_mapped_surface_scores.csv", index=False)

    crs = "EPSG:32632"
    tr = from_origin(440000, 5350600, 10, 10)
    def write(path, a):
        with rasterio.open(path, "w", driver="GTiff", height=a.shape[0],
                           width=a.shape[1], count=1, dtype="float32", crs=crs,
                           transform=tr, nodata=-9999) as ds:
            ds.write(a.astype("float32"), 1)
    write(tmp_path / "hs.tif", (np.indices((40, 40)).sum(axis=0) % 255))
    cls = np.zeros((40, 40)); cls[:, 15:20] = 1; cls[:, 20:30] = 2
    for m in METHODS:
        write(f"{stem}_{m}_surface_class.tif", cls)
    agree = np.zeros((40, 40)); agree[:, 15:18] = 3; agree[:, 18:20] = 1
    write(f"{stem}_method_agreement_floodplain.tif", agree)
    return stem


def test_method_list_and_stems(tmp_path):
    stem = _write_run(tmp_path)
    assert read_method_list(stem) == METHODS
    assert method_stem(stem, "ladder", METHODS) == f"{stem}_ladder"
    assert method_stem(stem, "ladder", ["ladder"]) == str(stem)


def test_method_list_single_run_reads_the_record(tmp_path):
    stem = tmp_path / "s"
    pd.DataFrame({"method": ["clubb2017"], "item": ["x"], "value": ["1"],
                  "source": ["preset"]}).to_csv(f"{stem}_method_record.csv", index=False)
    assert read_method_list(stem) == ["clubb2017"]


def test_agreement_combinations_follow_the_bits():
    c = agreement_combinations(["a", "b", "c"])
    assert len(c) == 7
    assert c[0] == (1, "a") and c[2] == (3, "a + b") and c[4] == (5, "a + c")
    assert c[-1] == (7, "all")


def test_levels_chart_one_column_per_method(tmp_path):
    stem = _write_run(tmp_path)
    out = tmp_path / "levels.png"
    fig = plot_method_levels(stem, str(out))
    ax = fig.axes[0]
    labels = [t.get_text() for t in ax.get_xticklabels()]
    assert labels[0] == "ladder\np = 0.02" and labels[1] == "clubb2017"
    # Only the cluster of >= 3 patches would be drawn; none here.
    assert len(ax.lines) == 0
    fig = plot_method_levels(stem, None, min_cluster_patches=2)
    assert len(fig.axes[0].lines) == 1       # the 12.25 m cluster of 2
    assert out.exists()


def test_scores_chart(tmp_path):
    stem = _write_run(tmp_path)
    out = tmp_path / "scores.png"
    fig = plot_mapped_scores(stem, str(out))
    assert [a.get_title() for a in fig.axes] == ["floodplain", "terrace"]
    assert out.exists()


def test_maps_write_pngs(tmp_path):
    pytest.importorskip("pygmt")
    from lsdviztools3.render.surface_methods import (
        render_method_agreement,
        render_method_classes,
    )
    stem = _write_run(tmp_path)
    fig = render_method_classes(stem, tmp_path / "hs.tif", panel_width="6c")
    fig.savefig(tmp_path / "classes.png")
    fig = render_method_agreement(stem, tmp_path / "hs.tif", "floodplain")
    fig.savefig(tmp_path / "agree.png")
    assert (tmp_path / "classes.png").stat().st_size > 0
    assert (tmp_path / "agree.png").stat().st_size > 0
