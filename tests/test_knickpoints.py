"""Knickpoint selection and plots (lsdtt3 _chi_knickpoints table)."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from lsdviztools3.charts.knickpoints import (
    plot_knickpoint_distribution,
    plot_knickpoints_along_chi,
    plot_knickpoints_on_profile,
    resolve_threshold,
    select_major_knickpoints,
)
from lsdviztools3.io.chi import KNICKPOINT_COLUMNS


def _knickpoints():
    """Basin 0: deltas +500, -400, +50, -30 (branches 0, 0, 1, 1) and a junction
    row of +900. Basin 1: +600 on branch 2."""
    rows = []

    def add(b, br, d, junction=0, chi=1.0):
        r = {c: 0 for c in KNICKPOINT_COLUMNS}
        r.update(basin_index=b, branch_id=br, delta_ksn=d, log_ratio_ksn=np.log10(2.0) * np.sign(d),
                 at_junction=junction, chi=chi, chi_normalised=chi / 10, elevation=100 + 10 * chi,
                 x=500000.0 + 100 * chi, y=4200000.0 + 10 * b, source_key=br)
        rows.append(r)

    add(0, 0, 500.0, chi=1.0)
    add(0, 0, -400.0, chi=2.0)
    add(0, 1, 50.0, chi=3.0)
    add(0, 1, -30.0, chi=4.0)
    add(0, 1, 900.0, junction=1, chi=5.0)
    add(1, 2, 600.0, chi=6.0)
    return pd.DataFrame(rows)


def _thresholds():
    return pd.DataFrame({
        "measure": ["abs_delta_ksn"] * 4 + ["abs_log_ratio_ksn"] * 4 + ["ksn"],
        "method": ["p90", "p95", "median_plus_3mad", "otsu"] * 2 + ["q3_plus_3iqr"],
        "threshold": [450.0, 550.0, 300.0, 100.0, 0.2, 0.25, 0.3, 0.1, 800.0],
        "n_above": 0, "n_total": 5})


def test_threshold_by_value_and_by_method():
    assert resolve_threshold(250, "delta")[0] == 250.0
    assert resolve_threshold("p95", "delta", _thresholds())[0] == 550.0
    assert resolve_threshold("otsu", "log_ratio", _thresholds())[0] == 0.1
    with pytest.raises(ValueError, match="thresholds table"):
        resolve_threshold("p95", "delta")
    with pytest.raises(ValueError, match="no abs_delta_ksn"):
        resolve_threshold("jenks", "delta", _thresholds())


def test_select_both_directions_without_junctions():
    rows, v, _ = select_major_knickpoints(_knickpoints(), 350)
    assert v == 350 and sorted(rows["delta_ksn"]) == [-400.0, 500.0, 600.0]
    assert (rows["at_junction"] == 0).all()
    assert "measure_value" in rows


def test_select_direction_and_junctions():
    up, _, _ = select_major_knickpoints(_knickpoints(), 350, direction="steepening")
    assert sorted(up["delta_ksn"]) == [500.0, 600.0]
    down, _, _ = select_major_knickpoints(_knickpoints(), 350, direction="shallowing")
    assert down["delta_ksn"].tolist() == [-400.0]
    withj, _, _ = select_major_knickpoints(_knickpoints(), 350, include_junctions=True)
    assert 900.0 in withj["delta_ksn"].tolist()
    with pytest.raises(ValueError, match="direction"):
        select_major_knickpoints(_knickpoints(), 350, direction="up")


def test_select_drops_removed_tributaries():
    removed = pd.DataFrame({"basin_index": [1], "branch_id": [2]})
    rows, _, _ = select_major_knickpoints(_knickpoints(), 350, removed=removed)
    assert 600.0 not in rows["delta_ksn"].tolist()


def test_select_named_method():
    rows, v, label = select_major_knickpoints(_knickpoints(), "p90", thresholds=_thresholds())
    assert v == 450.0 and "p90" in label
    assert sorted(rows["delta_ksn"]) == [500.0, 600.0]


def test_distribution_draws_thresholds_and_chosen():
    fig, (ah, ar) = plot_knickpoint_distribution(_knickpoints(), _thresholds(), chosen=350.0)
    # 4 methods x 2 sides + chosen x 2 on the histogram
    assert len(ah.get_lines()) == 10
    labels = [t.get_text() for t in ar.get_legend().get_texts()]
    assert any(l.startswith("chosen: 350") for l in labels)
    assert any("otsu: 100 (3 above)" in l for l in labels)
    plt.close(fig)


def _segments():
    rows = []
    for b in (0, 1):
        for i in range(3):
            rows.append(dict(basin_index=b, outlet_flat=0, segment_id=i, parent_segment_id=i - 1,
                             branch_id=0, downstream_flat=i, upstream_flat=i + 1,
                             downstream_chi=2.0 * i, upstream_chi=2.0 * (i + 1),
                             downstream_elevation=100 + 20.0 * i, upstream_elevation=100 + 20.0 * (i + 1),
                             k_sn=10.0 * (i + 1), max_abs_residual=0.0, n_nodes=5, source_key=0))
    return pd.DataFrame(rows)


def test_profile_markers_on_stack_positions():
    from lsdviztools3.charts.chi import stack_offsets
    major, _, _ = select_major_knickpoints(_knickpoints(), 350)
    fig, ax = plot_knickpoints_on_profile(_segments(), major)
    shifts = stack_offsets(_segments())
    offs = np.concatenate([c.get_offsets() for c in ax.collections if hasattr(c, "get_offsets")
                           and c.__class__.__name__ == "PathCollection"])
    xs = sorted(offs[:, 0].tolist())
    expected = sorted([1.0 + shifts[0], 2.0 + shifts[0], 6.0 + shifts[1]])
    assert xs == pytest.approx(expected)
    plt.close(fig)
    fig, ax = plot_knickpoints_on_profile(_segments(), major, basin=0)
    n = sum(len(c.get_offsets()) for c in ax.collections if c.__class__.__name__ == "PathCollection")
    assert n == 2
    plt.close(fig)


def test_along_chi_counts():
    major, _, _ = select_major_knickpoints(_knickpoints(), 350)
    fig, (ah, asr) = plot_knickpoints_along_chi(_knickpoints(), major, along="chi_normalised")
    labels = [t.get_text() for t in ah.get_legend().get_texts()]
    assert labels == ["all segment boundaries (5)", "major knickpoints (3)"]
    plt.close(fig)
    with pytest.raises(ValueError, match="along"):
        plot_knickpoints_along_chi(_knickpoints(), major, along="distance")


def test_map_writes_png(tmp_path):
    import geopandas as gpd
    import rasterio
    from rasterio.transform import from_origin
    from shapely.geometry import LineString
    from lsdviztools3.charts.knickpoints import render_knickpoint_map
    dem = tmp_path / "dem.tif"
    with rasterio.open(dem, "w", driver="GTiff", height=30, width=80, count=1, dtype="float32",
                       crs="EPSG:32633", transform=from_origin(499990, 4200100, 10, 10)) as d:
        d.write(np.add.outer(np.arange(30.0), np.arange(80.0)).astype("float32"), 1)
    lines = tmp_path / "seg.fgb"
    gpd.GeoDataFrame({"basin_index": [0], "branch_id": [0], "k_sn": [50.0]},
                     geometry=[LineString([(500000, 4200000), (500700, 4200010)])],
                     crs="EPSG:32633").to_file(lines, driver="FlatGeobuf")
    major, _, _ = select_major_knickpoints(_knickpoints(), 350)
    out = render_knickpoint_map(dem, lines, major, tmp_path / "kp.png", threshold_label="350")
    assert out.exists() and out.stat().st_size > 1000


def test_cli(tmp_path):
    from lsdviztools3.cli import plot_knickpoints
    d = tmp_path / "run"
    d.mkdir()
    _knickpoints().to_csv(d / "x_chi_knickpoints_mn_0p450.csv", index=False)
    _thresholds().to_csv(d / "x_chi_knickpoint_thresholds_mn_0p450.csv", index=False)
    _segments().to_csv(d / "x_chi_fit_segments_mn_0p450.csv", index=False)
    out = tmp_path / "out"
    assert plot_knickpoints.main(["--data-dir", str(d), "--prefix", "x", "--m-over-n", "0.45",
                                  "--out-dir", str(out), "--threshold", "p90"]) == 0
    m = pd.read_csv(out / "x_major_knickpoints.csv")
    assert len(m) == 2 and (m["threshold"] == 450.0).all()
    for f in ("knickpoint_distribution", "knickpoint_profile", "knickpoints_along_chi"):
        assert (out / f"x_{f}.png").exists()
    with pytest.raises(SystemExit):
        plot_knickpoints.main(["--data-dir", str(d), "--prefix", "x", "--m-over-n", "0.45",
                               "--out-dir", str(out), "--threshold", "jenks"])
