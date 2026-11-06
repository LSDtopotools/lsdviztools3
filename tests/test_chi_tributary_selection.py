"""Tributary selection for chi plots (lsdtt3 _chi_tributaries table)."""

from __future__ import annotations

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import LineCollection
from shapely.geometry import LineString

from lsdviztools3.io.chi import TRIBUTARY_COLUMNS, removed_mask, select_tributaries


def _tributaries():
    """Basin 0: main stem 0; tributaries 1 (big, steep), 2 (small), 3 (drains
    into 2), 4 (medium). Basin 1: main stem 0, tributary 1."""
    rows = []

    def add(b, br, parent, area, frac, rank, ksn_max):
        r = {c: 0 for c in TRIBUTARY_COLUMNS}
        r.update(basin_index=b, basin_key=b, source_key=10 * b + br, branch_id=br,
                 parent_branch_id=parent, junction_drainage_area_m2=area,
                 basin_drainage_area_m2=1e7, area_fraction=frac, area_rank_in_basin=rank,
                 ksn_max=ksn_max, junction_x=100.0 * br, junction_y=50.0 * b)
        rows.append(r)

    add(0, 0, -1, 1e7, 1.0, 0, 900.0)       # main stem: steep, but never removed
    add(0, 1, 0, 3e6, 0.3, 1, 2000.0)
    add(0, 2, 0, 5e4, 0.005, 3, 50.0)
    add(0, 3, 2, 2e4, 0.002, 4, 40.0)
    add(0, 4, 0, 1e6, 0.1, 2, 60.0)
    add(1, 0, -1, 1e7, 1.0, 0, 80.0)
    add(1, 1, 0, 4e5, 0.04, 1, 70.0)
    return pd.DataFrame(rows)


def _ids(df):
    return set(zip(df["basin_index"], df["branch_id"]))


def test_ksn_outlier_removes_tributary_not_main_stem():
    kept, removed = select_tributaries(_tributaries(), ksn_max=1000.0)
    assert _ids(removed) == {(0, 1)}
    assert removed["reason"].tolist() == ["ksn_outlier"]
    assert (0, 0) in _ids(kept)


def test_area_rules_and_descendants():
    _, removed = select_tributaries(_tributaries(), min_junction_area_m2=1e5)
    assert _ids(removed) == {(0, 2), (0, 3)}
    _, removed = select_tributaries(_tributaries(), min_area_fraction=0.05)
    assert _ids(removed) == {(0, 2), (0, 3), (1, 1)}
    r = dict(zip(zip(removed["basin_index"], removed["branch_id"]), removed["reason"]))
    assert r[(0, 2)] == "small_fraction" and r[(0, 3)] == "small_fraction"


def test_descendants_follow_unless_kept():
    t = _tributaries()
    t.loc[(t.basin_index == 0) & (t.branch_id == 3), "junction_drainage_area_m2"] = 5e5
    _, removed = select_tributaries(t, min_junction_area_m2=1e5)
    r = dict(zip(zip(removed["basin_index"], removed["branch_id"]), removed["reason"]))
    assert r == {(0, 2): "small_area", (0, 3): "drains_to_removed"}
    _, removed = select_tributaries(t, min_junction_area_m2=1e5, remove_descendants=False)
    assert _ids(removed) == {(0, 2)}


def test_top_n_by_area_and_reason_order():
    _, removed = select_tributaries(_tributaries(), top_n=2, ksn_max=1000.0)
    r = dict(zip(zip(removed["basin_index"], removed["branch_id"]), removed["reason"]))
    # (0,1) is rank 1 but a k_sn outlier: the first rule wins
    assert r == {(0, 1): "ksn_outlier", (0, 2): "rank", (0, 3): "rank"}


def test_named_ksn_threshold_from_lsdtt3_table():
    th = pd.DataFrame({"measure": ["ksn"], "method": ["q3_plus_3iqr"], "threshold": [500.0],
                       "n_above": [1], "n_total": [10]})
    _, removed = select_tributaries(_tributaries(), ksn_max="q3_plus_3iqr", thresholds=th)
    assert _ids(removed) == {(0, 1)}
    with pytest.raises(ValueError, match="thresholds"):
        select_tributaries(_tributaries(), ksn_max="q3_plus_3iqr")


def test_removed_mask_matches_on_basin_and_branch():
    _, removed = select_tributaries(_tributaries(), ksn_max=1000.0)
    frame = pd.DataFrame({"basin_index": [0, 0, 1], "branch_id": [1, 2, 1]})
    assert removed_mask(frame, removed).tolist() == [True, False, False]
    assert removed_mask(frame, None).tolist() == [False, False, False]


def _segments():
    rows = []
    for br, k in ((0, 50.0), (1, 5000.0), (2, 80.0)):
        rows.append(dict(basin_index=0, outlet_flat=0, segment_id=br, parent_segment_id=-1,
                         branch_id=br, downstream_flat=0, upstream_flat=1,
                         downstream_chi=float(br), upstream_chi=br + 1.0,
                         downstream_elevation=10.0 * br, upstream_elevation=10.0 * br + k / 100,
                         k_sn=k, max_abs_residual=0.0, n_nodes=5, source_key=br))
    return pd.DataFrame(rows)


def test_profile_ghosts_removed_and_scales_on_kept():
    from lsdviztools3.charts import plot_chi_profile
    removed = pd.DataFrame({"basin_index": [0], "branch_id": [1], "reason": ["ksn_outlier"]})
    fig, ax = plot_chi_profile(_segments(), basin=0, removed=removed)
    lcs = [c for c in ax.collections if isinstance(c, LineCollection)]
    coloured = [c for c in lcs if c.get_array() is not None]
    ghost = [c for c in lcs if c.get_array() is None]
    assert len(coloured) == 1 and len(coloured[0].get_segments()) == 2
    assert len(ghost) == 1 and len(ghost[0].get_segments()) == 1
    # the colour range ignores the removed 5000
    assert coloured[0].norm.vmax == pytest.approx(np.log10(80.0))
    fig2, ax2 = plot_chi_profile(_segments(), basin=0, removed=removed, show_removed="hide")
    assert len([c for c in ax2.collections if isinstance(c, LineCollection)]) == 1
    with pytest.raises(ValueError, match="show_removed"):
        plot_chi_profile(_segments(), basin=0, removed=removed, show_removed="fade")
    import matplotlib.pyplot as plt
    plt.close("all")


def test_stacked_hides_removed():
    from lsdviztools3.charts import plot_chi_stacked
    removed = pd.DataFrame({"basin_index": [0], "branch_id": [1]})
    fig, ax = plot_chi_stacked(_segments(), removed=removed, show_removed="hide")
    lc = next(c for c in ax.collections if isinstance(c, LineCollection))
    assert len(lc.get_segments()) == 2
    import matplotlib.pyplot as plt
    plt.close(fig)


def _lines_and_dem(tmp_path):
    import rasterio
    from rasterio.transform import from_origin
    dem = tmp_path / "dem.tif"
    z = np.add.outer(np.arange(20.0), np.arange(20.0)).astype("float32")
    with rasterio.open(dem, "w", driver="GTiff", height=20, width=20, count=1, dtype="float32",
                       crs="EPSG:32633", transform=from_origin(500000, 4200200, 10, 10)) as d:
        d.write(z, 1)
    gdf = gpd.GeoDataFrame(
        {"basin_index": [0, 0, 0], "branch_id": [0, 1, 2], "k_sn": [50.0, 5000.0, 80.0]},
        geometry=[LineString([(500010 + 50 * i, 4200010), (500040 + 50 * i, 4200150)]) for i in range(3)],
        crs="EPSG:32633")
    lines = tmp_path / "seg.fgb"
    gdf.to_file(lines, driver="FlatGeobuf")
    return dem, lines


def test_removed_map_writes_png(tmp_path):
    from lsdviztools3.workflows import render_chi_ksn_map, render_chi_removed_map
    dem, lines = _lines_and_dem(tmp_path)
    removed = pd.DataFrame({"basin_index": [0], "branch_id": [1], "reason": ["ksn_outlier"],
                            "junction_x": [500060.0], "junction_y": [4200010.0]})
    out = render_chi_removed_map(dem, lines, removed, tmp_path / "removed.png")
    assert out.exists() and out.stat().st_size > 1000
    out2 = render_chi_ksn_map(dem, lines, tmp_path / "ksn.png", removed=removed)
    assert out2.exists()


def test_cli_selection_writes_removed_csv(tmp_path):
    from lsdviztools3.cli import plot_chi
    d = tmp_path / "run"
    d.mkdir()
    seg = _segments()
    seg.to_csv(d / "x_chi_fit_segments_mn_0p450.csv", index=False)
    nodes = pd.DataFrame({"basin_index": [0, 0, 0], "outlet_flat": 0, "flat": [0, 1, 2],
                          "receiver_flat": [-2, 0, 1], "chi": [0.0, 1.0, 2.0],
                          "observed_elevation": [0.0, 1.0, 2.0], "fitted_elevation": [0.0, 1.0, 2.0],
                          "residual_elevation": 0.0, "abs_residual_elevation": 0.0,
                          "branch_id": [0, 1, 2], "is_main_stem": [1, 0, 0], "is_source": 0,
                          "is_junction": 0, "is_split_point": 0, "source_key": [0, 1, 2]})
    nodes.to_csv(d / "x_chi_fit_nodes_mn_0p450.csv", index=False)
    _tributaries().to_csv(d / "x_chi_tributaries_mn_0p450.csv", index=False)
    out = tmp_path / "out"
    assert plot_chi.main(["--data-dir", str(d), "--prefix", "x", "--m-over-n", "0.45",
                          "--out-dir", str(out), "--remove-ksn-above", "1000"]) == 0
    r = pd.read_csv(out / "x_removed_tributaries.csv")
    assert r["reason"].tolist() == ["ksn_outlier"]
    with pytest.raises(SystemExit):
        plot_chi.main(["--data-dir", str(d), "--prefix", "x", "--m-over-n", "0.45",
                       "--out-dir", str(out), "--remove-ksn-above", "q3_plus_3iqr"])
