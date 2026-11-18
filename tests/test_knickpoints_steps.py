"""Split-hierarchy measures, knickzones, steps and the step-aware (stage 2)
files in the knickpoint plots (lsdtt3 run_step_change_detection)."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from lsdviztools3.charts.knickpoints import (
    knickzone_reach_lines,
    plot_knickpoint_distribution,
    plot_knickpoints_on_profile,
    plot_steps_on_long_profile,
    resolve_threshold,
    select_major_knickpoints,
)
from lsdviztools3.io.chi import (
    KNICKPOINT_COLUMNS,
    KNICKPOINT_SPLIT_COLUMNS,
    KNICKZONE_COLUMNS,
    STEP_COLUMNS,
    attach_flow_distance,
    chi_fit_files,
    filter_chi_steps,
    load_chi_steps,
)


def _knickpoints():
    """Basin 0 branch 0: normalised splits +5 (major), -4 (major), +0.5;
    a junction row +9 (major in lsdtt3 would not be set at a junction);
    basin 1 branch 2: -6 (major)."""
    rows = []

    def add(b, br, sn, major, junction=0, chi=1.0):
        r = {c: 0 for c in KNICKPOINT_COLUMNS + KNICKPOINT_SPLIT_COLUMNS}
        r.update(basin_index=b, branch_id=br, split_delta_ksn_normalised=sn,
                 split_delta_ksn_metres=10 * sn, delta_ksn=100 * sn, is_major=major,
                 at_junction=junction, chi=chi, elevation=100 + 10 * chi, x=500000.0 + 100 * chi,
                 y=4200000.0, source_key=br)
        rows.append(r)

    add(0, 0, 5.0, 1, chi=1.0)
    add(0, 0, -4.0, 1, chi=2.0)
    add(0, 0, 0.5, 0, chi=3.0)
    add(0, 0, 9.0, 0, junction=1, chi=4.0)
    add(1, 2, -6.0, 1, chi=5.0)
    return pd.DataFrame(rows)


def _thresholds():
    return pd.DataFrame({"measure": ["abs_split_delta_ksn_normalised", "abs_split_delta_ksn_metres"],
                         "method": ["p95", "p95"], "threshold": [3.0, 45.0], "n_above": 0, "n_total": 4})


def _steps():
    rows = []
    for i, (b, br, drop, mouth) in enumerate([(0, 0, 30.0, 0), (0, 0, 8.0, 1), (1, 2, 120.0, 0)]):
        r = {c: 0 for c in STEP_COLUMNS}
        r.update(step_id=i, basin_index=b, branch_id=br, drop_m=drop, at_tributary_mouth=mouth,
                 bottom_chi=1.0 + i, top_chi=1.2 + i, bottom_elevation=100.0 + 10 * i,
                 top_elevation=100.0 + 10 * i + drop, bottom_flow_distance_m=1000.0 * (i + 1),
                 flow_length_m=50.0, bottom_x=500000.0 + 50 * i, bottom_y=4200000.0)
        rows.append(r)
    return pd.DataFrame(rows)


def _knickzones():
    rows = []
    for i, (b, br, ft, foot, lip) in enumerate([(0, 0, "split", 1.0, 2.0), (1, 2, "basin_outlet", 0.0, 5.0)]):
        r = {c: 0 for c in KNICKZONE_COLUMNS}
        r.update(knickzone_id=i, basin_index=b, branch_id=br, foot_type=ft, foot_chi=foot, lip_chi=lip,
                 foot_elevation=100 + 10 * foot, lip_elevation=100 + 10 * lip)
        rows.append(r)
    return pd.DataFrame(rows)


def _segments(is_step=False):
    rows = []
    for b, br in ((0, 0), (1, 2)):
        for i in range(3):
            rows.append(dict(basin_index=b, outlet_flat=0, segment_id=i, parent_segment_id=i - 1,
                             branch_id=br, downstream_flat=i, upstream_flat=i + 1,
                             downstream_chi=2.0 * i, upstream_chi=2.0 * (i + 1),
                             downstream_elevation=100 + 20.0 * i, upstream_elevation=100 + 20.0 * (i + 1),
                             k_sn=10.0 * (i + 1), max_abs_residual=0.0, n_nodes=5, source_key=br))
    df = pd.DataFrame(rows)
    if is_step:
        df["is_step"] = (df["segment_id"] == 1).astype(int)
    return df


def test_major_uses_lsdtt3_is_major_without_junctions():
    rows, v, label = select_major_knickpoints(_knickpoints(), "major", measure="split_normalised")
    assert sorted(rows["split_delta_ksn_normalised"]) == [-6.0, -4.0, 5.0]
    assert v == 4.0 and "is_major" in label
    up, _, _ = select_major_knickpoints(_knickpoints(), "major", measure="split_normalised",
                                        direction="steepening")
    assert up["split_delta_ksn_normalised"].tolist() == [5.0]


def test_major_needs_the_ranking_measure_and_columns():
    with pytest.raises(ValueError, match="split_normalised"):
        select_major_knickpoints(_knickpoints(), "major", measure="delta")
    old = _knickpoints().drop(columns=list(KNICKPOINT_SPLIT_COLUMNS))
    with pytest.raises(ValueError, match="no 'split_delta_ksn_normalised'"):
        select_major_knickpoints(old, "major", measure="split_normalised")
    with pytest.raises(ValueError, match="is_major"):
        resolve_threshold("major", "split_normalised", _thresholds())


def test_split_measure_by_method():
    rows, v, _ = select_major_knickpoints(_knickpoints(), "p95", measure="split", thresholds=_thresholds())
    assert v == 45.0 and sorted(rows["split_delta_ksn_metres"]) == [-60.0, 50.0]


def test_split_distribution_on_log_axes():
    fig, (ah, ar) = plot_knickpoint_distribution(_knickpoints(), _thresholds(), measure="split_normalised")
    assert ah.get_xscale() == "symlog" and ar.get_yscale() == "log"
    plt.close(fig)


def test_filter_steps_by_drop_flags_and_removed():
    s = filter_chi_steps(_steps())
    assert s["flagged"].tolist() == [False, True, False]
    assert filter_chi_steps(_steps(), min_drop=10)["step_id"].tolist() == [0, 2]
    assert filter_chi_steps(_steps(), hide_flagged=True)["step_id"].tolist() == [0, 2]
    removed = pd.DataFrame({"basin_index": [1], "branch_id": [2]})
    assert filter_chi_steps(_steps(), removed=removed)["step_id"].tolist() == [0, 1]


def test_load_steps_checks_columns(tmp_path):
    p = tmp_path / "s.csv"
    _steps().drop(columns=["drop_m"]).to_csv(p, index=False)
    with pytest.raises(ValueError, match="drop_m"):
        load_chi_steps(p)


def test_chi_fit_files_prefers_stage2(tmp_path):
    base = tmp_path / "x_chi"
    (tmp_path / "x_chi_fit_segments_mn_0p450.csv").write_text("a\n")
    f = chi_fit_files(base, "0p450")
    assert f["stage"] == "stage1" and f["knickpoints"].name == "x_chi_knickpoints_mn_0p450.csv"
    (tmp_path / "x_chi_fit_segments_stepaware_mn_0p450.csv").write_text("a\n")
    f = chi_fit_files(base, "0p450")
    assert f["stage"] == "stepaware"
    assert f["fit_segments"].name == "x_chi_fit_segments_stepaware_mn_0p450.csv"
    # knickpoints are ranked on stage 1 only
    assert f["knickpoints"].name == "x_chi_knickpoints_mn_0p450.csv"
    assert f["knickzones"].name == "x_chi_knickzones_mn_0p450.csv"
    assert f["steps"].name == "x_chi_steps_mn_0p450.csv"
    assert chi_fit_files(base, "0p450", stage="stage1")["stage"] == "stage1"
    with pytest.raises(ValueError, match="stage"):
        chi_fit_files(base, "0p450", stage="two")


def test_attach_flow_distance_by_source_and_nearest_chi():
    nodes = pd.DataFrame({"source_key": [3, 3, 4, -1], "chi": [0.10001, 0.5, 0.2, 0.3]})
    points = pd.DataFrame({"source_key": [3, 3, 3, 4], "chi": [0.0, 0.1, 0.5, 0.2],
                           "flow_distance": [0.0, 10.0, 50.0, 20.0]})
    out = attach_flow_distance(nodes, points)
    assert out["flow_distance"].tolist()[:3] == [10.0, 50.0, 20.0]
    assert np.isnan(out["flow_distance"].iloc[3])


def test_profile_draws_knickzones_and_steps():
    major, _, _ = select_major_knickpoints(_knickpoints(), "major", measure="split_normalised")
    steps = filter_chi_steps(_steps())
    fig, ax = plot_knickpoints_on_profile(_segments(), major, basin=0, knickzones=_knickzones(),
                                          steps=steps)
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "knickzone (foot to lip)" in labels and "step (square sized by drop)" in labels
    assert any("outlet / mouth" in l for l in labels)
    assert "1 knickzones, 2 steps" in ax.get_title()
    plt.close(fig)
    fig, ax = plot_knickpoints_on_profile(_segments(), major, knickzones=_knickzones(), steps=steps)
    plt.close(fig)


def test_long_profile_draws_steps_as_segments_green_when_flagged():
    from lsdviztools3.charts.knickpoints import STEP_COLOUR, STEP_FLAGGED_COLOUR
    from matplotlib.colors import to_hex
    nodes = pd.DataFrame({"basin_index": [1] * 3, "branch_id": [2] * 3, "is_main_stem": [1] * 3,
                          "flow_distance": [0.0, 3000.0, 4000.0], "observed_elevation": [0.0, 120.0, 260.0]})
    fig, axes = plot_steps_on_long_profile(nodes, filter_chi_steps(_steps()))
    ax0, ax1 = axes.ravel()[:2]
    assert ax0.get_title(loc="left") == "basin 0: 2 steps"
    # no squares: the steps are line segments only
    assert not [c for c in ax1.collections if c.__class__.__name__ == "PathCollection"]
    seg = [l for l in ax1.get_lines() if to_hex(l.get_color()) == STEP_COLOUR]
    assert len(seg) == 1
    assert seg[0].get_xdata().tolist() == pytest.approx([3.0, 3.05])     # bottom, top flow distance (km)
    assert seg[0].get_ydata().tolist() == pytest.approx([120.0, 240.0])
    # basin 0: step 1 is at a tributary mouth -> green
    colours = sorted(to_hex(l.get_color()) for l in ax0.get_lines() if to_hex(l.get_color()) in
                     (STEP_COLOUR, STEP_FLAGGED_COLOUR))
    assert colours == sorted([STEP_COLOUR, STEP_FLAGGED_COLOUR])
    labels = [t.get_text() for t in ax0.get_legend().get_texts()]
    assert labels == ["step", "step at basin outlet or tributary mouth"]
    plt.close(fig)
    with pytest.raises(ValueError, match="flow_distance"):
        plot_steps_on_long_profile(nodes.drop(columns="flow_distance"), _steps())


def test_selection_order_per_channel_and_direction():
    from lsdviztools3.charts.knickpoints import knickpoint_selection_order
    o = knickpoint_selection_order(_knickpoints())
    b0 = o[(o.basin_index == 0) & (o.branch_id == 0)].sort_values(["direction", "selection_order"])
    # steepening 5 > 0.5, shallowing -4 alone; the junction row is left out
    assert list(zip(b0.direction, b0.selection_order, b0.split_delta_ksn_normalised)) == [
        ("shallowing", 1, -4.0), ("steepening", 1, 5.0), ("steepening", 2, 0.5)]
    assert o[o.basin_index == 1][["direction", "selection_order"]].values.tolist() == [["shallowing", 1]]
    assert len(knickpoint_selection_order(_knickpoints(), majors_only=True)) == 3
    assert set(knickpoint_selection_order(_knickpoints(), basins=[1]).basin_index) == {1}


def test_order_kde_rows_by_direction():
    from lsdviztools3.charts.knickpoints import knickpoint_selection_order, plot_knickpoint_order_kde
    rows = []
    rng = np.random.default_rng(1)
    for br in range(30):
        # per channel: a lip near 1000 m, a weaker one near 600 m, a foot near 300 m
        for v, z0, c0 in [(-10.0, 1000, 8), (-5.0, 600, 4), (3.0, 300, 2)]:
            r = {c: 0 for c in KNICKPOINT_COLUMNS + KNICKPOINT_SPLIT_COLUMNS}
            r.update(basin_index=br % 3, branch_id=br, split_delta_ksn_normalised=v,
                     elevation=z0 + rng.normal(0, 30), chi=c0 + rng.normal(0, 0.3))
            rows.append(r)
    o = knickpoint_selection_order(pd.DataFrame(rows))
    fig, axes = plot_knickpoint_order_kde(o, max_order=3)
    assert axes.shape == (2, 2)
    assert [t.get_text() for t in axes[0, 1].get_legend().get_texts()] == [
        "1st on its channel (n = 30)", "2nd on its channel (n = 30)"]
    assert [t.get_text() for t in axes[1, 1].get_legend().get_texts()] == ["1st on its channel (n = 30)"]
    # 1st shallowing peaks near 1000 m, 1st steepening near 300 m
    for ax, z in ((axes[0, 0], 1000), (axes[1, 0], 300)):
        line = ax.get_lines()[0]
        assert line.get_ydata()[np.nanargmax(line.get_xdata())] == pytest.approx(z, abs=40)
    plt.close(fig)
    with pytest.raises(ValueError, match="max_order"):
        plot_knickpoint_order_kde(o, max_order=9)


def test_knickzone_reach_lines_select_segments_in_chi_range():
    seg = _segments()
    reach = knickzone_reach_lines(seg, _knickzones())
    got = sorted(zip(reach["basin_index"], reach["segment_id"]))
    # basin 0 chi 1-2: no whole segment; basin 1 chi 0-5: segments 0 (0-2) and 1 (2-4)
    assert got == [(1, 0), (1, 1)]


def test_cli_knickpoints_from_stage1_with_stepaware_ksn(tmp_path):
    from lsdviztools3.cli import plot_knickpoints
    d = tmp_path / "run"
    d.mkdir()
    _knickpoints().to_csv(d / "x_chi_knickpoints_mn_0p450.csv", index=False)
    _thresholds().to_csv(d / "x_chi_knickpoint_thresholds_mn_0p450.csv", index=False)
    _knickzones().to_csv(d / "x_chi_knickzones_mn_0p450.csv", index=False)
    _segments().to_csv(d / "x_chi_fit_segments_mn_0p450.csv", index=False)
    _segments(is_step=True).to_csv(d / "x_chi_fit_segments_stepaware_mn_0p450.csv", index=False)
    _steps().to_csv(d / "x_chi_steps_mn_0p450.csv", index=False)
    out = tmp_path / "out"
    assert plot_knickpoints.main(["--data-dir", str(d), "--prefix", "x", "--m-over-n", "0.45",
                                  "--out-dir", str(out), "--min-step-drop", "10"]) == 0
    m = pd.read_csv(out / "x_major_knickpoints.csv")
    assert len(m) == 3 and (m["threshold_measure"] == "split_normalised").all()
    assert pd.read_csv(out / "x_plotted_steps.csv")["step_id"].tolist() == [0, 2]
    for f in ("knickpoint_profile", "knickpoints_by_elevation", "knickpoint_order_kde"):
        assert (out / f"x_{f}.png").exists()
    assert plot_knickpoints.main(["--data-dir", str(d), "--prefix", "x", "--m-over-n", "0.45",
                                  "--out-dir", str(out), "--stepaware-ksn", "--no-steps", "--no-knickzones"]) == 0


def test_cli_basins_limit_tables_and_knickzone_flags(tmp_path):
    from lsdviztools3.cli import plot_knickpoints
    d = tmp_path / "run"
    d.mkdir()
    _knickpoints().to_csv(d / "x_chi_knickpoints_mn_0p450.csv", index=False)
    _thresholds().to_csv(d / "x_chi_knickpoint_thresholds_mn_0p450.csv", index=False)
    _knickzones().to_csv(d / "x_chi_knickzones_mn_0p450.csv", index=False)
    _segments().to_csv(d / "x_chi_fit_segments_mn_0p450.csv", index=False)
    _steps().to_csv(d / "x_chi_steps_mn_0p450.csv", index=False)
    out = tmp_path / "out"
    assert plot_knickpoints.main(["--data-dir", str(d), "--prefix", "x", "--m-over-n", "0.45",
                                  "--out-dir", str(out), "--basins", "1", "--profile-knickzones"]) == 0
    assert set(pd.read_csv(out / "x_major_knickpoints.csv")["basin_index"]) == {1}
    assert pd.read_csv(out / "x_plotted_steps.csv")["step_id"].tolist() == [2]
    p = plot_knickpoints.build_parser().parse_args(["--data-dir", "a", "--prefix", "x", "--m-over-n", "0.45",
                                                     "--out-dir", "o"])
    assert p.basins_per_row == 5 and not p.profile_knickzones and not p.no_map_knickzones


def test_by_elevation_counts_and_basin_columns():
    from lsdviztools3.charts.knickpoints import plot_knickpoints_by_elevation
    major, _, _ = select_major_knickpoints(_knickpoints(), "major", measure="split_normalised")
    fig, (ah, asr) = plot_knickpoints_by_elevation(_knickpoints(), major)
    labels = [t.get_text() for t in ah.get_legend().get_texts()]
    assert labels == ["all segment boundaries (4)", "major, shallowing upstream (2)",
                      "major, steepening upstream (1)"]
    assert [t.get_text() for t in asr.get_xticklabels()] == ["0", "1"]
    offs = np.concatenate([c.get_offsets() for c in asr.collections if c.__class__.__name__ == "PathCollection"])
    assert sorted(offs[:, 1].tolist()) == [110.0, 120.0, 150.0]   # elevations of the 3 majors
    assert ah.get_xscale() == "log"
    plt.close(fig)
