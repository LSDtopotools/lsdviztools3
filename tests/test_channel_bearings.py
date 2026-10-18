"""Tests for the grouped channel-bearing / junction-angle readers and charts.

Synthetic CSVs follow the exact lsdtt-channel-network-analysis headers (see
io.channel_bearings): rose bins centred on 0 (bin 0 lower edge -w/2), both
angular types per group, pooled "4+" order group, text "nan".
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lsdviztools3.io import channel_bearings as CB

mpl = pytest.importorskip("matplotlib")
mpl.use("Agg")

GROUPS = ("1", "2", "3", "4+")


def _rose(width=30.0) -> pd.DataFrame:
    rows = []
    for g in GROUPS:
        for atype, period in (("directional", 360.0), ("axial", 180.0)):
            for k in range(int(period / width)):
                centre = k * width
                count = 10 + (5 if (g == "1" and k == 0) else 0)
                rows.append({"measurand": "link_bearing", "group_key": "strahler_order",
                             "group_value": g, "stratum_key": "none", "stratum_value": "all",
                             "angular_type": atype, "bin_centre_deg": centre,
                             "bin_lower_deg": centre - width / 2,
                             "bin_upper_deg": centre + width / 2, "count": count,
                             "fraction_count": count / 100.0, "length_m": 30.0 * count,
                             "fraction_length": count / 100.0})
    return pd.DataFrame(rows)


def _bearing_summary() -> pd.DataFrame:
    rows = []
    for i, g in enumerate(GROUPS):
        for atype, period in (("directional", 360), ("axial", 180)):
            rows.append({"measurand": "link_bearing", "group_key": "strahler_order",
                         "group_value": g, "stratum_key": "none", "stratum_value": "all",
                         "angular_type": atype, "period_deg": period, "n": 100 + i,
                         "n_blocks": 10 + i, "total_length_m": 1e4,
                         "circular_mean_deg_count": 5.0 * i,
                         "mean_resultant_length_count": 0.1,
                         "circular_variance_count": 0.9, "circular_sd_deg_count": 120,
                         "circular_mean_deg_length": 10.0 + 20 * i,
                         "mean_resultant_length_length": 0.2, "n_eff_length": 80,
                         "rayleigh_z": 2.0, "rayleigh_p_iid_approx": 0.02,
                         "fraction_on_d8_multiples": 0.05,
                         "meets_min_group_size": 0 if g == "4+" else 1})
    return pd.DataFrame(rows)


def _junction_summary(with_ci=False) -> pd.DataFrame:
    rows = []
    for m, base in (("alpha", 60.0), ("beta_min", 120.0), ("beta_max", 170.0)):
        for g, n, nb in (("2", 611, 100), ("3", 262, 80), ("4+", 80, 0)):
            rows.append({"measurand": m, "group_key": "strahler_order", "group_value": g,
                         "stratum_key": "none", "stratum_value": "all", "n": n,
                         "n_blocks": nb, "mean_deg": base + 1, "median_deg": base,
                         "sd_deg": 20.0, "q25_deg": base - 15, "q75_deg": base + 15,
                         "iqr_deg": 30.0, "meets_min_group_size": 1})
    df = pd.DataFrame(rows)
    if with_ci:
        df["median_ci_low_deg"] = df["median_deg"] - 3
        df["median_ci_high_deg"] = df["median_deg"] + 4
        df["n_bootstrap"] = 1000
        df.loc[df.group_value == "4+", ["median_ci_low_deg", "median_ci_high_deg"]] = np.nan
    return df


def _junctions() -> pd.DataFrame:
    rng = np.random.default_rng(1)
    orders = np.repeat([2, 3, 4, 5], [61, 26, 6, 2])
    return pd.DataFrame({"junction": range(len(orders)), "junction_order": orders,
                         "receiver_order": orders + 1,
                         "alpha_deg": rng.uniform(10, 120, len(orders)),
                         "beta_1_deg": rng.uniform(90, 180, len(orders)),
                         "beta_2_deg": rng.uniform(90, 180, len(orders))})


def _comparisons() -> pd.DataFrame:
    return pd.DataFrame({
        "measurand": "link_bearing", "stratum_key": "none", "stratum_value": "all",
        "angular_type": ["directional", "axial"], "comparison_mode": "pairwise",
        "group_a": ["1", "1"], "group_b": ["4+", "4+"], "n_a": [100, 100], "n_b": [103, 103],
        "n_blocks_a": [10, 10], "n_blocks_b": [0, 0],
        "circular_wasserstein_deg": [2.6, 1.9], "rotation_invariant": 0,
        "optimal_rotation_deg": np.nan, "kuiper_V": [0.05, 0.04],
        "kuiper_p_iid_approx": [0.0004, 0.5], "kuiper_p_holm": [0.001, 0.5],
        "n_comparisons_in_family": 6, "permutation_p": np.nan, "note": np.nan})


@pytest.fixture
def run_dir(tmp_path):
    _rose().to_csv(tmp_path / "wf_channel_bearing_rose.csv", index=False, na_rep="nan")
    _bearing_summary().to_csv(tmp_path / "wf_channel_bearing_group_summary.csv", index=False)
    _comparisons().to_csv(tmp_path / "wf_channel_bearing_comparisons.csv", index=False,
                          na_rep="nan")
    _junction_summary().to_csv(tmp_path / "wf_junction_angle_group_summary.csv", index=False)
    _junctions().to_csv(tmp_path / "wf_junction_angles.csv", index=False)
    return tmp_path


def test_load_outputs_by_header(run_dir):
    out = CB.load_channel_geometry_outputs(run_dir, "wf")
    assert out["rose"] is not None and out["junction_comparisons"] is None
    sel = CB.select_groups(out["rose"], group_key="strahler_order", angular_type="axial")
    assert set(sel.group_value.astype(str)) == set(GROUPS)
    assert CB.select_groups(out["rose"], group_key="lithology").empty


def test_rose_missing_columns_raise(tmp_path):
    p = tmp_path / "bad.csv"
    pd.DataFrame({"group_key": ["x"]}).to_csv(p, index=False)
    with pytest.raises(KeyError):
        CB.load_bearing_rose(p)


def test_order_group_mask_pools_4plus():
    o = pd.Series([2, 3, 4, 5, 6])
    assert CB.order_group_mask(o, "4+").tolist() == [False, False, True, True, True]
    assert CB.order_group_mask(o, "3").tolist() == [False, True, False, False, False]


def test_petals_wrap_bin_centred_on_zero():
    from lsdviztools3.charts.channel_bearings import _petals

    r = _rose()
    rows = r[(r.group_value == "1") & (r.angular_type == "directional")]
    theta, width, h, axial = _petals(rows, "count")
    assert not axial
    assert np.allclose(width, np.deg2rad(30.0))  # [-15, 15) is 30 deg wide
    assert theta[0] == 0.0 and h[0] == 15


def test_petals_refuse_mixed_angular_types():
    from lsdviztools3.charts.channel_bearings import _petals

    r = _rose()
    with pytest.raises(ValueError, match="mix"):
        _petals(r[r.group_value == "1"], "count")


def test_axial_rose_is_mirrored():
    from lsdviztools3.charts.channel_bearings import _petals

    r = _rose()
    rows = r[(r.group_value == "1") & (r.angular_type == "axial")]
    theta, _w, h, axial = _petals(rows, "count")
    assert axial and len(theta) == 12
    assert np.allclose(np.rad2deg(theta[6:]) - np.rad2deg(theta[:6]), 180.0)
    assert np.array_equal(h[:6], h[6:])


def test_binned_rose_draws_summary_mean_and_labels_p_with_blocks():
    from lsdviztools3.charts.channel_bearings import plot_binned_rose

    fig, ax = plot_binned_rose(_rose(), _bearing_summary(), weight="fraction_length")
    assert ax.get_theta_direction() == -1
    labels = ax.get_legend_handles_labels()[1]
    assert labels[0].startswith("strahler order 1: n=100 (10 blocks), R=0.200")
    assert all("p_iid_approx" in lbl for lbl in labels)
    assert "below min group size" in labels[3]
    means = sorted(round(float(np.rad2deg(ln.get_xdata()[0])), 6) for ln in ax.lines)
    assert means == [10.0, 30.0, 50.0, 70.0]  # length-weighted means from the summary
    mpl.pyplot.close(fig)


def test_rose_comparison_grid_two_runs():
    from lsdviztools3.charts.channel_bearings import plot_rose_comparison

    fig, axes = plot_rose_comparison({"wv": (_rose(), _bearing_summary()),
                                      "null": (_rose(), None)}, angular_type="axial",
                                     weight="fraction_count")
    assert len(axes) == 8
    assert axes[0].get_title().startswith("wv: strahler order 1")
    assert {round(a.get_rmax(), 6) for a in axes} == {round(0.15 * 1.05, 6)}  # shared radius
    mpl.pyplot.close(fig)


@pytest.mark.parametrize("with_ci", [False, True])
def test_junction_summary_boxes(with_ci):
    from lsdviztools3.charts.channel_bearings import plot_junction_angle_summary

    fig, ax = plot_junction_angle_summary(_junction_summary(with_ci), measurand="alpha")
    texts = [t.get_text() for t in ax.texts]
    assert texts == ["n=611\n100 blocks", "n=262\n80 blocks", "n=80\n0 blocks"]
    assert [t.get_text() for t in ax.get_xticklabels()] == ["2", "3", "4+"]
    has_bar = any(c.__class__.__name__ == "LineCollection" for c in ax.collections)
    assert has_bar == with_ci
    assert ("median CI" in ax.get_title()) == with_ci
    mpl.pyplot.close(fig)


def test_junction_summary_two_runs():
    from lsdviztools3.charts.channel_bearings import plot_junction_angle_summary

    fig, ax = plot_junction_angle_summary({"wv": _junction_summary(),
                                           "null": _junction_summary()})
    assert ax.get_legend_handles_labels() is not None
    assert len(ax.texts) == 6
    mpl.pyplot.close(fig)


def test_junction_distributions_use_summary_groups():
    from lsdviztools3.charts.channel_bearings import plot_junction_angle_distributions

    fig, ax = plot_junction_angle_distributions(_junctions(), _junction_summary())
    labels = ax.get_legend_handles_labels()[1]
    assert labels == ["order 2 (n=611, 100 blocks)", "order 3 (n=262, 80 blocks)",
                      "order 4+ (n=80, 0 blocks)"]
    mpl.pyplot.close(fig)


def test_comparison_table_keeps_iid_label(run_dir):
    from lsdviztools3.charts.channel_bearings import plot_comparison_table

    fig, ax = plot_comparison_table(CB.load_bearing_comparisons(
        run_dir / "wf_channel_bearing_comparisons.csv"), kind="bearing",
        angular_type="directional")
    cells = {k: c.get_text().get_text() for k, c in ax.tables[0].get_celld().items()}
    assert cells[(0, 6)] == "kuiper_p_iid_approx"
    assert cells[(1, 3)] == "103 (0)"
    assert cells[(1, 6)] == "4.0e-04"
    mpl.pyplot.close(fig)


def test_channel_bearings_cli(run_dir, tmp_path):
    from lsdviztools3.cli import plot_channel_bearings

    out = tmp_path / "figs"
    assert plot_channel_bearings.main(["--data-dir", str(run_dir), "--prefix", "wf",
                                       "--out-dir", str(out), "--dpi", "50",
                                       "--compare-prefix", "wf"]) == 0
    names = {p.name for p in out.glob("*.png")}
    for atype in ("axial", "directional"):
        for w in ("count", "length"):
            assert f"wf_rose_{atype}_{w}.png" in names
            assert f"wf_vs_wf_rose_{atype}_{w}.png" in names
    assert {"wf_junction_alpha_summary.png", "wf_junction_alpha_distributions.png",
            "wf_bearing_comparisons_directional.png",
            "wf_vs_wf_junction_alpha_summary.png"} <= names
