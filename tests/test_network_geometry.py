"""Tests for the network-geometry sweep readers and charts.

The readers all return NaNs for missing files rather than raising, which is the
right behaviour for a sweep but makes silent misreads easy. Most of these tests
therefore check that a value actually arrives when the file IS present, not
just that nothing blew up.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lsdviztools3.charts import (
    plot_delta_L_cdf, plot_junction_angle_cdf, plot_metric_panels,
    plot_metric_vs_temperature,
)
from lsdviztools3.io import network_geometry as ng


# --------------------------------------------------------------------------
# Stem parsing
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "stem, expected",
    [
        ("ocn_T00_ocn_dem", 0.0),
        ("ocn_T001_ocn_dem", 0.01),
        ("ocn_T025_ocn_dem", 0.25),
        ("ocn_T05_ocn_dem", 0.5),
    ],
)
def test_temperature_of_tag_only_handles_leading_zero_forms(stem, expected):
    # The sweep writes the temperature with the decimal point deleted. For the
    # 0.xx forms that is still invertible, so the tag alone suffices.
    assert ng.temperature_of(stem) == pytest.approx(expected)


@pytest.mark.parametrize("stem", ["ocn_T10_ocn_dem", "ocn_T20_ocn_dem",
                                  "ocn_T100_ocn_dem", "ocn_T500_ocn_dem"])
def test_temperature_of_refuses_ambiguous_tags(stem):
    # Deleting the decimal point is NOT invertible once the value is >= 1:
    # "10" is 1.0, "100" is 10.0, "20" is 2.0. Guessing here silently rescaled
    # the x-axis of every temperature plot by 10x, so the loader now returns
    # NaN and makes the caller supply the parameter file.
    assert np.isnan(ng.temperature_of(stem))


def test_temperature_of_reads_the_parameter_file(tmp_path):
    # The parameter file is authoritative, and resolves the ambiguous cases.
    (tmp_path / "ocn_T20.param").write_text(
        "write_fname: ocn_T20\ntemperature_scale: 2.0\nn_rows: 1000\n")
    assert ng.temperature_of("ocn_T20_ocn_dem", tmp_path) == pytest.approx(2.0)

    (tmp_path / "ocn_T100.param").write_text("temperature_scale: 10.0\n")
    assert ng.temperature_of("ocn_T100_ocn_dem", tmp_path) == pytest.approx(10.0)


def test_temperature_of_falls_back_when_param_absent(tmp_path):
    # Directory given but no param file: the unambiguous tag still resolves.
    assert ng.temperature_of("ocn_T025_ocn_dem", tmp_path) == pytest.approx(0.25)


def test_temperature_of_is_nan_for_non_ocn():
    assert np.isnan(ng.temperature_of("wn_uniform_random"))


@pytest.mark.parametrize(
    "stem, family",
    [
        ("ocn_T05_ocn_dem", "Optimal Channel Network"),
        ("ds_r07_diamond_square_ss", "Diamond-square fractal surface"),
        ("rn_red_noise_ss", "Red-noise surface"),
        ("wn_uniform_random_ss", "White-noise surface"),
        ("something_else", "unknown"),
    ],
)
def test_family_of(stem, family):
    assert ng.family_of(stem) == family


# --------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------

def test_load_area_scaling_missing_file_gives_nans(tmp_path):
    out = ng.load_area_scaling(tmp_path, "nope")
    assert set(out) >= {"beta", "beta_n_decades"}
    assert all(np.isnan(v) for v in out.values())


def test_load_area_scaling_reads_values(tmp_path):
    pd.DataFrame([{
        "beta": 0.47, "beta_stderr": 0.01, "beta_ols": 0.52,
        "r_squared": 0.99, "n_used": 1234, "n_decades": 1.8,
    }]).to_csv(tmp_path / "s_area_scaling_summary.csv", index=False)

    out = ng.load_area_scaling(tmp_path, "s")
    assert out["beta"] == pytest.approx(0.47)
    assert out["beta_n_decades"] == pytest.approx(1.8)
    # beta_ols is carried separately: the MLE and OLS fits disagree in general
    # and collapsing them would hide that.
    assert out["beta_ols"] == pytest.approx(0.52)


def test_load_paired_flowpaths_keeps_rejection_counts(tmp_path):
    # The counts are what distinguish "no complexity" from "no pairs", so they
    # must survive the read.
    pd.DataFrame([{
        "n_pairs": 40, "median_delta_L": 0.13, "mean_delta_L": 0.2,
        "median_delta_chi": 0.05, "n_edges_examined": 157008,
        "n_rejected_distance": 156968, "n_rejected_no_junction": 0,
    }]).to_csv(tmp_path / "s_paired_flowpaths_summary.csv", index=False)

    out = ng.load_paired_flowpaths(tmp_path, "s")
    assert out["n_pairs"] == 40
    assert out["n_rejected_distance"] == 156968
    # Derived: what fraction of examined divide edges yielded a usable pair.
    assert out["internal_divide_fraction"] == pytest.approx(40 / 157008)


def test_internal_divide_fraction_is_nan_without_edges(tmp_path):
    # A summary with no examined edges must not divide by zero and must not
    # report 0.0, which would read as "no internal divides" rather than "not
    # measured".
    pd.DataFrame([{"n_pairs": 0, "n_edges_examined": 0}]).to_csv(
        tmp_path / "s_paired_flowpaths_summary.csv", index=False)
    assert np.isnan(ng.load_paired_flowpaths(tmp_path, "s")["internal_divide_fraction"])


def test_internal_divide_fraction_zero_when_nothing_rejoins(tmp_path):
    # This case is real: an OCN at T=0 examined 824851 edges and retained no
    # pairs at all. That is a genuine 0, not a missing measurement.
    pd.DataFrame([{"n_pairs": 0, "n_edges_examined": 824851}]).to_csv(
        tmp_path / "s_paired_flowpaths_summary.csv", index=False)
    assert ng.load_paired_flowpaths(tmp_path, "s")["internal_divide_fraction"] == 0.0


def test_collect_network_sweep_orders_families(tmp_path):
    for stem in ("wn_uniform_random", "ocn_T05_ocn_dem"):
        pd.DataFrame([{"beta": 0.5, "n_decades": 2.0}]).to_csv(
            tmp_path / f"{stem}_area_scaling_summary.csv", index=False)

    table = ng.collect_network_sweep(
        tmp_path, ["wn_uniform_random", "ocn_T05_ocn_dem"],
        include_junction_angles=False,
    )
    assert list(table["stem"]) == ["wn_uniform_random", "ocn_T05_ocn_dem"]
    # family is categorical and ordered, so sorting puts OCN first even though
    # it was supplied second.
    ordered = table.sort_values("family")["stem"].tolist()
    assert ordered[0] == "ocn_T05_ocn_dem"


def test_load_junction_angle_sets_omits_missing(tmp_path):
    pd.DataFrame({"donor_opening_angle_deg": [60.0, 70.0, np.nan]}).to_csv(
        tmp_path / "good_junction_angles.csv", index=False)

    out = ng.load_junction_angle_sets(tmp_path, ["good", "absent"])
    # A stem with no file is omitted entirely rather than mapped to an empty
    # array, so a caller cannot plot an invisible curve and think it worked.
    assert set(out) == {"good"}
    assert out["good"].tolist() == [60.0, 70.0]


# --------------------------------------------------------------------------
# Charts
# --------------------------------------------------------------------------

def test_junction_angle_cdf_is_monotone_and_reaches_one():
    rng = np.random.default_rng(0)
    sets = {
        "ocn_T05_ocn_dem": rng.normal(72, 10, 300),
        "wn_uniform_random": rng.normal(95, 25, 300),
    }
    fig, ax = plot_junction_angle_cdf(sets, reference_angle=72.0)

    assert ax.get_ylim() == (0.0, 1.0)
    # axvline also lands in get_lines(), with only two points; the CDFs are the
    # long ones.
    lines = [l for l in ax.get_lines() if np.asarray(l.get_ydata()).size > 2]
    assert len(lines) == 2
    for line in lines:
        y = np.asarray(line.get_ydata())
        assert np.all(np.diff(y) >= 0)          # a CDF cannot decrease
        assert y[-1] == pytest.approx(1.0)      # and must reach 1


def test_junction_angle_cdf_marks_sparse_curves():
    fig, ax = plot_junction_angle_cdf(
        {"ocn_T05_ocn_dem": np.array([60.0, 70.0, 80.0])},
        min_count=30, reference_angle=None,
    )
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any("sparse" in t for t in labels)


def test_delta_L_cdf_runs():
    rng = np.random.default_rng(1)
    fig, ax = plot_delta_L_cdf({"ocn_T05_ocn_dem": rng.uniform(0, 2, 200)})
    assert ax.get_ylim() == (0.0, 1.0)


def test_metric_panels_skips_absent_columns():
    table = pd.DataFrame({
        "stem": ["ocn_T00_ocn_dem", "wn_uniform_random"],
        "family": ["Optimal Channel Network", "White-noise surface"],
        "beta": [0.55, 0.70],
        "median_delta_L": [0.2, 0.6],
    })
    fig, axes = plot_metric_panels(table)
    # Only the two present metrics get panels; the other defaults are dropped
    # rather than drawn empty.
    assert len(axes) == 2


def test_metric_panels_raises_when_nothing_plottable():
    table = pd.DataFrame({
        "stem": ["a"], "family": ["Optimal Channel Network"],
        "beta": [np.nan],
    })
    with pytest.raises(ValueError, match="none of the requested metric"):
        plot_metric_panels(table)


def test_metric_vs_temperature_plots_ocn_series():
    table = pd.DataFrame({
        "stem": ["ocn_T00_ocn_dem", "ocn_T05_ocn_dem", "wn_uniform_random"],
        "family": pd.Categorical(
            ["Optimal Channel Network", "Optimal Channel Network",
             "White-noise surface"],
            categories=list(ng.FAMILY_ORDER) + ["unknown"], ordered=True,
        ),
        "temperature": [0.0, 0.5, np.nan],
        "median_delta_L": [0.2, 0.35, 0.6],
    })
    fig, ax = plot_metric_vs_temperature(table, metric="median_delta_L")
    ocn_line = ax.get_lines()[0]
    assert ocn_line.get_xdata().tolist() == [0.0, 0.5]


def test_metric_vs_temperature_rejects_unknown_metric():
    table = pd.DataFrame({"stem": ["a"], "family": ["Optimal Channel Network"],
                          "temperature": [0.0]})
    with pytest.raises(KeyError):
        plot_metric_vs_temperature(table, metric="not_a_column")


def test_paired_flowpaths_found_under_pf_prefix(tmp_path):
    # The delta_L step runs under a `pf_` write prefix so its other outputs
    # cannot collide with the sweep's. The loader must find that spelling
    # WITHOUT anyone copying files into place, or a rerun leaves stale copies
    # being plotted silently.
    pd.DataFrame([{
        "n_pairs": 100, "median_delta_L": 0.42, "n_edges_examined": 1000,
    }]).to_csv(tmp_path / "pf_s_paired_flowpaths_summary.csv", index=False)

    out = ng.load_paired_flowpaths(tmp_path, "s")
    assert out["median_delta_L"] == pytest.approx(0.42)
    assert out["internal_divide_fraction"] == pytest.approx(0.1)


def test_unprefixed_paired_flowpaths_wins_over_pf(tmp_path):
    pd.DataFrame([{"n_pairs": 1, "median_delta_L": 0.11,
                   "n_edges_examined": 10}]).to_csv(
        tmp_path / "s_paired_flowpaths_summary.csv", index=False)
    pd.DataFrame([{"n_pairs": 2, "median_delta_L": 0.99,
                   "n_edges_examined": 10}]).to_csv(
        tmp_path / "pf_s_paired_flowpaths_summary.csv", index=False)

    assert ng.load_paired_flowpaths(tmp_path, "s")["median_delta_L"] == pytest.approx(0.11)


def test_summarise_edge_basins_area_weighted_median():
    from lsdviztools3.io.network_geometry import summarise_edge_basins
    # Two basins of 0.4 and one of 0.2: half the drained area lies in basins of
    # 0.4 or bigger, so the area-weighted median is 0.4 -- while the PLAIN
    # median is also 0.4 here only by coincidence of the small sample.
    s = summarise_edge_basins([0.4, 0.4, 0.2])
    assert s["n_edge_basins"] == 3
    assert s["max_fraction"] == pytest.approx(0.4)
    assert s["area_weighted_median_fraction"] == pytest.approx(0.4)
    assert s["fraction_sum"] == pytest.approx(1.0)


def test_summarise_edge_basins_ignores_the_tiny_majority():
    from lsdviztools3.io.network_geometry import summarise_edge_basins
    # One large basin plus a thousand negligible ones -- the situation in every
    # north_south run. The plain median is ~0 and says nothing; the
    # area-weighted median correctly reports that one basin holds the tile.
    f = [0.5] + [0.0005] * 1000
    s = summarise_edge_basins(f)
    assert s["median_fraction"] == pytest.approx(0.0005)
    assert s["area_weighted_median_fraction"] == pytest.approx(0.5)


def test_summarise_edge_basins_empty():
    from lsdviztools3.io.network_geometry import summarise_edge_basins
    s = summarise_edge_basins([])
    assert s["n_edge_basins"] == 0
    assert np.isnan(s["area_weighted_median_fraction"])


def test_edge_basin_cdf_is_area_weighted_by_default():
    from lsdviztools3.charts import plot_edge_basin_cdf
    fig, ax = plot_edge_basin_cdf({"ocn_T05_ocn_dem": np.array([0.5, 0.25, 0.25])})
    line = ax.get_lines()[0]
    y = np.asarray(line.get_ydata())
    # Area weighting must reach 1 and be monotone.
    assert y[-1] == pytest.approx(1.0)
    assert np.all(np.diff(y) >= 0)
    assert ax.get_xscale() == "log"
