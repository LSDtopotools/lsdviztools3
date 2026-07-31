"""Tests for the lsdtt-divide-migration IO, charts and CLI.

The emphasis is on the three things a divide-migration figure can silently get
wrong: the A − B sign convention, NaN leaking in as a value, and treating the
``ttest`` p-value column as if it were a dispersion.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lsdviztools3.charts.divide_migration import (
    aggressor_annotation, plot_head_metric_vs_distance, plot_pair_comparison,
)
from lsdviztools3.io.divide_migration import (
    DIVIDE_METRICS, UNCERTAINTY_METHODS, find_outputs, load_pair_summary,
    resolve_metric, uncertainty_method,
)


def _summary(method: str = "std_err", **overrides) -> pd.DataFrame:
    """A two-pair summary; pair 2|3 is the strong one under every method."""
    rows = [
        {
            "pair_id": "1|2", "unit_a": 1, "unit_b": 2,
            "n_heads_a": 25, "n_heads_b": 17,
            "mean_head_elevation_a": 243.4, "mean_head_elevation_b": 251.1,
            "delta_head_elevation": -7.6,
            "n_head_elevation_a": 25, "n_head_elevation_b": 17,
            "mean_gradient_a": 0.0416, "mean_gradient_b": 0.0439,
            "delta_gradient": -0.0023, "n_gradient_a": 25, "n_gradient_b": 17,
            "mean_relief_a": 70.1, "mean_relief_b": 55.0, "delta_relief": 15.1,
            "n_relief_a": 22, "n_relief_b": 15,
            "mean_chi_a": 5.98, "mean_chi_b": 5.73, "delta_chi": 0.25,
            "n_chi_a": 25, "n_chi_b": 17,
            "uncertainty_method": method,
            "uncertainty_head_elevation": 39.2, "uncertainty_gradient": 0.0056,
            "uncertainty_relief": 13.8, "uncertainty_chi": 0.85,
        },
        {
            "pair_id": "2|3", "unit_a": 2, "unit_b": 3,
            "n_heads_a": 13, "n_heads_b": 35,
            "mean_head_elevation_a": 245.7, "mean_head_elevation_b": 252.5,
            "delta_head_elevation": -6.8,
            "n_head_elevation_a": 13, "n_head_elevation_b": 35,
            "mean_gradient_a": 0.0525, "mean_gradient_b": 0.0338,
            "delta_gradient": 0.0187, "n_gradient_a": 13, "n_gradient_b": 35,
            "mean_relief_a": 72.1, "mean_relief_b": 57.0, "delta_relief": 15.1,
            "n_relief_a": 13, "n_relief_b": 35,
            "mean_chi_a": 6.42, "mean_chi_b": 6.06, "delta_chi": 0.36,
            "n_chi_a": 13, "n_chi_b": 35,
            "uncertainty_method": method,
            "uncertainty_head_elevation": 40.5, "uncertainty_gradient": 0.0054,
            "uncertainty_relief": 14.9, "uncertainty_chi": 1.21,
        },
    ]
    df = pd.DataFrame(rows)
    for key, value in overrides.items():
        df[key] = value
    return df


def _heads() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    n = 40
    df = pd.DataFrame({
        "head_id": np.arange(1, n + 1),
        "unit_id": np.repeat([1, 2], n // 2),
        "pair_id_nearest": ["1|2"] * n,
        "head_elevation_m": rng.normal(250, 40, n),
        "chi": rng.normal(6, 1, n),
        "distance_to_divide_m": rng.uniform(100, 3000, n),
        "relief_m": rng.normal(60, 20, n),
        "gradient_to_divide": rng.normal(0.04, 0.01, n),
    })
    # The NaN cases the tool genuinely produces: nearest divide point on
    # nodata, and a head sitting on its divide (zero distance -> NaN gradient).
    df.loc[0:2, "relief_m"] = np.nan
    df.loc[3, "distance_to_divide_m"] = np.nan
    return df


# --------------------------------------------------------------------------
# IO / metadata
# --------------------------------------------------------------------------

def test_metric_aliases_resolve_to_the_same_metric():
    assert resolve_metric("relief") is resolve_metric("relief_m")
    assert resolve_metric("gradient_to_divide") is DIVIDE_METRICS["gradient"]
    assert resolve_metric("head_elevation_m").head_column == "head_elevation_m"
    with pytest.raises(KeyError):
        resolve_metric("not_a_metric")


def test_aggressor_sign_convention_matches_the_cli_reference():
    # Lower heads / lower chi implicate A; higher gradient / relief implicate A.
    assert DIVIDE_METRICS["head_elevation"].aggressor_delta_sign == -1
    assert DIVIDE_METRICS["chi"].aggressor_delta_sign == -1
    assert DIVIDE_METRICS["gradient"].aggressor_delta_sign == +1
    assert DIVIDE_METRICS["relief"].aggressor_delta_sign == +1


def test_uncertainty_method_read_from_the_file():
    assert uncertainty_method(_summary("bootstrap")).key == "bootstrap"
    # ttest is not a dispersion and supports no error bar.
    ttest = UNCERTAINTY_METHODS["ttest"]
    assert ttest.is_dispersion is False
    # std_dev is a dispersion but carries no significance test.
    assert UNCERTAINTY_METHODS["std_dev"].significance_multiplier is None
    # An unrecorded method claims nothing.
    unknown = uncertainty_method(pd.DataFrame({"pair_id": ["1|2"]}))
    assert unknown.significance_multiplier is None


def test_load_pair_summary_coerces_empty_fields_to_nan(tmp_path):
    csv = tmp_path / "x_divide_pair_summary.csv"
    csv.write_text(
        "pair_id,unit_a,unit_b,mean_chi_a,mean_chi_b,delta_chi,"
        "uncertainty_chi,uncertainty_method\n"
        "1|2,1,2,5.0,4.0,1.0,0.2,std_err\n"
        "2|3,2,3,,,,,std_err\n"
    )
    df = load_pair_summary(csv)
    assert df["pair_id"].tolist() == ["1|2", "2|3"]
    assert np.isnan(df.loc[1, "delta_chi"])
    assert df["delta_chi"].dtype.kind == "f"


def test_find_outputs_from_stem_and_directory(tmp_path):
    for name in (
        "sc_divides.fgb", "sc_channel_heads.fgb",
        "sc_divide_migration_metrics.csv", "sc_divide_pair_summary.csv",
    ):
        (tmp_path / name).write_text("")
    by_stem = find_outputs(tmp_path / "sc")
    by_dir = find_outputs(tmp_path)
    assert by_stem == by_dir
    assert by_stem["divides"].name == "sc_divides.fgb"
    assert by_stem["pair_summary"].name == "sc_divide_pair_summary.csv"

    with pytest.raises(FileNotFoundError):
        find_outputs(tmp_path / "not_a_run")


def test_find_outputs_reports_missing_optional_outputs(tmp_path):
    (tmp_path / "sc_divides.fgb").write_text("")
    outputs = find_outputs(tmp_path / "sc")
    assert outputs["divides"] is not None
    assert outputs["pair_summary"] is None
    assert outputs["metrics_csv"] is None


# --------------------------------------------------------------------------
# Aggressor annotation: sign + significance semantics
# --------------------------------------------------------------------------

def test_aggressor_side_follows_the_delta_sign_per_metric():
    row = _summary().iloc[1]  # pair 2|3, delta_gradient > 0 -> A (unit 2)
    text, state = aggressor_annotation(
        row, DIVIDE_METRICS["gradient"], UNCERTAINTY_METHODS["std_err"]
    )
    assert "unit 2 (A)" in text
    assert state == "strong"  # |0.0187| > 2 * 0.0054

    # Same pair, chi delta is positive -> B is the aggressor for chi.
    text, _ = aggressor_annotation(
        row, DIVIDE_METRICS["chi"], UNCERTAINTY_METHODS["std_err"]
    )
    assert "unit 3 (B)" in text


def test_weak_delta_is_flagged_not_silently_asserted():
    row = _summary().iloc[0]  # delta_chi 0.25 vs std_err 0.85
    text, state = aggressor_annotation(
        row, DIVIDE_METRICS["chi"], UNCERTAINTY_METHODS["std_err"]
    )
    assert state == "weak"
    assert "n.s." in text


def test_std_dev_claims_no_significance():
    row = _summary("std_dev").iloc[1]
    text, state = aggressor_annotation(
        row, DIVIDE_METRICS["gradient"], UNCERTAINTY_METHODS["std_dev"]
    )
    assert state == "untested"
    assert "not tested" in text


def test_ttest_column_is_read_as_a_p_value():
    ttest = UNCERTAINTY_METHODS["ttest"]
    row = _summary("ttest", uncertainty_gradient=0.0017).iloc[1]
    text, state = aggressor_annotation(row, DIVIDE_METRICS["gradient"], ttest)
    assert state == "strong"
    assert "p = 0.0017" in text

    row = _summary("ttest", uncertainty_gradient=0.42).iloc[1]
    _, state = aggressor_annotation(row, DIVIDE_METRICS["gradient"], ttest)
    assert state == "weak"


def test_nan_delta_is_reported_not_treated_as_zero():
    row = _summary(delta_chi=np.nan).iloc[0]
    text, state = aggressor_annotation(
        row, DIVIDE_METRICS["chi"], UNCERTAINTY_METHODS["std_err"]
    )
    assert state == "nodata"
    assert "NaN" in text


# --------------------------------------------------------------------------
# Charts
# --------------------------------------------------------------------------

def test_plot_pair_comparison_labels_both_sides_and_the_sign():
    fig, (ax_means, ax_delta) = plot_pair_comparison(_summary(), metric="chi")
    labels = [t.get_text() for t in ax_means.get_yticklabels()]
    assert any("A=1" in label and "B=2" in label for label in labels)
    assert "A − B" in ax_delta.get_xlabel()
    assert "standard error" in ax_delta.get_xlabel()
    # Both side means are drawn, joined by a bar, for both pairs.
    assert len(ax_means.lines) >= 2 * 3


def test_plot_pair_comparison_draws_error_bars_for_dispersions():
    fig, (_, ax_delta) = plot_pair_comparison(_summary("bootstrap"), metric="relief")
    assert ax_delta.containers, "expected an errorbar container for a dispersion"
    assert "95% CI half-width" in ax_delta.get_xlabel()


def test_plot_pair_comparison_never_draws_an_error_bar_for_ttest():
    fig, (_, ax_delta) = plot_pair_comparison(_summary("ttest"), metric="gradient")
    assert not ax_delta.containers, "a p-value must not be drawn as an error bar"
    assert "no error bar" in ax_delta.get_xlabel()
    assert "p-value" in ax_delta.get_xlabel()
    annotations = [t.get_text() for t in ax_delta.texts]
    assert any("p = " in text for text in annotations)


def test_plot_pair_comparison_drops_all_nan_pairs_and_says_so():
    summary = _summary()
    for column in ("mean_relief_a", "mean_relief_b", "delta_relief"):
        summary.loc[0, column] = np.nan
    fig, (ax_means, _) = plot_pair_comparison(summary, metric="relief")
    labels = [t.get_text() for t in ax_means.get_yticklabels()]
    assert len(labels) == 1 and "2|3" in labels[0]
    caption = " ".join(t.get_text() for t in fig.texts)
    assert "dropped" in caption


def test_plot_pair_comparison_keeps_a_one_sided_pair():
    summary = _summary()
    summary.loc[0, ["mean_chi_b", "n_chi_b", "delta_chi", "uncertainty_chi"]] = np.nan
    fig, (ax_means, _) = plot_pair_comparison(summary, metric="chi")
    assert len(ax_means.get_yticklabels()) == 2
    assert any("no valid values" in t.get_text() for t in ax_means.texts)


def test_plot_pair_comparison_all_nan_raises():
    summary = _summary()
    for column in ("mean_chi_a", "mean_chi_b", "delta_chi"):
        summary[column] = np.nan
    with pytest.raises(ValueError, match="finite"):
        plot_pair_comparison(summary, metric="chi")


def test_plot_pair_comparison_missing_columns_raises():
    with pytest.raises(KeyError):
        plot_pair_comparison(pd.DataFrame({"pair_id": ["1|2"]}), metric="chi")


def test_head_scatter_drops_nan_and_reports_the_count():
    fig, ax = plot_head_metric_vs_distance(_heads(), metric="relief")
    legend = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any("NaN dropped" in text for text in legend)
    plotted = sum(len(line.get_xdata()) for line in ax.lines)
    assert plotted == 40 - 4  # 3 NaN relief + 1 NaN distance
    caption = " ".join(t.get_text() for t in fig.texts)
    assert "dropped" in caption


def test_head_scatter_can_restrict_to_one_pair():
    heads = _heads()
    heads.loc[heads["unit_id"] == 2, "pair_id_nearest"] = "2|3"
    fig, ax = plot_head_metric_vs_distance(heads, metric="chi", pair_id="1|2")
    assert len(ax.lines) == 1


def test_head_scatter_all_nan_raises():
    heads = _heads()
    heads["chi"] = np.nan
    with pytest.raises(ValueError, match="nothing to plot"):
        plot_head_metric_vs_distance(heads, metric="chi")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def test_cli_help_exits_cleanly():
    from lsdviztools3.cli import plot_divide_migration

    with pytest.raises(SystemExit) as exc:
        plot_divide_migration.build_parser().parse_args(["--help"])
    assert exc.value.code == 0


def test_cli_main_smoke(monkeypatch, tmp_path):
    from lsdviztools3.cli import plot_divide_migration

    calls = {}

    def fake_render(prefix, out_dir, **kwargs):
        calls["prefix"] = prefix
        calls["out_dir"] = out_dir
        calls["kwargs"] = kwargs
        return {"map": "map.png"}

    monkeypatch.setattr(
        plot_divide_migration, "render_divide_migration_figures", fake_render
    )
    code = plot_divide_migration.main([
        "--prefix", "run/sc",
        "--outdir", str(tmp_path / "figs"),
        "--base", "dem.tif",
        "--metric", "relief_m",
        "--point-size", "0.1c",
    ])
    assert code == 0
    assert calls["prefix"] == "run/sc"
    assert calls["kwargs"]["base_raster"] == "dem.tif"
    assert calls["kwargs"]["metric"] == "relief_m"
    assert calls["kwargs"]["point_size"] == "0.1c"
    # The shared --cmap default ("geo") is overridden for a metric palette.
    assert calls["kwargs"]["cmap"] == "roma"
