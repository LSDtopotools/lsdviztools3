"""Tests for the charts (matplotlib x-y / distribution) plots."""

from __future__ import annotations

import numpy as np
import pytest

from lsdviztools3.charts import plot_distributions


def test_plot_distributions_returns_fig_ax():
    series = {
        "a": np.random.default_rng(0).normal(80, 15, 200),
        "b": np.random.default_rng(1).normal(95, 20, 200),
    }
    fig, ax = plot_distributions(
        series, bins=20, value_range=(0, 180),
        styles={"a": {"color": "crimson"}, "b": {"linestyle": "--"}},
        xlabel="angle", title="t",
    )
    assert ax.get_xlabel() == "angle"
    assert len(ax.get_legend().get_texts()) == 2
    assert ax.get_xlim() == (0.0, 180.0)


def test_plot_distributions_log_x_drops_nonpositive():
    # Mixed positive/zero/negative; log_x keeps only positives.
    series = {"r": np.array([0.001, 0.01, 0.1, 1.0, 0.0, -1.0])}
    fig, ax = plot_distributions(series, bins=10, value_range=(-3, 0), log_x=True)
    # 4 positive values -> histogram has area; axis spans the log range.
    assert ax.get_xlim() == (-3.0, 0.0)


def test_plot_distributions_skips_empty_series():
    series = {"empty": np.array([]), "good": np.array([1.0, 2.0, 3.0])}
    fig, ax = plot_distributions(series, bins=5)
    # Only the non-empty series is drawn / legended.
    assert len(ax.get_legend().get_texts()) == 1


def test_plot_distributions_all_empty_raises():
    with pytest.raises(ValueError, match="No plottable values"):
        plot_distributions({"a": np.array([]), "b": np.array([np.nan])}, bins=5)


# --- rose_diagram -----------------------------------------------------------

def test_rose_diagram_returns_north_up_clockwise_polar():
    from lsdviztools3.charts import rose_diagram

    rng = np.random.default_rng(0)
    fig, ax = rose_diagram(rng.uniform(0, 360, 300), bin_width_deg=15)
    assert ax.name == "polar"
    assert ax.get_theta_direction() == -1               # clockwise
    assert ax.get_theta_offset() == pytest.approx(np.pi / 2)  # North at top


def test_rose_diagram_grouped_has_legend_per_group():
    from lsdviztools3.charts import rose_diagram

    rng = np.random.default_rng(1)
    groups = {f"order {o}": rng.uniform(0, 360, 50) for o in (1, 2, 3)}
    fig, ax = rose_diagram(groups=groups, bin_width_deg=20, label_fmt="{}")
    assert len(ax.get_legend().get_texts()) == 3


def test_rose_diagram_no_data_raises():
    from lsdviztools3.charts import rose_diagram

    with pytest.raises(ValueError, match="No finite bearings"):
        rose_diagram(np.array([np.nan, np.inf]))


# --- circular statistics ----------------------------------------------------

def test_circular_wasserstein_identical_is_zero():
    from lsdviztools3.charts import circular_wasserstein

    x = np.array([10.0, 80.0, 200.0, 350.0])
    assert circular_wasserstein(x, x) == pytest.approx(0.0, abs=1e-9)


def test_circular_wasserstein_wraparound_is_small():
    from lsdviztools3.charts import circular_wasserstein

    # Samples near 350 vs near 10 are ~20 deg apart across the wrap, not ~340.
    a = np.array([350.0, 351.0, 349.0, 352.0])
    b = np.array([10.0, 11.0, 9.0, 12.0])
    d = circular_wasserstein(a, b, period=360.0)
    assert d < 40.0


def test_circular_wasserstein_symmetric():
    from lsdviztools3.charts import circular_wasserstein

    rng = np.random.default_rng(2)
    a = rng.uniform(0, 360, 100)
    b = rng.uniform(0, 360, 120)
    assert circular_wasserstein(a, b) == pytest.approx(
        circular_wasserstein(b, a), rel=1e-9)


def test_kuiper_identical_vs_disjoint():
    from lsdviztools3.charts import kuiper_two_sample

    rng = np.random.default_rng(3)
    base = rng.uniform(0, 360, 400)
    v_same, p_same = kuiper_two_sample(base, base)
    assert v_same == pytest.approx(0.0, abs=1e-9) and p_same > 0.5
    # Tightly clustered but at different bearings -> large V, tiny p.
    a = rng.normal(30, 3, 300) % 360
    b = rng.normal(200, 3, 300) % 360
    v_diff, p_diff = kuiper_two_sample(a, b)
    assert v_diff > 0.8 and p_diff < 0.01


def test_plot_circular_comparison_both_returns_two_axes():
    from lsdviztools3.charts import plot_circular_comparison

    rng = np.random.default_rng(4)
    fig, axes = plot_circular_comparison(
        rng.uniform(0, 360, 200), rng.uniform(0, 360, 200),
        labels=("X", "Y"), kind="both", annotate=True)
    assert len(axes) == 2


def _hanc_frame(n=300, seed=5):
    import pandas as pd

    rng = np.random.default_rng(seed)
    d = rng.uniform(0, 20_000, n)
    cls = rng.choice([1, 2], n)
    hanc = np.where(cls == 1, rng.uniform(0, 2, n), rng.uniform(5, 15, n))
    return pd.DataFrame({
        "channel_flow_distance_m": d, "hanc_m": hanc,
        "class": cls, "source_id": rng.choice([0, 1], n),
    })


def test_plot_hanc_distance_two_classes_legended():
    from lsdviztools3.charts import plot_hanc_distance

    fig, ax = plot_hanc_distance(_hanc_frame())
    texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert texts == ["floodplain", "terrace"]
    assert ax.get_xlabel() == "Distance along channel (km)"
    assert ax.get_ylabel() == "Height above channel (m)"


def test_plot_hanc_distance_source_filter_and_breaks():
    import pandas as pd
    from lsdviztools3.charts import plot_hanc_distance

    breaks = pd.DataFrame({
        "source_id": [0, 0], "window_centre_m": [5_000.0, 15_000.0],
        "breaks_m": ["3.5", "2.5;8.0"],
    })
    fig, ax = plot_hanc_distance(
        _hanc_frame(), source_id=0, breaks=breaks, distance_units="m")
    texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert texts == ["floodplain", "terrace", "band break"]
    assert ax.get_xlabel() == "Distance along channel (m)"


def test_plot_hanc_distance_subsamples_to_cap():
    from lsdviztools3.charts import plot_hanc_distance

    fig, ax = plot_hanc_distance(_hanc_frame(2_000), max_points=100)
    # Two scatter collections (one per class), each capped at 100 points.
    sizes = [len(c.get_offsets()) for c in ax.collections]
    assert sizes == [100, 100]


def test_plot_hanc_distance_empty_raises():
    from lsdviztools3.charts import plot_hanc_distance

    df = _hanc_frame(50)
    with pytest.raises(ValueError, match="No plottable pixels"):
        plot_hanc_distance(df[df["class"] == 99])
