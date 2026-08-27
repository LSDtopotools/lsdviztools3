"""Tests for the charts (matplotlib x-y / distribution) plots."""

from __future__ import annotations

import numpy as np
import pandas as pd
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


def _cosmo_sweep_frame():
    rows = []
    for lat in (-60.0, -30.0, 0.0, 30.0, 60.0):
        for nuclide in ("Be10", "Al26", "C14"):
            for scaling in ("stone", "lsdn"):
                for muon in ("braucherborchers", "tunedexp"):
                    rows.append(
                        {
                            "run_id": f"R_{lat}_{nuclide}_{scaling}_{muon}",
                            "lat_deg": lat,
                            "pressure_hpa": 900.0,
                            "erate_gcm2yr": 0.01,
                            "nuclide": nuclide,
                            "scaling_scheme": scaling,
                            "muon_scheme": muon,
                            "pct_diff_vs_cronus_v3": 0.05 * lat + (1.2 if scaling == "lsdn" else 0.0),
                            "pct_diff_vs_cronuscalc": -0.03 * lat,
                        }
                    )
    return pd.DataFrame(rows)


def test_plot_cosmo_sweep_panels_returns_grid():
    from lsdviztools3.charts import plot_cosmo_sweep_panels

    fig, axes = plot_cosmo_sweep_panels(_cosmo_sweep_frame(), reference="cronus_v3")
    assert axes.shape == (3, 2)
    assert "Residual vs cronus_v3" in axes[0, 0].get_ylabel()


def test_plot_cosmo_sample_residuals_grouped_boxplot():
    from lsdviztools3.charts import plot_cosmo_sample_residuals

    df = pd.DataFrame(
        {
            "dataset": ["dethier", "dethier", "clapp", "clapp"],
            "scaling_scheme": ["stone", "lsdn", "stone", "lsdn"],
            "pct_diff_lsdtt_minus_reference": [-5.2, -1.1, -3.8, 0.6],
        }
    )
    fig, ax = plot_cosmo_sample_residuals(df)
    assert ax.get_ylabel() == "Residual (%)"
    assert len(ax.get_xticklabels()) == 2


def test_plot_ksn_vs_erate_with_errorbars():
    from lsdviztools3.charts import plot_ksn_vs_erate

    df = pd.DataFrame(
        {
            "basin_index": [1, 2, 3],
            "inferred_erate_mm_yr": [0.2, 0.5, 1.0],
            "erate_err_low_mm_yr": [0.03, 0.07, 0.1],
            "erate_err_high_mm_yr": [0.04, 0.08, 0.12],
            "ksn_median": [110.0, 180.0, 260.0],
            "ksn_q25": [90.0, 140.0, 210.0],
            "ksn_q75": [135.0, 230.0, 310.0],
        }
    )
    fig, ax = plot_ksn_vs_erate(df)
    assert ax.get_xscale() == "log"
    assert ax.get_yscale() == "log"
    assert "erosion rate" in ax.get_xlabel().lower()


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


def _hand_frame(n=300, seed=5):
    import pandas as pd

    rng = np.random.default_rng(seed)
    d = rng.uniform(0, 20_000, n)
    cls = rng.choice([1, 2], n)
    hand = np.where(cls == 1, rng.uniform(0, 2, n), rng.uniform(5, 15, n))
    return pd.DataFrame({
        "channel_flow_distance_m": d, "hand_m": hand,
        "class": cls, "source_id": rng.choice([0, 1], n),
    })


def test_plot_hand_distance_two_classes_legended():
    from lsdviztools3.charts import plot_hand_distance

    fig, ax = plot_hand_distance(_hand_frame())
    texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert texts == ["floodplain", "terrace"]
    assert ax.get_xlabel() == "Distance along channel (km)"
    assert ax.get_ylabel() == "Height above channel (m)"


def test_plot_hand_distance_source_filter_and_breaks():
    import pandas as pd
    from lsdviztools3.charts import plot_hand_distance

    breaks = pd.DataFrame({
        "source_id": [0, 0], "window_centre_m": [5_000.0, 15_000.0],
        "breaks_m": ["3.5", "2.5;8.0"],
    })
    fig, ax = plot_hand_distance(
        _hand_frame(), source_id=0, breaks=breaks, distance_units="m")
    texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert texts == ["floodplain", "terrace", "band break"]
    assert ax.get_xlabel() == "Distance along channel (m)"


def test_plot_hand_distance_subsamples_to_cap():
    from lsdviztools3.charts import plot_hand_distance

    fig, ax = plot_hand_distance(_hand_frame(2_000), max_points=100)
    # Two scatter collections (one per class), each capped at 100 points.
    sizes = [len(c.get_offsets()) for c in ax.collections]
    assert sizes == [100, 100]


def test_plot_hand_distance_empty_raises():
    from lsdviztools3.charts import plot_hand_distance

    df = _hand_frame(50)
    with pytest.raises(ValueError, match="No plottable pixels"):
        plot_hand_distance(df[df["class"] == 99])


def _window_hist_frame():
    import pandas as pd

    rows = []
    rng = np.random.default_rng(7)
    for w in (1000.0, 3000.0, 5000.0):
        bins = np.logspace(-1, 2, 60)
        counts = (1000 * np.exp(-0.5 * ((np.log10(bins) - 0.5) / 0.2) ** 2)
                  + rng.integers(0, 5, 60))
        for b, c in zip(bins, counts):
            rows.append({"source_id": 0, "window_centre_m": w,
                         "hand_m": b, "n_pixels": int(c)})
    return pd.DataFrame(rows)


def test_plot_hand_window_populations_panels_and_breaks():
    import pandas as pd
    from lsdviztools3.charts import plot_hand_window_populations

    breaks = pd.DataFrame({
        "source_id": [0, 0, 0],
        "window_centre_m": [1000.0, 3000.0, 5000.0],
        "break_m": [5.0, 7.5, 6.0],
        "own_break": [True, False, True],
    })
    fig, axes = plot_hand_window_populations(
        _window_hist_frame(), breaks_df=breaks, source_id=0)
    assert len(axes) == 3
    # Each panel carries a break line legend; inherited window is marked.
    legends = [ax.get_legend() for ax in axes]
    assert all(lg is not None for lg in legends)
    assert "(inherited)" in legends[1].get_texts()[0].get_text()


def test_plot_hand_window_populations_missing_source_raises():
    from lsdviztools3.charts import plot_hand_window_populations

    with pytest.raises(ValueError, match="No HAND window histograms"):
        plot_hand_window_populations(_window_hist_frame(), source_id=9)


# --- plot_valley_cross_sections ---------------------------------------------

def _cross_section_frame(n_stations=8, n_samples=41):
    import pandas as pd

    rows = []
    offs = np.linspace(-500.0, 500.0, n_samples)
    for st in range(n_stations):
        lo, ro = -200.0, 220.0            # valley sides (offset metres)
        for off in offs:
            in_valley = 1 if lo <= off <= ro else 0
            elev = 100.0 + 0.02 * abs(off) + (0.0 if in_valley else 30.0)
            rows.append({
                "source_id": 0, "station_seq": st,
                "along_distance_m": 500.0 * st, "offset_m": off,
                "x": off, "y": 500.0 * st, "elevation": elev,
                "in_valley": in_valley,
                "left_offset_m": lo, "right_offset_m": ro,
                "valley_width_m": ro - lo,
            })
    return pd.DataFrame(rows)


def test_plot_valley_cross_sections_panel_per_station():
    from lsdviztools3.charts import plot_valley_cross_sections

    fig, axes = plot_valley_cross_sections(_cross_section_frame(6), source_id=0)
    assert len(axes) == 6
    # Panel title reports the along-distance and the valley width (ro - lo).
    assert "w=420 m" in axes[0].get_title()


def test_plot_valley_cross_sections_subsamples_to_cap():
    from lsdviztools3.charts import plot_valley_cross_sections

    fig, axes = plot_valley_cross_sections(
        _cross_section_frame(50), source_id=0, max_panels=12)
    assert len(axes) <= 12


def test_plot_valley_cross_sections_missing_source_raises():
    from lsdviztools3.charts import plot_valley_cross_sections

    with pytest.raises(ValueError, match="No cross-section stations"):
        plot_valley_cross_sections(_cross_section_frame(4), source_id=9)


def test_combine_width_steps_long_wall_not_oversegmented():
    from lsdviztools3.charts import combine_width_steps

    # Smooth widening wall: old cumulative logic split every ~2 m; we want one
    # merged step unless there are prominent widening events.
    stage = np.arange(312.5, 330.0, 0.1)
    width = 1000.0 + 6.0 * (stage - stage.min())
    boundaries, mids = combine_width_steps(
        stage, width, ds=30.0, min_step_height_m=2.0, min_step_width_px=3)

    assert len(boundaries) == 2
    assert len(mids) == 1


def test_combine_width_steps_keeps_major_widening_jumps():
    from lsdviztools3.charts import combine_width_steps

    stage = np.arange(300.0, 330.0, 0.1)
    # Two sharp widening events separated by broad near-plateaux.
    jump1 = 550.0 / (1.0 + np.exp(-(stage - 308.0) / 0.09))
    jump2 = 650.0 / (1.0 + np.exp(-(stage - 320.0) / 0.09))
    width = 700.0 + jump1 + jump2

    boundaries, _ = combine_width_steps(
        stage, width, ds=30.0, min_step_height_m=2.0, min_step_width_px=3)

    assert len(boundaries) >= 3
    assert any(abs(b - 308.0) < 1.0 for b in boundaries)
    assert any(abs(b - 320.0) < 1.0 for b in boundaries)


def test_optimize_step_midpoint_reduces_crossings_within_bounds():
    from lsdviztools3.charts.cross_sections import (
        _water_line_crossings, _optimize_step_midpoint,
    )

    off = np.linspace(-200.0, 200.0, 2001)
    # Concave section with a central hump: some levels give 4 crossings,
    # nearby levels give 2.
    elev = 100.0 + 0.0008 * off**2 + 1.8 * np.exp(-(off / 20.0) ** 2)
    zlo, zhi = 100.95, 102.0
    zmid = 0.5 * (zlo + zhi)

    n0 = len(_water_line_crossings(off, elev, zmid))
    zopt = _optimize_step_midpoint(off, elev, zmid, zlo, zhi)
    n1 = len(_water_line_crossings(off, elev, zopt))

    assert zlo < zopt < zhi
    assert n1 < n0


def test_adjust_step_midpoints_for_crossings_preserves_count_and_order():
    from lsdviztools3.charts.cross_sections import (
        _adjust_step_midpoints_for_crossings,
    )

    off = np.linspace(-200.0, 200.0, 2001)
    elev = 100.0 + 0.0008 * off**2 + 1.8 * np.exp(-(off / 20.0) ** 2)
    boundaries = [100.9, 102.0, 104.0]
    mids = [101.45, 103.0]

    adj = _adjust_step_midpoints_for_crossings(off, elev, boundaries, mids)

    assert len(adj) == 2
    assert boundaries[0] < adj[0] < boundaries[1]
    assert boundaries[1] < adj[1] < boundaries[2]


# --- plot_valley_break_metric -----------------------------------------------

def _break_metric_frame():
    import pandas as pd

    rows = []
    bins = np.logspace(-1.3, 2, 128)

    def add(scope, wc, break_bin, boundary_bin, peak_bins, raw=True):
        # Bimodal smoothed metric: a small valley-floor mode + a big wall mode.
        sm = (300 * np.exp(-0.5 * ((np.log10(bins) - 0.2) / 0.25) ** 2)
              + 1000 * np.exp(-0.5 * ((np.log10(bins) - 1.3) / 0.3) ** 2))
        for i, b in enumerate(bins):
            rows.append({
                "source_id": 0, "scope": scope, "window_centre_m": wc,
                "hand_m": b,
                "raw_count": int(sm[i]) if raw else 0,
                "smoothed": float(sm[i]),
                "is_peak": 1 if i in peak_bins else 0,
                "is_boundary": 1 if i == boundary_bin else 0,
                "is_break": 1 if i == break_bin else 0,
            })

    # Global: no clean break (reach-adaptive used).
    add("global", -1.0, break_bin=-1, boundary_bin=-1, peak_bins=[])
    # Three windows with a dip at bin 60 or none: the applied break in the test
    # then makes reach 1000 an accepted "own" dip, 3000 "inherited" (no own
    # dip), and 5000 a "median" (own dip replaced by a different applied break).
    add("window", 1000.0, break_bin=60, boundary_bin=60, peak_bins=[30, 100],
        raw=False)
    add("window", 3000.0, break_bin=-1, boundary_bin=-1, peak_bins=[], raw=False)
    add("window", 5000.0, break_bin=60, boundary_bin=60, peak_bins=[30, 100],
        raw=False)
    return pd.DataFrame(rows)


def test_plot_valley_break_metric_panels_and_applied_break():
    import pandas as pd
    from lsdviztools3.charts import plot_valley_break_metric

    dip60 = float(np.logspace(-1.3, 2, 128)[60])   # the dip HAND at bin 60
    breaks = pd.DataFrame({
        "source_id": [0, 0, 0],
        "window_centre_m": [1000.0, 3000.0, 5000.0],
        # own accepted (== dip), inherited, median (!= dip).
        "break_m": [dip60, 2.0, dip60 + 3.0],
        "own_break": [True, False, True],
    })
    fig, axes = plot_valley_break_metric(
        _break_metric_frame(), breaks_df=breaks, source_id=0)
    # Global panel + 3 window panels.
    assert len(axes) == 4
    # Global panel notes the reach-adaptive fallback.
    assert "reach-adaptive" in axes[0].get_title()
    # Window panels are tagged own / inherited / median.
    assert "(own)" in axes[1].get_title()
    assert "(inherited)" in axes[2].get_title()
    assert "(median)" in axes[3].get_title()

    # The descriptive break line appears under each reach title.
    def texts(ax):
        return " ".join(t.get_text() for t in ax.texts)
    assert "accepted dip at" in texts(axes[1])
    assert "inherited break at" in texts(axes[2])
    assert "accepted median at" in texts(axes[3])

    # Axis is LINEAR (Simon: log space was misleading) and cropped from 0 to
    # the decision region (not the full HAND range).
    assert axes[1].get_xscale() == "linear"
    xlo, xhi = axes[1].get_xlim()
    assert xlo == 0.0
    assert xhi <= 60.5


def test_plot_valley_break_metric_missing_source_raises():
    from lsdviztools3.charts import plot_valley_break_metric

    with pytest.raises(ValueError, match="No break-metric rows"):
        plot_valley_break_metric(_break_metric_frame(), source_id=9)


# --- plot_valley_break_fraction (relief-adaptive) ---------------------------

def _relief_breaks_frame():
    import pandas as pd
    return pd.DataFrame({
        "source_id": [0, 0, 0, 0],
        "window_centre_m": [1000.0, 2000.0, 3000.0, 4000.0],
        "break_m": [2.0, 3.0, 5.0, 2.0],
        "own_break": ["true", "false", "true", "true"],
        "break_fraction": [0.55, 0.60, 0.66, 0.58],
        "relief_m": [5.0, 8.0, 15.0, 6.0],
    })


def test_plot_valley_break_fraction_panels():
    from lsdviztools3.charts import plot_valley_break_fraction

    fig, axs = plot_valley_break_fraction(
        _relief_breaks_frame(), source_id=0, frac_min=0.5, frac_max=0.7)
    # Three stacked panels: fraction, relief, break.
    assert len(axs) == 3
    assert "acceptance" in axs[0].get_ylabel()
    assert "relief" in axs[1].get_ylabel()


def test_plot_valley_break_fraction_requires_column():
    import pandas as pd
    from lsdviztools3.charts import plot_valley_break_fraction

    # A fixed-method breaks table has no break_fraction column.
    fixed = pd.DataFrame({
        "source_id": [0, 0],
        "window_centre_m": [1000.0, 2000.0],
        "break_m": [2.0, 3.0],
        "own_break": ["true", "true"],
    })
    with pytest.raises(ValueError, match="break_fraction"):
        plot_valley_break_fraction(fixed, source_id=0)


# --- plot_valley_cross_section_map ------------------------------------------

def _nodes_frame(n=10):
    import pandas as pd
    return pd.DataFrame({
        "source_id": [0] * n,
        "node_seq": list(range(n)),
        "x": np.linspace(-400.0, 400.0, n),
        "y": np.linspace(0.0, 4500.0, n),
    })


def _tiny_hillshade_tif(tmp_path):
    """A small GeoTIFF whose extent covers the synthetic cross-section frame."""
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    # Data spans x in [-500,500], y in [0,4500]; give a generous margin.
    res = 20.0
    x0, y_top = -700.0, 5200.0
    width, height = 70, 280
    transform = from_origin(x0, y_top, res, res)
    path = tmp_path / "hs.tif"
    data = np.linspace(0, 255, width * height).reshape(height, width).astype(
        "float32")
    with rasterio.open(
        path, "w", driver="GTiff", height=height, width=width, count=1,
        dtype="float32", transform=transform, nodata=-9999.0,
    ) as dst:
        dst.write(data, 1)
    return str(path)


def test_plot_valley_cross_section_map_returns_fig_ax(tmp_path):
    from lsdviztools3.charts import plot_valley_cross_section_map

    hs = _tiny_hillshade_tif(tmp_path)
    fig, ax = plot_valley_cross_section_map(
        _cross_section_frame(8), _nodes_frame(), hs, source_id=0)
    # Hillshade + transect lines + centreline are drawn; numbers are labelled.
    assert len(ax.lines) >= 8
    texts = [t.get_text() for t in ax.texts]
    assert "0" in texts and "7" in texts


def test_plot_valley_cross_section_map_label_every_thins_labels(tmp_path):
    from lsdviztools3.charts import plot_valley_cross_section_map

    hs = _tiny_hillshade_tif(tmp_path)
    fig, ax = plot_valley_cross_section_map(
        _cross_section_frame(8), _nodes_frame(), hs, source_id=0,
        label_every=4)
    # Stations 0 and 4 labelled; 1,2,3 not.
    texts = [t.get_text() for t in ax.texts]
    assert "0" in texts and "4" in texts
    assert "1" not in texts


def test_plot_valley_cross_section_map_missing_source_raises(tmp_path):
    from lsdviztools3.charts import plot_valley_cross_section_map

    hs = _tiny_hillshade_tif(tmp_path)
    with pytest.raises(ValueError, match="No cross-section stations"):
        plot_valley_cross_section_map(
            _cross_section_frame(4), _nodes_frame(), hs, source_id=9)


def test_plot_valley_cross_section_features_map_with_mask_and_boundaries(tmp_path):
    import pandas as pd
    pytest.importorskip("rasterio")
    from lsdviztools3.charts import plot_valley_cross_section_features_map

    hs = _tiny_hillshade_tif(tmp_path)
    # Valley mask encoded as source_id + 1.
    vm = np.zeros((280, 70), dtype=float)
    vm[80:220, 20:50] = 1.0
    _write_geotiff(tmp_path / "valley_mask.tif", vm)

    # One channel feature + two hydraulic boundaries (with step labels).
    feats = pd.DataFrame({
        "source_id": [0],
        "feature_type": ["channel"],
        "x": [0.0],
        "y": [2000.0],
    })
    bounds = pd.DataFrame({
        "source_id": [0, 0],
        "station_seq": [2, 3],
        "boundary_type": ["hydraulic_edge", "hydraulic_edge"],
        "step": [0, 1],
        "x": [-60.0, 65.0],
        "y": [1000.0, 1500.0],
    })

    fig, ax = plot_valley_cross_section_features_map(
        feats,
        _cross_section_frame(8),
        _nodes_frame(),
        hs,
        source_id=0,
        boundaries_df=bounds,
        valley_mask_path=str(tmp_path / "valley_mask.tif"),
        hydraulic_label_every=1,
        label_every=0,
    )
    assert "features" in ax.get_title()
    assert len(ax.texts) >= 2


# --- plot_valley_reach_map --------------------------------------------------

def _write_geotiff(path, arr, nodata=None):
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin
    arr = np.asarray(arr, dtype="float32")
    tr = from_origin(0.0, arr.shape[0], 1.0, 1.0)  # 1 m cells, origin top-left
    with rasterio.open(
        path, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
        count=1, dtype="float32", transform=tr, nodata=nodata,
    ) as dst:
        dst.write(arr, 1)


def test_plot_valley_reach_map_returns_fig_ax(tmp_path):
    pytest.importorskip("rasterio")
    import pandas as pd
    from lsdviztools3.charts import plot_valley_reach_map
    # 8x8 reach raster: three reaches in bands, negative = outside corridor.
    reach = np.full((8, 8), -1.0)
    reach[1:7, 1:3] = 0
    reach[1:7, 3:5] = 1
    reach[1:7, 5:7] = 2
    _write_geotiff(tmp_path / "reach.tif", reach)
    _write_geotiff(tmp_path / "hs.tif", np.random.rand(8, 8) * 255)
    nodes = pd.DataFrame({"source_id": [0, 0, 0],
                          "x": [2.0, 4.0, 6.0], "y": [4.0, 4.0, 4.0]})
    fig, ax = plot_valley_reach_map(
        str(tmp_path / "reach.tif"), str(tmp_path / "hs.tif"),
        nodes_df=nodes, source_id=0)
    assert "channel 0" in ax.get_title()
    # Three reach-number labels (0,1,2) + the centreline line drawn.
    labels = {t.get_text() for t in ax.texts}
    assert {"0", "1", "2"} <= labels
    assert len(ax.get_lines()) >= 1


def test_plot_valley_reach_map_all_nodata_raises(tmp_path):
    pytest.importorskip("rasterio")
    from lsdviztools3.charts import plot_valley_reach_map
    _write_geotiff(tmp_path / "reach.tif", np.full((6, 6), -1.0))
    _write_geotiff(tmp_path / "hs.tif", np.ones((6, 6)))
    with pytest.raises(ValueError, match="no valid"):
        plot_valley_reach_map(str(tmp_path / "reach.tif"),
                              str(tmp_path / "hs.tif"))


# --- plot_valley_floor_map --------------------------------------------------

def test_plot_valley_floor_map_returns_fig_ax(tmp_path):
    pytest.importorskip("rasterio")
    import pandas as pd
    from lsdviztools3.charts import plot_valley_floor_map
    # 8x8 valley mask: source 0 (=1) and source 1 (=2) ribbons; 0 = not valley.
    vm = np.zeros((8, 8))
    vm[2:6, 1:3] = 1
    vm[2:6, 5:7] = 2
    _write_geotiff(tmp_path / "vm.tif", vm)
    _write_geotiff(tmp_path / "hs.tif", np.random.rand(8, 8) * 255)
    nodes = pd.DataFrame({"source_id": [0, 0, 1, 1],
                          "x": [2.0, 2.0, 6.0, 6.0],
                          "y": [3.0, 5.0, 3.0, 5.0]})
    fig, ax = plot_valley_floor_map(
        str(tmp_path / "vm.tif"), str(tmp_path / "hs.tif"),
        nodes_df=nodes, title="t")
    assert ax.get_title() == "t"
    # Two hillshade+overlay images and the two centrelines drawn.
    assert len(ax.get_images()) == 2
    assert len(ax.get_lines()) >= 2


def test_plot_valley_floor_map_empty_raises(tmp_path):
    pytest.importorskip("rasterio")
    from lsdviztools3.charts import plot_valley_floor_map
    _write_geotiff(tmp_path / "vm.tif", np.zeros((6, 6)))
    _write_geotiff(tmp_path / "hs.tif", np.ones((6, 6)))
    with pytest.raises(ValueError, match="no valley"):
        plot_valley_floor_map(str(tmp_path / "vm.tif"), str(tmp_path / "hs.tif"))


# --- hydraulic breaks -------------------------------------------------------

def _hydraulic_breaks_frame(n_stations=30, bimodal=True):
    """Synthetic hydraulic-break table: break 1 near HAND 2 m, and (when
    bimodal) break 2 near HAND 9 m, so a 10-section bin has two clear levels."""
    import pandas as pd
    rng = np.random.default_rng(0)
    rows = []
    for st in range(n_stations):
        along = 24.0 * st
        rows.append({"source_id": 0, "station_seq": st,
                     "along_distance_m": along, "break_index": 0,
                     "stage_m": 100.0, "hand_m": 0.0})
        rows.append({"source_id": 0, "station_seq": st,
                     "along_distance_m": along, "break_index": 1,
                     "stage_m": 102.0, "hand_m": 2.0 + rng.normal(0, 0.2)})
        if bimodal and st % 3 != 0:
            rows.append({"source_id": 0, "station_seq": st,
                         "along_distance_m": along, "break_index": 2,
                         "stage_m": 109.0, "hand_m": 9.0 + rng.normal(0, 0.3)})
    return pd.DataFrame(rows)


def test_plot_hydraulic_breaks_downstream_series_per_rank():
    from lsdviztools3.charts import plot_hydraulic_breaks_downstream
    fig, ax = plot_hydraulic_breaks_downstream(
        _hydraulic_breaks_frame(), source_id=0)
    # A points series + a connecting line per break rank (ranks 1 and 2).
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "break 1" in labels and "break 2" in labels
    assert ax.get_ylabel().startswith("break height")


def test_plot_hydraulic_breaks_downstream_missing_source_raises():
    from lsdviztools3.charts import plot_hydraulic_breaks_downstream
    with pytest.raises(ValueError, match="No hydraulic breaks"):
        plot_hydraulic_breaks_downstream(
            _hydraulic_breaks_frame(), source_id=9)


def test_hydraulic_break_clusters_bimodal_finds_two():
    from lsdviztools3.charts import hydraulic_break_clusters
    cl = hydraulic_break_clusters(
        _hydraulic_breaks_frame(bimodal=True), source_id=0, bin_size=10)
    assert (cl["n_clusters"] == 2).all()
    # Cluster centres bracket the two synthetic levels (~2 m and ~9 m).
    for centres in cl["cluster_centres"]:
        assert min(centres) < 4.0 and max(centres) > 7.0


def test_hydraulic_break_clusters_unimodal_finds_one():
    from lsdviztools3.charts import hydraulic_break_clusters
    df = _hydraulic_breaks_frame(bimodal=False)
    cl = hydraulic_break_clusters(df, source_id=0, bin_size=10)
    assert (cl["n_clusters"] == 1).all()


def test_hydraulic_break_clusters_columns_and_bins():
    from lsdviztools3.charts import hydraulic_break_clusters
    cl = hydraulic_break_clusters(
        _hydraulic_breaks_frame(n_stations=25), source_id=0, bin_size=10)
    assert list(cl["n_sections"]) == [10, 10, 5]      # 25 sections -> 3 bins
    for col in ("bin_index", "bin_centre_along_m", "along_lo_m", "along_hi_m",
                "n_values", "n_clusters", "cluster_centres", "silhouette"):
        assert col in cl.columns


def test_plot_hydraulic_break_clusters_overlays_levels():
    from lsdviztools3.charts import plot_hydraulic_break_clusters
    fig, ax = plot_hydraulic_break_clusters(
        _hydraulic_breaks_frame(), source_id=0, bin_size=10)
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "reach cluster level" in labels


def test_plot_cross_section_informed_valley_mask(tmp_path):
    import pandas as pd
    pytest.importorskip("rasterio")
    from lsdviztools3.charts import plot_cross_section_informed_valley_mask
    # Binary informed floodplain (value 1); original partly overlaps it, so the
    # figure has original-only, informed-only and overlap regions.
    inf = np.zeros((10, 10)); inf[3:7, 2:8] = 1
    orig = np.zeros((10, 10)); orig[2:6, 4:9] = 1
    _write_geotiff(tmp_path / "inf.tif", inf)
    _write_geotiff(tmp_path / "orig.tif", orig)
    _write_geotiff(tmp_path / "hs.tif", np.random.rand(10, 10) * 255)
    nodes = pd.DataFrame({"source_id": [0, 0, 0],
                          "x": [3.0, 5.0, 7.0], "y": [5.0, 5.0, 5.0]})
    fig, ax = plot_cross_section_informed_valley_mask(
        str(tmp_path / "inf.tif"), str(tmp_path / "orig.tif"),
        str(tmp_path / "hs.tif"), nodes_df=nodes)
    # hillshade + original overlay + informed overlay = 3 images; centreline drawn.
    assert len(ax.get_images()) == 3
    assert len(ax.get_lines()) >= 1


def test_plot_cross_section_informed_valley_mask_empty_raises(tmp_path):
    pytest.importorskip("rasterio")
    from lsdviztools3.charts import plot_cross_section_informed_valley_mask
    _write_geotiff(tmp_path / "inf.tif", np.zeros((6, 6)))
    _write_geotiff(tmp_path / "orig.tif", np.zeros((6, 6)))
    _write_geotiff(tmp_path / "hs.tif", np.ones((6, 6)))
    with pytest.raises(ValueError, match="no valley"):
        plot_cross_section_informed_valley_mask(
            str(tmp_path / "inf.tif"), str(tmp_path / "orig.tif"),
            str(tmp_path / "hs.tif"))


# --- chi profiles -----------------------------------------------------------
#
# These pin the lsdtt-chi-analysis column contract and the stacked-plot
# conventions documented in docs/documentation.md §5.1. Frames are built with
# the exact headers lsdtt3 writes, so a header change on the C++ side (or a
# drifting loader) fails here rather than silently producing a wrong plot.

#: Header of ``*_chi_fit_nodes_mn_*.csv``, in the order lsdtt3 writes it.
CHI_NODE_HEADER = (
    "basin_index", "outlet_flat", "flat", "receiver_flat", "chi",
    "observed_elevation", "fitted_elevation", "residual_elevation",
    "abs_residual_elevation", "branch_id", "is_main_stem", "is_source",
    "is_junction", "is_split_point",
)

#: Header of ``*_chi_fit_segments_mn_*.csv``, in the order lsdtt3 writes it.
CHI_SEGMENT_HEADER = (
    "basin_index", "outlet_flat", "segment_id", "parent_segment_id",
    "branch_id", "downstream_flat", "upstream_flat", "downstream_chi",
    "upstream_chi", "downstream_elevation", "upstream_elevation", "k_sn",
    "max_abs_residual", "n_nodes",
)


def _chi_nodes_frame():
    """Two basins with known chi spans, carrying the full lsdtt3 node header.

    Basin 0 spans chi 0..10 (one branch, two segments at k_sn 20 then 5);
    basin 1 spans chi 0..5 (one segment at k_sn 50). Elevation is the
    piecewise-linear integral of those k_sn, so observed == fitted.
    """
    import pandas as pd

    rows = []
    for basin, outlet, chi_max in ((0, 100, 10), (1, 200, 5)):
        for i in range(chi_max + 1):
            chi = float(i)
            if basin == 0:
                z = 20.0 * chi if chi <= 5 else 100.0 + 5.0 * (chi - 5.0)
            else:
                z = 50.0 * chi
            rows.append({
                "basin_index": basin, "outlet_flat": outlet,
                "flat": outlet + i, "receiver_flat": outlet + max(i - 1, 0),
                "chi": chi, "observed_elevation": z, "fitted_elevation": z,
                "residual_elevation": 0.0, "abs_residual_elevation": 0.0,
                "branch_id": 0, "is_main_stem": 1,
                "is_source": int(i == chi_max), "is_junction": 0,
                "is_split_point": 0,
            })
    return pd.DataFrame(rows, columns=list(CHI_NODE_HEADER))


def _chi_segments_frame():
    """Segments matching :func:`_chi_nodes_frame`, full lsdtt3 segment header."""
    import pandas as pd

    specs = [
        # basin, outlet, seg, d_chi, u_chi, d_z,   u_z,   k_sn
        (0, 100, 0, 0.0, 5.0, 0.0, 100.0, 20.0),
        (0, 100, 1, 5.0, 10.0, 100.0, 125.0, 5.0),
        (1, 200, 0, 0.0, 5.0, 0.0, 250.0, 50.0),
    ]
    rows = [{
        "basin_index": b, "outlet_flat": o, "segment_id": s,
        "parent_segment_id": s - 1, "branch_id": 0,
        "downstream_flat": o, "upstream_flat": o + 1,
        "downstream_chi": dc, "upstream_chi": uc,
        "downstream_elevation": dz, "upstream_elevation": uz,
        "k_sn": k, "max_abs_residual": 0.0, "n_nodes": 6,
    } for b, o, s, dc, uc, dz, uz, k in specs]
    return pd.DataFrame(rows, columns=list(CHI_SEGMENT_HEADER))


def _chi_data_map_frame():
    """Chi data-map points: note ``basin_key``, not ``basin_index``."""
    import pandas as pd

    rows = []
    for key, chi_max in ((0, 10), (1, 5)):
        for i in range(chi_max + 1):
            rows.append({
                "chi": float(i), "drainage_area": 1e6 - 1e4 * i,
                "flow_distance": 30.0 * i, "elevation": 10.0 * i,
                "source_key": key, "basin_key": key,
            })
    return pd.DataFrame(rows)


def test_chi_io_column_constants_match_lsdtt3_header():
    """The loader's declared columns are the contract documented in §5.1.1."""
    from lsdviztools3.io.chi import FIT_NODE_COLUMNS, FIT_SEGMENT_COLUMNS

    assert FIT_NODE_COLUMNS == CHI_NODE_HEADER
    assert FIT_SEGMENT_COLUMNS == CHI_SEGMENT_HEADER
    # k_sn lives on the segments only; it is grafted onto nodes by
    # assign_node_ksn. Documented because it surprises people.
    assert "k_sn" in FIT_SEGMENT_COLUMNS
    assert "k_sn" not in FIT_NODE_COLUMNS


def test_assign_node_ksn_uses_the_containing_segment():
    from lsdviztools3.io.chi import assign_node_ksn

    out = assign_node_ksn(_chi_nodes_frame(), _chi_segments_frame())
    b0 = out[out["basin_index"] == 0].set_index("chi")["k_sn"]
    assert b0[1.0] == pytest.approx(20.0)     # inside segment 0
    assert b0[7.0] == pytest.approx(5.0)      # inside segment 1
    assert out[out["basin_index"] == 1]["k_sn"].eq(50.0).all()
    assert out["k_sn"].notna().all()


def test_assign_node_ksn_falls_back_to_nearest_segment():
    """A node outside every interval takes the nearest segment in its branch."""
    import pandas as pd
    from lsdviztools3.io.chi import assign_node_ksn

    nodes = _chi_nodes_frame()
    stray = nodes[nodes["basin_index"] == 0].iloc[[0]].copy()
    stray["chi"] = 12.0                        # beyond the upstream end
    nodes = pd.concat([nodes, stray], ignore_index=True)

    out = assign_node_ksn(nodes, _chi_segments_frame())
    assert out.loc[out["chi"] == 12.0, "k_sn"].iloc[0] == pytest.approx(5.0)


def test_assign_node_ksn_leaves_nan_when_branch_has_no_segments():
    from lsdviztools3.io.chi import assign_node_ksn

    nodes = _chi_nodes_frame()
    nodes.loc[nodes["basin_index"] == 1, "branch_id"] = 7   # no such segment
    out = assign_node_ksn(nodes, _chi_segments_frame())
    assert out.loc[out["basin_index"] == 1, "k_sn"].isna().all()
    assert out.loc[out["basin_index"] == 0, "k_sn"].notna().all()


def test_load_chi_fit_reads_csvs_and_attaches_ksn(tmp_path):
    from lsdviztools3.io.chi import load_chi_fit

    nodes_csv = tmp_path / "x_chi_fit_nodes_mn_0p450.csv"
    segs_csv = tmp_path / "x_chi_fit_segments_mn_0p450.csv"
    _chi_nodes_frame().to_csv(nodes_csv, index=False)
    _chi_segments_frame().to_csv(segs_csv, index=False)

    out = load_chi_fit(nodes_csv, segs_csv)
    assert list(out.columns)[:len(CHI_NODE_HEADER)] == list(CHI_NODE_HEADER)
    assert out["k_sn"].notna().all()


def test_plot_chi_stacked_nodes_scatters_observed_with_fit():
    from lsdviztools3.charts import plot_chi_stacked
    from lsdviztools3.io.chi import assign_node_ksn

    nodes = assign_node_ksn(_chi_nodes_frame(), _chi_segments_frame())
    fig, ax = plot_chi_stacked(nodes)
    assert len(ax.collections) == 2                  # one scatter per basin
    assert len(ax.lines) == 2                        # one fitted branch each
    assert ax.get_xlabel().endswith("offset per basin)")


def test_plot_chi_stacked_segments_draws_one_line_collection():
    from matplotlib.collections import LineCollection
    from lsdviztools3.charts import plot_chi_stacked

    fig, ax = plot_chi_stacked(_chi_segments_frame())
    lcs = [c for c in ax.collections if isinstance(c, LineCollection)]
    assert len(lcs) == 1
    assert len(lcs[0].get_segments()) == 3           # all segments, both basins
    assert not ax.lines                              # segments *are* the fit


def test_plot_chi_stacked_default_offset_is_20_percent_of_widest_span():
    """Documented rule: gap = 0.2 * widest basin chi span; stacks are contiguous."""
    from matplotlib.collections import LineCollection
    from lsdviztools3.charts import plot_chi_stacked

    fig, ax = plot_chi_stacked(_chi_segments_frame())
    lc = [c for c in ax.collections if isinstance(c, LineCollection)][0]
    xs = np.array([[p[0] for p in seg] for seg in lc.get_segments()])

    # Basin 0 spans 0..10 (the widest), so the gap is 2.0 and basin 1 starts
    # at 10 + 2 = 12 and ends at 17.
    assert xs[:2].min() == pytest.approx(0.0)
    assert xs[:2].max() == pytest.approx(10.0)
    assert xs[2].min() == pytest.approx(12.0)
    assert xs[2].max() == pytest.approx(17.0)


def test_plot_chi_stacked_basins_argument_selects_and_orders():
    from matplotlib.collections import LineCollection
    from lsdviztools3.charts import plot_chi_stacked

    fig, ax = plot_chi_stacked(_chi_segments_frame(), basins=[1])
    lc = [c for c in ax.collections if isinstance(c, LineCollection)][0]
    xs = np.array([[p[0] for p in seg] for seg in lc.get_segments()])
    # Only basin 1, and it is normalised back to chi = 0.
    assert len(lc.get_segments()) == 1
    assert xs.min() == pytest.approx(0.0) and xs.max() == pytest.approx(5.0)


def test_plot_chi_stacked_labels_carry_basin_index():
    from lsdviztools3.charts import plot_chi_stacked

    fig, ax = plot_chi_stacked(_chi_segments_frame())
    assert {t.get_text() for t in ax.texts} == {"0", "1"}

    fig, ax = plot_chi_stacked(_chi_segments_frame(), label_basins=False)
    assert not ax.texts


def test_colour_values_floor_small_ksn_in_log_space():
    """log10 colouring floors k_sn < 0.1 to 0.0 rather than dropping it."""
    from lsdviztools3.charts.chi import _colour_values

    out = _colour_values(np.array([0.05, 0.0, -1.0, 10.0, 100.0]), True)
    assert out[0] == pytest.approx(0.0)      # floored, NOT log10(0.05)
    assert out[1] == pytest.approx(0.0)
    assert out[2] == pytest.approx(0.0)
    assert out[3] == pytest.approx(1.0)
    assert out[4] == pytest.approx(2.0)
    # Without log scaling the values pass through untouched.
    raw = _colour_values(np.array([0.05, 10.0]), False)
    assert raw.tolist() == [0.05, 10.0]


def test_plot_chi_stacked_rejects_a_data_map_frame_as_data():
    """Data maps carry basin_key, so they are an `observed=` input only."""
    from lsdviztools3.charts import plot_chi_stacked

    with pytest.raises(KeyError, match="basin_index"):
        plot_chi_stacked(_chi_data_map_frame())


def test_plot_chi_stacked_observed_overlay_uses_the_same_offset():
    from lsdviztools3.charts import plot_chi_stacked

    fig, ax = plot_chi_stacked(_chi_segments_frame(),
                               observed=_chi_data_map_frame())
    # Observed scatters are drawn per basin, under the segments (zorder 0).
    scatters = [c for c in ax.collections if c.get_zorder() == 0]
    assert len(scatters) == 2
    x_basin1 = scatters[1].get_offsets()[:, 0]
    # Basin 1's observed chi 0..5 is shifted by the same 12.0 as its segments.
    assert x_basin1.min() == pytest.approx(12.0)
    assert x_basin1.max() == pytest.approx(17.0)


def test_plot_chi_profile_selects_the_requested_basin():
    from lsdviztools3.charts import plot_chi_profile
    from lsdviztools3.io.chi import assign_node_ksn

    nodes = assign_node_ksn(_chi_nodes_frame(), _chi_segments_frame())
    fig, ax = plot_chi_profile(nodes, basin=1)
    assert ax.get_title() == "Basin 1"
    xs = ax.collections[0].get_offsets()[:, 0]
    # True chi, no offset applied.
    assert xs.min() == pytest.approx(0.0) and xs.max() == pytest.approx(5.0)


def test_plot_chi_profile_defaults_to_the_first_basin():
    from lsdviztools3.charts import plot_chi_profile

    fig, ax = plot_chi_profile(_chi_segments_frame())
    assert ax.get_title() == "Basin 0"


def test_plot_chi_profile_unknown_basin_raises():
    from lsdviztools3.charts import plot_chi_profile

    with pytest.raises(ValueError, match="basin_index=9"):
        plot_chi_profile(_chi_segments_frame(), basin=9)


# --- cosmo transient --------------------------------------------------------


def _default_cycle_colours(n):
    """The first n colours of the active (tab10) property cycle, as RGB."""
    import matplotlib as mpl

    cycle = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
    return [mpl.colors.to_rgb(c) for c in cycle[:n]]


def _line_colours(ax):
    import matplotlib as mpl

    return [mpl.colors.to_rgb(line.get_color()) for line in ax.lines]


def _cosmo_transient_frame(n=25):
    t = np.linspace(0.0, 50000.0, n)
    # Synthetic stand-ins for the lsdtt3-computed columns; the charts only draw.
    frac = {"Be10": 1.0 - np.exp(-t / 8000.0),
            "Al26": 1.0 - np.exp(-t / 6000.0),
            "C14": 1.0 - np.exp(-t / 1500.0)}
    df = pd.DataFrame({
        "run_label": ["s5_5x_final0p1"] * n,
        "production_scaling": 5.0,
        "erosion_ratio": 5.0,
        "eff_erate_old_g_cm2_yr": 0.02,
        "eff_erate_new_g_cm2_yr": 0.1,
        "time_since_step_yr": t,
        "app_erate_ratio_C14_Be10": 1.0 + 0.8 * np.exp(-t / 12000.0) * (t / 3000.0),
        "app_erate_ratio_C14_Al26": 1.0 + 0.5 * np.exp(-t / 12000.0) * (t / 3000.0),
    })
    for nuc, f in frac.items():
        df[f"approach_frac_{nuc}"] = f
        df[f"norm_conc_{nuc}"] = 1.0 - 0.8 * f
        df[f"norm_full_{nuc}"] = 1.0 - f
    return df


def test_cosmo_transient_title_exact_mathtext():
    from lsdviztools3.charts import cosmo_transient_title

    assert cosmo_transient_title(s_t=5, erosion_ratio=5, final_eps=0.1) == (
        r"$S_t$ = 5, erosion ratio = 5, "
        r"final $\epsilon$ = 0.1 g cm$^{-2}$yr$^{-1}$"
    )


def test_plot_cosmo_transient_concentration_three_series_and_style():
    from lsdviztools3.charts import plot_cosmo_transient_concentration

    fig, ax = plot_cosmo_transient_concentration(_cosmo_transient_frame())
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels == [r"$^{10}Be$", r"$^{26}Al$", r"$^{14}C$"]
    assert ax.get_xlim() == (0.0, 10000.0)
    assert ax.get_xlabel() == "Time since step change (yrs)"
    assert ax.get_ylabel() == "Normalised CRN\nconcentration (dimensionless)"
    assert fig.get_size_inches().tolist() == [5.0, 3.0]
    assert all(line.get_linewidth() == 2.0 for line in ax.lines)
    # Default tab10 cycle taken in series order, colours never set explicitly.
    assert _line_colours(ax) == _default_cycle_colours(3)
    # Deliberate deviation from the house convention: no gridlines.
    assert not any(g.get_visible() for g in ax.xaxis.get_gridlines())
    assert not any(g.get_visible() for g in ax.yaxis.get_gridlines())
    assert all(sp.get_visible() for sp in ax.spines.values())
    assert ax.get_title() == (
        r"$S_t$ = 5, erosion ratio = 5, "
        r"final $\epsilon$ = 0.1 g cm$^{-2}$yr$^{-1}$"
    )


def test_plot_cosmo_transient_concentration_defaults_to_approach_frac():
    from lsdviztools3.charts import plot_cosmo_transient_concentration

    df = _cosmo_transient_frame()
    fig, ax = plot_cosmo_transient_concentration(df)
    y = ax.lines[0].get_ydata()
    # approach_frac_ rises from 0 towards 1 (the other prefixes fall from 1).
    assert y[0] == pytest.approx(0.0)
    assert y[-1] > y[0]
    np.testing.assert_allclose(y, df["approach_frac_Be10"].to_numpy())


def test_plot_cosmo_transient_concentration_prefix_is_selectable():
    from lsdviztools3.charts import plot_cosmo_transient_concentration

    df = _cosmo_transient_frame()
    fig, ax = plot_cosmo_transient_concentration(df, value_prefix="norm_full_")
    np.testing.assert_allclose(ax.lines[0].get_ydata(),
                               df["norm_full_Be10"].to_numpy())


def test_plot_cosmo_transient_concentration_adopts_an_axes():
    import matplotlib.pyplot as plt
    from lsdviztools3.charts import plot_cosmo_transient_concentration

    fig0, ax0 = plt.subplots(figsize=(4.0, 2.0))
    fig, ax = plot_cosmo_transient_concentration(_cosmo_transient_frame(), ax=ax0)
    assert ax is ax0 and fig is fig0
    assert fig.get_size_inches().tolist() == [4.0, 2.0]


def test_plot_cosmo_transient_erate_ratio_two_series():
    from lsdviztools3.charts import plot_cosmo_transient_erate_ratio

    fig, ax = plot_cosmo_transient_erate_ratio(_cosmo_transient_frame())
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels == [r"$^{14}C$ to $^{10}Be$", r"$^{14}C$ to $^{26}Al$"]
    assert ax.get_xlim() == (0.0, 50000.0)
    assert ax.get_ylabel() == "Ratio of apparent\nerosion rates (dimensionless)"
    assert _line_colours(ax) == _default_cycle_colours(2)
    assert ax.get_legend().get_texts()[0].get_fontsize() == 16


def test_plot_cosmo_transient_show_title_false():
    from lsdviztools3.charts import plot_cosmo_transient_erate_ratio

    fig, ax = plot_cosmo_transient_erate_ratio(_cosmo_transient_frame(),
                                               show_title=False)
    assert ax.get_title() == ""


def test_plot_cosmo_transient_missing_columns_raises():
    from lsdviztools3.charts import plot_cosmo_transient_concentration

    df = _cosmo_transient_frame().drop(columns=["approach_frac_C14"])
    with pytest.raises(KeyError, match="approach_frac_C14"):
        plot_cosmo_transient_concentration(df)


def test_plot_cosmo_transient_erate_ratio_missing_columns_raises():
    from lsdviztools3.charts import plot_cosmo_transient_erate_ratio

    df = _cosmo_transient_frame().drop(columns=["app_erate_ratio_C14_Al26"])
    with pytest.raises(KeyError, match="app_erate_ratio_C14_Al26"):
        plot_cosmo_transient_erate_ratio(df)


def test_plot_cosmo_transient_all_nan_raises():
    from lsdviztools3.charts import plot_cosmo_transient_erate_ratio

    df = _cosmo_transient_frame()
    df["app_erate_ratio_C14_Be10"] = np.nan
    df["app_erate_ratio_C14_Al26"] = "NotAttempted"
    with pytest.raises(ValueError, match="No finite apparent erosion-rate"):
        plot_cosmo_transient_erate_ratio(df)


def test_plot_cosmo_transient_mixed_run_labels_raises():
    from lsdviztools3.charts import plot_cosmo_transient_concentration

    df = pd.concat([_cosmo_transient_frame(), _cosmo_transient_frame()],
                   ignore_index=True)
    df.loc[df.index[-5:], "run_label"] = "s12_10x_final0p2"
    with pytest.raises(ValueError, match="single run_label"):
        plot_cosmo_transient_concentration(df)


def test_load_cosmo_transient_table_roundtrip(tmp_path):
    from lsdviztools3.charts import load_cosmo_transient_table

    csv = tmp_path / "cosmo_transient_sweep_s5_5x_final0p1.csv"
    _cosmo_transient_frame().to_csv(csv, index=False)
    df = load_cosmo_transient_table(csv)
    assert "approach_frac_Be10" in df.columns
    assert len(df) == 25


# ---------------------------------------------------------------------------
# Spatial (per-sample) transient cosmogenic scatter charts
# ---------------------------------------------------------------------------

def _joined_samples_frame(n=30):
    """A joined sample table with every column the spatial charts need."""
    rng = np.random.default_rng(11)
    steptime_yr = np.linspace(0.0, 500000.0, n)
    return pd.DataFrame({
        "latitude": np.linspace(25.1, 25.3, n),
        "longitude": np.linspace(92.0, 92.2, n),
        "sc_steptime_yr": steptime_yr,
        "sc_steptime_kyr": steptime_yr / 1000.0,
        "upstream_mean_eff_erate_g_cm2_yr": np.linspace(0.0265, 0.13, n),
        "upstream_mean_erate_mm_yr": np.linspace(0.1, 0.49, n),
        "app_erate_ratio_C14_Be10": 1.0 + 0.3 * rng.random(n),
        "app_erate_ratio_C14_Al26": 1.0 + 0.2 * rng.random(n),
    })


def test_plot_ratio_vs_steptime_encodes_erate_as_colour_and_size():
    from lsdviztools3.charts import plot_ratio_vs_steptime
    import matplotlib.pyplot as plt

    df = _joined_samples_frame()
    fig, ax = plot_ratio_vs_steptime(df)
    try:
        coll = ax.collections[0]
        # Colour and size must both track the true erosion rate -- size is
        # redundant with colour by design, so a flat size array is a bug.
        assert np.allclose(coll.get_array(),
                           df["upstream_mean_erate_mm_yr"].to_numpy())
        sizes = coll.get_sizes()
        assert sizes.size == len(df)
        assert sizes.max() > sizes.min()
        assert "knickpoint" in ax.get_xlabel()
        assert "^{14}C" in ax.get_ylabel()
    finally:
        plt.close(fig)


def test_plot_ratio_vs_mean_erate_uses_kyr_for_colour():
    from lsdviztools3.charts import plot_ratio_vs_mean_erate
    import matplotlib.pyplot as plt

    df = _joined_samples_frame()
    fig, ax = plot_ratio_vs_mean_erate(df)
    try:
        coll = ax.collections[0]
        assert np.allclose(coll.get_array(), df["sc_steptime_kyr"].to_numpy())
        assert np.allclose(coll.get_offsets()[:, 0],
                           df["upstream_mean_erate_mm_yr"].to_numpy())
        assert "Average upstream erosion" in ax.get_xlabel()
    finally:
        plt.close(fig)


def test_spatial_charts_drop_only_non_finite_rows():
    from lsdviztools3.charts import plot_ratio_vs_steptime
    import matplotlib.pyplot as plt

    df = _joined_samples_frame()
    df.loc[df.index[3], "app_erate_ratio_C14_Be10"] = np.nan
    df.loc[df.index[7], "upstream_mean_erate_mm_yr"] = np.nan
    fig, ax = plot_ratio_vs_steptime(df)
    try:
        # A row missing any encoded value cannot be drawn faithfully; it is
        # dropped rather than given an invented colour or size.
        assert ax.collections[0].get_offsets().shape[0] == len(df) - 2
    finally:
        plt.close(fig)


def test_spatial_charts_reject_unknown_ratio_column():
    from lsdviztools3.charts import plot_ratio_vs_steptime

    with pytest.raises(KeyError, match="Unknown ratio column"):
        plot_ratio_vs_steptime(_joined_samples_frame(),
                               ratio_col="app_erate_ratio_Be10_C14")


def test_spatial_charts_report_missing_columns():
    from lsdviztools3.charts import plot_ratio_vs_mean_erate

    df = _joined_samples_frame().drop(columns=["sc_steptime_kyr"])
    with pytest.raises(KeyError, match="sc_steptime_kyr"):
        plot_ratio_vs_mean_erate(df)


def test_load_joined_samples_roundtrip(tmp_path):
    from lsdviztools3.charts import load_joined_samples

    csv = tmp_path / "mega_sc5x_joined_samples.csv"
    _joined_samples_frame().to_csv(csv, index=False)
    df = load_joined_samples(csv)
    assert "app_erate_ratio_C14_Be10" in df.columns
    assert len(df) == 30


# ---------------------------------------------------------------------------
# Family D: spatial samples overlaid with the 1-D column sweep
# ---------------------------------------------------------------------------

def _sweep_frame(n=200):
    """A 1-D step-change sweep table, in the column tool's own schema."""
    t = np.linspace(0.0, 500000.0, n)
    # Shape is irrelevant to these tests; only the columns and ordering are.
    return pd.DataFrame({
        "run_label": ["matched_5x"] * n,
        "time_since_step_yr": t,
        "app_erate_ratio_C14_Be10": 1.0 + 0.5 * np.exp(-t / 20000.0),
        "app_erate_ratio_C14_Al26": 1.0 + 0.4 * np.exp(-t / 20000.0),
    })


def test_overlay_draws_samples_and_sweep_on_one_axes():
    from lsdviztools3.charts import plot_ratio_overlay_1d
    import matplotlib.pyplot as plt

    df = _joined_samples_frame()
    sweep = _sweep_frame()
    fig, ax = plot_ratio_overlay_1d(df, [("matched", sweep)])
    try:
        # The samples keep Family B's encoding ...
        coll = ax.collections[0]
        assert np.allclose(coll.get_array(),
                           df["upstream_mean_erate_mm_yr"].to_numpy())
        # ... and the sweep is a line, not a second scatter.
        assert len(ax.lines) == 1
        xs, ys = ax.lines[0].get_xdata(), ax.lines[0].get_ydata()
        assert np.allclose(xs, sweep["time_since_step_yr"].to_numpy())
        assert np.allclose(ys, sweep["app_erate_ratio_C14_Be10"].to_numpy())
    finally:
        plt.close(fig)


def test_overlay_reads_the_same_ratio_from_both_tables():
    """Crossing the nuclide pairs between scatter and curve would be silent."""
    from lsdviztools3.charts import plot_ratio_overlay_1d
    import matplotlib.pyplot as plt

    df = _joined_samples_frame()
    sweep = _sweep_frame()
    fig, ax = plot_ratio_overlay_1d(
        df, [("matched", sweep)], ratio_col="app_erate_ratio_C14_Al26",
    )
    try:
        assert np.allclose(ax.collections[0].get_offsets()[:, 1],
                           df["app_erate_ratio_C14_Al26"].to_numpy())
        assert np.allclose(ax.lines[0].get_ydata(),
                           sweep["app_erate_ratio_C14_Al26"].to_numpy())
    finally:
        plt.close(fig)


def test_overlay_keeps_the_t_zero_anchor_samples():
    """The t = 0 anchors are the experiment's strongest check; a log x axis
    would drop them silently, so symlog must draw every one."""
    from lsdviztools3.charts import plot_ratio_overlay_1d
    import matplotlib.pyplot as plt

    df = _joined_samples_frame()
    assert (df["sc_steptime_yr"] == 0.0).sum() >= 1
    fig, ax = plot_ratio_overlay_1d(df, [("matched", _sweep_frame())])
    try:
        assert ax.get_xscale() == "symlog"
        assert ax.collections[0].get_offsets().shape[0] == len(df)
    finally:
        plt.close(fig)


def test_overlay_refuses_a_log_x_axis():
    from lsdviztools3.charts import plot_ratio_overlay_1d

    with pytest.raises(ValueError, match="anchor"):
        plot_ratio_overlay_1d(
            _joined_samples_frame(), [("matched", _sweep_frame())],
            xscale="log",
        )


def test_overlay_sweeps_take_distinct_styles():
    from lsdviztools3.charts import plot_ratio_overlay_1d
    import matplotlib.pyplot as plt

    sweeps = [("a", _sweep_frame()), ("b", _sweep_frame()),
              ("c", _sweep_frame())]
    fig, ax = plot_ratio_overlay_1d(_joined_samples_frame(), sweeps)
    try:
        assert len(ax.lines) == 3
        styles = {(ln.get_color(), ln.get_linestyle()) for ln in ax.lines}
        assert len(styles) == 3
    finally:
        plt.close(fig)


def test_overlay_reports_a_sweep_missing_its_columns():
    from lsdviztools3.charts import plot_ratio_overlay_1d

    bad = _sweep_frame().drop(columns=["time_since_step_yr"])
    with pytest.raises(KeyError, match="time_since_step_yr"):
        plot_ratio_overlay_1d(_joined_samples_frame(), [("bad", bad)])


def test_overlay_legend_headroom_clears_the_curve_peak():
    from lsdviztools3.charts import plot_ratio_overlay_1d
    import matplotlib.pyplot as plt

    sweep = _sweep_frame()
    fig, ax = plot_ratio_overlay_1d(_joined_samples_frame(),
                                    [("matched", sweep)])
    try:
        # Headroom is added above the data, so nothing drawn is clipped.
        assert ax.get_ylim()[1] > sweep["app_erate_ratio_C14_Be10"].max()
    finally:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Valley spacing
#
# The charts exist to make two failure modes VISIBLE, so the tests check that
# they stay visible rather than merely that a figure is produced.
# ---------------------------------------------------------------------------

def _spacing_frame():
    import pandas as pd
    return pd.DataFrame({
        "stem": ["a", "b", "c"],
        "channelization_index": [31.0, 78.0, 780.0],
        "convex_wavelength": [412.0, 245.0, 58.0],
        "spectral_wavelength": [500.0, 262.0, 287.0],
        "spectral_peak_wavelength": [500.0, 250.0, 250.0],
        "alpha": [-2.3, -3.7, -3.4],
        "alpha_predicted": [-2.1, -2.1, -2.1],
        "valleys_in_domain": [2.4, 4.1, 17.2],
        "under_sampled": [True, True, False],
    })


def test_spacing_vs_ci_shows_the_critical_band():
    """The Bonetti band must be INSIDE the axes, not clipped to the spine.

    The lowest run sits at CI = 31 and the band spans 32-37, so without an
    explicit x-limit the shading lands against the left edge and the legend
    carries an entry with no visible swatch -- a legend that promises a
    comparison the reader cannot actually make.
    """
    from lsdviztools3.charts import plot_spacing_vs_ci
    from lsdviztools3.io.valley_spacing import CI_CRITICAL_HIGH, CI_CRITICAL_LOW

    fig = plot_spacing_vs_ci(_spacing_frame())
    ax = fig.axes[0]
    lo, hi = ax.get_xlim()
    assert lo < CI_CRITICAL_LOW, "critical band clipped at the left spine"
    assert hi > CI_CRITICAL_HIGH
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any("Bonetti" in t for t in labels)


def test_spacing_vs_ci_marks_under_sampled_points_separately():
    """Under-sampled spectral points are drawn, but not joined by a trend line."""
    from lsdviztools3.charts import plot_spacing_vs_ci

    fig = plot_spacing_vs_ci(_spacing_frame())
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert any("too few valleys" in t for t in labels)


def test_spacing_vs_ci_optional_bonetti_reference():
    """The published lambda_c is opt-in: absent unless the caller asks for it."""
    from lsdviztools3.charts import plot_spacing_vs_ci

    without = plot_spacing_vs_ci(_spacing_frame())
    assert not any("lambda_c" in t.get_text() or "\\lambda_c" in t.get_text()
                   for t in without.axes[0].get_legend().get_texts())

    with_ref = plot_spacing_vs_ci(_spacing_frame(), bonetti_lambda_c=42.0)
    assert len(with_ref.axes[0].get_legend().get_texts()) > \
           len(without.axes[0].get_legend().get_texts())


def test_landscape_panel_one_axis_per_raster(tmp_path):
    """A sniff-test panel that silently drops a landscape is worse than none."""
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin
    from lsdviztools3.charts import plot_landscape_panel

    paths = []
    for i in range(3):
        p = tmp_path / f"dem{i}.tif"
        z = np.tile(np.sin(np.linspace(0, 6 * np.pi, 40)) * (i + 1), (40, 1))
        with rasterio.open(p, "w", driver="GTiff", height=40, width=40, count=1,
                           dtype="float32",
                           transform=from_origin(0, 200, 5, 5)) as dst:
            dst.write(z.astype("float32"), 1)
        paths.append(p)

    fig = plot_landscape_panel(paths, labels=["a", "b", "c"], ncols=3,
                               spacings=[100.0, 100.0, 100.0])
    drawn = [ax for ax in fig.axes if ax.images]
    assert len(drawn) == 3
