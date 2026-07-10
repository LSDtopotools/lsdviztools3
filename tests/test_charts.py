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
                         "hanc_m": b, "n_pixels": int(c)})
    return pd.DataFrame(rows)


def test_plot_hanc_window_populations_panels_and_breaks():
    import pandas as pd
    from lsdviztools3.charts import plot_hanc_window_populations

    breaks = pd.DataFrame({
        "source_id": [0, 0, 0],
        "window_centre_m": [1000.0, 3000.0, 5000.0],
        "break_m": [5.0, 7.5, 6.0],
        "own_break": [True, False, True],
    })
    fig, axes = plot_hanc_window_populations(
        _window_hist_frame(), breaks_df=breaks, source_id=0)
    assert len(axes) == 3
    # Each panel carries a break line legend; inherited window is marked.
    legends = [ax.get_legend() for ax in axes]
    assert all(lg is not None for lg in legends)
    assert "(inherited)" in legends[1].get_texts()[0].get_text()


def test_plot_hanc_window_populations_missing_source_raises():
    from lsdviztools3.charts import plot_hanc_window_populations

    with pytest.raises(ValueError, match="No HANC window histograms"):
        plot_hanc_window_populations(_window_hist_frame(), source_id=9)


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
                "hanc_m": b,
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

    dip60 = float(np.logspace(-1.3, 2, 128)[60])   # the dip HANC at bin 60
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
    # the decision region (not the full HANC range).
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
