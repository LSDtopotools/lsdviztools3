"""k_sn distributions (per channel pixel / per segment) and classes from them."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import LineCollection

from lsdviztools3.charts.chi import plot_chi_profile, plot_chi_stacked
from lsdviztools3.charts.ksn_distribution import (ksn_colour_classes, ksn_values, mpl_colour,
                                                  plot_ksn_distributions)


def _segments():
    """Basin 0: two long gentle segments (k_sn 10 and 20, 45 nodes each) and
    nine short steep ones (k_sn 100..900, 2 nodes each) on branch 1."""
    rows = [dict(basin_index=0, branch_id=0, segment_id=j, k_sn=10.0 * (j + 1), n_nodes=45,
                 downstream_chi=4.5 * j, upstream_chi=4.5 * (j + 1), downstream_elevation=45.0 * j,
                 upstream_elevation=45.0 * (j + 1)) for j in range(2)]
    for i in range(9):
        rows.append(dict(basin_index=0, branch_id=1, segment_id=i + 2, k_sn=100.0 * (i + 1),
                         n_nodes=2, downstream_chi=5.0 + i * 0.2, upstream_chi=5.2 + i * 0.2,
                         downstream_elevation=50.0 + i * 20, upstream_elevation=70.0 + i * 20))
    return pd.DataFrame(rows)


def _nodes(seg):
    rows = []
    for _, s in seg.iterrows():
        for j in range(int(s.n_nodes)):
            rows.append(dict(basin_index=s.basin_index, branch_id=s.branch_id, k_sn=s.k_sn,
                             chi=s.downstream_chi + j * 0.01, observed_elevation=s.downstream_elevation))
    return pd.DataFrame(rows)


def test_two_bases_sample_differently():
    seg = _segments()
    nodes = _nodes(seg)
    assert ksn_values(seg).size == 11 and ksn_values(nodes).size == 108
    q_seg = ksn_colour_classes("quantile", segments=seg, basis="segments", n_classes=2)
    q_pix = ksn_colour_classes("quantile", nodes=nodes, basis="pixels", n_classes=2)
    # per pixel the long gentle segment dominates; per segment the steep ones do
    assert q_pix.breaks[0] < q_seg.breaks[0]
    assert q_pix.scheme == "quantile (pixels)" and q_seg.scheme == "quantile (segments)"
    removed = pd.DataFrame({"basin_index": [0], "branch_id": [1]})
    assert ksn_values(seg, removed).tolist() == [10.0, 20.0]
    with pytest.raises(ValueError, match="nodes"):
        ksn_colour_classes("quantile", segments=seg, basis="pixels")
    with pytest.raises(ValueError, match="basis"):
        ksn_colour_classes("quantile", segments=seg, basis="length")


def test_distribution_plot_shows_classes():
    seg = _segments()
    nodes = _nodes(seg)
    cc = ksn_colour_classes([50.0, 300.0], segments=seg, basis="segments")
    fig, (ap, aseg) = plot_ksn_distributions(seg, nodes, classes=cc)
    for ax in (ap, aseg):
        # two class breaks + the median
        assert len(ax.get_lines()) == 3
        faces = {tuple(np.round(p.get_facecolor()[:3], 3)) for p in ax.patches if p.get_height() > 0}
        assert faces <= {tuple(np.round(mpl_colour(c), 3)) for c in cc.colors}
    assert [t.get_text() for t in fig.legends[0].get_texts()] == list(cc.labels)
    plt.close(fig)
    fig, (ap, _) = plot_ksn_distributions(seg, nodes, classes=cc, class_basis="pixels")
    assert "classes found here" in ap.get_title(loc="left")
    plt.close(fig)
    fig, (ap, _) = plot_ksn_distributions(seg, nodes)
    assert len(ap.get_lines()) == 1 and not fig.legends
    plt.close(fig)


def test_profiles_coloured_by_class():
    seg = _segments()
    cc = ksn_colour_classes([50.0, 300.0], segments=seg, basis="segments")
    for fig, ax in (plot_chi_profile(seg, basin=0, classes=cc), plot_chi_stacked(seg, classes=cc)):
        lc = [c for c in ax.collections if isinstance(c, LineCollection)][-1]
        assert sorted(set(lc.get_array().tolist())) == [0.0, 1.0, 2.0]
        assert fig.axes[-1].get_yticklabels()[0].get_text() == cc.labels[0]
        plt.close(fig)
    fig, ax = plot_chi_profile(_nodes(seg), basin=0, classes=cc)
    assert set(ax.collections[-1].get_array().tolist()) == {0.0, 1.0, 2.0}
    plt.close(fig)


def _many_basins(n=9):
    rows = []
    for b in range(n):
        for j in range(3):
            rows.append(dict(basin_index=b, branch_id=0, segment_id=j, k_sn=10.0 * (j + 1),
                             n_nodes=5, downstream_chi=float(j), upstream_chi=float(j + 1),
                             downstream_elevation=100.0 * b + 10 * j, upstream_elevation=100.0 * b + 10 * (j + 1)))
    return pd.DataFrame(rows)


def test_default_basis_is_pixels():
    seg = _segments()
    cc = ksn_colour_classes("quantile", nodes=_nodes(seg), n_classes=2)
    assert cc.scheme.endswith("(pixels)")


def test_stack_rows_wrap_and_order():
    from lsdviztools3.charts.chi import stack_rows
    rows = stack_rows(_many_basins(), basins_per_row=5)
    assert [len(r) for r in rows] == [5, 4]
    rows = stack_rows(_many_basins(), basins=[8, 2, 5], basins_per_row=2)
    assert [list(r) for r in rows] == [[8, 2], [5]]
    # same gap in every row: each row starts at 0 with the same spacing
    assert rows[0][8] == 0 and rows[1][5] == 0
    assert len(stack_rows(_many_basins())) == 1
    with pytest.raises(ValueError):
        stack_rows(_many_basins(), basins_per_row=0)


def test_stacked_rows_share_scales():
    fig, axes = plot_chi_stacked(_many_basins(), basins_per_row=5)
    assert len(axes) == 2
    assert axes[0].get_ylim() == axes[1].get_ylim() and axes[0].get_xlim() == axes[1].get_xlim()
    lo, hi = axes[1].get_ylim()
    assert lo < 0 and hi > 830          # the union of both rows' elevations
    plt.close(fig)
    with pytest.raises(ValueError, match="ax=None"):
        _, ax = plt.subplots()
        plot_chi_stacked(_many_basins(), basins_per_row=5, ax=ax)
    plt.close("all")


def test_gradient_inset_matches_plot_scale():
    from lsdviztools3.charts.chi import add_gradient_inset, reference_ksn
    assert reference_ksn(np.arange(1, 101)) == [11.0, 50.0, 90.0]
    fig, ax = plot_chi_stacked(_many_basins(3))
    ins = add_gradient_inset(ax, [5, 20], loc="lower right")
    # data units per axes-fraction are the same in the inset and the plot
    bx, bi = ax.get_position(), ins.get_position()
    (xa, xb), (ya, yb) = ax.get_xlim(), ax.get_ylim()
    (ia, ib), (ja, jb) = ins.get_xlim(), ins.get_ylim()
    assert np.isclose((ib - ia) / bi.width, (xb - xa) / bx.width)
    assert np.isclose((jb - ja) / bi.height, (yb - ya) / bx.height)
    # each line's slope is its k_sn
    slopes = sorted(np.diff(l.get_ydata())[0] / np.diff(l.get_xdata())[0] for l in ins.lines)
    assert np.allclose(slopes, [5, 20])
    plt.close(fig)


def test_gradient_inset_auto_avoids_data():
    from lsdviztools3.charts.chi import best_inset_spot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(np.linspace(0, 1, 200), np.linspace(0, 1, 200))   # fills lower left and upper right
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    _, loc = best_inset_spot([ax])
    assert loc in ("upper left", "lower right")
    fig2, axes = plot_chi_stacked(_many_basins(), basins_per_row=5, gradient_inset=True)
    assert sum(len(a.child_axes) for a in axes) == 1
    plt.close("all")


def test_class_limits_and_end_classes():
    from lsdviztools3.colour_classes import make_colour_classes
    rng = np.random.default_rng(0)
    inside = 10 ** rng.uniform(np.log10(5), np.log10(250), 10000)
    low = np.full(50, 1.0)       # 0.5%: an end class
    high = np.full(5, 900.0)     # 0.05%: joins the top class
    v = np.concatenate([inside, low, high])
    cc = make_colour_classes(v, "quantile", n_classes=4, log_scale=True, class_min=5, class_max=250)
    assert cc.breaks[0] == 5.0 and 250.0 not in cc.breaks
    assert cc.n_classes == 5                     # 4 between the limits + "below 5"
    assert all(5 < b < 250 for b in cc.breaks[1:])
    assert cc.labels[0] == "below 5"
    # both ends present
    v2 = np.concatenate([v, np.full(100, 900.0)])
    cc2 = make_colour_classes(v2, "quantile", n_classes=4, log_scale=True, class_min=5, class_max=250)
    assert cc2.breaks[0] == 5.0 and cc2.breaks[-1] == 250.0 and cc2.n_classes == 6
    # neither end: no data beyond the limits
    cc3 = make_colour_classes(inside, "equal_interval", n_classes=3, log_scale=True,
                              class_min=5, class_max=250)
    assert cc3.n_classes == 3
    assert np.allclose(np.log10(cc3.breaks), np.log10(5) + np.array([1, 2]) * np.log10(50) / 3)
    # breaks given by hand are untouched
    cc4 = make_colour_classes(v, [20, 40], log_scale=True, class_min=5, class_max=250)
    assert cc4.breaks == (20.0, 40.0)


def test_ksn_classes_default_limits():
    seg = _segments()
    nodes = _nodes(seg)
    cc = ksn_colour_classes("quantile", nodes=nodes, n_classes=2)
    # 10, 20 are inside [5, 250); 100..900 has 300..900 above 250 (12 of 108 pixels)
    assert cc.breaks[-1] == 250.0 and 5.0 not in cc.breaks
    cc = ksn_colour_classes("quantile", nodes=nodes, n_classes=2, class_min=None, class_max=None)
    assert 250.0 not in cc.breaks
