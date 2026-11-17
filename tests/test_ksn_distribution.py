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


def test_reference_ksn_rounds_to_one_significant_figure():
    from lsdviztools3.charts.chi import reference_ksn, round_ksn
    assert [round_ksn(x) for x in (13, 87, 217, 2.2, 0.47)] == [10, 90, 200, 2, 0.5]
    assert reference_ksn(np.arange(1, 101)) == [10, 50, 90]


def test_gradient_key_is_an_extra_profile_at_plot_scale():
    from matplotlib.collections import LineCollection
    fig, ax = plot_chi_stacked(_many_basins(3), gradient_key=[5, 20])
    key = [l for l in ax.lines if l.get_color() == "#3d3d3a"]
    assert len(key) == 2
    slopes = sorted(np.diff(l.get_ydata())[0] / np.diff(l.get_xdata())[0] for l in key)
    assert np.allclose(slopes, [5, 20])              # slope in data units = k_sn
    data_xmax = max(np.asarray(seg)[:, 0].max() for c in ax.collections
                    if isinstance(c, LineCollection) for seg in c.get_segments())
    assert min(l.get_xdata()[0] for l in key) > data_xmax   # after the last basin
    xl = ax.get_xlim()
    assert max(l.get_xdata()[1] for l in key) < xl[1]       # inside the axes
    plt.close(fig)


def test_gradient_key_on_last_row_and_profile():
    fig, axes = plot_chi_stacked(_many_basins(), basins_per_row=5, gradient_key=True)
    assert not [l for l in axes[0].lines if l.get_color() == "#3d3d3a"]
    assert [l for l in axes[-1].lines if l.get_color() == "#3d3d3a"]
    assert axes[0].get_xlim() == axes[-1].get_xlim()
    plt.close(fig)
    fig, ax = plot_chi_profile(_segments(), basin=0, gradient_key=[50])
    (key,) = [l for l in ax.lines if l.get_color() == "#3d3d3a"]
    assert np.isclose(np.diff(key.get_ydata())[0] / np.diff(key.get_xdata())[0], 50)
    plt.close(fig)
