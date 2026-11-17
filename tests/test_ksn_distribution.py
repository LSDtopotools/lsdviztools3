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
    cc = ksn_colour_classes([50.0, 300.0], segments=seg)
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
    cc = ksn_colour_classes([50.0, 300.0], segments=seg)
    for fig, ax in (plot_chi_profile(seg, basin=0, classes=cc), plot_chi_stacked(seg, classes=cc)):
        lc = [c for c in ax.collections if isinstance(c, LineCollection)][-1]
        assert sorted(set(lc.get_array().tolist())) == [0.0, 1.0, 2.0]
        assert fig.axes[-1].get_yticklabels()[0].get_text() == cc.labels[0]
        plt.close(fig)
    fig, ax = plot_chi_profile(_nodes(seg), basin=0, classes=cc)
    assert set(ax.collections[-1].get_array().tolist()) == {0.0, 1.0, 2.0}
    plt.close(fig)
