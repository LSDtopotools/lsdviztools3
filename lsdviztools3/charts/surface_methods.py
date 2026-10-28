"""Charts comparing the floodplain / terrace methods of one lsdtt3 run.

Read from lsdtt3 ``lsdtt-fluvial-surfaces`` outputs (see
``render.surface_methods`` for the file set):

  plot_method_levels   per method, every classified patch at its median HAND
                       (``<stem>_<method>_surface_patches.csv``), sized by area,
                       coloured by class; the height clusters and the
                       clustering p-value from ``_patch_level_clusters.csv``.
  plot_mapped_scores   per method, floodplain and terrace precision, recall and
                       intersection-over-union against the mapped surfaces
                       (``<stem>_mapped_surface_scores.csv``).

Visualisation only: nothing is recomputed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..render.surface_methods import method_stem, read_method_list

CLASS_COLOURS = {"floodplain": "#1f6fd1", "terrace": "#e8890c"}


def plot_method_levels(stem, out_png: Optional[str] = None, *,
                       methods: Optional[Sequence[str]] = None,
                       log_y: bool = False, min_cluster_patches: int = 3,
                       figsize=None):
    """Patch median HAND per method: the levels each method found.

    Each classified patch is a dot at its median HAND, area as dot size, colour
    the class; the dots are spread sideways within a method's column so they do
    not hide one another (the spread means nothing). Cluster centres from
    ``<method stem>_patch_level_clusters.csv`` (scope ``all``) holding at least
    ``min_cluster_patches`` patches are short black bars, with the clustering
    p-value under the method name.

    Returns:
        The matplotlib ``Figure``.
    """
    all_methods = read_method_list(stem)
    methods = list(methods) if methods is not None else all_methods
    fig, ax = plt.subplots(figsize=figsize or (2.2 * len(methods) + 2.0, 5.0))
    rng = np.random.default_rng(0)   # sideways jitter only
    labels = []
    for x, m in enumerate(methods):
        ms = method_stem(stem, m, all_methods)
        label = m
        pf = Path(f"{ms}_surface_patches.csv")
        if pf.exists():
            p = pd.read_csv(pf)
            p = p[p["class"].isin(list(CLASS_COLOURS))]
            if len(p):
                a = p["area_m2"].to_numpy(dtype=float)
                size = 8 + 120 * np.sqrt(a / a.max())
                ax.scatter(x + rng.uniform(-0.3, 0.3, len(p)), p["median_field_m"],
                           s=size, c=[CLASS_COLOURS[c] for c in p["class"]],
                           alpha=0.6, edgecolors="none")
        cf = Path(f"{ms}_patch_level_clusters.csv")
        if cf.exists():
            c = pd.read_csv(cf)
            c = c[c["scope"] == "all"]
            if len(c):
                r = c.iloc[0]
                counts = str(r["cluster_n_patches"]).split()
                for v, nc in zip(str(r["cluster_centres_m"]).split(), counts):
                    if int(nc) >= min_cluster_patches:
                        ax.plot([x - 0.38, x + 0.38], [float(v)] * 2,
                                color="black", lw=1)
                label += f"\np = {float(r['p_value']):.3g}"
        labels.append(label)
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlim(-0.6, len(methods) - 0.4)
    if log_y:
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.1)
    ax.set_ylabel("patch median HAND (m)")
    handles = [plt.Line2D([], [], marker="o", ls="", color=v, label=k)
               for k, v in CLASS_COLOURS.items()]
    handles.append(plt.Line2D([], [], color="black", lw=1,
                              label=f"height cluster (>= {min_cluster_patches} patches)"))
    ax.legend(handles=handles, fontsize=8, loc="upper left")
    ax.set_title("Levels found by each method (dot size: patch area)", fontsize=10)
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=200)
    return fig


def plot_mapped_scores(stem, out_png: Optional[str] = None, *, figsize=(7.0, 3.6)):
    """Precision, recall and IoU of each method against the mapped surfaces.

    Reads ``<stem>_mapped_surface_scores.csv``; one panel per class, grouped
    bars per method. A NaN (no cells of the class from the method) is left out.
    """
    s = pd.read_csv(f"{stem}_mapped_surface_scores.csv")
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    metrics = ("precision", "recall", "iou")
    width = 0.8 / len(metrics)
    x = np.arange(len(s))
    for ax, cls in zip(axes, ("floodplain", "terrace")):
        for i, metric in enumerate(metrics):
            ax.bar(x + (i - 1) * width, s[f"{cls}_{metric}"], width, label=metric)
        ax.set_xticks(x)
        ax.set_xticklabels(s["method"], fontsize=8, rotation=20)
        ax.set_title(cls, fontsize=10)
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("score")
    axes[0].legend(fontsize=8)
    fig.suptitle(f"Against mapped surfaces ({int(s['n_scored_cells'].iloc[0])} cells)",
                 fontsize=10)
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=200)
    return fig
