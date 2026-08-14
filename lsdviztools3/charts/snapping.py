"""Charts for the synthetic-network snapping experiment (matplotlib).

The experiment varies OCN annealing temperature to produce a range of branching
geometry, then snaps each landscape to a self-consistent fluvial + hillslope
form. Two questions:

  (a) Does temperature alone produce a real range of branching geometry?
  (b) Does snapping REDUCE that range?

:func:`plot_branching_range` is the figure question (b) turns on. Junction-angle
*distributions* are plotted with the existing
:func:`lsdviztools3.charts.distributions.plot_distributions`; this module only
adds the experiment-specific summary views.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np

_MPL_MSG = (
    "Snapping charts need matplotlib. Install the optional 'charts' extra, "
    "e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def plot_metric_by_temperature(
    table,
    *,
    metric: str = "median_deg",
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    ax=None,
):
    """One line per snapping variant, across annealing temperature.

    Reads directly off the tidy table from
    :func:`lsdviztools3.io.snapping.collect_experiment`.

    The visual claim is the SLOPE and SPREAD of each line, not its absolute
    level: a variant whose line is flatter than ``no snapping`` has lost
    sensitivity to the annealing temperature, which is what question (b)
    predicts snapping should do.

    Args:
        table: Tidy table with ``tag``, ``variant_label`` and ``metric`` columns.
        metric: Column to plot.
        ylabel: Y-axis label; defaults to the metric name.
        title: Optional title.
        ax: Existing matplotlib ``Axes`` or ``None``.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    import pandas as pd

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 4.4))
    else:
        fig = ax.figure

    tags = list(dict.fromkeys(table["tag"]))
    x = np.arange(len(tags))

    for label, grp in table.groupby("variant_label", sort=False):
        grp = grp.set_index("tag").reindex(tags)
        y = pd.to_numeric(grp[metric], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(y).any():
            continue
        emphasis = label.startswith("no snapping")
        ax.plot(
            x,
            y,
            marker="o",
            linewidth=2.4 if emphasis else 1.6,
            linestyle="-" if emphasis else "--",
            color="black" if emphasis else None,
            zorder=3 if emphasis else 2,
            label=label,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(tags)
    ax.set_xlabel("OCN annealing temperature (fixed n_iter, fixed seed)")
    ax.set_ylabel(ylabel or metric)
    if title:
        ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    return fig, ax


def plot_branching_range(
    range_table,
    *,
    metric_label: str = "junction angle",
    title: Optional[str] = None,
    ax=None,
):
    """Range of a branching metric across temperatures, one bar per variant.

    This is the direct test of question (b). The ``no snapping`` bar is the
    range the OCN annealing produced; a snapped variant with a SHORTER bar has
    had geometric freedom removed by the requirement that the surface be
    self-consistent.

    Args:
        range_table: Output of :func:`lsdviztools3.io.snapping.branching_range`.
        metric_label: Human name of the metric, for the axis label.
        title: Optional title.
        ax: Existing matplotlib ``Axes`` or ``None``.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()

    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
    else:
        fig = ax.figure

    labels = list(range_table["variant_label"])
    values = np.asarray(range_table["range"], dtype=float)
    colors = ["0.25" if lab.startswith("no snapping") else "tab:blue"
              for lab in labels]

    bars = ax.bar(labels, values, color=colors)

    # Annotate each bar with its value relative to the un-snapped case, which is
    # the number the hypothesis is about.
    if "range_relative_to_pre" in range_table.columns:
        rel = np.asarray(range_table["range_relative_to_pre"], dtype=float)
        for bar, r in zip(bars, rel):
            if not np.isfinite(r):
                continue
            ax.annotate(
                f"{r:.2f}x",
                (bar.get_x() + bar.get_width() / 2.0, bar.get_height()),
                textcoords="offset points",
                xytext=(0, 3),
                ha="center",
                fontsize=9,
            )

    ax.set_ylabel(f"range of {metric_label} across temperatures")
    if title:
        ax.set_title(title)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    return fig, ax


def plot_flow_change_fraction(
    table,
    *,
    title: Optional[str] = None,
    ax=None,
):
    """Fraction of pixels whose D8 flow direction changed during snapping.

    Grouped bars, one group per temperature. This quantifies how much
    rearrangement each hillslope treatment demanded of the OCN -- a large value
    means the annealed network was far from anything the snapped physics would
    accept.

    Args:
        table: Tidy table from :func:`lsdviztools3.io.snapping.collect_experiment`
            (rows with ``variant == "pre"`` carry no value and are skipped).
        title: Optional title.
        ax: Existing matplotlib ``Axes`` or ``None``.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    import pandas as pd

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 4.0))
    else:
        fig = ax.figure

    sub = table[table["variant"] != "pre"]
    if "fraction_changed" not in sub.columns:
        raise KeyError(
            "table has no 'fraction_changed' column -- were the flow-direction "
            "rasters written for both the pre- and post-snap surfaces?"
        )

    tags = list(dict.fromkeys(sub["tag"]))
    variants = list(dict.fromkeys(sub["variant_label"]))
    width = 0.8 / max(len(variants), 1)
    x = np.arange(len(tags))

    for i, label in enumerate(variants):
        grp = sub[sub["variant_label"] == label].set_index("tag").reindex(tags)
        y = pd.to_numeric(grp["fraction_changed"], errors="coerce").to_numpy(float)
        ax.bar(x + i * width - 0.4 + width / 2.0, y, width=width, label=label)

    ax.set_xticks(x)
    ax.set_xticklabels(tags)
    ax.set_xlabel("OCN annealing temperature")
    ax.set_ylabel("fraction of pixels with changed D8 direction")
    if title:
        ax.set_title(title)
    ax.grid(alpha=0.3, axis="y")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    return fig, ax
