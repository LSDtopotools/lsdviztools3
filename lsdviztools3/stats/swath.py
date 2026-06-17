"""Swath-profile plotting (matplotlib).

:func:`plot_swath` draws the classic swath envelope from an
``lsdtt-swath-profile`` per-bin table (see :func:`lsdviztools3.io.swath.load_swath`):
value (elevation by default) against distance along the baseline, with a
min/max envelope, an inter-quartile band, and median + mean lines.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

_MPL_MSG = (
    "Swath plotting needs matplotlib. Install the optional 'stats' extra, "
    "e.g. `pip install lsdviztools3[stats]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def plot_swath(
    df: pd.DataFrame,
    *,
    bands: Sequence[Tuple[int, int]] = ((0, 100), (25, 75)),
    show_median: bool = True,
    show_mean: bool = True,
    distance_units: str = "km",
    value_label: str = "Elevation (m)",
    color: str = "tab:blue",
    ax=None,
):
    """Plot a swath envelope: value vs distance along the baseline.

    Args:
        df: Swath per-bin table with ``along_axis_centre`` plus percentile
            columns ``p0``..``p100`` (and ``mean``). Empty bins (NaN) are
            dropped.
        bands: Percentile pairs to shade, low-to-high opacity by order. Default
            is the full min/max envelope ``(0, 100)`` plus the IQR ``(25, 75)``.
        show_median: Draw the p50 line.
        show_mean: Draw the mean line (dashed).
        distance_units: ``"km"`` (divide along-axis metres by 1000) or ``"m"``.
        value_label: Y-axis label (the swathed quantity).
        color: Base colour for bands and lines.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()

    data = df.dropna(subset=["p50"]) if "p50" in df.columns else df.dropna()
    if data.empty:
        raise ValueError("Swath table has no populated bins to plot.")

    scale = 1000.0 if distance_units == "km" else 1.0
    x = data["along_axis_centre"].to_numpy(dtype=float) / scale

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))
    else:
        fig = ax.figure

    # Shaded percentile bands: outermost (widest) drawn faintest.
    n = max(len(bands), 1)
    for i, (lo, hi) in enumerate(bands):
        lo_col, hi_col = f"p{lo}", f"p{hi}"
        if lo_col not in data.columns or hi_col not in data.columns:
            continue
        alpha = 0.15 + 0.25 * (i / max(n - 1, 1)) if n > 1 else 0.2
        ax.fill_between(
            x, data[lo_col], data[hi_col],
            color=color, alpha=alpha, linewidth=0,
            label=f"p{lo}–p{hi}",
        )

    if show_median and "p50" in data.columns:
        ax.plot(x, data["p50"], color=color, linewidth=1.2, label="median")
    if show_mean and "mean" in data.columns:
        ax.plot(x, data["mean"], color="0.2", linewidth=1.0,
                linestyle="--", label="mean")

    ax.set_xlabel(f"Distance along swath ({distance_units})")
    ax.set_ylabel(value_label)
    ax.margins(x=0)
    ax.legend(fontsize=8, loc="best", frameon=False)
    return fig, ax
