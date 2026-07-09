"""Valley cross-section plotting (matplotlib).

:func:`plot_valley_cross_sections` draws a grid of regularly-spaced transects
orthogonal to a valley centreline, from the per-sample table
(``*_valley_cross_section_samples.csv``) written by ``lsdtt-valley-centreline``
with ``write_valley_cross_sections``. Each panel plots elevation against the
signed offset from the centreline, shades the in-valley run, and marks the
detected left / right valley sides. Every transect of a channel shares the same
length (a fixed multiple of the widest valley on that channel), so the panels
are directly comparable.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

_MPL_MSG = (
    "Cross-section plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def plot_valley_cross_sections(
    df: pd.DataFrame,
    *,
    source_id: int = 0,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    offset_col: str = "offset_m",
    elevation_col: str = "elevation",
    in_valley_col: str = "in_valley",
    along_col: str = "along_distance_m",
    left_col: str = "left_offset_m",
    right_col: str = "right_offset_m",
    max_panels: int = 30,
    n_cols: int = 5,
    elevation_color: str = "0.35",
    valley_color: str = "tab:blue",
    left_color: str = "tab:green",
    right_color: str = "tab:red",
):
    """Grid of valley cross sections orthogonal to the centreline.

    Args:
        df: Per-sample table (``*_valley_cross_section_samples.csv``) with the
            source, station, offset, elevation, in-valley, along-distance and
            left/right-offset columns. Non-finite elevations are dropped from
            the elevation line but kept for the in-valley shading.
        source_id: Channel source to plot.
        source_col, station_col, offset_col, elevation_col, in_valley_col,
        along_col, left_col, right_col: Column names in ``df``.
        max_panels: Cap on the number of transect panels; when a channel has
            more stations they are evenly subsampled (the vector/CSV outputs
            keep every station — this is a visual summary). 0 or negative =
            draw all stations.
        n_cols: Panels per row.
        elevation_color: Colour of the elevation profile line.
        valley_color: Fill colour for the in-valley run.
        left_color, right_color: Colours of the left / right valley-side marks.

    Returns:
        ``(fig, axes)`` — axes as a flat list, one per transect drawn.

    Raises:
        ValueError: if the table holds no stations for ``source_id``.
    """
    plt = _import_mpl()
    from matplotlib.lines import Line2D

    ch = df[df[source_col] == source_id]
    stations = sorted(ch[station_col].unique())
    if not stations:
        raise ValueError(
            f"No cross-section stations for {source_col}={source_id}.")
    if max_panels and max_panels > 0 and len(stations) > max_panels:
        idx = np.linspace(0, len(stations) - 1, max_panels).round().astype(int)
        stations = [stations[i] for i in sorted(set(idx))]

    n = len(stations)
    n_cols = max(1, min(n_cols, n))
    n_rows = (n + n_cols - 1) // n_cols
    fig, axarr = plt.subplots(
        n_rows, n_cols, figsize=(3.1 * n_cols, 2.3 * n_rows), squeeze=False)

    ev = ch[elevation_col].replace([np.inf, -np.inf], np.nan).dropna()
    ymin, ymax = (float(ev.min()), float(ev.max())) if len(ev) else (0.0, 1.0)
    pad = 0.05 * (ymax - ymin + 1e-6)
    xmax = float(ch[offset_col].abs().max())

    axes = []
    for k, sseq in enumerate(stations):
        ax = axarr[k // n_cols][k % n_cols]
        axes.append(ax)
        s = ch[ch[station_col] == sseq].sort_values(offset_col)
        off = s[offset_col].to_numpy(dtype=float)
        elev = s[elevation_col].to_numpy(dtype=float)
        inval = s[in_valley_col].to_numpy(dtype=float)
        lo = float(s[left_col].iloc[0])
        ro = float(s[right_col].iloc[0])
        along = float(s[along_col].iloc[0])

        ax.plot(off, elev, "-", color=elevation_color, lw=1.0)
        ax.fill_between(off, ymin - pad, elev, where=(inval > 0),
                        color=valley_color, alpha=0.25, step="mid")
        if np.isfinite(lo):
            ax.axvline(lo, color=left_color, lw=1.1, ls="--")
        if np.isfinite(ro):
            ax.axvline(ro, color=right_color, lw=1.1, ls="--")
        ax.axvline(0.0, color="k", lw=0.7, ls=":")   # centreline
        ax.set_ylim(ymin - pad, ymax + pad)
        ax.set_xlim(-xmax, xmax)
        ax.set_title(f"{along / 1000:.2f} km  (w={ro - lo:.0f} m)", fontsize=8)
        ax.tick_params(labelsize=7)
        if k % n_cols == 0:
            ax.set_ylabel("elev (m)", fontsize=8)
        if k // n_cols == n_rows - 1:
            ax.set_xlabel("offset (m)", fontsize=8)

    for k in range(n, n_rows * n_cols):
        axarr[k // n_cols][k % n_cols].set_axis_off()

    handles = [
        Line2D([0], [0], color=elevation_color, lw=1.0, label="elevation"),
        Line2D([0], [0], color=left_color, lw=1.1, ls="--", label="left side"),
        Line2D([0], [0], color=right_color, lw=1.1, ls="--", label="right side"),
        Line2D([0], [0], color="k", lw=0.7, ls=":", label="centreline"),
    ]
    fig.legend(handles=handles, loc="upper right", fontsize=8, frameon=False,
               ncol=4, bbox_to_anchor=(0.995, 0.998))
    fig.suptitle(
        f"channel {source_id}: valley cross sections (orthogonal to centreline)",
        fontsize=11, y=0.999)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig, axes
