"""Height-above-channel (HANC) profile plotting (matplotlib).

:func:`plot_hanc_distance` draws HANC against distance along the channel for
classified floodplain / terrace pixels, from an ``lsdtt-flat-surface-finder``
per-pixel table (``*_hanc_profile.csv``, written by ``write_hanc_profile_csv``).
Optionally overlays the reach-relative band break levels from the companion
``*_hanc_band_breaks.csv`` table.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Mapping, Optional

import numpy as np
import pandas as pd

_MPL_MSG = (
    "HANC-profile plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Surface-class codes in the lsdtt3 ``_surface_class`` convention.
_DEFAULT_CLASS_LABELS = {1: "floodplain", 2: "terrace"}
_DEFAULT_CLASS_COLORS = {1: "tab:blue", 2: "tab:orange"}


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def plot_hanc_distance(
    df: pd.DataFrame,
    *,
    distance_col: str = "channel_flow_distance_m",
    hanc_col: str = "hanc_m",
    class_col: str = "class",
    source_col: Optional[str] = "source_id",
    source_id: Optional[int] = None,
    class_labels: Mapping[int, str] = _DEFAULT_CLASS_LABELS,
    class_colors: Mapping[int, str] = _DEFAULT_CLASS_COLORS,
    breaks: Optional[pd.DataFrame] = None,
    distance_units: str = "km",
    point_size: float = 2.0,
    alpha: float = 0.25,
    max_points: int = 200_000,
    seed: int = 0,
    log_y: bool = False,
    ax=None,
):
    """Plot HANC vs distance along the channel, coloured by surface class.

    Args:
        df: Per-pixel HANC profile table (``*_hanc_profile.csv``) with at
            least the distance, HANC, and class columns.
        distance_col: Along-channel flow distance column (metres, measured
            from the DEM outlet).
        hanc_col: Height-above-channel column (metres).
        class_col: Integer surface-class column (1 floodplain, 2 terrace).
        source_col: Channel source-id column, used only when ``source_id``
            filters to one channel. ``None`` pools all sources.
        source_id: Plot only pixels attributed to this channel source.
        class_labels: Legend text per class code.
        class_colors: Point colour per class code.
        breaks: Optional band-breaks table (``*_hanc_band_breaks.csv`` with
            ``window_centre_m`` and semicolon-separated ``breaks_m``) to
            overlay the reach-relative break levels.
        distance_units: ``"km"`` (divide metres by 1000) or ``"m"``.
        point_size: Scatter marker size.
        alpha: Scatter marker opacity.
        max_points: Random subsample cap per class (pixel tables can run to
            millions of rows). 0 or negative = no cap.
        seed: Subsample RNG seed (reproducible figures).
        log_y: Logarithmic HANC axis (drops non-positive values).
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``.

    Raises:
        ValueError: if no plottable pixels remain after filtering.
    """
    plt = _import_mpl()

    data = df
    if source_id is not None and source_col is not None:
        data = data[data[source_col] == source_id]

    scale = 1000.0 if distance_units == "km" else 1.0

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.figure

    rng = np.random.default_rng(seed)
    n_drawn = 0
    for code, label in class_labels.items():
        cls = data[data[class_col] == code]
        x = cls[distance_col].to_numpy(dtype=float)
        y = cls[hanc_col].to_numpy(dtype=float)
        keep = np.isfinite(x) & np.isfinite(y)
        if log_y:
            keep &= y > 0
        x, y = x[keep], y[keep]
        if x.size == 0:
            continue
        if 0 < max_points < x.size:
            idx = rng.choice(x.size, size=max_points, replace=False)
            x, y = x[idx], y[idx]
        ax.scatter(
            x / scale, y,
            s=point_size, alpha=alpha, linewidths=0,
            color=class_colors.get(code, "0.5"), label=label,
            rasterized=True,
        )
        n_drawn += x.size

    if n_drawn == 0:
        raise ValueError("No plottable pixels in the HANC profile table.")

    if breaks is not None and not breaks.empty:
        bx, by = [], []
        brk = breaks
        if source_id is not None and "source_id" in brk.columns:
            brk = brk[brk["source_id"] == source_id]
        for _, row in brk.iterrows():
            levels = str(row.get("breaks_m", "") or "")
            for lvl in levels.split(";"):
                if not lvl:
                    continue
                bx.append(float(row["window_centre_m"]) / scale)
                by.append(float(lvl))
        if bx:
            ax.scatter(bx, by, marker="_", s=120, color="0.15",
                       linewidths=1.2, label="band break", zorder=3)

    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel(f"Distance along channel ({distance_units})")
    ax.set_ylabel("Height above channel (m)")
    ax.margins(x=0)
    legend = ax.legend(fontsize=8, loc="best", frameon=False,
                       scatterpoints=1, markerscale=4)
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)
    return fig, ax


def plot_hanc_window_populations(
    hist_df: pd.DataFrame,
    *,
    breaks_df: Optional[pd.DataFrame] = None,
    source_id: int = 0,
    hanc_col: str = "hanc_m",
    count_col: str = "n_pixels",
    window_col: str = "window_centre_m",
    smooth_bins: int = 5,
    max_windows: int = 12,
    n_cols: int = 3,
):
    """Per-window HANC populations with the selected valley break.

    Draws one panel per along-channel window from the
    ``*_valley_hanc_window_hist.csv`` table written by
    ``lsdtt-valley-centreline``, so users can see the population each
    reach-adaptive threshold was selected from. When ``breaks_df`` (the
    companion ``*_valley_breaks.csv``) is given, the chosen break is drawn
    on each panel, marked inherited when the window did not produce its own.

    Args:
        hist_df: Per-window histogram table (``source_id``, window centre,
            bin-centre HANC, pixel count).
        breaks_df: Optional per-window break table (``source_id``,
            ``window_centre_m``, ``break_m``, ``own_break``).
        source_id: Channel to plot.
        hanc_col: Bin-centre column (m).
        count_col: Pixel-count column.
        window_col: Window-centre column (m along the channel).
        smooth_bins: Moving-average width (bins) for the drawn density.
        max_windows: Cap on the number of panels.
        n_cols: Panels per row.

    Returns:
        ``(fig, axes)`` — axes as a flat list, one per window drawn.

    Raises:
        ValueError: if the table holds no windows for ``source_id``.
    """
    plt = _import_mpl()

    data = hist_df[hist_df["source_id"] == source_id]
    windows = sorted(data[window_col].unique())[:max_windows]
    if not windows:
        raise ValueError(
            f"No HANC window histograms for source_id={source_id}.")

    n_rows = (len(windows) + n_cols - 1) // n_cols
    fig, axarr = plt.subplots(
        n_rows, n_cols, figsize=(4.6 * n_cols, 2.9 * n_rows),
        sharex=True, squeeze=False)
    axes = []
    for k, w in enumerate(windows):
        ax = axarr[k // n_cols][k % n_cols]
        axes.append(ax)
        g = data[data[window_col] == w].sort_values(hanc_col)
        counts = g[count_col].to_numpy(dtype=float)
        if smooth_bins > 1 and counts.size > smooth_bins:
            kernel = np.ones(smooth_bins) / smooth_bins
            pad = smooth_bins // 2
            padded = np.pad(counts, pad, mode="edge")
            counts = np.convolve(padded, kernel, mode="valid")[:len(g)]
        ax.plot(g[hanc_col].to_numpy(), counts, color="0.15", lw=1.0)
        ax.set_xscale("log")
        ax.set_title(f"window {w / 1000:.1f} km", fontsize=9)

        if breaks_df is not None:
            row = breaks_df[
                (breaks_df["source_id"] == source_id)
                & np.isclose(breaks_df[window_col], w)]
            if len(row):
                own = bool(row["own_break"].iloc[0])
                ax.axvline(
                    row["break_m"].iloc[0], color="crimson", ls="--", lw=1.2,
                    label=f"break {row['break_m'].iloc[0]:.1f} m"
                          + ("" if own else " (inherited)"))
                ax.legend(fontsize=7, frameon=False)

    for k in range(len(windows), n_rows * n_cols):
        axarr[k // n_cols][k % n_cols].set_axis_off()
    for ax in axes[-n_cols:]:
        ax.set_xlabel("HANC (m)")
    for r in range(n_rows):
        axarr[r][0].set_ylabel("pixel count")
    fig.suptitle(
        f"channel {source_id}: corridor HANC populations per reach window",
        fontsize=11)
    fig.tight_layout()
    return fig, axes


#: Colours for the break-metric figure.
_BREAK_METRIC_COLORS = {
    "valley": "#8ecae6",    # valley band shade (HANC <= break)
    "curve": "#1d3557",     # smoothed decision metric
    "raw": "0.75",          # raw histogram bars
    "dip": "#d1495b",       # raw detected dip on the metric
    "applied": "#5f0f40",   # applied break (post clamp / median / inheritance)
    "peak": "tab:orange",   # surviving mode peaks
    "boundary": "0.5",      # surviving mode boundaries
}


def plot_valley_break_metric(
    metric_df: pd.DataFrame,
    *,
    breaks_df: Optional[pd.DataFrame] = None,
    source_id: int = 0,
    n_cols: int = 4,
    colors: Mapping[str, str] = _BREAK_METRIC_COLORS,
):
    """Reviewer figure: how the valley HANC break is chosen.

    Shows the ACTUAL metric the break is decided on — the SMOOTHED log10-HANC
    density (a 2-pass moving average), NOT the raw, bumpy histogram — from the
    ``*_valley_break_metric.csv`` table written by ``lsdtt-valley-centreline``
    (``write_valley_break_metric``). The valley is everything at or below the
    first surviving inter-mode dip (the break). Each along-channel window has
    its own pooled (window +/- 1) decision curve, so the break adapts from a
    confined canyon to a broad valley floor.

    Per channel the figure draws a "whole corridor" panel (global smoothed
    metric with the raw histogram faint behind it, and mode peaks / surviving
    boundaries / the chosen break marked when the global population yields a
    break) followed by a grid of per-window panels. When ``breaks_df`` (the
    companion ``*_valley_breaks.csv``) is supplied, each window panel shows the
    RAW detected dip on the metric AND the APPLIED break (after the min-break
    clamp, median-of-three smoothing and inheritance), so the metric and the
    post-processing are both visible and the window is labelled own/inherited.

    Args:
        metric_df: Break-metric table with columns ``source_id``, ``scope``
            (``"global"`` / ``"window"``), ``window_centre_m``, ``hanc_m``,
            ``raw_count``, ``smoothed``, ``is_peak``, ``is_boundary``,
            ``is_break``.
        breaks_df: Optional per-window applied-break table (``source_id``,
            ``window_centre_m``, ``break_m``, ``own_break``).
        source_id: Channel to plot.
        n_cols: Panels per row.
        colors: Override any of the figure colours (see
            ``_BREAK_METRIC_COLORS`` for the keys).

    Returns:
        ``(fig, axes)`` — axes as a flat list (global panel first).

    Raises:
        ValueError: if the table holds no rows for ``source_id``.
    """
    plt = _import_mpl()
    from matplotlib.lines import Line2D

    ch = metric_df[metric_df["source_id"] == source_id]
    if ch.empty:
        raise ValueError(
            f"No break-metric rows for source_id={source_id}.")
    gl = ch[ch["scope"] == "global"]
    wins = ch[ch["scope"] == "window"]
    wcs = sorted(wins["window_centre_m"].unique())

    ab = None if breaks_df is None else breaks_df[
        breaks_df["source_id"] == source_id]

    def applied_at(wc):
        """(break_m, own) for the reach whose centre is nearest ``wc``."""
        if ab is None or len(ab) == 0:
            return None, None
        j = (ab["window_centre_m"] - wc).abs().idxmin()
        return float(ab.loc[j, "break_m"]), bool(ab.loc[j, "own_break"])

    def break_hanc(sub):
        b = sub[sub["is_break"] == 1]
        return float(b["hanc_m"].iloc[0]) if len(b) else None

    def draw_panel(ax, sub, title, show_raw, applied_break=None):
        sub = sub.sort_values("hanc_m")
        x = sub["hanc_m"].to_numpy(dtype=float)
        sm = sub["smoothed"].to_numpy(dtype=float)
        smn = sm / sm.max() if sm.max() > 0 else sm
        if show_raw and "raw_count" in sub and sub["raw_count"].max() > 0:
            raw = sub["raw_count"].to_numpy(dtype=float)
            ax.bar(x, raw / raw.max(), width=x * 0.08, color=colors["raw"],
                   align="center", zorder=1)
        ax.plot(x, smn, color=colors["curve"], lw=1.6, zorder=3)
        band = applied_break if applied_break is not None else break_hanc(sub)
        if band is not None and np.isfinite(band):
            ax.axvspan(x.min(), band, color=colors["valley"], alpha=0.35,
                       zorder=0)
        bk = break_hanc(sub)
        if bk is not None:
            ax.axvline(bk, color=colors["dip"], lw=1.6, ls="--", zorder=4)
        if applied_break is not None and np.isfinite(applied_break):
            ax.axvline(applied_break, color=colors["applied"], lw=1.8, zorder=5)
            ax.text(applied_break, 0.92, f" {applied_break:.1f} m",
                    color=colors["applied"], fontsize=7, ha="left", va="top",
                    zorder=6, transform=ax.get_xaxis_transform())
        for _, r in sub[sub["is_peak"] == 1].iterrows():
            ax.plot(r["hanc_m"], np.interp(r["hanc_m"], x, smn), "v", ms=5,
                    color=colors["peak"], zorder=6)
        for _, r in sub[sub["is_boundary"] == 1].iterrows():
            ax.axvline(r["hanc_m"], color=colors["boundary"], lw=0.7, ls=":",
                       zorder=2)
        ax.set_xscale("log")
        ax.set_xlim(max(x.min(), 0.05), x.max())
        ax.set_ylim(0, 1.2)
        ax.set_title(title, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_yticks([])

    npan = 1 + len(wcs)
    n_cols = max(1, min(n_cols, npan))
    n_rows = (npan + n_cols - 1) // n_cols
    fig, axarr = plt.subplots(
        n_rows, n_cols, figsize=(3.4 * n_cols, 2.4 * n_rows), squeeze=False)
    axf = [axarr[i // n_cols][i % n_cols] for i in range(n_rows * n_cols)]

    gbk = break_hanc(gl) if len(gl) else None
    gtitle = ("whole corridor"
              + (f": break {gbk:.1f} m" if gbk is not None
                 else ": no single break\n-> reach-adaptive used"))
    draw_panel(axf[0], gl, gtitle, show_raw=True)

    for i, wc in enumerate(wcs, start=1):
        sub = wins[wins["window_centre_m"] == wc]
        abrk, own = applied_at(wc)
        tag = "" if own is None else ("  (own)" if own else "  (inherited)")
        draw_panel(axf[i], sub, f"reach @ {wc / 1000:.1f} km{tag}",
                   show_raw=False, applied_break=abrk)

    for j in range(npan, n_rows * n_cols):
        axf[j].set_axis_off()

    handles = [
        Line2D([0], [0], color=colors["curve"], lw=1.6,
               label="smoothed decision metric"),
        Line2D([0], [0], color=colors["dip"], lw=1.6, ls="--",
               label="raw detected dip"),
        Line2D([0], [0], color=colors["applied"], lw=1.8,
               label="applied break (clamp/median/inherit)"),
        plt.Rectangle((0, 0), 1, 1, color=colors["valley"], alpha=0.5,
                      label="valley (HANC <= break)"),
        Line2D([0], [0], marker="v", color=colors["peak"], lw=0, ms=6,
               label="mode peak"),
        Line2D([0], [0], color=colors["boundary"], lw=0.7, ls=":",
               label="mode boundary"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(
        f"channel {source_id}: how the valley HANC break is chosen\n"
        f"(x = height above nearest channel, log scale; valley = at or below "
        f"the first inter-mode dip)", fontsize=11)
    # X labels on the bottom row of drawn panels.
    for ax in axf[max(0, npan - n_cols):npan]:
        ax.set_xlabel("HANC (m)", fontsize=8)
    fig.tight_layout(rect=[0, 0.03, 1, 0.94])
    return fig, axf[:npan]
