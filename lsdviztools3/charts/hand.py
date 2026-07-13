"""Height-above-channel (HAND) profile plotting (matplotlib).

:func:`plot_hand_distance` draws HAND against distance along the channel for
classified floodplain / terrace pixels, from an ``lsdtt-flat-surface-finder``
per-pixel table (``*_hand_profile.csv``, written by ``write_hand_profile_csv``).
Optionally overlays the reach-relative band break levels from the companion
``*_hand_band_breaks.csv`` table.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Mapping, Optional

import numpy as np
import pandas as pd

_MPL_MSG = (
    "HAND-profile plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Surface-class codes in the lsdtt3 ``_surface_class`` convention.
_DEFAULT_CLASS_LABELS = {1: "floodplain", 2: "terrace"}
_DEFAULT_CLASS_COLORS = {1: "tab:blue", 2: "tab:orange"}


_RASTERIO_MSG = (
    "The reach map needs rasterio to read the hillshade and reach-id rasters. "
    "Install it e.g. `pip install rasterio`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def _import_rasterio():
    try:
        import rasterio  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without rio
        raise ImportError(_RASTERIO_MSG) from exc
    return rasterio


def plot_hand_distance(
    df: pd.DataFrame,
    *,
    distance_col: str = "channel_flow_distance_m",
    hand_col: str = "hand_m",
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
    """Plot HAND vs distance along the channel, coloured by surface class.

    Args:
        df: Per-pixel HAND profile table (``*_hand_profile.csv``) with at
            least the distance, HAND, and class columns.
        distance_col: Along-channel flow distance column (metres, measured
            from the DEM outlet).
        hand_col: Height-above-channel column (metres).
        class_col: Integer surface-class column (1 floodplain, 2 terrace).
        source_col: Channel source-id column, used only when ``source_id``
            filters to one channel. ``None`` pools all sources.
        source_id: Plot only pixels attributed to this channel source.
        class_labels: Legend text per class code.
        class_colors: Point colour per class code.
        breaks: Optional band-breaks table (``*_hand_band_breaks.csv`` with
            ``window_centre_m`` and semicolon-separated ``breaks_m``) to
            overlay the reach-relative break levels.
        distance_units: ``"km"`` (divide metres by 1000) or ``"m"``.
        point_size: Scatter marker size.
        alpha: Scatter marker opacity.
        max_points: Random subsample cap per class (pixel tables can run to
            millions of rows). 0 or negative = no cap.
        seed: Subsample RNG seed (reproducible figures).
        log_y: Logarithmic HAND axis (drops non-positive values).
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
        y = cls[hand_col].to_numpy(dtype=float)
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
        raise ValueError("No plottable pixels in the HAND profile table.")

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


def plot_hand_window_populations(
    hist_df: pd.DataFrame,
    *,
    breaks_df: Optional[pd.DataFrame] = None,
    source_id: int = 0,
    hand_col: str = "hand_m",
    count_col: str = "n_pixels",
    window_col: str = "window_centre_m",
    smooth_bins: int = 5,
    max_windows: int = 12,
    n_cols: int = 3,
):
    """Per-window HAND populations with the selected valley break.

    Draws one panel per along-channel window from the
    ``*_valley_hand_window_hist.csv`` table written by
    ``lsdtt-valley-centreline``, so users can see the population each
    reach-adaptive threshold was selected from. When ``breaks_df`` (the
    companion ``*_valley_breaks.csv``) is given, the chosen break is drawn
    on each panel, marked inherited when the window did not produce its own.

    Args:
        hist_df: Per-window histogram table (``source_id``, window centre,
            bin-centre HAND, pixel count).
        breaks_df: Optional per-window break table (``source_id``,
            ``window_centre_m``, ``break_m``, ``own_break``).
        source_id: Channel to plot.
        hand_col: Bin-centre column (m).
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
            f"No HAND window histograms for source_id={source_id}.")

    n_rows = (len(windows) + n_cols - 1) // n_cols
    fig, axarr = plt.subplots(
        n_rows, n_cols, figsize=(4.6 * n_cols, 2.9 * n_rows),
        sharex=True, squeeze=False)
    axes = []
    for k, w in enumerate(windows):
        ax = axarr[k // n_cols][k % n_cols]
        axes.append(ax)
        g = data[data[window_col] == w].sort_values(hand_col)
        counts = g[count_col].to_numpy(dtype=float)
        if smooth_bins > 1 and counts.size > smooth_bins:
            kernel = np.ones(smooth_bins) / smooth_bins
            pad = smooth_bins // 2
            padded = np.pad(counts, pad, mode="edge")
            counts = np.convolve(padded, kernel, mode="valid")[:len(g)]
        ax.plot(g[hand_col].to_numpy(), counts, color="0.15", lw=1.0)
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
        ax.set_xlabel("HAND (m)")
    for r in range(n_rows):
        axarr[r][0].set_ylabel("pixel count")
    fig.suptitle(
        f"channel {source_id}: corridor HAND populations per reach window",
        fontsize=11)
    fig.tight_layout()
    return fig, axes


#: Colours for the break-metric figure.
_BREAK_METRIC_COLORS = {
    "valley": "#8ecae6",    # valley band shade (HAND <= break)
    "curve": "#1d3557",     # smoothed decision metric
    "raw": "0.75",          # raw histogram bars
    "dip": "#d1495b",       # raw detected dip on the metric
    "threshold": "#7b2cbf", # dip acceptance threshold (0.6x smaller flank)
    "applied": "#5f0f40",   # applied break (post clamp / median / inheritance)
    "peak": "tab:orange",   # surviving mode peaks
    "raw_peak": "0.55",     # candidate peaks before merging (merged-away)
    "boundary": "0.5",      # surviving mode boundaries
}


def plot_valley_break_metric(
    metric_df: pd.DataFrame,
    *,
    breaks_df: Optional[pd.DataFrame] = None,
    source_id: int = 0,
    n_cols: int = 4,
    valley_fraction: float = 0.6,
    colors: Mapping[str, str] = _BREAK_METRIC_COLORS,
):
    """Reviewer figure: how the valley HAND break is chosen.

    Shows the ACTUAL metric the break is decided on — the SMOOTHED HAND density
    (a 2-pass moving average of the log-binned height-above-nearest-channel
    histogram), NOT the raw, bumpy histogram — from the
    ``*_valley_break_metric.csv`` table written by ``lsdtt-valley-centreline``
    (``write_valley_break_metric``). Each along-channel window has its own
    pooled (window +/- 1) decision curve, so the break adapts from a confined
    canyon to a broad valley floor.

    How the break is defined (mirrors ``find_valley_break`` in the C++):

    * A **mode** is a bump in the smoothed density. The curve is split at its
      **local minima** — a bin ``b`` with ``sm[b] <= sm[b-1]`` and
      ``sm[b] < sm[b+1]`` — into segments; each segment is one mode, with a
      peak (its highest smoothed value) and a mass (its raw pixel count). The
      two lowest bins are excluded (the on-channel log-floor clamp spike).
    * The **inter-mode dip** is NOT a gradient. It is the smoothed density
      VALUE at the local minimum between two modes (the depth of the trough
      between the bumps).
    * A dip is **kept** (the two modes stay separate) only when it falls below
      ``valley_fraction`` x the SMALLER of its two flanking mode peaks;
      shallower dips are noise and the modes are merged. Modes whose mass is a
      negligible fraction of the population are also absorbed (except a
      protected low near-channel valley-floor mode).
    * The **break** is the first surviving dip — the minimum between the
      valley-floor mode and the first wall mode. The valley is every pixel at
      or below it.

    Design (Simon's critique): the x axis is **LINEAR** (0.1 m above the
    channel is obviously in the valley — log space hid that) and **cropped** to
    the decision region (the valley-floor mode, the dip, and the near flank of
    the first wall mode); pixels tens of metres above the channel, obviously
    not valley, are not shown. Each panel marks the mode peaks, the inter-mode
    minima, the chosen dip, and the dip acceptance threshold
    (``valley_fraction`` x smaller flanking peak, drawn in purple across its
    flanks), and — when ``breaks_df`` is supplied — the raw detected dip vs the
    APPLIED break (after the min-break clamp, median-of-three smoothing and
    inheritance). Each reach is tagged and captioned by how its break arose:
    ``(own)`` "accepted dip at X m" when the applied break equals the reach's
    own dip; ``(median)`` "dip at X m accepted median at Y m" when the applied
    break differs from the dip (replaced by the median-of-three across
    neighbours, or raised by the min-break clamp); ``(inherited)`` "inherited
    break at X m" when the reach produced no own dip.

    Args:
        metric_df: Break-metric table with columns ``source_id``, ``scope``
            (``"global"`` / ``"window"``), ``window_centre_m``, ``hand_m``,
            ``raw_count``, ``smoothed``, ``is_peak``, ``is_boundary``,
            ``is_break``.
        breaks_df: Optional per-window applied-break table (``source_id``,
            ``window_centre_m``, ``break_m``, ``own_break``).
        source_id: Channel to plot.
        n_cols: Panels per row.
        valley_fraction: The dip-acceptance fraction drawn on the figure; must
            match the run's ``break_valley_fraction`` (default 0.6, which is
            not CLI-exposed, so 0.6 is faithful for standard runs).
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

    def fraction_at(wc):
        """Acceptance fraction for the reach nearest ``wc``. Uses the
        per-reach ``break_fraction`` column (written by the relief-adaptive
        method) when present, else the scalar ``valley_fraction``."""
        if (ab is None or len(ab) == 0
                or "break_fraction" not in getattr(ab, "columns", [])):
            return valley_fraction
        j = (ab["window_centre_m"] - wc).abs().idxmin()
        v = ab.loc[j, "break_fraction"]
        try:
            v = float(v)
        except (TypeError, ValueError):
            return valley_fraction
        return v if np.isfinite(v) and v > 0 else valley_fraction

    def break_hand(sub):
        b = sub[sub["is_break"] == 1]
        return float(b["hand_m"].iloc[0]) if len(b) else None

    def peaks_hand(sub):
        """Surviving mode-peak HAND values, ascending."""
        return sorted(float(v) for v in sub[sub["is_peak"] == 1]["hand_m"])

    # --- Data-driven x crop (LINEAR axis) -----------------------------------
    # Simon: log space is misleading and pixels >~100 m above the channel are
    # obviously not valley. Crop to the decision region: the valley-floor mode,
    # the dip, and the near flank of the first wall mode. Upper limit = a
    # margin past the first WALL peak (the 2nd surviving peak) of the reach
    # with the highest such peak; fall back to a multiple of the applied break
    # where no reach segments (e.g. a confined all-wall gorge).
    wall_peaks = []
    for wc in wcs:
        p = peaks_hand(wins[wins["window_centre_m"] == wc])
        if len(p) >= 2:
            wall_peaks.append(p[1])
    applied_vals = []
    for wc in wcs:
        a, _ = applied_at(wc)
        if a is not None and np.isfinite(a):
            applied_vals.append(a)
    if wall_peaks:
        x_hi = 1.4 * max(wall_peaks)
    elif applied_vals:
        x_hi = 4.0 * max(applied_vals)
    else:
        x_hi = 40.0
    x_hi = float(min(max(x_hi, 15.0), 60.0))

    def draw_panel(ax, sub, title, show_raw, applied_break=None,
                   annotate_dip=True, subtitle=None, fraction=None):
        frac = fraction if fraction is not None else valley_fraction
        sub = sub.sort_values("hand_m")
        x = sub["hand_m"].to_numpy(dtype=float)
        sm = sub["smoothed"].to_numpy(dtype=float)
        in_crop = x <= x_hi
        pk = peaks_hand(sub)
        bk = break_hand(sub)

        # Y SCALING: the valley-floor mode and the dip that decide the break
        # are typically ~1-2% of the wall-mode height, so a wall-normalised
        # axis flattens them onto the baseline (exactly why the break looked
        # arbitrary). When the reach has a real break, scale y to its FLANKING
        # valley-side modes (the taller of the two peaks bracketing the dip)
        # so the dip and its threshold are legible; the wall mode then shoots
        # off the top and is clipped (a "wall mode ->" tag marks it). Otherwise
        # (no break / inherited / global) normalise to the max in the crop.
        def raw_at(hx):
            return float(np.interp(hx, x, sm))

        y_ref = None
        if bk is not None and np.isfinite(bk) and len(pk) >= 2:
            lpk = max([p for p in pk if p <= bk], default=pk[0])
            rpk = min([p for p in pk if p >= bk], default=pk[-1])
            # Scale to the SMALLER flank — the peak the dip threshold references
            # — so the dip and threshold are legible even when the other flank
            # is the towering wall mode (which then clips off the top).
            y_ref = min(raw_at(lpk), raw_at(rpk))
        if y_ref and y_ref > 0:
            norm = y_ref / 0.55         # smaller flanking mode sits at 0.55
        else:
            norm = sm[in_crop].max() if in_crop.any() and sm[in_crop].max() > 0 \
                else (sm.max() if sm.max() > 0 else 1.0)
        smn = sm / norm

        def y_at(hx):
            return float(np.interp(hx, x, smn))

        if show_raw and "raw_count" in sub and sub["raw_count"].max() > 0:
            raw = sub["raw_count"].to_numpy(dtype=float)
            rn = raw[in_crop].max() if in_crop.any() and raw[in_crop].max() > 0 \
                else raw.max()
            ax.fill_between(x, 0, raw / rn, color=colors["raw"], zorder=1,
                            step="mid", alpha=0.6)
        ax.plot(x, smn, color=colors["curve"], lw=1.8, zorder=3)
        # Tag the wall mode where the curve clips off the top.
        if y_ref:
            above = np.where(smn[in_crop] > 1.3)[0]
            if above.size:
                ax.annotate("wall mode", xy=(x[in_crop][above[0]], 1.3),
                            xytext=(0.62, 0.9), textcoords="axes fraction",
                            fontsize=6.5, color=colors["curve"], va="top",
                            arrowprops=dict(arrowstyle="->",
                                            color=colors["curve"], lw=0.8))

        band = applied_break if applied_break is not None else break_hand(sub)
        if band is not None and np.isfinite(band):
            ax.axvspan(0, band, color=colors["valley"], alpha=0.35, zorder=0)

        # Mode structure: surviving boundaries (the inter-mode minima) and the
        # mode peaks the segmentation kept (within the crop; the wall peak may
        # sit off the top of the zoomed axis).
        for hb in sub[sub["is_boundary"] == 1]["hand_m"]:
            if hb <= x_hi:
                ax.axvline(hb, color=colors["boundary"], lw=0.7, ls=":",
                           zorder=2)
        # Candidate peaks detected BEFORE merging: draw the merged-away ones as
        # faint hollow triangles so the merge decisions are visible (a bump you
        # can see but that is not a surviving mode was fused into a neighbour).
        if "is_raw_peak" in sub.columns:
            raw_pk = sorted(float(v)
                            for v in sub[sub["is_raw_peak"] == 1]["hand_m"])
            surviving = set(pk)
            for hp in raw_pk:
                if hp in surviving:
                    continue
                if hp <= x_hi and y_at(hp) <= 1.3:
                    ax.plot(hp, y_at(hp), "v", ms=6, mfc="none",
                            mec=colors["raw_peak"], mew=1.0, zorder=5)
        for hp in pk:
            if hp <= x_hi and y_at(hp) <= 1.3:
                ax.plot(hp, y_at(hp), "v", ms=6, color=colors["peak"], zorder=6)

        # Dip detection: the break is the first inter-mode minimum; it is
        # accepted because its density falls below valley_fraction x the
        # SMALLER of its two flanking mode peaks. Draw that threshold across
        # the flanks so a reviewer sees WHY the dip qualifies.
        if bk is not None and np.isfinite(bk):
            ax.axvline(bk, color=colors["dip"], lw=1.4, ls="--", zorder=4)
            if annotate_dip and len(pk) >= 2:
                left_pk = max([p for p in pk if p <= bk], default=pk[0])
                right_pk = min([p for p in pk if p >= bk], default=pk[-1])
                thr = frac * min(y_at(left_pk), y_at(right_pk))
                ax.plot([left_pk, min(right_pk, x_hi)], [thr, thr],
                        color=colors["threshold"], lw=1.2, ls=(0, (2, 1.5)),
                        zorder=5)
                ax.plot(bk, y_at(bk), "o", ms=5, color=colors["dip"], zorder=7)
                ax.annotate(
                    f"dip < {frac:g}x\nsmaller flank",
                    xy=(bk, y_at(bk)), xytext=(bk + 0.06 * x_hi, thr + 0.32),
                    fontsize=6.5, color=colors["threshold"], ha="left",
                    va="bottom", zorder=8,
                    arrowprops=dict(arrowstyle="->", color=colors["threshold"],
                                    lw=0.8))

        if applied_break is not None and np.isfinite(applied_break):
            ax.axvline(applied_break, color=colors["applied"], lw=2.0, zorder=6)
            ax.text(applied_break, 0.97, f" {applied_break:.1f} m",
                    color=colors["applied"], fontsize=7, ha="left", va="top",
                    zorder=9, transform=ax.get_xaxis_transform())

        ax.set_xlim(0, x_hi)
        ax.set_ylim(0, 1.35)
        # Reserve title pad for the descriptive break line (inherited break /
        # accepted dip / accepted median), placed just below the reach title.
        ax.set_title(title, fontsize=8, pad=14 if subtitle else 6)
        if subtitle:
            ax.annotate(subtitle, xy=(0.5, 1.0), xytext=(0, 2),
                        xycoords="axes fraction", textcoords="offset points",
                        ha="center", va="bottom", fontsize=6.6, color="0.25",
                        annotation_clip=False)
        ax.tick_params(labelsize=7)
        ax.set_yticks([])

    npan = 1 + len(wcs)
    n_cols = max(1, min(n_cols, npan))
    n_rows = (npan + n_cols - 1) // n_cols
    fig, axarr = plt.subplots(
        n_rows, n_cols, figsize=(3.4 * n_cols, 2.4 * n_rows), squeeze=False)
    axf = [axarr[i // n_cols][i % n_cols] for i in range(n_rows * n_cols)]

    gbk = break_hand(gl) if len(gl) else None
    gtitle = ("whole corridor"
              + (f": break {gbk:.1f} m" if gbk is not None
                 else ": no single break\n-> reach-adaptive used"))
    draw_panel(axf[0], gl, gtitle, show_raw=True)

    for i, wc in enumerate(wcs, start=1):
        sub = wins[wins["window_centre_m"] == wc]
        abrk, own = applied_at(wc)
        dip = break_hand(sub)
        # Three cases: a window with no own break INHERITS; a window with an
        # own dip that survives post-processing unchanged is ACCEPTED; a window
        # whose applied break differs from its own dip was replaced by the
        # median-of-3 across neighbouring reaches (or the min-break clamp).
        differs = (dip is not None and abrk is not None
                   and abs(dip - abrk) > 0.15)
        if own is None or not own:
            tag = "inherited"
            subtitle = (f"inherited break at {abrk:.1f} metres"
                        if abrk is not None else "inherited break")
        elif differs:
            tag = "median"
            subtitle = (f"dip at {dip:.1f} metres accepted "
                        f"median at {abrk:.1f} metres")
        else:
            tag = "own"
            acc = dip if dip is not None else abrk
            subtitle = f"accepted dip at {acc:.1f} metres"
        fr = fraction_at(wc)
        if abs(fr - valley_fraction) > 1e-6:
            subtitle = (subtitle or "") + f"  (fraction {fr:.2f})"
        draw_panel(axf[i], sub, f"reach @ {wc / 1000:.1f} km  ({tag})",
                   show_raw=False, applied_break=abrk, subtitle=subtitle,
                   fraction=fr)

    for j in range(npan, n_rows * n_cols):
        axf[j].set_axis_off()

    handles = [
        Line2D([0], [0], color=colors["curve"], lw=1.8,
               label="smoothed decision metric"),
        Line2D([0], [0], marker="v", color=colors["peak"], lw=0, ms=6,
               label="mode peak (surviving)"),
        Line2D([0], [0], marker="v", mfc="none", mec=colors["raw_peak"],
               mew=1.0, lw=0, ms=6, label="candidate peak (merged away)"),
        Line2D([0], [0], color=colors["boundary"], lw=0.7, ls=":",
               label="mode boundary (inter-mode minimum)"),
        Line2D([0], [0], color=colors["dip"], lw=1.4, ls="--",
               label="chosen dip"),
        Line2D([0], [0], color=colors["threshold"], lw=1.2, ls=(0, (2, 1.5)),
               label="dip acceptance threshold"),
        Line2D([0], [0], color=colors["applied"], lw=2.0,
               label="applied break (clamp/median/inherit)"),
        plt.Rectangle((0, 0), 1, 1, color=colors["valley"], alpha=0.5,
                      label="valley (HAND <= break)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=8, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(
        f"channel {source_id}: how the valley HAND break is chosen  "
        f"(break_valley_fraction = {valley_fraction:g})\n"
        f"(x = height above nearest channel, LINEAR & cropped to the decision "
        f"region; valley = at or below the first inter-mode dip)", fontsize=11)
    # X labels on the bottom row of drawn panels.
    for ax in axf[max(0, npan - n_cols):npan]:
        ax.set_xlabel("HAND (m)", fontsize=8)
    fig.tight_layout(rect=[0, 0.03, 1, 0.94])
    return fig, axf[:npan]


def plot_valley_reach_map(
    reach_id_path: str,
    hillshade_path: str,
    *,
    nodes_df: Optional[pd.DataFrame] = None,
    source_id: int = 0,
    source_col: str = "source_id",
    x_col: str = "x",
    y_col: str = "y",
    cmap_name: str = "tab20",
    alpha: float = 0.5,
    centreline_color: str = "crimson",
    label_reaches: bool = True,
    crop_margin_frac: float = 0.06,
    ax=None,
):
    """Map of the along-channel reaches over a hillshade.

    Draws the per-channel reach-id raster
    (``*_valley_reach_id_ch<sid>.tif`` from ``lsdtt-valley-centreline`` with
    ``write_valley_reach_id_raster``) semi-transparently over the hillshade,
    so every DEM pixel is coloured by the reach (window) of its nearest
    channel node — i.e. exactly the population each per-reach HAND break
    histogram is built from. Reaches are labelled with their number at the
    centroid of their pixels, and the valley centreline is overlaid.

    Each corridor pixel is attributed to its nearest channel node (Euclidean),
    the channel is segmented into evenly spaced reaches by along-channel
    distance, and the pixel takes that node's reach number. Non-corridor
    pixels are nodata (any negative value) and not drawn.

    Args:
        reach_id_path: Path to the per-channel reach-id GeoTIFF. Values >= 0
            are reach ids; negatives are outside the corridor (not drawn).
        hillshade_path: Path to the hillshade GeoTIFF (drawn in grey behind).
        nodes_df: Optional valley-centreline nodes table
            (``*_valley_centreline_nodes.csv``) to overlay the centreline.
        source_id: Channel source (used only to filter ``nodes_df`` and title).
        source_col, x_col, y_col: Column names in ``nodes_df``.
        cmap_name: Qualitative colormap cycled across reach ids so adjacent
            reaches differ.
        alpha: Opacity of the reach overlay.
        centreline_color: Colour of the overlaid centreline.
        label_reaches: Label each reach with its number at its pixel centroid.
        crop_margin_frac: Fraction of the reach extent added as a margin
            around the crop.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``.

    Raises:
        ValueError: if the reach raster has no valid (>= 0) pixels.
    """
    plt = _import_mpl()
    rasterio = _import_rasterio()
    import numpy as np

    with rasterio.open(reach_id_path) as src:
        reach = src.read(1).astype(float)
        rtr = src.transform
        rb = src.bounds
        nd = src.nodata
    if nd is not None:
        reach[reach == nd] = np.nan
    reach[reach < 0] = np.nan          # negatives = outside the corridor
    valid = np.isfinite(reach)
    if not valid.any():
        raise ValueError(
            f"Reach raster {reach_id_path} has no valid (>= 0) pixels.")

    with rasterio.open(hillshade_path) as src:
        hs = src.read(1).astype(float)
        hb = src.bounds
        hnd = src.nodata
    if hnd is not None:
        hs[hs == hnd] = np.nan

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 9))
    else:
        fig = ax.figure

    ax.imshow(hs, cmap="gray", extent=[hb.left, hb.right, hb.bottom, hb.top],
              origin="upper", interpolation="nearest", zorder=0)

    # Reach overlay: colour each reach id from a cycled qualitative colormap.
    ids = np.unique(reach[valid]).astype(int)
    cmap = plt.get_cmap(cmap_name)
    ncol = cmap.N
    rgba = np.zeros((reach.shape[0], reach.shape[1], 4), dtype=float)
    for rid in ids:
        m = reach == rid
        c = cmap(int(rid) % ncol)
        rgba[m] = (c[0], c[1], c[2], alpha)
    ax.imshow(rgba, extent=[rb.left, rb.right, rb.bottom, rb.top],
              origin="upper", interpolation="nearest", zorder=1)

    # Centreline overlay.
    if nodes_df is not None:
        g = nodes_df[nodes_df[source_col] == source_id]
        if len(g):
            ax.plot(g[x_col].to_numpy(), g[y_col].to_numpy(), "-",
                    lw=1.6, color=centreline_color, zorder=3)

    # Reach-number labels at each reach's pixel centroid (world coords).
    if label_reaches:
        for rid in ids:
            rows, cols = np.where(reach == rid)
            cr, cc = rows.mean(), cols.mean()
            wx, wy = rtr * (cc + 0.5, cr + 0.5)
            ax.text(wx, wy, str(int(rid)), fontsize=9, fontweight="bold",
                    ha="center", va="center", color="black", zorder=4,
                    path_effects=None,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white",
                              ec="none", alpha=0.7))

    # Crop to the reach extent (+ margin), in world coordinates.
    rows, cols = np.where(valid)
    xs, ys = rtr * (cols + 0.5, rows + 0.5)
    xmin, xmax = float(np.min(xs)), float(np.max(xs))
    ymin, ymax = float(np.min(ys)), float(np.max(ys))
    mx = crop_margin_frac * (xmax - xmin + 1e-6)
    my = crop_margin_frac * (ymax - ymin + 1e-6)
    ax.set_xlim(xmin - mx, xmax + mx)
    ax.set_ylim(ymin - my, ymax + my)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title(
        f"channel {source_id}: along-channel reaches "
        f"(pixels coloured by nearest-node reach; numbers = reach id)",
        fontsize=11)
    return fig, ax


def plot_valley_floor_map(
    valley_mask_path: str,
    hillshade_path: str,
    *,
    nodes_df: Optional[pd.DataFrame] = None,
    source_col: str = "source_id",
    x_col: str = "x",
    y_col: str = "y",
    valley_color: str = "#3a86ff",
    by_source: bool = True,
    cmap_name: str = "tab10",
    alpha: float = 0.5,
    centreline_color: str = "crimson",
    crop_margin_frac: float = 0.06,
    title: Optional[str] = None,
    ax=None,
):
    """Map of the extracted valley floor over a hillshade.

    Drapes the valley-mask raster (``*_valley_mask.tif`` from
    ``lsdtt-valley-centreline``; 0 = not valley, ``source_id + 1`` inside each
    channel's valley) semi-transparently over the hillshade and overlays the
    centreline(s). Handy for comparing how the extracted valley floor changes
    with a parameter such as ``break_valley_fraction``.

    Args:
        valley_mask_path: Path to the combined valley-mask GeoTIFF.
        hillshade_path: Path to the hillshade GeoTIFF (drawn in grey behind).
        nodes_df: Optional valley-centreline nodes table
            (``*_valley_centreline_nodes.csv``) to overlay the centrelines.
        source_col, x_col, y_col: Column names in ``nodes_df``.
        valley_color: Fill colour when ``by_source`` is False (single colour).
        by_source: Colour each channel's valley by its source id (cycled
            ``cmap_name``); otherwise a single ``valley_color``.
        cmap_name: Qualitative colormap used when ``by_source``.
        alpha: Opacity of the valley overlay.
        centreline_color: Colour of the overlaid centrelines.
        crop_margin_frac: Fraction of the valley extent added as a margin.
        title: Figure title; a default is used when None.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``.

    Raises:
        ValueError: if the mask has no valley (>0) pixels.
    """
    plt = _import_mpl()
    rasterio = _import_rasterio()
    import numpy as np
    from matplotlib.colors import to_rgb

    with rasterio.open(valley_mask_path) as src:
        vm = src.read(1).astype(float)
        vtr = src.transform
        vb = src.bounds
        nd = src.nodata
    if nd is not None:
        vm[vm == nd] = 0.0
    valley = vm > 0
    if not valley.any():
        raise ValueError(
            f"Valley mask {valley_mask_path} has no valley (>0) pixels.")

    with rasterio.open(hillshade_path) as src:
        hs = src.read(1).astype(float)
        hb = src.bounds
        hnd = src.nodata
    if hnd is not None:
        hs[hs == hnd] = np.nan

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 9))
    else:
        fig = ax.figure

    ax.imshow(hs, cmap="gray", extent=[hb.left, hb.right, hb.bottom, hb.top],
              origin="upper", interpolation="nearest", zorder=0)

    rgba = np.zeros((vm.shape[0], vm.shape[1], 4), dtype=float)
    if by_source:
        cmap = plt.get_cmap(cmap_name)
        for sid in np.unique(vm[valley]).astype(int):
            c = cmap((sid - 1) % cmap.N)
            m = vm == sid
            rgba[m] = (c[0], c[1], c[2], alpha)
    else:
        c = to_rgb(valley_color)
        rgba[valley] = (c[0], c[1], c[2], alpha)
    ax.imshow(rgba, extent=[vb.left, vb.right, vb.bottom, vb.top],
              origin="upper", interpolation="nearest", zorder=1)

    if nodes_df is not None and len(nodes_df):
        for sid, g in nodes_df.groupby(source_col):
            ax.plot(g[x_col].to_numpy(), g[y_col].to_numpy(), "-",
                    lw=1.4, color=centreline_color, zorder=3)

    rows, cols = np.where(valley)
    xs, ys = vtr * (cols + 0.5, rows + 0.5)
    xmin, xmax = float(np.min(xs)), float(np.max(xs))
    ymin, ymax = float(np.min(ys)), float(np.max(ys))
    mx = crop_margin_frac * (xmax - xmin + 1e-6)
    my = crop_margin_frac * (ymax - ymin + 1e-6)
    ax.set_xlim(xmin - mx, xmax + mx)
    ax.set_ylim(ymin - my, ymax + my)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title(title if title is not None else "valley floor", fontsize=11)
    return fig, ax


def plot_valley_break_fraction(
    breaks_df: pd.DataFrame,
    *,
    source_id: int = 0,
    source_col: str = "source_id",
    frac_min: Optional[float] = None,
    frac_max: Optional[float] = None,
    ax=None,
):
    """Per-reach relief-adaptive acceptance fraction along the channel.

    Shows, for one channel, how the relief-adaptive method varied the dip
    acceptance fraction with along-channel distance, alongside the local
    relief that drove it and the resulting valley HAND break. Reads the
    ``*_valley_breaks.csv`` table (written by ``lsdtt-valley-centreline``),
    which for the ``relief_adaptive`` method carries ``break_fraction`` and
    ``relief_m`` columns in addition to ``break_m`` / ``own_break``.

    Three stacked panels vs distance (km): (1) the acceptance fraction with
    own-break reaches marked and the [frac_min, frac_max] band shaded;
    (2) the median local relief per reach; (3) the applied HAND break.
    """
    import matplotlib.pyplot as plt

    df = breaks_df[breaks_df[source_col] == source_id].copy()
    if df.empty:
        raise ValueError(f"No break rows for source_id={source_id}.")
    if "break_fraction" not in df.columns:
        raise ValueError(
            "breaks_df has no 'break_fraction' column — run "
            "lsdtt-valley-centreline with valley_break_method=relief_adaptive.")
    df = df.sort_values("window_centre_m")
    d = df["window_centre_m"].to_numpy(dtype=float) / 1000.0
    frac = df["break_fraction"].to_numpy(dtype=float)
    own = (df["own_break"].astype(str).str.lower() == "true").to_numpy()
    relief = (df["relief_m"].to_numpy(dtype=float)
              if "relief_m" in df.columns else np.full(d.shape, np.nan))
    brk = df["break_m"].to_numpy(dtype=float)

    if ax is None:
        fig, axs = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    else:
        fig = ax.figure
        axs = [ax, ax.twinx(), ax.twinx()]

    fm = frac_min if frac_min is not None else float(np.nanmin(frac))
    fM = frac_max if frac_max is not None else float(np.nanmax(frac))
    axs[0].axhspan(fm, fM, color="#7b2cbf", alpha=0.10, zorder=0,
                   label=f"[{fm:g}, {fM:g}] fraction band")
    axs[0].plot(d, frac, "-", color="#7b2cbf", lw=1.6, zorder=2)
    axs[0].plot(d[own], frac[own], "o", ms=6, color="#7b2cbf",
                label="own break", zorder=3)
    axs[0].plot(d[~own], frac[~own], "o", ms=5, mfc="none",
                mec="#7b2cbf", label="inherited", zorder=3)
    axs[0].set_ylabel("acceptance\nfraction")
    axs[0].legend(fontsize=7, loc="best")
    axs[0].set_title(f"channel {source_id}: relief-adaptive break fraction",
                     fontsize=11)

    axs[1].plot(d, relief, "-", color="#2a7f3f", lw=1.6)
    axs[1].fill_between(d, 0, relief, color="#2a7f3f", alpha=0.15)
    axs[1].set_ylabel("median local\nrelief (m)")

    axs[2].plot(d, brk, "-", color="crimson", lw=1.6)
    axs[2].plot(d[own], brk[own], "o", ms=5, color="crimson")
    axs[2].set_ylabel("applied HAND\nbreak (m)")
    axs[2].set_xlabel("distance along channel (km)")

    for a in axs:
        a.grid(True, alpha=0.25)
    fig.tight_layout()
    return fig, axs
