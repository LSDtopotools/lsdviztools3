"""Slope-area charts: the "boomerang" that separates hillslopes from channels.

Plot local gradient against upslope drainage area on log-log axes and a
process signature appears, because the two transport laws scale differently
with area:

* **Hillslopes.** Diffusive transport carries a flux proportional to gradient,
  so at steady state gradient must GROW with distance from the divide -- and
  distance grows with area. The small-area limb therefore RISES to the right.
* **Channels.** Detachment-limited incision balances uplift when
  ``U = K A^m S^n``, so ``S = (U/K)^(1/n) A^(-m/n)``: a straight line FALLING
  to the right with slope ``-m/n``.

The join between them is the **channel head**, and the peak of the resulting
"boomerang" is the clearest picture of where it sits. A landscape with no
hillslope process shows no rising limb at all, only the fluvial line continued
to the smallest areas -- which is exactly what a fluvial-only model produces,
and is the quickest way to see whether a hillslope solver did anything.

Two cautions that decide whether the picture means anything:

* **Bin in log space and plot a median, not a mean.** Slope-area clouds span
  orders of magnitude in area and are heavily skewed at fixed area; an
  arithmetic mean of a skewed distribution sits above its own data and moves the
  apparent rollover.
* **The rollover is not the channel head of any particular pixel.** It is where
  the population changes character. Reading a threshold area off it is a
  reasonable estimate and not a measurement.

Howard, A. D. and Kerby, G. (1983). Channel changes in badlands. *Geological
Society of America Bulletin* 94:739-752.
Tarboton, D. G., Bras, R. L. and Rodriguez-Iturbe, I. (1991). On the extraction
of channel networks from digital elevation data. *Hydrological Processes*
5:81-100.
Montgomery, D. R. and Foufoula-Georgiou, E. (1993). Channel network source
representation using digital elevation models. *Water Resources Research*
29:3925-3934.

Nothing here computes a landscape quantity: slope and drainage area are read
from rasters lsdtt3 wrote. The second half of the module draws the binned
tables, low-concavity points and segments that ``lsdtt-slope-area`` writes
(read with :mod:`lsdviztools3.io.slope_area`); the bins, gradients and
concavities in them are lsdtt3's.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np

_MPL_MSG = (
    "The slope-area charts need matplotlib. Install it e.g. "
    "`pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def _read(path):
    """Read one raster band, with nodata as NaN."""
    try:
        import rasterio
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Reading rasters needs rasterio.") from exc

    with rasterio.open(path) as src:
        band = src.read(1).astype(float)
        nodata = src.nodata
    if nodata is not None:
        band = np.where(band == nodata, np.nan, band)
    # lsdtt3's usual sentinel, in case the file does not declare one.
    return np.where(band <= -9990.0, np.nan, band)


def binned_median(area, slope, n_bins: int = 40, min_count: int = 8):
    """Median slope in log-spaced area bins.

    Args:
        area: Upslope area per pixel, any consistent unit.
        slope: Gradient per pixel (dimensionless, not degrees).
        n_bins: Number of log-spaced bins across the area range.
        min_count: Bins with fewer pixels than this are dropped rather than
            plotted, because a median of three points in the tail reads as
            structure.

    Returns:
        ``(centres, median, q25, q75, count)`` as arrays over the surviving
        bins.
    """
    area = np.asarray(area, dtype=float)
    slope = np.asarray(slope, dtype=float)

    good = np.isfinite(area) & np.isfinite(slope) & (area > 0) & (slope > 0)
    area, slope = area[good], slope[good]
    if area.size == 0:
        raise ValueError("no pixels with positive, finite area and slope")

    edges = np.logspace(np.log10(area.min()), np.log10(area.max()), n_bins + 1)
    idx = np.digitize(area, edges) - 1

    centres, med, q25, q75, count = [], [], [], [], []
    for b in range(n_bins):
        sel = idx == b
        n = int(sel.sum())
        if n < min_count:
            continue
        s = slope[sel]
        centres.append(np.sqrt(edges[b] * edges[b + 1]))   # geometric centre
        med.append(np.median(s))
        q25.append(np.percentile(s, 25))
        q75.append(np.percentile(s, 75))
        count.append(n)

    return (np.array(centres), np.array(med), np.array(q25),
            np.array(q75), np.array(count))


def rollover_area(centres, median, tol: float = 0.10):
    """The hillslope/channel join, as an INTERVAL rather than a point.

    The naive answer is ``argmax`` of the binned median. That is fragile in
    exactly the case that matters: a real slope-area curve does not peak
    sharply, it flattens into a broad plateau where diffusive and fluvial
    transport are comparable, and ``argmax`` then picks an arbitrary bin off the
    noise. On the test landscape it chose 678 where the plateau ran to ~3000 --
    a factor of four, reported as if it were a measurement.

    So this returns the whole span of bins whose median is within ``tol`` of the
    peak. The width of that span is the honest uncertainty on the transition,
    and a wide one is a statement about the landscape (the two processes overlap
    over a decade of area) rather than a failure of the estimate.

    THE INTERVAL DEPENDS ON ``tol`` AND MUST BE QUOTED WITH IT. That is why the
    figure labels the tolerance rather than printing a bare pair of numbers: the
    plateau has no intrinsic edge, so any interval is a choice about how flat
    counts as flat.

    Args:
        centres: Geometric bin centres.
        median: Binned median slope.
        tol: Relative tolerance on the peak defining the plateau.

    Returns:
        ``(low, high, peak)`` areas, or None when the median never turns over --
        which is itself the answer: no rising limb means no hillslopes to
        separate.
    """
    if centres.size < 3:
        return None
    k = int(np.argmax(median))
    if k == 0 or k == centres.size - 1:
        return None          # the peak is at an end: no interior turnover

    on_plateau = median >= median[k] * (1.0 - tol)
    idx = np.flatnonzero(on_plateau)
    # Only the run containing the peak, so a second bump elsewhere is excluded.
    lo = k
    while lo - 1 in idx:
        lo -= 1
    hi = k
    while hi + 1 in idx:
        hi += 1
    return float(centres[lo]), float(centres[hi]), float(centres[k])


def plot_slope_area(slope_raster,
                    area_raster,
                    *,
                    channel_mask_raster=None,
                    channel_area_threshold: Optional[float] = None,
                    n_bins: int = 40,
                    m_over_n: Optional[float] = None,
                    mark_area: Optional[float] = None,
                    max_scatter: int = 40000,
                    title: Optional[str] = None,
                    figsize: Tuple[float, float] = (8.0, 6.5),
                    seed: int = 0):
    """The slope-area "boomerang".

    Args:
        slope_raster: Path to a gradient raster (dimensionless).
        area_raster: Path to an upslope drainage-area raster.
        channel_area_threshold: Alternative to a mask raster -- select channel
            pixels as ``area >= threshold``, in the units of ``area_raster``.
            This is a DISPLAY-time selection, not a computation: it reproduces
            the definition the model was run with so the two limbs can be
            coloured, and it is the same comparison lsdtt3 makes internally. Use
            the mask raster instead whenever one exists.
        channel_mask_raster: Optional channel mask; when given, hillslope and
            channel pixels are drawn in different colours so the reader can see
            whether the extracted network agrees with where the scaling breaks.
            The two are independent -- the mask comes from an extraction
            algorithm, the rollover from the topography -- so disagreement is
            informative rather than an error.
        n_bins: Log-spaced area bins for the median.
        mark_area: Optional area to mark with a vertical line -- typically the
            channel threshold the model was RUN with, so the reader can compare
            it against the transition the topography actually produced. The two
            are independent and need not agree.
        m_over_n: When given, overlay a reference line of slope ``-m/n``
            anchored on the fluvial limb. This is a REFERENCE, not a fit: it
            shows what concavity the model was asked for, so a mismatch is
            visible without anyone fitting anything.
        max_scatter: Cap on plotted points, for legibility and file size.
        seed: Sampling seed, so the same figure is drawn on a rerun.

    Returns:
        The matplotlib Figure.
    """
    plt = _import_mpl()

    slope = _read(slope_raster).ravel()
    area = _read(area_raster).ravel()
    if slope.size != area.size:
        raise ValueError(
            f"slope raster has {slope.size} cells but area raster has "
            f"{area.size}; they must be on the same grid")

    chan = None
    if channel_mask_raster is not None:
        chan = _read(channel_mask_raster).ravel()
        if chan.size != slope.size:
            raise ValueError("the channel mask is not on the same grid")
    elif channel_area_threshold is not None:
        chan = (area >= channel_area_threshold).astype(float)

    good = np.isfinite(slope) & np.isfinite(area) & (area > 0) & (slope > 0)

    fig, ax = plt.subplots(figsize=figsize)

    rng = np.random.default_rng(seed)

    def scatter(sel, colour, label):
        n = int(sel.sum())
        if n == 0:
            return
        i = np.flatnonzero(sel)
        if n > max_scatter:
            i = rng.choice(i, max_scatter, replace=False)
        ax.scatter(area[i], slope[i], s=2, alpha=0.12, color=colour,
                   linewidths=0, label=f"{label} ({n:,} px)")

    if chan is not None:
        scatter(good & (chan > 0.5), "#1f77b4", "channel")
        scatter(good & ~(chan > 0.5), "#8c8c8c", "hillslope")
    else:
        scatter(good, "#8c8c8c", "pixels")

    centres, med, q25, q75, _ = binned_median(area[good], slope[good], n_bins)
    ax.fill_between(centres, q25, q75, color="#d62728", alpha=0.18,
                    label="binned IQR")
    ax.plot(centres, med, "-o", color="#d62728", ms=4, lw=1.8,
            label="binned median")

    roll = rollover_area(centres, med)
    roll_hi = None
    if roll is not None:
        lo, hi, peak = roll
        roll_hi = hi
        ax.axvspan(lo, hi, color="#d62728", alpha=0.10)
        ax.axvline(peak, color="#d62728", ls=":", lw=1.4)
        ax.annotate(f"transition {lo:,.0f}-{hi:,.0f}\n(within 10% of peak)",
                    xy=(peak, med.max()), xytext=(6, 6),
                    textcoords="offset points", fontsize=9, color="#d62728")

    if mark_area is not None:
        ax.axvline(mark_area, color="#2a7f3f", ls="-.", lw=1.6)
        ax.annotate(f"imposed threshold {mark_area:,.0f}",
                    xy=(mark_area, np.nanmin(q25)), xytext=(6, 10),
                    textcoords="offset points", fontsize=9, color="#2a7f3f")

    if m_over_n is not None and roll_hi is not None:
        # Anchor the reference line on the binned median one decade above the
        # rollover, i.e. squarely on the fluvial limb.
        anchor = centres > roll_hi * 3.0
        if anchor.any():
            a0 = centres[anchor][0]
            s0 = med[anchor][0]
            xs = centres[centres >= roll_hi]
            ax.plot(xs, s0 * (xs / a0) ** (-m_over_n), "--", color="k", lw=1.4,
                    label=f"reference slope $-m/n$ = {-m_over_n:g}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("upslope drainage area")
    ax.set_ylabel("gradient  $S$")
    ax.grid(True, which="both", alpha=0.25)

    if title is None:
        title = "Slope-area: hillslopes rise to the right, channels fall"
    ax.set_title(title)

    leg = ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    for handle in leg.legend_handles:
        try:
            handle.set_alpha(1.0)
            handle.set_sizes([18])
        except Exception:  # noqa: BLE001 - line handles have no set_sizes
            pass

    fig.tight_layout()
    return fig


# =============================================================================
# lsdtt-slope-area outputs: binned tables, low-concavity points and segments
# =============================================================================
#
# These read what lsdtt-slope-area wrote (lsdviztools3.io.slope_area) and draw
# it. The bins, gradients and concavities are lsdtt3's; nothing is recomputed.
#
# ``plot_binned_slope_area``             d50 with a d25-d75 band, one line per
#     polyfit lambda / D8 window
# ``plot_slope_area_method_comparison``  the binned d50 of two runs (e.g.
#     polyfit and D8 flow path) on one plot
# ``plot_theta_distribution``            histogram of theta for channel pixels
#     with a valid fit, with the +/- threshold lines
# ``plot_low_concavity_slope_area``      the channel-pixel slope-area cloud,
#     low-concavity pixels in a second colour
# ``render_low_concavity_map``           the segments on a hillshade, coloured
#     by ``low_concavity``

#: Line colours of the binned tables, in fixed order (as charts.knickpoints).
SERIES_COLOURS = ("#2a78d6", "#1baf7a", "#9a5bd6", "#eb6834", "#d6a01b", "#d6457a")
#: low_concavity classes: 1, 0 and -1 (no valid fit).
LOW_CONCAVITY_COLOUR = "#eb6834"
NOT_LOW_CONCAVITY_COLOUR = "#2a78d6"
NO_FIT_COLOUR = "#9a9a9a"
#: The same colours as GMT r/g/b, for the map.
_GMT_CLASS = {1: ("235/104/52", "low concavity"), 0: ("42/120/214", "not low concavity"),
              -1: ("150/150/150", "no valid fit")}


def _draw_binned(ax, tables, *, colours, linestyle="-", band_alpha=0.18, label_prefix=""):
    """One d50 line and d25-d75 band per table. Returns the d50 lines."""
    lines = []
    for k, t in enumerate(tables):
        df = t["table"]
        colour = colours[k % len(colours)]
        a = df["drainage_area"].to_numpy(dtype=float)
        keep = a > 0
        a = a[keep]
        d25 = df["d25_slope"].to_numpy(dtype=float)[keep]
        d50 = df["d50_slope"].to_numpy(dtype=float)[keep]
        d75 = df["d75_slope"].to_numpy(dtype=float)[keep]
        pos = (d25 > 0) & (d75 > 0)
        if band_alpha > 0 and pos.any():
            ax.fill_between(a[pos], d25[pos], d75[pos], color=colour, alpha=band_alpha, lw=0)
        good = d50 > 0
        (line,) = ax.plot(a[good], d50[good], linestyle, marker="o", ms=3.5, lw=1.8, color=colour,
                          label=f"{label_prefix}{t['label']}")
        lines.append(line)
    return lines


def _finish_log_axes(ax, title):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("drainage area (m$^2$)")
    ax.set_ylabel("gradient $S$ (m/m)")
    ax.grid(True, which="both", alpha=0.25)
    if title:
        ax.set_title(title)


def plot_binned_slope_area(tables, *, title: Optional[str] = None,
                           figsize: Tuple[float, float] = (7.5, 5.5)):
    """Binned median gradient against drainage area, one line per scale.

    Args:
        tables: ``[{"label", "table"}, ...]`` from
            :func:`lsdviztools3.io.slope_area.find_binned_slope_area`; each
            ``table`` has ``drainage_area``, ``d25_slope``, ``d50_slope``,
            ``d75_slope`` as lsdtt-slope-area wrote them.
        title: Axes title (default: none).

    Returns:
        ``(fig, ax)``. The d50 is a line with markers; the shaded band is the
        d25-d75 range of the same bins.
    """
    plt = _import_mpl()
    if not tables:
        raise ValueError("no binned slope-area tables to plot")
    fig, ax = plt.subplots(figsize=figsize)
    _draw_binned(ax, tables, colours=SERIES_COLOURS)
    _finish_log_axes(ax, title)
    ax.legend(title="binned d50, band d25-d75", fontsize=8, title_fontsize=8,
              loc="best", framealpha=0.9)
    fig.tight_layout()
    return fig, ax


def plot_slope_area_method_comparison(tables, compare_tables, *, name: str = "run",
                                      compare_name: str = "comparison",
                                      title: Optional[str] = None, bands: bool = False,
                                      figsize: Tuple[float, float] = (7.5, 5.5)):
    """The binned d50 of two runs on one log-log plot.

    The first run's lines are solid, the comparison run's dashed. Colours
    run in one fixed order over both sets. The d25-d75 bands are off by
    default (``bands=True`` draws them, lighter than in
    :func:`plot_binned_slope_area`): four or more overlapping bands hide the
    lines.

    Args:
        tables, compare_tables: As for :func:`plot_binned_slope_area`.
        name, compare_name: Run names, put in front of each label.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    if not tables or not compare_tables:
        raise ValueError("the method comparison needs binned tables from both runs")
    fig, ax = plt.subplots(figsize=figsize)
    n = len(tables)
    alpha = 0.1 if bands else 0.0
    _draw_binned(ax, tables, colours=SERIES_COLOURS, band_alpha=alpha, label_prefix=f"{name}: ")
    _draw_binned(ax, compare_tables, colours=SERIES_COLOURS[n:] + SERIES_COLOURS[:n],
                 linestyle="--", band_alpha=alpha, label_prefix=f"{compare_name}: ")
    _finish_log_axes(ax, title)
    ax.legend(title="binned d50, band d25-d75" if bands else "binned d50", fontsize=8, title_fontsize=8,
              loc="best", framealpha=0.9)
    fig.tight_layout()
    return fig, ax


def _class_counts(flag):
    flag = np.asarray(flag)
    return int((flag == 1).sum()), int((flag == 0).sum()), int((flag == -1).sum())


def plot_theta_distribution(theta, flag, threshold: float, *, theta_range=None, bins: int = 60,
                            title: Optional[str] = None,
                            figsize: Tuple[float, float] = (7.5, 4.8)):
    """Histogram of the fitted concavity theta of the channel pixels.

    Only pixels with a valid fit (``flag >= 0``) are counted. Vertical lines
    mark ``-threshold`` and ``+threshold`` (lsdtt3's
    ``low_concavity_max_abs_concavity``); the bars inside them are the
    low-concavity pixels.

    Args:
        theta: Theta per channel pixel, as lsdtt3 wrote it.
        flag: ``low_concavity`` per pixel: 1, 0, -1 (no valid fit).
        threshold: ``low_concavity_max_abs_concavity`` of the run.
        theta_range: ``(low, high)`` of the x axis. Default: the 2nd to 98th
            percentile of the valid theta (taken at data values, so fewer
            than ~50 pixels are all shown), widened to include 2 x threshold
            either side. Values outside are not drawn, and their number is
            given in the legend.
        bins: Number of histogram bins across ``theta_range``.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    theta = np.asarray(theta, dtype=float)
    flag = np.asarray(flag)
    if theta.shape != flag.shape:
        raise ValueError("theta and flag must have the same length")
    n_low, n_not, n_none = _class_counts(flag)
    n_valid = n_low + n_not
    valid = (flag >= 0) & np.isfinite(theta)
    t = theta[valid]

    if theta_range is None:
        if t.size:
            # nearest data values, so a small sample is shown whole
            lo = np.percentile(t, 2, method="lower")
            hi = np.percentile(t, 98, method="higher")
        else:
            lo, hi = -threshold, threshold
        lo, hi = min(lo, -2 * threshold), max(hi, 2 * threshold)
    else:
        lo, hi = (float(v) for v in theta_range)
        if not hi > lo:
            raise ValueError("theta_range must be (low, high) with high > low")
    inside = (t >= lo) & (t <= hi)
    n_out = int((~inside).sum())

    fig, ax = plt.subplots(figsize=figsize)
    edges = np.linspace(lo, hi, bins + 1)
    if t.size:
        low = flag[valid] == 1          # lsdtt3's flag, not re-derived from theta
        ax.hist([t[inside & ~low], t[inside & low]], bins=edges, stacked=True,
                color=[NOT_LOW_CONCAVITY_COLOUR, LOW_CONCAVITY_COLOUR],
                label=[f"not low concavity ({n_not:,})", f"low concavity ({n_low:,})"],
                edgecolor="white", linewidth=0.4)
    for sign in (-1, 1):
        ax.axvline(sign * threshold, color="k", ls="--", lw=1.2,
                   label=f"|$\\theta$| = {threshold:g}" if sign == 1 else None)
    # Counts with no bar of their own, as legend entries.
    from matplotlib.patches import Patch
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Patch(color=NO_FIT_COLOUR))
    labels.append(f"no valid fit, not drawn ({n_none:,})")
    if n_out:
        handles.append(Patch(color="none"))
        labels.append(f"valid, outside the x range ({n_out:,})")
    ax.legend(handles, labels, fontsize=8, loc="upper right", framealpha=0.9)
    ax.set_xlim(lo, hi)
    ax.set_xlabel(r"concavity $\theta$ ($S \propto A^{-\theta}$)")
    ax.set_ylabel("channel pixels")
    ax.grid(True, axis="y", alpha=0.25)
    share = f"{n_low / n_valid:.1%}" if n_valid else "n/a"
    ax.set_title(title if title is not None else
                 f"{n_low:,} of {n_valid:,} valid channel pixels low concavity ({share}); "
                 f"{n_none:,} no valid fit", fontsize=10)
    fig.tight_layout()
    return fig, ax


def plot_low_concavity_slope_area(points: "pd.DataFrame", *, title: Optional[str] = None,
                                  figsize: Tuple[float, float] = (7.5, 5.5)):
    """The channel-pixel slope-area cloud, coloured by ``low_concavity``.

    Each point is one channel pixel of ``<prefix>_low_concavity_points.csv``
    (``drainage_area_m2`` against ``slope``, the D8 flow-path gradient of
    the low-concavity fit). Low-concavity pixels are drawn last, in a second
    colour; pixels with no valid fit are muted. Pixels without a positive
    area and gradient cannot go on log axes and are counted in the title.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    a = points["drainage_area_m2"].to_numpy(dtype=float)
    s = points["slope"].to_numpy(dtype=float)
    f = points["low_concavity"].to_numpy()
    drawable = np.isfinite(a) & np.isfinite(s) & (a > 0) & (s > 0)
    fig, ax = plt.subplots(figsize=figsize)
    for value, colour, label, size, alpha in (
            (-1, NO_FIT_COLOUR, "no valid fit", 5, 0.35),
            (0, NOT_LOW_CONCAVITY_COLOUR, "not low concavity", 7, 0.55),
            (1, LOW_CONCAVITY_COLOUR, "low concavity", 10, 0.9)):
        sel = drawable & (f == value)
        n_all = int((f == value).sum())
        ax.scatter(a[sel], s[sel], s=size, color=colour, alpha=alpha, linewidths=0,
                   label=f"{label} ({int(sel.sum()):,} of {n_all:,})")
    _finish_log_axes(ax, None)
    ax.set_ylabel("D8 flow-path gradient $S$ (m/m)")
    n_hidden = int((~drawable).sum())
    ax.set_title(title if title is not None else
                 f"channel pixels; {n_hidden:,} without a positive gradient not drawn",
                 fontsize=10)
    leg = ax.legend(fontsize=8, loc="best", framealpha=0.9, title="pixels drawn",
                    title_fontsize=8)
    for handle in leg.legend_handles:
        try:
            handle.set_alpha(1.0)
        except Exception:  # noqa: BLE001
            pass
    fig.tight_layout()
    return fig, ax


def stream_order_pen_width(order, base: float = 0.8, step: float = 0.6,
                           maximum: float = 3.2) -> float:
    """Map line width (points) for a Strahler order: ``base + step*(order-1)``,
    capped at ``maximum``."""
    return float(min(maximum, base + step * max(int(order) - 1, 0)))


def render_low_concavity_map(dem_path, segments, png_path, *, width_by_order: bool = True,
                             line_width: float = 1.2, coordinates: str = "local",
                             high_res: bool = False,
                             legend_position: str = "JTR+jTR+o0.2c",
                             title: Optional[str] = None):
    """The low-concavity segments on a hillshade of ``dem_path``.

    Segments are coloured by ``low_concavity``: low concavity (orange), not
    (blue), no valid fit (grey), drawn in that order from the bottom so the
    low-concavity reaches are on top. With ``width_by_order`` the line width
    grows with ``stream_order`` (:func:`stream_order_pen_width`), else every
    line is ``line_width`` points.

    Args:
        dem_path: DEM for the hillshade (the run's input DEM).
        segments: ``GeoDataFrame`` or path of ``<prefix>_low_concavity_segments.fgb``.
        png_path: Output PNG.

    Returns:
        ``Path`` to the written PNG.
    """
    import io as _io

    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from ..io.slope_area import load_low_concavity_segments
    from ..render.base import grid_to_display, save_figure
    from ..render.channels import render_channels
    from ..render.decorations import add_scalebar
    from ..render.hillshade import render_hillshade
    from ..render.style import MapStyle

    width = MapStyle.MAX_WIDTH
    gdf = segments if not isinstance(segments, (str, Path)) else load_low_concavity_segments(segments)
    if gdf.empty:
        raise ValueError("no low-concavity segments to draw")
    da = to_coordinates(load_raster(dem_path), coordinates)
    crs = get_crs(da)
    _, extent = grid_to_display(to_gmt_grid(da), crs)
    fig = render_hillshade(str(dem_path), coordinates=coordinates, shade=True,
                           style=MapStyle(cmap="gray", figure_size=width, transparent=False))
    has_order = width_by_order and "stream_order" in gdf.columns
    rows = []
    for value in (-1, 0, 1):
        colour, label = _GMT_CLASS[value]
        part = gdf[gdf["low_concavity"] == value]
        if part.empty:
            continue
        groups = (part.groupby("stream_order") if has_order else [(None, part)])
        for order, g in groups:
            w = stream_order_pen_width(order) if order is not None else line_width
            fig = render_channels(g, coordinates=coordinates,
                                  style=MapStyle(line_width=f"{w:g}p", line_color=colour,
                                                 figure_size=width), fig=fig)
        rows.append(f"S 0.3c - 0.5c - 2p,{colour} 0.9c {label} ({len(part):,} segments)")
    if has_order:
        lo, hi = int(gdf["stream_order"].min()), int(gdf["stream_order"].max())
        rows.append(f"L 8p - L line width: stream order {lo}-{hi}" if hi > lo
                    else f"L 8p - L stream order {lo}")
    if title:
        rows.insert(0, f"L 9p - L {title}")
    fig.legend(spec=_io.StringIO("\n".join(rows)), position=legend_position,
               box="+gwhite@20+p0.5p,black")
    add_scalebar(fig, extent, crs)
    save_figure(fig, png_path, style=MapStyle(figure_size=width, transparent=False), high_res=high_res)
    return Path(png_path)
