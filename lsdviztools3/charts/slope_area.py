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
from rasters lsdtt3 wrote.
"""

from __future__ import annotations

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
