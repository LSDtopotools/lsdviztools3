"""Valley-spacing charts.

Two figures, and they answer different questions:

* :func:`plot_valley_spectrum` -- ONE landscape. Shows the averaged power
  spectrum with the detected peaks, the power-weighted mean wavelength and the
  convex half-width marked on it. The point is to let the reader see *why* the
  two measures agree or disagree, rather than trusting a pair of scalars.

* :func:`plot_spacing_vs_ci` -- A SWEEP. Spacing against the Bonetti
  channelization index, with the critical band shaded.

**Why both measures are always drawn.** They fail differently, and the
disagreement is the diagnostic. On a domain holding only a few valleys the
spectral peak can land only on domain harmonics (extent/n), so it jumps between
discrete values and looks like noise; the convex half-width uses thousands of
samples and keeps working. At high channelization the reverse risk appears: the
dominant spectral power sits at the domain-scale drainage structure rather than
the valley scale, so the spectral wavelength saturates while the half-width
continues to fall. A chart showing only one of them would hide both failures.

Measures come from lsdtt3 (``lsdtt_valley_spacing_cpu``); nothing is recomputed
here.

References
----------
Perron, J. T., Dietrich, W. E., & Kirchner, J. W. (2008). Controls on the
spacing of first-order valleys. *JGR Earth Surface* 113, F04016.

Bonetti, S., Hooshyar, M., Camporeale, C., & Porporato, A. (2020).
Channelization cascade in landscape evolution. *PNAS* 117(3), 1375-1382.

Hooshyar, M., et al. (2021). Spectral signature of landscape channelization.
*Geophysical Research Letters*.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ..io.valley_spacing import CI_CRITICAL_HIGH, CI_CRITICAL_LOW


def plot_valley_spectrum(vs, ax=None, title=None):
    """Plot one landscape's spectrum with both spacing measures marked.

    Args:
        vs: a :class:`~lsdviztools3.io.valley_spacing.ValleySpacing`.
        ax: existing axes, or None to create a figure.
        title: figure title; defaults to the output stem.

    Returns:
        The matplotlib Figure.

    Raises:
        ValueError: if the set carries no spectrum.
    """
    if vs.spectrum is None or not len(vs.spectrum):
        raise ValueError("no spectrum in this output set "
                         "(was compute_valley_spacing enabled?)")

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
    else:
        fig = ax.figure

    lam = vs.spectrum["wavelength_m"].to_numpy()
    P = vs.spectrum["power"].to_numpy()
    good = (lam > 0) & (P > 0)
    ax.loglog(lam[good], P[good], color="0.35", lw=1.0, label="mean power spectrum")

    # Mark the peaks the weighted mean was taken over.
    if good.sum() > 3:
        p = P[good]
        is_peak = np.r_[False, (p[1:-1] > p[:-2]) & (p[1:-1] > p[2:]), False]
        ax.loglog(lam[good][is_peak], p[is_peak], ".", ms=5,
                  color="#1f77b4", label=f"peaks ({is_peak.sum()})")

    if np.isfinite(vs.spectral_wavelength) and vs.spectral_wavelength > 0:
        ax.axvline(vs.spectral_wavelength, color="#1f77b4", lw=1.8,
                   label=f"spectral, power-weighted: {vs.spectral_wavelength:.0f} m")
    if np.isfinite(vs.convex_wavelength) and vs.convex_wavelength > 0:
        ax.axvline(vs.convex_wavelength, color="#d62728", lw=1.8, ls="--",
                   label=f"convex half-width x2: {vs.convex_wavelength:.0f} m")

    ax.set_xlabel("wavelength (m)")
    ax.set_ylabel("mean power")
    ax.set_title(title or f"Valley spacing: {vs.stem}")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", alpha=0.25)

    # State the caveat on the figure, not just in the docs -- a reader looking at
    # the image alone should know when the spectral measure cannot be trusted.
    if vs.under_sampled:
        ax.text(0.02, 0.02,
                f"only ~{vs.valleys_in_domain:.0f} valleys across the domain:\n"
                "the spectral peak can only land on domain harmonics.\n"
                "Use the convex half-width here.",
                transform=ax.transAxes, fontsize=8, va="bottom",
                bbox=dict(boxstyle="round", fc="#fff3cd", ec="#856404", alpha=0.9))

    fig.tight_layout()
    return fig


def plot_spacing_vs_ci(frame, ax=None, title=None, bonetti_lambda_c=None):
    """Plot spacing against the channelization index for a sweep.

    Args:
        frame: from
            :func:`~lsdviztools3.io.valley_spacing.load_valley_spacing_sweep`.
        ax: existing axes, or None to create a figure.
        title: figure title.
        bonetti_lambda_c: valley spacing Bonetti et al. report at the critical
            index (~42 m in their configuration). Drawn as a reference line when
            given. It is only comparable if the domain and exponents match
            theirs, so it is opt-in rather than assumed.

    Returns:
        The matplotlib Figure.

    Raises:
        ValueError: if the frame is empty.
    """
    if frame is None or not len(frame):
        raise ValueError("empty sweep frame")

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
    else:
        fig = ax.figure

    ci = frame["channelization_index"].to_numpy()

    # The published transition, so a reader can see which runs are even in the
    # channelized regime.
    #
    # The axis limits are set FIRST and forced to bracket the critical band. On a
    # log axis a band from 32 to 37 is narrow, and if the lowest run sits at
    # CI = 31 the band lands hard against the left spine -- the legend then shows
    # an entry with no visible swatch, which is worse than not drawing it, because
    # the reader cannot check the comparison the legend claims to offer.
    lo = min(float(np.nanmin(ci)), CI_CRITICAL_LOW) / 1.6
    hi = max(float(np.nanmax(ci)), CI_CRITICAL_HIGH) * 1.6
    ax.set_xlim(lo, hi)

    ax.axvspan(CI_CRITICAL_LOW, CI_CRITICAL_HIGH, color="#f0c674", alpha=0.75,
               zorder=0,
               label=f"CI$_c$ = {CI_CRITICAL_LOW:.0f}–{CI_CRITICAL_HIGH:.0f}, first "
                     "channelization\n(Bonetti et al. 2020)")
    ax.axvline(CI_CRITICAL_LOW, color="#b58900", lw=0.8, zorder=0)
    ax.axvline(CI_CRITICAL_HIGH, color="#b58900", lw=0.8, zorder=0)

    # Bonetti's own reported spacing at the transition, so the comparison the
    # shaded band invites can actually be made rather than merely suggested.
    if bonetti_lambda_c is not None:
        ax.axhline(bonetti_lambda_c, color="#b58900", ls=":", lw=1.4, zorder=0,
                   label=f"Bonetti $\\lambda_c$ ≈ {bonetti_lambda_c:.0f} m at CI$_c$")

    # The distance-transform hillslope length is the primary measure: it is the
    # only one of the three that survives branching (see the module docstring).
    if "hillslope_length" in frame.columns and frame["hillslope_length"].notna().any():
        ax.loglog(ci, frame["hillslope_length"], "o-", color="#d62728",
                  label="hillslope length (distance transform)")
    else:
        ax.loglog(ci, frame["convex_wavelength"], "o-", color="#d62728",
                  label="convex half-width x2")

    # Distinguish the spectral points that are under-sampled: they are plotted,
    # because hiding them would hide the failure, but they are not joined by a
    # line that implies a trend.
    us = frame["under_sampled"].to_numpy(dtype=bool)
    lam_s = frame["spectral_wavelength"].to_numpy()
    ax.loglog(ci[~us], lam_s[~us], "s-", color="#1f77b4",
              label="spectral, power-weighted")
    if us.any():
        ax.loglog(ci[us], lam_s[us], "s", mfc="none", color="#1f77b4",
                  label="spectral (too few valleys)")

    ax.set_xlabel("channelization index  CI = $K\\,l^{m+n}/(D^n U^{1-n})$")
    ax.set_ylabel("valley spacing (m)")
    ax.set_title(title or "Valley spacing vs channelization index")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    return fig


def plot_landscape_panel(rasters, labels=None, ncols=None, title=None,
                         spacings=None, spacing_label="L$_h$",
                         azimuth=315.0, altitude=45.0):
    """Hillshade every landscape in a sweep, side by side, for a sniff test.

    A scaling law fitted to numbers can look convincing while the landscapes it
    came from are nonsense -- all channel, or a single valley, or a grid
    artefact. This exists so that judgement can be made by eye before any number
    is believed.

    Args:
        rasters: paths to the DEMs, in the order to display.
        labels: per-panel captions; defaults to file stems.
        ncols: panels per row; defaults to a near-square layout.
        title: figure title.
        spacings: optional per-panel length in metres. When given, a horizontal
            bar of that length is drawn TO SCALE on each panel and labelled, so
            the measured value can be held against the ridges actually visible.
            That comparison is the whole point of the panel.
        spacing_label: what the bar represents, used in its label. Defaults to
            "L_h"; pass e.g. "spacing" if the bar is a valley spacing rather than
            a hillslope length.
        azimuth, altitude: illumination, presentation only.

    Returns:
        The matplotlib Figure.

    Raises:
        ValueError: if no rasters were given.
    """
    import rasterio

    paths = list(rasters)
    if not paths:
        raise ValueError("no rasters given")

    if labels is None:
        labels = [str(p).split("/")[-1] for p in paths]
    n = len(paths)
    if ncols is None:
        ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.6 * ncols, 3.8 * nrows),
                             squeeze=False)

    for idx, path in enumerate(paths):
        ax = axes[idx // ncols][idx % ncols]
        with rasterio.open(path) as src:
            z = src.read(1).astype(float)
            cell = abs(src.transform.a)
        z[z < -1e30] = np.nan

        ax.imshow(_hillshade_local(z, cell, azimuth, altitude),
                  cmap="gray", vmin=0, vmax=1, interpolation="nearest")
        ax.set_title(labels[idx], fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])

        if spacings is not None and idx < len(spacings):
            lam = spacings[idx]
            if lam and np.isfinite(lam):
                # Drawn TO SCALE against the hillshade, so the reader can hold the
                # measured number against the ridges actually visible. That
                # comparison is the entire point of the panel, and an unlabelled
                # bar invites the reader to guess what it represents.
                px = lam / cell
                y = z.shape[0] * 0.90
                x0 = z.shape[1] * 0.06
                ax.plot([x0, x0 + px], [y, y], color="#d62728", lw=3.0,
                        solid_capstyle="butt", zorder=5)
                for xe in (x0, x0 + px):                    # end ticks
                    ax.plot([xe, xe], [y - z.shape[0] * 0.018,
                                       y + z.shape[0] * 0.018],
                            color="#d62728", lw=1.5, zorder=5)
                ax.text(x0, y - z.shape[0] * 0.030,
                        f"{spacing_label} = {lam:.0f} m",
                        color="#d62728", fontsize=7.5, va="bottom", zorder=5,
                        bbox=dict(boxstyle="square,pad=0.15", fc="white",
                                  ec="none", alpha=0.75))

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    if title:
        # tight_layout does not know about suptitle, so the panel titles ride up
        # into it. Lay out first, then reserve a strip at the top.
        fig.tight_layout()
        fig.subplots_adjust(top=1.0 - 0.4 / max(1.0, nrows * 3.8))
        fig.suptitle(title, fontsize=11)
    else:
        fig.tight_layout()
    return fig


def _hillshade_local(z, cell, azimuth=315.0, altitude=45.0):
    """Horn hillshade, GDAL's convention.

    aspect = atan2(dz/drow, -dz/dcol). The sign is easy to get wrong because row
    increases southward in a north-up raster; getting it wrong inverts the image
    so ridges read as valleys, which is exactly the error this panel exists to
    catch elsewhere.
    """
    dy, dx = np.gradient(z, cell)
    slope = np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(dy, -dx)
    az = np.radians(360.0 - azimuth + 90.0)
    alt = np.radians(altitude)
    shade = (np.sin(alt) * np.cos(slope)
             + np.cos(alt) * np.sin(slope) * np.cos(az - aspect))
    return np.clip(shade, 0.0, 1.0)

