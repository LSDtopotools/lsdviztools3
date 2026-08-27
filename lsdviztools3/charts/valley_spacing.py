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


def plot_spacing_vs_ci(frame, ax=None, title=None):
    """Plot spacing against the channelization index for a sweep.

    Args:
        frame: from
            :func:`~lsdviztools3.io.valley_spacing.load_valley_spacing_sweep`.
        ax: existing axes, or None to create a figure.
        title: figure title.

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
    ax.axvspan(CI_CRITICAL_LOW, CI_CRITICAL_HIGH, color="#cccccc", alpha=0.55,
               label=f"CI$_c$ {CI_CRITICAL_LOW:.0f}-{CI_CRITICAL_HIGH:.0f} "
                     "(Bonetti 2020)")

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
