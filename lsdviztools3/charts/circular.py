"""Circular (directional) statistics and comparison plots.

Compare two angular samples (e.g. channel bearings from two DEMs):

- :func:`circular_wasserstein` — 1-Wasserstein (earth-mover) distance on a
  circle, the natural distance for bearing data (handles 359 deg ~ 1 deg).
- :func:`kuiper_two_sample` — Kuiper's test, the rotation-invariant circular
  analogue of the two-sample Kolmogorov-Smirnov test (statistic + p-value).
- :func:`plot_circular_comparison` — overlay two samples' PDFs/CDFs and annotate
  the distance.

The numerics use only numpy (a core dependency); only the plot needs matplotlib
(the optional ``charts`` extra). Mirrors the C++ circular tests in lsdtt3's
``lsdtt_statistical_tests``.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np

_MPL_MSG = (
    "Circular comparison plots need matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def _prep(values, period: float) -> np.ndarray:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        raise ValueError("Angular sample has no finite values.")
    return np.sort(np.mod(v, period))


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cumw = np.cumsum(w)
    return float(v[np.searchsorted(cumw, 0.5 * cumw[-1])])


def _circular_w1(a_sorted: np.ndarray, b_sorted: np.ndarray, period: float) -> float:
    """Circular 1-Wasserstein for two empirical samples (sorted, in [0, period))."""
    support = np.unique(np.concatenate([a_sorted, b_sorted]))
    gaps = np.empty(support.size)
    gaps[:-1] = np.diff(support)
    gaps[-1] = period - support[-1] + support[0]   # wraparound interval
    f1 = np.searchsorted(a_sorted, support, side="right") / a_sorted.size
    f2 = np.searchsorted(b_sorted, support, side="right") / b_sorted.size
    d = f1 - f2
    # Optimal transport offset on the circle = length-weighted median of (F1-F2).
    mu = _weighted_median(d, gaps)
    return float(np.sum(np.abs(d - mu) * gaps))


def circular_wasserstein(
    a: Sequence[float],
    b: Sequence[float],
    *,
    period: float = 360.0,
    rotation_invariant: bool = False,
    rotation_step: float = 1.0,
) -> float:
    """1-Wasserstein distance between two angular samples on a circle.

    Args:
        a, b: Angular samples (degrees by default).
        period: Circle period (360 for full bearings, 180 for axial/orientation).
        rotation_invariant: Also minimise over a rigid rotation of ``b`` (use
            when only the *shape* matters, not the absolute orientation). Default
            ``False`` keeps the geographic frame.
        rotation_step: Rotation grid step (degrees) when ``rotation_invariant``.

    Returns:
        Distance in the same units as the inputs (degrees). 0 for identical
        distributions; bounded by ``period / 2``.
    """
    a_s, b_s = _prep(a, period), _prep(b, period)
    base = _circular_w1(a_s, b_s, period)
    if not rotation_invariant:
        return base
    best = base
    for rot in np.arange(rotation_step, period, rotation_step):
        shifted = np.sort(np.mod(b_s + rot, period))
        best = min(best, _circular_w1(a_s, shifted, period))
    return best


def kuiper_two_sample(
    a: Sequence[float],
    b: Sequence[float],
    *,
    period: float = 360.0,
) -> Tuple[float, float]:
    """Two-sample Kuiper test (rotation-invariant circular KS analogue).

    Args:
        a, b: Angular samples.
        period: Circle period (unused in the statistic, which is origin-free;
            kept for API symmetry).

    Returns:
        ``(V, p_value)`` where ``V = D+ + D-`` and ``p_value`` is the asymptotic
        Kuiper significance (small ``p`` => distributions differ).
    """
    a_s, b_s = _prep(a, period), _prep(b, period)
    support = np.unique(np.concatenate([a_s, b_s]))
    f1 = np.searchsorted(a_s, support, side="right") / a_s.size
    f2 = np.searchsorted(b_s, support, side="right") / b_s.size
    d = f1 - f2
    v = float(d.max() - d.min())   # D+ + D-

    n1, n2 = a_s.size, b_s.size
    ne = n1 * n2 / (n1 + n2)
    sne = np.sqrt(ne)
    lam = (sne + 0.155 + 0.24 / sne) * v
    if lam <= 0:
        return v, 1.0
    j = np.arange(1, 101)
    terms = (4.0 * j * j * lam * lam - 1.0) * np.exp(-2.0 * j * j * lam * lam)
    p = float(np.clip(2.0 * terms.sum(), 0.0, 1.0))
    return v, p


def plot_circular_comparison(
    a: Sequence[float],
    b: Sequence[float],
    *,
    labels: Tuple[str, str] = ("A", "B"),
    bin_width_deg: float = 15.0,
    kind: str = "both",
    period: float = 360.0,
    annotate: bool = True,
    rotation_invariant: bool = False,
    styles: Optional[dict] = None,
    title: Optional[str] = None,
    ax=None,
):
    """Overlay two angular samples' PDFs and/or CDFs; annotate the distance.

    Args:
        a, b: Angular samples (degrees).
        labels: Legend labels for ``a`` and ``b``.
        bin_width_deg: PDF histogram bin width.
        kind: ``"pdf"``, ``"cdf"``, or ``"both"`` (a 1x2 figure).
        period: Circle period.
        annotate: Write the circular-Wasserstein distance (and Kuiper V/p) on the
            plot.
        rotation_invariant: Passed to :func:`circular_wasserstein`.
        styles: Optional ``{label: matplotlib kwargs}`` for each series' line.
        title: Figure suptitle (``both``) or axes title.
        ax: Existing ``Axes`` (for single ``kind``) or ``None``.

    Returns:
        ``(fig, ax)`` for a single kind, or ``(fig, (ax_pdf, ax_cdf))`` for
        ``"both"``.
    """
    plt = _import_mpl()
    if kind not in ("pdf", "cdf", "both"):
        raise ValueError(f"kind must be 'pdf', 'cdf' or 'both', got {kind!r}.")
    styles = styles or {}

    series = {labels[0]: _prep(a, period), labels[1]: _prep(b, period)}
    dist = circular_wasserstein(a, b, period=period,
                                rotation_invariant=rotation_invariant)
    vstat, pval = kuiper_two_sample(a, b, period=period)
    note = (f"circular $W_1$ = {dist:.1f}$\\degree$\n"
            f"Kuiper V = {vstat:.3f}, p = {pval:.3g}")

    def _draw_pdf(axp):
        edges = np.arange(0.0, period + bin_width_deg, bin_width_deg)
        for lab, v in series.items():
            axp.hist(v, bins=edges, density=True, histtype="step",
                     linewidth=1.6, label=lab, **styles.get(lab, {}))
        axp.set_xlabel("Bearing (deg)")
        axp.set_ylabel("Density")
        axp.set_xlim(0, period)
        axp.legend(fontsize=8)

    def _draw_cdf(axc):
        for lab, v in series.items():
            x = np.concatenate([[0.0], v, [period]])
            y = np.concatenate([[0.0], np.arange(1, v.size + 1) / v.size, [1.0]])
            axc.step(x, y, where="post", linewidth=1.6, label=lab,
                     **styles.get(lab, {}))
        axc.set_xlabel("Bearing (deg)")
        axc.set_ylabel("Cumulative fraction")
        axc.set_xlim(0, period)
        axc.legend(fontsize=8)

    if kind == "both":
        fig, (ax_pdf, ax_cdf) = plt.subplots(1, 2, figsize=(11, 4.5))
        _draw_pdf(ax_pdf)
        _draw_cdf(ax_cdf)
        if annotate:
            ax_cdf.text(0.04, 0.96, note, transform=ax_cdf.transAxes,
                        va="top", ha="left", fontsize=9,
                        bbox=dict(boxstyle="round", fc="white", alpha=0.8))
        if title:
            fig.suptitle(title)
        return fig, (ax_pdf, ax_cdf)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4.5))
    else:
        fig = ax.figure
    (_draw_pdf if kind == "pdf" else _draw_cdf)(ax)
    if annotate:
        ax.text(0.04, 0.96, note, transform=ax.transAxes, va="top", ha="left",
                fontsize=9, bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    if title:
        ax.set_title(title)
    return fig, ax
