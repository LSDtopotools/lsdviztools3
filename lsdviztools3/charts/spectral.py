"""Polar plot of significant 2D spectral frequencies.

:func:`plot_significant_frequencies_polar` draws each significant cell of the
2D power spectrum (``lsdtt-spectral-analysis``,
``_significant_frequencies_2d.csv``, or the legacy ``_significant_freq.csv``)
on a North-up, clockwise polar axis:

* angle  -- the bearing, either of the wave vector (``orientation="wave"``, as
  written by the tool) or of the crests (``orientation="crest"``, wave + 90).
  For regularly spaced bands of vegetation the crest bearing is the direction
  the bands run;
* radius -- frequency (cycles/m), with rings labelled by wavelength so the
  spacing can be read off directly;
* colour -- normalised PSD (observed / background), log scale.

This is the spectrum seen from above: a characteristic spacing shows as a ring
segment at one radius, and a preferred orientation as a pair of opposed
clusters. The spectrum of a real surface is point-symmetric, so every cell has
a twin at bearing + 180; ``half=True`` keeps one of each pair.

Significance is decided by lsdtt3; nothing is recomputed here.
"""

from __future__ import annotations

import numpy as np

from ..io.spectral import SignificantFrequencies2D

_MPL_MSG = (
    "Spectral charts need matplotlib. Install the optional 'charts' extra, e.g. "
    "`pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Candidate ring wavelengths (m); the ones inside the plotted radius are used.
_RING_WAVELENGTHS = (0.1, 0.2, 0.25, 0.5, 1.0, 2.0, 2.5, 5.0, 10.0, 20.0,
                     25.0, 50.0, 100.0, 200.0, 250.0, 500.0, 1000.0, 2000.0,
                     5000.0)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt, LogNorm


def _format_wavelength(w: float) -> str:
    return f"{w:g} m"


def _ring_ticks(r_max: float, radius: str, n_rings: int = 4):
    """Ring positions and labels (always in wavelength) inside r_max.

    Rings are spread evenly in the plotted radius: for each of n_rings target
    radii, the nearest round wavelength is used. Long wavelengths all sit near
    the centre of a frequency plot, so choosing by wavelength would pile the
    rings up there.
    """
    to_r = (lambda w: 1.0 / w) if radius == "frequency" else (lambda w: w)
    candidates = [(to_r(w), w) for w in _RING_WAVELENGTHS if 0 < to_r(w) <= r_max]
    chosen = []
    for target in r_max * np.arange(1, n_rings + 1) / n_rings:
        if not candidates:
            break
        r, w = min(candidates, key=lambda rw: abs(rw[0] - target))
        if w not in chosen and r >= 0.1 * r_max:
            chosen.append(w)
    chosen.sort(key=to_r)
    return [to_r(w) for w in chosen], [_format_wavelength(w) for w in chosen]


def plot_significant_frequencies_polar(
    sig,
    *,
    orientation: str = "wave",
    radius: str = "frequency",
    half: bool = False,
    r_max: float | None = None,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    marker_size: float = 14.0,
    colorbar: bool = True,
    title: str | None = None,
    ax=None,
):
    """Polar scatter of significant 2D spectrum cells.

    Args:
        sig: a :class:`~lsdviztools3.io.spectral.SignificantFrequencies2D`, or a
            DataFrame with ``bearing_deg``, ``frequency_cycles_per_m``,
            ``wavelength_m`` and ``normalised_psd`` columns (e.g. one stem of
            :func:`~lsdviztools3.io.spectral.load_significant_frequencies_2d_sweep`).
        orientation: ``"wave"`` plots the wave-vector bearing as written by
            the tool; ``"crest"`` adds 90 deg, giving the direction the crests
            (bands, ridges) run.
        radius: ``"frequency"`` (cycles/m; short wavelengths at the rim, as in
            the 2D PSD raster) or ``"wavelength"`` (m; long wavelengths at the
            rim). Rings are labelled in wavelength either way.
        half: keep only bearings in [0, 180), one cell of each symmetric pair.
        r_max: outer radius in the units of ``radius``. Defaults to 5% beyond
            the largest plotted value. For ``radius="frequency"`` passing the
            Nyquist frequency, 1 / (2 * cell size), shows how close the cells
            sit to the resolution limit.
        cmap: colormap for normalised PSD.
        vmin, vmax: colour limits for normalised PSD; default to the data range
            (``vmin`` defaults to the threshold when it is known). Fix them to
            compare several runs on one scale.
        marker_size: scatter marker area (points^2).
        colorbar: draw a colorbar.
        title: axis title; defaults to the stem, the cell count and threshold.
        ax: an existing polar ``Axes``, or None to create one.

    Returns:
        ``(fig, ax)`` with ``ax`` a North-up, clockwise ``PolarAxes``. A run
        with no significant cells gives empty axes carrying a note -- that is a
        result, not an error.

    Raises:
        ValueError: on an unknown ``orientation`` or ``radius``, or a table
            missing the required columns.
    """
    plt, LogNorm = _import_mpl()
    if orientation not in ("wave", "crest"):
        raise ValueError(f"orientation must be 'wave' or 'crest', got {orientation!r}.")
    if radius not in ("frequency", "wavelength"):
        raise ValueError(f"radius must be 'frequency' or 'wavelength', got {radius!r}.")

    if isinstance(sig, SignificantFrequencies2D):
        table, stem, threshold = sig.table, sig.stem, sig.threshold
    else:
        table, stem, threshold = sig, "", float("nan")
    need = {"bearing_deg", "frequency_cycles_per_m", "wavelength_m", "normalised_psd"}
    if not need.issubset(table.columns):
        raise ValueError(f"table is missing columns {sorted(need - set(table.columns))}")

    bearing = np.asarray(table["bearing_deg"], dtype=float)
    if orientation == "crest":
        bearing = bearing + 90.0
    bearing = np.mod(bearing, 360.0)
    r_col = "frequency_cycles_per_m" if radius == "frequency" else "wavelength_m"
    r = np.asarray(table[r_col], dtype=float)
    power = np.asarray(table["normalised_psd"], dtype=float)

    keep = np.isfinite(bearing) & np.isfinite(r) & np.isfinite(power) & (power > 0)
    if half:
        keep &= bearing < 180.0
    bearing, r, power = bearing[keep], r[keep], power[keep]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6.0, 5.5), subplot_kw={"projection": "polar"})
    else:
        fig = ax.figure
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    if r_max is None:
        r_max = 1.05 * float(r.max()) if r.size else 1.0
    ax.set_rlim(0.0, r_max)
    ticks, labels = _ring_ticks(r_max, radius)
    ax.set_rticks(ticks)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_rlabel_position(22.5)

    if r.size:
        order = np.argsort(power)  # strongest drawn last, on top
        lo = vmin if vmin is not None else (
            threshold if np.isfinite(threshold) and threshold < power.min()
            else power.min())
        hi = vmax if vmax is not None else max(power.max(), lo * 1.0001)
        pts = ax.scatter(np.deg2rad(bearing[order]), r[order], c=power[order],
                         s=marker_size, cmap=cmap, norm=LogNorm(vmin=lo, vmax=hi),
                         edgecolors="none")
        if colorbar:
            cb = fig.colorbar(pts, ax=ax, pad=0.1, shrink=0.8)
            cb.set_label("normalised PSD (observed / background)")
    else:
        note = "no cells above threshold"
        if np.isfinite(threshold):
            note += f" {threshold:.3g}"
        ax.text(0.5, 0.5, note, transform=ax.transAxes, ha="center",
                va="center", fontsize=9, color="0.4")

    if title is None:
        parts = [stem] if stem else []
        parts.append(f"{int(keep.sum())} cells")
        if np.isfinite(threshold):
            parts.append(f"threshold {threshold:.3g}")
        parts.append("crest bearing" if orientation == "crest" else "wave-vector bearing")
        title = "  |  ".join(parts)
    if title:
        ax.set_title(title, fontsize=9, pad=14)
    return fig, ax
