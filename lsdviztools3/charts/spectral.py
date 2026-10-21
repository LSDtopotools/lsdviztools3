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

:func:`plot_spectral_tile_summary` puts one run on a 2 x 2 figure: the input
raster, the normalised 2D spectrum (:func:`plot_normalised_spectrum_2d`), the
radial spectrum (:func:`plot_radial_spectrum`) and the polar plot above. Load
the run with :func:`~lsdviztools3.io.spectral.load_spectral_analysis`.

Significance is decided by lsdtt3; nothing is recomputed here.
"""

from __future__ import annotations

import numpy as np

from ..io.spectral import BINNED_THRESHOLD, SignificantFrequencies2D

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


def _threshold_2d(run) -> float:
    sig = run.significant_2d
    if sig is not None and np.isfinite(sig.threshold):
        return sig.threshold
    return -np.log(1.0 - 0.9999)


def plot_normalised_spectrum_2d(run, *, ax=None, threshold=None, vmin=0.1,
                                vmax=None, cmap="magma", colorbar=True):
    """Normalised 2D spectrum (observed / background), north up, in cycles/m.

    Cells above the significance threshold are outlined in cyan. The rasters
    are placed by their own georeferencing, which lsdtt3 writes north-up with
    the zero-frequency cell centred on (0, 0).

    Args:
        run: a :class:`~lsdviztools3.io.spectral.SpectralAnalysisRun` with both
            2D rasters.
        ax: existing axes, or None to create a figure.
        threshold: outline level; defaults to the run's 2D threshold (9.21 at
            the CLI default).
        vmin, vmax: log colour limits; vmax defaults to max(data, 10).
        cmap: colormap.
        colorbar: draw a colorbar.

    Returns:
        ``(fig, ax)``.

    Raises:
        ValueError: if the run lacks either 2D raster
            (``write_2d_spectrum_rasters: false``).
    """
    plt, LogNorm = _import_mpl()
    from ..io.raster import get_extent

    norm = run.normalised_psd_2d
    if norm is None:
        raise ValueError("this run has no 2D PSD rasters "
                         "(was write_2d_spectrum_rasters switched off?)")
    if threshold is None:
        threshold = _threshold_2d(run)
    z = np.asarray(norm.values, dtype=float).copy()
    ny, nx = z.shape
    z[ny // 2, nx // 2] = np.nan  # zero frequency: detrended away, not data

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.5, 4.8))
    else:
        fig = ax.figure
    top = max(np.nanmax(z), 10.0) if vmax is None else vmax
    im = ax.imshow(z, extent=get_extent(norm), origin="upper", cmap=cmap,
                   norm=LogNorm(vmin=vmin, vmax=top), interpolation="nearest")
    ax.contour(np.asarray(norm.x), np.asarray(norm.y), np.nan_to_num(z) > threshold,
               levels=[0.5], colors="cyan", linewidths=0.6)
    ax.set_xlabel("east frequency (cycles/m)")
    ax.set_ylabel("north frequency (cycles/m)")
    ax.set_aspect("equal")
    ax.set_title(f"normalised 2D spectrum (cyan: > {threshold:.3g})", fontsize=9)
    if colorbar:
        fig.colorbar(im, ax=ax, shrink=0.8, label="observed / background")
    return fig, ax


def plot_radial_spectrum(run, *, ax=None, show_normalised=True):
    """Log-binned radial spectrum: observed and background PSD, and their ratio.

    Observed and background PSD go on log axes; the normalised PSD (green, right
    axis) is drawn with the 2.9955 line the CLI uses for
    ``_significant_frequencies.csv``. The rollover frequency is marked.

    Args:
        run: a :class:`~lsdviztools3.io.spectral.SpectralAnalysisRun` with
            ``radial`` (``write_radial_spectrum: true``).
        ax: existing axes, or None to create a figure.
        show_normalised: draw the normalised PSD on a twin axis.

    Returns:
        ``(fig, ax)``; ``ax`` is the PSD axis.

    Raises:
        ValueError: if the run has no radial spectrum.
    """
    plt, _ = _import_mpl()
    rad = run.radial
    if rad is None or not len(rad):
        raise ValueError("this run has no radial spectrum "
                         "(set write_radial_spectrum: true)")
    if ax is None:
        fig, ax = plt.subplots(figsize=(5.5, 4.5))
    else:
        fig = ax.figure

    f = rad["frequency_cycles_per_m"]
    ax.loglog(f, rad["psd"], "o-", ms=3, color="tab:blue", label="observed")
    ax.loglog(f, rad["background_psd"], "-", color="0.4", label="background")
    rollover = run.rollover.get("rollover_frequency_cycles_per_m", float("nan"))
    if isinstance(rollover, float) and np.isfinite(rollover) and rollover > 0:
        ax.axvline(rollover, ls=":", color="tab:red",
                   label=f"rollover {1.0 / rollover:.3g} m")
    ax.set_xlabel("frequency (cycles/m)")
    ax.set_ylabel("PSD (m$^2$ / (cycles/m))")
    ax.legend(fontsize=7, loc="lower left")

    n_sig = int((rad["normalised_psd"] > BINNED_THRESHOLD).sum())
    title = f"radial spectrum: {n_sig} bins > {BINNED_THRESHOLD}"
    beta = run.rollover.get("rollover_beta")
    if isinstance(beta, float) and np.isfinite(beta):
        title += f", beta {beta:.2f}"
    ax.set_title(title, fontsize=9)

    if show_normalised:
        ax2 = ax.twinx()
        ax2.semilogx(f, rad["normalised_psd"], color="tab:green", lw=1)
        ax2.axhline(BINNED_THRESHOLD, color="tab:green", ls="--", lw=0.8)
        ax2.set_ylabel("normalised PSD (binned)", color="tab:green")
    return fig, ax


def _plot_input(run, ax, cmap, label):
    plt, _ = _import_mpl()
    from ..io.raster import get_extent

    da = run.input_raster
    if da is None:
        ax.text(0.5, 0.5, "input raster not found", transform=ax.transAxes,
                ha="center", va="center", color="0.4")
        ax.set_axis_off()
        return
    w, e, s, n = get_extent(da)
    im = ax.imshow(np.asarray(da.values, dtype=float), cmap=cmap,
                   extent=[0.0, e - w, 0.0, n - s], origin="upper",
                   interpolation="nearest")
    ax.set_xlabel("m from west edge")
    ax.set_ylabel("m from south edge")
    ax.set_title(label, fontsize=9)
    ax.figure.colorbar(im, ax=ax, shrink=0.8)


def plot_spectral_tile_summary(run, *, input_cmap="YlGn", input_label="input",
                               polar_r_max=None, polar_vmin=None, polar_vmax=None,
                               orientation="wave", title=None, figsize=(12, 10)):
    """One run on a 2 x 2 figure.

    Panels: the input raster; the normalised 2D spectrum with significant cells
    outlined; the radial spectrum with rollover; and the significant 2D cells
    on a polar plot. Panels whose outputs are missing show a note instead.

    Args:
        run: a :class:`~lsdviztools3.io.spectral.SpectralAnalysisRun`, or an
            output stem to load with
            :func:`~lsdviztools3.io.spectral.load_spectral_analysis`.
        input_cmap: colormap for the input raster (default suits a
            soil-to-vegetation density).
        input_label: title of the input panel.
        polar_r_max: outer radius of the polar plot in cycles/m. Defaults to
            the corner of the 2D spectrum (sqrt(2) x Nyquist) when the rasters
            are present.
        polar_vmin, polar_vmax: colour limits of the polar plot; fix them to
            compare several tiles on one scale.
        orientation: ``"wave"`` or ``"crest"``, passed to
            :func:`plot_significant_frequencies_polar`.
        title: figure title; defaults to the stem.
        figsize: figure size in inches.

    Returns:
        The matplotlib Figure.
    """
    plt, _ = _import_mpl()
    from ..io.spectral import load_spectral_analysis

    if not hasattr(run, "radial"):
        run = load_spectral_analysis(run)

    fig = plt.figure(figsize=figsize, layout="constrained")
    _plot_input(run, fig.add_subplot(2, 2, 1), input_cmap, input_label)

    ax = fig.add_subplot(2, 2, 2)
    if run.normalised_psd_2d is not None:
        plot_normalised_spectrum_2d(run, ax=ax)
    else:
        ax.text(0.5, 0.5, "no 2D spectrum rasters", transform=ax.transAxes,
                ha="center", va="center", color="0.4")
        ax.set_axis_off()

    ax = fig.add_subplot(2, 2, 3)
    if run.radial is not None and len(run.radial):
        plot_radial_spectrum(run, ax=ax)
    else:
        ax.text(0.5, 0.5, "no radial spectrum\n(write_radial_spectrum: true)",
                transform=ax.transAxes, ha="center", va="center", color="0.4")
        ax.set_axis_off()

    ax = fig.add_subplot(2, 2, 4, projection="polar")
    if run.significant_2d is not None:
        r_max = polar_r_max
        if r_max is None and run.psd_2d is not None:
            from ..io.raster import get_extent
            r_max = np.sqrt(2.0) * max(abs(v) for v in get_extent(run.psd_2d))
        threshold = run.significant_2d.threshold
        plot_significant_frequencies_polar(
            run.significant_2d, ax=ax, r_max=r_max, orientation=orientation,
            vmin=polar_vmin if polar_vmin is not None else (
                threshold if np.isfinite(threshold) else None),
            vmax=polar_vmax,
            title=f"significant 2D cells: {run.significant_2d.n_cells}")
    else:
        ax.set_axis_off()
        ax.text(0.5, 0.5, "no _significant_frequencies_2d.csv",
                transform=ax.transAxes, ha="center", va="center", color="0.4")

    fig.suptitle(title if title is not None else run.stem, fontsize=11)
    return fig
