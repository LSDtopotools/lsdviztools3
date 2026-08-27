"""Gradient-curvature charts: does the flux law explain the topography?

After Roering, J. J. (2008). How well can hillslope evolution models 'explain'
topography? Simulating soil transport and production with high-resolution
topographic data. *Geological Society of America Bulletin* 120:1248-1262.
doi:10.1130/B26283.1

**Why this plot decides something.** On a hillslope at steady state, the two
candidate transport laws make *qualitatively* different predictions about how
curvature varies with gradient -- not merely different parameter values. So the
shape of the measured cloud discriminates between them without fitting anything.

Take a one-dimensional hillslope eroding uniformly at rate ``E``. Flux at
distance ``x`` from the divide must be ``q = E x``. Then:

**Linear diffusion**, ``q = D S``:

    E x = D S   =>   E = D dS/dx   =>   C = E / D

Curvature is CONSTANT, independent of gradient -- a horizontal line. A linearly
diffusing hillslope is a parabola, and a parabola has one curvature everywhere.

**Nonlinear (Roering) transport**, ``q = D S / (1 - (S/S_c)^2)``:

    E = D (dS/dx) (1 + (S/S_c)^2) / (1 - (S/S_c)^2)^2

    =>   C = (E/D) * (1 - (S/S_c)^2)^2 / (1 + (S/S_c)^2)

Curvature FALLS with gradient, from ``E/D`` at ``S = 0`` to zero as ``S -> S_c``.
Physically: near the critical gradient the flux diverges, so a vanishing
curvature is enough to carry the sediment, and the hillslope goes planar. That
is the observation the nonlinear law was introduced to explain -- real steep
hillslopes are planar, not parabolic, and linear diffusion cannot produce that
at any ``D``.

**Convexity, not curvature.** ``C`` above is ``-∇²z``, positive on a convex
hilltop. lsdtt3 uses the strict mathematical convention, in which a convex-up
hilltop has NEGATIVE Laplacian, so the sign is flipped on reading. Plotting the
raw signed value instead puts the whole diagnostic upside down and it still
looks like a plausible figure, which is why the flip is done here rather than
left to the caller.

**Use the Laplacian.** The derivation is in terms of ``∇²z``. Differential-
geometry mean curvature agrees only at small gradient, and this diagnostic is
most interesting exactly where gradients are large.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

_MPL_MSG = (
    "The gradient-curvature charts need matplotlib. Install it e.g. "
    "`pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def _read(path):
    try:
        import rasterio
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Reading rasters needs rasterio.") from exc
    with rasterio.open(path) as src:
        band = src.read(1).astype(float)
        nodata = src.nodata
    if nodata is not None:
        band = np.where(band == nodata, np.nan, band)
    return np.where(band <= -9990.0, np.nan, band)


def nonlinear_convexity(gradient, E_over_D: float, Sc: float):
    """Steady-state convexity predicted by the Roering (1999) flux law.

    ``C = (E/D) (1 - (S/Sc)^2)^2 / (1 + (S/Sc)^2)``, which is ``E/D`` at zero
    gradient and zero at ``S = Sc``. Verified against numerical differentiation
    of the flux law itself to 6e-9 relative.

    Args:
        gradient: Gradient ``S`` (dimensionless).
        E_over_D: Erosion rate divided by diffusivity, in m^-1. This sets the
            whole curve's height; it is not a free shape parameter.
        Sc: Critical gradient.

    Returns:
        Predicted convexity ``-∇²z``, NaN at or above ``Sc`` where the steady
        state does not exist.
    """
    s = np.asarray(gradient, dtype=float)
    r = np.where(s < Sc, (s / Sc) ** 2, np.nan)
    return E_over_D * (1.0 - r) ** 2 / (1.0 + r)


def linear_convexity(gradient, E_over_D: float):
    """Steady-state convexity predicted by linear diffusion: constant ``E/D``."""
    return np.full_like(np.asarray(gradient, dtype=float), float(E_over_D))


def plot_gradient_curvature(slope_raster,
                            curvature_raster,
                            *,
                            E_over_D: Optional[float] = None,
                            Sc: Optional[float] = None,
                            hillslope_mask_raster=None,
                            area_raster=None,
                            channel_area_threshold: Optional[float] = None,
                            n_bins: int = 30,
                            min_count: int = 25,
                            max_scatter: int = 40000,
                            convex_only: bool = False,
                            curvature_is_laplacian: bool = True,
                            title: Optional[str] = None,
                            figsize: Tuple[float, float] = (8.0, 6.5),
                            seed: int = 0):
    """Convexity against gradient, with the two flux laws overlaid.

    Args:
        slope_raster: Gradient raster (dimensionless).
        curvature_raster: Curvature raster. With ``curvature_is_laplacian``
            (the default) this is ``∇²z`` and is negated to give convexity.
        E_over_D: Erosion rate over diffusivity, m^-1. When given, both theory
            curves are drawn. This is an INPUT from the model or from field
            constraint, not a fit -- the figure tests the laws, it does not
            calibrate them.
        Sc: Critical gradient, for the nonlinear curve.
        hillslope_mask_raster: Optional mask of cells to include.
        area_raster, channel_area_threshold: Alternative to a mask -- exclude
            channel cells by drainage area. **Channels must be excluded**: the
            derivation is a hillslope sediment balance and says nothing about
            fluvial cells, which sit in concave ground and would fold a second
            population into the cloud.
        n_bins: Gradient bins for the median.
        min_count: Bins with fewer cells are dropped.
        convex_only: Restrict the BINNED MEDIAN to convex cells.

            This is not cherry-picking, but the reasoning has to be stated. The
            derivation is a one-dimensional balance for transport away from a
            divide, and it describes convex and planar ground. Convergent
            hollows accumulate sediment and are concave; the 1-D derivation
            makes no claim about them. Including both puts two populations at
            the same gradient, and a median over a bimodal population sits
            between two modes and represents neither.

            The scatter always shows every cell either way, so the excluded
            population stays visible and the reader can judge the choice.
        curvature_is_laplacian: Set False if the raster is already convexity.
        seed: Sampling seed, so a rerun draws the same figure.

    Returns:
        The matplotlib Figure.
    """
    plt = _import_mpl()

    s = _read(slope_raster).ravel()
    c = _read(curvature_raster).ravel()
    if s.size != c.size:
        raise ValueError("slope and curvature rasters are not on the same grid")

    convexity = -c if curvature_is_laplacian else c.copy()

    keep = np.isfinite(s) & np.isfinite(convexity)

    if hillslope_mask_raster is not None:
        m = _read(hillslope_mask_raster).ravel()
        keep &= m > 0.5
    elif area_raster is not None and channel_area_threshold is not None:
        a = _read(area_raster).ravel()
        keep &= np.isfinite(a) & (a < channel_area_threshold)

    if not keep.any():
        raise ValueError("no cells left after masking")

    s_k, cx_k = s[keep], convexity[keep]

    fig, ax = plt.subplots(figsize=figsize)

    rng = np.random.default_rng(seed)
    idx = np.flatnonzero(keep)
    if idx.size > max_scatter:
        idx = rng.choice(idx, max_scatter, replace=False)
    ax.scatter(s[idx], convexity[idx], s=2, alpha=0.10, color="#8c8c8c",
               linewidths=0, label=f"hillslope cells ({s_k.size:,})")

    # Binned median in gradient. Median again: convexity at fixed gradient is
    # skewed, and a mean would be pulled by the tail.
    if convex_only:
        conv = cx_k > 0.0
        s_k, cx_k = s_k[conv], cx_k[conv]
        if s_k.size == 0:
            raise ValueError("no convex cells left after masking")
    hi = float(np.nanpercentile(s_k, 99.5))
    edges = np.linspace(0.0, hi, n_bins + 1)
    which = np.digitize(s_k, edges) - 1
    xs, med, q25, q75 = [], [], [], []
    for b in range(n_bins):
        sel = which == b
        if int(sel.sum()) < min_count:
            continue
        xs.append(0.5 * (edges[b] + edges[b + 1]))
        med.append(np.median(cx_k[sel]))
        q25.append(np.percentile(cx_k[sel], 25))
        q75.append(np.percentile(cx_k[sel], 75))
    xs, med = np.array(xs), np.array(med)

    if xs.size:
        ax.fill_between(xs, q25, q75, color="#d62728", alpha=0.18,
                        label="binned IQR")
        ax.plot(xs, med, "-o", color="#d62728", ms=4, lw=1.8,
                label="binned median" + (" (convex cells)" if convex_only
                                         else " (all cells)"))

    if E_over_D is not None:
        grid = np.linspace(0.0, hi, 400)
        ax.plot(grid, linear_convexity(grid, E_over_D), "--", color="#1f77b4",
                lw=1.8, label=f"linear diffusion:  $C = E/D$ = {E_over_D:.4g}")
        if Sc is not None:
            ax.plot(grid, nonlinear_convexity(grid, E_over_D, Sc), "-",
                    color="#2a7f3f", lw=2.0,
                    label=f"Roering nonlinear, $S_c$ = {Sc:g}")
            ax.axvline(Sc, color="#2a7f3f", ls=":", lw=1.2)
            ax.annotate("$S_c$", xy=(Sc, 0), xytext=(4, 6),
                        textcoords="offset points", color="#2a7f3f")

    ax.axhline(0.0, color="k", lw=0.8, alpha=0.5)
    ax.set_xlabel("gradient  $S$")
    ax.set_ylabel(r"convexity  $-\nabla^2 z$   (m$^{-1}$)")
    ax.set_xlim(0.0, hi)
    ax.grid(True, alpha=0.25)

    if title is None:
        title = ("Gradient against curvature: linear diffusion is flat,\n"
                 "nonlinear transport falls to zero at $S_c$")
    ax.set_title(title)

    leg = ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    for h in leg.legend_handles:
        try:
            h.set_alpha(1.0)
            h.set_sizes([18])
        except Exception:  # noqa: BLE001
            pass

    fig.tight_layout()
    return fig
