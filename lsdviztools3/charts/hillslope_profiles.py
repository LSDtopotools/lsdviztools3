"""Matplotlib diagnostics for one-dimensional hillslope comparisons.

The functions accept arrays so they are useful for synthetic strips, while
``load_profile_raster`` provides the small rasterio bridge needed for lsdtt3
GeoTIFF/COG outputs.  Numerical derivatives use centred differences in the
same coordinate system as the supplied distance array.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np


def _mpl():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Hillslope plots need matplotlib.") from exc
    return plt


def load_profile_raster(path, *, row: Optional[int] = None):
    """Return ``(distance, elevation)`` for a raster centre row.

    GeoTIFF and COG are deliberately the only raster formats implied here;
    rasterio may additionally read legacy formats supported by the installed
    GDAL, but this function does not invent or write a format.
    """
    try:
        import rasterio
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Reading profile rasters needs rasterio.") from exc
    with rasterio.open(path) as src:
        z = src.read(1, masked=True)
        r = z.shape[0] // 2 if row is None else int(row)
        if not 0 <= r < z.shape[0]:
            raise IndexError("profile row is outside the raster")
        values = np.asarray(z[r].filled(np.nan), dtype=float)
        x = np.asarray([src.xy(r, c)[0] for c in range(z.shape[1])])
    return x, values


def profile_derivatives(distance, elevation):
    """Return ``(gradient, curvature)`` using centred finite differences."""
    x = np.asarray(distance, dtype=float)
    z = np.asarray(elevation, dtype=float)
    if x.ndim != 1 or z.ndim != 1 or x.size != z.size or x.size < 3:
        raise ValueError("distance and elevation must be 1-D arrays of length >= 3")
    if not np.all(np.isfinite(x)) or not np.all(np.diff(x) > 0):
        raise ValueError("distance must be finite and strictly increasing")
    if not np.all(np.isfinite(z)):
        raise ValueError("elevation contains non-finite values")
    gradient = np.gradient(z, x, edge_order=1)
    curvature = np.gradient(gradient, x, edge_order=1)
    return gradient, curvature


def profile_residuals(reference: Mapping[str, Sequence[float]],
                      comparison: Mapping[str, Sequence[float]]):
    """Return comparison-minus-reference arrays for shared profile quantities."""
    required = ("elevation", "gradient", "curvature")
    missing = [key for key in required if key not in reference or key not in comparison]
    if missing:
        raise KeyError(f"Missing profile quantities: {missing}")
    return {key: np.asarray(comparison[key], float) - np.asarray(reference[key], float)
            for key in required}


def plot_hillslope_profiles(
    distance,
    profiles: Mapping[str, Mapping[str, Sequence[float]] | Sequence[float]],
    *,
    analytic: Optional[Mapping[str, Mapping[str, Sequence[float]] | Sequence[float]]] = None,
    channel_indices: Optional[Sequence[int]] = None,
    figsize=(12, 10),
    title: Optional[str] = None,
):
    """Plot elevation, gradient, curvature and residual panels.

    ``profiles`` maps labels such as ``"Newton"`` and ``"cellular"`` to either
    elevation arrays or mappings containing ``elevation``, ``gradient`` and
    ``curvature``.  Analytic quantities are supplied in the same mapping. The
    residual panel contains every pairwise residual against the analytic linear
    and Roering profiles, plus cellular-minus-Newton when both are present.
    """
    plt = _mpl()
    x = np.asarray(distance, float)
    if x.ndim != 1:
        raise ValueError("distance must be one-dimensional")
    prepared = {}
    for label, values in profiles.items():
        if isinstance(values, Mapping):
            item = {k: np.asarray(v, float) for k, v in values.items()}
            if "elevation" not in item:
                raise KeyError(f"Profile {label!r} has no elevation")
            if "gradient" not in item or "curvature" not in item:
                item["gradient"], item["curvature"] = profile_derivatives(x, item["elevation"])
        else:
            item = {"elevation": np.asarray(values, float)}
            item["gradient"], item["curvature"] = profile_derivatives(x, item["elevation"])
        prepared[label] = item

    analytic_prepared = {}
    for label, values in (analytic or {}).items():
        if isinstance(values, Mapping):
            item = {k: np.asarray(v, float) for k, v in values.items()}
            if "gradient" not in item or "curvature" not in item:
                item["gradient"], item["curvature"] = profile_derivatives(
                    x, item["elevation"]
                )
        else:
            item = {"elevation": np.asarray(values, float)}
            item["gradient"], item["curvature"] = profile_derivatives(x, item["elevation"])
        analytic_prepared[label] = item

    fig, axes = plt.subplots(4, 1, figsize=figsize, sharex=True,
                             gridspec_kw={"height_ratios": (2, 1, 1, 1.4)})
    quantities = (("elevation", "Elevation (m)"),
                  ("gradient", "Gradient dz/dx"),
                  ("curvature", "Signed curvature d2z/dx2 (m-1)"))
    for key, ylabel in quantities:
        ax = axes[quantities.index((key, ylabel))]
        for label, item in prepared.items():
            ax.plot(x, item[key], label=label)
        for label, item in analytic_prepared.items():
            ax.plot(x, item[key], "--", lw=1.4, label=f"analytic {label}")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
    if channel_indices is not None:
        for ax in axes[:3]:
            for index in np.asarray(channel_indices, int):
                ax.axvline(x[index], color="0.3", ls=":", lw=0.8)

    residual_ax = axes[3]
    reference_names = []
    for ref_name in ("linear", "roering"):
        if ref_name in analytic_prepared:
            reference_names.append(ref_name)
            for label, item in prepared.items():
                residual_ax.plot(x, item["elevation"] - analytic_prepared[ref_name]["elevation"],
                                 label=f"{label} - {ref_name}")
    if "Newton" in prepared and "cellular" in prepared:
        residual_ax.plot(x, prepared["cellular"]["elevation"] - prepared["Newton"]["elevation"],
                         label="cellular - Newton")
    if not reference_names and not ("Newton" in prepared and "cellular" in prepared):
        residual_ax.text(0.5, 0.5, "No residual references supplied",
                         transform=residual_ax.transAxes, ha="center")
    residual_ax.axhline(0, color="k", lw=0.8)
    residual_ax.set_ylabel("Elevation residual (m)")
    residual_ax.set_xlabel("Distance from divide (m)")
    residual_ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8, ncol=2)
    residual_ax.legend(fontsize=7, ncol=2)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, axes


def plot_profile_residual_panels(distance, residuals, *, title=None, figsize=(12, 10)):
    """Plot one signed elevation-residual panel per named comparison."""
    plt = _mpl()
    x = np.asarray(distance, float)
    items = list(residuals.items())
    if not items:
        raise ValueError("residuals must not be empty")
    fig, axes = plt.subplots(len(items), 1, figsize=figsize, sharex=True, squeeze=False)
    for ax, (label, values) in zip(axes[:, 0], items):
        ax.plot(x, np.asarray(values, float), lw=1.2)
        ax.axhline(0.0, color="0.25", lw=0.8)
        ax.set_ylabel(f"{label}\n(m)")
        ax.grid(alpha=0.25)
    axes[-1, 0].set_xlabel("Distance from divide (m)")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, axes[:, 0]


plot_hillslope_profile = plot_hillslope_profiles
