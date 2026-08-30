"""Residual distributions and statistics for hillslope model comparisons."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


def residual_statistics(residual, *, mask=None) -> dict:
    """Return finite residual count, mean, RMS, Linf and requested percentiles."""
    values = np.asarray(residual, dtype=float)
    if mask is not None:
        values = values[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"n": 0, "mean": np.nan, "rms": np.nan, "linf": np.nan,
                "p05": np.nan, "median": np.nan, "p95": np.nan}
    return {"n": int(values.size), "mean": float(values.mean()),
            "rms": float(np.sqrt(np.mean(values ** 2))),
            "linf": float(np.max(np.abs(values))),
            "p05": float(np.percentile(values, 5)),
            "median": float(np.median(values)),
            "p95": float(np.percentile(values, 95))}


def summarise_residuals(residuals, *, hillslope_mask=None, labels=None) -> pd.DataFrame:
    """Summarise one or more named residual arrays, including hillslope-only rows."""
    if isinstance(residuals, np.ndarray):
        residuals = {"residual": residuals}
    rows = []
    for name, values in residuals.items():
        for scope, mask in (("full", None), ("hillslope", hillslope_mask)):
            row = {"comparison": name, "scope": scope}
            row.update(residual_statistics(values, mask=mask))
            rows.append(row)
    return pd.DataFrame(rows)


def plot_residual_density(residuals, *, hillslope_mask=None, bins=40,
                          fixed_xlim=None, figsize=(8, 5), title=None):
    """Plot residual density histograms, optionally using a hillslope mask."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Residual plots need matplotlib.") from exc
    if isinstance(residuals, np.ndarray):
        residuals = {"residual": residuals}
    fig, ax = plt.subplots(figsize=figsize)
    for name, values in residuals.items():
        values = np.asarray(values, float)
        if hillslope_mask is not None:
            values = values[np.asarray(hillslope_mask, bool)]
        values = values[np.isfinite(values)]
        if values.size:
            ax.hist(values, bins=bins, density=True, histtype="step", lw=1.8, label=name)
    if fixed_xlim is not None:
        ax.set_xlim(*fixed_xlim)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Elevation residual (m)")
    ax.set_ylabel("Density")
    ax.grid(alpha=0.25)
    ax.legend()
    if title:
        ax.set_title(title)
    fig.tight_layout()
    return fig, ax


def plot_residual_raster(residual, *, ax=None, vlim=None, cmap="RdBu_r",
                         title="Elevation residual (m)"):
    """Display an in-memory or rasterio-readable residual raster around zero."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Residual plots need matplotlib.") from exc
    if isinstance(residual, (str, bytes, Path)):
        try:
            import rasterio
        except ImportError as exc:  # pragma: no cover
            raise ImportError("Reading residual rasters needs rasterio.") from exc
        with rasterio.open(residual) as src:
            residual = src.read(1, masked=True).filled(np.nan)
    values = np.asarray(residual, float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError("residual raster contains no finite values")
    limit = float(np.max(np.abs(finite))) if vlim is None else float(vlim)
    limit = max(limit, np.finfo(float).eps)
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(values, cmap=cmap, vmin=-limit, vmax=limit)
    ax.set_title(title)
    ax.figure.colorbar(image, ax=ax, label="m")
    return ax.figure, ax


plot_hillslope_residual_density = plot_residual_density
plot_hillslope_residuals = plot_residual_density
summarize_residuals = summarise_residuals
