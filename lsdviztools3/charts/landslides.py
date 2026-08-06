"""Landslide/cosmo comparison charts (matplotlib).

Helpers for:
- realised landslide size (area-volume) against Larsen-style power laws
- sample-effect apparent-rate comparisons from ``*_LS_SAMPLE_EFFECT.csv``
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import pandas as pd

_MPL_MSG = (
    "Landslide plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def load_landslide_inventory(path: str) -> pd.DataFrame:
    """Load a ``*_LS_INVENTORY.csv`` table."""
    return pd.read_csv(path)


def plot_landslide_area_volume_larsen(
    inventory: pd.DataFrame,
    *,
    alpha: float = 0.2344,
    gamma: float = 1.410,
    title: str = "Realised landslide size vs Larsen power law",
    ax=None,
):
    """Scatter realised area-volume pairs with a Larsen-style reference curve.

    Args:
        inventory: DataFrame with ``area_m2`` and ``volume_m3`` columns.
        alpha: Prefactor in ``V = alpha * A**gamma``.
        gamma: Exponent in ``V = alpha * A**gamma``.
        title: Figure title.
        ax: Existing matplotlib Axes or ``None``.
    """
    plt = _import_mpl()
    req = {"area_m2", "volume_m3"}
    missing = req - set(inventory.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    a = pd.to_numeric(inventory["area_m2"], errors="coerce").to_numpy(dtype=float)
    v = pd.to_numeric(inventory["volume_m3"], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(a) & np.isfinite(v) & (a > 0.0) & (v > 0.0)
    a = a[keep]
    v = v[keep]
    if a.size == 0:
        raise ValueError("No positive finite area/volume pairs to plot.")

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 4.8))
    else:
        fig = ax.figure

    ax.scatter(a, v, s=18, alpha=0.7, label="Realised landslides")
    a_line = np.logspace(np.log10(float(a.min())), np.log10(float(a.max())), 200)
    v_line = alpha * np.power(a_line, gamma)
    ax.plot(a_line, v_line, linewidth=2.0,
            label=f"Reference: V={alpha:g} A^{gamma:g}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Area (m^2)")
    ax.set_ylabel("Volume (m^3)")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle=":", alpha=0.35)
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig, ax


def _finite_xy(df: pd.DataFrame, x_col: str, y_col: str) -> Tuple[np.ndarray, np.ndarray]:
    x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[y_col], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    return x[keep], y[keep]


def plot_sample_effect_apparent_vs_true(
    sample_effect: pd.DataFrame,
    *,
    control_col: str = "erate_control_mm_kyr",
    landslide_col: str = "erate_median_mm_kyr",
    true_rate_mm_kyr: Optional[float] = None,
    title: str = "Apparent erosion rates: control vs landslide",
    ax=None,
):
    """Plot apparent control rates against landslide-apparent rates."""
    plt = _import_mpl()
    req = {control_col, landslide_col}
    missing = req - set(sample_effect.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    x, y = _finite_xy(sample_effect, control_col, landslide_col)
    if x.size == 0:
        raise ValueError("No finite apparent-rate pairs to plot.")

    if ax is None:
        fig, ax = plt.subplots(figsize=(6.8, 4.8))
    else:
        fig = ax.figure

    ax.scatter(x, y, s=32, alpha=0.85)
    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    pad = 0.06 * (hi - lo if hi > lo else 1.0)
    lims = (lo - pad, hi + pad)
    ax.plot(lims, lims, linestyle="--", linewidth=1.1, color="0.25", label="1:1")
    if true_rate_mm_kyr is not None and np.isfinite(true_rate_mm_kyr):
        ax.axhline(float(true_rate_mm_kyr), linewidth=1.4, color="tab:red",
                   label=f"True rate: {true_rate_mm_kyr:g} mm/kyr")
    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_xlabel("Apparent control erosion rate (mm/kyr)")
    ax.set_ylabel("Apparent landslide erosion rate (mm/kyr)")
    ax.set_title(title)
    ax.grid(True, linestyle=":", alpha=0.35)
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig, ax


def plot_sample_effect_ratio(
    sample_effect: pd.DataFrame,
    *,
    ratio_col: str = "ratio_to_control_median",
    p05_col: str = "ratio_to_control_p05",
    p95_col: str = "ratio_to_control_p95",
    sample_name_col: str = "sample_name",
    title: str = "Landslide bias ratio by sample",
    ax=None,
):
    """Plot sample-wise landslide/control apparent-rate ratio with uncertainty."""
    plt = _import_mpl()
    req = {ratio_col}
    missing = req - set(sample_effect.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    ratio = pd.to_numeric(sample_effect[ratio_col], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(ratio)
    if not keep.any():
        raise ValueError("No finite ratio values to plot.")

    names = sample_effect[sample_name_col].astype(str).to_numpy() if sample_name_col in sample_effect.columns else None
    ratio = ratio[keep]
    if names is not None:
        names = names[keep]

    x = np.arange(ratio.size)

    has_ci = p05_col in sample_effect.columns and p95_col in sample_effect.columns
    if has_ci:
        p05 = pd.to_numeric(sample_effect[p05_col], errors="coerce").to_numpy(dtype=float)[keep]
        p95 = pd.to_numeric(sample_effect[p95_col], errors="coerce").to_numpy(dtype=float)[keep]
        ok = np.isfinite(p05) & np.isfinite(p95)
    else:
        p05 = p95 = None
        ok = np.zeros_like(ratio, dtype=bool)

    if ax is None:
        fig, ax = plt.subplots(figsize=(max(6.5, 0.75 * ratio.size), 4.2))
    else:
        fig = ax.figure

    if has_ci and ok.any():
        yerr_lo = np.where(ok, ratio - p05, 0.0)
        yerr_hi = np.where(ok, p95 - ratio, 0.0)
        ax.errorbar(x, ratio, yerr=[yerr_lo, yerr_hi], fmt="o", capsize=3)
    else:
        ax.plot(x, ratio, "o")

    ax.axhline(1.0, linestyle="--", linewidth=1.1, color="0.25")
    ax.set_ylabel("Apparent rate ratio (landslide/control)")
    ax.set_title(title)
    if names is not None:
        ax.set_xticks(x)
        ax.set_xticklabels(names)
    ax.grid(True, axis="y", linestyle=":", alpha=0.35)
    fig.tight_layout()
    return fig, ax
