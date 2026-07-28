"""Cosmogenic comparison charts (matplotlib).

Tools for plotting lsdtt3-vs-reference comparison CSVs from the lsdtt3 cosmo
test/workflow scripts, including single-point sweep outputs and real-sample
comparison tables.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

_MPL_MSG = (
    "Cosmo comparison plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def load_cosmo_table(path: str | Path) -> pd.DataFrame:
    """Load a cosmo comparison CSV, allowing leading `#` comment lines."""
    return pd.read_csv(path, comment="#")


def _residual_col(reference: str) -> str:
    ref = reference.strip().lower()
    if ref in {"cronus_v3", "cronus", "v3"}:
        return "pct_diff_vs_cronus_v3"
    if ref in {"cronuscalc", "cronus_calc"}:
        return "pct_diff_vs_cronuscalc"
    raise ValueError("reference must be one of: cronus_v3, cronuscalc")


def _nearest_filter(df: pd.DataFrame, col: str, target: Optional[float]) -> pd.DataFrame:
    if target is None:
        return df
    vals = pd.to_numeric(df[col], errors="coerce")
    if vals.isna().all():
        return df
    nearest = vals.iloc[(vals - target).abs().argmin()]
    return df[np.isclose(vals, nearest)]


def plot_cosmo_sweep_panels(
    df: pd.DataFrame,
    *,
    reference: str = "cronus_v3",
    representative_pressure_hpa: Optional[float] = None,
    representative_erate_gcm2yr: Optional[float] = None,
    aggregate: str = "median",
):
    """Plot sweep residuals vs latitude in nuclide x scaling panels.

    Expects a table like ``single_point_scan_comparison.csv`` with columns:
    ``nuclide, scaling_scheme, muon_scheme, lat_deg`` and
    ``pct_diff_vs_cronus_v3`` / ``pct_diff_vs_cronuscalc``.
    """
    plt = _import_mpl()
    resid_col = _residual_col(reference)

    req = {"nuclide", "scaling_scheme", "muon_scheme", "lat_deg", resid_col}
    missing = req - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns for sweep panels: {sorted(missing)}")

    work = df.copy()
    work = _nearest_filter(work, "pressure_hpa", representative_pressure_hpa)
    work = _nearest_filter(work, "erate_gcm2yr", representative_erate_gcm2yr)
    work[resid_col] = pd.to_numeric(work[resid_col], errors="coerce")
    work["lat_deg"] = pd.to_numeric(work["lat_deg"], errors="coerce")
    work = work[np.isfinite(work[resid_col]) & np.isfinite(work["lat_deg"])]
    if work.empty:
        raise ValueError("No plottable rows after filtering and numeric conversion.")

    if aggregate not in {"median", "mean"}:
        raise ValueError("aggregate must be 'median' or 'mean'")
    agg_name = aggregate

    nuclides = [n for n in ("Be10", "Al26", "C14") if n in set(work["nuclide"])]
    scalings = [s for s in ("stone", "lsdn") if s in set(work["scaling_scheme"])]
    if not nuclides or not scalings:
        raise ValueError("No nuclide/scaling rows available for panel plot.")

    fig, axes = plt.subplots(len(nuclides), len(scalings), figsize=(5.2 * len(scalings), 3.4 * len(nuclides)), sharex=True, sharey=True)
    if len(nuclides) == 1 and len(scalings) == 1:
        axes = np.array([[axes]])
    elif len(nuclides) == 1:
        axes = np.array([axes])
    elif len(scalings) == 1:
        axes = np.array([[ax] for ax in axes])

    for i, nuc in enumerate(nuclides):
        for j, scaling in enumerate(scalings):
            ax = axes[i, j]
            sub = work[(work["nuclide"] == nuc) & (work["scaling_scheme"] == scaling)]
            if sub.empty:
                ax.set_visible(False)
                continue
            for muon in sorted(sub["muon_scheme"].dropna().unique()):
                s2 = sub[sub["muon_scheme"] == muon]
                grp = s2.groupby("lat_deg", as_index=False)[resid_col].agg(agg_name).sort_values("lat_deg")
                ax.plot(grp["lat_deg"], grp[resid_col], marker="o", markersize=3.5, linewidth=1.2, label=str(muon))

            ax.axhline(0.0, color="0.5", linewidth=0.8, linestyle="--")
            ax.grid(True, alpha=0.25)
            ax.set_title(f"{nuc} | {scaling}")
            if i == len(nuclides) - 1:
                ax.set_xlabel("Latitude (deg)")
            if j == 0:
                ax.set_ylabel(f"Residual vs {reference} (%)")
            if i == 0 and j == len(scalings) - 1:
                ax.legend(title="Muon", fontsize=8)

    fig.suptitle("lsdtt3 comparison sweep: residual by latitude", y=0.995)
    fig.tight_layout()
    return fig, axes


def _first_existing(cols: Iterable[str], all_cols: Iterable[str]) -> Optional[str]:
    cset = set(all_cols)
    for c in cols:
        if c in cset:
            return c
    return None


def plot_cosmo_sample_residuals(
    df: pd.DataFrame,
    *,
    residual_col: Optional[str] = None,
    group_col: str = "dataset",
    split_col: str = "scaling_scheme",
):
    """Plot grouped sample residuals (boxplots) for real-sample comparison tables."""
    plt = _import_mpl()

    if residual_col is None:
        residual_col = _first_existing(
            ["pct_diff_lsdtt_minus_reference", "pct_diff_lsdtt_minus_octopus", "pct_diff_lsdtt_minus_published"],
            df.columns,
        )
        if residual_col is None:
            raise KeyError("Could not infer residual column; pass residual_col explicitly.")

    req = {residual_col, group_col}
    missing = req - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns for sample residual plot: {sorted(missing)}")

    work = df.copy()
    work[residual_col] = pd.to_numeric(work[residual_col], errors="coerce")
    work = work[np.isfinite(work[residual_col])]
    if work.empty:
        raise ValueError("No finite residual values to plot.")

    groups = list(pd.unique(work[group_col]))
    if not groups:
        raise ValueError("No groups available for boxplot.")

    has_split = split_col in work.columns and work[split_col].notna().any()
    fig, ax = plt.subplots(figsize=(max(7.0, 1.8 * len(groups)), 4.6))

    if has_split:
        splits = list(pd.unique(work[split_col]))
        width = 0.8 / max(1, len(splits))
        handles = []
        for i, sp in enumerate(splits):
            vals = []
            positions = []
            for g_idx, g in enumerate(groups):
                v = work[(work[group_col] == g) & (work[split_col] == sp)][residual_col].to_numpy()
                if v.size:
                    vals.append(v)
                    positions.append(g_idx + (i - (len(splits) - 1) / 2.0) * width)
            if not vals:
                continue
            bp = ax.boxplot(vals, positions=positions, widths=0.9 * width, patch_artist=True, manage_ticks=False)
            color = f"C{i % 10}"
            for box in bp["boxes"]:
                box.set(facecolor=color, alpha=0.5)
            for median in bp["medians"]:
                median.set(color="black", linewidth=1.1)
            handles.append((color, str(sp)))
        for color, label in handles:
            ax.plot([], [], color=color, linewidth=6, alpha=0.5, label=label)
        ax.legend(title=split_col, fontsize=8)
    else:
        vals = [work[work[group_col] == g][residual_col].to_numpy() for g in groups]
        ax.boxplot(vals, patch_artist=True)

    ax.axhline(0.0, color="0.5", linewidth=0.8, linestyle="--")
    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups, rotation=20, ha="right")
    ax.set_ylabel("Residual (%)")
    ax.set_title("lsdtt3 comparison residuals by group")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    return fig, ax
