"""Plots for the common-time transient residual CSV."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = ("time_yr", "mean", "rms", "linf", "p05", "median", "p95",
                    "hillslope_mean", "hillslope_rms")


def load_transient_residuals(path: str | Path) -> pd.DataFrame:
    """Load and validate a transient residual statistics CSV."""
    frame = pd.read_csv(path)
    missing = sorted(set(REQUIRED_COLUMNS) - set(frame.columns))
    if missing:
        raise KeyError(f"Missing transient residual columns: {missing}")
    frame = frame.sort_values("time_yr", kind="stable").reset_index(drop=True)
    if frame["time_yr"].duplicated().any():
        raise ValueError("transient residual CSV contains duplicate times")
    return frame


def plot_transient_residuals(data, *, figsize=(8, 4.5), title=None):
    """Plot RMS and maximum absolute residuals through time."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Transient residual plots need matplotlib.") from exc
    frame = load_transient_residuals(data) if isinstance(data, (str, Path)) else data.copy()
    missing = sorted(set(REQUIRED_COLUMNS) - set(frame.columns))
    if missing:
        raise KeyError(f"Missing transient residual columns: {missing}")
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(frame.time_yr, frame.rms, label="full RMS")
    ax.plot(frame.time_yr, frame.linf, label="maximum absolute residual")
    ax.plot(frame.time_yr, frame.hillslope_rms, "--", label="hillslope RMS")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Elevation residual (m)")
    ax.grid(alpha=0.25)
    ax.legend()
    if title:
        ax.set_title(title)
    fig.tight_layout()
    return fig, ax


def plot_strip_transient_metrics(data, *, reference_time=5000.0,
                                 figsize=(10, 9), title=None):
    """Plot changes in relief, maximum gradient and divide Laplacian by model."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Transient metric plots need matplotlib.") from exc
    frame = pd.read_csv(data) if isinstance(data, (str, Path)) else data.copy()
    required = {"time_yr", "model", "relief_m", "max_abs_gradient",
                "divide_laplacian_m-1"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise KeyError(f"Missing strip transient metric columns: {missing}")
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
    metrics = (("relief_m", "Change in relief (m)", False),
               ("max_abs_gradient", "Change in maximum |gradient|", False),
               ("divide_laplacian_m-1",
                "Change in |divide Laplacian| (m$^{-1}$)", True))
    for model, subset in frame.groupby("model", sort=False):
        subset = subset.sort_values("time_yr")
        reference_index = (subset.time_yr - reference_time).abs().idxmin()
        for ax, (column, ylabel, magnitude) in zip(axes, metrics):
            values = subset[column].abs() if magnitude else subset[column]
            reference = abs(subset.loc[reference_index, column]) if magnitude else \
                subset.loc[reference_index, column]
            ax.plot(subset.time_yr, values - reference,
                    label=model)
            ax.set_ylabel(ylabel); ax.grid(alpha=0.25)
    for ax in axes:
        ax.axvline(reference_time, color="0.3", ls=":", lw=1)
    axes[0].legend()
    axes[-1].set_xlabel("Time (yr)")
    fig.suptitle(title or "Strip response after uplift doubles")
    fig.tight_layout()
    return fig, axes


plot_hillslope_transient_residuals = plot_transient_residuals
