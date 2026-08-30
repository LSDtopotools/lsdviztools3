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
    """Plot full-domain and hillslope RMS/Linf residuals through time."""
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
    ax.plot(frame.time_yr, frame.linf, label="full Linf")
    ax.plot(frame.time_yr, frame.hillslope_rms, "--", label="hillslope RMS")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Elevation residual (m)")
    ax.grid(alpha=0.25)
    ax.legend()
    if title:
        ax.set_title(title)
    fig.tight_layout()
    return fig, ax


plot_hillslope_transient_residuals = plot_transient_residuals
