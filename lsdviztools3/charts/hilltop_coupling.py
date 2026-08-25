"""Charts for hilltop / channel coupling outputs.

Plan section 4.7 groups A and B: the diagnostics that answer *does the tracing
work?* and *is the ridge set right?*.  These are built during Phase 0/1 rather
than saved for the end, because looking at the traces over a hillshade is the
fastest way to find out that a trace kernel does something wrong, and it beats
reading the code.

Every quantity drawn here is read from a file lsdtt3 wrote.  Nothing about the
landscape is computed in this module — the one exception is the hillshade, which
is presentation only and is stated as such where it is used (the same latitude
``charts/hand.py`` takes for its reach maps).

Each function returns a ``matplotlib.figure.Figure``.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from ..io.hilltop_coupling import (
    OUTCOME_COLOUR,
    OUTCOME_NAME,
    SIDE_COLOUR,
    SIDE_NAME,
    HilltopCouplingOutputs,
)

_MPL_MSG = (
    "The hilltop-coupling charts need matplotlib. "
    "Install it e.g. `pip install matplotlib`."
)
_RASTERIO_MSG = (
    "The hilltop-coupling maps need rasterio to read the DEM and channel "
    "rasters. Install it e.g. `pip install rasterio`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def _import_rasterio():
    try:
        import rasterio  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_RASTERIO_MSG) from exc
    return rasterio


def _hillshade(z: np.ndarray, cell: float, azimuth: float = 315.0,
               altitude: float = 45.0) -> np.ndarray:
    """Horn hillshade of an elevation grid.

    Presentation only: nothing quantitative is read off it, so the choice of
    illumination is a styling decision rather than an analysis one.
    """
    dy, dx = np.gradient(z, cell)
    slope = np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(-dx, dy)
    az = np.radians(360.0 - azimuth + 90.0)
    alt = np.radians(altitude)
    shade = (np.sin(alt) * np.cos(slope)
             + np.cos(alt) * np.sin(slope) * np.cos(az - aspect))
    return np.clip(shade, 0.0, 1.0)


def _read_grid(path):
    rasterio = _import_rasterio()
    with rasterio.open(path) as src:
        return src.read(1), src.transform


def _grid_extent(n_rows: int, n_cols: int, cell: float):
    """Plot extent in metres with row 0 at the top."""
    return [0.0, n_cols * cell, n_rows * cell, 0.0]


def _rc_to_metres(row, col, cell):
    return (np.asarray(col) + 0.5) * cell, (np.asarray(row) + 0.5) * cell


def _draw_base(ax, dem, channels, cell, extent, channel_alpha=0.85):
    import matplotlib

    ax.imshow(_hillshade(dem, cell), cmap="Greys_r", extent=extent,
              vmin=0, vmax=1, interpolation="bilinear")
    if channels is not None:
        masked = np.ma.masked_where(channels < 0.5, channels)
        ax.imshow(masked,
                  cmap=matplotlib.colors.ListedColormap(["#1f6fb4"]),
                  extent=extent, alpha=channel_alpha, interpolation="nearest")
    ax.set_xlabel("easting (m, grid-relative)")
    ax.set_ylabel("northing (m, grid-relative)")


# --------------------------------------------------------------------------
# Group A: does the tracing work?
# --------------------------------------------------------------------------

def plot_traces_on_hillshade(outputs: HilltopCouplingOutputs, *,
                             cell: Optional[float] = None,
                             figsize=(9, 9)):
    """Traces from ridge pixels to the channel network, over a hillshade.

    The single most valuable diagnostic in the set: traces that cross a divide,
    run uphill, or die mid-slope are obvious here and invisible in any summary
    statistic.

    Args:
        outputs: A loaded output set; needs ``paths`` and the DEM raster.
        cell: Cell size in metres. Taken from ``run_summary`` when omitted.
        figsize: Figure size in inches.

    Returns:
        The figure.
    """
    plt = _import_mpl()

    cell = float(cell if cell is not None else outputs.summary.get("cell_size", 1.0))
    dem, transform = _read_grid(outputs.dem_path)
    channels = (_read_grid(outputs.channels_path)[0]
                if outputs.channels_path else None)

    n_rows, n_cols = dem.shape
    extent = _grid_extent(n_rows, n_cols, cell)

    fig, ax = plt.subplots(figsize=figsize)
    _draw_base(ax, dem, channels, cell, extent)

    inverse = ~transform
    n_drawn = 0
    if outputs.paths is not None and len(outputs.paths):
        for trace_id, group in outputs.paths.groupby("trace_id"):
            cols, rows = inverse * (group["x"].to_numpy(),
                                    group["y"].to_numpy())
            outcome = int(group["outcome"].iloc[0])
            ax.plot(np.asarray(cols) * cell, np.asarray(rows) * cell,
                    lw=0.7, alpha=0.85,
                    color=OUTCOME_COLOUR.get(outcome, "#333333"))
            n_drawn += 1

    handles = [plt.Line2D([], [], color="#1f6fb4", lw=3,
                          label="channel network")]
    if outputs.paths is not None and len(outputs.paths):
        for outcome in sorted(outputs.paths["outcome"].unique()):
            handles.append(plt.Line2D(
                [], [], color=OUTCOME_COLOUR[int(outcome)], lw=1.5,
                label=f"trace: {OUTCOME_NAME[int(outcome)]}"))
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.9)

    n_traces = int(outputs.summary.get("n_traces", 0))
    n_reached = int(outputs.summary.get("n_reached_sink", 0))
    pct = 100.0 * n_reached / n_traces if n_traces else 0.0
    ax.set_title(
        "Hilltop flow routing: traces from ridge pixels to the channel network\n"
        f"{n_traces} traces, {pct:.1f} % reach a channel; {n_drawn} drawn",
        fontsize=10)

    fig.tight_layout()
    return fig


def plot_trace_outcomes(outputs: HilltopCouplingOutputs, *,
                        cell: Optional[float] = None,
                        figsize=(16, 4.8)):
    """Outcome counts, where failures cluster, and the Lh / relief distributions.

    The border-drop count is written into the suptitle rather than left in a log:
    hilltop pixels on the DEM edge produce no trace at all, so they appear in
    none of these panels and the totals will not otherwise reconcile.

    Args:
        outputs: A loaded output set; needs ``traces``.
        cell: Cell size in metres. Taken from ``run_summary`` when omitted.
        figsize: Figure size in inches.

    Returns:
        The figure.
    """
    plt = _import_mpl()

    cell = float(cell if cell is not None else outputs.summary.get("cell_size", 1.0))
    traces = outputs.traces
    if traces is None or not len(traces):
        raise ValueError("plot_trace_outcomes: the output set carries no traces")

    fig, axes = plt.subplots(1, 3, figsize=figsize)

    counts = outputs.outcome_counts()
    axes[0].bar([OUTCOME_NAME[int(k)] for k in counts.index], counts.to_numpy(),
                color=[OUTCOME_COLOUR[int(k)] for k in counts.index])
    axes[0].set_ylabel("traces")
    axes[0].set_title("Trace outcomes")
    axes[0].tick_params(axis="x", rotation=20)
    for i, value in enumerate(counts.to_numpy()):
        axes[0].text(i, value, str(int(value)), ha="center", va="bottom",
                     fontsize=8)

    for outcome, group in traces.groupby("outcome"):
        x, y = _rc_to_metres(group["source_row"], group["source_col"], cell)
        axes[1].scatter(x, y, s=8, c=OUTCOME_COLOUR[int(outcome)],
                        label=OUTCOME_NAME[int(outcome)])
    axes[1].invert_yaxis()
    axes[1].set_aspect("equal")
    axes[1].legend(fontsize=7)
    axes[1].set_title("Where failures cluster")
    axes[1].set_xlabel("easting (m)")
    axes[1].set_ylabel("northing (m)")

    good = outputs.successful_traces()
    axes[2].hist(good["length_m"], bins=30, color="#2a7f3f", alpha=0.75,
                 label=f"Lh (median {good['length_m'].median():.0f} m)")
    twin = axes[2].twinx()
    twin.hist(good["relief_m"], bins=30, color="#8c564b", alpha=0.45,
              label=f"relief (median {good['relief_m'].median():.1f} m)")
    axes[2].set_xlabel("metres")
    axes[2].set_ylabel("count: hillslope length")
    twin.set_ylabel("count: relief")
    h1, l1 = axes[2].get_legend_handles_labels()
    h2, l2 = twin.get_legend_handles_labels()
    axes[2].legend(h1 + h2, l1 + l2, fontsize=8)
    axes[2].set_title("Hillslope length and relief (successful traces only)")

    dropped = outputs.border_dropped()
    n_hilltops = int(outputs.summary.get("n_hilltop_pixels", 0))
    if dropped and n_hilltops:
        fig.suptitle(
            f"Trace diagnostics. NOTE: {dropped} of {n_hilltops} hilltop pixels "
            f"({100.0 * dropped / n_hilltops:.0f} %) were dropped at the DEM "
            "border before tracing and appear nowhere here",
            fontsize=10, y=1.02)

    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# Group B: is the ridge set right?
# --------------------------------------------------------------------------

def plot_ridgelines_and_sides(outputs: HilltopCouplingOutputs, *,
                              cell: Optional[float] = None,
                              figsize=(16, 8)):
    """The ordered divides, and the hilltop pixels they separate.

    This is where the divide-grid indexing is checked by eye: the hilltop pixels
    must straddle each ridgeline, one colour on each flank. A panel where one
    colour dominates, or where the two interleave, means the side assignment has
    collapsed.

    Args:
        outputs: A loaded output set; needs ``ridgelines``, ``hilltops`` and the
            DEM raster.
        cell: Cell size in metres. Taken from ``run_summary`` when omitted.
        figsize: Figure size in inches.

    Returns:
        The figure.
    """
    plt = _import_mpl()

    cell = float(cell if cell is not None else outputs.summary.get("cell_size", 1.0))
    dem, transform = _read_grid(outputs.dem_path)
    channels = (_read_grid(outputs.channels_path)[0]
                if outputs.channels_path else None)

    n_rows, n_cols = dem.shape
    extent = _grid_extent(n_rows, n_cols, cell)

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    for ax in axes:
        _draw_base(ax, dem, channels, cell, extent, channel_alpha=0.7)

    inverse = ~transform
    cmap = plt.get_cmap("tab20")

    n_ridgelines = 0
    if outputs.ridgelines is not None and len(outputs.ridgelines):
        for rid, group in outputs.ridgelines.groupby("ridgeline_id"):
            ordered = group.sort_values("node_index")
            cols, rows = inverse * (ordered["x"].to_numpy(),
                                    ordered["y"].to_numpy())
            axes[0].plot(np.asarray(cols) * cell, np.asarray(rows) * cell,
                         lw=1.6, color=cmap(int(rid) % 20))
            n_ridgelines += 1

    axes[0].set_title(
        "Ordered ridgelines (the divides themselves)\n"
        f"{n_ridgelines} ridgelines on the pixel-corner divide grid, "
        "coloured by ridgeline_id", fontsize=10)

    n_hilltops = 0
    if outputs.hilltops is not None and len(outputs.hilltops):
        n_hilltops = len(outputs.hilltops)
        for side, group in outputs.hilltops.groupby("ridge_side"):
            x, y = _rc_to_metres(group["row"], group["col"], cell)
            axes[1].scatter(x, y, s=6, c=SIDE_COLOUR[int(side)],
                            label=f"side {int(side)} ({SIDE_NAME[int(side)]})")
        axes[1].legend(loc="upper right", fontsize=8, framealpha=0.9)

    axes[1].set_title(
        "Hilltop pixels: the pixels the ridgelines separate\n"
        f"{n_hilltops} pixels; the two colours are the opposing hillslopes "
        "of each divide", fontsize=10)

    fig.tight_layout()
    return fig


def plot_ridgeline_profiles(outputs: HilltopCouplingOutputs, *,
                            max_ridgelines: int = 8,
                            figsize=(10, 6)):
    """Hilltop-pixel count against distance along each divide.

    Uses the ordering directly, which is what makes the ridgeline representation
    worth having: a divide is a line with a coordinate along it, not a bag of
    pixels. Gaps here are places where the filters removed a stretch of divide.

    Args:
        outputs: A loaded output set; needs ``hilltops``.
        max_ridgelines: Draw at most this many, longest first.
        figsize: Figure size in inches.

    Returns:
        The figure.
    """
    plt = _import_mpl()

    hilltops = outputs.hilltops
    if hilltops is None or not len(hilltops):
        raise ValueError(
            "plot_ridgeline_profiles: the output set carries no hilltop pixels")

    fig, ax = plt.subplots(figsize=figsize)

    lengths = (hilltops.groupby("ridgeline_id")["ridgeline_total_length_m"]
               .max().sort_values(ascending=False))

    cmap = plt.get_cmap("tab10")
    for i, rid in enumerate(lengths.index[:max_ridgelines]):
        group = hilltops[hilltops["ridgeline_id"] == rid]
        for side, side_group in group.groupby("ridge_side"):
            distances = np.sort(side_group["distance_along_ridgeline_m"].to_numpy())
            ax.plot(distances, np.full_like(distances, i, dtype=float)
                    + (0.15 if int(side) else -0.15),
                    marker="|", ls="none", ms=6,
                    color=SIDE_COLOUR[int(side)])
        ax.axhline(i, color=cmap(i % 10), lw=0.4, alpha=0.4)

    ax.set_yticks(range(min(max_ridgelines, len(lengths))))
    ax.set_yticklabels([f"ridgeline {int(r)}"
                        for r in lengths.index[:max_ridgelines]], fontsize=8)
    ax.set_xlabel("distance along ridgeline (m)")
    ax.set_title(
        "Hilltop pixels along each divide, by side\n"
        "upper ticks = right of travel, lower = left; gaps are filtered stretches",
        fontsize=10)

    fig.tight_layout()
    return fig
