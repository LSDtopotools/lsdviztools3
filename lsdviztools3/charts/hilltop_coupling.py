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
    has_geometry,
    line_vertices,
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
    # np.gradient returns (d/drow, d/dcol); row increases SOUTHWARD in a
    # north-up raster, which is what makes the aspect sign easy to get wrong.
    # This is GDAL's convention (gdaldem hillshade, Horn): aspect =
    # atan2(dz/drow, -dz/dcol). The previous form, atan2(-dx, dy), was exactly
    # 180 degrees out, so every hillshade came back inverted -- ridges read as
    # valleys. Pinned by test_hillshade_illumination_matches_gdal.
    dy, dx = np.gradient(z, cell)
    slope = np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(dy, -dx)
    az = np.radians(360.0 - azimuth + 90.0)
    alt = np.radians(altitude)
    shade = (np.sin(alt) * np.cos(slope)
             + np.cos(alt) * np.sin(slope) * np.cos(az - aspect))
    return np.clip(shade, 0.0, 1.0)


def _load_channels(path):
    """Load a channel network from whatever it was written as.

    Vector (.fgb/.gpkg/.shp) is returned as a GeoDataFrame and drawn from its
    geometry; a raster is returned as an array and drawn as a mask. Returns
    None when no path was supplied.
    """
    if path is None:
        return None
    suffix = str(path).lower()
    if suffix.endswith((".fgb", ".gpkg", ".shp", ".geojson")):
        import geopandas as gpd
        return gpd.read_file(path)
    return _read_grid(path)[0]


def _read_grid(path):
    rasterio = _import_rasterio()
    with rasterio.open(path) as src:
        return src.read(1), src.transform


def _grid_extent(n_rows: int, n_cols: int, cell: float):
    """Plot extent in metres with row 0 at the top."""
    return [0.0, n_cols * cell, n_rows * cell, 0.0]


def _read_grid_extent(path, cell: Optional[float]):
    """Read a raster and return (array, transform, projected extent, cell size).

    Plotting in the raster's own projected CRS means vector features drawn from
    geometry need no transform, which is what lets the same chart serve both
    output shapes.
    """
    rasterio = _import_rasterio()
    with rasterio.open(path) as src:
        z = src.read(1)
        transform = src.transform
        b = src.bounds
        if cell is None:
            cell = float(abs(transform.a))
    return z, transform, [b.left, b.right, b.bottom, b.top], cell


def _hilltop_xy(hilltops, transform, cell):
    """Projected coordinates of hilltop pixels.

    Uses the geometry when the points came from a vector file, and falls back
    to row/col through the raster transform otherwise.
    """
    if has_geometry(hilltops):
        return (hilltops.geometry.x.to_numpy(),
                hilltops.geometry.y.to_numpy())
    cols = hilltops["col"].to_numpy() + 0.5
    rows = hilltops["row"].to_numpy() + 0.5
    xs, ys = rasterio_xy(transform, rows, cols)
    return xs, ys


def rasterio_xy(transform, rows, cols):
    """Vectorised (row, col) -> (x, y) through an affine transform."""
    xs = transform.c + transform.a * cols + transform.b * rows
    ys = transform.f + transform.d * cols + transform.e * rows
    return xs, ys


def _rc_to_metres(row, col, cell):
    return (np.asarray(col) + 0.5) * cell, (np.asarray(row) + 0.5) * cell


def _draw_channels(ax, channels, extent, alpha=0.85):
    """Draw the channel network, from vector geometry or from a raster mask.

    The VECTOR path is preferred and is tried first. A channel network is
    already a set of lines; rasterising it to draw it throws away the geometry,
    ties the drawing to the analysis grid, and needs an extra file that the
    extraction step did not have to write.
    """
    import matplotlib

    if channels is None:
        return

    # Vector: a GeoDataFrame of LineStrings (channel_network_lines) or Points
    # (channel_network_points).
    if has_geometry(channels):
        geom_types = set(channels.geometry.geom_type.unique())
        if geom_types & {"LineString", "MultiLineString"}:
            for _, xs, ys, _ in line_vertices(channels, channels.columns[0]):
                ax.plot(xs, ys, color="#1f6fb4", lw=1.0, alpha=alpha, zorder=3)
        else:
            ax.scatter(channels.geometry.x, channels.geometry.y,
                       s=0.5, c="#1f6fb4", alpha=alpha, zorder=3)
        return

    # Raster fallback: a 0/1 mask on the analysis grid.
    masked = np.ma.masked_where(np.asarray(channels) < 0.5, np.asarray(channels))
    ax.imshow(masked, cmap=matplotlib.colors.ListedColormap(["#1f6fb4"]),
              extent=extent, alpha=alpha, interpolation="nearest", zorder=3)


def _draw_base(ax, dem, channels, cell, extent, channel_alpha=0.85):
    masked_dem = np.ma.masked_invalid(np.where(dem < -1e30, np.nan, dem))
    filled = masked_dem.filled(np.nan)
    ax.imshow(_hillshade(np.nan_to_num(filled, nan=float(np.nanmin(filled))), cell),
              cmap="Greys_r", extent=extent,
              vmin=0, vmax=1, interpolation="bilinear")
    _draw_channels(ax, channels, extent, channel_alpha)
    ax.set_xlabel("easting (m)")
    ax.set_ylabel("northing (m)")


# --------------------------------------------------------------------------
# Group A: does the tracing work?
# --------------------------------------------------------------------------

def plot_traces_on_hillshade(outputs: HilltopCouplingOutputs, *,
                             cell: Optional[float] = None,
                             max_traces: Optional[int] = None,
                             figsize=(11, 9)):
    """Traces from ridge pixels to the channel network, over a hillshade.

    The single most valuable diagnostic in the set: traces that cross a divide,
    run uphill, or die mid-slope are obvious here and invisible in any summary
    statistic.

    Works from either output shape -- LineString features from the CLI or the
    development driver's long-form CSV -- and plots in the raster's own
    projected coordinates so no transform is needed.

    Args:
        outputs: A loaded output set; needs trace geometry and the DEM raster.
        cell: Cell size in metres. Read from the DEM when omitted.
        max_traces: Draw at most this many, for legibility on large runs.
        figsize: Figure size in inches.

    Returns:
        The figure.
    """
    plt = _import_mpl()

    dem, transform, extent, cell = _read_grid_extent(outputs.dem_path, cell)
    channels = _load_channels(outputs.channels_path)

    fig, ax = plt.subplots(figsize=figsize)
    _draw_base(ax, dem, channels, cell, extent)

    traces = outputs.any_traces()
    id_col = "trace_id" if (traces is not None and "trace_id" in traces.columns) \
        else "source_row"
    seen_outcomes = set()
    n_drawn = 0
    for _, xs, ys, attrs in line_vertices(traces, id_col):
        if max_traces is not None and n_drawn >= max_traces:
            break
        outcome = int(attrs.get("outcome", 0)) if hasattr(attrs, "get") else 0
        seen_outcomes.add(outcome)
        ax.plot(xs, ys, lw=0.6, alpha=0.85,
                color=OUTCOME_COLOUR.get(outcome, "#333333"))
        n_drawn += 1

    handles = [plt.Line2D([], [], color="#1f6fb4", lw=3, label="channel network")]
    for outcome in sorted(seen_outcomes):
        handles.append(plt.Line2D([], [], color=OUTCOME_COLOUR[outcome], lw=1.5,
                                  label=f"trace: {OUTCOME_NAME[outcome]}"))
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.9)

    n_traces = int(outputs.summary.get("n_traces", 0)) or n_drawn
    ax.set_title(
        "Hilltop flow routing: traces from ridge pixels to the channel network\n"
        f"{n_drawn} traces drawn (thinned on write)", fontsize=10)

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
    traces = outputs.trace_table()
    if traces is None or not len(traces):
        raise ValueError("plot_trace_outcomes: the output set carries no traces")
    # The product CLI keeps the outcome on the hilltop points as trace_status,
    # and names the source pixel row/col rather than source_row/source_col.
    outcome_col = "outcome" if "outcome" in traces.columns else "trace_status"
    row_col = "source_row" if "source_row" in traces.columns else "row"
    col_col = "source_col" if "source_col" in traces.columns else "col"
    lh_col = "length_m" if "length_m" in traces.columns else "hillslope_length_m"

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

    for outcome, group in traces.groupby(outcome_col):
        x, y = _rc_to_metres(group[row_col], group[col_col], cell)
        axes[1].scatter(x, y, s=8, c=OUTCOME_COLOUR[int(outcome)],
                        label=OUTCOME_NAME[int(outcome)])
    axes[1].invert_yaxis()
    axes[1].set_aspect("equal")
    axes[1].legend(fontsize=7)
    axes[1].set_title("Where failures cluster")
    axes[1].set_xlabel("easting (m)")
    axes[1].set_ylabel("northing (m)")

    good = outputs.successful_traces()
    axes[2].hist(good[lh_col], bins=30, color="#2a7f3f", alpha=0.75,
                 label=f"Lh (median {good[lh_col].median():.0f} m)")
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
                              figsize=(17, 8)):
    """The ordered divides, and the hilltop pixels they separate.

    This is where the divide-grid indexing is checked by eye: the hilltop pixels
    must straddle each ridgeline, one colour on each flank. A panel where one
    colour dominates, or where the two interleave, means the side assignment has
    collapsed.

    Args:
        outputs: A loaded output set; needs ``ridgelines``, ``hilltops`` and the
            DEM raster.
        cell: Cell size in metres. Read from the DEM when omitted.
        figsize: Figure size in inches.

    Returns:
        The figure.
    """
    plt = _import_mpl()

    dem, transform, extent, cell = _read_grid_extent(outputs.dem_path, cell)
    channels = _load_channels(outputs.channels_path)

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    for ax in axes:
        _draw_base(ax, dem, channels, cell, extent, channel_alpha=0.7)

    cmap = plt.get_cmap("tab20")
    n_ridgelines = 0
    for rid, xs, ys, _ in line_vertices(outputs.ridgelines, "ridgeline_id"):
        axes[0].plot(xs, ys, lw=1.0, color=cmap(int(rid) % 20))
        n_ridgelines += 1

    axes[0].set_title(
        "Ordered ridgelines (the divides themselves)\n"
        f"{n_ridgelines} ridgelines on the pixel-corner divide grid, "
        "coloured by ridgeline_id", fontsize=10)

    n_hilltops = 0
    hilltops = outputs.hilltops
    side_counts = {}
    if hilltops is not None and len(hilltops):
        n_hilltops = len(hilltops)
        xs_all, ys_all = _hilltop_xy(hilltops, transform, cell)
        for side in (0, 1):
            sel = (hilltops["ridge_side"] == side).to_numpy()
            side_counts[side] = int(sel.sum())
            if not sel.any():
                continue
            axes[1].scatter(xs_all[sel], ys_all[sel], s=1.0,
                            c=SIDE_COLOUR[side], linewidths=0,
                            label=f"side {side} ({SIDE_NAME[side]}), "
                                  f"{side_counts[side]}")
        axes[1].legend(loc="upper right", fontsize=8, framealpha=0.9, markerscale=8)

        # The two flanking pixels of a divide are ONE CELL apart, so at full
        # extent they overlap and whichever is drawn second appears to dominate.
        # The inset is the only place the pairing is actually visible, and the
        # counts are in the title so the figure states the balance even where
        # the plot cannot show it.
        if n_hilltops:
            import numpy as _np
            cx = float(_np.median(xs_all))
            cy = float(_np.median(ys_all))
            half = 60.0 * cell
            inset = axes[1].inset_axes([0.02, 0.02, 0.30, 0.34])
            near = ((_np.abs(xs_all - cx) < half) & (_np.abs(ys_all - cy) < half))
            for side in (0, 1):
                sel = near & (hilltops["ridge_side"] == side).to_numpy()
                if sel.any():
                    inset.scatter(xs_all[sel], ys_all[sel], s=6,
                                  c=SIDE_COLOUR[side], linewidths=0)
            for _, lxs, lys, _a in line_vertices(outputs.ridgelines, "ridgeline_id"):
                inset.plot(lxs, lys, color="0.35", lw=0.6)
            inset.set_xlim(cx - half, cx + half)
            inset.set_ylim(cy - half, cy + half)
            inset.set_xticks([]); inset.set_yticks([])
            inset.set_title(f"{int(2*half)} m zoom: pixels straddle the divide",
                            fontsize=7)
            for spine in inset.spines.values():
                spine.set_edgecolor("k")

    balance = ("" if not side_counts else
               f"  ({side_counts.get(0, 0)} left / {side_counts.get(1, 0)} right)")
    axes[1].set_title(
        "Hilltop pixels: the pixels the ridgelines separate\n"
        f"{n_hilltops} pixels; the two colours are the opposing hillslopes "
        f"of each divide{balance}", fontsize=10)

    fig.tight_layout()
    return fig


def plot_length_distributions(outputs: HilltopCouplingOutputs, *,
                              figsize=(15, 4.4)):
    """Hillslope length, relief and their ratio, over successful traces only.

    Every non-``reached channel`` outcome leaves ``Lh`` meaningless, so the
    distributions are taken over that subset and the discarded count is put in
    the title rather than left implicit -- a median hillslope length quoted
    without the fraction of hilltops it was computed from is not reusable.

    Args:
        outputs: A loaded output set; needs the hilltop table.
        figsize: Figure size in inches.

    Returns:
        The figure.
    """
    plt = _import_mpl()

    table = outputs.trace_table()
    if table is None or not len(table):
        raise ValueError("plot_length_distributions: no hilltop table")

    good = outputs.successful_traces()
    if not len(good):
        raise ValueError("plot_length_distributions: no successful traces")

    lh_col = "hillslope_length_m" if "hillslope_length_m" in good.columns else "length_m"
    lh = pd.to_numeric(good[lh_col], errors="coerce").dropna()
    relief = pd.to_numeric(good["relief_m"], errors="coerce").dropna()

    fig, axes = plt.subplots(1, 3, figsize=figsize)

    axes[0].hist(lh, bins=60, color="#2a7f3f", alpha=0.85)
    axes[0].axvline(lh.median(), color="k", ls="--", lw=1.5,
                    label=f"median {lh.median():.0f} m")
    axes[0].set_xlabel("hillslope length Lh (m)")
    axes[0].set_ylabel("hilltop pixels")
    axes[0].legend(fontsize=8)
    axes[0].set_title("Hillslope length")

    axes[1].hist(relief, bins=60, color="#8c564b", alpha=0.85)
    axes[1].axvline(relief.median(), color="k", ls="--", lw=1.5,
                    label=f"median {relief.median():.1f} m")
    axes[1].set_xlabel("hillslope relief R (m)")
    axes[1].legend(fontsize=8)
    axes[1].set_title("Hillslope relief")

    # Sinuosity needs no ground truth and diagnoses a wandering tracer.
    if "euclidean_distance_m" in good.columns:
        euc = pd.to_numeric(good["euclidean_distance_m"], errors="coerce")
        ratio = (pd.to_numeric(good[lh_col], errors="coerce") / euc).replace(
            [np.inf, -np.inf], np.nan).dropna()
        ratio = ratio[ratio < 5]
        axes[2].hist(ratio, bins=60, color="#1f77b4", alpha=0.85)
        axes[2].axvline(1.0, color="k", lw=1.5)
        axes[2].set_xlabel("Lh / straight-line distance")
        axes[2].set_title(f"Path sinuosity (median {ratio.median():.2f})")
    else:
        axes[2].axis("off")

    n_all = len(table)
    fig.suptitle(
        f"Hillslope length and relief, over the {len(good)} of {n_all} hilltop "
        f"pixels ({100.0 * len(good) / n_all:.0f} %) whose trace reached a channel",
        fontsize=11, y=1.03)
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

    # Derive the divide length from the along-divide coordinate rather than
    # requiring a total-length column: the product CLI does not carry one on
    # the points, and the maximum distance is the same number.
    length_col = ("ridgeline_total_length_m"
                  if "ridgeline_total_length_m" in hilltops.columns
                  else "distance_along_ridgeline_m")
    lengths = (hilltops.groupby("ridgeline_id")[length_col]
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

def plot_looped_trace_profiles(outputs: HilltopCouplingOutputs, *,
                               dem_path=None,
                               n_traces: int = 6,
                               near_channel_m: Optional[float] = None,
                               figsize=(13, 9),
                               seed: int = 0):
    """Elevation profiles of traces that terminated as ``loop_detected``.

    A loop outcome says only that the tracer revisited a cell; it does not say
    the path was wrong. This figure separates the two by drawing the elevation
    along each looping path and marking, in red, the stretch from the FIRST
    revisited cell onward -- the part the loop detector actually objected to.

    Read it as:

    * profile descends steadily, red stretch is a short flat wobble at the toe
      -> the path is fine and the loop is the tracer oscillating between two
      facets of a valley floor it has already reached;
    * red stretch climbs, or wanders over a long distance -> the path genuinely
      went somewhere it should not have.

    Args:
        outputs: The loaded output set; its ``trace_lines`` supply the geometry.
        dem_path: DEM to sample elevation from. Defaults to ``outputs.dem_path``.
        n_traces: How many looping traces to draw.
        near_channel_m: When set, keep only traces whose source is within this
            distance of a channel cell. Loops far from any channel are
            uninformative -- there is no channel for them to reach -- so this
            restricts the figure to the residual worth explaining.
        seed: Sampling seed, so the same traces are drawn on a rerun.

    Returns:
        The matplotlib Figure.
    """
    plt = _import_mpl()

    traces = outputs.trace_lines
    if traces is None or len(traces) == 0:
        raise ValueError("no trace geometry in this output set")
    if "outcome" not in traces.columns:
        raise ValueError("traces carry no outcome column")

    looped = traces[traces["outcome"] == 2]
    if len(looped) == 0:
        raise ValueError("no looping traces to draw")

    path = dem_path if dem_path is not None else outputs.dem_path
    if path is None:
        raise ValueError("a DEM is needed to sample elevation along the traces")
    grid, _transform, extent, cell = _read_grid_extent(path, None)
    grid = np.where(grid <= -9990, np.nan, grid.astype(float))

    if near_channel_m is not None and outputs.channels_path is not None:
        from scipy.ndimage import distance_transform_edt
        chan = _load_channels(outputs.channels_path)
        if isinstance(chan, np.ndarray):
            dist = distance_transform_edt(~(chan == 1)) * cell
            rr = looped["source_row"].to_numpy().astype(int)
            cc = looped["source_col"].to_numpy().astype(int)
            inb = (rr >= 0) & (rr < dist.shape[0]) & (cc >= 0) & (cc < dist.shape[1])
            keep = np.zeros(len(looped), bool)
            keep[inb] = dist[rr[inb], cc[inb]] <= near_channel_m
            if keep.any():
                looped = looped[keep]

    rng = np.random.default_rng(seed)
    take = min(n_traces, len(looped))
    picked = looped.iloc[rng.choice(len(looped), take, replace=False)]

    ncol = 2
    nrow = int(np.ceil(take / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)

    x0, x1, y0, y1 = extent
    for ax, (_, row) in zip(axes.ravel(), picked.iterrows()):
        geom = row["geometry"]
        parts = geom.geoms if geom.geom_type.startswith("Multi") else [geom]
        xs, ys = [], []
        for part in parts:
            px, py = part.xy
            xs.extend(px)
            ys.extend(py)
        xs = np.asarray(xs)
        ys = np.asarray(ys)

        col = np.clip(((xs - x0) / cell).astype(int), 0, grid.shape[1] - 1)
        r_ = np.clip(((y1 - ys) / cell).astype(int), 0, grid.shape[0] - 1)
        z = grid[r_, col]

        step = np.hypot(np.diff(xs), np.diff(ys))
        dist = np.concatenate([[0.0], np.cumsum(step)])

        # The first revisited CELL is what the tracer's own loop test fires on,
        # so the split is made on cell indices, not on coordinates.
        seen = {}
        first_repeat = len(xs)
        for i, key in enumerate(zip(r_.tolist(), col.tolist())):
            if key in seen:
                first_repeat = i
                break
            seen[key] = i

        ax.plot(dist[:first_repeat], z[:first_repeat], "-", color="#2a7f3f", lw=1.6,
                label="before first revisit")
        if first_repeat < len(xs):
            ax.plot(dist[first_repeat - 1:], z[first_repeat - 1:], "-",
                    color="#d62728", lw=1.6, label="after first revisit")
            ax.axvline(dist[first_repeat], color="#d62728", ls=":", lw=1)

        drop = float(z[0] - np.nanmin(z))
        ax.set_title(f"source ({int(row['source_row'])},{int(row['source_col'])})  "
                     f"path {dist[-1]:.0f} m, drop {drop:.1f} m",
                     fontsize=9)
        ax.set_xlabel("distance along path (m)", fontsize=8)
        ax.set_ylabel("elevation (m)", fontsize=8)
        ax.tick_params(labelsize=8)

    for ax in axes.ravel()[take:]:
        ax.axis("off")
    axes.ravel()[0].legend(fontsize=8)

    near = "" if near_channel_m is None else \
        f", sources within {near_channel_m:.0f} m of a channel"
    fig.suptitle("Traces that terminated as 'loop detected'"
                 f"{near}\nred = from the first revisited cell onward",
                 fontsize=11)
    fig.tight_layout()
    return fig
