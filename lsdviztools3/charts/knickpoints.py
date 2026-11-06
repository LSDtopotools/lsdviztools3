"""Knickpoint plots for lsdtt3 chi fits.

lsdtt-chi-analysis (``write_chi_knickpoint_table``) writes one row per node
where two fitted k_sn segments meet (``*_chi_knickpoints_mn_*.csv``), with
``delta_ksn`` (upstream minus downstream k_sn: positive where the channel
steepens upstream) and ``log_ratio_ksn``, and a table of suggested
thresholds for a "major" knickpoint (``*_chi_knickpoint_thresholds_mn_*.csv``:
p90, p95, median + 3 MAD, Otsu, over the absolute values). All the numbers
come from lsdtt3; this module selects rows and draws them.

``select_major_knickpoints``  rows beyond a threshold (a value, or a method name)
``plot_knickpoint_distribution``  the signed distribution with the suggested
    thresholds, and the ranked absolute values (for choosing by eye)
``plot_knickpoints_on_profile``  major knickpoints on a single-basin or the
    stacked chi profile
``plot_knickpoints_along_chi``  where major knickpoints sit along chi (or
    normalised chi, or elevation), against all candidates
``render_knickpoint_map``  major knickpoints on the hillshade and k_sn network
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

MEASURES = {"delta": ("delta_ksn", "abs_delta_ksn", r"$\Delta k_{sn}$ (upstream $-$ downstream)"),
            "log_ratio": ("log_ratio_ksn", "abs_log_ratio_ksn",
                          r"$\log_{10}(k_{sn,up}/k_{sn,down})$")}
METHODS = ("p90", "p95", "median_plus_3mad", "otsu")
# Validated together (light surface); every line is labelled in the legend.
METHOD_COLOURS = {"p90": "#2a78d6", "p95": "#1baf7a", "median_plus_3mad": "#9a5bd6", "otsu": "#eb6834"}
INK = "#3d3d3a"
ALL_COLOUR = "0.82"


def _measure(measure: str):
    if measure not in MEASURES:
        raise ValueError(f"measure must be one of {list(MEASURES)}, got {measure!r}")
    return MEASURES[measure]


def resolve_threshold(threshold, measure: str = "delta", thresholds: Optional[pd.DataFrame] = None
                      ) -> Tuple[float, str]:
    """``(value, label)`` for a threshold given as a number or a method name
    (``p90``, ``p95``, ``median_plus_3mad``, ``otsu``) from the lsdtt3
    thresholds table, on the absolute ``measure``."""
    _, abs_name, _ = _measure(measure)
    try:
        v = float(threshold)
        return v, f"{v:g} (user)"
    except (TypeError, ValueError):
        pass
    if thresholds is None:
        raise ValueError(f"threshold {threshold!r} is a method name: give the thresholds table.")
    hit = thresholds[(thresholds["measure"] == abs_name) & (thresholds["method"] == threshold)]
    if hit.empty or not np.isfinite(hit["threshold"].iloc[0]):
        raise ValueError(f"no {abs_name} threshold named {threshold!r} in the thresholds table")
    v = float(hit["threshold"].iloc[0])
    return v, f"{v:.3g} ({threshold})"


def select_major_knickpoints(
    knickpoints: pd.DataFrame,
    threshold,
    *,
    measure: str = "delta",
    direction: str = "both",
    include_junctions: bool = False,
    thresholds: Optional[pd.DataFrame] = None,
    removed: Optional[pd.DataFrame] = None,
) -> Tuple[pd.DataFrame, float, str]:
    """Select the major knickpoints.

    Args:
        knickpoints: The lsdtt3 knickpoint table.
        threshold: A value (in the measure's units) or a method name read from
            ``thresholds``.
        measure: ``"delta"`` (``delta_ksn``, default) or ``"log_ratio"``.
        direction: ``"both"`` (|value| above), ``"steepening"`` (value above
            the threshold) or ``"shallowing"`` (value below minus it).
        include_junctions: Keep rows at tributary mouths (``at_junction``).
        removed: Tributaries removed from the plots
            (:func:`lsdviztools3.io.chi.select_tributaries`); their knickpoints
            are dropped too.

    Returns:
        ``(rows, value, label)``: the selected rows (with a ``measure_value``
        column), the threshold value and a label for it.
    """
    col, _, _ = _measure(measure)
    if direction not in ("both", "steepening", "shallowing"):
        raise ValueError(f"direction must be both, steepening or shallowing, got {direction!r}")
    value, label = resolve_threshold(threshold, measure, thresholds)
    k = candidate_knickpoints(knickpoints, include_junctions=include_junctions, removed=removed)
    v = k[col].to_numpy(dtype=float)
    if direction == "both":
        keep = np.abs(v) > value
    elif direction == "steepening":
        keep = v > value
    else:
        keep = v < -value
    out = k[keep & np.isfinite(v)].copy()
    out["measure_value"] = out[col]
    return out, value, label


def candidate_knickpoints(knickpoints, *, include_junctions=False, removed=None) -> pd.DataFrame:
    """The rows the statistics and selections use: without junction rows
    (unless asked) and without removed tributaries."""
    k = knickpoints
    if not include_junctions:
        k = k[k["at_junction"] == 0]
    if removed is not None and len(removed):
        from ..io.chi import removed_mask
        k = k[~removed_mask(k, removed)]
    return k


def plot_knickpoint_distribution(
    knickpoints: pd.DataFrame,
    thresholds: Optional[pd.DataFrame] = None,
    *,
    measure: str = "delta",
    chosen: Optional[float] = None,
    chosen_label: Optional[str] = None,
    include_junctions: bool = False,
    removed: Optional[pd.DataFrame] = None,
    bins: int = 60,
):
    """The distribution of k_sn changes, for choosing a threshold.

    Left: a histogram of the signed measure (steepening right of 0, shallowing
    left), counts on a log axis so the tails show. Right: the absolute values
    ranked largest first, where a break in slope can be picked by eye. Both
    carry the suggested thresholds from lsdtt3 (one colour per method, its n
    above in the legend) and the chosen threshold (black, ±).

    Returns:
        ``(fig, (ax_hist, ax_ranked))``.
    """
    import matplotlib.pyplot as plt

    col, abs_name, xlabel = _measure(measure)
    k = candidate_knickpoints(knickpoints, include_junctions=include_junctions, removed=removed)
    v = k[col].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    fig, (ah, ar) = plt.subplots(1, 2, figsize=(11, 4), gridspec_kw={"width_ratios": [1.4, 1]})
    if v.size == 0:
        raise ValueError("no knickpoint candidates to plot")
    ah.hist(v, bins=bins, color=ALL_COLOUR, edgecolor="white", linewidth=0.5)
    ah.set_yscale("log")
    ah.set_xlabel(xlabel)
    ah.set_ylabel("number of segment boundaries")
    a = np.sort(np.abs(v))[::-1]
    ar.plot(np.arange(1, a.size + 1), a, color=INK, lw=1.2)
    ar.set_xscale("log")
    ar.set_xlabel("rank (largest first)")
    ar.set_ylabel(f"|{xlabel}|")
    handles = []
    if thresholds is not None:
        for method in METHODS:
            hit = thresholds[(thresholds["measure"] == abs_name) & (thresholds["method"] == method)]
            if hit.empty or not np.isfinite(hit["threshold"].iloc[0]):
                continue
            t = float(hit["threshold"].iloc[0])
            n = int((np.abs(v) > t).sum())
            c = METHOD_COLOURS[method]
            for x in (-t, t):
                ah.axvline(x, color=c, lw=1.0, ls="--")
            ar.axhline(t, color=c, lw=1.0, ls="--")
            handles.append(plt.Line2D([], [], color=c, ls="--", label=f"{method}: {t:.3g} ({n} above)"))
    if chosen is not None:
        n = int((np.abs(v) > chosen).sum())
        for x in (-chosen, chosen):
            ah.axvline(x, color="black", lw=1.8)
        ar.axhline(chosen, color="black", lw=1.8)
        handles.append(plt.Line2D([], [], color="black", lw=1.8,
                                  label=f"chosen: {chosen_label or f'{chosen:g}'} ({n} above)"))
    for ax in (ah, ar):
        ax.grid(alpha=0.25, lw=0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    ah.set_title(f"{v.size} segment boundaries" + ("" if include_junctions else " (tributary mouths left out)"),
                 fontsize=9, loc="left")
    if handles:
        ar.legend(handles=handles, fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    return fig, (ah, ar)


def _marker_sizes(values: np.ndarray, lo: float = 20.0, hi: float = 120.0) -> np.ndarray:
    a = np.abs(np.asarray(values, dtype=float))
    if a.size == 0:
        return a
    amin, amax = float(np.nanmin(a)), float(np.nanmax(a))
    if not amax > amin:
        return np.full(a.size, 0.5 * (lo + hi))
    return lo + (hi - lo) * (a - amin) / (amax - amin)


def _draw_markers(ax, x, y, values, zorder=6):
    v = np.asarray(values, dtype=float)
    s = _marker_sizes(v)
    up = v > 0
    h = []
    if up.any():
        ax.scatter(np.asarray(x)[up], np.asarray(y)[up], s=s[up], marker="^", facecolor="white",
                   edgecolor="black", linewidth=0.9, zorder=zorder)
        h.append(("^", "steepening upstream"))
    if (~up).any():
        ax.scatter(np.asarray(x)[~up], np.asarray(y)[~up], s=s[~up], marker="v", facecolor="0.25",
                   edgecolor="black", linewidth=0.9, zorder=zorder)
        h.append(("v", "shallowing upstream"))
    return h


def _marker_legend(ax, kinds, label_extra=""):
    import matplotlib.pyplot as plt
    hs = [plt.Line2D([], [], marker=m, ls="none", markerfacecolor=("white" if m == "^" else "0.25"),
                     markeredgecolor="black", markersize=7, label=lab + label_extra) for m, lab in kinds]
    if hs:
        ax.legend(handles=hs, fontsize=8, frameon=False, loc="upper left")


def plot_knickpoints_on_profile(
    segments: pd.DataFrame,
    major: pd.DataFrame,
    *,
    basin: Optional[int] = None,
    basins: Optional[Sequence[int]] = None,
    removed: Optional[pd.DataFrame] = None,
    show_removed: str = "ghost",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    threshold_label: str = "",
):
    """Major knickpoints on the chi profile: one basin (``basin=``) or the
    stacked plot of ``basins`` (default all). Markers: open upward triangles
    where the channel steepens upstream, filled downward where it shallows,
    sized by the size of the change.

    Returns:
        ``(fig, ax)``.
    """
    from .chi import plot_chi_profile, plot_chi_stacked, stack_offsets

    if basin is not None:
        fig, ax = plot_chi_profile(segments, basin=basin, removed=removed, show_removed=show_removed,
                                   vmin=vmin, vmax=vmax)
        m = major[major["basin_index"] == basin]
        kinds = _draw_markers(ax, m["chi"], m["elevation"], m["measure_value"])
    else:
        sel_basins = sorted(segments["basin_index"].unique()) if basins is None else list(basins)
        fig, ax = plot_chi_stacked(segments, basins=sel_basins, removed=removed,
                                   show_removed=show_removed, vmin=vmin, vmax=vmax)
        shifts = stack_offsets(segments, sel_basins)
        m = major[major["basin_index"].isin(list(shifts))]
        x = m["chi"].to_numpy(dtype=float) + m["basin_index"].map(shifts).to_numpy(dtype=float)
        kinds = _draw_markers(ax, x, m["elevation"], m["measure_value"])
    _marker_legend(ax, kinds)
    ax.set_title((ax.get_title() + "  " if basin is not None else "")
                 + f"{len(m)} major knickpoints" + (f", threshold {threshold_label}" if threshold_label else ""),
                 fontsize=9)
    return fig, ax


def plot_knickpoints_along_chi(
    knickpoints: pd.DataFrame,
    major: pd.DataFrame,
    *,
    along: str = "chi",
    include_junctions: bool = False,
    removed: Optional[pd.DataFrame] = None,
    bins: int = 40,
):
    """Where the major knickpoints sit: a histogram of ``along`` (``chi``,
    ``chi_normalised`` or ``elevation``) of the major knickpoints over that
    of all candidates (grey), and below, one row per basin with each major
    knickpoint as a marker. A cluster at one chi value across tributaries or
    basins is what a migrating wave of incision would leave.

    Returns:
        ``(fig, (ax_hist, ax_strip))``.
    """
    import matplotlib.pyplot as plt

    if along not in ("chi", "chi_normalised", "elevation"):
        raise ValueError(f"along must be chi, chi_normalised or elevation, got {along!r}")
    k = candidate_knickpoints(knickpoints, include_junctions=include_junctions, removed=removed)
    allv = k[along].to_numpy(dtype=float)
    allv = allv[np.isfinite(allv)]
    mv = major[along].to_numpy(dtype=float)
    edges = np.histogram_bin_edges(allv if allv.size else mv, bins=bins)
    fig, (ah, asr) = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                                  gridspec_kw={"height_ratios": [1.3, 1]})
    ah.hist(allv, bins=edges, color=ALL_COLOUR, edgecolor="white", linewidth=0.5,
            label=f"all segment boundaries ({allv.size})")
    ah.hist(mv[np.isfinite(mv)], bins=edges, color=INK, edgecolor="white", linewidth=0.5,
            label=f"major knickpoints ({np.isfinite(mv).sum()})")
    ah.set_ylabel("count")
    ah.legend(fontsize=8, frameon=False)
    basins_sorted = sorted(k["basin_index"].unique())
    pos = {b: i for i, b in enumerate(basins_sorted)}
    y = major["basin_index"].map(pos).to_numpy(dtype=float)
    kinds = _draw_markers(asr, major[along], y, major["measure_value"])
    asr.set_yticks(range(len(basins_sorted)))
    asr.set_yticklabels([str(int(b)) for b in basins_sorted])
    asr.set_ylabel("basin")
    asr.set_ylim(-0.7, len(basins_sorted) - 0.3)
    asr.set_xlabel({"chi": r"$\chi$ (m)", "chi_normalised": r"$\chi$ / basin maximum $\chi$",
                    "elevation": "elevation (m)"}[along])
    _marker_legend(asr, kinds)
    for ax in (ah, asr):
        ax.grid(alpha=0.25, lw=0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig, (ah, asr)


def render_knickpoint_map(
    dem_path,
    segment_lines,
    major: pd.DataFrame,
    png_path,
    *,
    removed: Optional[pd.DataFrame] = None,
    coordinates: str = "local",
    cmap: str = "viridis",
    ksn_vmin: Optional[float] = None,
    ksn_vmax: Optional[float] = None,
    threshold_label: str = "",
    high_res: bool = False,
    legend_position: str = "JTR+jTR+o0.2c",
):
    """Major knickpoints on the hillshade and the k_sn network: open upward
    triangles where the channel steepens upstream, filled downward ones where
    it shallows, sized by the size of the change (log-scaled). Removed
    tributaries are left off the network and their knickpoints are not shown.

    Returns:
        ``Path`` to the written PNG.
    """
    import io as _io

    import geopandas as gpd

    from ..io.chi import removed_mask
    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from ..io.vector import load_vector
    from ..render.base import grid_to_display, save_figure
    from ..render.channels import render_channels
    from ..render.decorations import add_scalebar
    from ..render.hillshade import render_hillshade
    from ..render.points import render_points
    from ..render.style import MapStyle

    width = MapStyle.MAX_WIDTH
    da = to_coordinates(load_raster(dem_path), coordinates)
    crs = get_crs(da)
    _, extent = grid_to_display(to_gmt_grid(da), crs)
    gdf = load_vector(str(segment_lines))
    if removed is not None and len(removed):
        gdf = gdf[~removed_mask(gdf, removed)]
    fig = render_hillshade(str(dem_path), coordinates=coordinates, shade=True,
                           style=MapStyle(cmap="gray", figure_size=width, transparent=False))
    fig = render_channels(gdf, color_by="k_sn", log_color=True, cmap=cmap, colorbar=True,
                          colorbar_label="log@-10@-k@-sn@-", coordinates=coordinates,
                          style=MapStyle(line_width="0.9p", cmap=cmap, figure_size=width),
                          fig=fig, vmin=ksn_vmin, vmax=ksn_vmax)
    rows = []
    m = major.dropna(subset=["x", "y"]).copy()
    m["abs_value"] = np.abs(m["measure_value"].to_numpy(dtype=float))
    for sign, symbol, fill, label in ((1, "t", "white", "steepening upstream"),
                                      (-1, "i", "gray25", "shallowing upstream")):
        part = m[np.sign(m["measure_value"]) == sign]
        if part.empty:
            continue
        pts = gpd.GeoDataFrame(part[["abs_value"]].copy(),
                               geometry=gpd.points_from_xy(part["x"], part["y"]), crs=gdf.crs)
        fig = render_points(pts, size_by="abs_value", size_range=(5.0, 12.0), symbol=symbol,
                            pen="0.8p,black", coordinates=coordinates,
                            style=MapStyle(line_color=fill, figure_size=width), fig=fig)
        rows.append(f"S 0.3c {symbol} 0.25c {fill} 0.8p,black 0.9c {label} ({len(part)})")
    if rows:
        if threshold_label:
            rows.insert(0, f"L 9p - L major knickpoints, threshold {threshold_label}")
        fig.legend(spec=_io.StringIO("\n".join(rows)), position=legend_position,
                   box="+gwhite@20+p0.5p,black")
    add_scalebar(fig, extent, crs)
    save_figure(fig, png_path, style=MapStyle(figure_size=width, transparent=False), high_res=high_res)
    return Path(png_path)
