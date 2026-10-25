"""Ridgeline profile plotting (matplotlib).

Charts of ``lsdtt-ridgeline-analysis`` profiles, ported from
R. Vervoordeldonk's ``ridgeline_plotting`` package:

* :func:`plot_ridgeline_detrending_curvature` -- four stacked panels against
  distance along the ridge crest: elevation with its fitted trend,
  moving-window detrended elevations, maximal curvature and the Laplacian at
  chosen polyfit scales.
* :func:`plot_ridgeline_profile_summary` -- native and topological triple
  divide junctions by threshold, the ridgeline series with a maximal-curvature
  sign barcode at every polyfit scale, and a key to the barcode. The map that
  went with it is :func:`lsdviztools3.render.ridgeline.render_ridgeline_map`.
* :func:`plot_mdd_ridgeline_map_profile` -- a main-drainage-divide profile
  above a drainage-area image of an east-west periodic model, with guides at
  fixed distances through both.
* :func:`plot_zero_spacing_slope_area` -- per site, the spread of maximal
  curvature zero-crossing spacing against polyfit scale beside the slope-area
  distribution, with the characteristic length L.

Inputs come from :mod:`lsdviztools3.io.ridgeline`. Every function returns
``(fig, axes)``; :func:`save_figure_formats` writes PNG, SVG and PDF.
matplotlib is imported lazily (DESIGN decision #1).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from ..io.ridgeline import (
    characteristic_length,
    lambda_column,
    lambda_columns,
    moving_window_column,
    node_spacing,
    periodic_median_crossing_spacing,
    slope_area_peak,
    tdj_display_positions,
)

_MPL_MSG = (
    "Ridgeline plotting needs matplotlib. Install the optional 'charts' extra, "
    "e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Colours shared by the ridgeline charts.
RIDGE_GREEN = "#178a3b"
TDJ_RED = "#d62728"
NEGATIVE_CURVATURE = "#006d2c"
POSITIVE_CURVATURE = "#74c476"
ZERO_CURVATURE = "#777777"
GUIDE_WIDTH = 0.6


def _import_mpl():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def save_figure_formats(fig, write_path, write_name: str, *, png_dpi: int = 300,
                        vector_raster_dpi: Optional[float] = None) -> list:
    """Write a figure as ``<write_path>/{png,svg,pdf}/<write_name>.<ext>``.

    Args:
        fig: The matplotlib figure.
        write_path: Output directory; the three subdirectories are created.
        write_name: File stem shared by the three files.
        png_dpi: PNG resolution.
        vector_raster_dpi: Resolution of raster images embedded in SVG/PDF
            (None: matplotlib's default).

    Returns:
        The written paths.
    """
    write_name = str(write_name).strip()
    if not write_name or Path(write_name).name != write_name:
        raise ValueError("write_name must be a file stem without path separators")
    root = Path(write_path)
    outputs = []
    for ext in ("png", "svg", "pdf"):
        (root / ext).mkdir(parents=True, exist_ok=True)
        dest = root / ext / f"{write_name}.{ext}"
        options = {"dpi": png_dpi} if ext == "png" else (
            {"dpi": vector_raster_dpi} if vector_raster_dpi else {})
        fig.savefig(dest, **options)
        outputs.append(dest)
    return outputs


def _style_axis(axis) -> None:
    axis.tick_params(width=GUIDE_WIDTH)
    for spine in axis.spines.values():
        spine.set_color("black")
        spine.set_linewidth(GUIDE_WIDTH)


def _km_ticks(distance) -> np.ndarray:
    return np.arange(math.floor(np.min(distance) / 1000) * 1000, np.max(distance) + 1, 1000)


def add_curvature_barcode(ax, distance, curvature, *, scale_labels: Sequence = (),
                          label_every: float = 25.0, alpha: float = 0.28):
    """Draw maximal-curvature signs as a barcode behind a profile.

    One row per polyfit scale (rows top to bottom), dark green where the
    curvature is negative (peaks), light green where positive (saddles), grey
    where zero, on a twin y axis.

    Args:
        ax: The profile axes.
        distance: Display distances (m), increasing.
        curvature: Array ``(n_scales, n_nodes)`` (or one row).
        scale_labels: Lambda of each row, for the twin-axis ticks.
        label_every: Only lambdas that are multiples of this are labelled.
        alpha: Barcode opacity.

    Returns:
        The twin axes holding the barcode.
    """
    _import_mpl()
    from matplotlib.colors import BoundaryNorm, ListedColormap

    curvature = np.atleast_2d(np.asarray(curvature, dtype=float))
    signs = np.where(curvature > 0, 1, np.where(curvature < 0, -1, 0))
    cmap = ListedColormap([NEGATIVE_CURVATURE, ZERO_CURVATURE, POSITIVE_CURVATURE])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    twin = ax.twinx()
    twin.set_zorder(0)
    ax.set_zorder(2)
    ax.patch.set_visible(False)
    n = signs.shape[0]
    twin.imshow(signs, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest",
                alpha=alpha, extent=(float(np.min(distance)), float(np.max(distance)),
                                     n - 0.5, -0.5))
    twin.set_ylim(n - 0.5, -0.5)
    labelled = [i for i, lam in enumerate(scale_labels) if lam and lam % label_every == 0]
    twin.set_yticks(labelled, [f"{scale_labels[i]:g}" for i in labelled])
    twin.set_yticks([], minor=True)
    twin.tick_params(axis="y", colors=NEGATIVE_CURVATURE, width=GUIDE_WIDTH, length=3)
    return twin


# -----------------------------------------------------------------------------
# Detrending and curvature
# -----------------------------------------------------------------------------
def plot_ridgeline_detrending_curvature(profile: pd.DataFrame,
                                        windows_pixels: Sequence[int] = (25, 49, 99, 149),
                                        lambdas_m: Sequence[float] = (25, 49, 99, 149),
                                        site_name: Optional[str] = None,
                                        font_size: float = 7):
    """Four stacked panels against distance along the ridge crest.

    1. Elevation (black over a green halo) and its fitted elevation (dashed).
    2. Moving-window detrended elevations, light to dark green with window
       size (legend in metres: window x node spacing).
    3. Maximal curvature at each ``lambdas_m`` (y limits +- the 99th percentile
       of |curvature|).
    4. The Laplacian at the same scales.

    Args:
        profile: One ridgeline after :func:`~lsdviztools3.io.ridgeline.select_profile`
            (uses ``display_distance``).
        windows_pixels: Four moving-window sizes; each needs a
            ``moving_window_detrended_elevation_<N>px`` column
            (``ridgeline_detrend_moving_average_window_pixels`` in
            lsdtt-ridgeline-analysis).
        lambdas_m: Four polyfit scales; each needs ``maximum_curvature_lambda_<L>m``
            and ``laplacian_lambda_<L>m``.
        site_name: Optional label above the top panel (red).
        font_size: Font size of every text element.

    Returns:
        ``(fig, axes)`` with the four axes top to bottom.
    """
    plt = _import_mpl()
    windows = sorted(int(w) for w in windows_pixels)
    lambdas = sorted(float(v) for v in lambdas_m)
    if len(windows) != 4 or len(lambdas) != 4:
        raise ValueError("Exactly four moving-window sizes and four lambdas are required")
    if len(set(windows)) != 4 or len(set(lambdas)) != 4:
        raise ValueError("Moving-window sizes and lambdas must not repeat")
    needed = ["elevation", "fitted_elevation"] + [moving_window_column(w) for w in windows]
    needed += [lambda_column(m, lam) for lam in lambdas for m in ("maximum_curvature", "laplacian")]
    missing = [c for c in needed if c not in profile.columns]
    if missing:
        raise ValueError(f"Profile is missing column(s): {', '.join(missing)}")

    d = profile["display_distance"].to_numpy(dtype=float)
    spacing = node_spacing(d)
    cmap = plt.get_cmap("Greens")
    window_colours = [cmap(v) for v in np.linspace(0.40, 0.85, 4)]
    lambda_colours = [cmap(v) for v in np.linspace(0.40, 0.85, 4)]
    fig, axes = plt.subplots(4, 1, figsize=(7.4, 4 * 3.6 / 2.54), sharex=True,
                             gridspec_kw={"hspace": 0.12})
    z = profile["elevation"].to_numpy(dtype=float)
    axes[0].plot(d, z, color=RIDGE_GREEN, linewidth=2.4, alpha=0.35)
    axes[0].plot(d, z, color="black", linewidth=GUIDE_WIDTH)
    axes[0].plot(d, profile["fitted_elevation"].to_numpy(dtype=float), color="black",
                 linewidth=GUIDE_WIDTH, linestyle="--")
    for colour, w in zip(window_colours, windows):
        axes[1].plot(d, profile[moving_window_column(w)].to_numpy(dtype=float), color=colour,
                     linewidth=0.8, label=f"{w * spacing:g} m")
    axes[1].axhline(0, color="black", linewidth=GUIDE_WIDTH)
    axes[1].legend(loc="upper right", frameon=False, ncol=4, fontsize=font_size * 0.85)
    for colour, lam in zip(lambda_colours, lambdas):
        axes[2].plot(d, profile[lambda_column("maximum_curvature", lam)].to_numpy(dtype=float),
                     color=colour, linewidth=0.7, label=f"{lam:g} m")
        axes[3].plot(d, profile[lambda_column("laplacian", lam)].to_numpy(dtype=float),
                     color=colour, linewidth=0.7, label=f"{lam:g} m")
    k = np.concatenate([profile[lambda_column("maximum_curvature", lam)].to_numpy(dtype=float)
                        for lam in lambdas])
    finite = np.abs(k[np.isfinite(k)])
    limit = float(np.percentile(finite, 99)) if finite.size else 1.0
    if not math.isfinite(limit) or limit == 0:
        limit = 1.0
    axes[2].set_ylim(-limit, limit)
    axes[2].axhline(0, color="black", linewidth=GUIDE_WIDTH)
    for axis in axes[2:]:
        axis.legend(loc="upper right", frameon=False, ncol=4, fontsize=font_size * 0.85)
    axes[3].set_xlabel("Distance along ridge-crest [m]")
    ticks = _km_ticks(d)
    for axis in axes:
        axis.set_xlim(d.min(), d.max())
        axis.set_xticks(ticks)
        axis.grid(axis="x", color="black", linewidth=GUIDE_WIDTH)
        axis.tick_params(labelsize=font_size)
        _style_axis(axis)
    for axis in axes[:-1]:
        axis.tick_params(axis="x", labelbottom=False)

    fig.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.10)
    fig.canvas.draw()
    _, height_px = fig.canvas.get_width_height()
    labels = ("Ridge-crest\nelevation [m]", "Elevation detrended\nwith smoothing window [m]",
              "Maximal\ncurvature [m$^{-1}$]", "Ridgetop curvature,\nLaplacian [m$^{-1}$]")
    for axis, label in zip(axes, labels):
        fig.text(0.035, (axis.bbox.y0 + axis.bbox.y1) / (2 * height_px), label,
                 rotation=90, ha="center", va="center")
    if site_name and str(site_name).strip():
        width_px = fig.canvas.get_width_height()[0]
        fig.text(axes[0].bbox.x0 / width_px, axes[0].bbox.y1 / height_px + 0.008,
                 str(site_name).strip(), ha="left", va="bottom", color=TDJ_RED)
    for text in fig.findobj(match=plt.Text):
        text.set_fontsize(font_size)
    return fig, axes


# -----------------------------------------------------------------------------
# Profile summary
# -----------------------------------------------------------------------------
def plot_ridgeline_profile_summary(profile: pd.DataFrame, native_tdj: Sequence[dict],
                                   topological_tdj: Sequence[dict], *,
                                   series: str = "detrended_elevation",
                                   y_axis_label: str = "Detrended\nridge-crest\nelevation [m]",
                                   highlight: tuple = ("native", 1000.0),
                                   reverse: bool = False,
                                   distance_range: Optional[tuple] = None,
                                   profile_distance=None,
                                   font_size: float = 7):
    """Triple divide junctions, the ridgeline series and its curvature barcode.

    Top to bottom: native TDJs (one row per drainage-area threshold, in km²),
    topological TDJs (one row per distance threshold, m), and the selected
    series (with a maximal-curvature sign barcode at every polyfit scale
    behind it). The highlighted TDJ set is drawn as open circles (topological)
    or omitted from its row (native, as in the original). A key below explains
    the barcode.

    Args:
        profile: One ridgeline after
            :func:`~lsdviztools3.io.ridgeline.select_profile`.
        native_tdj, topological_tdj: Three TDJ sets each
            (:func:`~lsdviztools3.io.ridgeline.load_tdj`); sorted by threshold.
        series: Profile column plotted (``detrended_elevation``,
            ``smoothed_detrended_elevation``, ``elevation``, ...).
        y_axis_label: Label of the series.
        highlight: ``(kind, threshold)`` of the TDJ set to highlight.
        reverse: Whether the profile was reversed (positions of the TDJs).
        distance_range: Displayed distance range (positions of the TDJs).
        profile_distance: ``distance`` of the whole, unclipped ridgeline, the
            anchor of the reverse rule (default: ``profile["distance"]``, right
            unless ``distance_range`` clipped it).
        font_size: Font size of every text element.

    Returns:
        ``(fig, axes)``: ``axes`` is ``{"native", "topological", "profile",
        "barcode", "key"}``.
    """
    plt = _import_mpl()
    from matplotlib.patches import FancyArrowPatch

    native = sorted(native_tdj, key=lambda t: t["threshold"])
    topo = sorted(topological_tdj, key=lambda t: t["threshold"])
    if len(native) != 3 or len(topo) != 3:
        raise ValueError("Three native and three topological TDJ sets are required")
    kind, threshold = highlight[0], float(highlight[1])
    if kind not in ("native", "topological"):
        raise ValueError("highlight kind must be 'native' or 'topological'")
    chosen = native if kind == "native" else topo
    thresholds = [t["threshold"] for t in chosen]
    if threshold not in thresholds:
        raise ValueError(f"Highlighted {kind} threshold {threshold:g} is not among "
                         f"{', '.join(f'{t:g}' for t in thresholds)}")
    selected_index = thresholds.index(threshold)
    if series not in profile.columns:
        raise ValueError(f"Profile has no column {series!r}")
    scales = lambda_columns(profile.columns)
    if not scales:
        raise ValueError("Profile must contain maximum_curvature or "
                         "maximum_curvature_lambda_<N>m columns")

    d = profile["display_distance"].to_numpy(dtype=float)
    y = profile[series].to_numpy(dtype=float)
    if not np.all(np.isfinite(d)) or not np.all(np.isfinite(y)):
        raise ValueError("The selected profile has non-finite distance or series values")
    k = np.vstack([profile[col].to_numpy(dtype=float) for _, col in scales])
    source = np.asarray(profile_distance if profile_distance is not None
                        else profile["distance"], dtype=float)
    shown = distance_range if distance_range is not None else (float(d.min()), float(d.max()))
    native_pos = [tdj_display_positions(t, source, reverse, shown) for t in native]
    topo_pos = [tdj_display_positions(t, source, reverse, shown) for t in topo]

    # The original's layout, with the map row (now a separate PyGMT map) left
    # to the key alone.
    top_height, key_row_height = 6.0 / 2.54, 3.6 / 2.54
    fig = plt.figure(figsize=(6.9, top_height + key_row_height))
    grid = fig.add_gridspec(2, 1, height_ratios=(top_height, key_row_height), hspace=0.3)
    top_grid = grid[0].subgridspec(3, 1, height_ratios=(0.7, 0.7, 1.8), hspace=0.42)
    native_ax = fig.add_subplot(top_grid[0])
    topo_ax = fig.add_subplot(top_grid[1], sharex=native_ax)
    profile_ax = fig.add_subplot(top_grid[2], sharex=native_ax)
    fig.subplots_adjust(left=0.115, right=0.89, top=0.94, bottom=0.04)
    key_row = grid[1].get_position(fig)

    for i in range(3):
        if not (kind == "native" and i == selected_index):
            native_ax.scatter(native_pos[i], np.full(native_pos[i].size, i), color=TDJ_RED,
                              s=13, edgecolors="none", zorder=3)
        if not (kind == "topological" and i == selected_index):
            topo_ax.scatter(topo_pos[i], np.full(topo_pos[i].size, i), color=TDJ_RED,
                            s=13, edgecolors="none", zorder=3)
    native_ax.set_yticks(range(3), [f"{t['threshold'] / 1e6:g}" for t in native])
    topo_ax.set_yticks(range(3), [f"{t['threshold']:g}" for t in topo])
    for axis in (native_ax, topo_ax):
        axis.set_ylim(2.5, -0.5)
        for spine in axis.spines.values():
            spine.set_visible(False)
        for row in range(3):
            axis.axhline(row, color="black", linewidth=GUIDE_WIDTH, zorder=1)
        axis.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    selected_pos = (native_pos if kind == "native" else topo_pos)[selected_index]
    if kind == "topological":
        topo_ax.scatter(selected_pos, np.full(selected_pos.size, selected_index),
                        facecolors="none", edgecolors=TDJ_RED, s=40, linewidths=1.0, zorder=4)
        if selected_pos.size:
            profile_ax.scatter(selected_pos, np.interp(selected_pos, d, y), facecolors="none",
                               edgecolors=TDJ_RED, s=40, linewidths=1.0, zorder=5)

    profile_ax.plot(d, y, color=RIDGE_GREEN, linewidth=2.4, alpha=0.35)
    profile_ax.plot(d, y, color="black", linewidth=GUIDE_WIDTH)
    if series in ("detrended_elevation", "smoothed_detrended_elevation"):
        lim = float(np.nanmax(np.abs(y)))
        if not math.isfinite(lim) or lim == 0:
            lim = 1.0
        profile_ax.set_ylim(-lim, lim)
        profile_ax.axhline(0, color="black", linewidth=GUIDE_WIDTH, zorder=1)
    profile_ax.set_xlabel("Distance along ridge-crest [m]")
    barcode_ax = add_curvature_barcode(profile_ax, d, k,
                                       scale_labels=[lam for lam, _ in scales])
    ticks = _km_ticks(d)
    for axis in (native_ax, topo_ax, profile_ax):
        axis.set_xticks(ticks)
        axis.set_xlim(d.min(), d.max())
        axis.grid(axis="x", color="black", linewidth=GUIDE_WIDTH)
        if axis is profile_ax:
            _style_axis(axis)

    # Key: alternating saddle / peak blocks on a lambda row, where the
    # original put it (right of a map 0.72 of the row wide at most).
    key_x = key_row.x0 + 0.72 * key_row.width + 0.012
    key_h = (0.5 / 2.54) / fig.get_figheight()
    key_y = key_row.y1 - key_h
    key_w = key_row.x1 - key_x
    key_ax = fig.add_axes([key_x, key_y, key_w, key_h])
    key_ax.set_xlim(0.0, 1.0)
    key_ax.set_ylim(0.0, 1.0)
    key_ax.set_xticks([])
    for spine in key_ax.spines.values():
        spine.set_visible(False)
    edges = np.linspace(0.06, 0.94, 6)
    for a, b, colour in zip(edges[:-1], edges[1:], (POSITIVE_CURVATURE, NEGATIVE_CURVATURE,
                                                    POSITIVE_CURVATURE, NEGATIVE_CURVATURE,
                                                    POSITIVE_CURVATURE)):
        key_ax.axvspan(a, b, ymin=0.0, ymax=1.0, facecolor=colour, edgecolor="none", alpha=0.28)
    for x, label in zip((edges[:3] + edges[1:4]) / 2, ("S", "P", "S")):
        key_ax.text(x, -0.20, label, transform=key_ax.transAxes, ha="center", va="top",
                    clip_on=False, color=NEGATIVE_CURVATURE)
    key_ax.plot([0.04, 1.0], [0.50, 0.50], color=NEGATIVE_CURVATURE, linewidth=GUIDE_WIDTH,
                solid_capstyle="butt")
    key_ax.set_yticks([0.50], [r"$\lambda$ value"])
    key_ax.yaxis.tick_right()
    key_ax.tick_params(axis="y", which="major", right=True, left=False, labelright=True,
                       labelleft=False, length=GUIDE_WIDTH, width=GUIDE_WIDTH,
                       colors=NEGATIVE_CURVATURE)
    key_ax.minorticks_off()
    # Offsets of the original (0.012 and 0.13 of an 11 cm figure), in cm.
    cm = 1 / 2.54 / fig.get_figheight()
    positive_text = fig.text(key_x + key_w + 0.008, key_y - 0.132 * cm,
                             "Positive\nmaximal\ncurvature:\nsaddles", ha="left", va="top",
                             color=NEGATIVE_CURVATURE)
    negative_text = fig.text(key_x + key_w + 0.008, key_y - 1.43 * cm,
                             "Negative\nmaximal\ncurvature:\npeaks", ha="left", va="top",
                             color=NEGATIVE_CURVATURE)

    fig.canvas.draw()
    width_px, height_px = fig.canvas.get_width_height()
    for axis, label in ((native_ax, "$A_{CH}$\nthreshold\n$[\\mathrm{km^2}]$"),
                        (topo_ax, "$D_{topo}$\nthreshold\n$[\\mathrm{m}]$"),
                        (profile_ax, y_axis_label)):
        fig.text(0.027, (axis.bbox.y0 + axis.bbox.y1) / (2 * height_px), label,
                 rotation=90, ha="center", va="center")
    fig.text(0.945, (profile_ax.bbox.y0 + profile_ax.bbox.y1) / (2 * height_px),
             "More topographic\nsmoothing, $\\lambda$ [m]", rotation=90, ha="center",
             va="center", color=NEGATIVE_CURVATURE)
    fig.text(native_ax.bbox.x1 / width_px + 0.01,
             (native_ax.bbox.y0 + native_ax.bbox.y1) / (2 * height_px),
             "Sparser\ndivide\nnetwork", ha="left", va="center")
    fig.text(topo_ax.bbox.x1 / width_px + 0.01,
             (topo_ax.bbox.y0 + topo_ax.bbox.y1) / (2 * height_px),
             "More\nrestrictive\nthreshold", ha="left", va="center")
    for axis, colour in ((native_ax, "black"), (topo_ax, "black"),
                         (profile_ax, NEGATIVE_CURVATURE)):
        fig.add_artist(FancyArrowPatch((0.985, axis.bbox.y1 / height_px),
                                       (0.985, axis.bbox.y0 / height_px),
                                       transform=fig.transFigure, arrowstyle="-|>",
                                       mutation_scale=6, linewidth=GUIDE_WIDTH, color=colour))
    fig.text(native_ax.bbox.x0 / width_px, native_ax.bbox.y1 / height_px + 0.008,
             "Triple divide junctions natively defined by divide network, by position "
             "along ridge-crest", ha="left", va="bottom", color=TDJ_RED)
    fig.text(topo_ax.bbox.x0 / width_px, topo_ax.bbox.y1 / height_px + 0.008,
             "Topologically defined triple divide junctions, by position along ridge-crest",
             ha="left", va="bottom", color=TDJ_RED)
    for axis in (native_ax, topo_ax, profile_ax, barcode_ax):
        axis.set_ylabel("")
    for text in fig.findobj(match=plt.Text):
        text.set_fontsize(font_size)
    fig.canvas.draw()
    to_fig = fig.transFigure.inverted()
    for text, block, radius in ((positive_text, 4, 0.25), (negative_text, 3, -0.25)):
        box = text.get_window_extent()
        start = to_fig.transform((box.x0, box.y1))
        target = to_fig.transform(key_ax.transAxes.transform(
            ((edges[block] + edges[block + 1]) / 2, 0.0)))
        fig.add_artist(FancyArrowPatch(start, target, transform=fig.transFigure,
                                       arrowstyle="-|>", connectionstyle=f"arc3,rad={radius}",
                                       mutation_scale=6, linewidth=GUIDE_WIDTH,
                                       color=NEGATIVE_CURVATURE))
    return fig, {"native": native_ax, "topological": topo_ax, "profile": profile_ax,
                 "barcode": barcode_ax, "key": key_ax}


# -----------------------------------------------------------------------------
# Main drainage divide: map image + profile
# -----------------------------------------------------------------------------
def _nice_amplitude(values) -> tuple:
    amplitude = max(float(np.nanmax(np.abs(values))) * 1.05, 1e-9)
    tick = 10.0 ** np.floor(np.log10(amplitude))
    tick *= 1 if amplitude / tick < 2 else 2 if amplitude / tick < 5 else 5
    return amplitude, float(tick)


def plot_mdd_ridgeline_map_profile(profile: pd.DataFrame, image: np.ndarray,
                                   grid_spacing: float, model_name: str, *,
                                   guide_interval: float = 200.0,
                                   curvature_column: Optional[str] = None):
    """A main-drainage-divide profile above the model's drainage-area image.

    The profile (detrended elevation, with a maximal-curvature sign barcode and
    the median zero-crossing spacing) and the image share vertical guides every
    ``guide_interval`` metres of profile distance. The ridge is drawn on the
    image from its divide-grid ``row``/``col``, shifted half a cell onto pixel
    centres. Meant for an east-west divide across a periodic model domain.

    Args:
        profile: One ridgeline sorted by distance, with ``row``, ``col``,
            ``distance``, ``detrended_elevation``, ``fitted_elevation`` and a
            maximal-curvature column.
        image: The drainage-area image (rows x cols [x bands]), e.g. the
            ``lsdtt-lem`` log10 drainage-area PNG.
        grid_spacing: Model cell size (m), for the domain-height label
            (:func:`~lsdviztools3.io.ridgeline.divide_grid_spacing`).
        model_name: Label at the top left.
        guide_interval: Distance between guides (m).
        curvature_column: Maximal-curvature column (default: the only one, or
            the smallest lambda).

    Returns:
        ``(fig, (profile_ax, map_ax))``.
    """
    plt = _import_mpl()
    from matplotlib import patheffects
    from matplotlib.lines import Line2D

    for column in ("row", "col", "distance", "detrended_elevation", "fitted_elevation"):
        if column not in profile.columns:
            raise ValueError(f"Profile is missing column {column!r}")
    if curvature_column is None:
        scales = lambda_columns(profile.columns)
        if not scales:
            raise ValueError("Profile needs a maximal-curvature column; enable a polyfit metric")
        curvature_column = scales[0][1]
    image = np.asarray(image)
    h, w = image.shape[:2]
    rows = profile["row"].to_numpy(dtype=float)
    cols = profile["col"].to_numpy(dtype=float)
    if rows.min() < 0 or rows.max() > h or cols.min() < 0 or cols.max() > w:
        raise ValueError("Ridgeline row/col fall outside the divide grid implied by the image")
    map_cols, map_rows = cols - 0.5, rows - 0.5
    d = profile["distance"].to_numpy(dtype=float)
    detrended = profile["detrended_elevation"].to_numpy(dtype=float)
    curvature = profile[curvature_column].to_numpy(dtype=float)
    x_max = float(np.nanmax(d))
    ticks = np.arange(0.0, x_max + 1.0, guide_interval)
    if not ticks.size:
        ticks = np.array([0.0])
    amplitude, y_tick = _nice_amplitude(detrended)
    spacing = periodic_median_crossing_spacing(d, curvature)

    width_in = 7.0 / 2.54
    left_in, right_in = 0.19 * width_in, 0.02 * width_in
    axis_w = width_in - left_in - right_in
    prof_h, gap, bottom, top = 0.60, 0.10, 0.45, 0.35
    map_h = axis_w * h / w
    fig_h = bottom + map_h + gap + prof_h + top
    fig = plt.figure(figsize=(width_in, fig_h))
    prof_ax = fig.add_axes([left_in / width_in, (bottom + map_h + gap) / fig_h,
                            axis_w / width_in, prof_h / fig_h])
    map_ax = fig.add_axes([left_in / width_in, bottom / fig_h, axis_w / width_in, map_h / fig_h])

    prof_ax.axhline(0, color="black", linewidth=0.5, zorder=0)
    prof_ax.plot(d, detrended, color="black", linewidth=0.6, zorder=2)
    prof_ax.plot(d, detrended, color=RIDGE_GREEN, linewidth=2.4, alpha=0.35, zorder=1)
    barcode = add_curvature_barcode(prof_ax, d, curvature)
    barcode.set_yticks([])
    for spine in barcode.spines.values():
        spine.set_visible(False)
    prof_ax.set_xlim(0, x_max)
    prof_ax.set_xticks(ticks)
    prof_ax.set_ylim(-amplitude, amplitude)
    prof_ax.set_yticks((-y_tick, 0, y_tick))
    prof_ax.set_yticklabels([f"-{y_tick:g} m",
                             f"{profile['fitted_elevation'].iloc[0]:.1f}", f"+{y_tick:g} m"])
    prof_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    prof_ax.tick_params(axis="y", labelsize=6, width=0.5)
    label = f"Median spacing: {spacing:.1f} m" if np.isfinite(spacing) else "Median spacing: n/a"
    prof_ax.text(0.02, 0.05, label, transform=prof_ax.transAxes, ha="left", va="bottom",
                 fontsize=6, path_effects=[patheffects.Stroke(linewidth=1.5, foreground="white"),
                                           patheffects.Normal()])
    for t in ticks:
        prof_ax.axvline(t, color="black", linewidth=0.5, alpha=0.75, zorder=0)
    for spine in prof_ax.spines.values():
        spine.set_visible(spine.spine_type == "left")
        spine.set_linewidth(0.5)

    map_ax.imshow(image, cmap="Greys_r", aspect="auto")
    map_ax.set_xlim(0, w - 1)
    map_ax.set_ylim(h - 1, 0)
    map_ax.plot(map_cols, map_rows, color="black", linewidth=0.6, zorder=2)
    map_ax.plot(map_cols, map_rows, color=RIDGE_GREEN, linewidth=2.4, alpha=0.35, zorder=3)
    map_x = {float(t): float(np.interp(t, d, map_cols)) for t in ticks}
    for x in map_x.values():
        map_ax.axvline(x, color="black", linewidth=0.5, alpha=0.75, zorder=1)
    map_ax.yaxis.set_visible(False)
    map_ax.set_xticks([map_x[float(t)] for t in ticks])
    labels = [f"{t:g}" for t in ticks]
    labels[-1] = f"{ticks[-1]:g} m"
    map_ax.set_xticklabels(labels)
    map_ax.tick_params(axis="x", labelsize=6, width=0.5, pad=2)
    map_ax.text(-0.10, 0.5, f"{h * grid_spacing:g} m", rotation=90, transform=map_ax.transAxes,
                ha="center", va="center", fontsize=6,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 2.0}, clip_on=False)
    map_ax.annotate("", xy=(-0.055, 1.03), xytext=(-0.055, -0.03), xycoords=map_ax.transAxes,
                    arrowprops={"arrowstyle": "<->", "color": "black", "linewidth": 0.5},
                    annotation_clip=False)
    for spine in map_ax.spines.values():
        spine.set_visible(spine.spine_type in ("left", "top", "bottom"))
        spine.set_linewidth(0.5)
    fig.patch.set_facecolor("white")
    fig.text(0.02, 0.98, model_name, ha="left", va="top", fontsize=7)
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    for t in ticks:
        px = inv.transform((prof_ax.transData.transform((t, 0))[0], 0))[0]
        mx = inv.transform((map_ax.transData.transform((map_x[float(t)], 0))[0], 0))[0]
        fig.lines.append(Line2D((px, mx), (prof_ax.get_position().y0, map_ax.get_position().y1),
                                transform=fig.transFigure, color="black", linewidth=0.5,
                                alpha=0.75, zorder=0))
    return fig, (prof_ax, map_ax)


# -----------------------------------------------------------------------------
# Zero-crossing spacing and slope-area
# -----------------------------------------------------------------------------
def _stepped_edges(values, logarithmic: bool = False) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 1:
        width = max(abs(values[0]) * 0.05, 1.0)
        return np.asarray([values[0] - width, values[0] + width])
    if logarithmic:
        inner = np.sqrt(values[:-1] * values[1:])
        first, last = values[0] ** 2 / inner[0], values[-1] ** 2 / inner[-1]
    else:
        inner = 0.5 * (values[:-1] + values[1:])
        first, last = values[0] - (inner[0] - values[0]), values[-1] + (values[-1] - inner[-1])
    return np.concatenate(([first], inner, [last]))


def _fill_quantiles(ax, x, low, high, colour, logarithmic: bool = False) -> None:
    if not np.size(x):
        return
    edges = _stepped_edges(x, logarithmic)
    ax.fill_between(edges, np.r_[low[0], low], np.r_[high[0], high], step="post",
                    color=colour, linewidth=0)


def plot_zero_spacing_slope_area(sites: Sequence[dict], *, log_spacing: bool = False,
                                 fixed_spacing_ylim: bool = True, spacing_y_max: float = 200.0,
                                 absolute_curvature_label: bool = False,
                                 figure_width: float = 2.75, row_height: float = 1.1,
                                 xlabel_space: float = 0.45, font_size: float = 6.0):
    """Zero-crossing spacing against polyfit scale beside slope-area, per site.

    Left column: 5-95 % and 25-75 % bands and the median (+) of the spacing of
    maximal-curvature sign changes against lambda, with a dashed line at 2L.
    Right column: the slope-area distribution (bands and median, log area),
    with the area of the largest median slope and L = sqrt((A / dx)^2).

    Args:
        sites: ``[{"name", "spacing", "slope_area"}, ...]`` with ``spacing`` from
            :func:`~lsdviztools3.io.ridgeline.curvature_spacing_quantiles` and
            ``slope_area`` from :func:`~lsdviztools3.io.ridgeline.load_slope_area`.
        log_spacing: Logarithmic spacing axis.
        fixed_spacing_ylim: Spacing axis 0..``spacing_y_max`` (ignored with
            ``log_spacing``).
        spacing_y_max: Upper spacing limit (m).
        absolute_curvature_label: Label the left axis as absolute curvature.
        figure_width: Figure width (inches), independent of the site count.
        row_height: Height per site (inches).
        xlabel_space: Extra height for the bottom x labels (inches).
        font_size: Font size.

    Returns:
        ``(fig, axes)`` with ``axes`` of shape ``(n_sites, 2)``.
    """
    plt = _import_mpl()
    if not sites:
        raise ValueError("At least one site is required")
    if spacing_y_max <= 0 or figure_width <= 0 or row_height <= 0 or xlabel_space < 0:
        raise ValueError("Spacing limit, figure width and row height must be positive")
    n = len(sites)
    fig_h = row_height * n + xlabel_space
    fig, axes = plt.subplots(n, 2, figsize=(figure_width, fig_h), squeeze=False)
    labels = []
    for i, site in enumerate(sites):
        spacing_ax, slope_ax = axes[i]
        spacing = site["spacing"]
        sa = site["slope_area"]
        site_label = spacing_ax.text(0.03, 0.95, site["name"], transform=spacing_ax.transAxes,
                                     ha="left", va="top", fontsize=font_size, color=TDJ_RED)
        x = spacing["polyfit_lambda"]
        if x.size:
            _fill_quantiles(spacing_ax, x, spacing["d05"], spacing["d95"], "#bbbbbb")
            _fill_quantiles(spacing_ax, x, spacing["d25"], spacing["d75"], "#777777")
            spacing_ax.plot(x, spacing["d50"], linestyle="None", marker="+", color="black",
                            markersize=2, markeredgewidth=0.35)
        positive = sa["drainage_area"].to_numpy(dtype=float) > 0
        if not np.any(positive):
            raise ValueError(f"Slope-area data for {site['name']} has no positive drainage areas")
        area = sa["drainage_area"].to_numpy(dtype=float)[positive]
        q = {c: sa[c].to_numpy(dtype=float)[positive]
             for c in ("d05_slope", "d25_slope", "d50_slope", "d75_slope", "d95_slope")}
        _fill_quantiles(slope_ax, area, q["d05_slope"], q["d95_slope"], "#bbbbbb", True)
        _fill_quantiles(slope_ax, area, q["d25_slope"], q["d75_slope"], "#777777", True)
        slope_ax.plot(area, q["d50_slope"], linestyle="None", marker="+", color="black",
                      markersize=2, markeredgewidth=0.35)
        peak = slope_area_peak(sa)
        length = characteristic_length(peak, float(spacing["grid_resolution"]))
        slope_ax.axvline(peak, color="green", linewidth=1.0)
        slope_ax.text(0.97, 0.95, f"L $\\approx$ {length:.0f} m", transform=slope_ax.transAxes,
                      ha="right", va="top", fontsize=font_size, color="green", fontweight="bold")
        if x.size:
            spacing_ax.text(1.02 * length, 0.97 * spacing_y_max, "2 x L", color="green",
                            fontsize=font_size, ha="left", va="top", fontweight="bold")
            line = spacing_ax.axvline(2 * length, color="green", linewidth=1.0, linestyle="--")
            labels.append((spacing_ax, line, site_label))
        slope_ax.set_xscale("log")
        if log_spacing:
            spacing_ax.set_yscale("log")
        elif fixed_spacing_ylim:
            spacing_ax.set_ylim(0, spacing_y_max)
        else:
            spacing_ax.autoscale(axis="y", tight=True)
        for axis in (spacing_ax, slope_ax):
            axis.tick_params(direction="in", labelsize=font_size)
            axis.margins(x=0, y=0)
            for spine in axis.spines.values():
                spine.set_linewidth(0.6)
        slope_ax.yaxis.tick_right()
        slope_ax.yaxis.set_label_position("right")
        if i < n - 1:
            spacing_ax.tick_params(labelbottom=False)
            slope_ax.tick_params(labelbottom=False)
    axes[-1, 0].set_xlabel("Degree of topographic\nsmoothing, $\\lambda$ [m]", fontsize=font_size)
    axes[-1, 1].set_xlabel("D$\\infty$ drainage area [m$^2$]", fontsize=font_size)
    fig.subplots_adjust(left=0.16, right=0.84, top=0.98, bottom=xlabel_space / fig_h,
                        wspace=0.22, hspace=0.14)
    fig.text(0.055, 0.5, "Absolute maximal curvature [m$^{-1}$]" if absolute_curvature_label
             else "Peak-and-saddle spacing from maximal curvature sign flips [m]",
             rotation=90, ha="center", va="center", fontsize=font_size)
    fig.text(0.975, 0.5, "Pixel slope [-]", rotation=90, ha="center", va="center",
             fontsize=font_size)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for axis, line, site_label in labels:
        bottom_px = site_label.get_window_extent(renderer=renderer).y0
        bottom = axis.transAxes.inverted().transform((0, bottom_px))[1]
        line.set_ydata((0.0, max(0.0, bottom - 0.015)))
    return fig, axes
