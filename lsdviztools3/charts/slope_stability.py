"""Charts and maps of an ``lsdtt-slope-stability`` run.

Everything drawn is read from the run's files (:mod:`lsdviztools3.io.slope_stability`);
nothing is recomputed. Histogram bin edges and point subsampling are drawing
choices, not results.

* SINMAP's slope - specific catchment area plot: the cells coloured by
  stability class, with the class boundaries lsdtt3 wrote
  (``<prefix>_sinmap_slope_area_curves.csv``).
* Catalogue charts for the rupture-plane landslides and the slip bodies: area
  against volume and against mean depth (log-log), the distribution of area,
  and the distribution of length / width.
* The rupture-plane depth against the slip surface's deepest body thickness,
  per crown: how much deeper the integrated or curved surfaces go.
* Maps on a hillshade: the SINMAP class map; rupture depth; probabilities of
  failure; slip-surface factor of safety with the screened-out crowns shown;
  slip-body depth.

Pack, R. T., Tarboton, D. G., Goodwin, C. N. and Prasad, A. (2005). SINMAP 2,
user's manual. Utah State University.
Steer, P., Pelascini, L. and Longuevergne, L. (2026). EGUsphere preprint,
https://doi.org/10.5194/egusphere-2026-4142.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np

from .slope_area import SERIES_COLOURS

_MPL_MSG = "The slope-stability charts need matplotlib. Install it e.g. `pip install matplotlib`."

#: SINMAP classes 1 (stable) to 6 (defended): an ordinal blue-to-red ramp
#: (ColorBrewer RdBu, its six strongest steps), readable with red-green colour
#: blindness and with no class near white.
SINMAP_CLASS_COLOURS = {1: "#2166ac", 2: "#4393c3", 3: "#92c5de", 4: "#f4a582", 5: "#d6604d",
                        6: "#b2182b"}

#: Line styles of SINMAP's class boundaries, all in dark ink.
_CURVE_STYLE = {
    "always_saturated_above": ("-", "always saturated above"),
    "never_saturated_below": (":", "never saturated below"),
    "fs_min_1.5": ("--", "FS_min = 1.5"),
    "fs_min_1.0": ("-.", "FS_min = 1.0"),
    "fs_max_1.0": ((0, (5, 1, 1, 1, 1, 1)), "FS_max = 1.0"),
}

#: Colour maps of the maps.
SAFETY_FACTOR_CMAP = "RdBu"       # diverging about F = 1: red unstable, blue stable
PROBABILITY_CMAP = "Reds"         # sequential, 0 to 1
DEPTH_CMAP = "YlOrBr"             # sequential, light (shallow) to dark (deep)
#: Slip crowns not evaluated. Slightly off pure grey: GMT paints the NaN cells of
#: a pure-grey CPT red despite nan_transparent (see render.colormaps.resolve_cmap).
SCREENED_COLOUR = "150/150/152"


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def _finite_positive(*arrays):
    keep = np.ones(arrays[0].shape, dtype=bool)
    for a in arrays:
        keep &= np.isfinite(a) & (a > 0)
    return [a[keep] for a in arrays]


def _subsample(n: int, max_points: int) -> np.ndarray:
    """Every k-th index, deterministic, so at most max_points are drawn."""
    if n <= max_points:
        return np.arange(n)
    return np.arange(0, n, int(np.ceil(n / max_points)))


# -----------------------------------------------------------------------------
# SINMAP
# -----------------------------------------------------------------------------
def plot_sinmap_slope_area(slope, sca, stability_class, curves, *, region_id: int = 0,
                           max_points: int = 200_000, title: Optional[str] = None, ax=None):
    """SINMAP's slope - specific catchment area plot.

    Args:
        slope: D-infinity slope raster values (m/m), any shape.
        sca: Specific catchment area raster values (m), same shape.
        stability_class: SINMAP class raster values (1-6), same shape.
        curves: The ``_sinmap_slope_area_curves.csv`` table.
        region_id: Which region's curves to draw.
        max_points: At most this many cells are drawn (every k-th cell).

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    s, a, c = (np.asarray(x, dtype=float).ravel() for x in (slope, sca, stability_class))
    keep = np.isfinite(s) & np.isfinite(a) & np.isfinite(c) & (a > 0) & (s > 0)
    s, a, c = s[keep], a[keep], c[keep].astype(int)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
    else:
        fig = ax.figure
    idx = _subsample(s.size, max_points)
    for k in sorted(SINMAP_CLASS_COLOURS):
        sel = idx[c[idx] == k]
        if sel.size == 0:
            continue
        ax.scatter(a[sel], s[sel], s=2, color=SINMAP_CLASS_COLOURS[k], alpha=0.5, linewidths=0,
                   label=f"{k} {_class_name(k)} ({int((c == k).sum()):,} cells)", rasterized=True)
    reg = curves[curves["region_id"] == region_id]
    for name, (ls, label) in _CURVE_STYLE.items():
        part = reg[reg["curve"] == name]
        if part.empty:
            continue
        ax.plot(part["specific_catchment_area_m"], part["slope"], ls=ls, color="#222222", lw=1.3,
                label=label)
    ax.set_xscale("log")
    ax.set_xlabel("specific catchment area a (m)")
    ax.set_ylabel("slope (m/m)")
    ax.set_ylim(bottom=0)
    ax.grid(True, which="major", color="#dddddd", lw=0.6)
    ax.set_axisbelow(True)
    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, markerscale=4,
                    framealpha=0.9, borderaxespad=0.0)
    for h in leg.legend_handles:
        if hasattr(h, "set_alpha"):
            h.set_alpha(1.0)
    ax.set_title(title or f"SINMAP slope - specific catchment area (region {region_id} boundaries)")
    return fig, ax


def _class_name(k: int) -> str:
    from ..io.slope_stability import SINMAP_CLASSES
    return SINMAP_CLASSES.get(int(k), str(k))


# -----------------------------------------------------------------------------
# Catalogues
# -----------------------------------------------------------------------------
def _catalogue_axes(ax, figsize=(6.0, 4.8)):
    plt = _import_mpl()
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    return fig, ax


def plot_catalogue_area_volume(catalogues: Sequence[Tuple[str, "pd.DataFrame"]], *,
                               title: Optional[str] = None, ax=None):
    """Volume against area, log-log, one series per (label, catalogue)."""
    fig, ax = _catalogue_axes(ax)
    for k, (label, df) in enumerate(catalogues):
        a, v = _finite_positive(df["area_m2"].to_numpy(float), df["volume_m3"].to_numpy(float))
        ax.scatter(a, v, s=18 if k == 0 else 7, alpha=0.6, color=SERIES_COLOURS[k % len(SERIES_COLOURS)],
                   linewidths=0, label=f"{label} ({a.size:,})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("area (m²)")
    ax.set_ylabel("volume (m³)")
    ax.grid(True, which="major", color="#dddddd", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8)
    ax.set_title(title or "area - volume")
    return fig, ax


def plot_catalogue_area_depth(catalogues, *, title: Optional[str] = None, ax=None):
    """Mean depth against area, log-log, one series per (label, catalogue)."""
    fig, ax = _catalogue_axes(ax)
    for k, (label, df) in enumerate(catalogues):
        a, d = _finite_positive(df["area_m2"].to_numpy(float), df["mean_depth_m"].to_numpy(float))
        ax.scatter(a, d, s=18 if k == 0 else 7, alpha=0.6, color=SERIES_COLOURS[k % len(SERIES_COLOURS)],
                   linewidths=0, label=f"{label} ({a.size:,})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("area (m²)")
    ax.set_ylabel("mean depth (m)")
    ax.grid(True, which="major", color="#dddddd", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8)
    ax.set_title(title or "area - mean depth")
    return fig, ax


def _log_edges(values, per_decade: int = 8):
    lo, hi = np.log10(values.min()), np.log10(values.max())
    if hi - lo < 1e-9:
        lo, hi = lo - 0.5, hi + 0.5
    n = max(1, int(np.ceil((hi - lo) * per_decade)))
    return np.logspace(lo, hi, n + 1)


def plot_catalogue_area_distribution(catalogues, *, per_decade: int = 8, title: Optional[str] = None,
                                     ax=None):
    """Counts of area in logarithmic bins, one step line per catalogue (log-log)."""
    fig, ax = _catalogue_axes(ax)
    allv = np.concatenate([_finite_positive(df["area_m2"].to_numpy(float))[0] for _, df in catalogues])
    if allv.size == 0:
        raise ValueError("no positive areas to draw")
    edges = _log_edges(allv, per_decade)
    for k, (label, df) in enumerate(catalogues):
        a = _finite_positive(df["area_m2"].to_numpy(float))[0]
        counts, _ = np.histogram(a, bins=edges)
        ax.stairs(np.where(counts > 0, counts, np.nan), edges, lw=1.8,
                  color=SERIES_COLOURS[k % len(SERIES_COLOURS)], label=f"{label} ({a.size:,})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("area (m²)")
    ax.set_ylabel(f"number per bin ({per_decade} bins per decade)")
    ax.grid(True, which="major", color="#dddddd", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8)
    ax.set_title(title or "area distribution")
    return fig, ax


#: Readable names of the catalogues in lsdtt3's frequency - area tables.
CATALOGUE_LABELS = {"by_daylight": "by daylight point", "by_toe": "by toe", "amalgamated": "amalgamated"}


def plot_frequency_area(density, fit, *, title: Optional[str] = None, ax=None):
    """Frequency - area: the log-binned probability density p(A) with lsdtt3's power-law tail.

    ``density`` and ``fit`` are the tables of
    :func:`lsdviztools3.io.slope_stability.load_frequency_area`. For each
    catalogue the binned density is drawn as points, and the fitted tail as a
    line from A_min: p(A) = (n_tail / n) (alpha - 1) / A_min (A / A_min)^(-alpha),
    alpha = -eta, the maximum-likelihood fit lsdtt3 wrote (its normalisation over
    the whole catalogue). The legend gives eta with its standard error and A_min;
    a dotted line marks A_min.
    """
    fig, ax = _catalogue_axes(ax, figsize=(6.2, 5.0))
    for k, cat in enumerate(fit["catalogue"].tolist()):
        colour = SERIES_COLOURS[k % len(SERIES_COLOURS)]
        part = density[(density["catalogue"] == cat) & (density["n_landslides"] > 0)]
        label = CATALOGUE_LABELS.get(cat, cat)
        ax.scatter(part["area_centre_m2"], part["probability_density_per_m2"], s=22, color=colour, linewidths=0,
                   alpha=0.85, label=f"{label} ({int(fit.loc[fit['catalogue'] == cat, 'n_landslides'].iloc[0]):,})")
        row = fit[fit["catalogue"] == cat].iloc[0]
        eta, se, amin = float(row["eta"]), float(row["eta_standard_error"]), float(row["area_min_m2"])
        if not (np.isfinite(eta) and np.isfinite(amin) and amin > 0):
            continue
        alpha = -eta
        frac = float(row["n_tail"]) / float(row["n_landslides"])
        amax = float(part["area_upper_m2"].max()) if "area_upper_m2" in part else float(part["area_centre_m2"].max())
        a = np.logspace(np.log10(amin), np.log10(max(amax, amin * 1.01)), 50)
        p = frac * (alpha - 1.0) / amin * (a / amin) ** (-alpha)
        ax.plot(a, p, color=colour, lw=1.8,
                label=f"  \u03b7 = {eta:.2f} \u00b1 {se:.2f}, A \u2265 {amin:,.0f} m\u00b2 ({int(row['n_tail']):,})")
        ax.axvline(amin, color=colour, lw=0.9, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("area A (m²)")
    ax.set_ylabel("probability density p(A) (m⁻²)")
    ax.grid(True, which="major", color="#dddddd", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, loc="lower left")
    ax.set_title(title or "frequency - area")
    return fig, ax


def plot_catalogue_length_width(catalogues, *, per_decade: int = 10, title: Optional[str] = None,
                                ax=None):
    """Counts of length / width (logarithmic bins from 1), one step line per catalogue."""
    fig, ax = _catalogue_axes(ax)
    allv = np.concatenate([_finite_positive(df["length_to_width"].to_numpy(float))[0] for _, df in catalogues])
    if allv.size == 0:
        raise ValueError("no finite length / width ratios to draw")
    edges = _log_edges(np.maximum(allv, 1.0), per_decade)
    for k, (label, df) in enumerate(catalogues):
        r = _finite_positive(df["length_to_width"].to_numpy(float))[0]
        counts, _ = np.histogram(np.maximum(r, 1.0), bins=edges)
        ax.stairs(counts, edges, lw=1.8, color=SERIES_COLOURS[k % len(SERIES_COLOURS)],
                  label=f"{label} ({r.size:,})")
    ax.set_xscale("log")
    ax.set_xlabel("length / width")
    ax.set_ylabel("number per bin")
    ax.grid(True, which="major", color="#dddddd", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8)
    ax.set_title(title or "length / width")
    return fig, ax


def plot_depth_comparison(rupture_depth, slip_max_depth, *, slip_label: str = "slip surface",
                          max_points: int = 100_000, title: Optional[str] = None, ax=None):
    """Per crown: the slip surface's deepest body thickness against the rupture-plane depth.

    Both rasters are on the same grid; only cells with both values are drawn,
    with the 1:1 line. Log-log.
    """
    plt = _import_mpl()
    r = np.asarray(rupture_depth, dtype=float).ravel()
    s = np.asarray(slip_max_depth, dtype=float).ravel()
    r, s = _finite_positive(r, s)
    if r.size == 0:
        raise ValueError("no crown has both a rupture depth and a slip depth")
    if ax is None:
        fig, ax = plt.subplots(figsize=(5.4, 5.0))
    else:
        fig = ax.figure
    idx = _subsample(r.size, max_points)
    ax.scatter(r[idx], s[idx], s=5, alpha=0.4, color=SERIES_COLOURS[0], linewidths=0, rasterized=True,
               label=f"crowns ({r.size:,})")
    lo, hi = min(r.min(), s.min()), max(r.max(), s.max())
    ax.plot([lo, hi], [lo, hi], color="#222222", lw=1.0, ls="--", label="1:1")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("rupture-plane depth h (m)")
    ax.set_ylabel(f"{slip_label}: deepest body thickness (m)")
    ax.grid(True, which="major", color="#dddddd", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8)
    ax.set_title(title or "depth per crown: slip surface against rupture plane")
    return fig, ax


# -----------------------------------------------------------------------------
# Maps (pygmt)
# -----------------------------------------------------------------------------
def _map_style():
    from ..render.style import MapStyle
    return MapStyle(figure_size=MapStyle.MAX_WIDTH, transparent=False)


def _finish_map(fig, dem_path, png_path, coordinates, title, extra_rows=(), legend_position="JTR+jTR+o0.2c",
                high_res=False):
    import io as _io

    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from ..render.base import grid_to_display, save_figure
    from ..render.decorations import add_scalebar

    da = to_coordinates(load_raster(dem_path), coordinates)
    crs = get_crs(da)
    _, extent = grid_to_display(to_gmt_grid(da), crs)
    rows = []
    if title:
        rows.append(f"L 9p - L {title}")
    rows.extend(extra_rows)
    if rows:
        fig.legend(spec=_io.StringIO("\n".join(rows)), position=legend_position,
                   box="+gwhite@20+p0.5p,black")
    add_scalebar(fig, extent, crs)
    save_figure(fig, png_path, style=_map_style(), high_res=high_res)
    return Path(png_path)


def render_sinmap_class_map(dem_path, class_raster, png_path, *, class_areas=None,
                            coordinates: str = "local", high_res: bool = False):
    """The SINMAP stability classes on a hillshade, with each class's share of the area."""
    from ..render.categorical import render_categorical_drape

    shares = {}
    if class_areas is not None and not class_areas.empty:
        tot = class_areas.groupby("stability_class")["area_m2"].sum()
        frac = tot / tot.sum()
        shares = {int(k): float(v) for k, v in frac.items()}
    categories = {k: (f"{k} {_class_name(k)}" + (f" ({100 * shares[k]:.1f} %)" if k in shares else ""),
                      SINMAP_CLASS_COLOURS[k]) for k in SINMAP_CLASS_COLOURS}
    fig = render_categorical_drape(dem_path, class_raster, categories=categories, coordinates=coordinates,
                                   style=_map_style(), transparency=25.0, legend_header="SINMAP stability class",
                                   require_labels=False, decorate=False)
    return _finish_map(fig, dem_path, png_path, coordinates, None, high_res=high_res)


def render_value_map(dem_path, value_raster, png_path, *, cmap: str, label: str,
                     vmin: Optional[float] = None, vmax: Optional[float] = None,
                     title: Optional[str] = None, overlay_mask=None, overlay_label: Optional[str] = None,
                     coordinates: str = "local", transparency: float = 25.0, high_res: bool = False):
    """A continuous raster on a hillshade, with a colour bar.

    ``overlay_mask`` (a raster: non-zero = drawn) is painted grey on top, e.g.
    the slip-surface crowns that were screened out and not evaluated.
    """
    from ..render.channels import render_drape

    fig = render_drape(dem_path, value_raster, coordinates=coordinates, base_is_hillshade=False,
                       cmap=cmap, vmin=vmin, vmax=vmax, transparency=transparency, legend=False,
                       colorbar=True, colorbar_label=label, style=_map_style())
    rows = []
    if overlay_mask is not None:
        n = _overlay_mask(fig, overlay_mask, coordinates)
        rows.append(f"S 0.3c s 0.35c {SCREENED_COLOUR} - 0.9c {overlay_label or 'masked'} ({n:,} cells)")
    return _finish_map(fig, dem_path, png_path, coordinates, title, extra_rows=rows, high_res=high_res)


def _overlay_mask(fig, mask_raster, coordinates) -> int:
    """Paint the non-zero cells of mask_raster grey over fig; returns their count."""
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from ..render.base import auto_projection, grid_to_display, set_grid_gtype

    mask = mask_raster if hasattr(mask_raster, "rio") else load_raster(mask_raster)
    mask = to_gmt_grid(to_coordinates(mask, coordinates))
    crs = get_crs(mask)
    mask, extent = grid_to_display(mask, crs)
    set_grid_gtype(mask, crs)
    values = np.asarray(mask.values, dtype=float)
    drape = np.where(np.isfinite(values) & (values != 0), 1.0, np.nan)
    n = int(np.isfinite(drape).sum())
    if n == 0:
        return 0
    mask = mask.copy(data=drape)
    style = _map_style()
    projection = style.projection or auto_projection(extent, crs, style.figure_size)
    pygmt.makecpt(cmap=f"{SCREENED_COLOUR},{SCREENED_COLOUR}", series=[0.5, 1.5, 1.0])
    fig.grdimage(grid=mask, region=extent, projection=projection, cmap=True, shading=False,
                 nan_transparent=True, transparency=55)
    return n


def screened_mask(evaluated_raster):
    """In memory: 1 where a slip crown was screened out (``_slip_evaluated`` == 0), else NaN."""
    from ..io.raster import load_raster

    da = evaluated_raster if hasattr(evaluated_raster, "rio") else load_raster(evaluated_raster)
    v = np.asarray(da.values, dtype=float)
    return da.copy(data=np.where(v == 0.0, 1.0, np.nan))
