"""Perimeter valley sources: where big valleys cross the edge of a DEM.

Reads the outputs of one lsdtt3 ``lsdtt-channel-extraction`` run made with
``run_perimeter_sources: true``, by its stem (write_path + write_fname):

  <stem>_perimeter_crossings.fgb   every gentle (flat) run along the perimeter,
                                   as a line, with its verdict
  <stem>_perimeter_outlets.fgb     one point per outlet run, ranked by summed
                                   drainage area, with ``selected``
  <stem>_perimeter_sources.fgb     the source points
  <stem>_perimeter_profile.csv     one row per perimeter visit

Two figures:

``plot_perimeter_sources_map``: the hillshade with the traced perimeter, the
flat runs along it coloured by verdict, the outlets (the selected ones filled,
the next largest hollow, numbered by rank) and the sources.

``plot_perimeter_profile``: elevation round the perimeter from its highest
cell, flat runs coloured by verdict, outlets and sources marked, and the
perimeter gradient in a panel below with the max_gradient threshold.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Fixed order, most interesting first. Rejected runs are dashed so identity
# does not rest on colour alone; the two size rejections are grey.
VERDICT_COLOURS = {
    "source": "#1baf7a",
    "outlet_crossing": "#2a78d6",
    "failed_inward": "#eb6834",
    "drains_elsewhere": "#8e5bd0",
    "too_close_to_outlet": "#d4a20a",
    "not_low_spot": "#b4b2a9",
    "too_narrow": "#8a8a85",
}
VERDICT_LABELS = {
    "source": "source crossing",
    "outlet_crossing": "outlet crossing (flow leaves)",
    "failed_inward": "rejected: too steep inward",
    "drains_elsewhere": "rejected: drains to an unselected outlet",
    "too_close_to_outlet": "rejected: too close to the outlet",
    "not_low_spot": "rejected: not a low spot",
    "too_narrow": "rejected: narrower than min width",
}
PERIMETER_COLOUR = "#3d3d3a"
OUTLET_COLOUR = "#2a78d6"
SOURCE_COLOUR = "#1baf7a"


def _read(stem, suffix):
    p = Path(f"{stem}{suffix}")
    if not p.exists():
        raise FileNotFoundError(
            f"{p} (run lsdtt-channel-extraction with run_perimeter_sources: true).")
    return p


def _ingested_value(stem, key, default):
    """A parameter from <stem>_ingested_params.param, else ``default``."""
    p = Path(f"{stem}_ingested_params.param")
    if p.exists():
        for line in p.read_text().splitlines():
            k, _, v = line.partition(":")
            if k.strip() == key:
                try:
                    return float(v.strip())
                except ValueError:
                    break
    return default


def _hillshade_from(hillshade=None, dem=None, max_pixels: int = 3000):
    """(image, extent) from a hillshade raster, or computed from a DEM.

    Large rasters are read decimated to at most ``max_pixels`` on a side; the
    hillshade is only a backdrop.
    """
    import rasterio
    from rasterio.plot import plotting_extent

    path = hillshade if hillshade is not None else dem
    if path is None:
        raise ValueError("give a hillshade or a DEM")
    with rasterio.open(str(path)) as src:
        step = max(1, int(np.ceil(max(src.width, src.height) / max_pixels)))
        shape = (max(1, src.height // step), max(1, src.width // step))
        z = src.read(1, out_shape=shape, masked=True).astype("float64").filled(np.nan)
        extent = plotting_extent(src)
        cell = abs(src.transform.a) * src.width / shape[1]
    if hillshade is None:
        from .network_maps import hillshade as _hs
        z = np.where(np.isfinite(z), _hs(np.nan_to_num(z, nan=np.nanmean(z)),
                                         cell_size=cell), np.nan)
    return z, extent


def _ranked_outlets(outlets, n_unselected: int):
    sel = outlets[outlets["selected"] == 1]
    rest = outlets[outlets["selected"] != 1].sort_values("rank").head(n_unselected)
    return sel, rest


def plot_perimeter_sources_map(
    stem,
    out_png,
    hillshade=None,
    dem=None,
    bbox=None,
    n_unselected_outlets: int = 3,
    title: str | None = None,
    figsize=(10.0, 11.0),
    dpi: int = 200,
) -> Path:
    """Flat perimeter runs by verdict, outlets and sources on the hillshade.

    Args:
        stem: lsdtt-channel-extraction output stem.
        out_png: Output PNG.
        hillshade: Hillshade raster on the run's grid; or give ``dem``.
        dem: DEM to compute a hillshade from when no hillshade is given.
        bbox: Optional (xmin, ymin, xmax, ymax); default the raster extent.
        n_unselected_outlets: How many of the largest unselected outlet runs
            to draw (hollow), besides the selected ones.
        title: Figure title.
        figsize: Figure size (inches).
        dpi: Output resolution.

    Returns:
        The path written.
    """
    import geopandas as gpd
    from matplotlib.collections import LineCollection

    hs, extent = _hillshade_from(hillshade, dem)
    crossings = gpd.read_file(str(_read(stem, "_perimeter_crossings.fgb")))
    outlets = gpd.read_file(str(_read(stem, "_perimeter_outlets.fgb")))
    sources = gpd.read_file(str(_read(stem, "_perimeter_sources.fgb")))
    profile = pd.read_csv(_read(stem, "_perimeter_profile.csv"))

    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(hs, cmap="gray", extent=extent, origin="upper", zorder=0)
    # Closed loop: the trace starts and ends at the highest cell.
    px = np.r_[profile["x"].to_numpy(), profile["x"].iloc[0]]
    py = np.r_[profile["y"].to_numpy(), profile["y"].iloc[0]]
    ax.plot(px, py, color=PERIMETER_COLOUR, lw=0.6, zorder=1)

    handles = [plt.Line2D([], [], color=PERIMETER_COLOUR, lw=0.8,
                          label="traced perimeter")]
    order = list(VERDICT_COLOURS)
    for z, verdict in enumerate(reversed(order)):
        sub = crossings[crossings["verdict"] == verdict]
        if sub.empty:
            continue
        good = verdict in ("source", "outlet_crossing")
        thin = verdict in ("too_narrow", "not_low_spot")
        lw = 4.0 if good else (2.0 if thin else 3.0)
        ls = "-" if good else (0, (3, 1.5))
        lines = [np.asarray(part.coords)[:, :2] for g in sub.geometry
                 for part in getattr(g, "geoms", [g])]
        ax.add_collection(LineCollection(lines, colors=VERDICT_COLOURS[verdict],
                                         linewidths=lw, linestyles=ls,
                                         capstyle="round", zorder=2 + z))
        # A dot at each run's middle, so runs shorter than a screen pixel
        # still show at full extent.
        mid = np.array([ln[len(ln) // 2] for ln in lines])
        ax.scatter(mid[:, 0], mid[:, 1], s=10 if thin else 28,
                   color=VERDICT_COLOURS[verdict], edgecolor="black", lw=0.4,
                   zorder=2 + z, clip_on=False)
    for verdict in order:
        n = int((crossings["verdict"] == verdict).sum())
        if n:
            good = verdict in ("source", "outlet_crossing")
            handles.append(plt.Line2D([], [], color=VERDICT_COLOURS[verdict],
                                      lw=3.0 if good else 2.0, ls="-" if good else "--",
                                      marker="o", ms=4, mec="black", mew=0.4,
                                      label=f"flat run: {VERDICT_LABELS[verdict]} ({n})"))

    fail = crossings[(crossings["verdict"] == "failed_inward")
                     & crossings["failure_x"].notna()]
    if not fail.empty:
        ax.scatter(fail["failure_x"], fail["failure_y"], marker="x", s=40,
                   color=VERDICT_COLOURS["failed_inward"], lw=1.5, zorder=12)
        handles.append(plt.Line2D([], [], color=VERDICT_COLOURS["failed_inward"],
                                  marker="x", ls="none", ms=7,
                                  label="where the inward path got too steep"))

    sel, rest = _ranked_outlets(outlets, n_unselected_outlets)
    for sub, filled in ((rest, False), (sel, True)):
        if sub.empty:
            continue
        ax.scatter(sub.geometry.x, sub.geometry.y, marker="v", s=150,
                   facecolor=OUTLET_COLOUR if filled else "white",
                   edgecolor=OUTLET_COLOUR if not filled else "black",
                   lw=1.8 if not filled else 1.0, zorder=14, clip_on=False)
        for r in sub.itertuples():
            ax.annotate(str(int(r.rank)), (r.geometry.x, r.geometry.y),
                        xytext=(7, 5), textcoords="offset points", fontsize=8,
                        weight="bold", color="black",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                  ec="none", alpha=0.8), zorder=16)
    if not sel.empty:
        handles.append(plt.Line2D([], [], marker="v", ls="none", ms=10,
                                  mfc=OUTLET_COLOUR, mec="black",
                                  label=f"selected outlet, numbered by rank ({len(sel)})"))
    if not rest.empty:
        handles.append(plt.Line2D([], [], marker="v", ls="none", ms=10, mfc="white",
                                  mec=OUTLET_COLOUR, mew=1.8,
                                  label=f"next largest outlets, not selected "
                                        f"(of {len(outlets)} outlet runs)"))

    if not sources.empty:
        ax.scatter(sources.geometry.x, sources.geometry.y, marker="*", s=320,
                   facecolor=SOURCE_COLOUR, edgecolor="black", lw=1.0, zorder=15,
                   clip_on=False)
        for r in sources.itertuples():
            ax.annotate(f"S{int(r.source_id)}", (r.geometry.x, r.geometry.y),
                        xytext=(8, -12), textcoords="offset points", fontsize=8,
                        weight="bold", color="black",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                  ec="none", alpha=0.8), zorder=16)
    handles.append(plt.Line2D([], [], marker="*", ls="none", ms=14,
                              mfc=SOURCE_COLOUR, mec="black",
                              label=f"source ({len(sources)})"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.07),
              ncol=2, fontsize=8, frameon=False)

    if bbox is not None:
        x0, y0, x1, y1 = bbox
    else:
        x0, x1, y0, y1 = extent
        # Outlets and sources sit on the edge; leave room for their markers.
        m = 0.02 * max(x1 - x0, y1 - y0)
        x0, x1, y0, y1 = x0 - m, x1 + m, y0 - m, y1 + m
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xlabel("easting (m)")
    ax.set_ylabel("northing (m)")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set_title(title or "Perimeter valley sources: flat runs, outlets and sources")

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png


def plot_perimeter_profile(
    stem,
    out_png,
    max_gradient: float | None = None,
    n_unselected_outlets: int = 3,
    title: str | None = None,
    figsize=(11.0, 6.5),
    dpi: int = 200,
) -> Path:
    """Elevation and gradient round the perimeter, flat runs by verdict.

    Args:
        stem: lsdtt-channel-extraction output stem.
        out_png: Output PNG.
        max_gradient: Threshold line on the gradient panel; default the run's
            perimeter_max_gradient from <stem>_ingested_params.param, else 0.05.
        n_unselected_outlets: How many of the largest unselected outlet runs
            to mark besides the selected ones.
        title: Figure title.
        figsize: Figure size (inches).
        dpi: Output resolution.

    Returns:
        The path written.
    """
    import geopandas as gpd

    profile = pd.read_csv(_read(stem, "_perimeter_profile.csv"))
    crossings = gpd.read_file(str(_read(stem, "_perimeter_crossings.fgb")))
    outlets = gpd.read_file(str(_read(stem, "_perimeter_outlets.fgb")))
    sources = gpd.read_file(str(_read(stem, "_perimeter_sources.fgb")))
    if max_gradient is None:
        max_gradient = _ingested_value(stem, "perimeter_max_gradient", 0.05)

    s_km = profile["s_m"].to_numpy() / 1e3
    fig, (ax, axg) = plt.subplots(2, 1, figsize=figsize, sharex=True,
                                  gridspec_kw={"height_ratios": [2.2, 1]})
    ax.plot(s_km, profile["elevation"], color=PERIMETER_COLOUR, lw=0.7, zorder=1)
    axg.plot(s_km, profile["gradient"].clip(lower=1e-4), color=PERIMETER_COLOUR,
             lw=0.5, zorder=1)

    verdict_of = dict(zip(crossings["run_id"].astype(int), crossings["verdict"]))
    handles = [plt.Line2D([], [], color=PERIMETER_COLOUR, lw=0.8,
                          label="perimeter (from its highest cell)")]
    gentle = profile[profile["run_id"] >= 0]
    seen = []
    for run_id, g in gentle.groupby("run_id"):
        verdict = verdict_of.get(int(run_id), "too_narrow")
        # A run that wraps past the end of the trace is split at the seam.
        breaks = np.flatnonzero(np.diff(g["visit"].to_numpy()) > 1) + 1
        for part in np.split(np.arange(len(g)), breaks):
            gp = g.iloc[part]
            good = verdict in ("source", "outlet_crossing")
            ax.plot(gp["s_m"] / 1e3, gp["elevation"], color=VERDICT_COLOURS[verdict],
                    lw=3.5 if good else 2.0, ls="-" if good else (0, (3, 1.5)),
                    solid_capstyle="round", zorder=3 if good else 2)
        if verdict not in seen:
            seen.append(verdict)
    for verdict in VERDICT_COLOURS:
        if verdict in seen:
            n = int((crossings["verdict"] == verdict).sum())
            good = verdict in ("source", "outlet_crossing")
            handles.append(plt.Line2D([], [], color=VERDICT_COLOURS[verdict],
                                      lw=3.0 if good else 2.0, ls="-" if good else "--",
                                      label=f"flat run: {VERDICT_LABELS[verdict]} ({n})"))

    # Outlet points sit on the perimeter cell the outlet run is attached to:
    # the visit nearest each outlet point.
    vx, vy = profile["x"].to_numpy(), profile["y"].to_numpy()

    def nearest_visit(x, y):
        return int(np.argmin((vx - x) ** 2 + (vy - y) ** 2))

    sel, rest = _ranked_outlets(outlets, n_unselected_outlets)
    for sub, filled in ((rest, False), (sel, True)):
        for r in sub.itertuples():
            i = nearest_visit(r.geometry.x, r.geometry.y)
            ax.scatter(s_km[i], profile["elevation"].iloc[i], marker="v", s=110,
                       facecolor=OUTLET_COLOUR if filled else "white",
                       edgecolor="black" if filled else OUTLET_COLOUR,
                       lw=1.0 if filled else 1.6, zorder=6)
            ax.annotate(str(int(r.rank)), (s_km[i], profile["elevation"].iloc[i]),
                        xytext=(0, -16), textcoords="offset points", ha="center",
                        fontsize=8, weight="bold")
    if not sel.empty:
        handles.append(plt.Line2D([], [], marker="v", ls="none", ms=9,
                                  mfc=OUTLET_COLOUR, mec="black",
                                  label="selected outlet (rank)"))
    if not rest.empty:
        handles.append(plt.Line2D([], [], marker="v", ls="none", ms=9, mfc="white",
                                  mec=OUTLET_COLOUR, mew=1.6,
                                  label="next largest outlets, not selected"))
    for r in sources.itertuples():
        i = nearest_visit(r.geometry.x, r.geometry.y)
        ax.scatter(s_km[i], profile["elevation"].iloc[i], marker="*", s=260,
                   facecolor=SOURCE_COLOUR, edgecolor="black", lw=1.0, zorder=7)
        ax.annotate(f"S{int(r.source_id)}", (s_km[i], profile["elevation"].iloc[i]),
                    xytext=(0, 10), textcoords="offset points", ha="center",
                    fontsize=8, weight="bold")
    if not sources.empty:
        handles.append(plt.Line2D([], [], marker="*", ls="none", ms=13,
                                  mfc=SOURCE_COLOUR, mec="black", label="source"))

    axg.axhline(max_gradient, color=VERDICT_COLOURS["failed_inward"], lw=1.0,
                ls="--", zorder=2)
    axg.text(0.995, max_gradient, f" max_gradient {max_gradient:g}",
             transform=axg.get_yaxis_transform(), ha="right", va="bottom",
             fontsize=8)
    axg.set_yscale("log")
    axg.set_ylabel("gradient (m/m)")
    axg.set_xlabel("distance round the perimeter (km)")
    ax.set_ylabel("elevation (m)")
    ax.set_xlim(0, s_km[-1])
    axg.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.32),
               ncol=2, fontsize=8, frameon=False)
    ax.set_title(title or "Perimeter profile: flat runs, outlets and sources")

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png
