"""Smoothed channel profiles and along-channel slopes from lsdtt3.

Reads ``<write_fname>_channel_slope_points.<ext>`` from lsdtt-channel-extraction
(``write_channel_slope_points``): one point per channel cell with
``drainage_area_m2``, ``flow_distance_m``, ``elevation_m`` (pit-filled DEM),
``smoothed_elevation_m`` and ``slope_m_per_m`` (the least-squares line of
elevation against along-channel distance over +/- the
``slope_window_half_length_metres`` window: its value and its gradient),
``slope_valid``. Optionally field points snapped to the channel by the same
tool (``run_channel_point_attributes``), which carry the same
``drainage_area_m2`` / ``slope_m_per_m`` and their own measured attributes.
lsdtt3 computes every value; this module only draws.

``load_channel_points``      a slope-points or field-points file (csv, fgb, gpkg)
``plot_channel_profiles``    elevation against flow distance, raw and smoothed
``plot_slope_distance``      smoothed slope against flow distance
``plot_slope_area``          smoothed slope against drainage area (log-log)
``render_channel_slope_map`` channel cells coloured by log10 slope over a hillshade
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

RAW_COLOUR = "0.78"          # pit-filled elevation, every channel cell
SMOOTH_COLOUR = "#2a78d6"    # smoothed elevation
FIELD_COLOUR = "#d6266e"     # field points when not coloured by an attribute
INK = "#3d3d3a"

REQUIRED = ("drainage_area_m2", "slope_m_per_m")


def load_channel_points(path):
    """Load a channel slope-points or field-points file as a GeoDataFrame.

    CSV files carry ``latitude`` / ``longitude`` (WGS84, as lsdtt3 writes
    them); FlatGeobuf and GeoPackage carry their own geometry and CRS.
    """
    import geopandas as gpd

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    if p.suffix.lower() == ".csv":
        df = pd.read_csv(p)
        if not {"latitude", "longitude"} <= set(df.columns):
            raise ValueError(f"{p}: a CSV needs latitude and longitude columns")
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["longitude"], df["latitude"]),
                               crs="EPSG:4326")
    else:
        gdf = gpd.read_file(p)
    missing = [c for c in REQUIRED if c not in gdf.columns]
    if missing:
        raise ValueError(f"{p}: missing {missing}; is it from lsdtt-channel-extraction "
                         "write_channel_slope_points or run_channel_point_attributes?")
    return gdf


def _valid(df: pd.DataFrame) -> pd.DataFrame:
    if "slope_valid" in df.columns:
        df = df[df["slope_valid"].astype(float) == 1]
    s = pd.to_numeric(df["slope_m_per_m"], errors="coerce")
    return df[np.isfinite(s) & (s > 0)]


def _area_filter(df: pd.DataFrame, min_area_m2: Optional[float]) -> pd.DataFrame:
    if min_area_m2 is None:
        return df
    return df[df["drainage_area_m2"] >= min_area_m2]


def _field_scatter(ax, field: Optional[pd.DataFrame], x: str, y: str, color_by: Optional[str],
                   cmap: str, label: str):
    """Field points on top: one colour, or coloured by an attribute (with a colourbar)."""
    if field is None or len(field) == 0 or x not in field.columns or y not in field.columns:
        return None
    f = field[np.isfinite(pd.to_numeric(field[x], errors="coerce")) &
              np.isfinite(pd.to_numeric(field[y], errors="coerce"))]
    if color_by:
        if color_by not in field.columns:
            raise ValueError(f"field points have no column {color_by!r}")
        sc = ax.scatter(f[x], f[y], c=pd.to_numeric(f[color_by], errors="coerce"), cmap=cmap, s=46,
                        marker="o", edgecolors="black", linewidths=0.8, zorder=6, label=label)
        ax.figure.colorbar(sc, ax=ax, label=color_by)
        return sc
    return ax.scatter(f[x], f[y], c=FIELD_COLOUR, s=46, marker="o", edgecolors="black",
                      linewidths=0.8, zorder=6, label=label)


def _legend(ax, loc):
    """Legend with fixed-size proxy handles: the channel cells are drawn as
    1.5-point dots, the field points as 46-point circles, so one markerscale
    cannot suit both."""
    from matplotlib.lines import Line2D

    handles = []
    for h, lab in zip(*ax.get_legend_handles_labels()):
        if lab == "field points":
            colour = FIELD_COLOUR if getattr(h, "get_array", lambda: None)() is None else "0.6"
            handles.append(Line2D([], [], ls="none", marker="o", ms=7, mfc=colour, mec="black",
                                  mew=0.8, label=lab))
        else:
            fc = h.get_facecolor()
            handles.append(Line2D([], [], ls="none", marker="o", ms=5,
                                  mfc=fc[0] if len(fc) else "0.5", mec="none", label=lab))
    ax.legend(handles=handles, fontsize=8, frameon=False, loc=loc)


def _style(ax):
    ax.grid(alpha=0.25, lw=0.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def plot_channel_profiles(points: pd.DataFrame, *, field: Optional[pd.DataFrame] = None,
                          field_color_by: Optional[str] = None, min_area_m2: Optional[float] = None,
                          cmap: str = "viridis", ax=None):
    """Elevation against flow distance to the outlet for every channel cell:
    the pit-filled DEM (grey) and the smoothed profile (blue, cells with a
    valid fit), with the field points (their smoothed elevation, or elevation)
    on top. ``min_area_m2`` keeps only cells draining at least that area."""
    import matplotlib.pyplot as plt

    df = _area_filter(points, min_area_m2)
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4.5))
    else:
        fig = ax.figure
    km = df["flow_distance_m"] / 1000.0
    ax.scatter(km, df["elevation_m"], s=1.5, c=RAW_COLOUR, edgecolors="none", zorder=1,
               label="pit-filled DEM")
    v = _valid(df)
    ax.scatter(v["flow_distance_m"] / 1000.0, v["smoothed_elevation_m"], s=1.5, c=SMOOTH_COLOUR,
               edgecolors="none", zorder=2, label="smoothed profile")
    if field is not None and "flow_distance_m" in field.columns:
        f = field.assign(flow_distance_km=field["flow_distance_m"] / 1000.0)
        zcol = "smoothed_elevation_m" if "smoothed_elevation_m" in f.columns else "elevation_m"
        _field_scatter(ax, f, "flow_distance_km", zcol, field_color_by, cmap, "field points")
    ax.set_xlabel("flow distance to the outlet (km)")
    ax.set_ylabel("elevation (m)")
    _legend(ax, "upper left")
    _style(ax)
    return fig, ax


def plot_slope_distance(points: pd.DataFrame, *, field: Optional[pd.DataFrame] = None,
                        field_color_by: Optional[str] = None, min_area_m2: Optional[float] = None,
                        cmap: str = "viridis", ax=None):
    """Smoothed slope against flow distance (log slope axis), channel cells in
    grey, field points on top."""
    import matplotlib.pyplot as plt

    v = _valid(_area_filter(points, min_area_m2))
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4.5))
    else:
        fig = ax.figure
    ax.scatter(v["flow_distance_m"] / 1000.0, v["slope_m_per_m"], s=1.5, c=RAW_COLOUR,
               edgecolors="none", zorder=1, label="channel cells")
    if field is not None and "flow_distance_m" in field.columns:
        f = _valid(field).assign(flow_distance_km=lambda d: d["flow_distance_m"] / 1000.0)
        _field_scatter(ax, f, "flow_distance_km", "slope_m_per_m", field_color_by, cmap, "field points")
    ax.set_yscale("log")
    ax.set_xlabel("flow distance to the outlet (km)")
    ax.set_ylabel("smoothed channel slope (m/m)")
    _legend(ax, "upper right")
    _style(ax)
    return fig, ax


def plot_slope_area(points: pd.DataFrame, *, field: Optional[pd.DataFrame] = None,
                    field_color_by: Optional[str] = None, min_area_m2: Optional[float] = None,
                    cmap: str = "viridis", ax=None):
    """Smoothed slope against drainage area, both log axes: channel cells in
    grey, field points on top (coloured by an attribute if asked)."""
    import matplotlib.pyplot as plt

    v = _valid(_area_filter(points, min_area_m2))
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5.5))
    else:
        fig = ax.figure
    ax.scatter(v["drainage_area_m2"], v["slope_m_per_m"], s=1.5, c=RAW_COLOUR, edgecolors="none",
               zorder=1, label="channel cells")
    if field is not None:
        _field_scatter(ax, _valid(field), "drainage_area_m2", "slope_m_per_m", field_color_by, cmap,
                       "field points")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"drainage area (m$^2$)")
    ax.set_ylabel("smoothed channel slope (m/m)")
    _legend(ax, "upper right")
    _style(ax)
    return fig, ax


def render_channel_slope_map(dem_path, points, png_path, *, field=None, coordinates: str = "local",
                             cmap: str = "viridis", slope_vmin: Optional[float] = None,
                             slope_vmax: Optional[float] = None, min_area_m2: Optional[float] = None,
                             high_res: bool = False):
    """Channel cells coloured by log10 smoothed slope over a hillshade, with the
    field points as open circles. ``slope_vmin`` / ``slope_vmax`` are in m/m."""
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    from ..colour_limits import colour_limits
    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from ..render.base import grid_to_display, save_figure
    from ..render.channels import _prepare_overlay
    from ..render.decorations import add_colourbar, add_scalebar
    from ..render.hillshade import render_hillshade
    from ..render.points import render_points
    from ..render.style import MapStyle

    width = MapStyle.MAX_WIDTH
    raw = load_raster(dem_path)
    dem_crs = get_crs(raw)
    # Points from a CSV are in WGS84: put both layers in the DEM's CRS first.
    points = points.to_crs(dem_crs) if points.crs is not None else points
    if field is not None and len(field) and field.crs is not None:
        field = field.to_crs(dem_crs)
    da = to_coordinates(raw, coordinates)
    crs = get_crs(da)
    _, extent = grid_to_display(to_gmt_grid(da), crs)
    fig = render_hillshade(str(dem_path), coordinates=coordinates, shade=True,
                           style=MapStyle(cmap="gray", figure_size=width, transparent=False))
    v = _valid(_area_filter(points, min_area_m2))
    disp, _ = _prepare_overlay(v, coordinates, MapStyle(figure_size=width))
    z = np.log10(disp["slope_m_per_m"].to_numpy(dtype=float))
    lo, hi, extend = colour_limits(z, slope_vmin, slope_vmax, log_scale=True)
    pygmt.makecpt(cmap=cmap, series=[lo, hi], background=True)
    fig.plot(x=disp.geometry.x.to_numpy(), y=disp.geometry.y.to_numpy(), style="c0.06c",
             fill=z, cmap=True)
    add_colourbar(fig, label="log@-10@- slope (m/m)", extend=extend)
    if field is not None and len(field):
        fig = render_points(field, coordinates=coordinates, symbol="c", pen="1p,black",
                            style=MapStyle(point_size="0.22c", line_color="white", figure_size=width),
                            fig=fig)
    add_scalebar(fig, extent, crs)
    save_figure(fig, png_path, style=MapStyle(figure_size=width, transparent=False), high_res=high_res)
    return Path(png_path)
