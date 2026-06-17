"""Orthographic locator basemap rendering."""

from __future__ import annotations

import math
from typing import Optional

from ..io.raster import load_raster, to_coordinates, to_gmt_grid
from .base import new_figure
from .style import MapStyle

_EARTH_RADIUS_KM = 6371.0


def _diag_km(west: float, south: float, east: float, north: float) -> float:
    """Approximate bbox diagonal in kilometers."""
    mid_lat = 0.5 * (south + north)
    dlon_km = abs(east - west) * 111.32 * math.cos(math.radians(mid_lat))
    dlat_km = abs(north - south) * 110.57
    return math.hypot(dlon_km, dlat_km)


def _horizon_angle_deg(elevation_km: float) -> float:
    """Angular distance to horizon for observer elevation above Earth."""
    ratio = _EARTH_RADIUS_KM / (_EARTH_RADIUS_KM + max(elevation_km, 0.0))
    ratio = min(1.0, max(-1.0, ratio))
    return math.degrees(math.acos(ratio))


def _auto_elevation_km(diag_km: float) -> float:
    """Scaled default elevation from DEM size, clamped to sensible range."""
    return max(30.0, min(2000.0, 4.0 * diag_km))


def render_orthographic_basemap(
    raster,
    *,
    style: Optional[MapStyle] = None,
    fig=None,
    mode: str = "auto",
    global_area_threshold_deg2: float = 25.0,
    elevation_km: Optional[float] = None,
    view_limit_deg: Optional[float] = None,
    look_bearing: float = 0.0,
    look_vertical_angle: float = 0.0,
    view_twist: float = 0.0,
    global_horizon_deg: float = 90.0,
    bbox_pen: str = "3p,0/255/255",
    coast_pen: str = "0.35p,black",
    land: str = "gray90",
    water: str = "aliceblue",
    grid_spacing_deg: Optional[float] = None,
    show_north_arrow: bool = True,
    north_arrow: str = "jTR+o0.35c/0.35c+w0.9c+f1+l",
    show_center: bool = True,
):
    """Render a global orthographic basemap with the raster footprint outlined.

    Args:
        raster: Raster ``xarray.DataArray`` or path to raster file.
        style: Optional map style. ``figure_size`` controls globe diameter.
        fig: Existing ``pygmt.Figure`` to draw onto, or ``None``.
        mode: ``"global"``, ``"aerial"``, or ``"auto"`` (default).
        global_area_threshold_deg2: In ``auto`` mode, bboxes above this area
            use the global hemisphere style.
        elevation_km: Camera elevation above Earth in km (aerial mode).
        view_limit_deg: Half view angle in degrees for aerial mode; if omitted,
            this is derived from elevation and bbox size.
        look_bearing: Camera look bearing in degrees clockwise from north.
        look_vertical_angle: GMT perspective tilt angle (0 = nadir/straight
            down, 60 = 30 degrees above horizon). Default is 0.
        view_twist: Clockwise image twist in degrees for aerial mode.
        global_horizon_deg: Orthographic horizon limit (<= 90) in global mode.
        bbox_pen: GMT pen string for raster bounding-box outline.
        coast_pen: GMT pen string for coastlines.
        land: GMT land fill color.
        water: GMT water fill color.
        grid_spacing_deg: Graticule spacing in degrees. Defaults to 30 for
            global mode and 5 for aerial mode.
        show_north_arrow: Draw a north arrow/rose.
        north_arrow: GMT rose specification for north arrow placement/style.
        show_center: Plot a small center marker at the raster midpoint.

    Returns:
        The ``pygmt.Figure`` with an orthographic locator map.
    """
    style = style or MapStyle()
    da = raster if hasattr(raster, "rio") else load_raster(raster)
    geo = to_gmt_grid(to_coordinates(da, "geographic"))
    west, south, east, north = geo.rio.bounds()

    center_lon = 0.5 * (west + east)
    center_lat = 0.5 * (south + north)
    global_horizon_deg = max(1.0, min(90.0, float(global_horizon_deg)))

    bbox_area_deg2 = abs((east - west) * (north - south))
    mode = "aerial" if mode == "local" else mode
    if mode not in {"auto", "global", "aerial"}:
        raise ValueError("mode must be one of: 'auto', 'global', 'aerial'.")
    if mode == "auto":
        mode = "global" if bbox_area_deg2 >= global_area_threshold_deg2 else "aerial"

    if mode == "global":
        projection = f"G{center_lon}/{center_lat}/{global_horizon_deg}/{style.figure_size}"
        spacing = 30.0 if grid_spacing_deg is None else float(grid_spacing_deg)
        frame = f"g{spacing}"
    else:
        diag_km = _diag_km(west, south, east, north)
        auto_elev = _auto_elevation_km(diag_km)
        elev = auto_elev if elevation_km is None else float(elevation_km)
        horizon_deg = _horizon_angle_deg(elev)
        half_diag_deg = 0.5 * math.hypot(
            abs(north - south),
            abs(east - west) * math.cos(math.radians(center_lat)),
        )
        auto_radius = max(10.0, min(60.0, horizon_deg))
        radius_deg = float(view_limit_deg) if view_limit_deg is not None else auto_radius
        radius_deg = max(radius_deg, 1.2 * half_diag_deg)
        radius_deg = min(radius_deg, 85.0)
        vwidth = 2.0 * radius_deg
        vheight = 2.0 * radius_deg
        spacing = 5.0 if grid_spacing_deg is None else float(grid_spacing_deg)
        frame = f"g{spacing}"
        projection = (
            f"G{center_lon}/{center_lat}/{style.figure_size}"
            f"+a{look_bearing}+t{look_vertical_angle}"
            f"+v{vwidth}/{vheight}+w{view_twist}+z{elev}"
        )

    fig = fig or new_figure(style)
    fig.coast(
        region="g",
        projection=projection,
        frame=frame,
        land=land,
        water=water,
        shorelines=coast_pen,
    )

    if show_north_arrow:
        fig.basemap(region="g", projection=projection, rose=north_arrow)

    fig.plot(
        x=[west, east, east, west],
        y=[south, south, north, north],
        pen=bbox_pen,
        close=True,
    )
    fig.plot(
        x=[west, east, east, west],
        y=[south, south, north, north],
        style="c0.2c",
        fill="0/255/255",
        pen="0.8p,black",
    )

    if show_center:
        fig.plot(
            x=[center_lon],
            y=[center_lat],
            style="c0.08c",
            fill="red",
            pen="0.2p,black",
        )
    return fig
