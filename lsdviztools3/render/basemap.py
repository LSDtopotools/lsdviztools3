"""Locator maps: where a raster sits, on a globe or in its region.

Three kinds, from :func:`render_orthographic_basemap` (``mode``) or directly:

- ``"global"``: a hemisphere centred on the raster, with coastlines and a
  graticule; the raster's footprint outlined.
- ``"regional"`` (:func:`render_regional_basemap`): a flat map several times
  larger than the raster, with coastlines, country and state borders, major
  rivers and lakes, optional shaded relief, a scale bar and the footprint.
- ``"aerial"``: a perspective view from a height above the raster.

``"auto"`` chooses ``global`` for large rasters and ``regional`` otherwise.
:func:`add_locator_inset` draws either of the first two as an inset in a
corner of a map that is already drawn.
"""

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
    extent_km: Optional[float] = None,
    relief: bool = False,
    regional_bbox_pen: str = "1.6p,red",
):
    """Render a locator map of a raster: a globe, its region, or an aerial view.

    Args:
        raster: Raster ``xarray.DataArray`` or path to raster file.
        style: Optional map style. ``figure_size`` controls the map width
            (the globe's diameter).
        fig: Existing ``pygmt.Figure`` to draw onto, or ``None``.
        mode: ``"global"``, ``"regional"`` (see
            :func:`render_regional_basemap`), ``"aerial"``, or ``"auto"``
            (default: global for a footprint of ``global_area_threshold_deg2``
            or more, regional otherwise).
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
        show_north_arrow: Draw a north arrow (aerial and regional maps; a
            globe has none).
        north_arrow: GMT rose specification for north arrow placement/style.
        show_center: Plot a small center marker at the raster midpoint.
        extent_km: Regional mode: width of the map (km). Default ten times
            the footprint's diagonal, at least 250 km.
        relief: Regional mode: shade the land with GMT's global relief data
            (downloaded from the GMT server on first use).
        regional_bbox_pen: Regional mode: pen of the footprint outline.

    Returns:
        The ``pygmt.Figure``.
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
    if mode not in {"auto", "global", "regional", "aerial"}:
        raise ValueError("mode must be one of: 'auto', 'global', 'regional', 'aerial'.")
    if mode == "auto":
        mode = "global" if bbox_area_deg2 >= global_area_threshold_deg2 else "regional"
    if mode == "regional":
        return render_regional_basemap(
            geo, style=style, fig=fig, bbox_pen=regional_bbox_pen, coast_pen=coast_pen,
            land=land, water=water, extent_km=extent_km, relief=relief,
            show_north_arrow=show_north_arrow)

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

    # No north arrow on the globe: north points a different way at every place
    # on a hemisphere (the graticule shows it), and GMT's rose on a whole-globe
    # orthographic map writes PostScript that Ghostscript cannot read (the PNG
    # export failed).
    if show_north_arrow and mode == "aerial":
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


def _nice_length_km(target_km: float) -> float:
    """1, 2 or 5 times a power of ten, the largest not above ``target_km``."""
    if target_km <= 0:
        return 1.0
    power = 10.0 ** math.floor(math.log10(target_km))
    for step in (5.0, 2.0, 1.0):
        if step * power <= target_km:
            return step * power
    return power


def regional_extent(west: float, south: float, east: float, north: float,
                    extent_km: Optional[float] = None):
    """Region ``[w, e, s, n]`` (degrees) of a square regional map around a footprint.

    Args:
        west, south, east, north: The footprint, WGS84 degrees.
        extent_km: Width (and height) of the map in km. Default: ten times the
            footprint's diagonal, and at least 250 km.

    Returns:
        ``([w, e, s, n], extent_km)``.
    """
    if extent_km is None:
        extent_km = max(10.0 * _diag_km(west, south, east, north), 250.0)
    extent_km = float(extent_km)
    if extent_km <= 0:
        raise ValueError(f"extent_km must be > 0, got {extent_km}.")
    clon = 0.5 * (west + east)
    clat = 0.5 * (south + north)
    half = 0.5 * extent_km
    dlat = half / 110.57
    dlon = half / (111.32 * max(math.cos(math.radians(clat)), 0.01))
    s = max(clat - dlat, -85.0)
    n = min(clat + dlat, 85.0)
    return [clon - dlon, clon + dlon, s, n], extent_km


def _relief_resolution(extent_km: float) -> str:
    """GMT earth_relief resolution for a map ``extent_km`` wide (a few hundred cells)."""
    for limit, res in ((200.0, "15s"), (800.0, "30s"), (2500.0, "01m"), (8000.0, "05m")):
        if extent_km <= limit:
            return res
    return "10m"


def render_regional_basemap(
    raster,
    *,
    style: Optional[MapStyle] = None,
    fig=None,
    extent_km: Optional[float] = None,
    relief: bool = False,
    borders: bool = True,
    rivers: bool = True,
    land: str = "gray90",
    water: str = "aliceblue",
    coast_pen: str = "0.35p,black",
    bbox_pen: str = "1.6p,red",
    show_north_arrow: bool = True,
    north_arrow: str = "jTR+o0.3c/0.3c+w0.8c+f1+l",
    scale_bar: bool = True,
    frame: bool = True,
):
    """A flat locator map of the region around a raster, its footprint outlined.

    Mercator, square, ``extent_km`` wide, centred on the raster: coastlines,
    country borders (solid) and state or province borders (dashed), permanent
    rivers and lakes, optionally shaded relief, the raster's bounding box, a
    scale bar and a north arrow. Coastline detail is chosen by GMT from the
    map's size.

    Args:
        raster: Raster ``xarray.DataArray`` or path (any CRS).
        style: Map style; ``figure_size`` is the map width.
        fig: Existing ``pygmt.Figure`` to draw onto (e.g. inside an inset), or
            ``None`` for a new one.
        extent_km: Width of the map in km. Default: ten times the footprint's
            diagonal, at least 250 km.
        relief: Shade the land with GMT's global relief grid
            (``earth_relief``, resolution chosen from ``extent_km``),
            downloaded from the GMT server and cached on first use.
        borders: Draw country and state/province borders.
        rivers: Draw permanent rivers.
        land, water: Fills (land is not filled when ``relief`` is on).
        coast_pen: Coastline pen.
        bbox_pen: Pen of the footprint outline.
        show_north_arrow: Draw a north arrow.
        north_arrow: GMT rose specification.
        scale_bar: Draw a scale bar (bottom left).
        frame: Annotate the frame with longitude and latitude; ``False`` draws
            a plain border (insets).

    Returns:
        The ``pygmt.Figure``.
    """
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    style = style or MapStyle()
    da = raster if hasattr(raster, "rio") else load_raster(raster)
    geo = to_gmt_grid(to_coordinates(da, "geographic"))
    west, south, east, north = geo.rio.bounds()
    region, extent_km = regional_extent(west, south, east, north, extent_km)
    projection = f"M{style.figure_size}"
    frame_spec = ["WSne", "af"] if frame else ["wsne", "f0"]

    fig = fig or new_figure(style)
    with pygmt.config(MAP_FRAME_TYPE="plain", FONT_ANNOT_PRIMARY="7p", FONT_LABEL="7p",
                      MAP_TICK_LENGTH_PRIMARY="2p", MAP_SCALE_HEIGHT="4p"):
        if relief:
            grid = pygmt.datasets.load_earth_relief(
                resolution=_relief_resolution(extent_km), region=region)
            shade = pygmt.grdgradient(grid=grid, radiance=[315, 45], normalize="t1")
            pygmt.makecpt(cmap="gray", series=[-1.2, 0.6])
            fig.grdimage(grid=shade, region=region, projection=projection, cmap=True,
                         transparency=55, frame=frame_spec)
            land_fill = None
        else:
            fig.basemap(region=region, projection=projection, frame=frame_spec)
            land_fill = land
        coast_kw = dict(region=region, projection=projection, water=water, lakes=water,
                        shorelines=coast_pen, resolution="auto")
        if land_fill:
            coast_kw["land"] = land_fill
        if borders:
            coast_kw["borders"] = ["1/0.9p,gray20", "2/0.5p,gray35,-"]
        if rivers:
            coast_kw["rivers"] = ["r/0.4p,79/128/181"]
        fig.coast(**coast_kw)
        fig.plot(x=[west, east, east, west], y=[south, south, north, north],
                 close=True, pen=bbox_pen)
        if scale_bar or show_north_arrow:
            basemap_kw = {}
            if scale_bar:
                length = _nice_length_km(extent_km / 4.0)
                basemap_kw["map_scale"] = f"jBL+o0.3c/0.3c+w{length:g}k+f+l"
            if show_north_arrow:
                basemap_kw["rose"] = north_arrow
            fig.basemap(region=region, projection=projection, **basemap_kw)
    return fig


def _width_cm(width: str) -> float:
    w = str(width).strip()
    for unit, factor in (("c", 1.0), ("i", 2.54), ("p", 2.54 / 72.0)):
        if w.endswith(unit):
            return float(w[:-1]) * factor
    return float(w)


def add_locator_inset(
    fig,
    raster,
    *,
    mode: str = "global",
    position: str = "TR",
    width: str = "4c",
    offset: str = "0.2c",
    box: str = "+gwhite+p0.5p,black",
    extent_km: Optional[float] = None,
    relief: bool = False,
):
    """Draw a small locator map in a corner of a map that is already drawn.

    The inset is a globe (``mode="global"``) or the region
    (``mode="regional"``, without frame annotations, scale bar or north
    arrow), with the raster's footprint outlined. Draw the main map first:
    the inset is placed on its frame.

    Args:
        fig: The ``pygmt.Figure`` holding the main map.
        raster: The raster whose footprint is shown (usually the main map's).
        mode: ``"global"`` or ``"regional"``.
        position: Corner of the main map: ``"TR"`` (default), ``"TL"``,
            ``"BR"`` or ``"BL"``, or a full GMT inset position (starting
            ``j``, ``g``, ``n`` or ``x``).
        width: Inset width (and height), e.g. ``"4c"``.
        offset: Distance from the corner.
        box: GMT box around the inset (fill and pen).
        extent_km: Regional mode: width of the region shown (km).
        relief: Regional mode: shaded relief (downloads GMT relief data).

    Returns:
        The same ``pygmt.Figure``.
    """
    if mode not in ("global", "regional"):
        raise ValueError(f"Inset mode must be 'global' or 'regional', got {mode!r}.")
    pos = str(position)
    if pos[:1] not in ("j", "g", "n", "x"):
        if pos.upper() not in ("TR", "TL", "BR", "BL"):
            raise ValueError(f"position must be TR, TL, BR, BL or a GMT inset position, got {position!r}.")
        pos = f"j{pos.upper()}+w{width}+o{offset}"
    # The map inside the box, a little smaller than it so the border shows.
    inner = f"{max(_width_cm(width) - 0.2, 0.5):.2f}c"
    inner_style = MapStyle(figure_size=inner, transparent=False)
    with fig.inset(position=pos, box=box, clearance="0.1c"):
        if mode == "global":
            render_orthographic_basemap(raster, style=inner_style, fig=fig, mode="global",
                                        show_center=False, bbox_pen="1.5p,red",
                                        coast_pen="0.2p,black", grid_spacing_deg=30.0)
        else:
            render_regional_basemap(raster, style=inner_style, fig=fig, extent_km=extent_km,
                                    relief=relief, show_north_arrow=False, scale_bar=False,
                                    frame=False, coast_pen="0.25p,black", bbox_pen="1.2p,red")
    return fig

