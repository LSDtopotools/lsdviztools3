"""3D perspective relief rendering via PyGMT ``grdview``.

A rayshader-style 3D terrain view: the DEM is drawn as an illuminated surface
seen from a viewing azimuth/elevation, coloured by height (or by an optional
drape raster). Stays in the PyGMT stack and is headless-friendly.

Like the 2D maps, projected DEMs are shown in km (axes via ``grid_to_display``)
and geographic DEMs in degrees; elevation (z) stays in metres.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
from .base import auto_projection, grid_to_display, is_projected, new_figure, set_grid_gtype
from .style import MapStyle


def _relief_frame(crs) -> list:
    """3D frame: labelled x/y/z axes and a back-wall box."""
    if is_projected(crs):
        return ["xaf+lEasting (km)", "yaf+lNorthing (km)",
                "zaf+lElevation (m)", "wSEnZ"]
    return ["xaf+lLongitude", "yaf+lLatitude", "zaf+lElevation (m)", "wSEnZ"]


def _figure_width_cm(figure_size: str) -> float:
    """Parse the leading number from a GMT width string like ``"12c"`` -> 12.0."""
    digits = "".join(ch for ch in str(figure_size) if (ch.isdigit() or ch == "."))
    return float(digits) if digits else 12.0


def _zsize_from_ve(ve, extent, zrange_m, figure_size, crs) -> str:
    """GMT z-height (``-JZ``) giving a true vertical:horizontal scale ratio ``ve``.

    Only meaningful for projected (metric) CRS, where the x axis is in km. With
    horizontal scale ``width_cm / (x_range_km * 1000)`` cm per ground-metre and
    vertical scale ``zsize_cm / z_range_m``, ``ve`` is their ratio, so
    ``zsize_cm = ve * z_range_m * width_cm / (x_range_km * 1000)``.
    """
    width_cm = _figure_width_cm(figure_size)
    x_range_km = float(extent[1]) - float(extent[0])
    if x_range_km <= 0 or zrange_m <= 0:
        raise ValueError("Cannot derive z-size: non-positive x-range or z-range.")
    zsize_cm = ve * zrange_m * width_cm / (x_range_km * 1000.0)
    return f"{zsize_cm:g}c"


def _finite_range(values) -> tuple:
    """(min, max) over finite values, guarding against an empty/flat grid."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError("Raster has no finite cells to render.")
    vmin, vmax = float(finite.min()), float(finite.max())
    if vmin == vmax:
        vmax = vmin + 1.0
    return vmin, vmax


def render_relief3d(
    raster,
    *,
    coordinates: str = "local",
    azimuth: float = 135.0,
    elevation: float = 30.0,
    vertical_exaggeration: Optional[float] = None,
    zsize: str = "4c",
    shade: bool = True,
    drape=None,
    cmap: Optional[str] = None,
    style: Optional[MapStyle] = None,
    fig=None,
    zlimits: Optional[tuple] = None,
):
    """Render a DEM as a 3D perspective relief surface (``grdview``).

    Args:
        raster: DEM — ``xarray.DataArray`` or path.
        coordinates: ``"local"`` (UTM, km axes) / ``"geographic"`` / explicit CRS.
        azimuth: Viewing azimuth in degrees (compass bearing of the camera).
        elevation: Viewing elevation angle in degrees above the horizon.
        vertical_exaggeration: True vertical:horizontal scale ratio. When set
            (projected CRS only) it sets ``zsize`` from the data ranges and
            overrides the ``zsize`` argument; ignored for geographic CRS.
        zsize: Plotted height of the full elevation range (GMT ``-JZ``); larger
            values exaggerate relief. Ignored when ``vertical_exaggeration`` set.
        shade: Illuminate the surface (hillshading from the topography).
        drape: Optional raster to colour the surface by instead of height; it is
            reprojected and resampled onto the DEM grid. ``None`` colours by
            elevation.
        cmap: CPT name (defaults to ``style.cmap``; ``"geo"`` suits elevation).
        style: Styling options (``cmap``, ``figure_size``, fonts, ...).
        fig: Existing ``pygmt.Figure`` to draw onto, or ``None`` for a new one.
        zlimits: Optional fixed ``(zmin, zmax)`` for the z axis, the vertical
            scale and the colour range. ``None`` uses this raster's own range.
            An animation must pin it, or every frame is rescaled and relief
            growth is invisible.

    Returns:
        The ``pygmt.Figure``.
    """
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    style = style or MapStyle()

    da = to_coordinates(load_raster(raster), coordinates)
    crs = get_crs(da)
    grid, extent = grid_to_display(to_gmt_grid(da), crs)
    set_grid_gtype(grid, crs)

    if zlimits is not None:
        zmin, zmax = float(zlimits[0]), float(zlimits[1])
        if zmax <= zmin:
            zmax = zmin + 1.0
    else:
        zmin, zmax = _finite_range(grid.values)
    region = [extent[0], extent[1], extent[2], extent[3], zmin, zmax]

    if vertical_exaggeration is not None and is_projected(crs):
        zsize = _zsize_from_ve(vertical_exaggeration, extent, zmax - zmin,
                               style.figure_size, crs)
    projection = style.projection or auto_projection(extent, crs, style.figure_size)
    perspective = [azimuth, elevation]
    frame = style.frame or _relief_frame(crs)

    if fig is None:
        fig = new_figure(style)

    common = dict(
        grid=grid, region=region, projection=projection,
        perspective=perspective, zsize=zsize, frame=frame,
        shading=shade,
    )

    if drape is not None:
        # Colour the surface by the drape raster, resampled onto the DEM grid.
        drape_da = to_coordinates(load_raster(drape), coordinates)
        drape_grid, _ = grid_to_display(to_gmt_grid(drape_da), crs)
        drape_grid = drape_grid.interp(x=grid.x, y=grid.y)
        set_grid_gtype(drape_grid, crs)
        dmin, dmax = _finite_range(drape_grid.values)
        pygmt.makecpt(cmap=cmap or style.cmap, series=[dmin, dmax])
        fig.grdview(drapegrid=drape_grid, surftype="i", cmap=True, **common)
    else:
        pygmt.makecpt(cmap=cmap or style.cmap, series=[zmin, zmax])
        fig.grdview(surftype="s", cmap=True, **common)

    return fig
