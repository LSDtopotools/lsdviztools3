"""Hillshade rendering."""

from __future__ import annotations

import numpy as np
from typing import Optional

from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
from .base import (
    add_crs_annotations,
    auto_projection,
    display_frame,
    grid_to_display,
    new_figure,
    set_grid_gtype,
)
from .style import MapStyle


def render_hillshade(
    raster,
    *,
    coordinates: str = "local",
    shade: bool = True,
    style: Optional[MapStyle] = None,
    fig=None,
    decorate: bool = True,
    clim: Optional[tuple] = None,
):
    """Render an elevation grid as shaded relief, or a ready-made hillshade.

    With ``shade=True`` (default) the grid is treated as elevation: it is
    coloured with the style colour map and given automatic GMT hillshade
    illumination. With ``shade=False`` the grid is drawn flat with the colour
    map only — use this for a pre-computed hillshade raster (e.g. lsdtt3's
    ``_hillshade.tif``) with ``MapStyle(cmap="gray")``.

    Args:
        raster: Grid (``xarray.DataArray``) or a path to one.
        coordinates: ``"local"`` (native projected/UTM, the default) or
            ``"geographic"`` (EPSG:4326), or an explicit CRS. See
            ``io.raster.to_coordinates``.
        shade: Apply hillshade illumination (elevation input) or not (the input
            is already a hillshade).
        style: Styling options (colour map, figure size, projection/frame). With
            ``style.clean`` the map furniture is stripped and a solid-rectangle
            scale bar is drawn (see ``style.clean``).
        fig: Existing ``pygmt.Figure`` to draw onto, or ``None`` for a new one.
        decorate: Add the "clean" scale bar when ``style.clean`` (default True).
            Layering callers (e.g. ``render_drape``) set this ``False`` for the
            base layer and add the scale bar themselves on top.
        clim: Optional ``(zmin, zmax)`` colour limits. When given, the colour map
            is built against this fixed range instead of the per-raster min/max.
            Use this to keep one colour scale across a sequence of frames (e.g.
            an animation of a transient model run).

    Returns:
        The ``pygmt.Figure``.
    """
    import pygmt

    style = style or MapStyle()
    grid = raster if hasattr(raster, "rio") else load_raster(raster)
    grid = to_gmt_grid(to_coordinates(grid, coordinates))
    crs = get_crs(grid)

    # Display units: km for projected CRS (avoids 6-7 digit axis labels).
    grid, extent = grid_to_display(grid, crs)
    set_grid_gtype(grid, crs)

    projection = style.projection or auto_projection(extent, crs, style.figure_size)
    frame = display_frame(crs, style)

    fig = fig or new_figure(style)

    if clim is not None:
        zmin, zmax = float(clim[0]), float(clim[1])
    else:
        values = np.asarray(grid.values, dtype=float)
        if not np.isfinite(values).any():
            raise ValueError("Hillshade raster has no finite values to render.")
        zmin = float(np.nanmin(values))
        zmax = float(np.nanmax(values))
    if zmin == zmax:
        zmax = zmin + 1.0
    # Route through resolve_cmap so matplotlib names work. For grey maps, only
    # apply the anti-red-NaN grey nudge when transparency is enabled.
    from .colormaps import resolve_cmap
    pygmt.makecpt(
        cmap=resolve_cmap(
            style.cmap,
            avoid_pure_gray_bug=style.transparent,
        )
        or style.cmap,
        series=[zmin, zmax],
    )

    fig.grdimage(
        grid=grid,
        region=extent,
        projection=projection,
        frame=frame,
        cmap=True,
        shading=shade,
        nan_transparent=style.transparent,
    )
    if style.clean:
        if decorate:
            from .decorations import add_bar_scalebar

            add_bar_scalebar(fig, extent, crs)
    else:
        add_crs_annotations(fig, extent, crs, style)
    return fig
