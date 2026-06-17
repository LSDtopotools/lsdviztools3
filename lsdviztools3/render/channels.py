"""Drape and channel-network rendering."""

from __future__ import annotations

from typing import Optional

import numpy as np

from ..io.vector import load_vector, to_coordinates
from .base import (
    add_crs_label,
    auto_projection,
    display_frame,
    gdf_to_display,
    new_figure,
)
from .decorations import add_colourbar
from .style import MapStyle


def _as_gdf(data):
    """Accept a GeoDataFrame or a path; return a GeoDataFrame."""
    return data if hasattr(data, "geometry") else load_vector(data)


def _iter_line_coords(geom):
    """Yield (xs, ys) coordinate lists for a (Multi)LineString geometry."""
    if geom is None or geom.is_empty:
        return
    if geom.geom_type == "LineString":
        xs, ys = geom.xy
        yield list(xs), list(ys)
    elif geom.geom_type == "MultiLineString":
        for part in geom.geoms:
            xs, ys = part.xy
            yield list(xs), list(ys)


def _colour_z(values, log_color: bool):
    """Return z-values for colouring, optionally log10-scaled (k_sn-friendly)."""
    v = np.asarray(values, dtype=float)
    if not log_color:
        return v
    out = np.full_like(v, np.nan)
    pos = v > 0
    out[pos] = np.log10(v[pos])
    return out


def _prepare_overlay(data, coordinates, style):
    """Reproject and scale a vector to display units. Returns the display gdf."""
    gdf = to_coordinates(_as_gdf(data), coordinates)
    return gdf_to_display(gdf, gdf.crs), gdf.crs


def _ensure_overlay_figure(gdf_display, crs, style, fig):
    """If fig is None, build a publication-framed figure around the gdf bounds."""
    if fig is not None:
        return fig
    if gdf_display.empty:
        raise ValueError("Cannot render overlay: input vector has no features.")
    fig = new_figure(style)
    west, south, east, north = gdf_display.total_bounds
    if not all(v == v for v in (west, south, east, north)):
        raise ValueError("Cannot render overlay: vector bounds are not finite.")
    region = [west, east, south, north]
    projection = style.projection or auto_projection(region, crs, style.figure_size)
    fig.basemap(region=region, projection=projection, frame=display_frame(crs, style))
    add_crs_label(fig, region, crs)
    return fig


def render_drape(base_raster, drape_raster, *, style: MapStyle, fig=None):
    """Drape a coloured raster over a base (hillshade) raster.

    Args:
        base_raster: Base grid (e.g. hillshade) — ``DataArray`` or path.
        drape_raster: Grid to colour-drape on top — ``DataArray`` or path.
        style: Styling options (drape ``cmap`` etc.).
        fig: Existing ``pygmt.Figure`` to compose onto, or ``None``.

    Returns:
        The ``pygmt.Figure``.
    """
    # TODO: grdimage base (shading) then grdimage drape with transparency/cmap.
    raise NotImplementedError


def render_channels(
    lines_or_points,
    *,
    color_by=None,
    log_color: bool = True,
    cmap: Optional[str] = None,
    colorbar: bool = False,
    colorbar_label: Optional[str] = None,
    coordinates: str = "local",
    style: Optional[MapStyle] = None,
    fig=None,
):
    """Render a channel network as lines, typically over a hillshade.

    With ``color_by`` (e.g. ``"k_sn"``) the lines are coloured through a CPT;
    ``log_color`` (default ``True``) colours by log10 of the value, the natural
    scale for k_sn. Without it the lines use ``style.line_color``.

    Args:
        lines_or_points: Channel geometry — ``GeoDataFrame`` or a path to a
            vector file (e.g. ``*_chi_fit_segment_lines.fgb`` carrying ``k_sn``).
        color_by: Column to colour lines by (e.g. ``"k_sn"``). ``None`` = uniform.
        log_color: Colour by log10 of the value (default True).
        cmap: CPT name for the colour scale (defaults to ``style.cmap``).
        colorbar: Add a colourbar for the colour scale.
        colorbar_label: Colourbar label (e.g. ``"log@-10@-k@-sn@-"``).
        coordinates: Coordinate system to reproject into so the network lines up
            with the base map (``"local"`` / ``"geographic"`` / explicit CRS).
        style: Styling options (``line_width``, ``line_color``, ``cmap``).
        fig: Existing ``pygmt.Figure`` to compose onto (e.g. from
            ``render_hillshade``). If ``None``, a new framed figure is created.

    Returns:
        The ``pygmt.Figure``.
    """
    import pygmt

    style = style or MapStyle()
    gdf_display, crs = _prepare_overlay(lines_or_points, coordinates, style)
    if gdf_display.empty:
        raise ValueError("Cannot render channels: input vector has no features.")
    fig = _ensure_overlay_figure(gdf_display, crs, style, fig)

    if color_by is None:
        fig.plot(data=gdf_display, pen=f"{style.line_width},{style.line_color}")
        return fig

    if color_by not in gdf_display.columns:
        raise KeyError(
            f"color_by={color_by!r} not in channel attributes: "
            f"{[c for c in gdf_display.columns if c != gdf_display.geometry.name]}"
        )

    z = _colour_z(gdf_display[color_by].to_numpy(), log_color)
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        raise ValueError(f"No finite values to colour by in {color_by!r}.")
    vmin, vmax = float(finite.min()), float(finite.max())
    if vmin == vmax:
        vmax = vmin + 1.0
    pygmt.makecpt(cmap=cmap or style.cmap, series=[vmin, vmax])

    # Colour each line segment by its (log10) value via a per-feature z.
    for geom, zval in zip(gdf_display.geometry, z):
        if not np.isfinite(zval):
            continue
        for xs, ys in _iter_line_coords(geom):
            fig.plot(x=xs, y=ys, zvalue=zval, cmap=True,
                     pen=f"{style.line_width}+z")

    if colorbar:
        add_colourbar(fig, label=colorbar_label)
    return fig
