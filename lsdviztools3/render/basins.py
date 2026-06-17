"""Basin polygon rendering."""

from __future__ import annotations

from typing import Optional

from .channels import _ensure_overlay_figure, _prepare_overlay
from .style import MapStyle


def render_basins(
    basins_gdf,
    *,
    label_col: Optional[str] = None,
    coordinates: str = "local",
    fill: str = "red",
    transparency: float = 60,
    style: Optional[MapStyle] = None,
    fig=None,
):
    """Render basin polygons semi-transparently, optionally labelled.

    Args:
        basins_gdf: Basin polygons as a ``GeoDataFrame`` or a path to a vector
            file (e.g. an lsdtt3 ``_basin_outlines.fgb`` written as polygons).
        label_col: Optional attribute column for basin labels, placed at each
            polygon's ``representative_point()`` (decision §8.4).
        coordinates: Coordinate system to reproject into so the polygons line up
            with the base map (``"local"`` / ``"geographic"`` / explicit CRS).
        fill: GMT fill colour for the polygons.
        transparency: Percent transparency for the fill (0 opaque, 100 invisible).
        style: Styling options (outline ``line_width``/``line_color``).
        fig: Existing ``pygmt.Figure`` to compose onto, or ``None`` for a new one.

    Returns:
        The ``pygmt.Figure``.
    """
    style = style or MapStyle()
    gdf_display, crs = _prepare_overlay(basins_gdf, coordinates, style)
    if gdf_display.empty:
        raise ValueError("Cannot render basins: input vector has no features.")
    fig = _ensure_overlay_figure(gdf_display, crs, style, fig)

    fig.plot(
        data=gdf_display,
        fill=fill,
        transparency=transparency,
        pen=f"{style.line_width},{style.line_color}",
    )

    if label_col and label_col in gdf_display.columns:
        for geom, label in zip(gdf_display.geometry, gdf_display[label_col]):
            point = geom.representative_point()
            fig.text(x=point.x, y=point.y, text=str(label))

    return fig
