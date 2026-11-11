"""Basin polygon rendering."""

from __future__ import annotations

from typing import Optional

import numpy as np

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
    color_by: Optional[str] = None,
    cmap: Optional[str] = None,
    colorbar: bool = False,
    colorbar_label: Optional[str] = None,
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
        color_by: Optional numeric column used to colour polygons by value.
        cmap: Optional CPT name for ``color_by`` mode (defaults to ``style.cmap``).
        colorbar: Add a colourbar for ``color_by`` mode.
        colorbar_label: Optional colourbar label.

    Returns:
        The ``pygmt.Figure``.
    """
    style = style or MapStyle()
    gdf_display, crs = _prepare_overlay(basins_gdf, coordinates, style)
    if gdf_display.empty:
        raise ValueError("Cannot render basins: input vector has no features.")
    fig = _ensure_overlay_figure(gdf_display, crs, style, fig)

    if color_by is None:
        fig.plot(
            data=gdf_display,
            fill=fill,
            transparency=transparency,
            pen=f"{style.line_width},{style.line_color}",
        )
    else:
        from lsdviztools3._gmt import require_pygmt
        pygmt = require_pygmt()

        from .colormaps import resolve_cmap
        from .decorations import add_colourbar

        if color_by not in gdf_display.columns:
            raise KeyError(
                f"color_by={color_by!r} not in basin attributes: "
                f"{[c for c in gdf_display.columns if c != gdf_display.geometry.name]}"
            )

        z = np.asarray(gdf_display[color_by], dtype=float)
        finite = z[np.isfinite(z)]
        if finite.size == 0:
            raise ValueError(f"No finite values to colour by in {color_by!r}.")
        vmin, vmax = float(finite.min()), float(finite.max())
        if vmin == vmax:
            vmax = vmin + 1.0

        pygmt.makecpt(cmap=resolve_cmap(cmap or style.cmap), series=[vmin, vmax])

        for geom, zval in zip(gdf_display.geometry, z):
            if not np.isfinite(zval) or geom is None or geom.is_empty:
                continue

            polygons = [geom] if geom.geom_type == "Polygon" else list(getattr(geom, "geoms", []))
            for poly in polygons:
                if poly.is_empty:
                    continue
                xs, ys = poly.exterior.xy
                fig.plot(
                    x=list(xs),
                    y=list(ys),
                    close=True,
                    fill="+z",
                    zvalue=float(zval),
                    cmap=True,
                    pen=f"{style.line_width},{style.line_color}",
                    transparency=transparency,
                )

        if colorbar:
            add_colourbar(fig, label=colorbar_label)

    if label_col and label_col in gdf_display.columns:
        for geom, label in zip(gdf_display.geometry, gdf_display[label_col]):
            point = geom.representative_point()
            fig.text(x=point.x, y=point.y, text=str(label))

    return fig
