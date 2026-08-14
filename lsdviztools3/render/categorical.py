"""Categorical and binary drape rendering.

``render_binary_drape`` is implemented. ``render_categorical_drape`` remains a
deferred extension point in v1.
"""

from __future__ import annotations

from typing import Optional

from .style import MapStyle

_DEFERRED_MSG = (
    "lsdviztools3.render.categorical is deferred in v1 and not implemented yet."
)


def render_categorical_drape(base_raster, category_raster, *, categories=None,
                             style: MapStyle, fig=None):
    """Drape a categorical (class) raster over a base raster.

    Args:
        base_raster: Base grid (e.g. hillshade) — ``DataArray`` or path.
        category_raster: Integer-coded class grid — ``DataArray`` or path.
        categories: Optional mapping of class value -> label/colour.
        style: Styling options.
        fig: Existing ``pygmt.Figure`` to compose onto, or ``None``.

    Returns:
        The ``pygmt.Figure``.
    """
    raise NotImplementedError(_DEFERRED_MSG)


def render_binary_drape(
    base_raster,
    mask_raster,
    *,
    color="red",
    transparency: float = 50.0,
    coordinates: str = "local",
    style: Optional[MapStyle] = None,
    fig=None,
    shade: bool = True,
):
    """Drape a binary mask (true/false) over a base raster.

    The base is drawn as shaded relief and the mask's true cells are painted
    over it in a single flat colour at the requested transparency, so the
    underlying topography stays readable through the drape.

    Only the true cells are drawn. False cells and nodata are made transparent
    rather than painted in a second colour -- a two-colour drape hides the base
    everywhere and defeats the point of draping at all.

    Args:
        base_raster: Base grid (elevation, or a ready-made hillshade with
            ``shade=False``) — ``DataArray`` or path.
        mask_raster: Binary mask grid — ``DataArray`` or path. Anything
            non-zero and finite counts as true.
        color: Fill colour for the true cells (any GMT colour name or R/G/B).
        transparency: Drape transparency in percent; 0 is opaque, 100 invisible.
        coordinates: ``"local"`` (default) or ``"geographic"``.
        style: Styling options for the base layer.
        fig: Existing ``pygmt.Figure`` to compose onto, or ``None``.
        shade: Treat ``base_raster`` as elevation and illuminate it (default),
            or draw it flat if it is already a hillshade.

    Returns:
        The ``pygmt.Figure``.
    """
    import numpy as np
    import pygmt

    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from .base import auto_projection, grid_to_display, set_grid_gtype
    from .hillshade import render_hillshade

    style = style or MapStyle()

    fig = render_hillshade(
        base_raster,
        coordinates=coordinates,
        shade=shade,
        style=style,
        fig=fig,
        decorate=False,
    )

    mask = mask_raster if hasattr(mask_raster, "rio") else load_raster(mask_raster)
    mask = to_gmt_grid(to_coordinates(mask, coordinates))
    crs = get_crs(mask)
    mask, extent = grid_to_display(mask, crs)
    set_grid_gtype(mask, crs)

    # Keep only the true cells; everything else becomes NaN so that
    # nan_transparent leaves the base visible.
    values = np.asarray(mask.values, dtype=float)
    drape = np.where(np.isfinite(values) & (values != 0), 1.0, np.nan)
    if not np.isfinite(drape).any():
        # Nothing to draw. Returning the base unchanged is more useful than
        # raising: an empty mask is a legitimate result (nothing changed).
        return fig

    mask = mask.copy(data=drape)

    projection = style.projection or auto_projection(extent, crs, style.figure_size)

    # A one-value CPT: every true cell takes `color`, and the series is widened
    # slightly because GMT needs a non-degenerate range.
    pygmt.makecpt(cmap=f"{color},{color}", series=[0.5, 1.5, 1.0])

    fig.grdimage(
        grid=mask,
        region=extent,
        projection=projection,
        cmap=True,
        shading=False,
        nan_transparent=True,
        transparency=transparency,
    )

    if style.clean:
        from .decorations import add_bar_scalebar

        add_bar_scalebar(fig, extent, crs)

    return fig
