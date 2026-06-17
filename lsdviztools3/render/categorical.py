"""Categorical and binary drape rendering.

This module is a deferred extension point in v1.
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


def render_binary_drape(base_raster, mask_raster, *, color="red", style: MapStyle, fig=None):
    """Drape a binary mask (true/false) over a base raster.

    Args:
        base_raster: Base grid (e.g. hillshade) — ``DataArray`` or path.
        mask_raster: Binary mask grid — ``DataArray`` or path.
        color: Fill colour for the true cells.
        style: Styling options.
        fig: Existing ``pygmt.Figure`` to compose onto, or ``None``.

    Returns:
        The ``pygmt.Figure``.
    """
    raise NotImplementedError(_DEFERRED_MSG)
