"""Divide-migration map rendering.

``render_divide_migration`` composes the standard three layers for an
``lsdtt-divide-migration`` run — a grey hillshade base, the divide polylines,
and the channel heads coloured by one metric — by reusing
:func:`~lsdviztools3.render.hillshade.render_hillshade`,
:func:`~lsdviztools3.render.channels.render_channels` and
:func:`~lsdviztools3.render.points.render_points` rather than reimplementing
them.

Heads whose colour metric is NaN (``relief_m`` / ``gradient_to_divide`` are
legitimately NaN when the nearest divide point falls on nodata, or the head sits
on its divide) are removed before the CPT is built, so a NaN can never be
painted as if it were a value. The number removed is returned by
:func:`divide_migration_layers` and reported by the CLI.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Optional

import numpy as np

from ..io.divide_migration import load_channel_heads, load_divides, resolve_metric
from .channels import render_channels
from .hillshade import render_hillshade
from .points import render_points
from .style import MapStyle


def drop_nan_metric(heads, metric):
    """Drop heads whose colour metric is not finite.

    Args:
        heads: Channel-head ``GeoDataFrame``.
        metric: A :class:`~lsdviztools3.io.divide_migration.DivideMetric`.

    Returns:
        ``(kept, n_dropped)``.
    """
    if metric.head_column not in heads.columns:
        raise KeyError(
            f"Channel heads carry no {metric.head_column!r} column; available: "
            f"{[c for c in heads.columns if c != heads.geometry.name]}"
        )
    values = heads[metric.head_column].to_numpy(dtype=float)
    finite = np.isfinite(values)
    return heads[finite], int((~finite).sum())


def render_divide_migration(
    base_raster,
    divides,
    heads,
    *,
    metric: str = "chi",
    coordinates: str = "local",
    base_is_hillshade: bool = False,
    cmap: str = "roma",
    log_color: bool = False,
    colorbar: bool = True,
    divide_color: str = "black",
    divide_width: str = "1.6p",
    symbol: str = "c",
    point_size: Optional[str] = None,
    style: Optional[MapStyle] = None,
    fig=None,
):
    """Render divides over a hillshade with channel heads coloured by a metric.

    Args:
        base_raster: DEM or hillshade grid (``DataArray`` or path). Drawn grey.
        divides: ``<stem>_divides.fgb`` / ``.gpkg`` path or ``GeoDataFrame``.
        heads: ``<stem>_channel_heads.fgb`` / ``.gpkg`` path or ``GeoDataFrame``
            (the metrics CSV has no CRS, so use the vector here).
        metric: Head metric to colour by — ``chi`` / ``relief`` / ``gradient`` /
            ``head_elevation`` (channel-head column names accepted as aliases).
        coordinates: ``"local"`` (default), ``"geographic"``, or an explicit CRS.
            Applied to all three layers so they line up.
        base_is_hillshade: ``True`` when ``base_raster`` is already a hillshade
            (drawn flat); ``False`` (default) illuminates it as elevation.
        cmap: Colour map for the head metric.
        log_color: Colour by log10 of the metric. Default ``False`` — relief and
            gradient can be negative or zero, and chi is not a log quantity.
        colorbar: Draw a colourbar for the head colour scale.
        divide_color: Pen colour for the divide polylines.
        divide_width: Pen width for the divide polylines.
        symbol: GMT symbol code for the heads (default ``"c"``, circle).
        point_size: GMT symbol size for the heads (e.g. ``"0.12c"``). Defaults
            to ``style.point_size``.
        style: Styling options.
        fig: Existing ``pygmt.Figure`` to compose onto.

    Returns:
        The ``pygmt.Figure``.
    """
    fig, _ = divide_migration_layers(
        base_raster, divides, heads, metric=metric, coordinates=coordinates,
        base_is_hillshade=base_is_hillshade, cmap=cmap, log_color=log_color,
        colorbar=colorbar, divide_color=divide_color, divide_width=divide_width,
        symbol=symbol, point_size=point_size, style=style, fig=fig,
    )
    return fig


def divide_migration_layers(
    base_raster,
    divides,
    heads,
    *,
    metric: str = "chi",
    coordinates: str = "local",
    base_is_hillshade: bool = False,
    cmap: str = "roma",
    log_color: bool = False,
    colorbar: bool = True,
    divide_color: str = "black",
    divide_width: str = "1.6p",
    symbol: str = "c",
    point_size: Optional[str] = None,
    style: Optional[MapStyle] = None,
    fig=None,
):
    """As :func:`render_divide_migration`, also returning the NaN-drop count.

    Returns:
        ``(fig, n_dropped)`` where ``n_dropped`` is the number of channel heads
        removed because their colour metric was NaN.
    """
    style = style or MapStyle()
    if point_size is not None:
        style = replace(style, point_size=point_size)
    metric = resolve_metric(metric)

    divides_gdf = divides if hasattr(divides, "geometry") else load_divides(divides)
    heads_gdf = heads if hasattr(heads, "geometry") else load_channel_heads(heads)
    if not hasattr(heads_gdf, "geometry"):
        raise TypeError(
            "Channel heads must carry geometry for the map; pass the "
            "<stem>_channel_heads.fgb/.gpkg vector rather than the metrics CSV."
        )

    heads_gdf, n_dropped = drop_nan_metric(heads_gdf, metric)
    if heads_gdf.empty:
        raise ValueError(
            f"Every channel head has a NaN {metric.head_column!r}; nothing to "
            "colour."
        )

    import pygmt

    # The grey CPT is nudged off pure grey (see render.colormaps) so NaN stays
    # transparent, and GMT's illumination amplifies saturation on dark cells —
    # which turns a shaded grey base faintly blue. Pinning the HSV saturation
    # limits to 0 for the base layer keeps illumination on brightness only.
    with pygmt.config(COLOR_HSV_MIN_S=0, COLOR_HSV_MAX_S=0):
        fig = render_hillshade(
            base_raster,
            coordinates=coordinates,
            shade=not base_is_hillshade,
            style=replace(style, cmap="gray"),
            fig=fig,
        )
    fig = render_channels(
        divides_gdf,
        color_by=None,
        coordinates=coordinates,
        style=replace(style, line_color=divide_color, line_width=divide_width),
        fig=fig,
    )
    label = metric.axis_label.replace("Mean ", "")
    if log_color:
        label = f"log@-10@- {label}"
    fig = render_points(
        heads_gdf,
        color_by=metric.head_column,
        log_color=log_color,
        cmap=cmap,
        colorbar=colorbar and not style.clean,
        colorbar_label=label,
        symbol=symbol,
        coordinates=coordinates,
        style=style,
        fig=fig,
    )
    return fig, n_dropped
