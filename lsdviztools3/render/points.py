"""Point rendering (e.g. channel points, sample sites) over a base map.

``render_points`` plots point geometries with optional colour-by and size-by
columns. Size-by maps log10 of a column (e.g. ``drainage_area``) onto a point-
size range (default 1–3 pt), so points fatten downstream and a channel drawn as
points reads as a width-graded line.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .channels import _ensure_overlay_figure, _prepare_overlay
from .decorations import add_colourbar
from .style import MapStyle

#: GMT symbol size: centimetres per point.
_CM_PER_PT = 0.0352778


def _log_scaled_sizes(values, size_range: Tuple[float, float]) -> np.ndarray:
    """Map log10(values) linearly onto ``size_range`` (pt), returned in cm."""
    v = np.asarray(values, dtype=float)
    logv = np.full_like(v, np.nan)
    pos = v > 0
    logv[pos] = np.log10(v[pos])
    finite = logv[np.isfinite(logv)]
    lo_pt, hi_pt = size_range
    if finite.size == 0:
        return np.full_like(v, lo_pt * _CM_PER_PT)
    vmin, vmax = float(finite.min()), float(finite.max())
    if vmin == vmax:
        norm = np.zeros_like(logv)
    else:
        norm = (logv - vmin) / (vmax - vmin)
    sizes_pt = lo_pt + (hi_pt - lo_pt) * norm
    # Non-positive / NaN values fall back to the smallest size.
    sizes_pt = np.where(np.isfinite(sizes_pt), sizes_pt, lo_pt)
    return sizes_pt * _CM_PER_PT


def render_points(
    df,
    *,
    color_by: Optional[str] = None,
    size_by: Optional[str] = None,
    log_color: bool = True,
    size_range: Tuple[float, float] = (1.0, 3.0),
    cmap: Optional[str] = None,
    colorbar: bool = False,
    colorbar_label: Optional[str] = None,
    coordinates: str = "local",
    style: Optional[MapStyle] = None,
    fig=None,
):
    """Render point data, optionally coloured and/or size-scaled by a column.

    Args:
        df: Point geometry — a ``GeoDataFrame`` (e.g. from
            ``io.vector.load_vector`` on ``*_chi_data_maps.fgb`` or
            ``*_channel_network_points.fgb``) or a path to a point vector file.
        color_by: Column to colour points by (e.g. ``"k_sn"``, ``"stream_order"``).
        size_by: Column whose log10 sets point size (e.g. ``"drainage_area"``),
            so larger channels render as larger/thicker points.
        log_color: Colour by log10 of ``color_by`` (default True).
        size_range: ``(min_pt, max_pt)`` point sizes for ``size_by`` (default 1–3).
        cmap: CPT name (defaults to ``style.cmap``).
        colorbar: Add a colourbar for the colour scale.
        colorbar_label: Colourbar label.
        coordinates: Coordinate system to reproject into (``"local"`` /
            ``"geographic"`` / explicit CRS).
        style: Styling options (``point_size``, ``line_color``, ``cmap``).
        fig: Existing ``pygmt.Figure`` to compose onto, or ``None`` for a new one.

    Returns:
        The ``pygmt.Figure``.
    """
    import pygmt

    style = style or MapStyle()
    gdf_display, crs = _prepare_overlay(df, coordinates, style)
    if gdf_display.empty:
        raise ValueError("Cannot render points: input has no features.")
    fig = _ensure_overlay_figure(gdf_display, crs, style, fig)

    xs = gdf_display.geometry.x.to_numpy()
    ys = gdf_display.geometry.y.to_numpy()

    # Per-point size from log10(size_by), else the fixed style size.
    if size_by is not None:
        if size_by not in gdf_display.columns:
            raise KeyError(f"size_by={size_by!r} not in point attributes.")
        size = _log_scaled_sizes(gdf_display[size_by].to_numpy(), size_range)
    else:
        size = style.point_size

    plot_kwargs = dict(x=xs, y=ys, style="c", size=size)

    if color_by is not None:
        if color_by not in gdf_display.columns:
            raise KeyError(f"color_by={color_by!r} not in point attributes.")
        v = np.asarray(gdf_display[color_by].to_numpy(), dtype=float)
        z = v.copy()
        if log_color:
            z = np.full_like(v, np.nan)
            z[v > 0] = np.log10(v[v > 0])
        finite = z[np.isfinite(z)]
        if finite.size == 0:
            raise ValueError(f"No finite values to colour by in {color_by!r}.")
        vmin, vmax = float(finite.min()), float(finite.max())
        if vmin == vmax:
            vmax = vmin + 1.0
        pygmt.makecpt(cmap=cmap or style.cmap, series=[vmin, vmax])
        fig.plot(fill=z, cmap=True, **plot_kwargs)
        if colorbar:
            add_colourbar(fig, label=colorbar_label)
    else:
        fig.plot(fill=style.line_color, **plot_kwargs)

    return fig
