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


#: GMT symbols that draw an outline only (no fill area). These are styled with a
#: pen rather than a fill.
_OUTLINE_ONLY_SYMBOLS = frozenset({"x", "y", "+", "-"})


def render_points(
    df,
    *,
    color_by: Optional[str] = None,
    size_by: Optional[str] = None,
    log_color: bool = True,
    size_range: Tuple[float, float] = (1.0, 3.0),
    symbol: str = "c",
    pen: Optional[str] = None,
    cmap: Optional[str] = None,
    colorbar: bool = False,
    colorbar_label: Optional[str] = None,
    coordinates: str = "local",
    style: Optional[MapStyle] = None,
    fig=None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    classes=None,
    n_classes: int = 5,
    class_colors=None,
    class_color_mode: str = "auto",
    legend: bool = True,
    legend_title: Optional[str] = None,
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
        symbol: GMT point-symbol code — ``"c"`` circle (default), ``"s"`` square,
            ``"t"`` triangle, ``"x"`` cross, ``"+"`` plus, etc. Outline-only
            symbols (``x``, ``y``, ``+``, ``-``) are drawn with a pen.
        pen: Outline pen (e.g. ``"1p,black"``). Defaults to ``style.line_width``
            + ``style.line_color`` for outline-only symbols; ``None`` = no
            outline on filled symbols.
        cmap: CPT name (defaults to ``style.cmap``).
        colorbar: Add a colourbar for the colour scale.
        colorbar_label: Colourbar label.
        coordinates: Coordinate system to reproject into (``"local"`` /
            ``"geographic"`` / explicit CRS).
        style: Styling options (``point_size``, ``line_color``, ``cmap``).
        fig: Existing ``pygmt.Figure`` to compose onto, or ``None`` for a new one.
        vmin, vmax: Colour-scale limits in data units (``None``: the data's
            range). With ``classes`` they are the ends of the colour gradient
            and the range the schemes work in.
        classes: Colour the points in classes of ``color_by`` instead of a
            continuous scale: the breaks (a list, tuple or NumPy array of
            increasing values), or a scheme (``"equal_interval"``,
            ``"quantile"``, ``"natural_breaks"`` / ``"jenks"``). ``log_color``
            makes the schemes and the gradient work in log10. A legend of the
            classes replaces the colourbar. See
            :mod:`lsdviztools3.colour_classes`.
        n_classes: Number of classes for a scheme (default 5).
        class_colors: One colour per class instead of sampling ``cmap``.
        class_color_mode: ``"gradient"``, ``"distinct"`` or ``"auto"`` (see
            :func:`lsdviztools3.render.channels.render_drape`).
        legend: Draw the class legend (with ``classes``; default True).
        legend_title: Class legend title (GMT text; default ``colorbar_label``).

    Returns:
        The ``pygmt.Figure``.
    """
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    style = style or MapStyle()
    gdf_display, crs = _prepare_overlay(df, coordinates, style)
    if gdf_display.empty:
        raise ValueError("Cannot render points: input has no features.")
    if classes is not None and color_by is None:
        raise ValueError("classes needs color_by: the column whose values are classified.")
    fig = _ensure_overlay_figure(gdf_display, crs, style, fig)

    xs = gdf_display.geometry.x.to_numpy()
    ys = gdf_display.geometry.y.to_numpy()

    outline_only = symbol in _OUTLINE_ONLY_SYMBOLS

    # Per-point size from log10(size_by) (numeric cm array passed via size=),
    # else a fixed size baked into the symbol string. A scalar size= does not
    # render, so a fixed size MUST be part of the style token (e.g. "c0.2c").
    if size_by is not None:
        if size_by not in gdf_display.columns:
            raise KeyError(f"size_by={size_by!r} not in point attributes.")
        size = _log_scaled_sizes(gdf_display[size_by].to_numpy(), size_range)
        plot_kwargs = dict(x=xs, y=ys, style=symbol, size=size)
    else:
        plot_kwargs = dict(x=xs, y=ys, style=f"{symbol}{style.point_size}")

    # Outline-only symbols carry a pen rather than a fill.
    outline_pen = pen or (
        f"{style.line_width},{style.line_color}" if outline_only else None
    )

    if classes is not None:
        from ..colour_classes import make_colour_classes

        if color_by not in gdf_display.columns:
            raise KeyError(f"color_by={color_by!r} not in point attributes.")
        raw = np.asarray(gdf_display[color_by].to_numpy(), dtype=float)
        if not np.isfinite(raw).any():
            raise ValueError(f"No finite values to colour by in {color_by!r}.")
        colour_classes = make_colour_classes(
            raw, classes, n_classes=n_classes, cmap=cmap or style.cmap,
            class_colors=class_colors, vmin=vmin, vmax=vmax, log_scale=log_color,
            color_mode=class_color_mode)
        index = colour_classes.classify(raw)
        if log_color:
            index[~(raw > 0)] = np.nan
        for k, colour in enumerate(colour_classes.colors):
            sel = index == k
            if not sel.any():
                continue
            kw = dict(plot_kwargs, x=xs[sel], y=ys[sel])
            if "size" in kw:
                kw["size"] = np.asarray(kw["size"])[sel]
            if outline_only:
                fig.plot(pen=f"{style.line_width},{colour}", **kw)
            else:
                fig.plot(fill=colour, pen=pen, **kw)
        if legend:
            from .categorical import add_class_legend

            add_class_legend(fig, colour_classes,
                             symbol="square" if outline_only else symbol,
                             title=legend_title if legend_title is not None else colorbar_label)
        return fig

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
        from ..colour_limits import colour_limits

        lo, hi, extend = colour_limits(z, vmin, vmax, log_scale=log_color)
        # background=True: values beyond the limits take the end colours.
        pygmt.makecpt(cmap=cmap or style.cmap, series=[lo, hi], background=True)
        if outline_only:
            fig.plot(cmap=True, pen=outline_pen, **plot_kwargs)
        else:
            fig.plot(fill=z, cmap=True, **plot_kwargs)
        if colorbar:
            add_colourbar(fig, label=colorbar_label, extend=extend)
    elif outline_only:
        fig.plot(pen=outline_pen, **plot_kwargs)
    else:
        fig.plot(fill=style.line_color, pen=outline_pen, **plot_kwargs)

    return fig
