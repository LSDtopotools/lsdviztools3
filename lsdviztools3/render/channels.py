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


def render_drape(
    base_raster,
    drape_raster,
    *,
    coordinates: str = "local",
    base_is_hillshade: bool = True,
    log_color: bool = False,
    transparency: float = 40.0,
    cmap: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    percentile_min: Optional[float] = None,
    percentile_max: Optional[float] = None,
    out_of_bounds: str = "clamp",
    categorical_gt_breaks=None,
    categorical_colors=None,
    colorbar: bool = False,
    colorbar_label: Optional[str] = None,
    colorbar_position: str = "bottom",
    style: Optional[MapStyle] = None,
    fig=None,
):
    """Drape a coloured, semi-transparent raster over a gray hillshade base.

    Two layers: the base raster is drawn as a **gray** hillshade, then the drape
    raster is drawn on top through ``cmap`` with ``transparency`` so the relief
    texture shows through. The drape grid must share the base's CRS and grid (run
    them through the same ``coordinates`` mode); reproject/resample beforehand if
    they don't align.

    Args:
        base_raster: Base grid — ``DataArray`` or path. With
            ``base_is_hillshade=True`` (default) it is drawn flat through a gray
            CPT (it is already a hillshade); with ``False`` it is treated as
            elevation and given GMT hillshade illumination.
        drape_raster: Grid to colour-drape on top — ``DataArray`` or path.
        coordinates: Coordinate mode for both layers (``"local"`` /
            ``"geographic"`` / explicit CRS). See ``io.raster.to_coordinates``.
        base_is_hillshade: Treat ``base_raster`` as a precomputed hillshade
            (``True``) or as elevation to illuminate (``False``).
        log_color: Colour the drape by ``log10`` of its values; values ``<= 0``
            are clamped to NaN (drawn transparent), so the colourbar reads in
            log10 units.
        transparency: GMT layer transparency for the drape, 0 (opaque) to 100
            (invisible). ~30–50 keeps both colour and terrain readable.
        cmap: Colour map for the drape (defaults to ``style.cmap``). Accepts GMT
            CPTs (``geo``, ``roma``, ``seis``, ...), matplotlib colour maps
            (``cividis``, ``Spectral``, ...; converted on the fly), a single
            named colour (``"black"``) for a solid mask, a ``colour,colour``
            list, or a ``.cpt`` path. See ``render.colormaps.resolve_cmap``. The
            base is always gray regardless of this.
        vmin, vmax: Fixed lower/upper bounds for the colour scale. ``None`` uses
            the data extent (or the percentile, if given). Override either side.
        percentile_min, percentile_max: Derive the bound from a percentile of the
            finite drape values (0–100), e.g. ``percentile_min=90`` starts the
            scale at the 90th percentile. Ignored on a side where ``vmin``/
            ``vmax`` is set.
        out_of_bounds: How to treat values outside ``[vmin, vmax]``:
            ``"clamp"`` (default) paints them the nearest end colour;
            ``"clip"`` leaves them unrendered (transparent — the hillshade shows
            through). Combine ``percentile_min`` + ``"clip"`` + a solid ``cmap``
            to mask "values above the Nth percentile" in one colour.
        categorical_gt_breaks: Optional list of thresholds for discrete classes
            based on strict ``value > threshold`` tests. Provide with
            ``categorical_colors``. Values not exceeding the smallest threshold
            remain nodata/transparent.
        categorical_colors: Optional list of colours matching
            ``categorical_gt_breaks`` (same length).
        colorbar: Add a colourbar for the drape CPT.
        colorbar_label: Colourbar label (e.g. ``"log@-10@-curvature"``).
        colorbar_position: Colourbar placement: ``"bottom"`` (default),
            ``"top"``, ``"left"``, ``"right"``.
        style: Styling options (drape ``cmap``, figure size, transparency of the
            saved background, ...).
        fig: Existing ``pygmt.Figure`` to compose onto, or ``None`` for a new one.

    Returns:
        The ``pygmt.Figure``.
    """
    import pygmt
    import tempfile
    from dataclasses import replace

    from ..io.raster import (
        get_crs as _get_crs,
        load_raster as _load_raster,
        to_coordinates as _raster_to_coordinates,
        to_gmt_grid as _to_gmt_grid,
    )
    from .base import grid_to_display as _grid_to_display, set_grid_gtype as _set_gtype
    from .colormaps import resolve_cmap
    from .hillshade import render_hillshade

    if out_of_bounds not in ("clamp", "clip"):
        raise ValueError(
            f"out_of_bounds must be 'clamp' or 'clip', got {out_of_bounds!r}."
        )
    style = style or MapStyle()

    if categorical_gt_breaks is not None or categorical_colors is not None:
        if not categorical_gt_breaks or not categorical_colors:
            raise ValueError(
                "Categorical drape needs both categorical_gt_breaks and categorical_colors."
            )
        if len(categorical_gt_breaks) != len(categorical_colors):
            raise ValueError(
                "categorical_gt_breaks and categorical_colors must have the same length."
            )

    # 1) Gray hillshade base (flat if already a hillshade, else illuminated).
    fig = render_hillshade(
        base_raster,
        coordinates=coordinates,
        shade=not base_is_hillshade,
        style=replace(style, cmap="gray"),
        fig=fig,
        decorate=False,  # clean scale bar (if any) is added on top, after the drape
    )

    # 2) Coloured, semi-transparent drape on top, through the same pipeline so it
    #    lines up with the base.
    grid = drape_raster if hasattr(drape_raster, "rio") else _load_raster(drape_raster)
    grid = _to_gmt_grid(_raster_to_coordinates(grid, coordinates))
    crs = _get_crs(grid)
    grid, extent = _grid_to_display(grid, crs)

    if log_color:
        values = np.asarray(grid.values, dtype=float)
        grid = grid.copy(data=np.log10(np.where(values > 0, values, np.nan)))
    _set_gtype(grid, crs)

    values = np.asarray(grid.values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        clamp = " (after clamping values <= 0 for log scale)" if log_color else ""
        raise ValueError(f"Drape raster has no finite values to render{clamp}.")

    if categorical_gt_breaks is not None:
        breaks = [float(v) for v in categorical_gt_breaks]
        pairs = sorted(zip(breaks, categorical_colors), key=lambda p: p[0])
        sorted_breaks = [p[0] for p in pairs]
        sorted_colors = [p[1] for p in pairs]

        values = np.asarray(grid.values, dtype=float)
        classes = np.full(values.shape, np.nan, dtype=float)
        for idx, threshold in enumerate(sorted_breaks, start=1):
            classes[values > threshold] = float(idx)
        grid = grid.copy(data=classes)

        cpt_lines = ["# COLOR_MODEL = RGB", "# Discrete categories for > threshold classes"]
        for idx, color in enumerate(sorted_colors, start=1):
            lo = idx - 0.5
            hi = idx + 0.5
            cpt_lines.append(f"{lo:.1f} {color} {hi:.1f} {color}")
        cpt_lines.append(f"B {sorted_colors[0]}")
        cpt_lines.append(f"F {sorted_colors[-1]}")
        cpt_lines.append("N 255/255/255")

        with tempfile.NamedTemporaryFile("w", suffix=".cpt", delete=False) as tf:
            tf.write("\n".join(cpt_lines) + "\n")
            cpt_path = tf.name

        pygmt.makecpt(
            cmap=cpt_path,
            series=[0.5, len(sorted_colors) + 0.5, 1],
            continuous=False,
        )
    else:
        # Colour-scale bounds: data extent, overridden by percentile, then by vmin/vmax.
        lo, hi = float(finite.min()), float(finite.max())
        if percentile_min is not None:
            lo = float(np.percentile(finite, percentile_min))
        if percentile_max is not None:
            hi = float(np.percentile(finite, percentile_max))
        if vmin is not None:
            lo = float(vmin)
        if vmax is not None:
            hi = float(vmax)
        if lo == hi:
            hi = lo + 1.0

        # "clip": drop out-of-range values so they render transparent (hillshade
        # shows through). "clamp": keep them — GMT paints them the end colours.
        if out_of_bounds == "clip":
            grid = grid.where((grid >= lo) & (grid <= hi))

        pygmt.makecpt(cmap=resolve_cmap(cmap or style.cmap), series=[lo, hi])

    projection = style.projection or auto_projection(extent, crs, style.figure_size)
    fig.grdimage(
        grid=grid,
        region=extent,
        projection=projection,
        cmap=True,
        transparency=transparency,
        nan_transparent=True,
    )
    if style.clean:
        from .decorations import add_bar_scalebar

        add_bar_scalebar(fig, extent, crs)
    elif colorbar:
        add_colourbar(fig, label=colorbar_label, placement=colorbar_position)
    return fig


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
