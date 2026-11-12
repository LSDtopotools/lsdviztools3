"""Categorical and binary drape rendering.

``render_categorical_drape`` draws an integer-coded class raster (e.g. an lsdtt3
``_lithology_code.tif``) either as semi-transparent colours over a gray
hillshade or as flat, opaque colours on their own, with a legend keyed by class
**name** (never the bare code). ``render_binary_drape`` paints the true cells of
a mask over shaded relief.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

from .style import MapStyle

#: Colour for the "unmapped" class (code -1). Not pure grey: a grey CPT row with
#: a NaN slot trips GMT's red-NaN bug (see render.colormaps).
UNMAPPED_COLOUR = "200/200/202"


def _qualitative_palette(n: int) -> list:
    """``n`` distinct colours (matplotlib tab20, tab20b, tab20c, then cycled).

    Pure-grey entries are nudged off grey for the GMT red-NaN bug.
    """
    from matplotlib import colormaps

    rgb = []
    for name in ("tab20", "tab20b", "tab20c"):
        cmap = colormaps[name]
        rgb += [tuple(int(round(255 * c)) for c in cmap(i)[:3]) for i in range(cmap.N)]
    # tab20 alternates dark/light; take the dark half first so small legends
    # get the strongest colours.
    order = list(range(0, 20, 2)) + list(range(1, 20, 2)) + list(range(20, len(rgb)))
    rgb = [rgb[i] for i in order]
    out = []
    for i in range(n):
        r, g, b = rgb[i % len(rgb)]
        if r == g == b:
            b = b + 1 if b < 255 else b - 1
        out.append(f"{r}/{g}/{b}")
    return out


def categorical_legend_entries(
    codes: Sequence[int],
    categories: Optional[Mapping] = None,
    *,
    colors: Optional[Mapping[int, str]] = None,
    unmapped_code: int = -1,
    unmapped_label: str = "unmapped",
    require_labels: bool = True,
) -> list:
    """Pair each class code with its legend label and colour.

    Args:
        codes: Class codes to draw (typically those present in the raster).
        categories: ``code -> label``, or ``code -> (label, colour)``, or
            ``code -> {"label": ..., "color": ...}``.
        colors: Optional ``code -> colour`` overriding the default palette.
        unmapped_code: Code meaning "no unit" (drawn light grey).
        unmapped_label: Its legend label.
        require_labels: Raise if any code lacks a label (default), rather than
            falling back to the bare integer, which is meaningless to a reader.

    Returns:
        ``[(code, label, colour), ...]``: codes >= 0 ascending, then mixed (-3),
        arms disagree (-4), other negative codes, and unmapped (-1) last. Names
        starting "mixed" / "arms disagree" / "unmapped" get the reserved greys of
        :data:`render.lithology_colours.RESERVED_HEX` unless a colour is given.

    Raises:
        KeyError: a code has no label and ``require_labels`` is set.
    """
    categories = dict(categories or {})
    colors = dict(colors or {})
    from ..io.lithology import ARMS_DISAGREE_CODE, MIXED_CODE, reserved_category
    from .lithology_colours import RESERVED_HEX

    # Real classes (code >= 0) first, then the reserved lsdtt3 categories in a
    # fixed order: mixed (-3), arms disagree (-4), other negatives, unmapped.
    reserved_rank = {MIXED_CODE: 0, ARMS_DISAGREE_CODE: 1, unmapped_code: 3}

    def _key(c):
        if c == unmapped_code:
            return (2, reserved_rank[c], 0)
        if c >= 0:
            return (0, 0, c)
        return (1, reserved_rank.get(c, 2), -c)

    codes = sorted({int(c) for c in codes}, key=_key)
    named = [c for c in codes if c >= 0 or c not in reserved_rank]
    palette = dict(zip(named, _qualitative_palette(len(named))))

    entries = []
    missing = []
    for code in codes:
        spec = categories.get(code)
        label, colour = None, None
        if isinstance(spec, Mapping):
            label, colour = spec.get("label"), spec.get("color")
        elif isinstance(spec, (tuple, list)):
            label = spec[0]
            colour = spec[1] if len(spec) > 1 else None
        elif spec is not None:
            label = str(spec)
        if code == unmapped_code:
            label = label or unmapped_label
            colour = colour or UNMAPPED_COLOUR
        if label is None:
            missing.append(code)
            label = str(code)
        cat = reserved_category(label)
        if colour is None and cat is not None:
            colour = RESERVED_HEX[cat]
        colour = colors.get(code) or colour or palette.get(code, RESERVED_HEX["mixed"])
        entries.append((code, label, colour))
    if missing and require_labels:
        raise KeyError(
            f"No name for class code(s) {missing}; supply a code->name mapping "
            "(see io.lithology.resolve_code_names).")
    return entries


def neutral_shading():
    """``pygmt.config`` context keeping GMT illumination colour-neutral.

    GMT illuminates by moving each colour in HSV space toward
    ``COLOR_HSV_MIN_S``/``MAX_S`` (default 1 / 0.1). For a single near-grey
    base colour (whose hue is arbitrary, e.g. the slightly green 209/210/208 of
    the solid "lightgray" CPT) that pushes shaded slopes to a saturated hue.
    Holding saturation at 0 keeps the relief gray. Used only around the base
    layer of an illumination-only hillshade.
    """
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    return pygmt.config(COLOR_HSV_MIN_S=0, COLOR_HSV_MAX_S=0)


def _legend_text(label: str) -> str:
    """Escape GMT's text escape character."""
    return str(label).replace("@", "@@")


def legend_spec(entries, *, ncols: int = 1, header: Optional[str] = None,
                extra_rows: Sequence[str] = ()) -> str:
    """GMT legend spec: one filled square + name per ``(code, label, colour)``."""
    rows = []
    if header:
        rows.append(f"H 8p,Helvetica-Bold {_legend_text(header)}")
    if ncols > 1:
        rows.append(f"N {ncols}")
    for _code, label, colour in entries:
        rows.append(f"S 0.15c s 0.28c {colour} 0.25p,black 0.4c {_legend_text(label)}")
    if extra_rows:
        if ncols > 1:
            rows.append("N 1")
        rows += list(extra_rows)
    return "\n".join(rows) + "\n"


def add_categorical_legend(
    fig,
    entries,
    *,
    position: str = "JMR+jML+o0.3c/0c",
    ncols: Optional[int] = None,
    font_size: str = "7p",
    header: Optional[str] = None,
    extra_rows: Sequence[str] = (),
    box: str = "+gwhite+p0.5p,black",
):
    """Draw a legend of class names beside the map.

    Args:
        fig: The ``pygmt.Figure``.
        entries: ``[(code, label, colour)]`` from :func:`categorical_legend_entries`.
        position: GMT legend position (default: outside the right map edge).
        ncols: Legend columns; ``None`` = 1 up to 16 entries, else 2.
        font_size: Entry font size.
        header: Optional legend title.
        extra_rows: Extra raw GMT legend rows appended after the entries
            (e.g. a channel symbol).
        box: GMT legend box string.

    Returns:
        The same ``pygmt.Figure``.
    """
    import io as _io

    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    if not entries:
        return fig
    if ncols is None:
        ncols = 1 if len(entries) <= 16 else 2
    if ncols > 1 and "+w" not in position:
        # Multi-column legends need an explicit width: estimate from the longest
        # label (Helvetica averages ~0.55 em per character).
        pt = float(font_size.rstrip("p"))
        longest = max(len(str(lbl)) for _c, lbl, _col in entries)
        col_w = 0.6 + longest * 0.55 * pt * 0.0352778
        position = f"{position}+w{ncols * col_w:.2f}c"
    spec = legend_spec(entries, ncols=ncols, header=header, extra_rows=extra_rows)
    with pygmt.config(FONT_ANNOT_PRIMARY=font_size):
        fig.legend(spec=_io.StringIO(spec), position=position, box=box)
    return fig


def _write_class_cpt(colours: Sequence[str]) -> str:
    """Discrete CPT: index i (1-based) -> colours[i-1]; NaN transparent."""
    import tempfile

    lines = ["# COLOR_MODEL = RGB", "# Discrete categorical classes"]
    for idx, colour in enumerate(colours, start=1):
        lines.append(f"{idx - 0.5:.1f} {colour} {idx + 0.5:.1f} {colour}")
    lines += [f"B {colours[0]}", f"F {colours[-1]}", "N 255/255/255"]
    with tempfile.NamedTemporaryFile("w", suffix=".cpt", delete=False) as tf:
        tf.write("\n".join(lines) + "\n")
        return tf.name


def present_codes(category_raster, *, ignore: Sequence[int] = ()) -> list:
    """Integer class codes present (finite) in a class raster, ascending."""
    import numpy as np

    from ..io.raster import load_raster

    da = category_raster if hasattr(category_raster, "rio") else load_raster(category_raster)
    values = np.asarray(da.values, dtype=float)
    values = values[np.isfinite(values)]
    return [int(c) for c in np.unique(values.astype(np.int64)) if int(c) not in set(ignore)]


def render_categorical_drape(
    base_raster,
    category_raster,
    *,
    categories=None,
    style: Optional[MapStyle] = None,
    fig=None,
    colors: Optional[Mapping[int, str]] = None,
    coordinates: str = "local",
    hillshade: bool = True,
    shade: bool = True,
    transparency: float = 40.0,
    nodata_codes: Sequence[int] = (),
    legend: bool = True,
    legend_position: str = "JMR+jML+o0.3c/0c",
    legend_ncols: Optional[int] = None,
    legend_header: Optional[str] = None,
    legend_extra_rows: Sequence[str] = (),
    require_labels: bool = True,
    decorate: bool = True,
    base_cmap: str = "gray",
    flat_transparency: float = 0.0,
):
    """Drape a categorical (class) raster, legend keyed by class name.

    Two modes:

    * ``hillshade=True`` (default): ``base_raster`` is drawn as a **gray**
      hillshade and the class colours go on top at ``transparency`` so the
      relief reads through.
    * ``hillshade=False``: flat, opaque class colours with no hillshade;
      ``base_raster`` is ignored (may be ``None``).

    Only codes present in the raster appear in the legend. Codes in
    ``nodata_codes`` (and NaN/nodata cells) are left transparent.

    Args:
        base_raster: Base grid (elevation, or a ready-made hillshade with
            ``shade=False``) — ``DataArray`` or path; unused without hillshade.
        category_raster: Integer-coded class grid — ``DataArray`` or path.
        categories: ``code -> label`` (or ``(label, colour)`` / dict). Every
            code present must have a label unless ``require_labels=False``.
        style: Styling options (figure size, frame, clean mode).
        fig: Existing ``pygmt.Figure`` to compose onto, or ``None``.
        colors: Optional ``code -> colour`` overriding the default palette.
        coordinates: ``"local"`` (default) / ``"geographic"`` / explicit CRS.
        hillshade: Draw over a gray hillshade (semi-transparent) or flat/opaque.
        shade: With ``hillshade``: treat ``base_raster`` as elevation and
            illuminate it (default), or draw it flat if it is already a hillshade.
        transparency: Class-layer transparency (percent) over the hillshade.
            Ignored when ``hillshade=False`` (always opaque).
        nodata_codes: Codes to leave undrawn.
        legend: Draw the class legend.
        legend_position, legend_ncols, legend_header, legend_extra_rows: see
            :func:`add_categorical_legend`.
        require_labels: Refuse to draw a bare integer code in the legend.
        decorate: Add the clean solid scale bar when ``style.clean``. Callers
            that add their own furniture set this ``False``.
        base_cmap: Colour map of the hillshade base. ``"gray"`` (default, as
            ``render_drape``) maps *elevation* onto a black-to-white ramp and
            then illuminates it, so high ground is lighter than low ground
            whatever the relief. A single colour (e.g. ``"lightgray"``) gives
            an illumination-only base: brightness then varies with slope
            facing only, not with elevation.
        flat_transparency: Class-layer transparency (percent) when
            ``hillshade=False`` (default 0 = opaque); a non-zero value mutes
            the classes against the white page, e.g. under a coloured overlay.

    Returns:
        The ``pygmt.Figure``.
    """
    import numpy as np
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()
    from dataclasses import replace

    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from .base import (
        add_crs_annotations, auto_projection, display_frame, grid_to_display,
        new_figure, set_grid_gtype,
    )
    from .hillshade import render_hillshade

    style = style or MapStyle()

    grid = category_raster if hasattr(category_raster, "rio") else load_raster(category_raster)
    # Nearest-neighbour if a reprojection is needed: classes must not blend.
    if coordinates not in ("local", "utm", "projected", "native"):
        from rasterio.enums import Resampling

        from ..io.raster import resolve_target_crs

        grid = grid.rio.reproject(resolve_target_crs(grid, coordinates),
                                  resampling=Resampling.nearest)
    grid = to_gmt_grid(to_coordinates(grid, "local"))
    crs = get_crs(grid)
    grid, extent = grid_to_display(grid, crs)

    values = np.asarray(grid.values, dtype=float)
    finite = np.isfinite(values)
    codes = sorted({int(c) for c in np.unique(values[finite].astype(np.int64))}
                   - {int(c) for c in nodata_codes})
    if not codes:
        raise ValueError("Category raster has no classes to draw.")
    entries = categorical_legend_entries(
        codes, categories, colors=colors, require_labels=require_labels)

    # Codes -> 1..K so a discrete CPT colours each class exactly.
    index = np.full(values.shape, np.nan)
    for k, (code, _label, _colour) in enumerate(entries, start=1):
        index[finite & (values == code)] = float(k)
    grid = grid.copy(data=index)
    set_grid_gtype(grid, crs)

    projection = style.projection or auto_projection(extent, crs, style.figure_size)

    if hillshade:
        if base_raster is None:
            raise ValueError("hillshade=True needs a base_raster (DEM or hillshade).")
        import contextlib

        ctx = neutral_shading() if base_cmap != "gray" else contextlib.nullcontext()
        with ctx:
            fig = render_hillshade(
                base_raster, coordinates=coordinates, shade=shade,
                style=replace(style, cmap=base_cmap), fig=fig, decorate=False)
        layer_transparency = transparency
        frame = None
    else:
        fig = fig or new_figure(style)
        layer_transparency = flat_transparency
        frame = display_frame(crs, style)

    pygmt.makecpt(cmap=_write_class_cpt([c for _k, _l, c in entries]),
                  series=[0.5, len(entries) + 0.5], continuous=False)
    kwargs = {}
    if frame is not None and layer_transparency == 0:
        kwargs["frame"] = frame
    if layer_transparency > 0:
        kwargs["transparency"] = layer_transparency
    fig.grdimage(grid=grid, region=extent, projection=projection, cmap=True,
                 nan_transparent=True, **kwargs)
    if frame is not None and layer_transparency > 0:
        # Frame separately so it is not faded with a muted class layer.
        fig.basemap(region=extent, projection=projection, frame=frame)
    if not hillshade and not style.clean:
        add_crs_annotations(fig, extent, crs, style)

    if legend:
        add_categorical_legend(
            fig, entries, position=legend_position, ncols=legend_ncols,
            header=legend_header, extra_rows=legend_extra_rows)
    if style.clean and decorate:
        from .decorations import add_bar_scalebar

        add_bar_scalebar(fig, extent, crs)
    return fig


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
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

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


#: Approximate Helvetica character widths, as a fraction of the font size:
#: regular text (digits are 0.556 em) and bold headers.
_CHAR_EM = 0.58
_BOLD_CHAR_EM = 0.64
_CM_PER_PT = 2.54 / 72.0


def class_legend_width(labels, *, font_size: str = "7p", header: Optional[str] = None,
                       header_size: str = "8p", symbol_width_cm: float = 0.5) -> float:
    """Width (cm) of a one-column class legend that holds its longest text.

    GMT does not measure text when it sizes a legend box, so the width is
    estimated from the longest label (and the header) with generous Helvetica
    character widths, plus the symbol, the gaps and the box padding.
    """
    pt = float(str(font_size).rstrip("p"))
    longest = max((len(str(lbl)) for lbl in labels), default=0)
    width = 0.15 + symbol_width_cm + 0.2 + longest * _CHAR_EM * pt * _CM_PER_PT + 0.3
    if header:
        # GMT escapes such as @- (subscript) take no width.
        visible = len(str(header).replace("@-", "").replace("@+", "").replace("@@", "@"))
        hpt = float(str(header_size).rstrip("p"))
        width = max(width, 0.3 + visible * _BOLD_CHAR_EM * hpt * _CM_PER_PT + 0.3)
    return round(width, 2)


def add_class_legend(
    fig,
    colour_classes,
    *,
    symbol: str = "square",
    line_width: str = "2p",
    title: Optional[str] = None,
    position: str = "JMR+jML+o0.3c/0c",
    font_size: str = "7p",
    box: str = "+gwhite+p0.5p,black",
):
    """Legend for colour classes: one entry per class, lowest at the top.

    Drawn instead of a colourbar when a map is coloured in classes (the
    ``classes`` argument of ``render_drape``, ``render_channels`` and
    ``render_points``). By default it sits outside the right edge of the map,
    so it covers no data; the saved figure includes it.

    Args:
        fig: The ``pygmt.Figure``.
        colour_classes: A :class:`lsdviztools3.colour_classes.ColourClasses`.
        symbol: ``"square"`` (rasters, polygons), ``"line"`` (channels) or a
            GMT point symbol code (``"c"`` circle, ``"t"`` triangle, ...).
        line_width: Pen width of the ``"line"`` symbol.
        title: Legend title, in GMT text (``"k@-sn@-"`` gives a subscript).
        position: GMT legend position. Its width is always set from the text
            (``+w``), so labels and title stay inside the box.
        font_size: Entry font size.
        box: GMT legend box string.

    Returns:
        The same ``pygmt.Figure``.
    """
    import io as _io

    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    labels = list(colour_classes.labels)
    colours = list(colour_classes.colors)
    symbol_width = 0.6 if symbol == "line" else 0.4
    width = class_legend_width(labels, font_size=font_size, header=title,
                               symbol_width_cm=symbol_width)
    rows = []
    if title:
        rows.append(f"H 8p,Helvetica-Bold {title}")
        rows.append("D 0.1c 0p")          # a small gap, no line
    for label, colour in zip(labels, colours):
        text = _legend_text(label)
        if symbol == "line":
            rows.append(f"S 0.35c - 0.6c - {line_width},{colour} 0.85c {text}")
        elif symbol == "square":
            rows.append(f"S 0.2c s 0.3c {colour} 0.25p,black 0.6c {text}")
        else:
            rows.append(f"S 0.2c {symbol} 0.25c {colour} 0.25p,black 0.6c {text}")
    if "+w" not in position:
        position = f"{position}+w{width}c"
    with pygmt.config(FONT_ANNOT_PRIMARY=font_size):
        fig.legend(spec=_io.StringIO("\n".join(rows) + "\n"), position=position, box=box)
    return fig
