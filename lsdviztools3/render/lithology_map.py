"""Geologic / lithologic unit map: unit colours, names, credit, scale, north.

:func:`render_lithology_map` composes one map from lsdtt3 lithology outputs:

* the unit code raster (``<prefix>_lithology_code.tif``) drawn through
  :func:`render.categorical.render_categorical_drape` — semi-transparent over a
  gray hillshade (``hillshade=True``) or flat and opaque (``hillshade=False``);
* a legend of unit **names** (from :func:`io.lithology.resolve_code_names`);
* an optional channel-network overlay (a vector file, or the x/y nodes of the
  chi lithology node CSV), in one colour or coloured by stream order;
* a small credit line (geologic map + DEM source, from the fetch-lithology
  ``_ATTRIBUTION.txt``), a solid scale bar with its label below, and a plain
  north arrow.

All furniture is drawn in every combination of the two options.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd


from .style import MapStyle

#: Stream-order columns recognised in channel inputs, in preference order.
ORDER_COLUMNS = ("strahler_order", "stream_order", "order", "basin_order")

#: Vector file suffixes (everything else is read as a CSV of x/y nodes).
_VECTOR_SUFFIXES = {".fgb", ".gpkg", ".shp", ".geojson", ".json"}


def _order_styles(orders: Sequence[int], base_width: str = "0.6p") -> dict:
    """``{order: (colour, width)}``: gray-to-black and thicker with order.

    Deliberately achromatic, so stream order stays readable over (and is never
    confused with) the categorical unit colours.
    """
    orders = sorted(set(int(o) for o in orders))
    n = len(orders)
    w0 = float(str(base_width).rstrip("p"))
    out = {}
    for i, o in enumerate(orders):
        f = i / max(n - 1, 1)
        g = int(round(150 * (1.0 - f)))  # 150 (mid gray) -> 0 (black)
        out[o] = (f"{g}/{g}/{g + 1 if g < 255 else g}", f"{w0 * (0.7 + 1.1 * f):.2f}p")
    return out


#: Default continuous colour map for k_sn-coloured channels: sequential and
#: perceptually ordered, unlike the qualitative unit palette.
KSN_CMAP = "plasma"


def ksn_colour_limits(values, clip_percentiles=(2.0, 98.0), log_scale: bool = False):
    """Colour limits for k_sn: a percentile clip of the plotted values.

    Display scaling only (like ``render_drape``'s ``percentile_min/max``): the
    k_sn values themselves are lsdtt3's. With ``log_scale`` values <= 0 are
    dropped first and the limits are in log10 units.

    Returns:
        ``(z, lo, hi)`` — the values to colour by (log10 if requested; NaN
        where not drawable) and the limits.
    """
    v = np.asarray(values, dtype=float)
    if log_scale:
        v = np.where(v > 0, v, np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.log10(v)
    finite = v[np.isfinite(v)]
    if finite.size == 0:
        raise ValueError("No finite k_sn values to colour by"
                         + (" (after dropping values <= 0 for log scale)." if log_scale else "."))
    lo, hi = (float(x) for x in np.percentile(finite, list(clip_percentiles)))
    if lo == hi:
        hi = lo + 1.0
    return v, lo, hi


def ksn_colourbar_label(m_over_n=None, a0=None, log_scale: bool = False) -> str:
    """GMT colourbar label, e.g. ``k_sn (m^0.9; m/n = 0.45, A_0 = 1 m^2)``."""
    head = "log@-10@- k@-sn@-" if log_scale else "k@-sn@-"
    if m_over_n is None:
        return f"{head} (m@+2@~q@~@+)"
    parts = [f"m/n = {m_over_n:g}"]
    if a0 is not None:
        parts.append(f"A@-0@- = {a0:g} m@+2@+")
    return f"{head} (m@+{2 * m_over_n:g}@+; {', '.join(parts)})"


def _load_channels(channels, crs, coordinates):
    """Channels as ``("lines", GeoDataFrame)`` or ``("points", DataFrame x/y)``.

    Coordinates are converted to the map's display units (km if projected).
    A CSV must carry ``x``/``y`` in the map CRS (as the chi lithology node CSV
    does).
    """
    import geopandas as gpd
    import pandas as pd

    from ..io.lithology import NODE_COLUMNS
    from ..io.vector import load_vector, to_coordinates
    from .base import KM, gdf_to_display, is_projected

    if isinstance(channels, gpd.GeoDataFrame) or (
            isinstance(channels, (str, Path)) and Path(channels).suffix.lower() in _VECTOR_SUFFIXES):
        gdf = channels if isinstance(channels, gpd.GeoDataFrame) else load_vector(channels)
        gdf = gdf_to_display(to_coordinates(gdf, coordinates), crs)
        kinds = set(gdf.geometry.geom_type)
        if kinds <= {"Point", "MultiPoint"}:
            df = pd.DataFrame(gdf.drop(columns=gdf.geometry.name))
            df["x"], df["y"] = gdf.geometry.x.to_numpy(), gdf.geometry.y.to_numpy()
            return "points", df
        return "lines", gdf

    df = channels if isinstance(channels, pd.DataFrame) else pd.read_csv(channels)
    xc, yc = NODE_COLUMNS["x"], NODE_COLUMNS["y"]
    if xc not in df.columns or yc not in df.columns:
        raise KeyError(f"Channel CSV needs {xc!r}/{yc!r} columns; has {list(df.columns)}")
    df = df.copy()
    scale = 1.0 / KM if is_projected(crs) else 1.0
    df["x"], df["y"] = df[xc].to_numpy(float) * scale, df[yc].to_numpy(float) * scale
    return "points", df


def render_channel_overlay(
    fig,
    channels,
    crs,
    *,
    coordinates: str = "local",
    color: str = "black",
    width: str = "0.6p",
    order_column: Optional[str] = None,
    point_size: str = "0.035c",
    color_by: Optional[str] = None,
    cmap: str = KSN_CMAP,
    clip_percentiles=(2.0, 98.0),
    log_scale: bool = False,
    m_over_n: Optional[float] = None,
    a0: Optional[float] = None,
    colorbar_position: str = "JBC+o0c/0.5c+w7c/0.3c+h",
) -> list:
    """Overlay a channel network: one colour, by stream order, or by k_sn.

    Args:
        fig: The ``pygmt.Figure`` (already framed by the base map).
        channels: Path or table: a line/point vector (``.fgb``/``.gpkg``/...)
            or a CSV with ``x``/``y`` in the map CRS (e.g. the chi lithology
            node CSV), or a ``(Geo)DataFrame``.
        crs: The map CRS (for display-unit scaling).
        coordinates: Coordinate mode used for the map.
        color: Single channel colour (when not colouring by order).
        width: Line width for line inputs.
        order_column: Column holding a stream order to colour by. ``"auto"``
            picks the first of :data:`ORDER_COLUMNS` present; ``None`` = single
            colour. If the column is absent, the input is drawn in one colour.
        point_size: Symbol size for point inputs (nodes drawn as dots read as
            a line at map scale).
        color_by: A numeric column to colour by with a continuous colour map
            and its own colourbar — ``"k_sn"`` (the per-node k_sn of the chi
            lithology node CSV, or a ``k_sn`` attribute of segment lines).
            Overrides ``order_column``.
        cmap: Colour map for ``color_by`` (default :data:`KSN_CMAP`).
        clip_percentiles: Colour limits as percentiles of the drawn values
            (default 2-98); values beyond are drawn in the end colours.
        log_scale: Colour by log10 (values <= 0 are not drawn).
        m_over_n, a0: For the colourbar label; default to the data's
            ``m_over_n`` / ``A_0`` columns if present.
        colorbar_position: GMT ``-D`` position of the colourbar.

    Returns:
        Extra GMT legend rows describing the overlay.
    """
    kind, data = _load_channels(channels, crs, coordinates)
    if color_by is not None:
        return _draw_by_value(fig, kind, data, color_by, cmap=cmap, width=width,
                              point_size=point_size, clip_percentiles=clip_percentiles,
                              log_scale=log_scale, m_over_n=m_over_n, a0=a0,
                              colorbar_position=colorbar_position)
    if order_column == "auto":
        order_column = next((c for c in ORDER_COLUMNS if c in data.columns), None)
    elif order_column is not None and order_column not in data.columns:
        order_column = None

    def _draw(subset, colour, w=width):
        if kind == "lines":
            from .channels import _iter_line_coords

            for geom in subset.geometry:
                for xs, ys in _iter_line_coords(geom):
                    fig.plot(x=xs, y=ys, pen=f"{w},{colour}")
        else:
            size = point_size
            if w != width:  # scale dots with the order width
                ratio = float(w.rstrip("p")) / float(str(width).rstrip("p"))
                size = f"{float(point_size.rstrip('c')) * ratio:.4f}c"
            fig.plot(x=subset["x"].to_numpy(), y=subset["y"].to_numpy(),
                     style=f"c{size}", fill=colour)

    if order_column is None:
        _draw(data, color)
        return [f"S 0.15c - 0.28c - 1.2p,{color} 0.4c Channels"]

    orders = data[order_column].dropna().astype(int)
    styles = _order_styles(orders.unique(), width)
    rows = []
    for o in sorted(styles):
        colour, w = styles[o]
        _draw(data[data[order_column].astype("Int64") == o], colour, w)
        rows.append(f"S 0.15c - 0.28c - {w},{colour} 0.4c Order {o}")
    return rows


def render_network_status_overlay(
    fig,
    status: Mapping,
    crs,
    *,
    width: str = "0.9p",
    unanalysed_width_factor: float = 0.35,
    junction_size: str = "0.075c",
    header: Optional[str] = "Strahler order",
) -> list:
    """The whole channel network, analysed links and junctions told apart.

    ``status`` is :func:`lsdviztools3.io.channel_bearings.load_network_status`.
    Every link is drawn in its stream-order gray: analysed links at the order
    width, the rest thinner (``unanalysed_width_factor``). Confluences are
    circles: analysed ones filled black, the rest open.

    Returns:
        Extra GMT legend rows describing the overlay (under ``header``).
    """
    from ..io.channel_bearings import NETWORK_LINE_COLUMNS
    from .base import gdf_to_display

    lines = gdf_to_display(status["lines"].to_crs(crs), crs)
    order_col = NETWORK_LINE_COLUMNS["order"]
    styles = _order_styles(lines[order_col].dropna().astype(int).unique(), width)

    def _thin(w: str) -> str:
        return f"{float(w.rstrip('p')) * unanalysed_width_factor:.2f}p"

    # One GMT call per (analysed, order) group: thousands of links stay fast.
    geom_only = lines[["analysed", order_col, lines.geometry.name]]
    for analysed in (False, True):
        for o, (colour, w) in sorted(styles.items()):
            sel = geom_only[(geom_only["analysed"] == analysed) &
                            (geom_only[order_col].astype(int) == o)]
            if len(sel):
                fig.plot(data=sel[[sel.geometry.name]],
                         pen=f"{w if analysed else _thin(w)},{colour}")

    rows = ([f"H 8p,Helvetica-Bold {header}"] if header else []) + [f"S 0.15c - 0.28c - {w},{colour} 0.4c Order {o}"
            for o, (colour, w) in sorted(styles.items())]
    top_colour, top_w = styles[max(styles)]
    rows.append(f"S 0.15c - 0.28c - {_thin(top_w)},{top_colour} 0.4c "
                "Not analysed (thin)")

    junctions = status.get("junctions")
    if junctions is not None and len(junctions):
        j = gdf_to_display(junctions.to_crs(crs), crs)
        rest, done = j[~j["analysed"]], j[j["analysed"]]
        if len(rest):
            fig.plot(x=rest.geometry.x.to_numpy(), y=rest.geometry.y.to_numpy(),
                     style=f"c{junction_size}", pen="0.4p,60/60/60", fill="white")
        if len(done):
            fig.plot(x=done.geometry.x.to_numpy(), y=done.geometry.y.to_numpy(),
                     style=f"c{junction_size}", pen="0.25p,white", fill="black")
        rows.append(f"S 0.15c c {junction_size} black 0.25p,white 0.4c Analysed junction")
        rows.append(f"S 0.15c c {junction_size} white 0.4p,60/60/60 0.4c "
                    "Junction not analysed")
    return rows


#: Region fills / outlines (Okabe-Ito, colour-blind safe).
REGION_COLOURS = ("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2",
                  "#D55E00", "#F0E442")


def render_region_overlay(
    fig,
    regions,
    crs,
    extent,
    *,
    field: str,
    fill: bool = True,
    transparency: float = 72.0,
    outline_width: str = "1.4p",
    colors: Optional[Mapping[str, str]] = None,
    header: Optional[str] = None,
) -> list:
    """Region polygons (e.g. physiographic provinces) clipped to the map.

    ``regions`` is a polygon vector (path or GeoDataFrame) with a ``field`` of
    region names. With ``fill`` each region is a semi-transparent colour with a
    dark outline (for a plain hillshade); without, a coloured outline only (for
    a map that already carries unit colours).

    Returns:
        Legend rows, one per region present in the map.
    """
    import geopandas as gpd
    from shapely.geometry import box

    from ..io.vector import load_vector
    from .base import KM, gdf_to_display, is_projected

    gdf = regions if isinstance(regions, gpd.GeoDataFrame) else load_vector(regions)
    if field not in gdf.columns:
        raise KeyError(f"region field {field!r} not in {list(gdf.columns)}")
    gdf = gdf.to_crs(crs)
    scale = KM if is_projected(crs) else 1.0
    xmin, xmax, ymin, ymax = (float(v) * scale for v in extent[:4])
    gdf = gpd.clip(gdf[[field, gdf.geometry.name]], box(xmin, ymin, xmax, ymax))
    gdf = gdf[~gdf.geometry.is_empty]
    gdf = gdf_to_display(gdf.dissolve(by=field, as_index=False), crs)

    values = sorted(gdf[field].astype(str).unique())
    palette = {v: REGION_COLOURS[i % len(REGION_COLOURS)] for i, v in enumerate(values)}
    palette.update(dict(colors or {}))
    rows = [f"H 8p,Helvetica-Bold {header}"] if header else []
    for v in values:
        part = gdf[gdf[field].astype(str) == v]
        colour = palette[v]
        label = v.title() if v.isupper() else v
        if fill:
            fig.plot(data=part, fill=colour, transparency=transparency)
            fig.plot(data=part, pen="0.8p,40/40/40")
            rows.append(f"S 0.15c s 0.28c {colour}@{int(transparency)} 0.5p,40/40/40 0.4c {label}")
        else:
            fig.plot(data=part, pen=f"{outline_width},{colour}")
            rows.append(f"S 0.15c - 0.28c - {outline_width},{colour} 0.4c {label}")
    return rows


def _draw_by_value(fig, kind, data, column, *, cmap, width, point_size,
                   clip_percentiles, log_scale, m_over_n, a0, colorbar_position):
    """Draw channels coloured by a numeric column through a continuous CPT."""
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    from .colormaps import resolve_cmap
    from .decorations import add_colourbar

    if column not in data.columns:
        raise KeyError(f"color_by={column!r} not in channel columns {list(data.columns)}")
    z, lo, hi = ksn_colour_limits(data[column].to_numpy(), clip_percentiles, log_scale)
    pygmt.makecpt(cmap=resolve_cmap(cmap) or cmap, series=[lo, hi])
    ok = np.isfinite(z)
    if kind == "lines":
        from .channels import _iter_line_coords

        for geom, zval in zip(data.geometry[ok], z[ok]):
            for xs, ys in _iter_line_coords(geom):
                fig.plot(x=xs, y=ys, zvalue=float(zval), cmap=True, pen=f"{width}+z")
    else:
        fig.plot(x=data["x"].to_numpy()[ok], y=data["y"].to_numpy()[ok], fill=z[ok],
                 cmap=True, style=f"c{point_size}")
    if m_over_n is None and "m_over_n" in data.columns:
        v = pd.to_numeric(data["m_over_n"], errors="coerce").dropna()
        m_over_n = float(v.iloc[0]) if len(v) else None
    if a0 is None and "A_0" in data.columns:
        v = pd.to_numeric(data["A_0"], errors="coerce").dropna()
        a0 = float(v.iloc[0]) if len(v) else None
    p0, p1 = clip_percentiles
    label = ksn_colourbar_label(m_over_n, a0, log_scale) if column == "k_sn" else column
    label += f", colour limits {p0:g}-{p1:g}th pct"
    add_colourbar(fig, label=label, position=colorbar_position, frame=[f"xaf+l{label}"])
    return []


def render_lithology_map(
    code_raster,
    *,
    dem=None,
    names: Mapping[int, str],
    colors: Optional[Mapping[int, str]] = None,
    hillshade: bool = True,
    shade: bool = True,
    transparency: float = 35.0,
    relief: str = "illumination",
    channels=None,
    channel_color: str = "black",
    channel_width: str = "0.6p",
    channel_order_column: Optional[str] = None,
    channel_color_by: Optional[str] = None,
    ksn_cmap: str = KSN_CMAP,
    ksn_clip_percentiles=(2.0, 98.0),
    ksn_log: bool = False,
    ksn_m_over_n: Optional[float] = None,
    ksn_a0: Optional[float] = None,
    geology_mute_transparency: float = 65.0,
    network_status: Optional[Mapping] = None,
    regions=None,
    region_field: Optional[str] = None,
    region_header: Optional[str] = None,
    unit_colours: Optional[Mapping[str, str]] = None,
    credit: Optional[str] = None,
    legend_title: Optional[str] = None,
    legend_ncols: Optional[int] = None,
    legend_position: Optional[str] = None,
    legend_units: Optional[Sequence[int]] = None,
    show_unmapped: bool = True,
    scalebar_anchor: str = "BL",
    north_arrow_anchor: str = "TR",
    credit_anchor: str = "BR",
    coordinates: str = "local",
    style: Optional[MapStyle] = None,
    fig=None,
):
    """Render a geologic / lithologic unit map from an lsdtt3 code raster.

    Args:
        code_raster: ``<prefix>_lithology_code.tif`` (path or ``DataArray``).
        dem: Elevation (or a ready hillshade with ``shade=False``) for the
            hillshade; needed only when ``hillshade=True``.
        names: Unit ``code -> name`` (see
            :func:`lsdviztools3.io.lithology.resolve_code_names`). Every code
            in the raster must be named; bare codes are never shown.
        colors: Optional ``code -> colour`` overrides (take precedence over
            ``unit_colours``).
        network_status: :func:`lsdviztools3.io.channel_bearings.load_network_status`
            of a CNA run: the whole network by order, analysed links thick and
            the rest thin, junctions as circles (filled = analysed). Drawn
            instead of ``channels``.
        regions, region_field: Optional region polygons (e.g. physiographic
            provinces) drawn as coloured outlines over the units, with legend rows.
        unit_colours: Optional ``{unit name: colour}``. Unit colours otherwise
            come from :func:`render.lithology_colours.lithology_palette` on the
            names in ``names`` — the same rule the k_sn charts use, so a unit
            has one colour across a run's figures.
        hillshade: Semi-transparent unit colours over a gray hillshade (True,
            default) or flat, opaque colours and no hillshade (False).
        shade: Treat ``dem`` as elevation and illuminate it (default).
        transparency: Unit-colour transparency over the hillshade, percent.
        relief: Hillshade base: ``"illumination"`` (default; a uniform light
            gray illuminated by GMT, so one unit reads as one colour family
            everywhere) or ``"elevation"`` (the package's usual gray
            elevation ramp + illumination; high ground lighter, which makes a
            single unit look lighter on plateaus than in valleys).
        channels: Optional channel overlay (see :func:`render_channel_overlay`).
        channel_color, channel_width: Single-colour overlay styling.
        channel_order_column: Colour channels by this stream-order column
            (``"auto"`` to detect one); ``None`` = single colour. Orders are
            drawn gray-to-black and thicker with order.
        channel_color_by: ``"k_sn"`` colours channels by the per-node k_sn
            lsdtt3 wrote (continuous ``ksn_cmap`` + its own colourbar), and
            mutes the geology to ``geology_mute_transparency`` so the k_sn
            colours read first. Overrides ``channel_order_column``.
        ksn_cmap, ksn_clip_percentiles, ksn_log: k_sn colour map, percentile
            colour limits (default 2-98) and optional log10 scaling.
        ksn_m_over_n, ksn_a0: m/n and A_0 for the colourbar label (default:
            the channel table's ``m_over_n``/``A_0`` columns if present).
        geology_mute_transparency: Unit-colour transparency (percent) used
            when channels are coloured by k_sn, with or without hillshade.
        credit: Credit line (see
            :func:`lsdviztools3.io.lithology.lithology_map_credit`); ``None``
            or empty draws none.
        legend_title: Optional legend header (e.g. the lithology field name).
        legend_ncols: Legend columns (``None`` = automatic).
        legend_position: GMT legend position. ``None`` = right of the map for
            up to 16 units, else centred below it (two columns), so long
            map-unit lists do not shrink the map.
        legend_units: Restrict the drawn/legended units to these codes (others
            are left transparent). ``None`` = every unit present.
        show_unmapped: Draw/legend the ``-1`` unmapped class.
        scalebar_anchor, north_arrow_anchor, credit_anchor: Corners.
        coordinates: ``"local"`` (default) / ``"geographic"`` / explicit CRS.
        style: Map styling. ``style.clean`` strips axes and the CRS title.
        fig: Existing figure to compose onto.

    Returns:
        The ``pygmt.Figure``.
    """
    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from .base import grid_to_display
    from .categorical import present_codes, render_categorical_drape
    from .decorations import add_bar_scalebar, add_map_credit, add_north_arrow

    style = style or MapStyle()
    da = code_raster if hasattr(code_raster, "rio") else load_raster(code_raster)

    nodata_codes = []
    if not show_unmapped:
        nodata_codes.append(-1)
    if legend_units is not None:
        keep = {int(c) for c in legend_units}
        nodata_codes += [c for c in present_codes(da) if c not in keep]

    # Overlay rows are only known after drawing the channels, but the legend
    # is drawn by the drape; draw it here afterwards instead.
    if relief not in ("illumination", "elevation"):
        raise ValueError(f"relief must be 'illumination' or 'elevation', got {relief!r}.")
    from .lithology_colours import code_colours

    colors = {**code_colours(names, unit_colours), **dict(colors or {})}
    muted = channels is not None and channel_color_by is not None
    unit_transparency = geology_mute_transparency if muted else transparency
    fig = render_categorical_drape(
        dem, da, categories=names, colors=colors, coordinates=coordinates,
        hillshade=hillshade, shade=shade, transparency=unit_transparency,
        flat_transparency=geology_mute_transparency if muted else 0.0,
        nodata_codes=nodata_codes, legend=False, style=style, fig=fig,
        decorate=False,
        base_cmap="lightgray" if (relief == "illumination" and shade) else "gray")

    grid = to_gmt_grid(to_coordinates(da, coordinates))
    crs = get_crs(grid)
    _, extent = grid_to_display(grid, crs)

    from .categorical import add_categorical_legend, categorical_legend_entries

    codes = [c for c in present_codes(da) if c not in set(nodata_codes)]
    entries = categorical_legend_entries(codes, names, colors=colors)
    legend_below = len(entries) > 16

    region_rows = []
    if regions is not None:
        if not region_field:
            raise ValueError("region_field is required with regions")
        region_rows = render_region_overlay(fig, regions, crs, extent, field=region_field,
                                            fill=False, header=region_header or region_field)
    extra_rows = []
    if network_status is not None:
        extra_rows = render_network_status_overlay(fig, network_status, crs,
                                                   width=channel_width)
    elif channels is not None:
        extra_rows = render_channel_overlay(
            fig, channels, crs, coordinates=coordinates, color=channel_color,
            width=channel_width, order_column=channel_order_column,
            color_by=channel_color_by, cmap=ksn_cmap,
            clip_percentiles=ksn_clip_percentiles, log_scale=ksn_log,
            m_over_n=ksn_m_over_n, a0=ksn_a0)
        if channel_color_by is not None:
            extra_rows = [f"S 0.15c - 0.28c - 1.2p,black 0.4c Channels: colour = "
                          f"{channel_color_by} (see colour bar)"]
    extra_rows = region_rows + extra_rows

    if muted:
        # Legend swatches at the same transparency as the muted drape.
        t = int(round(unit_transparency))
        entries = [(c, lbl, f"{col}@{t}") for c, lbl, col in entries]
    if legend_position is None:
        if legend_below:
            # Leave room for a k_sn colourbar between the map and the legend.
            legend_position = "JBC+jTC+o0c/" + ("2.0c" if muted else "0.3c")
        else:
            legend_position = "JMR+jML+o0.3c/0c"
    add_categorical_legend(fig, entries, position=legend_position, ncols=legend_ncols,
                           header=legend_title, extra_rows=extra_rows)

    # Furniture last, so it sits on top of the channels; the white halo keeps
    # channel lines from running through the "2 km" and "N" glyphs.
    add_bar_scalebar(fig, extent, crs, anchor=scalebar_anchor, label_position="below",
                     label_halo="white")
    add_north_arrow(fig, extent, crs, projection=style.projection,
                    anchor=north_arrow_anchor, label_halo="white")
    add_map_credit(fig, extent, credit or "", anchor=credit_anchor)
    return fig


def render_channel_order_map(
    dem,
    channels=None,
    *,
    network_status: Optional[Mapping] = None,
    regions=None,
    region_field: Optional[str] = None,
    region_transparency: float = 72.0,
    region_header: Optional[str] = None,
    order_column: Optional[str] = "auto",
    channel_color: str = "black",
    channel_width: str = "0.6p",
    credit: Optional[str] = None,
    legend_title: Optional[str] = None,
    relief: str = "illumination",
    coordinates: str = "local",
    style: Optional[MapStyle] = None,
    fig=None,
):
    """A channel network coloured by stream order over a plain hillshade.

    The no-geology counterpart of :func:`render_lithology_map` (same overlay,
    legend and furniture), for areas without a lithology code raster.

    Args:
        dem: Elevation raster for the hillshade.
        channels: As :func:`render_channel_overlay` (e.g.
            ``*_channel_link_bearing_vectors.fgb`` with ``stream_order``).
        network_status: :func:`lsdviztools3.io.channel_bearings.load_network_status`
            of a CNA run, drawn instead of ``channels``: the whole network by
            order, analysed links thick and the rest thin, junctions as
            circles (filled = analysed).
        regions, region_field: Optional region polygons (e.g. physiographic
            provinces), filled semi-transparently under the channels.
        region_transparency: Region fill transparency, percent.
        order_column: Order column (``"auto"`` detects one; ``None`` = one colour).
        channel_color, channel_width: Single-colour styling / line width.
        credit: Optional credit line.
        legend_title: Optional legend header.
        relief: ``"illumination"`` (uniform light gray, illuminated) or
            ``"elevation"`` (gray elevation ramp + illumination).
        coordinates: Coordinate mode.
        style: Map styling.
        fig: Existing figure.

    Returns:
        The ``pygmt.Figure``.
    """
    from dataclasses import replace

    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from .base import grid_to_display
    from .decorations import add_bar_scalebar, add_map_credit, add_north_arrow
    from .hillshade import render_hillshade

    style = style or MapStyle()
    base_cmap = "lightgray" if relief == "illumination" else "gray"
    import contextlib

    from .categorical import neutral_shading

    with neutral_shading() if relief == "illumination" else contextlib.nullcontext():
        fig = render_hillshade(dem, coordinates=coordinates, shade=True,
                               style=replace(style, cmap=base_cmap), fig=fig,
                               decorate=False)
    da = dem if hasattr(dem, "rio") else load_raster(dem)
    grid = to_gmt_grid(to_coordinates(da, coordinates))
    crs = get_crs(grid)
    _, extent = grid_to_display(grid, crs)
    rows = []
    if regions is not None:
        if not region_field:
            raise ValueError("region_field is required with regions")
        rows += render_region_overlay(fig, regions, crs, extent, field=region_field,
                                      fill=True, transparency=region_transparency,
                                      header=region_header or region_field)
    if network_status is not None:
        # The legend title heads the channel rows.
        rows += render_network_status_overlay(fig, network_status, crs, width=channel_width,
                                              header=legend_title)
    elif channels is not None:
        rows += render_channel_overlay(fig, channels, crs, coordinates=coordinates,
                                       color=channel_color, width=channel_width,
                                       order_column=order_column)
    # Legend of overlay rows only (no unit entries).
    import io as _io

    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()

    title_rows = ([f"H 8p,Helvetica-Bold {legend_title}"]
                  if legend_title and network_status is None else [])
    spec = title_rows + rows
    with pygmt.config(FONT_ANNOT_PRIMARY="7p"):
        fig.legend(spec=_io.StringIO("\n".join(spec) + "\n"),
                   position="JMR+jML+o0.3c/0c", box="+gwhite+p0.5p,black")
    add_bar_scalebar(fig, extent, crs, label_position="below", label_halo="white")
    add_north_arrow(fig, extent, crs, projection=style.projection, label_halo="white")
    add_map_credit(fig, extent, credit or "")
    return fig
