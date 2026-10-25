"""Ridgeline map (PyGMT).

The map half of R. Vervoordeldonk's ridgeline profile summary: a hillshade
around one ``lsdtt-ridgeline-analysis`` ridgeline, the ridge (green halo under
a black line) with a perpendicular tick every kilometre of displayed distance,
and either the channel network (native triple divide junctions) or the paired
flowpaths of one topological threshold, coloured per junction, with the
highlighted junctions circled.

The original rotated its hillshade by one of eight fixed angles to lay the
ridge across the page. Here the ridge is drawn in whatever CRS the hillshade
has: run the analysis on an lsdtt3 oblique frame DEM
(``rotate_to_linear_feature``) and the map is laid along the ridge, with the
reference graticule, a true-north arrow and an ``Oblique Mercator`` title
added by the shared render layer.

The ridge is drawn from the profile's own ``latitude`` / ``longitude`` per
node, so the clipped distance range and the kilometre ticks follow the same
reverse / range rules as the profile charts
(:func:`lsdviztools3.io.ridgeline.select_profile`).
"""

from __future__ import annotations

import colorsys
from typing import Optional

import numpy as np

from ..io.raster import get_crs, load_raster, to_coordinates
from ..io.ridgeline import tdj_rows_in_range
from ..io.vector import load_vector
from .base import KM, is_projected
from .style import MapStyle

#: Ridge halo, as in the profile charts (alpha 0.35 is GMT transparency 65).
RIDGE_HALO_PEN = "2.4p,#178a3b@65"
RIDGE_PEN = "0.6p,black"
CHANNEL_PEN = "0.6p,#1f77b4@20"
TDJ_RED = "#d62728"


def _lonlat_to_crs(lon, lat, crs):
    from pyproj import Transformer

    transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = transformer.transform(np.asarray(lon, dtype=float), np.asarray(lat, dtype=float))
    return np.asarray(x, dtype=float), np.asarray(y, dtype=float)


def junction_colours(node_indices) -> dict:
    """One evenly spaced hue per junction (``node_index``), as hex strings."""
    unique = sorted({int(v) for v in node_indices})
    n = max(len(unique), 1)
    return {node: "#%02x%02x%02x" % tuple(int(round(255 * c))
                                          for c in colorsys.hsv_to_rgb(i / n, 1.0, 1.0))
            for i, node in enumerate(unique)}


def ridgeline_ticks(display_distance, x, y, *, interval: float = 1000.0,
                    half_length: float = 50.0, step: float = 10.0) -> list:
    """Perpendicular ticks along a ridge at multiples of ``interval``.

    Positions are interpolated along the displayed distance; the direction is
    taken over +-``step`` metres of it. Distances and coordinates share units
    (metres).

    Returns:
        ``[((x0, y0), (x1, y1)), ...]``.
    """
    d = np.asarray(display_distance, dtype=float)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ticks = []
    for t in np.arange(np.ceil(d.min() / interval) * interval, d.max() + 1.0, interval):
        px, py = np.interp(t, d, x), np.interp(t, d, y)
        a, b = max(d.min(), t - step), min(d.max(), t + step)
        dx = np.interp(b, d, x) - np.interp(a, d, x)
        dy = np.interp(b, d, y) - np.interp(a, d, y)
        length = float(np.hypot(dx, dy))
        if length == 0:
            continue
        nx, ny = -dy / length, dx / length
        ticks.append(((px - nx * half_length, py - ny * half_length),
                      (px + nx * half_length, py + ny * half_length)))
    return ticks


def _iter_parts(geom):
    if geom is None or geom.is_empty:
        return
    if geom.geom_type in ("LineString", "Point"):
        yield geom
    elif hasattr(geom, "geoms"):
        for part in geom.geoms:
            yield from _iter_parts(part)


def render_ridgeline_map(hillshade, profile, *, highlighted_tdj: Optional[dict] = None,
                         overlay=None, reverse: bool = False,
                         distance_range: Optional[tuple] = None, profile_distance=None,
                         shade: bool = False, pad_m: float = 100.0,
                         tick_interval: float = 1000.0, tick_half_length: float = 50.0,
                         coordinates: str = "local", style: Optional[MapStyle] = None,
                         north_arrow: bool = True, fig=None):
    """Map one ridgeline over a hillshade, with its junction overlay.

    Args:
        hillshade: Hillshade raster (path or DataArray), e.g. lsdtt3's
            ``_hillshade.tif``; or an elevation grid with ``shade=True``.
        profile: One ridgeline after
            :func:`~lsdviztools3.io.ridgeline.select_profile` (needs
            ``latitude``, ``longitude`` and ``display_distance``).
        highlighted_tdj: The highlighted TDJ set
            (:func:`~lsdviztools3.io.ridgeline.load_tdj`); its junctions within
            the displayed range are circled. Its ``kind`` picks the overlay
            style: ``native`` draws ``overlay`` as a blue channel network;
            ``topological`` keeps the overlay features with this set's
            ``threshold_m`` and visible ``node_index`` values and colours them
            per junction.
        overlay: Channel network (native) or paired flowpaths (topological)
            vector, path or GeoDataFrame; None for none.
        reverse, distance_range: As given to ``select_profile`` (placing the
            highlighted junctions).
        profile_distance: ``distance`` of the whole, unclipped ridgeline, the
            anchor of the reverse rule (default: ``profile["distance"]``, right
            unless ``distance_range`` clipped it).
        shade: Treat ``hillshade`` as elevation and shade it.
        pad_m: Margin around the ridge and junctions (m), clipped to the raster.
        tick_interval: Distance between perpendicular ridge ticks (m).
        tick_half_length: Half-length of a tick (m).
        coordinates: Map coordinates (``io.raster.to_coordinates``); must be
            projected.
        style: Map style. Defaults to a grey colour map.
        north_arrow: Draw the true-north arrow.
        fig: Existing figure, or None.

    Returns:
        The ``pygmt.Figure``.
    """
    from .decorations import add_bar_scalebar, add_north_arrow
    from .hillshade import render_hillshade

    style = style or MapStyle(cmap="gray")
    for column in ("latitude", "longitude", "display_distance"):
        if column not in profile.columns:
            raise ValueError(f"Profile is missing column {column!r}")
    grid = hillshade if hasattr(hillshade, "rio") else load_raster(hillshade)
    grid = to_coordinates(grid, coordinates)
    crs = get_crs(grid)
    if not is_projected(crs):
        raise ValueError("The ridgeline map needs projected coordinates (metres)")

    rx, ry = _lonlat_to_crs(profile["longitude"], profile["latitude"], crs)
    d = profile["display_distance"].to_numpy(dtype=float)
    xs, ys = [rx], [ry]

    points = None
    if highlighted_tdj is not None:
        full = profile_distance if profile_distance is not None else profile["distance"]
        points = tdj_rows_in_range(highlighted_tdj, full, reverse, distance_range)
        if len(points):
            px, py = _lonlat_to_crs(points["longitude"], points["latitude"], crs)
            xs.append(px)
            ys.append(py)

    gx = np.asarray(grid.x.values, dtype=float)
    gy = np.asarray(grid.y.values, dtype=float)
    allx, ally = np.concatenate(xs), np.concatenate(ys)
    west = max(float(allx.min()) - pad_m, float(gx.min()))
    east = min(float(allx.max()) + pad_m, float(gx.max()))
    south = max(float(ally.min()) - pad_m, float(gy.min()))
    north = min(float(ally.max()) + pad_m, float(gy.max()))
    if west >= east or south >= north:
        raise ValueError("The ridgeline does not overlap the hillshade raster")
    xsel = slice(west, east) if gx[0] < gx[-1] else slice(east, west)
    ysel = slice(south, north) if gy[0] < gy[-1] else slice(north, south)
    grid = grid.sel(x=xsel, y=ysel)

    fig = render_hillshade(grid, coordinates=coordinates, shade=shade, style=style, fig=fig,
                           decorate=False)
    extent = [west / KM, east / KM, south / KM, north / KM]

    if overlay is not None and highlighted_tdj is not None:
        gdf = overlay if hasattr(overlay, "geometry") else load_vector(overlay)
        gdf = gdf.to_crs(crs)
        if "ridgeline_number" in gdf.columns and "ridgeline_number" in profile.columns:
            gdf = gdf[gdf["ridgeline_number"].astype(str)
                      == str(int(profile["ridgeline_number"].iloc[0]))]
        topological = highlighted_tdj["kind"] == "topological"
        if topological:
            if "threshold_m" in gdf.columns:
                gdf = gdf[gdf["threshold_m"].astype(float) == float(highlighted_tdj["threshold"])]
            if "node_index" in gdf.columns and points is not None and "node_index" in points:
                gdf = gdf[gdf["node_index"].astype(int).isin(points["node_index"].astype(int))]
        else:
            gdf = gdf.cx[west:east, south:north]
        colours = (junction_colours(gdf["node_index"])
                   if topological and "node_index" in gdf.columns else {})
        for _, feature in gdf.iterrows():
            if topological:
                colour = colours.get(int(feature["node_index"]), TDJ_RED) if colours else TDJ_RED
                pen = f"0.6p,{colour}@15"
            else:
                pen = CHANNEL_PEN
            for part in _iter_parts(feature.geometry):
                px, py = (np.asarray(v, dtype=float) / KM for v in part.xy)
                if part.geom_type == "Point":
                    fig.plot(x=px, y=py, style="c0.05c", fill=pen.split(",", 1)[1])
                else:
                    fig.plot(x=px, y=py, pen=pen)

    for (x0, y0), (x1, y1) in ridgeline_ticks(d, rx, ry, interval=tick_interval,
                                              half_length=tick_half_length):
        fig.plot(x=[x0 / KM, x1 / KM], y=[y0 / KM, y1 / KM], pen=RIDGE_PEN)
    fig.plot(x=rx / KM, y=ry / KM, pen=RIDGE_HALO_PEN)
    fig.plot(x=rx / KM, y=ry / KM, pen=RIDGE_PEN)
    if points is not None and len(points):
        px, py = _lonlat_to_crs(points["longitude"], points["latitude"], crs)
        fig.plot(x=px / KM, y=py / KM, style="c0.25c", pen=f"1p,{TDJ_RED}")

    if north_arrow:
        add_north_arrow(fig, extent, crs, projection=style.projection)
    add_bar_scalebar(fig, extent, crs)
    return fig
