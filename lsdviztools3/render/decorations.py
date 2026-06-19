"""Map decorations: scalebar, colourbar, inset.

``add_colourbar`` wraps ``Figure.colorbar`` for the active CPT (e.g. the log10
k_sn palette set by ``render_channels``). ``add_scalebar`` draws a scale bar that
works on both the projected (km-linear) maps this package uses and geographic
(Mercator) maps — see :func:`add_scalebar` for why the two cases differ.
"""

from __future__ import annotations

import math
from typing import Optional

from .base import is_projected

_DEFERRED_MSG = (
    "lsdviztools3.render.decorations.add_inset is deferred in v1 and not "
    "implemented yet."
)

#: GMT symbol size: centimetres per point.
_CM_PER_PT = 0.0352778


#: Bar dimensions (length/thickness) shared by all placements.
_BAR_SIZE = "+w8c/0.35c"

#: Named colourbar placements -> GMT ``-D`` position strings. Offsets are chosen
#: to clear the map's axis labels (e.g. "Easting (km)") with a comfortable gap.
#: Horizontal placements add ``+h``; vertical ones omit it.
_PLACEMENTS = {
    "bottom": "JBC+o0c/1.3c" + _BAR_SIZE + "+h",
    "top": "JTC+o0c/1.4c" + _BAR_SIZE + "+h",
    "right": "JMR+o1.4c/0c" + _BAR_SIZE,
    "left": "JML+o1.4c/0c" + _BAR_SIZE,
}
#: Friendly aliases.
_PLACEMENTS["below"] = _PLACEMENTS["bottom"]
_PLACEMENTS["above"] = _PLACEMENTS["top"]


def colourbar_position(placement: str) -> str:
    """Resolve a named placement (or a raw GMT string) to a ``-D`` position.

    Args:
        placement: One of ``bottom``/``below``, ``top``/``above``, ``left``,
            ``right``, or a raw GMT ``-D`` position string (used as-is if it
            isn't a known name).

    Returns:
        A GMT ``-D`` position string.
    """
    return _PLACEMENTS.get(placement.lower(), placement)


def add_colourbar(
    fig,
    *,
    label: Optional[str] = None,
    placement: str = "bottom",
    position: Optional[str] = None,
    frame=None,
    cmap=True,
):
    """Add a colourbar for the active CPT.

    Call after a ``makecpt`` (e.g. ``render_channels(color_by=...)`` leaves the
    log10 k_sn palette active). Default placement is a horizontal bar centred
    below the map, offset to clear the x-axis label.

    Args:
        fig: The ``pygmt.Figure``.
        label: Axis label (e.g. ``"log@-10@-k@-sn@-"``). Becomes ``x+l<label>``.
        placement: Named placement: ``"bottom"`` (default), ``"top"``,
            ``"left"``, ``"right"`` (``"below"``/``"above"`` also accepted).
            Ignored if ``position`` is given.
        position: Explicit GMT ``-D`` position string; overrides ``placement``.
        frame: Explicit colourbar frame; overrides ``label`` when given.
        cmap: CPT to use (``True`` = the active one).

    Returns:
        The same ``pygmt.Figure``.
    """
    if position is None:
        position = colourbar_position(placement)
    if frame is None:
        frame = f'x+l{label}' if label else "af"
    fig.colorbar(cmap=cmap, position=position, frame=frame)
    return fig


#: Preferred scale-bar lengths in km (used by the "clean" solid scale bar).
NICE_SCALEBAR_KM = (
    0.1, 0.5, 1, 2, 3, 4, 5, 8, 10, 15, 20, 25, 40, 50, 75, 100, 150, 200,
)


def _nice_bar_km(max_km: float) -> float:
    """Largest preferred length (NICE_SCALEBAR_KM) not exceeding ``max_km``."""
    options = [v for v in NICE_SCALEBAR_KM if v <= max_km]
    return float(options[-1]) if options else float(NICE_SCALEBAR_KM[0])


def add_bar_scalebar(
    fig,
    extent,
    crs,
    *,
    max_fraction: float = 0.1,
    anchor: str = "BL",
    color: str = "black",
    font: str = "9p,Helvetica,black",
    height_frac: float = 0.018,
    pad_frac: float = 0.04,
):
    """Draw a small solid-rectangle scale bar (the "clean" map scale bar).

    The bar length is the largest :data:`NICE_SCALEBAR_KM` value not exceeding
    ``max_fraction`` of the map width (default 1/10), so it stays compact and
    lands on a round number (e.g. 5 km, 25 km). A ``"<n> km"`` label sits above
    the bar.

    Args:
        fig: The ``pygmt.Figure``.
        extent: ``[west, east, south, north]`` in the map's display units (km for
            projected, degrees for geographic).
        crs: The map CRS (selects the km<->display-unit conversion).
        max_fraction: Cap on bar length as a fraction of map width.
        anchor: Corner: ``"BL"`` (default), ``"BR"``, ``"TL"``, ``"TR"``.
        color: Fill/outline colour of the bar and label.
        font: Label font.
        height_frac: Bar height as a fraction of map height.
        pad_frac: Inset from the chosen corner, as a fraction of map size.

    Returns:
        The same ``pygmt.Figure``.
    """
    west, east, south, north = extent
    width = east - west
    height = north - south

    if is_projected(crs):
        width_km = width  # display units are already km
    else:
        lat = 0.5 * (south + north)
        width_km = 111.319 * width * math.cos(math.radians(lat))

    bar_km = _nice_bar_km(width_km * max_fraction)
    if is_projected(crs):
        bar_len = bar_km
    else:
        lat = 0.5 * (south + north)
        bar_len = bar_km / (111.319 * math.cos(math.radians(lat)))

    px, py = pad_frac * width, pad_frac * height
    bh = height_frac * height
    x0 = west + px if "L" in anchor else east - px - bar_len
    y0 = south + py if "B" in anchor else north - py - bh

    fig.plot(
        x=[x0, x0 + bar_len, x0 + bar_len, x0, x0],
        y=[y0, y0, y0 + bh, y0 + bh, y0],
        fill=color,
        pen=f"0.5p,{color}",
        no_clip=True,
    )
    label = f"{bar_km:g} km"
    fig.text(
        x=x0 + bar_len / 2.0,
        y=y0 + bh + 0.012 * height,
        text=label,
        justify="CB",
        font=font,
        no_clip=True,
    )
    return fig


def _nice_length(x: float) -> float:
    """Largest 'nice' number (1/2/5 x 10^n) not exceeding x."""
    if x <= 0:
        return 1.0
    exp = math.floor(math.log10(x))
    base = 10.0 ** exp
    for m in (5.0, 2.0, 1.0):
        if x >= m * base:
            return m * base
    return base


def add_scalebar(
    fig,
    extent,
    crs,
    *,
    length_km: Optional[float] = None,
    anchor: str = "BL",
    font: str = "9p,Helvetica,black",
    pen: str = "1.2p,black",
):
    """Add a scale bar, branching on whether the map CRS is projected.

    Projected maps here are drawn with a GMT *linear* projection in km
    (1 axis unit = 1 km), where GMT's geographic ``map_scale`` does not apply —
    so the bar is drawn directly as a ruler (a line + end ticks + label) in the
    km coordinate system. Geographic (Mercator) maps use the native
    ``Figure.basemap(map_scale=...)`` so GMT converts to true ground distance.

    Args:
        fig: The ``pygmt.Figure``.
        extent: ``[west, east, south, north]`` in the map's display units
            (km for projected, degrees for geographic).
        crs: The map CRS (decides which scale-bar mechanism is used).
        length_km: Bar length in km. ``None`` -> a 'nice' ~1/5 of map width.
        anchor: Corner: ``"BL"`` (default), ``"BR"``, ``"TL"``, ``"TR"``.
        font: Label font.
        pen: Bar/tick pen.

    Returns:
        The same ``pygmt.Figure``.
    """
    west, east, south, north = extent
    width = east - west
    height = north - south

    if is_projected(crs):
        if length_km is None:
            length_km = _nice_length(width / 5.0)
        mx, my = 0.06 * width, 0.06 * height
        x0 = west + mx if "L" in anchor else east - mx - length_km
        y0 = south + my if "B" in anchor else north - my
        tick = 0.015 * height

        # Horizontal bar.
        fig.plot(x=[x0, x0 + length_km], y=[y0, y0], pen=pen, no_clip=True)
        # End ticks.
        for xt in (x0, x0 + length_km):
            fig.plot(x=[xt, xt], y=[y0 - tick, y0 + tick], pen=pen, no_clip=True)
        # Label.
        label = (f"{length_km:g} km" if length_km >= 1
                 else f"{length_km * 1000:g} m")
        fig.text(x=x0 + length_km / 2.0, y=y0 + 1.8 * tick, text=label,
                 justify="CB", font=font, no_clip=True)
        return fig

    # Geographic (Mercator): native map scale (true ground distance).
    if length_km is None:
        length_km = _nice_length(111.0 * width / 5.0)
    just = {"BL": "LB", "BR": "RB", "TL": "LT", "TR": "RT"}.get(anchor, "LB")
    lon_ref = west if "L" in anchor else east
    lat_ref = south if "B" in anchor else north
    map_scale = f"g{lon_ref}/{lat_ref}+w{length_km:g}k+j{just}+o0.4c/0.4c+f+u"
    fig.basemap(map_scale=map_scale)
    return fig


def add_inset(fig, **kwargs):
    """Add a locator inset map to a figure. Deferred in v1."""
    raise NotImplementedError(_DEFERRED_MSG)
