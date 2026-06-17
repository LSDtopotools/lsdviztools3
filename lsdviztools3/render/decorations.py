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


def add_colourbar(
    fig,
    *,
    label: Optional[str] = None,
    position: str = "JBC+o0c/1.0c+w8c/0.35c+h",
    frame=None,
    cmap=True,
):
    """Add a colourbar for the active CPT.

    Call after a ``makecpt`` (e.g. ``render_channels(color_by=...)`` leaves the
    log10 k_sn palette active). Default placement is a horizontal bar centred
    below the map.

    Args:
        fig: The ``pygmt.Figure``.
        label: Axis label (e.g. ``"log@-10@-k@-sn@-"``). Becomes ``x+l<label>``.
        position: GMT ``-D`` position string.
        frame: Explicit colourbar frame; overrides ``label`` when given.
        cmap: CPT to use (``True`` = the active one).

    Returns:
        The same ``pygmt.Figure``.
    """
    if frame is None:
        frame = f'x+l{label}' if label else "af"
    fig.colorbar(cmap=cmap, position=position, frame=frame)
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
