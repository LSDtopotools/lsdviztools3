"""Map decorations: scalebar, north arrow, credit line, colourbar, inset.

``add_colourbar`` wraps ``Figure.colorbar`` for the active CPT (e.g. the log10
k_sn palette set by ``render_channels``). ``add_scalebar`` draws a scale bar that
works on both the projected (km-linear) maps this package uses and geographic
(Mercator) maps — see :func:`add_scalebar` for why the two cases differ.
``add_bar_scalebar`` / ``add_north_arrow`` / ``add_map_credit`` are the minimal
"clean" furniture: a solid bar, a plain arrowhead with an N, and a small data
credit, all drawn with ``fig.plot``/``fig.text`` in display coordinates.
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
    label_position: str = "above",
    label_halo: Optional[str] = None,
):
    """Draw a small solid-rectangle scale bar (the "clean" map scale bar).

    The bar length is the largest :data:`NICE_SCALEBAR_KM` value not exceeding
    ``max_fraction`` of the map width (default 1/10), so it stays compact and
    lands on a round number (e.g. 5 km, 25 km). A ``"<n> km"`` label sits above
    the bar (default) or below it (``label_position="below"``).

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
        label_position: ``"above"`` (default, the original layout) or
            ``"below"`` the bar. With ``"below"`` on a bottom anchor the bar is
            lifted so the label keeps the same inset from the map edge.
        label_halo: Colour of a thin outline drawn *behind* the label glyphs
            (GMT ``=~`` font outline), so linework under the label (e.g. a
            channel overlay) cannot run through the text. ``None`` (default)
            keeps the original plain label.

    Returns:
        The same ``pygmt.Figure``.
    """
    if label_position not in ("above", "below"):
        raise ValueError(
            f"label_position must be 'above' or 'below', got {label_position!r}.")
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
    gap = 0.012 * height
    # Room for a below-bar label (~9 pt text) when anchored at the bottom.
    label_room = 0.05 * height if (label_position == "below" and "B" in anchor) else 0.0
    x0 = west + px if "L" in anchor else east - px - bar_len
    y0 = south + py + label_room if "B" in anchor else north - py - bh

    fig.plot(
        x=[x0, x0 + bar_len, x0 + bar_len, x0, x0],
        y=[y0, y0, y0 + bh, y0 + bh, y0],
        fill=color,
        pen=f"0.5p,{color}",
        no_clip=True,
    )
    label = f"{bar_km:g} km"
    if label_position == "below":
        label_y, justify = y0 - gap, "CT"
    else:
        label_y, justify = y0 + bh + gap, "CB"
    fig.text(
        x=x0 + bar_len / 2.0,
        y=label_y,
        text=label,
        justify=justify,
        font=_halo_font(font, label_halo),
        no_clip=True,
    )
    return fig


def _halo_font(font: str, halo: Optional[str], width: str = "2p") -> str:
    """GMT font string with an outline drawn beneath the fill (``=~pen``)."""
    return f"{font}=~{width},{halo}" if halo else font


#: GMT projection codes whose +y axis is grid north (when not flipped).
_NORTH_UP_PROJECTIONS = ("X", "x", "M", "m")


def _check_north_up(crs, projection: Optional[str]) -> None:
    """Raise unless the map's +y direction is grid north.

    Two things can break "up = north": the GMT projection (a rotated, polar or
    y-flipped projection) and the CRS itself (some projected CRSs have southing
    or westing axes). Both are checked; anything not provably north-up raises.
    """
    if projection is not None:
        code = projection.strip()[:1]
        if code not in _NORTH_UP_PROJECTIONS:
            raise ValueError(
                f"add_north_arrow: projection {projection!r} is not a north-up "
                "linear (X) or Mercator (M) projection; grid north is not "
                "straight up, so no arrow is drawn.")
        # A negative y scale (e.g. "X12c/-8c") flips the y axis.
        if code in ("X", "x") and "/" in projection:
            if projection.split("/", 1)[1].strip().startswith("-"):
                raise ValueError(
                    f"add_north_arrow: projection {projection!r} flips the y "
                    "axis; grid north points down.")
    if crs is None or not is_projected(crs):
        return  # geographic: +y is latitude, i.e. north
    from pyproj import CRS as _PyprojCRS

    pcrs = _PyprojCRS.from_user_input(crs.to_wkt() if hasattr(crs, "to_wkt") else crs)
    dirs = [a.direction.lower() for a in pcrs.axis_info[:2]]
    # Polar stereographic axes are named for the meridian they run along
    # ("south" / "south" for EPSG:3413); the arrow's angle comes from the
    # convergence there, so only a reversed (westing/southing) CRS is refused.
    polar = dirs[0] == dirs[1] and dirs[0] in ("north", "south")
    if sorted(dirs) != ["east", "north"] and not polar:
        raise ValueError(
            f"add_north_arrow: CRS axes are {dirs}, not east/north; the "
            "map's +y direction is not grid north.")


def north_arrow_geometry(extent, crs=None, *, anchor: str = "TR",
                         size_frac: float = 0.07, pad_frac: float = 0.04,
                         angle_degrees: float = 0.0):
    """Vertices of the north arrowhead and the position of its ``N`` label.

    Returned in display units: ``(xs, ys, (label_x, label_y))`` where
    ``xs``/``ys`` close a triangle whose apex points ``angle_degrees``
    clockwise from +y (0: straight up, grid north), turned about its centre.
    The ``N`` sits centred beyond the base, on the arrow's axis.
    """
    west, east, south, north = extent
    width, height = east - west, north - south
    # Size from the shorter map side, so tall maps do not get a giant arrow.
    h = size_frac * min(width, height)
    w = 0.6 * h
    if not is_projected(crs) and crs is not None:
        # Degrees of longitude are shorter on the ground; keep the head's shape.
        lat = 0.5 * (south + north)
        w = w / max(math.cos(math.radians(lat)), 1e-6)
    label_h = 0.045 * height
    px, py = pad_frac * width, pad_frac * height
    cx = west + px + w / 2.0 if "L" in anchor else east - px - w / 2.0
    if "T" in anchor:
        apex_y = north - py
        base_y = apex_y - h
    else:
        base_y = south + py + label_h
        apex_y = base_y + h
    xs = [cx - w / 2.0, cx, cx + w / 2.0, cx - w / 2.0]
    ys = [base_y, apex_y, base_y, base_y]
    label = (cx, base_y - 0.01 * height)
    if angle_degrees:
        # Turn clockwise by angle about the arrow's centre (equal x/y scale:
        # projected maps are drawn in km on both axes).
        cy = 0.5 * (apex_y + base_y)
        a = math.radians(angle_degrees)
        ca, sa = math.cos(a), math.sin(a)

        def turn(x, y):
            dx, dy = x - cx, y - cy
            return cx + dx * ca + dy * sa, cy - dx * sa + dy * ca

        pts = [turn(x, y) for x, y in zip(xs, ys)]
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        label = turn(cx, cy - 0.5 * h - 0.01 * height)
    return xs, ys, label


def add_north_arrow(
    fig,
    extent,
    crs,
    *,
    projection: Optional[str] = None,
    anchor: str = "TR",
    color: str = "black",
    font: str = "9p,Helvetica-Bold,black",
    size_frac: float = 0.07,
    pad_frac: float = 0.04,
    label_halo: Optional[str] = None,
    true_north: bool = True,
):
    """Draw an unornamented north arrow: a solid arrowhead with an ``N`` under it.

    No compass rose, circle or other decoration. On a projected map the arrow
    points to **true north at the arrow's own position**: it is turned from
    grid up by the grid convergence there (``lsdviztools3.frames
    .grid_convergence_degrees``). That is under ~3 degrees within a UTM zone,
    up to 180 degrees in polar stereographic, and 90 minus the azimuth (60
    degrees for azimuth 30) in an lsdtt3 oblique frame, whose grid north is
    not north at all. GMT's own map rose points straight up on these Cartesian
    maps, so it is not used. Geographic (Mercator) maps are north-up.

    If the map is not north-up — a rotated/polar/y-flipped GMT projection, or a
    CRS whose axes are not east/north (e.g. a south-orientated CRS) — this
    raises rather than drawing a wrong arrow.

    Args:
        fig: The ``pygmt.Figure``.
        extent: ``[west, east, south, north]`` in the map's display units.
        crs: The map CRS.
        projection: The GMT projection the map was drawn with (``None`` = the
            package's auto projection, which is north-up).
        anchor: Corner: ``"TR"`` (default), ``"TL"``, ``"BR"``, ``"BL"``.
        color: Arrowhead fill.
        font: Font of the ``N``.
        size_frac: Arrowhead height as a fraction of map height.
        pad_frac: Inset from the chosen corner, as a fraction of map size
            (same convention as :func:`add_bar_scalebar`).
        label_halo: Optional outline colour behind the ``N`` (as for
            :func:`add_bar_scalebar`).
        true_north: Point to true north (default). False points straight up,
            to grid north.

    Returns:
        The same ``pygmt.Figure``.

    Raises:
        ValueError: if the map's +y direction is not grid north.
    """
    _check_north_up(crs, projection)
    angle = 0.0
    if true_north and crs is not None and is_projected(crs):
        from ..frames import grid_convergence_degrees

        xs0, ys0, _ = north_arrow_geometry(
            extent, crs, anchor=anchor, size_frac=size_frac, pad_frac=pad_frac)
        # Display units are km on projected maps.
        angle = grid_convergence_degrees(crs, 1000.0 * sum(xs0[:3]) / 3.0,
                                         1000.0 * sum(ys0[:3]) / 3.0)
    xs, ys, (lx, ly) = north_arrow_geometry(
        extent, crs, anchor=anchor, size_frac=size_frac, pad_frac=pad_frac,
        angle_degrees=angle)
    fig.plot(x=xs, y=ys, fill=color, pen=f"0.5p,{color}", no_clip=True)
    fig.text(x=lx, y=ly, text="N", justify="CT" if abs(angle) < 10.0 else "CM",
             font=_halo_font(font, label_halo), no_clip=True)
    return fig


def add_map_credit(
    fig,
    extent,
    text: str,
    *,
    anchor: str = "BR",
    font: str = "6p,Helvetica,gray20",
    box_fill: Optional[str] = "white@35",
    pad_frac: float = 0.01,
):
    """Write a small, unobtrusive data-credit line in a map corner.

    Used for the geologic-map / DEM citation (see
    :func:`lsdviztools3.io.lithology.lithology_map_credit`). The text is drawn
    inside the map on a faint white box so it stays legible over any colour.

    Args:
        fig: The ``pygmt.Figure``.
        extent: ``[west, east, south, north]`` in display units.
        text: The credit line (empty/None draws nothing). ``"\n"`` separates
            several lines, stacked away from the anchored map edge.
        anchor: Corner: ``"BR"`` (default), ``"BL"``, ``"TR"``, ``"TL"``.
        font: Font. 6 pt by default — deliberately below the 8 pt floor used
            for map annotation, because a full credit line at 8 pt spans most
            of a single-column map. Raise it if a journal requires.
        box_fill: Fill behind the text (``None`` for none).
        pad_frac: Inset from the corner as a fraction of map size.

    Returns:
        The same ``pygmt.Figure``.
    """
    if not text:
        return fig
    west, east, south, north = extent
    px, py = pad_frac * (east - west), pad_frac * (north - south)
    x = west + px if "L" in anchor else east - px
    y = south + py if "B" in anchor else north - py
    just = ("L" if "L" in anchor else "R") + ("B" if "B" in anchor else "T")
    kwargs = {"fill": box_fill, "clearance": "0.05c/0.03c"} if box_fill else {}
    # Several lines ("\n"-separated) stack away from the anchored edge.
    lines = [ln for ln in str(text).split("\n") if ln.strip()]
    pt = float(font.split(",")[0].rstrip("p")) if font.split(",")[0].rstrip("p") else 6.0
    step = 1.45 * pt * 0.0352778  # line spacing in cm
    order = lines[::-1] if "B" in anchor else lines
    for i, line in enumerate(order):
        dy = i * step if "B" in anchor else -i * step
        extra = {"offset": f"0c/{dy:.3f}c"} if i else {}
        # GMT text escapes: a literal "@" must be written "@@".
        fig.text(x=x, y=y, text=line.replace("@", "@@"), justify=just,
                 font=font, no_clip=True, **kwargs, **extra)
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
