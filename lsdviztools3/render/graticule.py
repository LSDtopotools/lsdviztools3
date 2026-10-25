"""Reference-grid graticule over an oblique frame map.

A raster in an lsdtt3 oblique frame is drawn as a rectangle in frame km, with
the feature along +x. To locate it, the easting/northing grid of a reference
CRS (the UTM, or polar stereographic, CRS of the frame centre; see
``lsdviztools3.frames.reference_crs_for``) is drawn across it: lines of
constant easting and of constant northing, densified, transformed into the
frame and clipped to the map, at a "nice" km interval, each labelled
horizontally where it meets the map edge. Maps already in UTM or polar
coordinates need no graticule; their own axes are that grid.
"""

from __future__ import annotations

import math
from typing import List, Optional

import numpy as np

from ..frames import as_pyproj, is_oblique_frame, reference_crs_for
from .base import KM
from .style import MapStyle

#: Points per graticule line before clipping.
LINE_POINTS = 200


def nice_interval(span: float, target_lines: int = 4) -> float:
    """A 1, 2 or 5 x 10^k interval giving about ``target_lines`` lines over ``span``."""
    if not (span > 0.0):
        raise ValueError("span must be positive")
    raw = span / max(target_lines, 1)
    power = 10.0 ** math.floor(math.log10(raw))
    for step in (1.0, 2.0, 5.0, 10.0):
        if raw <= step * power:
            return step * power
    return 10.0 * power


def graticule_lines(extent_m, crs, reference_crs=None, interval_m: Optional[float] = None,
                    target_lines: int = 4) -> List[dict]:
    """Reference-grid lines clipped to a frame map, in frame metres.

    Args:
        extent_m: ``[west, east, south, north]`` of the map in the frame CRS (m).
        crs: The frame CRS.
        reference_crs: Override of the reference CRS (default: the analysis CRS
            of the frame centre).
        interval_m: Line spacing in reference-CRS metres; default a 1-2-5 step
            giving about ``target_lines`` lines across the longer reference span.

    Returns:
        One dict per clipped piece: ``axis`` ("E" or "N"), ``value`` (the
        reference coordinate, m), ``x``/``y`` (vertex arrays, frame m), and
        ``label`` (``(x, y, edge)`` where the piece meets the map edge, edge in
        "bottom", "left", "top", "right"), or ``label`` None.
    """
    from pyproj import Transformer
    from shapely.geometry import LineString, box

    west, east, south, north = (float(v) for v in extent_m)
    ref = reference_crs_for(crs, reference_crs)
    if ref is None:
        raise ValueError("graticule_lines: the CRS is not an oblique frame; pass reference_crs")
    frame = as_pyproj(crs)
    to_ref = Transformer.from_crs(frame, ref, always_xy=True)
    to_frame = Transformer.from_crs(ref, frame, always_xy=True)

    # The map rectangle in the reference CRS, edges densified.
    t = np.linspace(0.0, 1.0, 50)
    bx = np.concatenate([west + (east - west) * t, np.full_like(t, east),
                         east - (east - west) * t, np.full_like(t, west)])
    by = np.concatenate([np.full_like(t, south), south + (north - south) * t,
                         np.full_like(t, north), north - (north - south) * t])
    re, rn = to_ref.transform(bx, by)
    emin, emax, nmin, nmax = float(re.min()), float(re.max()), float(rn.min()), float(rn.max())
    step = interval_m or nice_interval(max(emax - emin, nmax - nmin), target_lines)

    rect = box(west, south, east, north)
    tol = 1e-6 * max(east - west, north - south)
    lines = []
    for axis, lo, hi, across_lo, across_hi in (("E", emin, emax, nmin, nmax),
                                               ("N", nmin, nmax, emin, emax)):
        for value in np.arange(math.ceil(lo / step) * step, hi + 0.5 * step, step):
            if value > hi:
                break
            s = np.linspace(across_lo, across_hi, LINE_POINTS)
            v = np.full_like(s, value)
            fx, fy = to_frame.transform(v, s) if axis == "E" else to_frame.transform(s, v)
            piece = LineString(np.column_stack([fx, fy])).intersection(rect)
            parts = getattr(piece, "geoms", [piece])
            for part in parts:
                if part.is_empty or part.geom_type != "LineString" or part.length <= 0.0:
                    continue
                xy = np.asarray(part.coords)
                lines.append({"axis": axis, "value": float(value), "x": xy[:, 0], "y": xy[:, 1],
                              "label": _label_point(xy, west, east, south, north, tol, axis)})
    return lines


def _label_point(xy, west, east, south, north, tol, axis):
    """Where a clipped line meets the map edge, preferring bottom (E) or left (N)."""
    order = ("bottom", "left", "top", "right") if axis == "E" else ("left", "bottom", "right", "top")
    on_edge = {
        "bottom": lambda x, y: abs(y - south) <= tol,
        "left": lambda x, y: abs(x - west) <= tol,
        "top": lambda x, y: abs(y - north) <= tol,
        "right": lambda x, y: abs(x - east) <= tol,
    }
    ends = [tuple(xy[0]), tuple(xy[-1])]
    for edge in order:
        for x, y in ends:
            if on_edge[edge](x, y):
                return (float(x), float(y), edge)
    return None


def graticule_label(axis: str, value_m: float) -> str:
    """Label text: the reference coordinate in km, e.g. ``"512 km E"``."""
    return f"{value_m / KM:g} km {axis}"


#: GMT justification and offset of a label, by the edge it sits on (inside the map).
_LABEL_PLACEMENT = {
    "bottom": ("CB", "0c/0.08c"),
    "top": ("CT", "0c/-0.08c"),
    "left": ("LM", "0.08c/0c"),
    "right": ("RM", "-0.08c/0c"),
}


def add_reference_graticule(fig, extent, crs, style: Optional[MapStyle] = None,
                            reference_crs=None):
    """Draw the reference-CRS easting/northing grid over an oblique frame map.

    A no-op unless ``crs`` is an oblique frame and ``style.graticule`` is on.

    Args:
        fig: The ``pygmt.Figure`` (already holding the map).
        extent: ``[west, east, south, north]`` in display units (frame km).
        crs: The map CRS.
        style: ``graticule``, ``graticule_pen``, ``graticule_font``,
            ``graticule_interval_km`` and ``reference_crs`` are read from it.
        reference_crs: Override of ``style.reference_crs``.

    Returns:
        The same ``pygmt.Figure``.
    """
    style = style or MapStyle()
    if not style.graticule or not is_oblique_frame(crs):
        return fig
    extent_m = [float(v) * KM for v in extent]
    interval = style.graticule_interval_km * KM if style.graticule_interval_km else None
    lines = graticule_lines(extent_m, crs, reference_crs or style.reference_crs, interval)
    for line in lines:
        fig.plot(x=line["x"] / KM, y=line["y"] / KM, pen=style.graticule_pen)
    for line in lines:
        if line["label"] is None:
            continue
        x, y, edge = line["label"]
        justify, offset = _LABEL_PLACEMENT[edge]
        fig.text(x=x / KM, y=y / KM, text=graticule_label(line["axis"], line["value"]),
                 justify=justify, offset=offset, font=f"{style.graticule_font}=~1.5p,white",
                 no_clip=True)
    return fig
