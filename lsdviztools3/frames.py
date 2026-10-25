"""Oblique analysis frames (lsdtt3 rotated feature frames) and grid convergence.

lsdtt3 can analyse a DEM in an oblique frame: a Hotine Oblique Mercator CRS
whose central line runs along a fault, river or range, with the feature along
+x and the frame centre at (0, 0). Such a raster is an ordinary north-up
raster in that CRS, so it plots as a rectangle in the Cartesian ``X``
projection. These helpers recognise a frame, pick the reference CRS whose
easting/northing grid is drawn over it (``render.graticule``), and give the
grid convergence used to point the north arrow. No PyGMT imports here.

Conventions follow lsdtt3 (``src/core/lsdtt_true_north.hpp``): the
convergence is the angle of true north measured clockwise from grid up, so a
grid bearing ``g`` is the true bearing ``g - c``. pyproj's
``meridian_convergence`` is the opposite sign.
"""

from __future__ import annotations

import math
from typing import Optional, Tuple

#: Beyond this latitude the reference CRS is polar stereographic (lsdtt3 rule).
POLAR_LATITUDE = 75.0


def as_pyproj(crs):
    """A ``pyproj.CRS`` from a rasterio CRS, pyproj CRS, WKT, EPSG code, etc."""
    from pyproj import CRS

    if crs is None:
        return None
    if isinstance(crs, CRS):
        return crs
    if hasattr(crs, "to_wkt"):
        return CRS.from_user_input(crs.to_wkt())
    return CRS.from_user_input(crs)


def _operation(crs):
    pcrs = as_pyproj(crs)
    if pcrs is None or not pcrs.is_projected:
        return None
    return pcrs.coordinate_operation


def is_oblique_frame(crs) -> bool:
    """True if ``crs`` is an oblique Mercator frame (method name contains "Oblique Mercator")."""
    op = _operation(crs)
    return op is not None and "oblique mercator" in op.method_name.lower()


def _parameter(op, *names) -> Optional[float]:
    for p in op.params:
        if any(n in p.name.lower() for n in names):
            return float(p.value)
    return None


def frame_azimuth(crs) -> Optional[float]:
    """Azimuth of the frame's central line, degrees clockwise from true north; None if not a frame."""
    if not is_oblique_frame(crs):
        return None
    return _parameter(_operation(crs), "azimuth")


def frame_centre(crs) -> Optional[Tuple[float, float]]:
    """(longitude, latitude) of the frame centre; None if not a frame."""
    if not is_oblique_frame(crs):
        return None
    op = _operation(crs)
    return (_parameter(op, "longitude of projection centre", "longitude of origin"),
            _parameter(op, "latitude of projection centre", "latitude of origin"))


def analysis_epsg(longitude: float, latitude: float) -> int:
    """The CRS lsdtt3 analyses a point in: UTM, or NSIDC polar stereographic beyond 75 deg."""
    if latitude > POLAR_LATITUDE:
        return 3413
    if latitude < -POLAR_LATITUDE:
        return 3976
    lon = (longitude + 180.0) % 360.0 - 180.0
    zone = int(math.floor((lon + 180.0) / 6.0)) % 60 + 1
    return (32600 if latitude >= 0.0 else 32700) + zone


def reference_crs_for(crs, reference_crs=None):
    """The reference CRS whose grid is drawn over a frame.

    ``reference_crs`` (anything pyproj reads) wins; otherwise the analysis CRS
    of the frame centre (``analysis_epsg``). None when ``crs`` is not a frame
    and no override is given.
    """
    from pyproj import CRS

    if reference_crs is not None:
        return as_pyproj(reference_crs)
    centre = frame_centre(crs)
    if centre is None:
        return None
    return CRS.from_epsg(analysis_epsg(*centre))


def grid_convergence_degrees(crs, x: float, y: float) -> float:
    """Angle of true north, degrees clockwise from grid up, at (x, y) in ``crs`` units.

    0 for a geographic CRS. Computed by stepping 1e-4 deg north (south within
    a degree of the north pole) through the CRS's own geographic CRS.
    """
    from pyproj import Transformer

    pcrs = as_pyproj(crs)
    if pcrs is None or not pcrs.is_projected:
        return 0.0
    to_geo = Transformer.from_crs(pcrs, pcrs.geodetic_crs, always_xy=True)
    to_grid = Transformer.from_crs(pcrs.geodetic_crs, pcrs, always_xy=True)
    lon, lat = to_geo.transform(x, y)
    step = -1e-4 if lat > 89.0 else 1e-4
    x2, y2 = to_grid.transform(lon, lat + step)
    dx, dy = x2 - x, y2 - y
    if step < 0.0:
        dx, dy = -dx, -dy
    return math.degrees(math.atan2(dx, dy))
