"""Named demo regions (bounding boxes) for examples and tests.

Each ``bbox`` is ``(west, south, east, north)`` in WGS84 degrees. These drive
keyless Microsoft Planetary Computer (Copernicus GLO-30) downloads, so no DEM
data is committed to the repository.

Small single-tile regions (Eigg) are used for the fast test suite; the larger,
high-relief tropical islands are showcase examples that span several GLO-30
tiles (merged automatically on download).

It also reads boxes typed or pasted by a user: :func:`parse_lat_lon` reads one
``"latitude, longitude"`` pair as Google Maps gives it (right-click a place and
click the numbers at the top of the menu to copy them), and
:func:`box_from_corners` turns two opposite corners into a box. The notebooks
use these so a box takes two right-clicks and two pastes.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import NamedTuple, Tuple


@dataclass(frozen=True)
class Region:
    """A named demo region."""

    name: str
    bbox: Tuple[float, float, float, float]  # (west, south, east, north)
    note: str


REGIONS = {
    "eigg": Region(
        "Isle of Eigg, Scotland",
        (-6.20, 56.84, -6.08, 56.96),
        "Tiny single-tile island (An Sgurr, 393 m); the fast test fixture.",
    ),
    "rum": Region(
        "Isle of Rum, Scotland",
        (-6.42, 56.96, -6.26, 57.07),
        "Rugged relief (Rum Cuillin); spans the 57N tile boundary (2 tiles).",
    ),
    "viti_levu": Region(
        "Viti Levu, Fiji",
        (177.20, -18.30, 178.75, -17.25),
        "Largest Fijian island; high relief (Tomanivi 1324 m), complex drainage.",
    ),
    "guadalcanal": Region(
        "Guadalcanal, Solomon Islands",
        (159.55, -9.95, 160.90, -9.25),
        "High relief (Popomanaseu 2335 m); dense dendritic drainage.",
    ),
}

#: Region used by the (network-gated) test suite — must stay small.
DEFAULT_TEST_REGION = "eigg"


class Box(NamedTuple):
    """A longitude/latitude box in WGS84 degrees.

    The field order is west, east, south, north, so
    ``west, east, south, north = box_from_corners(...)`` unpacks it. (Note that
    :attr:`Region.bbox` uses the other common order, west, south, east, north;
    :attr:`bbox` gives that.)
    """

    west: float
    east: float
    south: float
    north: float

    @property
    def bbox(self) -> Tuple[float, float, float, float]:
        """The box as ``(west, south, east, north)``."""
        return (self.west, self.south, self.east, self.north)


_NUMBER = re.compile(r"[-+]?\d+(?:\.\d+)?")


def parse_lat_lon(text: str) -> Tuple[float, float]:
    """Read a ``"latitude, longitude"`` pair, as copied from Google Maps.

    Google Maps gives decimal degrees, latitude first, for example
    ``"55.9476273103734, -3.158220317286066"`` (Edinburgh). Anything other than
    the two numbers (spaces, brackets) is ignored. Degrees-minutes-seconds text
    such as ``55°56'51"N`` is not read: it holds more than two numbers, and the
    error message says how to copy the decimal form instead.

    Args:
        text: The pasted text.

    Returns:
        ``(latitude, longitude)`` in decimal degrees.

    Raises:
        ValueError: If the text does not hold exactly two numbers, or the
            latitude is outside -90 to 90, or the longitude outside -180 to 180.
    """
    numbers = _NUMBER.findall(text)
    if len(numbers) != 2:
        raise ValueError(
            f'Could not find a latitude and a longitude in "{text}". Right-click the place '
            'on Google Maps, click the two numbers at the top of the menu to copy them, and '
            'paste them between the quotes, for example "55.9476, -3.1582".')
    latitude = float(numbers[0])
    longitude = float(numbers[1])
    if not -90 <= latitude <= 90:
        raise ValueError(f'"{text}": the first number is the latitude and must be between '
                         '-90 and 90. Google Maps gives the latitude first.')
    if not -180 <= longitude <= 180:
        raise ValueError(f'"{text}": the second number is the longitude and must be '
                         'between -180 and 180.')
    return latitude, longitude


def box_from_corners(lower_left: str, upper_right: str, verbose: bool = True) -> Box:
    """The box with these two opposite corners, each a pasted ``"lat, lon"``.

    The corners are meant to be the lower left (south-west) and upper right
    (north-east), but either pair of opposite corners gives the same box: the
    edges are the minimum and maximum of the two. A box that crosses the
    180° meridian is not supported.

    Args:
        lower_left: One corner, as read by :func:`parse_lat_lon`.
        upper_right: The opposite corner.
        verbose: Print the box and its approximate size in km.

    Returns:
        A :class:`Box` (west, east, south, north).

    Raises:
        ValueError: If a corner cannot be read, or the two corners share a
            latitude or a longitude (the box would have no area).
    """
    lat_1, lon_1 = parse_lat_lon(lower_left)
    lat_2, lon_2 = parse_lat_lon(upper_right)
    if lat_1 == lat_2 or lon_1 == lon_2:
        raise ValueError("The two corners must have different latitudes and different longitudes.")
    box = Box(west=min(lon_1, lon_2), east=max(lon_1, lon_2),
              south=min(lat_1, lat_2), north=max(lat_1, lat_2))
    if verbose:
        # Approximate: degrees of latitude are ~110.57 km, of longitude
        # ~111.32 km times the cosine of the latitude.
        width_km = (box.east - box.west) * 111.32 * math.cos(math.radians((box.south + box.north) / 2))
        height_km = (box.north - box.south) * 110.57
        print(f"Box: west {box.west:.5f}, east {box.east:.5f}, "
              f"south {box.south:.5f}, north {box.north:.5f}")
        print(f"Size: about {width_km:.0f} km east-west by {height_km:.0f} km north-south")
    return box
