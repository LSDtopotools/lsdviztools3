"""Named demo regions (bounding boxes) for examples and tests.

Each ``bbox`` is ``(west, south, east, north)`` in WGS84 degrees. These drive
keyless Microsoft Planetary Computer (Copernicus GLO-30) downloads, so no DEM
data is committed to the repository.

Small single-tile regions (Eigg) are used for the fast test suite; the larger,
high-relief tropical islands are showcase examples that span several GLO-30
tiles (merged automatically on download).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple


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
