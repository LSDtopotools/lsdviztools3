"""The hillshade illumination convention.

This existed as a real defect: ``_hillshade`` used ``aspect = atan2(-dx, dy)``,
exactly 180 degrees from GDAL's convention, so every hillshade lsdviztools3 drew
was inverted -- ridges read as valleys. Nothing caught it because an inverted
hillshade is a perfectly plausible-looking image.

The test is phrased as a statement about light rather than about the formula:
under a sun in the north-west, the north-west flank of a hill is brighter than
the south-east flank. That holds for any correct implementation and fails for
the inverted one.
"""

import numpy as np
import pytest

from lsdviztools3.charts.hilltop_coupling import _hillshade


def _peak(n=101, sigma=20.0):
    """A single Gaussian peak. Row 0 is NORTH, as in a north-up raster."""
    y, x = np.mgrid[0:n, 0:n]
    return 100.0 * np.exp(-(((x - n // 2) ** 2 + (y - n // 2) ** 2)
                            / (2.0 * sigma ** 2)))


@pytest.mark.parametrize("azimuth,bright,dark", [
    (315.0, (slice(30, 45), slice(30, 45)), (slice(55, 70), slice(55, 70))),
    (135.0, (slice(55, 70), slice(55, 70)), (slice(30, 45), slice(30, 45))),
])
def test_hillshade_illuminates_the_flank_facing_the_sun(azimuth, bright, dark):
    shade = _hillshade(_peak(), cell=1.0, azimuth=azimuth)
    assert shade[bright].mean() > shade[dark].mean() + 0.2


def test_hillshade_is_not_inverted_relative_to_gdal():
    """Directly against GDAL's documented Horn formulation."""
    z = _peak()
    cell = 1.0
    dy, dx = np.gradient(z, cell)
    slope = np.arctan(np.hypot(dx, dy))
    gdal_aspect = np.arctan2(dy, -dx)
    az = np.radians(360.0 - 315.0 + 90.0)
    alt = np.radians(45.0)
    expected = np.clip(np.sin(alt) * np.cos(slope)
                       + np.cos(alt) * np.sin(slope) * np.cos(az - gdal_aspect),
                       0.0, 1.0)
    assert np.allclose(_hillshade(z, cell), expected)


def test_flat_ground_takes_the_sine_of_the_solar_altitude():
    shade = _hillshade(np.zeros((20, 20)), cell=1.0, altitude=30.0)
    assert np.allclose(shade, np.sin(np.radians(30.0)))
