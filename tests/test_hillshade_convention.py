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


# --- the GMT path the maps use (render_hillshade, grdimage shading) -----------
# The tests above cover the matplotlib helper only. Maps draw elevation through
# render_hillshade, where GMT's automatic illumination (shading=True) does the
# shading. Same statement about light: on a symmetric hill, the north-west
# flank renders brighter than the south-east flank. Run in two projected CRSs,
# so a north/south or east/west flip from the display transform would show.

@pytest.mark.parametrize("epsg,x0,y0", [(32630, 500000.0, 4200000.0), (27700, 400000.0, 400000.0)])
def test_gmt_hillshade_lights_the_north_west_flank(tmp_path, epsg, x0, y0):
    pytest.importorskip("pygmt")
    rasterio = pytest.importorskip("rasterio")
    from matplotlib.image import imread
    from rasterio.transform import from_origin

    from lsdviztools3.render.hillshade import render_hillshade
    from lsdviztools3.render.style import MapStyle

    n, cell = 101, 30.0
    dem = tmp_path / "hill.tif"
    with rasterio.open(dem, "w", driver="GTiff", height=n, width=n, count=1, dtype="float32",
                       crs=f"EPSG:{epsg}", transform=from_origin(x0, y0 + n * cell, cell, cell)) as d:
        d.write(_peak(n).astype("float32"), 1)
    style = MapStyle(cmap="gray", transparent=False, clean=True)
    fig = render_hillshade(str(dem), shade=True, style=style, decorate=False)
    png = tmp_path / "hill.png"
    fig.savefig(png, dpi=100, crop=True)

    img = imread(png)[..., :3].mean(axis=2)
    # the map area: the rows and columns that are not the white page
    rows = np.where((img < 0.99).any(axis=1))[0]
    cols = np.where((img < 0.99).any(axis=0))[0]
    map_ = img[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
    h, w = map_.shape

    def patch(fy, fx):
        r, c = int(fy * h), int(fx * w)
        return map_[r - 3:r + 4, c - 3:c + 4].mean()

    # the flanks at 0.35 / 0.65 of the map along the NW-SE diagonal; the hill
    # is symmetric, so elevation colour is the same at both
    north_west, south_east = patch(0.35, 0.35), patch(0.65, 0.65)
    assert north_west > south_east + 0.15, (north_west, south_east)
