"""Tests for lsdviztools3.io.raster."""

import pytest
import xarray as xr
import numpy as np

from lsdviztools3.io import raster


@pytest.mark.network
def test_load_raster_returns_dataarray(small_dem):
    da = raster.load_raster(small_dem)
    assert isinstance(da, xr.DataArray)
    assert da.ndim == 2  # single band squeezed
    assert raster.get_crs(da) is not None


@pytest.mark.network
def test_get_extent_ordering(small_dem):
    da = raster.load_raster(small_dem)
    west, east, south, north = raster.get_extent(da)
    assert west < east
    assert south < north


@pytest.mark.network
def test_to_coordinates_geographic(small_dem):
    da = raster.to_coordinates(raster.load_raster(small_dem), "geographic")
    crs = raster.get_crs(da)
    assert crs.is_geographic
    assert crs.to_epsg() == 4326


@pytest.mark.network
def test_to_coordinates_local_is_projected_utm(small_dem):
    da = raster.to_coordinates(raster.load_raster(small_dem), "local")
    crs = raster.get_crs(da)
    assert crs.is_projected
    # Isle of Eigg (~6.1 W) falls in UTM zone 29N (EPSG:32629).
    assert crs.to_epsg() == 32629


@pytest.mark.network
def test_to_coordinates_sets_nodata_to_nan_on_reproject(small_dem):
    da = raster.to_coordinates(raster.load_raster(small_dem), "local")
    assert np.isnan(da.values).any()


@pytest.mark.network
def test_to_coordinates_explicit_epsg(small_dem):
    da = raster.to_coordinates(raster.load_raster(small_dem), "EPSG:3857")
    assert raster.get_crs(da).to_epsg() == 3857


def test_to_coordinates_missing_crs_raises():
    da = xr.DataArray(np.ones((2, 2)), dims=("y", "x"), coords={"x": [0, 1], "y": [0, 1]})
    with pytest.raises(ValueError, match="no CRS metadata"):
        raster.to_coordinates(da, "local")


def test_resolve_target_crs_local_utm_unavailable_raises(monkeypatch):
    class DummyRio:
        def estimate_utm_crs(self):
            return None

    class DummyDa:
        rio = DummyRio()

    with pytest.raises(ValueError, match="Could not estimate a local/projected CRS"):
        raster.resolve_target_crs(DummyDa(), "local")
