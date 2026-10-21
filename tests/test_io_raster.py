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


def _write_small_tif(path):
    import rasterio
    from rasterio.transform import from_origin

    with rasterio.open(path, "w", driver="GTiff", height=8, width=8, count=1,
                       dtype="float32", crs="EPSG:32630",
                       transform=from_origin(500000, 6000080, 10, 10)) as dst:
        dst.write(np.arange(64, dtype="float32").reshape(8, 8), 1)


def _run_python(code):
    import subprocess
    import sys

    return subprocess.run([sys.executable, "-c", code], capture_output=True,
                          text=True, timeout=120)


def test_lazy_raster_left_open_exits_cleanly(tmp_path):
    """A script that loads a raster and exits must not print a shutdown error.

    rioxarray keeps the file open in xarray's global cache; closed during
    interpreter finalisation it printed "Error in sys.excepthook" at the end of
    every such script. Only a fresh interpreter's exit shows this, hence the
    subprocess.
    """
    tif = tmp_path / "small.tif"
    _write_small_tif(tif)
    proc = _run_python(
        "from lsdviztools3.io.raster import load_raster\n"
        f"da = load_raster({str(tif)!r})\n"
        "print(float(da.max()))\n")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "63.0"
    assert "Error in sys.excepthook" not in proc.stderr
    assert "Original exception was" not in proc.stderr


def test_exit_cleanup_is_registered_once():
    raster.register_exit_cleanup()
    raster.register_exit_cleanup()
    assert raster._EXIT_HANDLER_REGISTERED
