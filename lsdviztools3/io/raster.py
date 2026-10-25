"""Raster IO for lsdtt3 grids.

Reads COG/GeoTIFF rasters (and legacy ENVI ``.bil``, read-only) into
``xarray.DataArray`` via rioxarray, and exposes helpers to hand grids to PyGMT.

No PyGMT imports here — this is part of the backend-agnostic core.
"""

from __future__ import annotations

import atexit

import rioxarray  # noqa: F401  (registers the .rio accessor on xarray objects)
import xarray as xr
from rasterio.crs import CRS

_EXIT_HANDLER_REGISTERED = False


def _close_cached_files() -> None:
    """Close every file xarray still holds open in its global file cache."""
    try:
        from xarray.backends.file_manager import FILE_CACHE
    except ImportError:  # pragma: no cover - private API moved in a future xarray
        return
    try:
        FILE_CACHE.clear()  # evicting a file closes it
    except Exception:  # pragma: no cover - never let cleanup break an exit
        pass


def register_exit_cleanup() -> None:
    """Close rasters opened through rioxarray before the interpreter shuts down.

    ``rioxarray.open_rasterio`` is lazy: the rasterio dataset stays open in
    xarray's global file cache until something closes it. Left open, it is
    closed during interpreter finalisation, after Python has torn down what it
    needs to report errors, and every script that loaded a raster ends by
    printing ``Error in sys.excepthook: / Original exception was:`` with no
    traceback. Closing the cache from an ``atexit`` handler runs while the
    interpreter is still intact. Seen with rioxarray 0.22, xarray 2026.4,
    rasterio 1.5, Python 3.12. Safe to call more than once.
    """
    global _EXIT_HANDLER_REGISTERED
    if not _EXIT_HANDLER_REGISTERED:
        atexit.register(_close_cached_files)
        _EXIT_HANDLER_REGISTERED = True


register_exit_cleanup()

#: Accepted aliases for the two coordinate modes.
_GEOGRAPHIC_ALIASES = {"geographic", "geo", "4326", "epsg:4326", "wgs84"}
_LOCAL_ALIASES = {"local", "utm", "projected", "native"}


def load_raster(path) -> xr.DataArray:
    """Load a raster as an ``xarray.DataArray`` via rioxarray.

    Supports COG/GeoTIFF and legacy ENVI ``.bil`` (read-only). The returned
    array carries CRS and affine transform metadata (``.rio`` accessor).

    Args:
        path: Path to the raster file.

    Returns:
        The raster as a 2-D ``xarray.DataArray`` (a single band is squeezed out).
    """
    da = rioxarray.open_rasterio(path, masked=True)
    if "band" in da.dims and da.sizes["band"] == 1:
        da = da.squeeze("band", drop=True)
    return da


def to_gmt_grid(da: xr.DataArray) -> xr.DataArray:
    """Return a grid object suitable for ``pygmt.Figure.grdimage``.

    PyGMT accepts an ``xarray.DataArray`` directly; this hook centralises the
    normalisation it needs: 2-D with ``y``/``x`` dims.

    Args:
        da: Raster as an ``xarray.DataArray``.

    Returns:
        A grdimage-ready ``DataArray``.
    """
    if "band" in da.dims:
        da = da.squeeze("band", drop=True)
    if da.ndim != 2:
        raise ValueError(f"Expected a 2-D grid for PyGMT, got dims {da.dims}.")
    return da


def get_extent(da: xr.DataArray) -> list:
    """Return the bounding extent as ``[west, east, south, north]`` (GMT order).

    Args:
        da: Raster as an ``xarray.DataArray``.

    Returns:
        Extent ``[w, e, s, n]`` in the raster's CRS units.
    """
    west, south, east, north = da.rio.bounds()
    return [west, east, south, north]


def get_crs(da: xr.DataArray):
    """Return the raster CRS (``da.rio.crs``)."""
    return da.rio.crs


def reproject_raster(da: xr.DataArray, crs) -> xr.DataArray:
    """Reproject a raster to a target CRS (``da.rio.reproject``)."""
    return da.rio.reproject(crs)


def resolve_target_crs(da: xr.DataArray, coordinates="local") -> CRS:
    """Resolve a ``coordinates`` choice to a concrete CRS.

    Args:
        da: Raster (used to estimate the local UTM zone when needed).
        coordinates: ``"local"``/``"utm"`` for the native projected CRS (UTM
            estimated from the raster's extent; a raster in an lsdtt3 oblique
            frame keeps its frame), ``"geographic"`` for EPSG:4326,
            or any explicit CRS (EPSG code, proj string, or CRS object).

    Returns:
        A ``rasterio.crs.CRS``.
    """
    if isinstance(coordinates, str):
        key = coordinates.strip().lower()
        if key in _GEOGRAPHIC_ALIASES:
            return CRS.from_epsg(4326)
        if key in _LOCAL_ALIASES:
            # An lsdtt3 oblique frame is the local projected CRS of its own
            # raster: keep it rather than warping back to UTM.
            from ..frames import is_oblique_frame

            native = getattr(da.rio, "crs", None)
            if native is not None and is_oblique_frame(native):
                return CRS.from_user_input(native)
            utm = da.rio.estimate_utm_crs()
            if utm is None:
                raise ValueError(
                    "Could not estimate a local/projected CRS from raster extent. "
                    "Provide an explicit CRS (for example coordinates='EPSG:32629')."
                )
            return CRS.from_user_input(utm)
        return CRS.from_user_input(coordinates)
    return CRS.from_user_input(coordinates)


def to_coordinates(da: xr.DataArray, coordinates="local") -> xr.DataArray:
    """Reproject a raster to the chosen coordinate system if needed.

    ``"local"`` -> native UTM (estimated), ``"geographic"`` -> EPSG:4326, or an
    explicit CRS. A no-op when the raster is already in the target CRS.

    Args:
        da: Raster as an ``xarray.DataArray``.
        coordinates: see :func:`resolve_target_crs`.

    Returns:
        The (possibly reprojected) ``xarray.DataArray``.
    """
    if da.rio.crs is None:
        raise ValueError(
            "Input raster has no CRS metadata; cannot reproject. "
            "Assign a CRS first or pass data that includes georeferencing."
        )
    target = resolve_target_crs(da, coordinates)
    if da.rio.crs is not None and CRS.from_user_input(da.rio.crs) == target:
        return da
    return da.rio.reproject(target, nodata=float("nan"))
