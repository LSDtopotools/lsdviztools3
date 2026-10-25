"""Vector IO for lsdtt3 outputs.

Reads FlatGeobuf (``.fgb``), GeoPackage (``.gpkg``), and Shapefile (``.shp``)
into ``geopandas.GeoDataFrame``. No PyGMT imports here.
"""

from __future__ import annotations

import geopandas as gpd

_GEOGRAPHIC_ALIASES = {"geographic", "geo", "4326", "epsg:4326", "wgs84"}
_LOCAL_ALIASES = {"local", "utm", "projected", "native"}


def load_vector(path, layer=None) -> gpd.GeoDataFrame:
    """Load a vector dataset as a ``GeoDataFrame``.

    Handles ``.fgb`` / ``.gpkg`` / ``.shp``. ``layer`` selects a layer for
    multi-layer containers such as GeoPackage.

    Args:
        path: Path to the vector file.
        layer: Optional layer name (GeoPackage / multi-layer sources).

    Returns:
        The features as a ``geopandas.GeoDataFrame``.
    """
    if layer is not None:
        return gpd.read_file(path, layer=layer)
    return gpd.read_file(path)


def to_xy(gdf: gpd.GeoDataFrame):
    """Extract x/y coordinates (and attributes) for PyGMT ``plot``.

    For point layers, returns per-feature x/y. For line/polygon layers,
    returns the structure PyGMT expects (e.g. coordinate sequences).

    Args:
        gdf: Input ``GeoDataFrame``.

    Returns:
        x/y data suitable for ``pygmt.Figure.plot`` (e.g. a DataFrame or
        coordinate arrays). Exact shape TBD by geometry type.
    """
    # TODO: branch on geometry type; explode multi-geoms as needed.
    raise NotImplementedError


def to_coordinates(gdf: gpd.GeoDataFrame, coordinates="local") -> gpd.GeoDataFrame:
    """Reproject a ``GeoDataFrame`` to the chosen coordinate system if needed.

    Mirrors ``io.raster.to_coordinates`` so vector overlays land in the same CRS
    as the base raster: ``"local"`` -> native UTM (estimated; a vector in an
    lsdtt3 oblique frame stays in it), ``"geographic"`` -> EPSG:4326, or an
    explicit CRS. Over a frame raster, pass the raster's CRS so vectors in any
    CRS land on it. A no-op when already in the target CRS.
    """
    if gdf.crs is None:
        raise ValueError(
            "Input vector has no CRS metadata; cannot reproject. "
            "Set gdf.crs before calling to_coordinates()."
        )

    if isinstance(coordinates, str):
        key = coordinates.strip().lower()
        if key in _GEOGRAPHIC_ALIASES:
            target = "EPSG:4326"
        elif key in _LOCAL_ALIASES:
            from ..frames import is_oblique_frame

            # Vectors made in an lsdtt3 oblique frame stay in it.
            target = gdf.crs if is_oblique_frame(gdf.crs) else gdf.estimate_utm_crs()
        else:
            target = coordinates
    else:
        target = coordinates

    from pyproj import CRS

    target_crs = CRS.from_user_input(target)
    if gdf.crs is not None and CRS.from_user_input(gdf.crs) == target_crs:
        return gdf
    return gdf.to_crs(target_crs)
