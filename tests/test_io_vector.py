"""Tests for lsdviztools3.io.vector.

Uses a tiny synthetic GeoDataFrame written to FlatGeobuf — no heavy fixtures.
"""

import geopandas as gpd
import pytest
from shapely.geometry import Point

from lsdviztools3.io import vector


def _tiny_fgb(tmp_path):
    gdf = gpd.GeoDataFrame(
        {"source_id": [1, 2, 3]},
        geometry=[Point(-6.13, 56.90), Point(-6.12, 56.91), Point(-6.11, 56.89)],
        crs="EPSG:4326",
    )
    path = tmp_path / "points.fgb"
    gdf.to_file(path, driver="FlatGeobuf")
    return path


def test_load_vector_roundtrip(tmp_path):
    path = _tiny_fgb(tmp_path)
    gdf = vector.load_vector(path)
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert len(gdf) == 3
    assert gdf.crs is not None


def test_to_xy_not_implemented():
    """to_xy is still a stub."""
    with pytest.raises(NotImplementedError):
        vector.to_xy(None)


def test_to_coordinates_missing_crs_raises():
    gdf = gpd.GeoDataFrame({"id": [1]}, geometry=[Point(0, 0)])
    with pytest.raises(ValueError, match="no CRS metadata"):
        vector.to_coordinates(gdf, "geographic")


def test_local_keeps_a_projected_vector_in_its_own_crs():
    """A British National Grid layer stays in EPSG:27700 for "local"; "utm" moves it."""
    import geopandas as gpd
    from shapely.geometry import LineString
    from lsdviztools3.io import vector

    gdf = gpd.GeoDataFrame(geometry=[LineString([(172000, 633000), (172100, 633100)])],
                           crs="EPSG:27700")
    assert vector.to_coordinates(gdf, "local").crs.to_epsg() == 27700
    assert vector.to_coordinates(gdf, "utm").crs.to_epsg() == 32630
