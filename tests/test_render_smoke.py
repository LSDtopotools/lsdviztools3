"""Smoke tests for the render layer."""

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Point

from lsdviztools3.io.raster import load_raster
from lsdviztools3.render import channels, hillshade, points
from lsdviztools3.render.style import MapStyle


def test_mapstyle_defaults():
    """MapStyle should construct with defaults."""
    style = MapStyle()
    assert style.cmap
    assert style.figure_size


@pytest.mark.network
@pytest.mark.parametrize("coordinates", ["local", "geographic"])
def test_render_hillshade_writes_png(small_dem, tmp_path, coordinates):
    """render_hillshade produces a non-empty PNG in both coordinate modes."""
    dem = load_raster(small_dem)
    fig = hillshade.render_hillshade(dem, coordinates=coordinates, style=MapStyle())
    out = tmp_path / f"hillshade_{coordinates}.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_render_channels_standalone_writes_png(tmp_path):
    """render_channels builds its own framed figure from a synthetic network."""
    gdf = gpd.GeoDataFrame(
        {"id": [1, 2]},
        geometry=[
            LineString([(-6.14, 56.88), (-6.12, 56.90)]),
            LineString([(-6.12, 56.90), (-6.10, 56.92)]),
        ],
        crs="EPSG:4326",
    )
    fig = channels.render_channels(gdf, coordinates="geographic", style=MapStyle())
    out = tmp_path / "channels.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_render_channels_empty_raises():
    gdf = gpd.GeoDataFrame({"id": []}, geometry=[], crs="EPSG:4326")
    with pytest.raises(ValueError, match="no features"):
        channels.render_channels(gdf, coordinates="geographic", style=MapStyle())


def test_render_basins_empty_raises():
    from lsdviztools3.render import basins

    gdf = gpd.GeoDataFrame({"id": []}, geometry=[], crs="EPSG:4326")
    with pytest.raises(ValueError, match="no features"):
        basins.render_basins(gdf, coordinates="geographic", style=MapStyle())


def test_log_scaled_sizes_maps_to_point_range():
    """log10(area) maps to the requested pt range (smallest..largest), in cm."""
    from lsdviztools3.render.points import _CM_PER_PT, _log_scaled_sizes

    sizes = _log_scaled_sizes([10.0, 100.0, 1000.0], (1.0, 3.0))
    # Min value -> 1 pt, max value -> 3 pt; monotonic in between.
    assert sizes[0] == pytest.approx(1.0 * _CM_PER_PT)
    assert sizes[-1] == pytest.approx(3.0 * _CM_PER_PT)
    assert sizes[0] < sizes[1] < sizes[2]


def _point_gdf():
    return gpd.GeoDataFrame(
        {"drainage_area": [1e5, 1e6, 1e7], "stream_order": [1, 2, 3]},
        geometry=[Point(-6.14, 56.88), Point(-6.12, 56.90), Point(-6.10, 56.92)],
        crs="EPSG:4326",
    )


def test_render_points_empty_raises():
    gdf = gpd.GeoDataFrame({"id": []}, geometry=[], crs="EPSG:4326")
    with pytest.raises(ValueError, match="no features"):
        points.render_points(gdf, coordinates="geographic", style=MapStyle())


def test_render_points_size_and_colour_writes_png(tmp_path):
    """render_points scales size by log10(area) and colours by a column."""
    fig = points.render_points(
        _point_gdf(), size_by="drainage_area", color_by="stream_order",
        log_color=False, size_range=(1.0, 3.0), coordinates="geographic",
        style=MapStyle(),
    )
    out = tmp_path / "points.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_render_channels_colour_by_writes_png(tmp_path):
    """render_channels colours lines by a column through a CPT."""
    gdf = gpd.GeoDataFrame(
        {"k_sn": [5.0, 50.0]},
        geometry=[
            LineString([(-6.14, 56.88), (-6.12, 56.90)]),
            LineString([(-6.12, 56.90), (-6.10, 56.92)]),
        ],
        crs="EPSG:4326",
    )
    fig = channels.render_channels(
        gdf, color_by="k_sn", colorbar=True, cmap="viridis",
        coordinates="geographic", style=MapStyle(),
    )
    out = tmp_path / "channels_ksn.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


@pytest.mark.network
def test_render_hillshade_adds_crs_label_for_projected(monkeypatch, small_dem):
    calls = {"count": 0}

    def _spy(fig, extent, crs):
        calls["count"] += 1
        assert crs is not None and crs.is_projected
        return fig

    monkeypatch.setattr(hillshade, "add_crs_label", _spy)
    hillshade.render_hillshade(load_raster(small_dem), coordinates="local", style=MapStyle())
    assert calls["count"] == 1


@pytest.mark.network
def test_render_hillshade_transparent_pixels_are_not_red(small_dem, tmp_path):
    try:
        from PIL import Image
    except ImportError:
        pytest.skip("Pillow not available")

    fig = hillshade.render_hillshade(load_raster(small_dem), coordinates="local", style=MapStyle())
    out = tmp_path / "hillshade_transparency.png"
    fig.savefig(out, transparent=True)

    rgba = np.array(Image.open(out).convert("RGBA"))
    transparent = rgba[:, :, 3] == 0
    if not transparent.any():
        pytest.skip("No transparent pixels found in rendered output")

    transparent_rgb = rgba[:, :, :3][transparent]
    red_like = (
        (transparent_rgb[:, 0] > 200)
        & (transparent_rgb[:, 1] < 80)
        & (transparent_rgb[:, 2] < 80)
    )
    assert not red_like.any()
