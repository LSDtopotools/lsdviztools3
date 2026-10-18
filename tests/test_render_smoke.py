"""Smoke tests for the render layer."""

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Point

from lsdviztools3.io.raster import load_raster
from lsdviztools3.render import channels, hillshade, points, sti
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


def _fixed_size_points_are_visible(fig, tmp_path, name, rgb_match):
    """Return True if the rendered figure contains pixels close to rgb_match."""
    try:
        from PIL import Image
    except ImportError:
        pytest.skip("Pillow not available")
    out = tmp_path / name
    fig.savefig(out)
    arr = np.array(Image.open(out).convert("RGB")).reshape(-1, 3)
    r, g, b = rgb_match
    near = (
        (np.abs(arr[:, 0].astype(int) - r) < 60)
        & (np.abs(arr[:, 1].astype(int) - g) < 60)
        & (np.abs(arr[:, 2].astype(int) - b) < 60)
    )
    return bool(near.any())


def test_render_points_fixed_size_filled_visible(tmp_path):
    """A fixed-size filled symbol must actually render (regression: a scalar
    size= passed to fig.plot drew nothing)."""
    fig = points.render_points(
        _point_gdf(), coordinates="geographic",
        style=MapStyle(point_size="0.5c", line_color="yellow"),
    )
    assert _fixed_size_points_are_visible(
        fig, tmp_path, "pts_fixed_c.png", (255, 255, 0))


def test_render_points_outline_only_symbol_visible(tmp_path):
    """An outline-only symbol (x) renders with its pen, not a fill."""
    fig = points.render_points(
        _point_gdf(), coordinates="geographic", symbol="x",
        pen="3p,red", style=MapStyle(point_size="0.6c"),
    )
    assert _fixed_size_points_are_visible(
        fig, tmp_path, "pts_x.png", (255, 0, 0))


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


def _synthetic_raster(values, crs="EPSG:32630"):
    """Small projected DataArray with a CRS, for offline render tests."""
    import xarray as xr
    import rioxarray  # noqa: F401 - registers the .rio accessor

    values = np.asarray(values, dtype=float)
    ny, nx = values.shape
    x = 400000.0 + np.arange(nx) * 30.0
    y = 6300000.0 + np.arange(ny)[::-1] * 30.0
    da = xr.DataArray(values, dims=("y", "x"), coords={"y": y, "x": x})
    return da.rio.write_crs(crs)


def test_render_drape_writes_png(tmp_path):
    """render_drape layers a coloured drape over a gray hillshade base."""
    base = _synthetic_raster(np.outer(np.linspace(0, 1, 8), np.linspace(0, 1, 8)))
    drape = _synthetic_raster(np.outer(np.arange(1, 9), np.arange(1, 9)))
    fig = channels.render_drape(
        base, drape, coordinates="local", transparency=40,
        cmap="roma", colorbar=True, colorbar_label="Curvature",
        style=MapStyle(transparent=False),
    )
    out = tmp_path / "drape.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_render_drape_log_clamps_nonpositive(tmp_path):
    """log_color renders (zeros/negatives clamped to transparent, not crashing)."""
    base = _synthetic_raster(np.ones((8, 8)))
    # Mix of <=0 (clamped) and positive values.
    drape = _synthetic_raster(np.linspace(-5, 100, 64).reshape(8, 8))
    fig = channels.render_drape(
        base, drape, coordinates="local", log_color=True,
        style=MapStyle(transparent=False),
    )
    out = tmp_path / "drape_log.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_render_drape_log_all_nonpositive_raises():
    """All values <= 0 under log scale leaves nothing finite to render."""
    base = _synthetic_raster(np.ones((4, 4)))
    drape = _synthetic_raster(np.full((4, 4), -1.0))
    with pytest.raises(ValueError, match="clamping values <= 0"):
        channels.render_drape(base, drape, coordinates="local", log_color=True,
                              style=MapStyle())


def test_render_drape_bad_out_of_bounds_raises():
    base = _synthetic_raster(np.ones((4, 4)))
    drape = _synthetic_raster(np.arange(16.0).reshape(4, 4))
    with pytest.raises(ValueError, match="clamp.*clip"):
        channels.render_drape(base, drape, coordinates="local",
                              out_of_bounds="nope", style=MapStyle())


def test_render_drape_percentile_clip_threshold_mask(tmp_path):
    """percentile_min + clip + solid colour = mask above the Nth percentile."""
    base = _synthetic_raster(np.ones((16, 16)))
    drape = _synthetic_raster(np.arange(256.0).reshape(16, 16))
    fig = channels.render_drape(
        base, drape, coordinates="local",
        cmap="black", percentile_min=90, out_of_bounds="clip",
        transparency=0, style=MapStyle(transparent=False),
    )
    out = tmp_path / "mask.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_nice_bar_km_snaps_to_round_numbers():
    """The clean scale bar picks the largest preferred length <= width/10."""
    from lsdviztools3.render.decorations import _nice_bar_km

    assert _nice_bar_km(8.5) == 8.0       # ~85 km map
    assert _nice_bar_km(5.0) == 5.0
    assert _nice_bar_km(1.2) == 1.0
    assert _nice_bar_km(0.07) == 0.1      # below smallest -> smallest
    assert _nice_bar_km(90.0) == 75.0


def test_render_drape_clean_writes_png(tmp_path):
    """clean drape: no colourbar/crs furniture, just map + solid scale bar."""
    base = _synthetic_raster(np.outer(np.linspace(0, 1, 16), np.linspace(0, 1, 16)))
    drape = _synthetic_raster(np.arange(256.0).reshape(16, 16))
    fig = channels.render_drape(
        base, drape, coordinates="local", cmap="Spectral",
        colorbar=True,  # should be suppressed by clean
        style=MapStyle(clean=True, transparent=False),
    )
    out = tmp_path / "clean_drape.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_render_hillshade_clean_writes_png(tmp_path):
    base = _synthetic_raster(np.outer(np.linspace(0, 1, 16), np.linspace(0, 1, 16)))
    fig = hillshade.render_hillshade(
        base, coordinates="local", shade=False,
        style=MapStyle(cmap="gray", clean=True, transparent=False),
    )
    out = tmp_path / "clean_hs.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_render_sti_p95_drape_writes_png(tmp_path):
    """STI helper writes a non-empty STI P95 drape PNG."""
    base = _synthetic_raster(np.outer(np.linspace(0, 1, 12), np.linspace(0, 1, 12)))
    sti_vals = np.linspace(0.0, 1.0, 144).reshape(12, 12)
    sti_vals[sti_vals < 0.9] = np.nan
    drape = _synthetic_raster(sti_vals)
    fig = sti.render_sti_p95_drape(
        base,
        drape,
        coordinates="local",
        style=MapStyle(transparent=False),
    )
    out = tmp_path / "sti_p95_drape.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_resolve_cmap_matplotlib_and_solid(tmp_path):
    """Matplotlib-only names -> CPT file; GMT names pass through; colour -> CPT."""
    from lsdviztools3.render.colormaps import resolve_cmap

    # GMT-native names and colour-lists pass straight through.
    assert resolve_cmap("geo") == "geo"
    assert resolve_cmap("roma") == "roma"
    assert resolve_cmap("red,blue") == "red,blue"
    # Matplotlib-only map -> a generated .cpt file.
    spectral = resolve_cmap("Spectral")
    assert spectral.endswith(".cpt") and "Spectral" in spectral
    # A single named colour -> a solid CPT file (for threshold masks).
    black = resolve_cmap("black")
    assert black.endswith(".cpt")


def test_resolve_cmap_gray_nudge_can_be_disabled():
    """Gray CPT nudge is optional when transparency workaround is unnecessary."""
    from lsdviztools3.render.colormaps import resolve_cmap

    nudged = resolve_cmap("gray", avoid_pure_gray_bug=True)
    neutral = resolve_cmap("gray", avoid_pure_gray_bug=False)

    assert nudged.endswith(".cpt")
    assert neutral.endswith(".cpt")
    assert nudged != neutral


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
@pytest.mark.parametrize("cmap", ["geo", "gray"])
def test_render_hillshade_transparent_pixels_are_not_red(small_dem, tmp_path, cmap):
    try:
        from PIL import Image
    except ImportError:
        pytest.skip("Pillow not available")

    # cmap="gray" is the regression case: GMT's native grey CPT with
    # nan_transparent paints NaN cells red unless the CPT is made non-pure-grey.
    fig = hillshade.render_hillshade(
        load_raster(small_dem), coordinates="local", style=MapStyle(cmap=cmap))
    out = tmp_path / f"hillshade_transparency_{cmap}.png"
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


class _RecordingFig:
    """Stand-in figure recording plot/text calls (no GMT needed)."""

    def __init__(self):
        self.plots, self.texts = [], []

    def plot(self, **kwargs):
        self.plots.append(kwargs)

    def text(self, **kwargs):
        self.texts.append(kwargs)



def _utm():
    from rasterio.crs import CRS

    return CRS.from_epsg(32617)


def test_bar_scalebar_label_above_by_default():
    """Default layout is unchanged: label centred above the bar."""
    from lsdviztools3.render.decorations import add_bar_scalebar

    fig = _RecordingFig()
    add_bar_scalebar(fig, [600.0, 626.0, 4300.0, 4322.0], _utm())
    bar_top = max(fig.plots[0]["y"])
    assert fig.texts[0]["y"] > bar_top
    assert fig.texts[0]["justify"] == "CB"
    assert fig.texts[0]["text"] == "2 km"
    # Unchanged position: bottom-left inset of 4 % of the map.
    assert min(fig.plots[0]["y"]) == pytest.approx(4300.0 + 0.04 * 22.0)


def test_bar_scalebar_label_below():
    from lsdviztools3.render.decorations import add_bar_scalebar

    fig = _RecordingFig()
    add_bar_scalebar(fig, [600.0, 626.0, 4300.0, 4322.0], _utm(), label_position="below")
    bar_bottom = min(fig.plots[0]["y"])
    assert fig.texts[0]["y"] < bar_bottom
    assert fig.texts[0]["justify"] == "CT"
    assert fig.texts[0]["y"] > 4300.0  # label stays inside the map


def test_bar_scalebar_bad_label_position():
    from lsdviztools3.render.decorations import add_bar_scalebar

    with pytest.raises(ValueError, match="label_position"):
        add_bar_scalebar(_RecordingFig(), [0, 10, 0, 10], _utm(), label_position="left")


@pytest.mark.parametrize("anchor", ["TR", "TL", "BR", "BL"])
def test_north_arrow_points_up(anchor):
    """The arrowhead apex is directly above its base midpoint, N beneath it."""
    from lsdviztools3.render.decorations import add_north_arrow

    fig = _RecordingFig()
    extent = [600.0, 626.0, 4300.0, 4322.0]
    add_north_arrow(fig, extent, _utm(), anchor=anchor)
    xs, ys = fig.plots[0]["x"], fig.plots[0]["y"]
    apex = int(np.argmax(ys))
    base = [i for i in range(3) if i != apex]
    assert ys[base[0]] == ys[base[1]] < ys[apex]
    assert xs[apex] == pytest.approx(0.5 * (xs[base[0]] + xs[base[1]]))
    assert fig.plots[0]["fill"] == "black"
    assert fig.texts[0]["text"] == "N"
    assert fig.texts[0]["y"] < ys[base[0]]
    assert all(extent[0] < x < extent[1] for x in xs)
    assert all(extent[2] < y < extent[3] for y in ys + [fig.texts[0]["y"]])


def test_north_arrow_geographic_ok():
    from rasterio.crs import CRS

    from lsdviztools3.render.decorations import add_north_arrow

    fig = _RecordingFig()
    add_north_arrow(fig, [-79.4, -79.1, 38.8, 39.0], CRS.from_epsg(4326))
    assert fig.texts[0]["text"] == "N"


@pytest.mark.parametrize("projection", ["G-79/39/12c", "X12c/-10c", "S0/90/12c"])
def test_north_arrow_refuses_non_north_up_projection(projection):
    from lsdviztools3.render.decorations import add_north_arrow

    with pytest.raises(ValueError, match="north"):
        add_north_arrow(_RecordingFig(), [0, 10, 0, 10], _utm(), projection=projection)


def test_north_arrow_refuses_south_oriented_crs():
    from rasterio.crs import CRS

    from lsdviztools3.render.decorations import add_north_arrow

    # EPSG:2053 (Hartebeesthoek94 / Lo29) has westing/southing axes.
    with pytest.raises(ValueError, match="not east/north"):
        add_north_arrow(_RecordingFig(), [0, 10, 0, 10], CRS.from_epsg(2053))


def test_map_credit_corner_and_escape():
    from lsdviztools3.render.decorations import add_map_credit

    fig = _RecordingFig()
    add_map_credit(fig, [0.0, 10.0, 0.0, 10.0], "Data @ source")
    t = fig.texts[0]
    assert t["text"] == "Data @@ source" and t["justify"] == "RB"
    assert t["x"] < 10.0 and t["y"] > 0.0
    fig2 = _RecordingFig()
    add_map_credit(fig2, [0, 10, 0, 10], "")
    assert fig2.texts == []


def test_scalebar_and_arrow_label_halo():
    """Optional halo: GMT '=~' outline behind the glyphs; default font unchanged."""
    from lsdviztools3.render.decorations import add_bar_scalebar, add_north_arrow

    plain, halo = _RecordingFig(), _RecordingFig()
    add_bar_scalebar(plain, [0.0, 26.0, 0.0, 22.0], _utm())
    add_bar_scalebar(halo, [0.0, 26.0, 0.0, 22.0], _utm(), label_halo="white")
    assert plain.texts[0]["font"] == "9p,Helvetica,black"
    assert halo.texts[0]["font"] == "9p,Helvetica,black=~2p,white"
    arrow = _RecordingFig()
    add_north_arrow(arrow, [0.0, 26.0, 0.0, 22.0], _utm(), label_halo="white")
    assert arrow.texts[0]["font"].endswith("=~2p,white")
