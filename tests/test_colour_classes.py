"""Colour classes: breaks, schemes, limits, colours, labels, legends, and maps."""

from __future__ import annotations

import os

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Point

from lsdviztools3 import colour_classes as cc

# --------------------------------------------------------------------------
# Breaks and classify


def test_breaks_given_make_one_more_class_and_edges_go_up():
    c = cc.make_colour_classes(None, [50, 100, 200], cmap="viridis", vmin=0, vmax=400)
    assert c.n_classes == 4
    # A value equal to a break is in the class above it; nodata has no class.
    got = c.classify([0, 49.9, 50, 99.9, 100, 200, 1e9, np.nan, np.inf])
    np.testing.assert_array_equal(got[:7], [0, 0, 1, 1, 2, 3, 3])
    assert np.isnan(got[7]) and np.isnan(got[8])
    assert c.labels == ("below 50", "50 to 100", "100 to 200", "200 and above")


@pytest.mark.parametrize("breaks", [[1, 1], [3, 2], [1, np.nan], []])
def test_bad_breaks_refused(breaks):
    with pytest.raises(ValueError):
        cc.make_colour_classes([1, 2, 3], breaks)


def test_breaks_accept_list_tuple_and_array():
    a = cc.make_colour_classes([0, 10], [2, 5], cmap="viridis")
    b = cc.make_colour_classes([0, 10], (2, 5), cmap="viridis")
    d = cc.make_colour_classes([0, 10], np.array([2.0, 5.0]), cmap="viridis")
    assert a.breaks == b.breaks == d.breaks == (2.0, 5.0)


def test_natural_breaks_finds_the_clusters():
    assert cc.natural_breaks([1, 2, 3, 10, 11, 12, 20, 21, 22], 3).tolist() == [10, 20]
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(10, 1, 300), rng.normal(50, 1, 300), rng.normal(90, 1, 300)])
    b = cc.class_breaks(x, "jenks", 3)
    # Each break is the lowest value of the cluster above it.
    assert b[0] == pytest.approx(x[(x > 30) & (x < 70)].min())
    assert b[1] == pytest.approx(x[x > 70].min())
    c = cc.make_colour_classes(x, "natural_breaks", n_classes=3, cmap="viridis")
    counts = np.bincount(c.classify(x).astype(int))
    assert counts.tolist() == [300, 300, 300]


def test_natural_breaks_matches_brute_force_on_small_data():
    rng = np.random.default_rng(3)
    x = np.sort(rng.gamma(2.0, 5.0, 14))

    def cost(groups):
        return sum(((g - g.mean()) ** 2).sum() for g in groups)

    best = None
    for i in range(1, len(x) - 1):
        for j in range(i + 1, len(x)):
            c = cost([x[:i], x[i:j], x[j:]])
            if best is None or c < best[0] - 1e-12:
                best = (c, [x[i], x[j]])
    assert cc.natural_breaks(x, 3).tolist() == pytest.approx(best[1])


def test_natural_breaks_samples_large_data():
    x = np.random.default_rng(1).uniform(0, 1, 50_000)
    b = cc.natural_breaks(x, 4)
    assert len(b) == 3 and np.all(np.diff(b) > 0)


def test_equal_interval_and_quantile():
    x = np.arange(101.0)
    assert cc.class_breaks(x, "equal_interval", 4).tolist() == [25, 50, 75]
    assert cc.class_breaks(x, "quantile", 4).tolist() == [25, 50, 75]
    # Limits bound the scheme.
    assert cc.class_breaks(x, "equal_interval", 2, vmin=20, vmax=40).tolist() == [30]
    assert cc.class_breaks(x, "quantile", 2, vmin=0, vmax=50).tolist() == [25]


def test_log_scale_scheme_gives_equal_ratios():
    b = cc.class_breaks([1, 10, 100, 1000, 10000], "equal_interval", 4, log_scale=True)
    assert b.tolist() == pytest.approx([10, 100, 1000])


def test_too_few_distinct_values_merge_classes():
    b = cc.class_breaks([1, 1, 1, 2, 2], "quantile", 5)
    assert b.tolist() == [2.0]
    with pytest.raises(ValueError, match="single value"):
        cc.make_colour_classes([3, 3, 3], "jenks", n_classes=3)


@pytest.mark.parametrize("kwargs, match", [
    (dict(classes="nope"), "Unknown classification scheme"),
    (dict(classes="quantile", n_classes=1), "at least 2"),
    (dict(classes="quantile", vmin=5, vmax=5), "below vmax"),
    (dict(classes="quantile", log_scale=True, vmin=-1), "> 0"),
    (dict(classes=[1, 2], color_mode="rainbow"), "color_mode"),
    (dict(classes=[1, 2], class_colors=["red", "blue"]), "3 classes"),
])
def test_bad_arguments_refused(kwargs, match):
    with pytest.raises(ValueError, match=match):
        cc.make_colour_classes([1.0, 2.0, 3.0, 4.0], **kwargs)


def test_scheme_needs_values():
    with pytest.raises(ValueError, match="needs the data"):
        cc.make_colour_classes(None, "quantile")


# --------------------------------------------------------------------------
# Colours


def test_gradient_colours_sit_at_class_centres_between_limits():
    from matplotlib import colormaps

    c = cc.make_colour_classes(None, [20, 60], cmap="viridis", vmin=0, vmax=100)
    # Centres 10, 40, 80 of 0..100.
    expected = [cc._rgb_text(colormaps["viridis"](t)) for t in (0.1, 0.4, 0.8)]
    assert list(c.colors) == expected
    assert c.limits == (0.0, 100.0)


def test_uneven_breaks_give_uneven_colours():
    pos = cc.class_positions([10, 20, 90], 0, 100)
    assert pos.tolist() == pytest.approx([0.05, 0.15, 0.55, 0.95])
    # Log scale: centres in log10.
    pos = cc.class_positions([100], 10, 1000, log_scale=True)
    assert pos.tolist() == pytest.approx([0.25, 0.75])


def test_limits_default_to_the_data():
    c = cc.make_colour_classes([0.0, 100.0], [50], cmap="viridis")
    assert c.limits == (0.0, 100.0)


def test_distinct_colours_in_order():
    from matplotlib import colormaps

    c = cc.make_colour_classes(None, [1, 2, 3], cmap="tab10")      # auto -> distinct
    assert list(c.colors) == [cc._rgb_text(colormaps["tab10"](i)) for i in range(4)]
    assert c.limits is None
    with pytest.raises(ValueError, match="10 colours"):
        cc.make_colour_classes(None, list(range(1, 12)), cmap="tab10")
    with pytest.raises(ValueError, match="qualitative"):
        cc.make_colour_classes(None, [1, 2], cmap="viridis", color_mode="distinct")


def test_class_colors_given_and_grey_nudged():
    c = cc.make_colour_classes(None, [1], class_colors=["#808080", "10/20/30"])
    assert c.colors == ("128/128/129", "10/20/30")


def test_gmt_colour_map_is_sampled():
    pytest.importorskip("pygmt")
    c = cc.make_colour_classes(None, [20, 60], cmap="geo", vmin=0, vmax=100)
    assert len(c.colors) == 3
    for col in c.colors:
        assert all(0 <= int(v) <= 255 for v in col.split("/"))


def test_labels_keep_breaks_distinct_and_take_a_format():
    assert cc.class_labels([1.0001, 1.0002]) == ["below 1.0001", "1.0001 to 1.0002", "1.0002 and above"]
    # Written out in full, not 1.23e+03.
    assert cc.class_labels([1234.5678, 20000.25]) == ["below 1230", "1230 to 20000", "20000 and above"]
    assert cc.class_labels([10, 20], ".1f")[1] == "10.0 to 20.0"


def test_parse_classes_arg():
    assert cc.parse_classes_arg("Jenks") == "jenks"
    assert cc.parse_classes_arg("50, 100,200") == [50.0, 100.0, 200.0]
    assert cc.parse_classes_arg("") is None
    with pytest.raises(ValueError, match="neither"):
        cc.parse_classes_arg("50,abc")


def test_legend_width_holds_the_longest_text():
    from lsdviztools3.render.categorical import class_legend_width

    short = class_legend_width(["below 5"])
    long = class_legend_width(["123456.789 and above"])
    assert long > short
    # A long title widens the legend; GMT escapes take no width.
    assert class_legend_width(["a"], header="A rather long legend title") > class_legend_width(["a"])
    assert class_legend_width(["a"], header="k@-sn@-") == class_legend_width(["a"], header="ksn")


# --------------------------------------------------------------------------
# Maps (need GMT)

pygmt = pytest.importorskip("pygmt")
from lsdviztools3.render import channels, points  # noqa: E402
from lsdviztools3.render.style import MapStyle  # noqa: E402


def _raster(values):
    import rioxarray  # noqa: F401
    import xarray as xr

    values = np.asarray(values, dtype=float)
    ny, nx = values.shape
    da = xr.DataArray(values, dims=("y", "x"),
                      coords={"y": 6300000.0 + np.arange(ny)[::-1] * 30.0,
                              "x": 400000.0 + np.arange(nx) * 30.0})
    return da.rio.write_crs("EPSG:32630")


def _png(fig, path):
    fig.savefig(path)
    assert path.exists() and path.stat().st_size > 0
    return path


@pytest.mark.parametrize("classes", [[10, 30, 50], "equal_interval", "quantile", "natural_breaks"])
def test_render_drape_in_classes(tmp_path, classes):
    base = _raster(np.ones((16, 16)))
    drape = _raster(np.arange(256.0).reshape(16, 16) / 4)
    fig = channels.render_drape(base, drape, coordinates="local", classes=classes,
                                n_classes=4, cmap="viridis", legend_title="Gradient",
                                style=MapStyle(transparent=False))
    _png(fig, tmp_path / "drape_classes.png")


def test_render_drape_classes_and_gt_breaks_exclusive():
    base = _raster(np.ones((4, 4)))
    with pytest.raises(ValueError, match="not both"):
        channels.render_drape(base, base, coordinates="local", classes=[1],
                              categorical_gt_breaks=[0.5], categorical_colors=["red"])


def test_render_drape_classes_log_scale(tmp_path):
    base = _raster(np.ones((8, 8)))
    drape = _raster(np.logspace(-1, 3, 64).reshape(8, 8))
    fig = channels.render_drape(base, drape, coordinates="local", log_color=True,
                                classes="equal_interval", n_classes=4, vmin=1, vmax=1000,
                                cmap="magma", style=MapStyle(transparent=False))
    _png(fig, tmp_path / "drape_log_classes.png")


def _lines():
    return gpd.GeoDataFrame(
        {"k_sn": [5.0, 30.0, 80.0, 300.0]},
        geometry=[LineString([(-6.14, 56.88), (-6.12, 56.90)]),
                  LineString([(-6.12, 56.90), (-6.10, 56.92)]),
                  LineString([(-6.10, 56.92), (-6.08, 56.93)]),
                  LineString([(-6.08, 56.93), (-6.06, 56.95)])],
        crs="EPSG:4326")


@pytest.mark.parametrize("classes,cmap", [([10, 50, 100], "viridis"), ("jenks", "tab10"),
                                          ("quantile", "geo")])
def test_render_channels_in_classes(tmp_path, classes, cmap):
    fig = channels.render_channels(_lines(), color_by="k_sn", classes=classes, n_classes=3,
                                   cmap=cmap, coordinates="geographic", legend_title="k@-sn@-",
                                   style=MapStyle())
    _png(fig, tmp_path / "channels_classes.png")


def test_render_points_in_classes(tmp_path):
    gdf = gpd.GeoDataFrame(
        {"area": [1e5, 1e6, 1e7, 1e8]},
        geometry=[Point(-6.14, 56.88), Point(-6.12, 56.90), Point(-6.10, 56.92), Point(-6.08, 56.93)],
        crs="EPSG:4326")
    fig = points.render_points(gdf, color_by="area", size_by="area", classes="equal_interval",
                               n_classes=3, cmap="viridis", coordinates="geographic",
                               style=MapStyle())
    _png(fig, tmp_path / "points_classes.png")
    with pytest.raises(ValueError, match="color_by"):
        points.render_points(gdf, classes=[1], coordinates="geographic", style=MapStyle())


def test_render_points_continuous_uses_vmin_vmax(tmp_path):
    gdf = gpd.GeoDataFrame({"v": [1.0, 10.0, 100.0]},
                           geometry=[Point(0, 0), Point(0.01, 0.01), Point(0.02, 0.02)],
                           crs="EPSG:4326")
    fig = points.render_points(gdf, color_by="v", vmin=2, vmax=50, colorbar=True,
                               coordinates="geographic", style=MapStyle())
    _png(fig, tmp_path / "points_limits.png")


# --------------------------------------------------------------------------
# Command lines


def test_cli_drape_and_channels_take_classes(tmp_path):
    from lsdviztools3.cli import plot_channels, plot_drape

    base = tmp_path / "base.tif"
    drape = tmp_path / "drape.tif"
    _raster(np.ones((16, 16))).rio.to_raster(base)
    _raster(np.arange(256.0).reshape(16, 16)).rio.to_raster(drape)
    out = tmp_path / "drape.png"
    assert plot_drape.main(["--base", str(base), "--drape", str(drape), "--out", str(out),
                            "--classes", "50,100,200", "--cmap", "viridis"]) == 0
    assert out.exists()
    lines = tmp_path / "lines.fgb"
    gdf = _lines().to_crs("EPSG:32630")
    gdf.to_file(lines)
    out2 = tmp_path / "channels.png"
    assert plot_channels.main(["--base", str(base), "--lines", str(lines), "--out", str(out2),
                               "--color-by", "k_sn", "--classes", "jenks", "--n-classes", "3",
                               "--cmap", "viridis"]) == 0
    with pytest.raises(SystemExit, match="color-by"):
        plot_channels.main(["--base", str(base), "--lines", str(lines), "--out", str(out2),
                            "--classes", "jenks"])
    with pytest.raises(SystemExit, match="neither"):
        plot_drape.main(["--base", str(base), "--drape", str(drape), "--out", str(out),
                         "--classes", "1,x"])


def test_cli_chi_and_knickpoints_have_ksn_classes():
    from lsdviztools3.cli import plot_chi, plot_knickpoints

    for mod in (plot_chi, plot_knickpoints):
        a = mod.build_parser().parse_args(["--data-dir", ".", "--prefix", "p", "--m-over-n", "0.45",
                                           "--out-dir", ".", "--ksn-classes", "20,50,100",
                                           "--class-colors", "a,b,c,d"])
        kw = mod._ksn_class_kwargs(a)
        assert kw["ksn_classes"] == [20.0, 50.0, 100.0]
        assert kw["class_colors"] == ["a", "b", "c", "d"]


def test_basins_map_has_no_cosmogenic_text_unless_asked(tmp_path, monkeypatch):
    """lsdviz3-basins wrote 'scaling=stone, muon=..., atmosphere=ERA40' on every map."""
    from lsdviztools3.cli import plot_basins

    calls = []
    monkeypatch.setattr(pygmt.Figure, "text", lambda self, **kw: calls.append(kw))
    base = tmp_path / "base.tif"
    _raster(np.ones((16, 16))).rio.to_raster(base)
    basins = tmp_path / "basins.fgb"
    from shapely.geometry import box
    gpd.GeoDataFrame({"basin_key": [0]}, geometry=[box(400060, 6300060, 400300, 6300300)],
                     crs="EPSG:32630").to_file(basins)
    out = tmp_path / "basins.png"
    assert plot_basins.main(["--base", str(base), "--basins", str(basins), "--out", str(out)]) == 0
    def cosmo_text():
        return [c["text"] for c in calls if "scaling" in str(c.get("text", ""))]

    assert cosmo_text() == []          # (the CRS label is drawn with text too)
    assert plot_basins.main(["--base", str(base), "--basins", str(basins), "--out", str(out),
                             "--scaling-scheme", "stone"]) == 0
    assert cosmo_text() == ["scaling=stone"]
