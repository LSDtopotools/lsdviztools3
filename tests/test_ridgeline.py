"""Ridgeline profile plots (io, charts, map and CLIs) on synthetic outputs.

The synthetic ridgeline runs 3 km due east in UTM 30N at 2 m node spacing,
with elevation ``500 + 0.01 d + 5 sin(2 pi (d + 1) / 200)``. Its curvature
columns are the analytic second derivative of the sine, so they change sign
every 100 m (exactly; the 1 m phase keeps zero off the nodes). Junction
CSVs, paired flowpaths, a channel network, a hillshade and a slope-area
table (median slope peaking at 1000 m^2) are written in the lsdtt3 formats.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from lsdviztools3.io import ridgeline as rio_

EPSG = 32630
X0, Y0 = 500000.0, 6200000.0
LAMBDAS = (25, 49, 99, 149)
WINDOWS = (25, 49, 99, 149)
NATIVE = (1000, 10000, 100000)
TOPOLOGICAL = (100, 500, 1000)


def _profile_frame(length=3000.0, step=2.0, ridgeline=0, start=0.0):
    from pyproj import Transformer

    d = np.arange(start, start + length + step / 2, step)
    k = 2 * math.pi / 200.0
    wave = 5 * np.sin(k * (d + 1))
    x, y = X0 + d - start, np.full(d.size, Y0)
    lon, lat = Transformer.from_crs(EPSG, 4326, always_xy=True).transform(x, y)
    df = pd.DataFrame({
        "ridgeline_number": ridgeline, "row": 50, "col": (d - start) / step + 10,
        "latitude": lat, "longitude": lon, "ridgeline_type": "interfluvial",
        "distance": d, "elevation": 500 + 0.01 * d + wave, "fitted_elevation": 500 + 0.01 * d,
        "detrended_elevation": wave, "smoothed_detrended_elevation": wave,
        "smoothed_detrended_slope": 5 * k * np.cos(k * (d + 1)),
    })
    for w in WINDOWS:
        df[f"moving_window_fitted_elevation_{w}px"] = 500 + 0.01 * d
        df[f"moving_window_detrended_elevation_{w}px"] = wave / (1 + w / 100)
    for lam in LAMBDAS:
        df[f"maximum_curvature_lambda_{lam}m"] = -k * k * wave / (1 + lam / 100)
        df[f"laplacian_lambda_{lam}m"] = -k * k * wave / (1 + lam / 50)
    return df


def _tdj_frame(distances, ridgeline=0, start=0.0):
    from pyproj import Transformer

    d = np.asarray(distances, dtype=float)
    lon, lat = Transformer.from_crs(EPSG, 4326, always_xy=True).transform(
        X0 + d - start, np.full(d.size, Y0))
    return pd.DataFrame({
        "ridgeline_number": ridgeline, "ridgeline_type": "interfluvial",
        "triple_divide_mode": "x", "node_index": (d * 10).astype(int),
        "distance_along_ridgeline_m": d, "elevation": 500.0, "primary_junction": 1,
        "secondary_junction": 2, "latitude": lat, "longitude": lon})


@pytest.fixture
def outputs(tmp_path):
    """Write the synthetic lsdtt3 outputs; return their paths."""
    import geopandas as gpd
    import rasterio
    from rasterio.transform import from_origin
    from shapely.geometry import LineString

    paths = {"dir": tmp_path}
    profile = pd.concat([_profile_frame(), _profile_frame(ridgeline=1)])
    paths["profile"] = tmp_path / "x_direct_interfluvial_ridgelines_with_polyfit.csv"
    profile.to_csv(paths["profile"], index=False)
    for i, area in enumerate(NATIVE, 1):
        p = tmp_path / f"x_native_triple_divides_area_{area}m2.csv"
        _tdj_frame(np.arange(150.0, 3000.0, 300.0 * i)).to_csv(p, index=False)
        paths[f"native_tdj_csv_{i}"] = p
    for i, dist in enumerate(TOPOLOGICAL, 1):
        p = tmp_path / f"x_topological_distance_triple_divides_{dist}m.csv"
        _tdj_frame(np.arange(250.0, 3000.0, 400.0 * i)).to_csv(p, index=False)
        paths[f"topological_tdj_csv_{i}"] = p

    tdj = _tdj_frame(np.arange(250.0, 3000.0, 400.0))
    lines = []
    for _, row in tdj.iterrows():
        x = X0 + row["distance_along_ridgeline_m"]
        for side, sgn in (("left", 1), ("right", -1)):
            lines.append({"ridgeline_number": 0, "ridgeline_type": "interfluvial",
                          "node_index": int(row["node_index"]), "threshold_m": 100.0,
                          "side": side, "branch": 0, "topological_distance_m": 100.0,
                          "junction_flat": 0,
                          "geometry": LineString([(x, Y0), (x + 20, Y0 + sgn * 150)])})
    paths["paired"] = tmp_path / "x_topological_distance_paired_flowpaths_100m.fgb"
    gpd.GeoDataFrame(lines, crs=EPSG).to_file(paths["paired"])
    channels = [{"id": i, "geometry": LineString([(X0 + 100 * i, Y0 + 30),
                                                  (X0 + 100 * i + 50, Y0 + 300)])}
                for i in range(30)]
    paths["channels"] = tmp_path / "x_channel_network_lines.fgb"
    gpd.GeoDataFrame(channels, crs=EPSG).to_file(paths["channels"])

    cell = 5.0
    west, north = X0 - 500, Y0 + 500
    nx, ny = 800, 200
    yy, xx = np.mgrid[0:ny, 0:nx]
    shade = (127 + 100 * np.sin(xx / 15.0) * np.cos(yy / 11.0)).astype("uint8")
    paths["hillshade"] = tmp_path / "x_hillshade.tif"
    with rasterio.open(paths["hillshade"], "w", driver="GTiff", width=nx, height=ny, count=1,
                       dtype="uint8", crs=f"EPSG:{EPSG}",
                       transform=from_origin(west, north, cell, cell)) as dst:
        dst.write(shade, 1)

    area = np.logspace(1, 6, 26)
    peak = 1.0 - 0.1 * (np.log10(area) - 3.0) ** 2
    sa = pd.DataFrame({"drainage_area": area, "d05_slope": peak * 0.5,
                       "d25_slope": peak * 0.8, "d50_slope": peak, "d75_slope": peak * 1.2,
                       "d95_slope": peak * 1.5, "number_of_pixels": 10})
    paths["slope_area"] = tmp_path / "x_slope_area.csv"
    sa.to_csv(paths["slope_area"], index=False)

    from PIL import Image

    paths["png"] = tmp_path / "x_modelled_log10_drainage_area.png"
    Image.fromarray(np.tile(np.linspace(0, 255, 1600, dtype="uint8"), (80, 1))).save(paths["png"])
    return paths


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


# --- io ---------------------------------------------------------------------------

def test_profile_loading_and_column_helpers(outputs):
    df = rio_.load_ridgeline_profile(outputs["profile"], 1, required=("distance", "elevation"))
    assert (df["ridgeline_number"] == 1).all() and df["distance"].is_monotonic_increasing
    assert rio_.lambda_columns(df.columns) == [
        (float(lam), f"maximum_curvature_lambda_{lam}m") for lam in LAMBDAS]
    assert rio_.lambda_columns(["maximum_curvature"]) == [(0.0, "maximum_curvature")]
    assert rio_.moving_window_sizes(df.columns) == list(WINDOWS)
    assert rio_.ridgeline_id_column(["basin_junction_number", "distance"]) == \
        "basin_junction_number"
    with pytest.raises(ValueError):
        rio_.load_ridgeline_profile(outputs["profile"], 7)


def test_reverse_and_range_selection():
    df = _profile_frame(length=100.0, step=2.0, start=40.0)
    out = rio_.select_profile(df, reverse=True, distance_range=(10.0, 20.0))
    # Reversed: max - (d - min) - min = 140 - d, so 10..20 keeps d = 120..130.
    assert out["display_distance"].tolist() == pytest.approx(np.arange(10.0, 20.1, 2.0))
    assert out["distance"].tolist() == pytest.approx(np.arange(130.0, 119.9, -2.0))
    assert rio_.node_spacing(out["display_distance"]) == pytest.approx(2.0)
    with pytest.raises(ValueError):
        rio_.select_profile(df, distance_range=(20.0, 10.0))
    with pytest.raises(ValueError):
        rio_.node_spacing([0.0, 1.0, 2.0, 10.0])


def test_tdj_thresholds_and_positions_on_a_clipped_reversed_profile(outputs, tmp_path):
    assert rio_.tdj_threshold_from_filename(outputs["native_tdj_csv_2"], "native") == 10000
    assert rio_.tdj_threshold_from_filename(outputs["topological_tdj_csv_3"],
                                            "topological") == 1000
    with pytest.raises(ValueError):
        rio_.tdj_threshold_from_filename(tmp_path / "x_area_1m2_area_2m2.csv", "native")
    tdj = rio_.load_tdj(outputs["native_tdj_csv_1"], "native", 0)
    full = rio_.load_ridgeline_profile(outputs["profile"], 0)["distance"]
    pos = rio_.tdj_display_positions(tdj, full, reverse=True, distance_range=(0.0, 1000.0))
    # Junctions at 150, 450, ... 2850; reversed 3000 - d within 0..1000.
    assert sorted(pos) == pytest.approx([150.0, 450.0, 750.0])


def test_spacing_statistics_on_a_known_sine(outputs):
    df = rio_.load_ridgeline_profile(outputs["profile"], 0)
    spacing = rio_.curvature_spacing_quantiles(df)
    assert spacing["grid_resolution"] == pytest.approx(2.0)
    assert spacing["polyfit_lambda"].tolist() == list(map(float, LAMBDAS))
    # Sign changes every 100 m, located at the node after each crossing.
    for key in ("d05", "d50", "d95"):
        assert spacing[key] == pytest.approx(np.full(4, 100.0))
    absolute = rio_.curvature_spacing_quantiles(df, absolute_curvature=True)
    assert absolute["d95"][0] > absolute["d95"][-1] > 0
    d, c = df["distance"].to_numpy(), df["maximum_curvature_lambda_25m"].to_numpy()
    assert rio_.periodic_median_crossing_spacing(d, c) == pytest.approx(100.0, abs=1e-6)


def test_characteristic_length_and_summary(outputs):
    sa = rio_.load_slope_area(outputs["slope_area"])
    assert rio_.slope_area_peak(sa) == pytest.approx(1000.0)
    assert rio_.characteristic_length(1000.0, 2.0) == pytest.approx(500.0)
    spacing = rio_.curvature_spacing_quantiles(rio_.load_ridgeline_profile(outputs["profile"]))
    row = rio_.spacing_summary("A", spacing, sa)
    assert row["L"] == pytest.approx(500.0) and row["2xL"] == pytest.approx(1000.0)
    assert row["closest_lambda_to_2xL"] == 149.0
    assert row["median_zero_spacing_divided_by_L"] == pytest.approx(0.2)


# --- charts -----------------------------------------------------------------------

def test_curvature_barcode_signs():
    import matplotlib.pyplot as plt

    from lsdviztools3.charts.ridgeline import add_curvature_barcode

    fig, ax = plt.subplots()
    twin = add_curvature_barcode(ax, [0, 1, 2], [[-1.0, 0.0, 2.0], [3.0, -1.0, 0.0]],
                                 scale_labels=[25.0, 50.0])
    image = twin.get_images()[0].get_array()
    assert image.tolist() == [[-1, 0, 1], [1, -1, 0]]
    assert [t.get_text() for t in twin.get_yticklabels()] == ["25", "50"]


def test_charts_render(outputs):
    from lsdviztools3.charts import ridgeline as ch

    full = rio_.load_ridgeline_profile(outputs["profile"], 0)
    profile = rio_.select_profile(full)
    fig, axes = ch.plot_ridgeline_detrending_curvature(profile, WINDOWS, LAMBDAS,
                                                       site_name="Site A")
    assert len(axes) == 4 and len(axes[1].get_lines()) == 5
    native = [rio_.load_tdj(outputs[f"native_tdj_csv_{i}"], "native", 0) for i in (1, 2, 3)]
    topo = [rio_.load_tdj(outputs[f"topological_tdj_csv_{i}"], "topological", 0)
            for i in (1, 2, 3)]
    fig, axes = ch.plot_ridgeline_profile_summary(profile, native, topo,
                                                  highlight=("topological", 100))
    # Topological threshold 100 is highlighted: drawn as open circles, also on the profile.
    assert len(axes["profile"].collections) == 1
    with pytest.raises(ValueError):
        ch.plot_ridgeline_profile_summary(profile, native, topo, highlight=("native", 5))
    image = np.zeros((80, 1600))
    ch.plot_mdd_ridgeline_map_profile(full, image, rio_.divide_grid_spacing(full), "model")
    sa = rio_.load_slope_area(outputs["slope_area"])
    sites = [{"name": n, "spacing": rio_.curvature_spacing_quantiles(full), "slope_area": sa}
             for n in ("A", "B")]
    fig, axes = ch.plot_zero_spacing_slope_area(sites)
    assert axes.shape == (2, 2)


def test_save_figure_formats(tmp_path):
    import matplotlib.pyplot as plt

    from lsdviztools3.charts.ridgeline import save_figure_formats

    fig, _ = plt.subplots()
    written = save_figure_formats(fig, tmp_path / "out", "stem")
    assert [p.relative_to(tmp_path).as_posix() for p in written] == [
        "out/png/stem.png", "out/svg/stem.svg", "out/pdf/stem.pdf"]
    with pytest.raises(ValueError):
        save_figure_formats(fig, tmp_path, "a/b")


# --- map --------------------------------------------------------------------------

def test_ridgeline_ticks_are_perpendicular():
    from lsdviztools3.render.ridgeline import junction_colours, ridgeline_ticks

    d = np.arange(0.0, 2501.0, 2.0)
    ticks = ridgeline_ticks(d, 100 + d, np.full(d.size, 7.0))
    assert len(ticks) == 3
    (x0, y0), (x1, y1) = ticks[1]
    assert (x0, x1) == pytest.approx((1100.0, 1100.0))
    assert sorted((y0, y1)) == pytest.approx([-43.0, 57.0])
    colours = junction_colours([5, 3, 5, 9])
    assert list(colours) == [3, 5, 9] and colours[3] == "#ff0000"


@pytest.mark.parametrize("kind, threshold", [("native", 1000.0), ("topological", 100.0)])
def test_ridgeline_map_renders(outputs, tmp_path, kind, threshold):
    pytest.importorskip("pygmt")
    from lsdviztools3.render.ridgeline import render_ridgeline_map

    full = rio_.load_ridgeline_profile(outputs["profile"], 0)
    rng = (500.0, 2500.0)
    profile = rio_.select_profile(full, True, rng)
    index = (NATIVE if kind == "native" else TOPOLOGICAL).index(int(threshold)) + 1
    tdj = rio_.load_tdj(outputs[f"{kind}_tdj_csv_{index}"], kind, 0)
    overlay = outputs["channels"] if kind == "native" else outputs["paired"]
    fig = render_ridgeline_map(outputs["hillshade"], profile, highlighted_tdj=tdj,
                               overlay=overlay, reverse=True, distance_range=rng,
                               profile_distance=full["distance"])
    out = tmp_path / f"map_{kind}.png"
    fig.savefig(out)
    assert out.stat().st_size > 0


# --- CLIs -------------------------------------------------------------------------

def test_driver_file_parsing(tmp_path):
    from lsdviztools3.cli._ridgeline import read_driver_file

    p = tmp_path / "a.driver"
    p.write_text('# comment\nname: "a # b"  # trailing\nflag: True\nn: 3\nx: 2.5e1\n'
                 'y_axis_label: "Top\\nBottom"\nlist: 25,49\n')
    assert read_driver_file(p) == {"name": "a # b", "flag": True, "n": 3, "x": 25.0,
                                   "y_axis_label": "Top\nBottom", "list": "25,49"}
    p.write_text("a: 1\na: 2\n")
    with pytest.raises(ValueError, match="Duplicate"):
        read_driver_file(p)


def _driver(path, entries):
    path.write_text("".join(f"{k}: {v}\n" for k, v in entries.items()))
    return str(path)


def test_cli_curvature_and_mdd(outputs, capsys):
    from lsdviztools3.cli import plot_mdd_ridgeline, plot_ridgeline_curvature

    out = outputs["dir"] / "figs"
    driver = _driver(outputs["dir"] / "c.driver", {
        "ridgeline_profile_csv": outputs["profile"], "write_path": out,
        "moving_window_sizes_pixels": "25,49,99,149", "curvature_lambdas_m": "25,49,99,149",
        "reverse_ridgeline": "true", "use_distance_range": "true",
        "minimum_distance": 100, "maximum_distance": 2000})
    assert plot_ridgeline_curvature.main(["--driver", driver, "--write-name", "c"]) == 0
    assert (out / "pdf" / "c.pdf").is_file()
    assert plot_mdd_ridgeline.main([
        "--ridgeline-profile-csv", str(outputs["profile"]),
        "--drainage-area-map-png", str(outputs["png"]), "--model-name", "m",
        "--write-path", str(out)]) == 0
    assert (out / "svg" / "mdd_ridgeline_map_profile.svg").is_file()
    assert plot_ridgeline_curvature.main(["--driver", driver,
                                          "--curvature-lambdas-m", "1,2,3,4"]) == 2
    assert "missing column" in capsys.readouterr().err


def test_cli_rejects_unknown_driver_keys(outputs):
    from lsdviztools3.cli import plot_ridgeline_curvature

    driver = _driver(outputs["dir"] / "bad.driver", {"hillshade_orientation": "north_up"})
    with pytest.raises(SystemExit):
        plot_ridgeline_curvature.main(["--driver", driver])


def test_cli_summary(outputs):
    pytest.importorskip("pygmt")
    from lsdviztools3.cli import plot_ridgeline_summary

    entries = {k: v for k, v in outputs.items() if "tdj" in k}
    entries.update({"ridgeline_profile_csv": outputs["profile"],
                    "hillshade_raster_file": outputs["hillshade"],
                    "channel_network_vector_file": outputs["channels"],
                    "paired_flowpaths_vector_file": outputs["paired"],
                    "ridgeline_vector_file": "unused.fgb",
                    "write_path": outputs["dir"] / "summary",
                    "y_axis_label": '"Detrended\\nelevation [m]"',
                    "highlighted_tdj_type": "topological", "highlighted_tdj_threshold": 100})
    driver = _driver(outputs["dir"] / "s.driver", entries)
    assert plot_ridgeline_summary.main(["--driver", driver]) == 0
    root = outputs["dir"] / "summary"
    for rel in ("png/ridgeline_0_profile_summary.png", "svg/ridgeline_0_profile_summary.svg",
                "png/ridgeline_0_profile_summary_map.png",
                "pdf/ridgeline_0_profile_summary_map.pdf"):
        assert (root / rel).is_file(), rel


def test_cli_spacing(outputs):
    from lsdviztools3.cli import plot_ridgeline_spacing

    driver = _driver(outputs["dir"] / "z.driver", {
        "site_count": 2, "write_path": outputs["dir"] / "z",
        "site_1_display_name": "A", "site_1_ridgeline_profile_csv": outputs["profile"],
        "site_1_slope_area_csv": outputs["slope_area"],
        "site_2_display_name": "B", "site_2_ridgeline_profile_csv": outputs["profile"],
        "site_2_slope_area_csv": outputs["slope_area"], "site_2_ridgeline_number": 1,
        "site_2_use_distance_range": "true", "site_2_minimum_distance": 0,
        "site_2_maximum_distance": 1000})
    assert plot_ridgeline_spacing.main(["--driver", driver]) == 0
    summary = pd.read_csv(outputs["dir"] / "z" / "csv" / "spacing_summary_statistics.csv")
    assert summary["site"].tolist() == ["A", "B"]
    assert summary["L"].tolist() == pytest.approx([500.0, 500.0])
    assert (outputs["dir"] / "z" / "png" / "zero_spacing_slope_area.png").is_file()
    with pytest.raises(SystemExit):
        plot_ridgeline_spacing.main(["--write-path", str(outputs["dir"])])


@pytest.mark.parametrize("module", ["plot_ridgeline_curvature", "plot_ridgeline_summary",
                                    "plot_ridgeline_spacing", "plot_mdd_ridgeline"])
def test_cli_help(module):
    import importlib

    mod = importlib.import_module(f"lsdviztools3.cli.{module}")
    with pytest.raises(SystemExit) as exc:
        mod.build_parser().parse_args(["--help"])
    assert exc.value.code == 0
