"""Single-channel basins: basins draining into one trunk channel (Nujiang).

Demonstrates lsdtt3's `extract_basins_from_single_channel` mode: from an upstream
point (and an optional downstream point) it builds ONE main channel by following
D8 (upstream identified by reachability; truncated at the common-downstream point
where the downstream point's flow joins), then extracts the side-tributary basins
draining into that channel (complete + minimum-size filtered). All of that is
done in lsdtt3 (C++); this script orchestrates and plots.

Test site: a stretch of the Nujiang (Nu/Salween) gorge around
26.835, 98.880. Outputs under ../data/nujang. lsdtt3 via $LSDTT3_DIR.

Run:
    pixi run python examples/single_channel_basins.py
    pixi run python examples/single_channel_basins.py --skip-compute
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import geopandas as gpd
from shapely.geometry import Point

from lsdviztools3.driver import run_lsdtt
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.basins import render_basins
from lsdviztools3.render.channels import render_channels
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.render.points import render_points
from lsdviztools3.render.style import MapStyle

ROOT = Path(os.environ.get("NUJANG_DIR", "/mnt/c/LSDTopoTools/data/nujang")).resolve()
LSDTT_DIR = os.environ.get("LSDTT3_DIR", "/mnt/c/LSDTopoTools/lsdtt3/build")
STEM = "nujang_COP30"
MIN_BASIN_PIXELS = 2000     # min side-basin size (~1.8 km2 on 30 m COP30)

DEM = ROOT / f"{STEM}.tif"
FILL = ROOT / f"{STEM}_fill.tif"
HILLSHADE = ROOT / f"{STEM}_hillshade.tif"
POINTS = ROOT / "nujang_mainstem_points.gpkg"   # upstream + downstream points
CHANNEL = ROOT / "nujang_sc_side_tributary_main_channel.gpkg"
BASINS = ROOT / "nujang_sc_side_tributary_basin_outlines.gpkg"
OUTLETS = ROOT / "nujang_sc_side_tributary_outlets.gpkg"
PLOTS = ROOT / "plots"


def prepare_points() -> None:
    """Place upstream/downstream points on the Nujiang main stem (test-data prep).

    In normal use the user supplies the points file. Here, for reproducibility,
    we put two points on the main stem by drainage area (extract the channel
    network, then pick the highest-area point near two latitudes).
    """
    if POINTS.exists():
        return
    run_lsdtt("lsdtt-channel-extraction",
              {"read_path": str(ROOT), "write_path": str(ROOT),
               "read_fname": f"{STEM}_fill.tif", "write_fname": "nujang",
               "raster_is_filled": "true", "source_pixels_threshold": 5000,
               "write_channel_network_points": "true",
               "channel_network_points_output_format": "gpkg"},
              directory=ROOT, lsdtt_dir=LSDTT_DIR)
    pts = gpd.read_file(ROOT / "nujang_channel_network_points.gpkg").to_crs(4326)
    area = [c for c in pts.columns if "area" in c.lower()][0]

    def main_stem_near(lat, half=0.02):
        sub = pts[(pts.geometry.y > lat - half) & (pts.geometry.y < lat + half)]
        return sub.loc[sub[area].idxmax()].geometry

    up, dn = main_stem_near(26.90), main_stem_near(26.77)
    gpd.GeoDataFrame({"role": ["upstream", "downstream"]},
                     geometry=[Point(up.x, up.y), Point(dn.x, dn.y)],
                     crs=4326).to_file(POINTS, driver="GPKG")


def compute() -> None:
    """Preprocess the DEM and run the single-channel side-basin extraction."""
    print("Preprocessing (fill + hillshade) ...")
    run_lsdtt("lsdtt-raster-preprocessing",
              {"read_path": str(ROOT), "write_path": str(ROOT),
               "read_fname": f"{STEM}.tif", "write_fname": STEM,
               "write_fill_raster": "true", "write_hillshade_raster": "true"},
              directory=ROOT, lsdtt_dir=LSDTT_DIR)

    prepare_points()

    print("Extracting single channel + side basins ...")
    run_lsdtt("lsdtt-basin-extraction",
              {"read_path": str(ROOT), "write_path": str(ROOT),
               "read_fname": f"{STEM}_fill.tif", "write_fname": "nujang_sc",
               "raster_is_filled": "true",
               "extract_basins_from_single_channel": "true",
               "single_channel_points_fname": str(POINTS),
               "side_tributary_minimum_contributing_pixels": MIN_BASIN_PIXELS,
               "reject_incomplete_side_tributary_basins": "true",
               "write_side_tributary_main_channel": "true",
               "write_side_tributary_outlets": "true",
               "write_side_tributary_basin_outlines": "true",
               "side_tributary_vector_output_format": "gpkg"},
              directory=ROOT, lsdtt_dir=LSDTT_DIR)


def plot() -> None:
    """Render the single channel and the basins draining into it."""
    PLOTS.mkdir(parents=True, exist_ok=True)
    style = MapStyle(cmap="gray", clean=True, transparent=False)
    fig = render_hillshade(str(HILLSHADE), coordinates="local", shade=False, style=style)
    fig = render_basins(str(BASINS), coordinates="local", fill="orange", transparency=55,
                        style=MapStyle(line_width="0.7p", line_color="saddlebrown"), fig=fig)
    fig = render_channels(str(CHANNEL), coordinates="local",
                          style=MapStyle(line_width="2p", line_color="blue"), fig=fig)
    fig = render_points(str(OUTLETS), coordinates="local",
                        style=MapStyle(line_color="yellow", point_size="0.14c"), fig=fig)
    fig = render_points(str(POINTS), coordinates="local",
                        style=MapStyle(line_color="red", point_size="0.30c"), fig=fig)
    out = PLOTS / "single_channel_basins.png"
    save_figure(fig, out, style=MapStyle(transparent=False))
    print(f"wrote {out}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-compute", action="store_true",
                    help="Skip lsdtt3; plot from existing outputs.")
    args = ap.parse_args(argv)
    if not args.skip_compute:
        compute()
    plot()
    print(f"\nDone. Outputs under {ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
