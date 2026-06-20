"""Parkfield: basin geometry vs the San Andreas Fault (a baseline line).

End-to-end example of comparing basin geometry to a baseline structure (a fault).
All *analysis* is done in lsdtt3 (C++); this script only orchestrates and plots.

It:
  - downloads a COP30 tile around Parkfield, CA,
  - fills + hillshades it and clips the global GEM active-fault database to the
    DEM footprint (lsdtt-raster-preprocessing),
  - extracts basins whose channels CROSS the fault with drainage area above a
    threshold, keeping only the largest complete basin per crossing, and writes
    a per-basin analysis CSV with outlet, centre-of-mass, area, outlet->COM
    bearing, local fault azimuth, side of the fault, and the bearing relative to
    the fault-orthogonal — all computed in lsdtt-basin-extraction
    (select_basins_from_line_crossings; the fault is reprojected internally),
  - then plots:
      a) hillshade + fault trace,
      b) + basins (semi-transparent),
      c) + basin centre-of-mass points and outlet -> COM vectors,
      d) rose of outlet -> COM bearings with the fault direction marked,
      e) outlet -> COM bearing relative to the fault-orthogonal, split by side
         of the fault (does dextral motion bend the basins?).

Outputs under ../data/parkfield/plots. lsdtt3 via $LSDTT3_DIR.

Run:
    pixi run python examples/fault_basin_deformation.py
    pixi run python examples/fault_basin_deformation.py --skip-compute
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString

from lsdviztools3.charts import plot_distributions, rose_diagram
from lsdviztools3.download.ot_scraper import download_dem
from lsdviztools3.driver import run_lsdtt
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.basins import render_basins
from lsdviztools3.render.channels import render_channels
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.render.points import render_points
from lsdviztools3.render.style import MapStyle

ROOT = Path(os.environ.get("PARKFIELD_DIR",
                           "/mnt/c/LSDTopoTools/data/parkfield")).resolve()
FAULTS_GLOBAL = Path(os.environ.get(
    "GEM_FAULTS",
    "/mnt/c/LSDTopoTools/data/fault_vector/gem_active_faults_harmonized.gpkg")
).resolve()
LSDTT_DIR = os.environ.get("LSDTT3_DIR", "/mnt/c/LSDTopoTools/lsdtt3/build")

PARK_LAT, PARK_LON = 35.89865076816407, -120.44109502276305
HALF_DEG = 0.20
SRC_THRESH = 200        # channel source threshold (pixels)
CROSS_MIN = 500         # min contributing pixels at a fault crossing (~0.45 km2)
MAX_AREA_KM2 = 20.0     # cap: keep fault-transverse catchments, drop the
                        # along-fault trunk valley (the SAF runs *along* it)
STEM = "parkfield_COP30"

PLOTS = ROOT / "plots"
FAULT = ROOT / f"{STEM}_clipped_vector.gpkg"           # clipped faults (EPSG:4326)
BASINS = ROOT / "parkfield_junction_basin_polygons.fgb"
ANALYSIS = ROOT / "parkfield_line_crossing_basins.csv"  # all analysis (from C++)


def compute() -> None:
    """Download + preprocess + clip + extract crossing basins (all C++ analysis)."""
    ROOT.mkdir(parents=True, exist_ok=True)
    if not (ROOT / f"{STEM}.tif").exists():
        print("Downloading Parkfield COP30 ...")
        download_dem((PARK_LON - HALF_DEG, PARK_LAT - HALF_DEG,
                      PARK_LON + HALF_DEG, PARK_LAT + HALF_DEG),
                     source="planetary_computer", out_file=f"{STEM}.tif",
                     directory=ROOT)

    print("Preprocessing (fill + hillshade) and clipping faults ...")
    run_lsdtt("lsdtt-raster-preprocessing",
              {"read_path": str(ROOT), "write_path": str(ROOT),
               "read_fname": f"{STEM}.tif", "write_fname": STEM,
               "write_fill_raster": "true", "write_hillshade_raster": "true",
               "clip_vector_data_to_raster": "true",
               "vector_data_fname": str(FAULTS_GLOBAL),
               "vector_clip_output_format": "gpkg"},
              directory=ROOT, lsdtt_dir=LSDTT_DIR)

    print("Extracting basins crossing the fault (+ per-basin analysis) ...")
    # The fault stays in its native CRS; basin-extraction reprojects it
    # internally and computes every per-basin quantity in C++.
    run_lsdtt("lsdtt-basin-extraction",
              {"read_path": str(ROOT), "write_path": str(ROOT),
               "read_fname": f"{STEM}_fill.tif", "write_fname": "parkfield",
               "raster_is_filled": "true", "source_pixels_threshold": SRC_THRESH,
               "select_basins_from_line_crossings": "true",
               "crossing_line_fname": str(FAULT),
               "crossing_minimum_contributing_pixels": CROSS_MIN,
               "filter_basins_by_area": "true",
               "basins_max_area_pixels": int(MAX_AREA_KM2 * 1e6 / 900.0),
               "write_basin_outlines_as_polygon": "true"},
              directory=ROOT, lsdtt_dir=LSDTT_DIR)


def load_outputs():
    """Load the C++ analysis CSV + basin polygons; build GeoDataFrames to plot."""
    df = pd.read_csv(ANALYSIS)
    basins = gpd.read_file(BASINS)
    crs = basins.crs
    com_gdf = gpd.GeoDataFrame(
        df.copy(), geometry=gpd.points_from_xy(df.com_x, df.com_y), crs=crs)
    vec_gdf = gpd.GeoDataFrame(
        df.copy(),
        geometry=[LineString([(r.outlet_x, r.outlet_y), (r.com_x, r.com_y)])
                  for r in df.itertuples()],
        crs=crs)
    fault_az = float(df["local_fault_azimuth_deg"].median())   # from C++
    return basins, df, com_gdf, vec_gdf, fault_az


# --- plots (the only thing this script computes) ----------------------------

def _base_hillshade():
    style = MapStyle(cmap="gray", clean=True, transparent=False)
    return render_hillshade(str(ROOT / f"{STEM}_hillshade.tif"),
                            coordinates="local", shade=False, style=style), style


def make_plots(basins, df, com_gdf, vec_gdf, fault_az):
    import matplotlib.pyplot as plt
    PLOTS.mkdir(parents=True, exist_ok=True)
    fault_pen = MapStyle(line_width="2p", line_color="red")

    # a) hillshade + fault
    fig, style = _base_hillshade()
    fig = render_channels(str(FAULT), coordinates="local", style=fault_pen, fig=fig)
    save_figure(fig, PLOTS / "a_hillshade_fault.png", style=style)

    # b) + basins (semi-transparent)
    fig, _ = _base_hillshade()
    fig = render_basins(basins, coordinates="local", fill="skyblue", transparency=55,
                        style=MapStyle(line_width="0.8p", line_color="navy"), fig=fig)
    fig = render_channels(str(FAULT), coordinates="local", style=fault_pen, fig=fig)
    save_figure(fig, PLOTS / "b_hillshade_fault_basins.png",
                style=MapStyle(transparent=False))

    # c) + COM points and outlet -> COM vectors
    fig, _ = _base_hillshade()
    fig = render_basins(basins, coordinates="local", fill="skyblue", transparency=70,
                        style=MapStyle(line_width="0.6p", line_color="navy"), fig=fig)
    fig = render_channels(str(FAULT), coordinates="local", style=fault_pen, fig=fig)
    fig = render_channels(vec_gdf, coordinates="local",
                          style=MapStyle(line_width="1.4p", line_color="yellow"), fig=fig)
    fig = render_points(com_gdf, coordinates="local",
                        style=MapStyle(line_color="black", point_size="0.18c"), fig=fig)
    save_figure(fig, PLOTS / "c_outlet_com_vectors.png",
                style=MapStyle(transparent=False))

    # d) rose of outlet -> COM bearings with fault direction
    fig, ax = rose_diagram(df["outlet_to_com_bearing_deg"].to_numpy(),
                           bin_width_deg=15, mode="directional",
                           title="Outlet → COM bearings (with fault direction)")
    rmax = ax.get_ylim()[1]
    for a in (fault_az, fault_az + 180):
        ax.plot([math.radians(a), math.radians(a)], [0, rmax],
                color="red", lw=2.2, ls="--")
    ax.text(math.radians(fault_az), rmax, "  fault", color="red", fontsize=9)
    fig.savefig(PLOTS / "d_rose_bearings_fault.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # e) deflection from fault-orthogonal, by side (all computed in C++)
    col = "deflection_from_fault_orthogonal_deg"
    series = {f"{s} side (n={int((df.side == s).sum())})": df.loc[df.side == s, col].to_numpy()
              for s in ("east", "west")}
    fig, ax = plot_distributions(
        series, bins=18, value_range=(-90, 90),
        xlabel="Outlet→COM bearing minus fault-orthogonal (deg)",
        title="Basin axis vs fault-orthogonal (dextral bending?)")
    ax.axvline(0, color="0.3", ls="--", lw=1)
    ax.text(0, ax.get_ylim()[1] * 0.95, " orthogonal to fault", fontsize=8,
            rotation=90, va="top", ha="left", color="0.3")
    fig.savefig(PLOTS / "e_deflection_by_side.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    for n in ["a_hillshade_fault", "b_hillshade_fault_basins",
              "c_outlet_com_vectors", "d_rose_bearings_fault",
              "e_deflection_by_side"]:
        print(f"wrote {PLOTS / (n + '.png')}")
    print(f"\nfault azimuth (C++): {fault_az:.1f} deg ; "
          f"east basins: {(df.side == 'east').sum()}, "
          f"west basins: {(df.side == 'west').sum()}")
    print(f"median deflection  east: {df.loc[df.side=='east', col].median():.1f} deg, "
          f"west: {df.loc[df.side=='west', col].median():.1f} deg")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-compute", action="store_true",
                    help="Skip download/lsdtt3; plot from existing outputs.")
    args = ap.parse_args(argv)
    if not args.skip_compute:
        compute()
    make_plots(*load_outputs())
    print(f"\nDone. Plots under {PLOTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
