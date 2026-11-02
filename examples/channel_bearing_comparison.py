"""Channel-bearing comparison: dendritic vs trellis drainage (COP30 + lsdtt3).

Downloads two COP30 tiles — an Allegheny Plateau site (dendritic) and a Valley
and Ridge site (trellis, strongly structurally aligned) — fills + hillshades
each, extracts the channel network with junction/bearing analysis restricted to
*complete* basins, then plots:

  a) channel bearings as a rose diagram separated by Strahler stream order,
  b) junction opening angles (distribution, the two sites compared),
  c) outlet -> basin-centroid bearings as a rose diagram,
  and a circular (radial) comparison of the two sites' channel-bearing
  distributions.

The statistics of that comparison come from lsdtt3, not from this script: one
lsdtt-channel-network-analysis run over a DEM covering both sites, with each
site's box as a region (channel_bearing_group_by: region), writes the circular
Wasserstein distance, the Kuiper test and a block-permutation p-value between
the two regions (<stem>_channel_bearing_comparisons.csv). The comparison plot
draws the links of that same run (<stem>_channel_bearing_group_membership.csv).

Outputs under ../data/channel_bearing_comparison (one subdir per site, plus a
comparison/ dir). lsdtt3 binaries via $LSDTT3_DIR (default
/mnt/c/LSDTopoTools/lsdtt3/build).

Run:
    pixi run python examples/channel_bearing_comparison.py
    pixi run python examples/channel_bearing_comparison.py --skip-compute
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from lsdviztools3.charts import (
    plot_circular_comparison, plot_distributions, rose_diagram,
)
from lsdviztools3.download.ot_scraper import download_dem
from lsdviztools3.driver import run_lsdtt
from lsdviztools3.paths import ensure_dir
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.render.style import MapStyle

ROOT = Path(os.environ.get(
    "CHANNEL_BEARING_DIR", "/mnt/c/LSDTopoTools/data/channel_bearing_comparison"))
LSDTT_DIR = os.environ.get("LSDTT3_DIR", "/mnt/c/LSDTopoTools/lsdtt3/build")

HALF_DEG = 0.10        # half-width of the (square) bounding box in degrees
SRC_THRESH = 500       # threshold_contributing_pixels (~0.45 km^2 on 30 m COP30)


@dataclass
class Site:
    key: str
    name: str
    lat: float
    lon: float
    network: str


SITES = [
    Site("allegheny", "Allegheny Plateau", 39.603215, -80.526961, "dendritic"),
    Site("valley_ridge", "Valley and Ridge", 39.278199, -78.666881, "trellis"),
]


def _stem(site: Site) -> str:
    return f"{site.key}_COP30"


def _paths(site: Site) -> dict:
    d = ROOT / site.key
    s = _stem(site)
    return {
        "dir": d,
        "dem": d / f"{s}.tif",
        "fill": d / f"{s}_fill.tif",
        "hillshade": d / f"{s}_hillshade.tif",
        "links": d / f"{s}_channel_link_bearings.csv",
        "net_lines": d / f"{s}_channel_network_lines.fgb",
        "angles": d / f"{s}_junction_angles.csv",
        "com": d / f"{s}_junction_com_bearings.csv",
        "hillshade_png": d / f"{site.key}_hillshade.png",
        "channels_png": d / f"{site.key}_hillshade_channels.png",
    }


def compute_site(site: Site) -> dict:
    """Download + preprocess + channel-network analysis for one site."""
    p = _paths(site)
    site_dir = ensure_dir(p["dir"])
    stem = _stem(site)
    base = {"read_path": str(site_dir), "write_path": str(site_dir)}

    if not p["dem"].exists():
        print(f"[{site.key}] downloading COP30 ...")
        download_dem(
            (site.lon - HALF_DEG, site.lat - HALF_DEG,
             site.lon + HALF_DEG, site.lat + HALF_DEG),
            source="planetary_computer", out_file=f"{stem}.tif",
            directory=site_dir,
        )

    print(f"[{site.key}] preprocessing (fill + hillshade) ...")
    run_lsdtt(
        "lsdtt-raster-preprocessing",
        dict(base, read_fname=f"{stem}.tif", write_fname=stem,
             write_fill_raster="true", write_hillshade_raster="true"),
        directory=site_dir, lsdtt_dir=LSDTT_DIR,
    )

    print(f"[{site.key}] channel-network analysis (complete basins) ...")
    run_lsdtt(
        "lsdtt-channel-network-analysis",
        dict(base, read_fname=f"{stem}_fill.tif", write_fname=stem,
             raster_is_filled="true", threshold_contributing_pixels=SRC_THRESH,
             analyse_all_junctions="true", only_complete_basins="true",
             write_channel_link_bearings="true", write_junction_angles="true",
             write_junction_centre_of_mass_bearings="true",
             write_vector_format="csv"),
        directory=site_dir, lsdtt_dir=LSDTT_DIR,
    )

    print(f"[{site.key}] channel-extraction (network lines for plotting) ...")
    run_lsdtt(
        "lsdtt-channel-extraction",
        dict(base, read_fname=f"{stem}_fill.tif", write_fname=stem,
             raster_is_filled="true", threshold_contributing_pixels=SRC_THRESH,
             write_channel_network_lines="true"),
        directory=site_dir, lsdtt_dir=LSDTT_DIR,
    )
    return p


def render_site(site: Site, p: dict) -> None:
    """Clean grey hillshade PNG, and a hillshade + channel-network overlay."""
    from lsdviztools3.render.channels import render_channels

    print(f"[{site.key}] rendering hillshade + channel network ...")
    style = MapStyle(cmap="gray", clean=True, transparent=False)
    fig = render_hillshade(str(p["hillshade"]), coordinates="local",
                           shade=False, style=style)
    save_figure(fig, p["hillshade_png"], style=style)

    # Hillshade with the channel network drawn on top.
    fig = render_hillshade(str(p["hillshade"]), coordinates="local",
                           shade=False,
                           style=MapStyle(cmap="gray", clean=True, transparent=False))
    fig = render_channels(
        str(p["net_lines"]), coordinates="local",
        style=MapStyle(line_width="0.8p", line_color="dodgerblue4"), fig=fig)
    save_figure(fig, p["channels_png"], style=MapStyle(transparent=False))


def _savefig(fig, out: Path) -> None:
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")


def plot_site(site: Site, p: dict) -> None:
    """Per-site roses: channel bearings by stream order, and outlet->com bearings."""
    import matplotlib.pyplot as plt

    links = pd.read_csv(p["links"])
    groups = {f"order {o}": links.loc[links.stream_order == o, "bearing_deg"].to_numpy()
              for o in sorted(links.stream_order.unique())}
    fig, _ = rose_diagram(groups=groups, bin_width_deg=10, mode="directional",
                          stacked=True, label_fmt="{}",
                          title=f"{site.name}: channel bearings by order")
    _savefig(fig, p["dir"] / f"{site.key}_rose_channel_bearings.png")
    plt.close(fig)

    com = pd.read_csv(p["com"])
    fig, _ = rose_diagram(com["bearing_deg"].to_numpy(), bin_width_deg=15,
                          mode="directional", stacked=True,
                          title=f"{site.name}: outlet -> centroid bearings")
    _savefig(fig, p["dir"] / f"{site.key}_rose_com_bearings.png")
    plt.close(fig)


def plot_comparison(sites_paths) -> None:
    """Cross-site comparisons into comparison/."""
    import matplotlib.pyplot as plt

    comp = ensure_dir(ROOT / "comparison")

    # Junction opening angles (b).
    angle_series = {
        site.name: pd.read_csv(p["angles"])["donor_opening_angle_deg"].to_numpy()
        for site, p in sites_paths
    }
    fig, _ = plot_distributions(
        angle_series, bins=30, value_range=(0, 180),
        xlabel="Junction opening angle (deg)",
        title="Junction angles: dendritic vs trellis")
    _savefig(fig, comp / "junction_angle_comparison.png")
    plt.close(fig)

    # Channel-bearing distributions: circular comparison. The statistics are
    # lsdtt3's (one run over both sites, grouped by region); the plot draws the
    # links of that run and annotates lsdtt3's numbers.
    (site_a, _), (site_b, _) = sites_paths
    members = pd.read_csv(comp / f"{COMPARE_STEM}_channel_bearing_group_membership.csv")
    bearings_a = members.loc[members.group_value == site_a.name, "bearing_deg"].to_numpy()
    bearings_b = members.loc[members.group_value == site_b.name, "bearing_deg"].to_numpy()
    fig, _ = plot_circular_comparison(
        bearings_a, bearings_b, labels=(site_a.name, site_b.name),
        bin_width_deg=10, kind="both", annotation=_comparison_annotation(comp),
        title="Channel-bearing distributions")
    _savefig(fig, comp / "channel_bearing_circular_comparison.png")
    plt.close(fig)


COMPARE_STEM = "both_sites_COP30"


def _comparison_box() -> tuple:
    """A WGS84 box (west, south, east, north) covering every site's box."""
    pad = HALF_DEG + 0.02
    return (min(s.lon for s in SITES) - pad, min(s.lat for s in SITES) - pad,
            max(s.lon for s in SITES) + pad, max(s.lat for s in SITES) + pad)


def compute_comparison() -> None:
    """One lsdtt3 run over both sites, comparing their bearings as regions."""
    import geopandas as gpd
    from shapely.geometry import box

    comp = ensure_dir(ROOT / "comparison")
    dem = comp / f"{COMPARE_STEM}.tif"
    if not dem.exists():
        print("[comparison] downloading COP30 over both sites ...")
        download_dem(_comparison_box(), source="planetary_computer",
                     out_file=dem.name, directory=comp)

    # Each site's box is a region; lsdtt3 reprojects the polygons to the DEM.
    regions = comp / "site_regions.gpkg"
    gpd.GeoDataFrame(
        {"site": [s.name for s in SITES]},
        geometry=[box(s.lon - HALF_DEG, s.lat - HALF_DEG, s.lon + HALF_DEG, s.lat + HALF_DEG)
                  for s in SITES],
        crs="EPSG:4326").to_file(regions, driver="GPKG")

    base = {"read_path": str(comp), "write_path": str(comp)}
    print("[comparison] preprocessing ...")
    run_lsdtt(
        "lsdtt-raster-preprocessing",
        dict(base, read_fname=dem.name, write_fname=COMPARE_STEM, write_fill_raster="true"),
        directory=comp, lsdtt_dir=LSDTT_DIR,
    )
    print("[comparison] channel-network analysis, bearings compared by region ...")
    run_lsdtt(
        "lsdtt-channel-network-analysis",
        dict(base, read_fname=f"{COMPARE_STEM}_fill.tif", write_fname=COMPARE_STEM,
             raster_is_filled="true", threshold_contributing_pixels=SRC_THRESH,
             analyse_all_junctions="true", only_complete_basins="true",
             write_channel_link_bearings="true", write_vector_format="csv",
             write_channel_bearing_statistics="true",
             channel_bearing_group_by="region",
             channel_region_fname=regions.name, channel_region_field_name="site",
             channel_bearing_comparison_mode="pairwise"),
        directory=comp, lsdtt_dir=LSDTT_DIR,
    )


def _comparison_annotation(comp: Path) -> str:
    """lsdtt3's statistics for the directional bearing comparison of the sites."""
    rows = pd.read_csv(comp / f"{COMPARE_STEM}_channel_bearing_comparisons.csv")
    rows = rows[(rows.measurand == "link_bearing") & (rows.angular_type == "directional")]
    if "all" in set(rows.stratum_value.astype(str)):
        rows = rows[rows.stratum_value.astype(str) == "all"]
    row = rows.iloc[0]
    return (f"circular $W_1$ = {row.circular_wasserstein_deg:.1f}$\\degree$\n"
            f"Kuiper V = {row.kuiper_V:.3f}, p (iid) = {row.kuiper_p_iid_approx:.3g}\n"
            f"block permutation p = {row.permutation_p:.3g} "
            f"({int(row.n_blocks_a)} / {int(row.n_blocks_b)} blocks)\n"
            f"(lsdtt-channel-network-analysis)")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-compute", action="store_true",
                        help="Skip download/lsdtt3; replot from existing CSVs.")
    args = parser.parse_args(argv)

    ensure_dir(ROOT)
    sites_paths = []
    for site in SITES:
        p = _paths(site) if args.skip_compute else compute_site(site)
        if not args.skip_compute:
            render_site(site, p)
        plot_site(site, p)
        sites_paths.append((site, p))

    if not args.skip_compute:
        compute_comparison()
    plot_comparison(sites_paths)
    print(f"\nDone. Outputs under {ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
