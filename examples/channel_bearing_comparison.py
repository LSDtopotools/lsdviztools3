"""Channel-bearing comparison: dendritic vs trellis drainage (COP30 + lsdtt3).

Downloads two COP30 tiles — an Allegheny Plateau site (dendritic) and a Valley
and Ridge site (trellis, strongly structurally aligned) — fills + hillshades
each, extracts the channel network with junction/bearing analysis restricted to
*complete* basins, then plots:

  a) channel bearings as a rose diagram separated by Strahler stream order,
  b) junction opening angles (distribution, the two sites compared),
  c) outlet -> basin-centroid bearings as a rose diagram,
  and a circular (radial) comparison of the two sites' channel-bearing
  distributions with a circular-Wasserstein distance + Kuiper test.

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
import subprocess
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
SRC_THRESH = 500       # source_pixels_threshold (~0.45 km^2 on 30 m COP30)


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
        "angles": d / f"{s}_junction_angles.csv",
        "com": d / f"{s}_junction_com_bearings.csv",
        "hillshade_png": d / f"{site.key}_hillshade.png",
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
    analysis = dict(
        base, read_fname=f"{stem}_fill.tif", write_fname=stem,
        raster_is_filled="true", source_pixels_threshold=SRC_THRESH,
        analyse_all_junctions="true",
        write_channel_link_bearings="true", channel_link_bearings_output_format="csv",
        write_junction_angles="true", junction_angles_output_format="csv",
        write_junction_centre_of_mass_bearings="true",
        only_complete_basins="true",
    )
    try:
        run_lsdtt("lsdtt-channel-network-analysis", analysis,
                  directory=site_dir, lsdtt_dir=LSDTT_DIR)
    except subprocess.CalledProcessError:
        # Older binary without the complete-basin filter: retry without it.
        print(f"  ! only_complete_basins not supported; retrying without it")
        analysis.pop("only_complete_basins")
        run_lsdtt("lsdtt-channel-network-analysis", analysis,
                  directory=site_dir, lsdtt_dir=LSDTT_DIR)
    return p


def render_site(site: Site, p: dict) -> None:
    """Clean grey hillshade PNG from the precomputed hillshade raster."""
    print(f"[{site.key}] rendering hillshade ...")
    style = MapStyle(cmap="gray", clean=True, transparent=False)
    fig = render_hillshade(str(p["hillshade"]), coordinates="local",
                           shade=False, style=style)
    save_figure(fig, p["hillshade_png"], style=style)


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

    # Channel-bearing distributions: circular comparison + Wasserstein.
    (site_a, pa), (site_b, pb) = sites_paths
    bearings_a = pd.read_csv(pa["links"])["bearing_deg"].to_numpy()
    bearings_b = pd.read_csv(pb["links"])["bearing_deg"].to_numpy()
    fig, _ = plot_circular_comparison(
        bearings_a, bearings_b, labels=(site_a.name, site_b.name),
        bin_width_deg=10, kind="both", rotation_invariant=False, annotate=True,
        title="Channel-bearing distributions")
    _savefig(fig, comp / "channel_bearing_circular_comparison.png")
    plt.close(fig)


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

    plot_comparison(sites_paths)
    print(f"\nDone. Outputs under {ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
