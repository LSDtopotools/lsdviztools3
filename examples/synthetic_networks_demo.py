"""Synthetic drainage-network demo: generate, render, and compare networks.

Drives lsdtt3 to build synthetic Optimal Channel Networks (OCNs), renders each
as a clean hillshade with its channel network overlaid, and compares the
distributions of junction angles and junction area ratios across cases.

Four cases: {rectangle, island} x {more-branching, less-branching}, all on a
30 m grid. Branchiness is set by the simulated-annealing "heat loss" schedule of
``lsdtt-synthetic-networks`` (higher ``temperature_scale`` / faster
``cooling_rate`` / lower ``energy_exponent_gamma`` => more branching).

Outputs go under ``../data/synthetic_networks`` (one subdir per case, plus a
``comparison/`` dir). lsdtt3 binaries are located via ``$LSDTT3_DIR`` (default:
``/mnt/c/LSDTopoTools/lsdtt3/build``).

Run:
    pixi run python examples/synthetic_networks_demo.py
    pixi run python examples/synthetic_networks_demo.py --only-plots   # re-plot
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from lsdviztools3.driver import run_lsdtt
from lsdviztools3.paths import ensure_dir
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.channels import render_channels
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.render.style import MapStyle

ROOT = Path(os.environ.get(
    "SYNTHETIC_NETWORKS_DIR", "/mnt/c/LSDTopoTools/data/synthetic_networks"))
LSDTT_DIR = os.environ.get("LSDTT3_DIR", "/mnt/c/LSDTopoTools/lsdtt3/build")

GRID = 256            # n_rows = n_cols
CELL = 30             # cell size in metres
SRC_THRESH = 50       # threshold_contributing_pixels (same for extraction + analysis)
ISLAND_FILL = 0.92    # island_fill_fraction: bigger island => more channels


@dataclass
class Case:
    """One synthetic-network case: a domain shape and a branchiness setting."""

    key: str
    domain_type: str          # "rectangle" or "efd_island"
    temperature_scale: float
    cooling_rate: float
    gamma: float
    shape: str                # "rectangle" | "island" (for plot grouping)
    branching: str            # "more" | "less"


# More branching: hot, fast-cooling, weak aggregation. Less: cool, slow, strong.
_MORE = dict(temperature_scale=3.0, cooling_rate=2.0, gamma=0.4)
_LESS = dict(temperature_scale=0.5, cooling_rate=0.5, gamma=0.7)

CASES = [
    Case("rectangle_more", "rectangle", shape="rectangle", branching="more", **_MORE),
    Case("rectangle_less", "rectangle", shape="rectangle", branching="less", **_LESS),
    Case("island_more", "efd_island", shape="island", branching="more", **_MORE),
    Case("island_less", "efd_island", shape="island", branching="less", **_LESS),
]


def run_case(case: Case) -> dict:
    """Generate the OCN DEM, extract channels, and compute junction stats.

    Returns a dict of the produced file Paths.
    """
    case_dir = ensure_dir(ROOT / case.key)
    stem = case.key
    base = {"write_path": str(case_dir), "write_fname": stem}
    dem = f"{stem}_ocn_dem.tif"

    print(f"[{case.key}] synthetic-networks ({case.domain_type}, "
          f"T={case.temperature_scale}) ...")
    syn = dict(base,
               n_rows=GRID, n_cols=GRID, cell_size=CELL,
               domain_type=case.domain_type,
               temperature_scale=case.temperature_scale,
               cooling_rate=case.cooling_rate,
               energy_exponent_gamma=case.gamma,
               steady_state_mode="relief", target_relief=1000,
               m_exponent=0.45, n_exponent=1.0, baselevel_elevation=0.0,
               random_seed=1, write_dem_raster="true")
    if case.domain_type != "rectangle":
        syn["island_fill_fraction"] = ISLAND_FILL
    run_lsdtt("lsdtt-synthetic-networks", syn,
              directory=case_dir, lsdtt_dir=LSDTT_DIR)

    print(f"[{case.key}] channel-extraction ...")
    run_lsdtt(
        "lsdtt-channel-extraction",
        dict(base, read_path=str(case_dir), read_fname=dem,
             threshold_contributing_pixels=SRC_THRESH,
             write_channel_network_lines="true"),
        directory=case_dir, lsdtt_dir=LSDTT_DIR,
    )

    print(f"[{case.key}] channel-network-analysis (junction stats) ...")
    run_lsdtt(
        "lsdtt-channel-network-analysis",
        dict(base, read_path=str(case_dir), read_fname=dem,
             threshold_contributing_pixels=SRC_THRESH,
             analyse_all_junctions="true",
             write_junction_angles="true", write_vector_format="csv",
             write_junction_area_ratio_data="true"),
        directory=case_dir, lsdtt_dir=LSDTT_DIR,
    )

    return _case_paths(case)


def _case_paths(case: Case) -> dict:
    """The expected output paths for a case (whether or not it has been run)."""
    case_dir = ROOT / case.key
    stem = case.key
    return {
        "dir": case_dir,
        "dem": case_dir / f"{stem}_ocn_dem.tif",
        "lines": case_dir / f"{stem}_channel_network_lines.fgb",
        "angles": case_dir / f"{stem}_junction_angles.csv",
        "ratios": case_dir / f"{stem}_junction_area_ratios.csv",
        "hillshade_png": case_dir / f"{stem}_hillshade.png",
        "channels_png": case_dir / f"{stem}_hillshade_channels.png",
    }


def render_case(case: Case, paths: dict) -> None:
    """Render a clean shaded-relief PNG and a hillshade+channels PNG."""
    dem = str(paths["dem"])
    print(f"[{case.key}] rendering hillshade + channels ...")

    # 1) Clean shaded relief (grey, scale bar, no furniture).
    style = MapStyle(cmap="gray", clean=True, transparent=False)
    fig = render_hillshade(dem, coordinates="local", shade=True, style=style)
    save_figure(fig, paths["hillshade_png"], style=style)

    # 2) Hillshade + channel network overlaid.
    fig = render_hillshade(
        dem, coordinates="local", shade=True,
        style=MapStyle(cmap="gray", clean=True, transparent=False),
    )
    try:
        fig = render_channels(
            str(paths["lines"]), coordinates="local",
            style=MapStyle(line_width="0.8p", line_color="dodgerblue4"),
            fig=fig,
        )
    except (ValueError, FileNotFoundError) as exc:
        print(f"  ! skipping channel overlay for {case.key}: {exc}")
    save_figure(fig, paths["channels_png"], style=MapStyle(transparent=False))


# --- Comparison distribution plots (plain matplotlib) ---------------------

# Colour by branchiness, line style by shape, so the legend reads the 2x2 grid.
_COLOR = {"more": "crimson", "less": "steelblue"}
_LS = {"rectangle": "-", "island": "--"}


def _label(case: Case) -> str:
    return f"{case.shape} / {case.branching}-branching"


def _style_for(case: Case) -> dict:
    return {"color": _COLOR[case.branching], "linestyle": _LS[case.shape]}


def _series_by_label(by_case: dict, kind: str, column: str) -> tuple[dict, dict]:
    """Collect ``{label: values}`` and ``{label: style}`` for the four cases."""
    series, styles = {}, {}
    for case in CASES:
        df = by_case.get(case.key, {}).get(kind)
        if df is None or column not in df or df.empty:
            continue
        series[_label(case)] = df[column].to_numpy()
        styles[_label(case)] = _style_for(case)
    return series, styles


def plot_angle_comparison(by_case: dict, out: Path) -> None:
    """Overlaid junction-opening-angle distributions across the four cases."""
    from lsdviztools3.charts import plot_distributions

    series, styles = _series_by_label(by_case, "angles", "donor_opening_angle_deg")
    fig, _ = plot_distributions(
        series, bins=30, value_range=(0, 180), styles=styles,
        xlabel="Junction opening angle (deg)",
        title="Junction-angle distributions",
    )
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")


def plot_ratio_comparison(by_case: dict, out: Path) -> None:
    """Overlaid junction area-ratio distributions (log10) across the cases."""
    from lsdviztools3.charts import plot_distributions

    series, styles = _series_by_label(by_case, "ratios", "ratio")
    fig, _ = plot_distributions(
        series, bins=30, value_range=(-3.0, 0.0), log_x=True, styles=styles,
        xlabel="log$_{10}$(junction area ratio = tributary / main)",
        title="Junction area-ratio distributions",
    )
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")


def _load_csv(path: Path):
    try:
        return pd.read_csv(path)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only-plots", action="store_true",
                        help="Skip lsdtt3/rendering; rebuild comparison plots "
                             "from existing CSVs.")
    args = parser.parse_args(argv)

    ensure_dir(ROOT)
    by_case = {}
    for case in CASES:
        if args.only_plots:
            paths = _case_paths(case)
        else:
            paths = run_case(case)
            render_case(case, paths)
        by_case[case.key] = {
            "angles": _load_csv(paths["angles"]),
            "ratios": _load_csv(paths["ratios"]),
        }

    comp = ensure_dir(ROOT / "comparison")
    plot_angle_comparison(by_case, comp / "junction_angle_distributions.png")
    plot_ratio_comparison(by_case, comp / "junction_area_ratio_distributions.png")
    print(f"\nDone. Outputs under {ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
