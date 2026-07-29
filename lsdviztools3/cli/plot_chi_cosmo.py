"""CLI: lsdviz3-chi-cosmo — chi basin products + optional cosmo comparison."""

from __future__ import annotations

import argparse
from dataclasses import replace

from ..paths import resolve_in_dir
from . import _common
from ..workflows import render_chi_and_cosmo_basin_workflow


def build_parser() -> argparse.ArgumentParser:
    """Build parser for ``lsdviz3-chi-cosmo``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-chi-cosmo",
        description=(
            "Render one chi profile PNG per basin, chi stacked plot, basin k_sn "
            "summary/map, and optional k_sn-vs-cosmo erosion-rate comparison."
        ),
    )
    parser.add_argument("--fit-nodes-csv", required=True, help="Path to *_chi_fit_nodes_mn_*.csv")
    parser.add_argument("--fit-segments-csv", required=True, help="Path to *_chi_fit_segments_mn_*.csv")
    parser.add_argument("--outdir", required=True, help="Output directory for generated files.")
    parser.add_argument("--stem", default="chi", help="Output filename stem.")
    parser.add_argument("--dpi", type=int, default=300, help="Output PNG dpi.")
    parser.add_argument("--cmap", default="viridis", help="Colormap for basin overlays.")
    parser.add_argument("--theta", type=float, default=None, help="Concavity index theta (m/n) for labels.")
    parser.add_argument("--scaling-scheme", default="stone", help="Metadata label.")
    parser.add_argument("--atmosphere", default="ERA40", help="Metadata label.")
    parser.add_argument("--muon-scheme", default="BraucherBorchers", help="Metadata label.")

    parser.add_argument("--base", default=None, help="Optional base raster for basin k_sn map.")
    parser.add_argument("--hillshade", default=None, help="Optional lsdtt3 hillshade raster path.")
    parser.add_argument("--basins", default=None, help="Optional basin polygons for basin k_sn map.")
    parser.add_argument("--basin-id-col", default=None, help="Optional basin ID column in basin polygons.")

    parser.add_argument("--cosmo-csv", default=None, help="Optional cosmo erosion-rate CSV for k_sn-vs-erate plot.")
    parser.add_argument("--erate-col", default=None, help="Optional erosion-rate column in cosmo CSV.")
    parser.add_argument("--lon-col", default=None, help="Optional longitude column in cosmo CSV.")
    parser.add_argument("--lat-col", default=None, help="Optional latitude column in cosmo CSV.")
    parser.add_argument("--sample-col", default=None, help="Optional sample-name column in cosmo CSV.")
    parser.add_argument(
        "--erate-input-units",
        default="auto",
        choices=["auto", "mm_per_yr", "mm_per_kyr"],
        help="Units of erosion-rate values in --cosmo-csv.",
    )
    parser.add_argument(
        "--allow-incomplete-basins",
        action="store_true",
        help="Allow basin polygons that touch DEM edge/nodata buffers.",
    )

    parser.add_argument(
        "--dir",
        dest="directory",
        default=None,
        help="Target directory for outputs (created if absent). A relative --outdir is written inside it.",
    )
    return parser


def main(argv=None) -> int:
    """Entry point for chi/cosmo basin workflow plotting."""
    args = build_parser().parse_args(argv)
    outdir = resolve_in_dir(args.directory, args.outdir)
    style = _common.style_from_args(args)
    if style.cmap == "geo":
        style = replace(style, cmap="viridis")

    out = render_chi_and_cosmo_basin_workflow(
        {
            "fit_nodes": args.fit_nodes_csv,
            "fit_segments": args.fit_segments_csv,
        },
        outdir,
        stem=args.stem,
        dpi=args.dpi,
        basin_polygons=args.basins,
        base_raster=args.base,
        basin_id_col=args.basin_id_col,
        cosmo_csv=args.cosmo_csv,
        erate_col=args.erate_col,
        lon_col=args.lon_col,
        lat_col=args.lat_col,
        sample_col=args.sample_col,
        theta=args.theta,
        scaling_scheme=args.scaling_scheme,
        atmosphere=args.atmosphere,
        muon_scheme=args.muon_scheme,
        erate_input_units=args.erate_input_units,
        hillshade_raster=args.hillshade,
        require_complete_basins=not args.allow_incomplete_basins,
        style=style,
    )

    for key, val in out.items():
        if isinstance(val, list):
            print(f"{key}:")
            for p in val:
                print(f"  - {p}")
        else:
            print(f"{key}: {val}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
