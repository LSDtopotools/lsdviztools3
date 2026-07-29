"""CLI: lsdviz3-cosmo-basins — standard basin erosion-rate overlay plot."""

from __future__ import annotations

import argparse
from dataclasses import replace

from . import _common
from ..paths import resolve_in_dir
from ..workflows import render_cosmo_basin_erate_map


def build_parser() -> argparse.ArgumentParser:
    """Build parser for ``lsdviz3-cosmo-basins``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-cosmo-basins",
        description=(
            "Join cosmo erosion-rate CSV samples to basin polygons and render a "
            "hillshade + semi-transparent basin overlay colored by inferred rate."
        ),
    )
    parser.add_argument("--base", required=True, help="Base raster path (DEM or hillshade).")
    parser.add_argument(
        "--hillshade",
        default=None,
        help="Optional lsdtt3 hillshade raster path. When set, this is used as the base image.",
    )
    parser.add_argument("--basins", required=True, help="Basin polygon vector path.")
    parser.add_argument("--cosmo-csv", required=True, help="Cosmo erosion-rate CSV path.")
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument(
        "--joined-basins-out",
        default=None,
        help="Optional path to write joined basin polygons (e.g. .fgb/.gpkg).",
    )
    parser.add_argument("--basin-id-col", default=None, help="Optional basin ID column.")
    parser.add_argument("--erate-col", default=None, help="Optional erosion-rate CSV column.")
    parser.add_argument("--lon-col", default=None, help="Optional longitude CSV column.")
    parser.add_argument("--lat-col", default=None, help="Optional latitude CSV column.")
    parser.add_argument("--sample-col", default=None, help="Optional sample-name CSV column.")
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
        "--out-erate-col",
        default="inferred_erate_mm_yr",
        help="Joined basin attribute name for erosion rate.",
    )
    parser.add_argument(
        "--plot-title",
        default="",
        help="Optional title text drawn at the top-left.",
    )
    parser.add_argument("--scaling-scheme", default="stone", help="Metadata label.")
    parser.add_argument("--muon-scheme", default="BraucherBorchers", help="Metadata label.")
    parser.add_argument("--atmosphere", default="ERA40", help="Metadata label.")
    parser.add_argument("--fill", default="orange", help="Fallback basin fill colour.")
    parser.add_argument(
        "--transparency",
        type=float,
        default=65.0,
        help="Fill transparency percent (0 opaque, 100 invisible).",
    )
    parser.add_argument(
        "--colorbar",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Draw erosion-rate colorbar (default: true).",
    )
    parser.add_argument(
        "--colorbar-label",
        default="Inferred erosion rate (mm/yr)",
        help="Colorbar label.",
    )
    _common.add_style_args(parser)
    _common.add_target_dir_arg(parser)
    return parser


def main(argv=None) -> int:
    """Entry point for standardized basin-erate plot generation."""
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)
    if style.cmap == "geo":
        style = replace(style, cmap="viridis")
    out_png = resolve_in_dir(args.directory, args.out)
    joined_out = (
        resolve_in_dir(args.directory, args.joined_basins_out)
        if args.joined_basins_out
        else None
    )

    out = render_cosmo_basin_erate_map(
        base_raster=args.base,
        basin_polygons=args.basins,
        cosmo_csv=args.cosmo_csv,
        out_png=out_png,
        joined_basins_out=joined_out,
        basin_id_col=args.basin_id_col,
        erate_col=args.erate_col,
        lon_col=args.lon_col,
        lat_col=args.lat_col,
        sample_col=args.sample_col,
        out_erate_col=args.out_erate_col,
        erate_input_units=args.erate_input_units,
        hillshade_raster=args.hillshade,
        scaling_scheme=args.scaling_scheme,
        muon_scheme=args.muon_scheme,
        atmosphere=args.atmosphere,
        plot_title=args.plot_title,
        coordinates=args.coordinates,
        transparency=args.transparency,
        fill=args.fill,
        colorbar=args.colorbar,
        colorbar_label=args.colorbar_label,
        style=style,
        require_complete_basins=not args.allow_incomplete_basins,
    )

    for key, value in out.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
