"""CLI: lsdviz3-lithology-map — geologic / lithologic unit map to PNG.

Draws an lsdtt3 ``<prefix>_lithology_code.tif`` with a legend of unit names, a
small geologic-map + DEM credit, a solid scale bar (label below) and a plain
north arrow. Two independent options:

* ``--hillshade`` (default) / ``--no-hillshade``: semi-transparent unit colours
  over a gray hillshade of ``--dem``, or flat opaque colours with no hillshade.
* ``--channels PATH``: overlay a channel network (a line/point vector, or the
  chi lithology node CSV with x/y), in one colour or by stream order.

Unit names, the credit and the DEM source are discovered from the lsdtt3
outputs next to the code raster unless given explicitly.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from . import _common


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-lithology-map``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-lithology-map",
        description=(
            "Render an lsdtt3 lithology code raster as a geologic unit map with "
            "unit-name legend, credit line, scale bar and north arrow."
        ),
    )
    parser.add_argument("--codes", dest="code_raster", required=True,
                        help="Unit code raster (<prefix>_lithology_code.tif).")
    parser.add_argument("--dem", default=None,
                        help="DEM for the hillshade (required unless --no-hillshade).")
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument(
        "--prefix", default=None,
        help="lsdtt3 write_fname of the run, used to find unit names and the "
        "attribution (default: the code raster name minus _lithology_code.tif).")
    parser.add_argument(
        "--names", default=None,
        help="CSV giving code->name (code,name or lith_code,lith_name). "
        "Default: <prefix>_lithology_units.csv, then the chi k_sn CSVs, then "
        "--crosswalk, then <prefix>_lithology_legend.csv.")
    parser.add_argument("--crosswalk", default=None,
                        help="Crosswalk CSV with code,name columns (name fallback).")
    parser.add_argument(
        "--attribution", default=None,
        help="lsdtt-fetch-lithology _ATTRIBUTION.txt (default: discovered from "
        "<prefix>_ATTRIBUTION.txt or the run's lithology_fname).")
    parser.add_argument(
        "--credit", default=None,
        help="Explicit credit line; overrides --attribution/--dem-source.")
    parser.add_argument(
        "--dem-source", default=None,
        help="DEM product to credit (e.g. 'Copernicus GLO-30'). Default: read "
        "from the fetch run's <stem>_ingested_params.param beside --dem.")
    parser.add_argument(
        "--dem-params", default=None,
        help="lsdtt-fetch-raster parameter file recording dem_source (default: "
        "<stem>_ingested_params.param beside --dem, if it records one).")
    parser.add_argument("--no-credit", action="store_true",
                        help="Draw no credit line.")
    hs = parser.add_mutually_exclusive_group()
    hs.add_argument("--hillshade", dest="hillshade", action="store_true", default=True,
                    help="Semi-transparent unit colours over a hillshade (default).")
    hs.add_argument("--no-hillshade", dest="hillshade", action="store_false",
                    help="Flat, opaque unit colours with no hillshade.")
    parser.add_argument("--base-is-hillshade", action="store_true",
                        help="--dem is already a hillshade (draw it flat).")
    parser.add_argument("--transparency", type=float, default=35.0,
                        help="Unit-colour transparency over the hillshade, percent "
                        "(default 35).")
    parser.add_argument(
        "--relief", default="illumination", choices=["illumination", "elevation"],
        help="Hillshade base: illumination only (default; one unit = one colour "
        "family) or the gray elevation ramp + illumination (high ground lighter).")
    parser.add_argument(
        "--channels", default=None,
        help="Channel overlay: a line/point vector (.fgb/.gpkg/...) or a CSV with "
        "x,y in the map CRS (e.g. *_chi_ksn_lithology_nodes_mn_*.csv).")
    parser.add_argument(
        "--network-status", action="store_true",
        help="Draw the whole channel network of the CNA run <prefix> "
        "(<prefix>_channel_network_lines.fgb): analysed links thick, the rest thin, "
        "confluences as circles (filled = analysed). Replaces --channels.")
    parser.add_argument("--regions", default=None,
                        help="Region polygons (e.g. USGS physiographic divisions "
                        "physio.shp) drawn as coloured outlines.")
    parser.add_argument("--region-field", default="PROVINCE",
                        help="Region name field in --regions (default PROVINCE).")
    parser.add_argument("--region-title", default="Physiographic province",
                        help="Legend header for the regions (default "
                        "'Physiographic province').")
    parser.add_argument("--channel-color", default="black",
                        help="Channel colour (default black).")
    parser.add_argument("--channel-width", default="0.6p",
                        help="Channel line width for line inputs (default 0.6p).")
    parser.add_argument(
        "--channel-order-column", default=None,
        help="Colour channels by this stream-order column ('auto' to detect "
        "strahler_order/stream_order/order).")
    parser.add_argument(
        "--channel-color-by", default=None, choices=["k_sn"],
        help="Colour channels by k_sn (the k_sn column lsdtt3 wrote in the channel "
        "CSV/vector) with a continuous colour map and colour bar; mutes the geology.")
    parser.add_argument("--ksn-cmap", default="plasma",
                        help="Colour map for k_sn channels (default plasma).")
    parser.add_argument("--ksn-percentiles", default="2,98",
                        help="k_sn colour limits as percentiles 'lo,hi' (default 2,98).")
    parser.add_argument("--ksn-log", action="store_true",
                        help="Colour k_sn on a log10 scale (values <= 0 not drawn).")
    parser.add_argument("--geology-mute", type=float, default=65.0,
                        help="Unit-colour transparency (percent) when channels are "
                        "coloured by k_sn (default 65).")
    parser.add_argument(
        "--unit-colours", default=None,
        help="name,colour CSV or JSON {name: colour} fixing unit colours by NAME "
        "(use the same file for lsdviz3-ksn-by-lithology).")
    parser.add_argument("--legend-title", default=None,
                        help="Legend header (e.g. the lithology field name).")
    parser.add_argument("--legend-columns", type=int, default=None,
                        help="Legend columns (default: 1, or 2 above 16 units).")
    parser.add_argument("--hide-unmapped", action="store_true",
                        help="Leave unmapped (-1) cells undrawn and out of the legend.")
    _common.add_style_args(parser)
    _common.add_target_dir_arg(parser)
    return parser


def resolve_inputs(args) -> dict:
    """Resolve names, credit and prefix from the arguments + discovery.

    Returns:
        ``{"prefix", "names", "credit", "directory"}``.
    """
    from ..io import lithology as L

    code_path = Path(args.code_raster)
    directory = code_path.parent
    prefix = args.prefix
    if prefix is None:
        name = code_path.name
        suffix = L.SUFFIXES["code_raster"]
        prefix = name[: -len(suffix)] if name.endswith(suffix) else code_path.stem

    names = L.resolve_code_names(directory, prefix, names_csv=args.names,
                                 crosswalk=args.crosswalk)

    credit = ""
    if not args.no_credit:
        attribution = args.attribution or L.discover_attribution(directory, prefix)
        dem_params = getattr(args, "dem_params", None) or (
            L.discover_dem_params(args.dem) if args.dem else None)
        credit = L.lithology_map_credit(
            attribution=attribution, credit=args.credit,
            dem_source=args.dem_source, dem_params=dem_params)
    unit_colours = (L.read_unit_colours(args.unit_colours)
                    if getattr(args, "unit_colours", None) else None)

    # m/n and A_0 for a k_sn colour bar: the channel table's own columns are
    # used if present (render side); otherwise the run's summary file.
    ksn_mn = ksn_a0 = None
    if getattr(args, "channel_color_by", None) == "k_sn":
        found = L.find_mn_suffixes(directory, prefix)
        if found:
            summ = L.lithology_paths(directory, prefix, found[0])["summary"]
            if summ.exists():
                s = L.load_ksn_lithology_summary(summ)
                ksn_mn = L._first_value(s, L.SUMMARY_COLUMNS["m_over_n"])
                ksn_a0 = L._first_value(s, L.SUMMARY_COLUMNS["A_0"])
    network_status = None
    if getattr(args, "network_status", False):
        from ..io.channel_bearings import load_network_status

        network_status = load_network_status(directory, prefix)
        if network_status is None:
            raise SystemExit(f"--network-status: no {prefix}_channel_network_lines.* in "
                             f"{directory} (run CNA with write_channel_network_lines: true).")
    return {"prefix": prefix, "names": names, "credit": credit,
            "network_status": network_status,
            "directory": directory, "unit_colours": unit_colours,
            "ksn_m_over_n": ksn_mn, "ksn_a0": ksn_a0}


def main(argv=None) -> int:
    """Entry point: render the lithology map and save a PNG."""
    from ..render.lithology_map import render_lithology_map

    args = build_parser().parse_args(argv)
    if args.hillshade and not args.dem:
        raise SystemExit("--dem is required with --hillshade (or pass --no-hillshade).")
    style = _common.style_from_args(args)
    resolved = resolve_inputs(args)

    fig = render_lithology_map(
        args.code_raster,
        dem=args.dem,
        names=resolved["names"],
        hillshade=args.hillshade,
        shade=not args.base_is_hillshade,
        transparency=args.transparency,
        relief=args.relief,
        channels=args.channels,
        channel_color=args.channel_color,
        channel_width=args.channel_width,
        channel_order_column=args.channel_order_column,
        channel_color_by=args.channel_color_by,
        ksn_cmap=args.ksn_cmap,
        ksn_clip_percentiles=tuple(float(v) for v in args.ksn_percentiles.split(",")),
        ksn_log=args.ksn_log,
        ksn_m_over_n=resolved["ksn_m_over_n"],
        ksn_a0=resolved["ksn_a0"],
        geology_mute_transparency=args.geology_mute,
        network_status=resolved["network_status"],
        regions=args.regions,
        region_field=args.region_field,
        region_header=args.region_title,
        unit_colours=resolved["unit_colours"],
        credit=resolved["credit"],
        legend_title=args.legend_title,
        legend_ncols=args.legend_columns,
        show_unmapped=not args.hide_unmapped,
        coordinates=args.coordinates,
        style=style,
    )
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res,
                        directory=args.directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
