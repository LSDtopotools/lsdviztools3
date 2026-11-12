"""CLI: lsdviz3-basemap — a locator map of a raster (globe, region or aerial view) to PNG."""

from __future__ import annotations

import argparse

from . import _common
from ..render.basemap import render_orthographic_basemap, render_regional_basemap


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-basemap``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-basemap",
        description="A locator map showing where a raster is: a globe, its region (borders, "
                    "rivers, optional relief), or an aerial perspective view, with the raster's "
                    "footprint outlined.",
    )
    parser.add_argument("--in", dest="input_raster", required=True, help="Input raster path.")
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument(
        "--mode",
        choices=("auto", "global", "regional", "aerial"),
        default="auto",
        help="global: a hemisphere centred on the raster; regional: a flat map of the region "
             "around it; aerial: a perspective view from above. auto (default): global for "
             "large rasters, regional otherwise.",
    )
    parser.add_argument(
        "--elevation-km",
        type=float,
        default=None,
        help="Observer elevation in km for aerial mode.",
    )
    parser.add_argument(
        "--view-limit-deg",
        type=float,
        default=None,
        help="Half-angle view limit in degrees for aerial mode.",
    )
    parser.add_argument(
        "--look-bearing",
        type=float,
        default=0.0,
        help="Look bearing in degrees clockwise from north (default: 0).",
    )
    parser.add_argument(
        "--look-vertical-angle",
        type=float,
        default=0.0,
        help="Vertical look angle in degrees (default: 0; straight down).",
    )
    parser.add_argument(
        "--grid-spacing-deg",
        type=float,
        default=None,
        help="Graticule spacing in degrees.",
    )
    parser.add_argument(
        "--no-north-arrow",
        action="store_true",
        help="Disable north arrow (regional and aerial maps; a globe has none).",
    )
    parser.add_argument(
        "--extent-km", type=float, default=None,
        help="Regional mode: width of the map in km (default: ten times the raster's "
             "diagonal, at least 250 km).",
    )
    parser.add_argument(
        "--relief", action="store_true",
        help="Regional mode: shade the land with GMT's global relief data (downloaded from "
             "the GMT server on first use).",
    )
    parser.add_argument("--no-borders", action="store_true",
                        help="Regional mode: leave out country and state borders.")
    parser.add_argument("--no-rivers", action="store_true",
                        help="Regional mode: leave out rivers.")
    parser.add_argument("--no-scale-bar", action="store_true",
                        help="Regional mode: leave out the scale bar.")
    _common.add_style_args(parser)
    _common.add_target_dir_arg(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: render a basemap from raster footprint."""
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)
    if args.mode == "regional":
        fig = render_regional_basemap(
            args.input_raster,
            style=style,
            extent_km=args.extent_km,
            relief=args.relief,
            borders=not args.no_borders,
            rivers=not args.no_rivers,
            scale_bar=not args.no_scale_bar,
            show_north_arrow=not args.no_north_arrow,
        )
    else:
        fig = render_orthographic_basemap(
            args.input_raster,
            style=style,
            mode=args.mode,
            elevation_km=args.elevation_km,
            view_limit_deg=args.view_limit_deg,
            look_bearing=args.look_bearing,
            look_vertical_angle=args.look_vertical_angle,
            grid_spacing_deg=args.grid_spacing_deg,
            show_north_arrow=not args.no_north_arrow,
            extent_km=args.extent_km,
            relief=args.relief,
        )
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res,
                        directory=args.directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
