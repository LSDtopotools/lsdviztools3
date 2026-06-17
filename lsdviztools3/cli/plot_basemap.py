"""CLI: lsdviz3-basemap — render global/aerial basemap to PNG."""

from __future__ import annotations

import argparse

from . import _common
from ..render.basemap import render_orthographic_basemap


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-basemap``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-basemap",
        description="Render a global orthographic or aerial perspective basemap.",
    )
    parser.add_argument("--in", dest="input_raster", required=True, help="Input raster path.")
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument(
        "--mode",
        choices=("auto", "global", "aerial"),
        default="auto",
        help="Basemap mode (default: auto).",
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
        help="Disable north arrow.",
    )
    _common.add_style_args(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: render a basemap from raster footprint."""
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)
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
    )
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
