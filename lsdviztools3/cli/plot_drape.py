"""CLI: lsdviz3-drape — drape a coloured raster over a gray hillshade to PNG."""

from __future__ import annotations

import argparse

from . import _common
from ..render.channels import render_drape


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-drape``.

    Returns:
        The configured parser (``--base`` + ``--drape`` rasters, ``--out`` PNG,
        drape options, and the shared style flags).
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-drape",
        description=(
            "Drape a coloured, semi-transparent data raster over a gray "
            "hillshade base and render to PNG. The base is always gray; the "
            "drape uses --cmap."
        ),
    )
    parser.add_argument(
        "--base",
        dest="base_raster",
        required=True,
        help="Base raster path (a precomputed hillshade, or a DEM with --shade-base).",
    )
    parser.add_argument(
        "--drape",
        dest="drape_raster",
        required=True,
        help="Data raster to colour-drape on top.",
    )
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument(
        "--shade-base",
        action="store_true",
        help="Treat --base as elevation and illuminate it, rather than as a "
        "precomputed hillshade.",
    )
    parser.add_argument(
        "--log",
        action="store_true",
        help="Colour the drape by log10 of its values (values <= 0 are clamped "
        "to transparent).",
    )
    parser.add_argument(
        "--transparency",
        type=float,
        default=40.0,
        help="Drape layer transparency, 0 (opaque) to 100 (invisible); "
        "default 40.",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Fixed lower bound for the drape colour scale.",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Fixed upper bound for the drape colour scale.",
    )
    parser.add_argument(
        "--percentile-min",
        type=float,
        default=None,
        help="Lower bound as a percentile (0-100) of the drape values, e.g. 90.",
    )
    parser.add_argument(
        "--percentile-max",
        type=float,
        default=None,
        help="Upper bound as a percentile (0-100) of the drape values.",
    )
    parser.add_argument(
        "--out-of-bounds",
        default="clamp",
        choices=["clamp", "clip"],
        help="Out-of-range handling: 'clamp' to the end colour (default), or "
        "'clip' to leave unrendered (transparent).",
    )
    parser.add_argument(
        "--categorical-gt-breaks",
        default=None,
        help=(
            "Comma-separated thresholds for categorical rendering using strictly "
            "greater-than classes (e.g. '0.05,0.12,0.21'). Values <= first "
            "threshold remain nodata/transparent."
        ),
    )
    parser.add_argument(
        "--categorical-colors",
        default=None,
        help=(
            "Comma-separated colours matching --categorical-gt-breaks "
            "(e.g. 'yellow,orange,red')."
        ),
    )
    parser.add_argument(
        "--colorbar",
        action="store_true",
        help="Add a colourbar for the drape colour scale.",
    )
    parser.add_argument(
        "--colorbar-label",
        default=None,
        help="Colourbar label (e.g. 'Curvature').",
    )
    parser.add_argument(
        "--colorbar-position",
        default="bottom",
        choices=["bottom", "top", "left", "right"],
        help="Colourbar placement relative to the map (default: bottom).",
    )
    _common.add_class_args(parser, what="the drape")
    _common.add_locator_args(parser)
    _common.add_style_args(parser)
    _common.add_target_dir_arg(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: render the gray hillshade + coloured drape, save PNG.

    Args:
        argv: Optional argument list (defaults to ``sys.argv``).

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)
    cat_breaks = None
    cat_colors = None
    if args.categorical_gt_breaks:
        cat_breaks = [float(v.strip()) for v in args.categorical_gt_breaks.split(",") if v.strip()]
    if args.categorical_colors:
        cat_colors = [v.strip() for v in args.categorical_colors.split(",") if v.strip()]

    fig = render_drape(
        args.base_raster,
        args.drape_raster,
        coordinates=args.coordinates,
        base_is_hillshade=not args.shade_base,
        log_color=args.log,
        transparency=args.transparency,
        cmap=style.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        percentile_min=args.percentile_min,
        percentile_max=args.percentile_max,
        out_of_bounds=args.out_of_bounds,
        categorical_gt_breaks=cat_breaks,
        categorical_colors=cat_colors,
        **_common.class_kwargs(args),
        colorbar=args.colorbar,
        colorbar_label=args.colorbar_label,
        colorbar_position=args.colorbar_position,
        style=style,
    )
    _common.add_locator_from_args(fig, args, args.base_raster)
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res,
                        directory=args.directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
