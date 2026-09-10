"""CLI: lsdviz3-hillshade — render a hillshade raster to PNG."""

from __future__ import annotations

import argparse

from . import _common
from ..render.hillshade import render_hillshade


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-hillshade``.

    Returns:
        The configured parser (``--in`` raster, ``--out`` PNG, style flags).
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-hillshade",
        description="Render a hillshade/elevation raster to PNG.",
    )
    parser.add_argument(
        "--in",
        dest="input_raster",
        required=True,
        help="Input raster path.",
    )
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument(
        "--input-is-hillshade",
        action="store_true",
        help="Treat input as a pre-computed hillshade and disable extra shading.",
    )
    _common.add_style_args(parser)
    _common.add_target_dir_arg(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: load raster, render hillshade, save PNG.

    Args:
        argv: Optional argument list (defaults to ``sys.argv``).

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)
    # A pre-computed hillshade is illumination, not elevation: colouring it with
    # the elevation colour map paints slope aspect in browns and greens. Grey
    # unless the caller asked for something else.
    if args.input_is_hillshade:
        from dataclasses import replace
        from ..render.style import MapStyle
        if args.cmap == MapStyle().cmap:
            style = replace(style, cmap="gray")
    fig = render_hillshade(
        args.input_raster,
        coordinates=args.coordinates,
        shade=not args.input_is_hillshade,
        style=style,
    )
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res,
                        directory=args.directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
