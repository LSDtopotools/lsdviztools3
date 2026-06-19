"""CLI: lsdviz3-basins — hillshade base + basin polygons to PNG."""

from __future__ import annotations

import argparse
from dataclasses import replace

from . import _common
from ..render.basins import render_basins
from ..render.hillshade import render_hillshade


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-basins``.

    Returns:
        The configured parser (``--base`` raster, ``--basins`` vector,
        ``--label-col``, ``--out`` PNG, style flags).
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-basins",
        description="Render a hillshade base with basin polygon overlays.",
    )
    parser.add_argument("--base", required=True, help="Base raster path.")
    parser.add_argument("--basins", required=True, help="Basin vector path.")
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument(
        "--label-col",
        default=None,
        help="Optional column to label polygons.",
    )
    parser.add_argument("--fill", default="orange", help="Polygon fill colour.")
    parser.add_argument(
        "--transparency",
        type=float,
        default=65.0,
        help="Fill transparency percent (0 opaque, 100 invisible).",
    )
    _common.add_style_args(parser)
    _common.add_target_dir_arg(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: render hillshade base then basin polygons, save PNG.

    Args:
        argv: Optional argument list.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)
    base_style = replace(style, cmap="gray")

    fig = render_hillshade(
        args.base,
        coordinates=args.coordinates,
        shade=False,
        style=base_style,
    )
    fig = render_basins(
        args.basins,
        label_col=args.label_col,
        coordinates=args.coordinates,
        fill=args.fill,
        transparency=args.transparency,
        style=style,
        fig=fig,
    )
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res,
                        directory=args.directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
