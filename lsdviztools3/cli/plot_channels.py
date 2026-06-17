"""CLI: lsdviz3-channels — hillshade base + draped channel network to PNG."""

from __future__ import annotations

import argparse
from dataclasses import replace

from . import _common
from ..render.channels import render_channels
from ..render.hillshade import render_hillshade


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-channels``.

    Returns:
        The configured parser (``--base`` raster, ``--lines`` vector,
        ``--out`` PNG, ``--color-by``, style flags).
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-channels",
        description="Render a hillshade base with channel-network overlays.",
    )
    parser.add_argument("--base", required=True, help="Base raster path.")
    parser.add_argument("--lines", required=True, help="Channel vector path.")
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument(
        "--color-by",
        default=None,
        help="Optional attribute used to colour channels (API placeholder).",
    )
    _common.add_style_args(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: render hillshade base then drape channels, save PNG.

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
    fig = render_channels(
        args.lines,
        color_by=args.color_by,
        coordinates=args.coordinates,
        style=style,
        fig=fig,
    )
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
