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
        help="Attribute to colour the channels by (e.g. k_sn). Without it the lines are one colour.",
    )
    parser.add_argument("--linear", action="store_true",
                        help="Colour by the value itself rather than its log10 (the default, "
                             "suited to k_sn and drainage area).")
    parser.add_argument("--vmin", type=float, default=None,
                        help="Lower colour limit, in data units.")
    parser.add_argument("--vmax", type=float, default=None,
                        help="Upper colour limit, in data units.")
    parser.add_argument("--colorbar", action="store_true",
                        help="Add a colourbar (continuous colours only; classes get a legend).")
    parser.add_argument("--legend-title", default=None,
                        help="Title of the colourbar or class legend (GMT text).")
    _common.add_class_args(parser, what="the channels (by --color-by)")
    _common.add_style_args(parser)
    _common.add_target_dir_arg(parser)
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
    class_kw = _common.class_kwargs(args)
    if class_kw["classes"] is not None and args.color_by is None:
        raise SystemExit("--classes needs --color-by: the attribute whose values are classified.")
    fig = render_channels(
        args.lines,
        color_by=args.color_by,
        log_color=not args.linear,
        colorbar=args.colorbar,
        colorbar_label=args.legend_title,
        vmin=args.vmin,
        vmax=args.vmax,
        coordinates=args.coordinates,
        style=style,
        fig=fig,
        **class_kw,
    )
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res,
                        directory=args.directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
