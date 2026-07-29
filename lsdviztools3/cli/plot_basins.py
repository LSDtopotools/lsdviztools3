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
    parser.add_argument(
        "--color-by",
        default=None,
        help="Optional numeric basin column used to color polygons.",
    )
    parser.add_argument(
        "--colorbar",
        action="store_true",
        help="Draw a colorbar when --color-by is used.",
    )
    parser.add_argument(
        "--colorbar-label",
        default="Inferred erosion rate (mm/yr)",
        help="Colorbar label when --color-by is used.",
    )
    parser.add_argument(
        "--plot-title",
        default="",
        help="Optional title text drawn at the top-left.",
    )
    parser.add_argument(
        "--scaling-scheme",
        default="stone",
        help="Metadata label for scaling scheme (default: stone).",
    )
    parser.add_argument(
        "--muon-scheme",
        default="BraucherBorchers",
        help="Metadata label for muon scheme (default: BraucherBorchers).",
    )
    parser.add_argument(
        "--atmosphere",
        default="ERA40",
        help="Metadata label for atmosphere dataset (default: ERA40).",
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
        shade=True,
        style=base_style,
    )
    fig = render_basins(
        args.basins,
        label_col=args.label_col,
        coordinates=args.coordinates,
        fill=args.fill,
        transparency=args.transparency,
        color_by=args.color_by,
        colorbar=args.colorbar,
        colorbar_label=args.colorbar_label,
        cmap=style.cmap,
        style=style,
        fig=fig,
    )

    title_prefix = args.plot_title.strip()
    metadata = (
        f"scaling={args.scaling_scheme}, "
        f"muon={args.muon_scheme}, atmosphere={args.atmosphere}"
    )
    title_text = metadata if not title_prefix else f"{title_prefix} | {metadata}"
    if title_text:
        fig.text(
            x=0.01,
            y=0.99,
            text=title_text,
            font="14p,Helvetica-Bold,black",
            region="0/1/0/1",
            projection="X1c",
            justify="TL",
            no_clip=True,
        )

    _common.save_figure(fig, args.out, style=style, high_res=args.high_res,
                        directory=args.directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
