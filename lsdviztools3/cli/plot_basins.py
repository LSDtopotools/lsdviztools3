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
    # Cosmogenic labels, for erosion-rate maps. Off unless given: they used to
    # default to stone / BraucherBorchers / ERA40 and were written over every
    # basin map, cosmogenic or not.
    parser.add_argument(
        "--scaling-scheme",
        default=None,
        help="Cosmogenic scaling scheme to note on the map (e.g. stone). Not shown unless given.",
    )
    parser.add_argument(
        "--muon-scheme",
        default=None,
        help="Cosmogenic muon scheme to note on the map (e.g. BraucherBorchers). Not shown unless given.",
    )
    parser.add_argument(
        "--atmosphere",
        default=None,
        help="Atmosphere model to note on the map (e.g. ERA40). Not shown unless given.",
    )
    parser.add_argument("--fill", default="orange", help="Polygon fill colour.")
    parser.add_argument(
        "--transparency",
        type=float,
        default=65.0,
        help="Fill transparency percent (0 opaque, 100 invisible).",
    )
    _common.add_locator_args(parser)
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
    metadata = ", ".join(
        f"{name}={value}"
        for name, value in (("scaling", args.scaling_scheme), ("muon", args.muon_scheme),
                            ("atmosphere", args.atmosphere))
        if value
    )
    title_text = " | ".join(t for t in (title_prefix, metadata) if t)
    if title_text:
        # In the map's top-left corner, on a white box so it reads over the relief.
        fig.text(
            position="TL",
            text=title_text,
            font="9p,Helvetica-Bold,black",
            justify="TL",
            offset="0.15c/-0.15c",
            fill="white@15",
            clearance="0.1c/0.1c",
        )

    _common.add_locator_from_args(fig, args, args.base)
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res,
                        directory=args.directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
