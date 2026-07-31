"""CLI: lsdviz3-divide-migration — figures for an ``lsdtt-divide-migration`` run.

Reads one output set (divide polylines, channel heads, pair summary) and writes
a map plus the side-vs-side charts. Nothing is recomputed here: the deltas,
uncertainties and per-head metrics all come from the CLI's own output.
"""

from __future__ import annotations

import argparse

from . import _common
from ..io.divide_migration import METRIC_CHOICES
from ..workflows import render_divide_migration_figures


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-divide-migration``.

    Returns:
        The configured parser (``--prefix`` output set, optional ``--base``
        raster for the map, ``--outdir``, ``--metric``, style flags).
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-divide-migration",
        description=(
            "Render lsdtt-divide-migration outputs: a hillshade map with divide "
            "polylines and channel heads coloured by a metric, plus the "
            "side-vs-side pair-summary chart (delta = A - B, A the "
            "lower-numbered unit) and a per-head scatter against distance to "
            "the divide."
        ),
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Output stem (<write_path>/<write_fname>) or a directory holding "
        "one lsdtt-divide-migration output set.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for the figure PNGs.",
    )
    parser.add_argument(
        "--base",
        default=None,
        help="Base raster (DEM, or a hillshade with --base-is-hillshade) for "
        "the map. Without it only the charts are written.",
    )
    parser.add_argument(
        "--base-is-hillshade",
        action="store_true",
        help="Treat --base as a precomputed hillshade rather than elevation.",
    )
    parser.add_argument(
        "--metric",
        default="chi",
        choices=sorted(METRIC_CHOICES),
        help="Channel-head metric to map and compare (default: chi).",
    )
    parser.add_argument(
        "--point-size",
        default=None,
        help="GMT symbol size for channel heads (e.g. 0.12c).",
    )
    parser.add_argument(
        "--divide-color",
        default="black",
        help="Pen colour for the divide polylines (default: black).",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="p-value threshold, used only when the run's uncertainty method "
        "is ttest (default: 0.05).",
    )
    parser.add_argument(
        "--stem",
        default="divide_migration",
        help="Output filename stem (default: divide_migration).",
    )
    parser.add_argument(
        "--chart-dpi",
        type=int,
        default=220,
        help="DPI for the matplotlib charts (default: 220).",
    )
    _common.add_style_args(parser)
    _common.add_target_dir_arg(parser)
    # The shared --cmap default ("geo") is a terrain palette; a divide metric
    # wants a sequential/diverging scientific map instead. The hillshade base is
    # forced to grey by the renderer regardless.
    parser.set_defaults(cmap="roma")
    return parser


def main(argv=None) -> int:
    """Entry point: render the divide-migration figure set.

    Args:
        argv: Optional argument list.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)

    from ..paths import resolve_in_dir

    outdir = resolve_in_dir(args.directory, args.outdir)

    written = render_divide_migration_figures(
        args.prefix,
        outdir,
        base_raster=args.base,
        metric=args.metric,
        coordinates=args.coordinates,
        base_is_hillshade=args.base_is_hillshade,
        cmap=args.cmap,
        point_size=args.point_size,
        divide_color=args.divide_color,
        alpha=args.alpha,
        stem=args.stem,
        style=style,
        high_res=args.high_res,
        dpi=args.chart_dpi,
    )
    for key, value in written.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
