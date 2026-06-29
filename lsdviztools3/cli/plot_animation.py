"""CLI: lsdviz3-animate — animate a transient lsdtt-lem run as a GIF.

Globs the ``_model_iteration_<N>`` rasters in a directory (the lsdtt-lem
``output/`` folder), renders each on one shared colour scale with an iteration
and model-time label (read from ``_model_elevation_stats.csv``), and stitches
them into an animated GIF.
"""

from __future__ import annotations

import argparse

from . import _common
from ..render.animation import DEFAULT_FRAME_GLOB, animate_transient_rasters


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-animate``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-animate",
        description="Animate a transient lsdtt-lem raster sequence to a GIF.",
    )
    parser.add_argument(
        "--in",
        dest="frames_dir",
        required=True,
        help="Directory holding the transient frame rasters.",
    )
    parser.add_argument("--out", required=True, help="Output GIF path.")
    parser.add_argument(
        "--glob",
        default=DEFAULT_FRAME_GLOB,
        help=f"Frame filename pattern (default: {DEFAULT_FRAME_GLOB}).",
    )
    parser.add_argument(
        "--stats-csv",
        default=None,
        help="Optional _model_elevation_stats.csv for time labels "
        "(auto-detected in --in by default).",
    )
    parser.add_argument(
        "--fps", type=float, default=4.0, help="Frames per second (default: 4)."
    )
    parser.add_argument(
        "--loop", type=int, default=0, help="GIF loop count (0 = forever)."
    )
    parser.add_argument(
        "--no-annotate",
        action="store_true",
        help="Do not stamp iteration/time labels on the frames.",
    )
    parser.add_argument(
        "--keep-frames",
        action="store_true",
        help="Keep the intermediate PNG frames next to the GIF.",
    )
    _common.add_style_args(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: build frames, assemble GIF."""
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)
    out = animate_transient_rasters(
        args.frames_dir,
        args.out,
        glob=args.glob,
        stats_csv=args.stats_csv,
        coordinates=args.coordinates,
        style=style,
        annotate=not args.no_annotate,
        fps=args.fps,
        loop=args.loop,
        keep_frames=args.keep_frames,
    )
    print(f"Wrote animation: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
