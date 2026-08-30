"""CLI for paired Newton/cellular hillslope comparison animations."""

from __future__ import annotations

import argparse

from ..render.hillslope_comparison_animation import (
    render_hillslope_comparison_animation,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lsdviz3-hillslope-comparison",
        description="Animate paired Newton and cellular hillslope rasters.",
    )
    parser.add_argument("--newton", required=True, help="Newton frame directory.")
    parser.add_argument("--cellular", required=True, help="Cellular frame directory.")
    parser.add_argument("--out", required=True, help="Output GIF path.")
    parser.add_argument("--glob", default="*_model_iteration_*.tif",
                        help="Frame glob (default: *_model_iteration_*.tif).")
    parser.add_argument("--stats-csv", help="Optional common-time residual statistics CSV.")
    parser.add_argument("--fps", type=float, default=4.0)
    parser.add_argument("--loop", type=int, default=0)
    parser.add_argument("--keep-frames", action="store_true")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    output = render_hillslope_comparison_animation(
        args.newton, args.cellular, args.out, pattern=args.glob,
        fps=args.fps, loop=args.loop, stats_csv=args.stats_csv,
        keep_frames=args.keep_frames,
    )
    print(f"Wrote animation: {output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
