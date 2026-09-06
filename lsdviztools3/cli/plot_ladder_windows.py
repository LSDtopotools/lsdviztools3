"""CLI: lsdviz3-ladder-windows — how lsdtt-fluvial-surfaces picked its levels.

Reads a `<stem>_ladder_windows.csv` and plots, for each channel, where the
breaks between fluvial levels sat along the channel, where the levels themselves
sat, and which windows the DEM could not resolve.

Nothing is recomputed: every number comes from lsdtt3's own ladder output.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..charts.ladder_windows import plot_ladder_windows, summarise_ladder


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-ladder-windows``.

    Returns:
        The configured parser.
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-ladder-windows",
        description=(
            "Plot the reference-field ladder's level selection along each "
            "channel: the breaks between levels, the level centres, and the "
            "windows the DEM could not resolve."
        ),
    )
    parser.add_argument("--windows", required=True,
                        help="A <stem>_ladder_windows.csv from "
                             "lsdtt-fluvial-surfaces.")
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument("--summary", default=None,
                        help="Optional CSV path for the per-channel summary.")
    parser.add_argument("--title", default=None, help="Figure title.")
    parser.add_argument("--dpi", type=int, default=160)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Optional argument list (defaults to sys.argv).

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    png = plot_ladder_windows(args.windows, args.out, title=args.title,
                              dpi=args.dpi)
    print(f"wrote {png}")

    summary = summarise_ladder(args.windows)
    if args.summary:
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.summary, index=False)
        print(f"wrote {args.summary}")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
