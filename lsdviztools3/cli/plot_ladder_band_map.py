"""CLI: lsdviz3-ladder-bands — map the ladder band index over a hillshade."""

from __future__ import annotations

import argparse

from ..charts.ladder_band_map import plot_ladder_band_map


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-ladder-bands``.

    Returns:
        The configured parser.
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-ladder-bands",
        description="Map the ladder band index over a hillshade.")
    parser.add_argument("--hillshade", required=True)
    parser.add_argument("--bands", required=True,
                        help="<stem>_band_index.tif")
    parser.add_argument("--channel", default=None)
    parser.add_argument("--centreline", default=None)
    parser.add_argument("--max-band", type=int, default=None)
    parser.add_argument("--out", required=True)
    parser.add_argument("--title", default=None)
    parser.add_argument("--dpi", type=int, default=200)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Optional argument list.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    png = plot_ladder_band_map(
        args.hillshade, args.bands, args.out, channel_fgb=args.channel,
        centreline_fgb=args.centreline, max_band=args.max_band,
        title=args.title, dpi=args.dpi)
    print(f"wrote {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
