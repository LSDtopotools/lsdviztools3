"""CLI: lsdviz3-ladder-density — the HAND densities the ladder cut its breaks from.

Plots one panel per along-channel window: the reference-field distribution, the
breaks drawn on it, and how cleanly each break separates the populations either
side of it.

This is the diagnostic that says whether a break is real. The break value alone
cannot, and neither can the resolvability floor - see the module docstring in
lsdviztools3.charts.ladder_density.
"""

from __future__ import annotations

import argparse

from ..charts.ladder_density import plot_ladder_density


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-ladder-density``.

    Returns:
        The configured parser.
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-ladder-density",
        description=("Plot each window's reference-field density with the "
                     "ladder's breaks and their separation ratios."),
    )
    parser.add_argument("--histograms", required=True,
                        help="A <stem>_ladder_histograms.csv "
                             "(write_ladder_histograms_csv: true).")
    parser.add_argument("--windows", required=True,
                        help="The matching <stem>_ladder_windows.csv.")
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument("--source-id", type=int, default=None,
                        help="Restrict to one channel.")
    parser.add_argument("--max-panels", type=int, default=12,
                        help="Panels per channel, sampled evenly downstream.")
    parser.add_argument("--title", default=None)
    parser.add_argument("--dpi", type=int, default=160)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Optional argument list.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    png = plot_ladder_density(args.histograms, args.windows, args.out,
                              source_id=args.source_id,
                              max_panels=args.max_panels,
                              title=args.title, dpi=args.dpi)
    print(f"wrote {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
