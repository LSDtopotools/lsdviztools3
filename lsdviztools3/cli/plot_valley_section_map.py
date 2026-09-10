"""CLI: lsdviz3-valley-section-map — cross sections against centreline and channel.

Maps the numbered cross sections on the hillshade with the valley centreline and
the channel network, and marks any section that crosses no channel cell.
"""

from __future__ import annotations

import argparse

from ..charts.valley_section_map import plot_valley_section_map


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-valley-section-map``.

    Returns:
        The configured parser.
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-valley-section-map",
        description=("Map valley cross sections with the centreline and the "
                     "channel network."))
    parser.add_argument("--hillshade", required=True)
    parser.add_argument("--sections", required=True,
                        help="<stem>_valley_cross_sections.fgb")
    parser.add_argument("--centreline", default=None,
                        help="<stem>_valley_centreline.fgb")
    parser.add_argument("--channel", default=None,
                        help="<stem>_channel_network_lines.fgb")
    parser.add_argument("--samples", default=None,
                        help="<stem>_valley_cross_section_samples.csv, used to "
                             "flag sections crossing no channel cell.")
    parser.add_argument("--sides", default=None,
                        help="<stem>_valley_cross_section_sides.fgb, the "
                             "detected valley margins.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--source-id", type=int, default=0)
    parser.add_argument("--label-every", type=int, default=2)
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
    png = plot_valley_section_map(
        args.hillshade, args.sections, args.out,
        centreline_fgb=args.centreline, channel_fgb=args.channel,
        samples_csv=args.samples, sides_fgb=args.sides,
        source_id=args.source_id,
        label_every=args.label_every, title=args.title, dpi=args.dpi)
    print(f"wrote {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
