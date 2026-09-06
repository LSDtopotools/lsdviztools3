"""CLI: lsdviz3-valley-sections — valley cross sections with the picked floor.

One panel per station along a channel: the ground surface across the section,
the extent picked out as valley floor, the section's low point, and optionally
a field raster (normally HAND) sampled at every point.

Nothing is recomputed; the sections come from lsdtt-valley-centreline's own
sample table.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..charts.valley_sections import plot_valley_sections, summarise_sections


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-valley-sections``.

    Returns:
        The configured parser.
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-valley-sections",
        description=("Plot valley cross sections with the picked valley floor, "
                     "the section low point, and an optional height field."),
    )
    parser.add_argument("--samples", required=True,
                        help="A <stem>_valley_cross_section_samples.csv.")
    parser.add_argument("--out", required=True, help="Output PNG path.")
    parser.add_argument("--source-id", type=int, default=0)
    parser.add_argument("--field", default=None,
                        help="Optional raster sampled at each section point, "
                             "normally the HAND raster.")
    parser.add_argument("--field-label", default="HAND (m)")
    parser.add_argument("--max-panels", type=int, default=30)
    parser.add_argument("--n-cols", type=int, default=5)
    parser.add_argument("--summary", default=None,
                        help="Optional CSV of per-station width and "
                             "centreline-to-low-point offset.")
    parser.add_argument("--title", default=None)
    parser.add_argument("--dpi", type=int, default=150)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Optional argument list.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    png = plot_valley_sections(args.samples, args.out,
                               source_id=args.source_id,
                               field_raster=args.field,
                               field_label=args.field_label,
                               max_panels=args.max_panels,
                               n_cols=args.n_cols,
                               title=args.title, dpi=args.dpi)
    print(f"wrote {png}")
    if args.summary:
        s = summarise_sections(args.samples, source_id=args.source_id)
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        s.to_csv(args.summary, index=False)
        print(f"wrote {args.summary}")
        off = s["low_point_offset_m"].abs()
        print(f"  median |centreline -> low point| = {off.median():.0f} m "
              f"(max {off.max():.0f} m)")
        print(f"  median valley width = {s['valley_width_m'].median():.0f} m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
