"""CLI: lsdviz3-valley-widths - valley widths for every channel.

Reads one lsdtt3 ``lsdtt-valley-centreline`` run made with
``run_valley_width_network: true`` by its stem (write_path + write_fname) and
writes any of:

  --out-map      section width lines on the hillshade, coloured by
                 excluded_reason, with the segment centrelines
  --out-profile  valley width down one segment (default the trunk with the
                 largest drainage area), tributary confluences marked,
                 drainage area below
  --out-channels each segment's extracted D8 channel and its valley
                 centreline on the hillshade
  --out-centreline-profile
                 elevation of the DEM under one segment's centreline against
                 its baseline, and their difference
"""

from __future__ import annotations

import argparse


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-valley-widths``."""
    p = argparse.ArgumentParser(
        prog="lsdviz3-valley-widths",
        description="Map and profile of valley widths from an lsdtt-valley-"
                    "centreline run_valley_width_network run.")
    p.add_argument("--stem", required=True,
                   help="lsdtt-valley-centreline output stem (path + write_fname).")
    p.add_argument("--base", help="Hillshade raster on the run's grid (--out-map).")
    p.add_argument("--valley-mask", help="Valley mask raster to draw under the "
                                         "sections (the run's valley_mask_fname).")
    p.add_argument("--out-map", help="PNG: section width lines on the hillshade.")
    p.add_argument("--out-profile", help="PNG: valley width down one segment.")
    p.add_argument("--out-channels",
                   help="PNG: extracted channels and valley centrelines on the hillshade.")
    p.add_argument("--out-centreline-profile",
                   help="PNG: DEM elevation under one segment's centreline "
                        "against its baseline.")
    p.add_argument("--segment-id", type=int, default=None,
                   help="Segment for --out-profile and --out-centreline-profile "
                        "(default the trunk with the largest drainage area).")
    p.add_argument("--map-segments", default=None,
                   help="Comma-separated segment ids to draw on --out-map and "
                        "--out-channels (default all).")
    p.add_argument("--bbox", default=None,
                   help="Map extent xmin,ymin,xmax,ymax in map units "
                        "(default the sections' extent).")
    p.add_argument("--hide-excluded", action="store_true",
                   help="Draw only the kept sections / stations.")
    p.add_argument("--title", default=None)
    p.add_argument("--dpi", type=int, default=200)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if not (args.out_map or args.out_profile or args.out_channels or
            args.out_centreline_profile):
        raise SystemExit("Nothing to do: give --out-map, --out-profile, "
                         "--out-channels and/or --out-centreline-profile.")
    from ..charts.valley_width_network import (
        plot_centreline_profile,
        plot_channel_and_centreline,
        plot_valley_width_map,
        plot_valley_widths_along_channel,
    )
    if (args.out_map or args.out_channels) and not args.base:
        raise SystemExit("--out-map and --out-channels need --base.")
    bbox = tuple(float(v) for v in args.bbox.split(",")) if args.bbox else None
    if bbox is not None and len(bbox) != 4:
        raise SystemExit("--bbox needs four numbers: xmin,ymin,xmax,ymax.")
    segs = ([int(s) for s in args.map_segments.split(",")]
            if args.map_segments else None)
    if args.out_centreline_profile:
        plot_centreline_profile(args.stem, args.out_centreline_profile,
                                segment_id=args.segment_id, title=args.title,
                                dpi=args.dpi)
        print("wrote", args.out_centreline_profile)
    if args.out_channels:
        plot_channel_and_centreline(args.stem, args.base, args.out_channels,
                                    valley_mask=args.valley_mask, bbox=bbox,
                                    segment_ids=segs, title=args.title, dpi=args.dpi)
        print("wrote", args.out_channels)
    if args.out_map:
        plot_valley_width_map(args.stem, args.base, args.out_map,
                              valley_mask=args.valley_mask, bbox=bbox,
                              segment_ids=segs,
                              show_excluded=not args.hide_excluded,
                              title=args.title, dpi=args.dpi)
        print("wrote", args.out_map)
    if args.out_profile:
        plot_valley_widths_along_channel(args.stem, args.out_profile,
                                         segment_id=args.segment_id,
                                         show_excluded=not args.hide_excluded,
                                         title=args.title, dpi=args.dpi)
        print("wrote", args.out_profile)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
