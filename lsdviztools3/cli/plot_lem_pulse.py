"""CLI: lsdviz3-lem-pulse — animate a LEM run responding to a pulse of incision.

Takes the output directory of an ``lsdtt-lem`` run that was driven by a
``forcing_history_fname`` and wrote per-frame diagnostics
(``frequency_of_output_writing``, ``write_erosion_rate_raster``,
``write_channel_mask_raster``), and renders each frame as an illuminated
hillshade with the erosion rate draped on top, the timestep's channel network
overlaid, and a mean-erosion-rate history panel underneath marking where the
frame sits in the forcing. The frames are stitched into a GIF.

The erosion-rate colour scale is pinned across all frames, so a colour means the
same rate in every frame and the pulse is visible as a bright band of knickpoints
migrating up the network.
"""

from __future__ import annotations

import argparse

from . import _common
from ..render.lem_pulse import animate_lem_pulse


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-lem-pulse``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-lem-pulse",
        description="Animate an lsdtt-lem incision-pulse run (hillshade + erosion "
        "rate + channel network + forcing history) as a GIF.",
    )
    parser.add_argument(
        "--in",
        dest="frames_dir",
        required=True,
        help="lsdtt-lem output directory holding the _model_iteration_* rasters.",
    )
    parser.add_argument("--out", required=True, help="Output GIF path.")
    parser.add_argument(
        "--prefix",
        default=None,
        help="The run's write_fname. Default: any prefix. Set this when the "
        "directory holds more than one run.",
    )
    parser.add_argument(
        "--ext",
        default=".tif",
        help="Raster extension (default: .tif; use .bil for ENVI output).",
    )
    parser.add_argument(
        "--stats-csv",
        default=None,
        help="The run's _model_elevation_stats.csv (auto-detected in --in by default). "
        "Supplies frame times and the history panel.",
    )
    parser.add_argument(
        "--erosion-cmap",
        default="inferno",
        help="Colour map for the erosion-rate drape (default: inferno).",
    )
    parser.add_argument(
        "--erosion-vmin",
        type=float,
        default=None,
        help="Low end of the shared erosion-rate colour scale (default: 0).",
    )
    parser.add_argument(
        "--erosion-vmax",
        type=float,
        default=None,
        help="High end of the shared erosion-rate colour scale. Default: the "
        "--erosion-percentile of the pooled positive rates over all frames.",
    )
    parser.add_argument(
        "--erosion-percentile",
        type=float,
        default=99.0,
        help="Percentile setting the top of the derived colour scale (default: 99).",
    )
    parser.add_argument(
        "--erosion-transparency",
        type=float,
        default=30.0,
        help="Transparency of the erosion drape, so relief shows through (default: 30).",
    )
    parser.add_argument(
        "--channel-color",
        default="dodgerblue",
        help="Colour of the channel-network overlay (default: dodgerblue).",
    )
    parser.add_argument(
        "--no-timeseries",
        action="store_true",
        help="Omit the mean-erosion-rate history panel under each map.",
    )
    parser.add_argument(
        "--no-annotate",
        action="store_true",
        help="Do not stamp the model time on the frames.",
    )
    parser.add_argument(
        "--fps", type=float, default=5.0, help="Frames per second (default: 5)."
    )
    parser.add_argument(
        "--loop", type=int, default=0, help="GIF loop count (0 = forever)."
    )
    parser.add_argument(
        "--keep-frames",
        action="store_true",
        help="Keep the intermediate PNG frames next to the GIF.",
    )
    _common.add_style_args(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: build the pulse frames and assemble the GIF."""
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)

    # Only pin the colour scale explicitly if the user asked for a limit; otherwise
    # let the renderer derive one that spans every frame.
    erosion_clim = None
    if args.erosion_vmin is not None or args.erosion_vmax is not None:
        erosion_clim = (args.erosion_vmin or 0.0, args.erosion_vmax)

    out = animate_lem_pulse(
        args.frames_dir,
        args.out,
        prefix=args.prefix,
        ext=args.ext,
        stats_csv=args.stats_csv,
        coordinates=args.coordinates,
        erosion_clim=erosion_clim,
        erosion_percentile=args.erosion_percentile,
        erosion_cmap=args.erosion_cmap,
        erosion_transparency=args.erosion_transparency,
        channel_color=args.channel_color,
        style=style,
        show_timeseries=not args.no_timeseries,
        annotate=not args.no_annotate,
        fps=args.fps,
        loop=args.loop,
        keep_frames=args.keep_frames,
    )
    print(f"Wrote animation: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
