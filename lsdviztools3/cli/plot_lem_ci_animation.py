"""CLI: lsdviz3-lem-ci — animate an LEM run whose channelization index is walked.

Takes the output directory of an ``lsdtt-lem`` run driven by a forcing history
that steps ``erodibility`` and/or ``hillslope_diffusivity`` (and so the
channelization index), with ``frequency_of_output_writing`` set. Each frame has a
narrow CI time-series panel on the left -- time running down the vertical axis,
the full run drawn, a dot at the current time -- and the landscape on the right,
as a shaded-relief map (``--view map``) or a 3D oblique view (``--view oblique``).
The frames are stitched into a GIF; its duration is n_frames / fps.
"""

from __future__ import annotations

import argparse

from . import _common
from ..render.lem_ci_animation import animate_ci_run


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-lem-ci``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-lem-ci",
        description="Animate an lsdtt-lem channelization-index ladder run "
        "(CI history panel + map or oblique landscape) as a GIF.",
    )
    parser.add_argument("--in", dest="frames_dir", required=True,
                        help="lsdtt-lem output directory holding the _model_iteration_* rasters.")
    parser.add_argument("--out", required=True, help="Output GIF path.")
    parser.add_argument("--view", choices=["map", "oblique"], default="map",
                        help="Landscape panel: shaded-relief map or 3D oblique (default: map).")
    parser.add_argument("--prefix", default=None,
                        help="The run's write_fname (default: any prefix).")
    parser.add_argument("--ext", default=".tif", help="Raster extension (default: .tif).")
    parser.add_argument("--stats-csv", default=None,
                        help="The run's _model_elevation_stats.csv (auto-detected by default).")
    parser.add_argument("--no-initial", action="store_true",
                        help="Do not start from <prefix>_initial_topography.")
    parser.add_argument("--map-width", default="15c",
                        help="Width of the landscape panel (default: 15c).")
    parser.add_argument("--azimuth", type=float, default=200.0,
                        help="Oblique camera bearing in degrees (default: 200).")
    parser.add_argument("--elevation", type=float, default=35.0,
                        help="Oblique camera elevation angle in degrees (default: 35).")
    parser.add_argument("--vertical-exaggeration", type=float, default=2.0,
                        help="Oblique vertical exaggeration (default: 2).")
    parser.add_argument("--landscape-cmap", default="gray",
                        help="Colour map of the landscape panel (default: gray). The "
                        "shared --cmap flag is not used here: its default is 'geo'.")
    parser.add_argument("--fps", type=float, default=10.0,
                        help="Frames per second (default: 10).")
    parser.add_argument("--frame-dpi", type=int, default=110,
                        help="Frame resolution; GIF size grows quickly with it (default: 110).")
    parser.add_argument("--loop", type=int, default=0, help="GIF loop count (0 = forever).")
    parser.add_argument("--keep-frames", action="store_true",
                        help="Keep the intermediate PNG frames next to the GIF.")
    _common.add_style_args(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: render the frames and assemble the GIF."""
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)
    out = animate_ci_run(
        args.frames_dir,
        args.out,
        view=args.view,
        prefix=args.prefix,
        ext=args.ext,
        stats_csv=args.stats_csv,
        include_initial=not args.no_initial,
        coordinates=args.coordinates,
        map_width=args.map_width,
        oblique_azimuth=args.azimuth,
        oblique_elevation=args.elevation,
        vertical_exaggeration=args.vertical_exaggeration,
        cmap=args.landscape_cmap,
        fps=args.fps,
        dpi=args.frame_dpi,
        loop=args.loop,
        keep_frames=args.keep_frames,
        style=style,
    )
    print(f"Wrote animation: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
