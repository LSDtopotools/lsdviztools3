"""CLI: lsdviz3-ridgeline-summary — triple divide junctions, profile and map.

The chart: native and topological triple divide junctions at three
thresholds each, the ridgeline series with a maximal-curvature sign barcode
at every polyfit scale, and a key (``<write_path>/{png,svg,pdf}/<write_name>``).
The map: the ridgeline over the hillshade with the channel network (native)
or the paired flowpaths (topological) of the highlighted threshold
(``<write_path>/{png,pdf}/<write_name>_map``; GMT writes no SVG).

Run the analysis on an lsdtt3 oblique frame DEM to lay the map along the
ridge; the map follows whatever CRS the hillshade has.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ..charts.ridgeline import plot_ridgeline_profile_summary, save_figure_formats
from ..io.ridgeline import load_ridgeline_profile, load_tdj, select_profile
from ._ridgeline import (
    add_driver_arg, add_output_args, add_profile_selection_args, distance_range,
    parse_with_driver, require,
)

_TDJ_KEYS = [f"{kind}_tdj_csv_{i}" for kind in ("native", "topological") for i in (1, 2, 3)]


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-ridgeline-summary``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-ridgeline-summary",
        description="Plot triple divide junctions, a ridgeline series with its "
                    "curvature barcode, and a map of the ridgeline.")
    add_driver_arg(parser)
    parser.add_argument("--ridgeline-profile-csv", dest="ridgeline_profile_csv", default=None,
                        help="*_with_polyfit.csv from lsdtt-ridgeline-analysis.")
    parser.add_argument("--hillshade-raster-file", dest="hillshade_raster_file", default=None,
                        help="Hillshade raster for the map (omit for no map).")
    for key in _TDJ_KEYS:
        parser.add_argument("--" + key.replace("_", "-"), dest=key, default=None,
                            help="TDJ CSV; the threshold is read from the file name.")
    parser.add_argument("--channel-network-vector-file", dest="channel_network_vector_file",
                        default=None, help="Channel network vector (native highlight).")
    parser.add_argument("--paired-flowpaths-vector-file", dest="paired_flowpaths_vector_file",
                        default=None, help="Paired flowpaths vector (topological highlight).")
    parser.add_argument("--ridgeline-vector-file", dest="ridgeline_vector_file", default=None,
                        help="Accepted for old drivers and ignored: the map draws the ridge "
                             "from the profile's latitude/longitude.")
    add_profile_selection_args(parser)
    parser.add_argument("--ridgeline-series", dest="ridgeline_series",
                        default="detrended_elevation", help="Profile column to plot.")
    parser.add_argument("--y-axis-label", dest="y_axis_label",
                        default="Detrended\nridge-crest\nelevation [m]",
                        help="Label of the series (\\n for new lines).")
    parser.add_argument("--highlighted-tdj-type", dest="highlighted_tdj_type", default="native",
                        choices=("native", "topological"))
    parser.add_argument("--highlighted-tdj-threshold", dest="highlighted_tdj_threshold",
                        type=float, default=1000.0)
    parser.add_argument("--vector-raster-dpi", dest="vector_raster_dpi", type=float,
                        default=300.0, help="Resolution of rasters inside SVG/PDF.")
    add_output_args(parser, None)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    parser = build_parser()
    args, _ = parse_with_driver(parser, argv)
    args.y_axis_label = args.y_axis_label.replace("\\n", "\n")
    require(parser, args, "ridgeline_profile_csv", "write_path", *_TDJ_KEYS)
    name = args.write_name or f"ridgeline_{args.ridgeline_number}_profile_summary"
    try:
        full = load_ridgeline_profile(args.ridgeline_profile_csv, args.ridgeline_number,
                                      required=("distance", args.ridgeline_series))
        rng = distance_range(args)
        profile = select_profile(full, args.reverse_ridgeline, rng)
        native = [load_tdj(getattr(args, f"native_tdj_csv_{i}"), "native",
                           args.ridgeline_number) for i in (1, 2, 3)]
        topo = [load_tdj(getattr(args, f"topological_tdj_csv_{i}"), "topological",
                         args.ridgeline_number) for i in (1, 2, 3)]
        highlight = (args.highlighted_tdj_type, args.highlighted_tdj_threshold)
        fig, _ = plot_ridgeline_profile_summary(
            profile, native, topo, series=args.ridgeline_series, y_axis_label=args.y_axis_label,
            highlight=highlight, reverse=args.reverse_ridgeline, distance_range=rng,
            profile_distance=full["distance"])
        outputs = save_figure_formats(fig, args.write_path, name,
                                      vector_raster_dpi=args.vector_raster_dpi)
        if args.hillshade_raster_file:
            from ..render.ridgeline import render_ridgeline_map

            chosen = native if highlight[0] == "native" else topo
            tdj = next(t for t in chosen if t["threshold"] == highlight[1])
            overlay = (args.channel_network_vector_file if highlight[0] == "native"
                       else args.paired_flowpaths_vector_file)
            fig_map = render_ridgeline_map(
                args.hillshade_raster_file, profile, highlighted_tdj=tdj, overlay=overlay,
                reverse=args.reverse_ridgeline, distance_range=rng,
                profile_distance=full["distance"])
            for ext in ("png", "pdf"):
                dest = Path(args.write_path) / ext / f"{name}_map.{ext}"
                dest.parent.mkdir(parents=True, exist_ok=True)
                fig_map.savefig(dest, dpi=int(args.vector_raster_dpi))
                outputs.append(dest)
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
