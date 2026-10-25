"""CLI: lsdviz3-ridgeline-curvature — ridge-crest detrending and curvature.

Four stacked panels against distance along one ``lsdtt-ridgeline-analysis``
ridgeline: elevation with its fitted trend, moving-window detrended
elevations, maximal curvature and the Laplacian at four polyfit scales.
Writes ``<write_path>/{png,svg,pdf}/<write_name>``.

The profile CSV needs the moving-window columns
(``ridgeline_detrend_moving_average_window_pixels``) and
``maximum_curvature``/``laplacian`` polyfit metrics at each lambda.
"""

from __future__ import annotations

import argparse
import sys

from ..charts.ridgeline import plot_ridgeline_detrending_curvature, save_figure_formats
from ..io.ridgeline import load_ridgeline_profile, select_profile
from ._ridgeline import (
    add_driver_arg, add_output_args, add_profile_selection_args, distance_range,
    number_list, parse_with_driver, require,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-ridgeline-curvature``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-ridgeline-curvature",
        description="Plot ridge-crest elevation, moving-window detrending, maximal "
                    "curvature and the Laplacian along one ridgeline.")
    add_driver_arg(parser)
    parser.add_argument("--ridgeline-profile-csv", dest="ridgeline_profile_csv", default=None,
                        help="*_with_polyfit.csv from lsdtt-ridgeline-analysis.")
    add_profile_selection_args(parser)
    parser.add_argument("--moving-window-sizes-pixels", dest="moving_window_sizes_pixels",
                        type=number_list, default=[25, 49, 99, 149],
                        help="Four moving-window sizes, comma separated (default 25,49,99,149).")
    parser.add_argument("--curvature-lambdas-m", dest="curvature_lambdas_m", type=number_list,
                        default=[25, 49, 99, 149],
                        help="Four polyfit lambdas (m), comma separated (default 25,49,99,149).")
    parser.add_argument("--site-display-name", dest="site_display_name", default=None,
                        help="Optional site name above the top panel.")
    add_output_args(parser, None)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    parser = build_parser()
    args, _ = parse_with_driver(parser, argv)
    require(parser, args, "ridgeline_profile_csv", "write_path")
    windows = args.moving_window_sizes_pixels
    if any(w <= 0 or not float(w).is_integer() for w in windows):
        parser.error("moving_window_sizes_pixels must be positive integers")
    try:
        profile = load_ridgeline_profile(args.ridgeline_profile_csv, args.ridgeline_number,
                                         required=("distance", "elevation", "fitted_elevation"))
        profile = select_profile(profile, args.reverse_ridgeline, distance_range(args))
        fig, _ = plot_ridgeline_detrending_curvature(
            profile, [int(w) for w in windows], args.curvature_lambdas_m,
            site_name=args.site_display_name)
        name = args.write_name or f"ridgeline_{args.ridgeline_number}_profile_detrending_and_curvature"
        outputs = save_figure_formats(fig, args.write_path, name)
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
