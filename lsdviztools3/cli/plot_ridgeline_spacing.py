"""CLI: lsdviz3-ridgeline-spacing — curvature zero-crossing spacing and slope-area.

One row per site: the spread of the spacing between maximal-curvature sign
changes against polyfit scale, beside the slope-area distribution with the
characteristic length L and 2L. Writes ``<write_path>/{png,svg,pdf}/<write_name>``
and ``<write_path>/csv/spacing_summary_statistics.csv``.

Sites are given with repeated ``--site NAME PROFILE_CSV SLOPE_AREA_CSV`` or,
in a driver file, as ``site_count`` and numbered ``site_<N>_*`` keys
(``display_name``, ``ridgeline_profile_csv``, ``slope_area_csv``,
``ridgeline_number``, ``reverse_ridgeline``, ``use_distance_range``,
``minimum_distance``, ``maximum_distance``, ``absolute_curvature``).

The spacing statistics and L are computed here for now; they are due to move
into lsdtt3.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from ..charts.ridgeline import plot_zero_spacing_slope_area, save_figure_formats
from ..io.ridgeline import (
    curvature_spacing_quantiles, load_ridgeline_profile, load_slope_area, spacing_summary,
)
from ._ridgeline import add_driver_arg, add_output_args, boolean, parse_with_driver, require

SUMMARY_COLUMNS = (
    "site", "A_max_median_slope", "dx", "L", "2xL", "closest_lambda_to_2xL",
    "median_zero_crossing_spacing_at_closest_lambda", "median_zero_spacing_divided_by_L",
)


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-ridgeline-spacing``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-ridgeline-spacing",
        description="Plot maximal-curvature zero-crossing spacing against polyfit scale "
                    "beside slope-area, per site.")
    add_driver_arg(parser)
    parser.add_argument("--site", action="append", nargs=3, default=[],
                        metavar=("NAME", "PROFILE_CSV", "SLOPE_AREA_CSV"),
                        help="A site (ridgeline 0, whole profile); repeat per site.")
    parser.add_argument("--site-count", dest="site_count", type=int, default=None,
                        help="Number of site_<N>_* driver entries.")
    parser.add_argument("--log-spacing", dest="log_spacing", type=boolean, default=False)
    parser.add_argument("--fixed-spacing-ylim", dest="fixed_spacing_ylim", type=boolean,
                        default=True)
    parser.add_argument("--spacing-y-max-m", dest="spacing_y_max_m", type=float, default=200.0)
    parser.add_argument("--absolute-curvature-label", dest="absolute_curvature_label",
                        type=boolean, default=False)
    parser.add_argument("--figure-width-inches", dest="figure_width_inches", type=float,
                        default=2.75)
    parser.add_argument("--row-height-inches", dest="row_height_inches", type=float,
                        default=1.1)
    parser.add_argument("--xlabel-space-inches", dest="xlabel_space_inches", type=float,
                        default=0.45)
    parser.add_argument("--font-size", dest="font_size", type=float, default=6.0)
    add_output_args(parser, "zero_spacing_slope_area")
    return parser


def _driver_sites(parser, count, extra) -> list:
    sites = []
    for i in range(1, count + 1):
        prefix = f"site_{i}_"
        site = {}
        for key in ("display_name", "ridgeline_profile_csv", "slope_area_csv"):
            value = extra.get(prefix + key)
            if not isinstance(value, str) or not value.strip():
                parser.error(f"driver is missing {prefix + key!r}")
            site[key] = value.strip()
        number = extra.get(prefix + "ridgeline_number", 0)
        if isinstance(number, bool) or not isinstance(number, int):
            parser.error(f"{prefix}ridgeline_number must be an integer")
        site["ridgeline_number"] = number
        for key in ("reverse_ridgeline", "use_distance_range", "absolute_curvature"):
            value = extra.get(prefix + key, False)
            if not isinstance(value, bool):
                parser.error(f"{prefix}{key} must be true or false")
            site[key] = value
        lo = float(extra.get(prefix + "minimum_distance", 0.0))
        hi = float(extra.get(prefix + "maximum_distance", 0.0))
        if site["use_distance_range"] and lo > hi:
            parser.error(f"{prefix}minimum_distance must not exceed maximum_distance")
        site["distance_range"] = (lo, hi) if site["use_distance_range"] else None
        sites.append(site)
    unused = sorted(k for k in extra
                    if int(k.split("_")[1]) > count)
    if unused:
        parser.error(f"site entries beyond site_count: {', '.join(unused)}")
    return sites


def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    parser = build_parser()
    args, extra = parse_with_driver(parser, argv, extra_key_patterns=(r"site_\d+_\w+",))
    require(parser, args, "write_path")
    sites = [{"display_name": n, "ridgeline_profile_csv": p, "slope_area_csv": s,
              "ridgeline_number": 0, "reverse_ridgeline": False, "distance_range": None,
              "absolute_curvature": False} for n, p, s in args.site]
    if args.site_count is not None:
        if args.site_count < 1:
            parser.error("site_count must be a positive integer")
        sites += _driver_sites(parser, args.site_count, extra)
    elif extra:
        parser.error("site_<N>_* driver entries need site_count")
    if not sites:
        parser.error("give at least one --site, or site_count with site_<N>_* driver entries")
    try:
        plotted = []
        for site in sites:
            profile = load_ridgeline_profile(site["ridgeline_profile_csv"],
                                             site["ridgeline_number"])
            spacing = curvature_spacing_quantiles(profile, site["reverse_ridgeline"],
                                                  site["distance_range"],
                                                  site["absolute_curvature"])
            plotted.append({"name": site["display_name"], "spacing": spacing,
                            "slope_area": load_slope_area(site["slope_area_csv"])})
        fig, _ = plot_zero_spacing_slope_area(
            plotted, log_spacing=args.log_spacing, fixed_spacing_ylim=args.fixed_spacing_ylim,
            spacing_y_max=args.spacing_y_max_m,
            absolute_curvature_label=args.absolute_curvature_label,
            figure_width=args.figure_width_inches, row_height=args.row_height_inches,
            xlabel_space=args.xlabel_space_inches, font_size=args.font_size)
        summary = Path(args.write_path) / "csv" / "spacing_summary_statistics.csv"
        summary.parent.mkdir(parents=True, exist_ok=True)
        with summary.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=SUMMARY_COLUMNS)
            writer.writeheader()
            writer.writerows(spacing_summary(s["name"], s["spacing"], s["slope_area"])
                             for s in plotted)
        outputs = save_figure_formats(fig, args.write_path, args.write_name) + [summary]
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
