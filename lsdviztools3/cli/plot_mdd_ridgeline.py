"""CLI: lsdviz3-mdd-ridgeline — main drainage divide profile over its model.

The detrended main-drainage-divide profile, with a maximal-curvature sign
barcode and the median zero-crossing spacing, above the model's drainage-area
image (an ``lsdtt-lem`` log10 drainage-area PNG), with guides every 200 m
through both. For an east-west divide across a periodic model domain.
Writes ``<write_path>/{png,svg,pdf}/<write_name>``.
"""

from __future__ import annotations

import argparse
import sys

from ..charts.ridgeline import plot_mdd_ridgeline_map_profile, save_figure_formats
from ..io.ridgeline import divide_grid_spacing, load_ridgeline_profile
from ._ridgeline import add_driver_arg, add_output_args, parse_with_driver, require


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-mdd-ridgeline``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-mdd-ridgeline",
        description="Plot a main-drainage-divide ridgeline profile above the model's "
                    "drainage-area image.")
    add_driver_arg(parser)
    parser.add_argument("--ridgeline-profile-csv", dest="ridgeline_profile_csv", default=None,
                        help="Range-defined *_with_polyfit.csv with a maximum_curvature metric.")
    parser.add_argument("--drainage-area-map-png", dest="drainage_area_map_png", default=None,
                        help="lsdtt-lem log10 drainage-area PNG of the same model grid.")
    parser.add_argument("--model-name", dest="model_name", default=None,
                        help="Model label at the top left.")
    parser.add_argument("--ridgeline-number", dest="ridgeline_number", type=int, default=0,
                        help="Ridgeline to plot (default 0, the main divide).")
    parser.add_argument("--guide-interval", dest="guide_interval", type=float, default=200.0,
                        help="Distance between guides (m, default 200).")
    add_output_args(parser, "mdd_ridgeline_map_profile")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    parser = build_parser()
    args, _ = parse_with_driver(parser, argv)
    require(parser, args, "ridgeline_profile_csv", "drainage_area_map_png", "model_name",
            "write_path")
    try:
        import numpy as np
        from PIL import Image

        profile = load_ridgeline_profile(
            args.ridgeline_profile_csv, args.ridgeline_number,
            required=("row", "col", "distance", "detrended_elevation", "fitted_elevation"))
        with Image.open(args.drainage_area_map_png) as png:
            image = np.asarray(png.convert("RGB"))
        fig, _ = plot_mdd_ridgeline_map_profile(
            profile, image, divide_grid_spacing(profile), str(args.model_name),
            guide_interval=args.guide_interval)
        outputs = save_figure_formats(fig, args.write_path, args.write_name)
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
