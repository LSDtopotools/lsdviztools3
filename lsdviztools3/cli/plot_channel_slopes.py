"""CLI: lsdviz3-channel-slopes - smoothed channel profiles and slopes.

Reads <write_fname>_channel_slope_points.<ext> from lsdtt-channel-extraction
(write_channel_slope_points: every channel cell with its pit-filled and
smoothed elevation and its along-channel slope over the
slope_window_half_length_metres window) and, optionally, field points snapped
to the channel by the same tool (run_channel_point_attributes). Writes:

  <out-prefix>_channel_profiles.png   elevation against flow distance, raw and smoothed
  <out-prefix>_slope_distance.png     smoothed slope against flow distance
  <out-prefix>_slope_area.png         smoothed slope against drainage area
  <out-prefix>_channel_slope_map.png  with --dem: channel cells by log10 slope
"""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-channel-slopes``."""
    p = argparse.ArgumentParser(
        prog="lsdviz3-channel-slopes",
        description="Smoothed channel profiles and along-channel slopes from lsdtt-channel-extraction "
                    "(write_channel_slope_points), with optional field points.")
    p.add_argument("--slope-points", required=True,
                   help="<write_fname>_channel_slope_points.<csv|fgb|gpkg> from lsdtt-channel-extraction.")
    p.add_argument("--field-points", default=None,
                   help="Points on the channel (csv with latitude/longitude, fgb or gpkg) carrying drainage_area_m2, "
                        "slope_m_per_m and, for the profile plots, flow_distance_m and smoothed_elevation_m or "
                        "elevation_m.")
    p.add_argument("--field-color-by", default=None,
                   help="Colour the field points by this attribute (e.g. a grain size column).")
    p.add_argument("--min-area-m2", type=float, default=None,
                   help="Plot only channel cells draining at least this area (m^2).")
    p.add_argument("--dem", default=None, help="DEM for the map's hillshade. Without it no map is drawn.")
    p.add_argument("--slope-vmin", type=float, default=None, help="Map colour range lower end (m/m).")
    p.add_argument("--slope-vmax", type=float, default=None, help="Map colour range upper end (m/m).")
    p.add_argument("--cmap", default="viridis", help="Colour map (default viridis).")
    p.add_argument("--coordinates", default="local", help="Map coordinates: local/projected or geographic.")
    p.add_argument("--out-dir", required=True, help="Directory for the PNGs.")
    p.add_argument("--out-prefix", default=None,
                   help="Prefix of the PNG names (default: the slope-points file stem).")
    p.add_argument("--dpi", type=int, default=300)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ..charts.channel_slopes import (load_channel_points, plot_channel_profiles, plot_slope_area,
                                         plot_slope_distance, render_channel_slope_map)

    try:
        points = load_channel_points(args.slope_points)
        field = load_channel_points(args.field_points) if args.field_points else None
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc))
    if args.field_color_by and field is None:
        raise SystemExit("--field-color-by needs --field-points.")
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = args.out_prefix or Path(args.slope_points).stem.replace("_channel_slope_points", "")
    written = []
    for name, fn in (("channel_profiles", plot_channel_profiles), ("slope_distance", plot_slope_distance),
                     ("slope_area", plot_slope_area)):
        try:
            fig, _ = fn(points, field=field, field_color_by=args.field_color_by,
                        min_area_m2=args.min_area_m2, cmap=args.cmap)
        except ValueError as exc:
            raise SystemExit(str(exc))
        path = out / f"{stem}_{name}.png"
        fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
    if args.dem:
        written.append(render_channel_slope_map(
            args.dem, points, out / f"{stem}_channel_slope_map.png", field=field,
            coordinates=args.coordinates, cmap=args.cmap, slope_vmin=args.slope_vmin,
            slope_vmax=args.slope_vmax, min_area_m2=args.min_area_m2))
    for w in written:
        print(w)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
