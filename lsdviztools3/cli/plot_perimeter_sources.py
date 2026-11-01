"""CLI: lsdviz3-perimeter-sources - valleys crossing the edge of a DEM.

Reads one lsdtt3 ``lsdtt-channel-extraction`` run made with
``run_perimeter_sources: true`` by its stem (write_path + write_fname) and
writes any of:

  --out-map      the flat perimeter runs coloured by verdict, the outlets and
                 the sources on the hillshade
  --out-profile  elevation and gradient round the perimeter, flat runs by
                 verdict, outlets and sources marked
"""

from __future__ import annotations

import argparse


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-perimeter-sources``."""
    p = argparse.ArgumentParser(
        prog="lsdviz3-perimeter-sources",
        description="Map and profile of the perimeter valley sources from an "
                    "lsdtt-channel-extraction run_perimeter_sources run.")
    p.add_argument("--stem", required=True,
                   help="lsdtt-channel-extraction output stem (path + write_fname).")
    p.add_argument("--base", help="Hillshade raster on the run's grid (--out-map).")
    p.add_argument("--dem", help="DEM to compute the hillshade from when no "
                                 "--base is given (--out-map).")
    p.add_argument("--out-map", help="PNG: flat runs, outlets and sources on the hillshade.")
    p.add_argument("--out-profile", help="PNG: elevation and gradient round the perimeter.")
    p.add_argument("--n-unselected-outlets", type=int, default=3,
                   help="Also draw this many of the largest unselected outlet "
                        "runs (default 3).")
    p.add_argument("--max-gradient", type=float, default=None,
                   help="Threshold line on the profile's gradient panel "
                        "(default the run's perimeter_max_gradient).")
    p.add_argument("--bbox", default=None,
                   help="Map extent xmin,ymin,xmax,ymax in map units "
                        "(default the raster extent).")
    p.add_argument("--title", default=None)
    p.add_argument("--dpi", type=int, default=200)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if not (args.out_map or args.out_profile):
        raise SystemExit("Nothing to do: give --out-map and/or --out-profile.")
    from ..charts.perimeter_sources import (
        plot_perimeter_profile,
        plot_perimeter_sources_map,
    )
    if args.out_map and not (args.base or args.dem):
        raise SystemExit("--out-map needs --base or --dem.")
    bbox = tuple(float(v) for v in args.bbox.split(",")) if args.bbox else None
    if bbox is not None and len(bbox) != 4:
        raise SystemExit("--bbox needs four numbers: xmin,ymin,xmax,ymax.")
    if args.out_map:
        plot_perimeter_sources_map(args.stem, args.out_map, hillshade=args.base,
                                   dem=args.dem, bbox=bbox,
                                   n_unselected_outlets=args.n_unselected_outlets,
                                   title=args.title, dpi=args.dpi)
        print("wrote", args.out_map)
    if args.out_profile:
        plot_perimeter_profile(args.stem, args.out_profile,
                               max_gradient=args.max_gradient,
                               n_unselected_outlets=args.n_unselected_outlets,
                               title=args.title, dpi=args.dpi)
        print("wrote", args.out_profile)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
