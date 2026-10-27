"""CLI: lsdviz3-floodplain - floodplain polygons over a hillshade, and the HAND QQ.

Reads the outputs of one lsdtt3 ``lsdtt-fluvial-surfaces`` run by its stem
(read_path + write_fname):

  <stem>_floodplain_first_pass_polygons.<fgb|gpkg>  (second pass only)
  <stem>_floodplain_polygons.<fgb|gpkg>
  <stem>_floodplain_thresholds.csv and the QQ tables it names

and writes ``--out-map`` (polygons over ``--base``) and ``--out-qq`` (one QQ
panel per HAND pass). Either output may be omitted.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..render.style import MapStyle


def _find(stem: str, suffix: str):
    for ext in (".fgb", ".gpkg"):
        p = Path(stem + suffix + ext)
        if p.exists():
            return p
    return None


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-floodplain``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-floodplain",
        description="Floodplain polygons (first and second pass) over a "
                    "hillshade, and the HAND QQ plots the bounds came from.",
    )
    parser.add_argument("--stem", required=True,
                        help="lsdtt-fluvial-surfaces output stem (path + write_fname).")
    parser.add_argument("--base", help="Hillshade raster (or DEM with --base-is-dem).")
    parser.add_argument("--base-is-dem", action="store_true",
                        help="Illuminate --base instead of drawing it flat.")
    parser.add_argument("--out-map", help="Output PNG for the polygon map.")
    parser.add_argument("--out-qq", help="Output PNG for the HAND QQ plots.")
    parser.add_argument("--fill-transparency", type=float, default=85.0,
                        help="Percent transparency of the polygon fills (default 85).")
    parser.add_argument("--outline-width", default="0.35p",
                        help="GMT pen width of the outlines (default 0.35p).")
    parser.add_argument("--outline-transparency", type=float, default=40.0,
                        help="Percent transparency of the outlines (default 40).")
    parser.add_argument("--drape",
                        help="Raster coloured under the polygons, e.g. "
                             "<stem>_reference_field.tif (HAND).")
    parser.add_argument("--drape-label", default="HAND (m)",
                        help="Colourbar label for --drape (default 'HAND (m)').")
    parser.add_argument("--drape-vmin", type=float, default=None)
    parser.add_argument("--drape-vmax", type=float, default=None)
    parser.add_argument("--drape-cmap", default="viridis")
    parser.add_argument("--colours", default=None,
                        help="Comma-separated GMT colours, first pass then second "
                             "(default blue, orange).")
    parser.add_argument("--figure-size", default="24c", help="Map width (GMT size).")
    parser.add_argument("--dpi", type=int, default=200)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if not args.out_map and not args.out_qq:
        raise SystemExit("Nothing to do: give --out-map and/or --out-qq.")

    if args.out_map:
        if not args.base:
            raise SystemExit("--out-map needs --base.")
        first = _find(args.stem, "_floodplain_first_pass_polygons")
        final = _find(args.stem, "_floodplain_polygons")
        if final is None:
            raise SystemExit(f"No {args.stem}_floodplain_polygons.fgb/.gpkg "
                             "(run with write_floodplain_polygons: true).")
        sets, labels = ([first, final], ["first pass", "second pass"]) \
            if first is not None else ([final], ["floodplain"])
        from ..render.floodplain import render_floodplain_polygons
        fig = render_floodplain_polygons(
            args.base, sets, labels=labels,
            colours=args.colours.split(",") if args.colours else None,
            base_is_hillshade=not args.base_is_dem,
            fill_transparency=args.fill_transparency,
            outline_width=args.outline_width,
            outline_transparency=args.outline_transparency,
            drape_raster=args.drape,
            drape_cmap=args.drape_cmap,
            drape_vmin=args.drape_vmin, drape_vmax=args.drape_vmax,
            drape_label=args.drape_label if args.drape else None,
            style=MapStyle(figure_size=args.figure_size),
        )
        fig.savefig(args.out_map, dpi=args.dpi)
        print("wrote", args.out_map)

    if args.out_qq:
        from ..charts.floodplain_qq import plot_floodplain_hand_qq
        plot_floodplain_hand_qq(args.stem + "_floodplain_thresholds.csv", args.out_qq)
        print("wrote", args.out_qq)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
