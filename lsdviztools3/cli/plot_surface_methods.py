"""CLI: lsdviz3-surface-methods - compare the methods of one fluvial-surfaces run.

Reads the outputs of one lsdtt3 ``lsdtt-fluvial-surfaces`` run by its stem
(read_path + write_fname), single- or several-method, and writes any of:

  --out-classes     each method's floodplain / terrace, side by side
  --out-agreement   which methods agree, one map per class
                    (<out>_floodplain.png, <out>_terrace.png)
  --out-levels      patch median HAND per method, with the height clusters
  --out-scores      precision / recall / IoU against the mapped surfaces

``--mapped`` outlines mapped polygons on the maps.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..render.style import MapStyle


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-surface-methods``."""
    p = argparse.ArgumentParser(
        prog="lsdviz3-surface-methods",
        description="Compare the floodplain / terrace methods of one "
                    "lsdtt-fluvial-surfaces run.")
    p.add_argument("--stem", required=True,
                   help="lsdtt-fluvial-surfaces output stem (path + write_fname).")
    p.add_argument("--base", help="Hillshade raster on the run's grid (maps only).")
    p.add_argument("--mapped", help="Mapped-surface polygons to outline on the maps.")
    p.add_argument("--methods", help="Comma-separated subset of methods (default all).")
    p.add_argument("--out-classes", help="PNG: classes per method, side by side.")
    p.add_argument("--out-agreement",
                   help="PNG stem: writes <stem>_floodplain.png and <stem>_terrace.png.")
    p.add_argument("--out-levels", help="PNG: patch median HAND per method.")
    p.add_argument("--levels-log-y", action="store_true",
                   help="Log HAND axis (from 0.1 m) on --out-levels.")
    p.add_argument("--out-scores", help="PNG: scores against the mapped surfaces.")
    p.add_argument("--panel-width", default="10c", help="GMT width of each map panel.")
    p.add_argument("--dpi", type=int, default=200)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if not any([args.out_classes, args.out_agreement, args.out_levels, args.out_scores]):
        raise SystemExit("Nothing to do: give --out-classes, --out-agreement, "
                         "--out-levels and/or --out-scores.")
    methods = args.methods.split(",") if args.methods else None
    if (args.out_classes or args.out_agreement) and not args.base:
        raise SystemExit("The maps need --base (a hillshade).")

    if args.out_classes:
        from ..render.surface_methods import render_method_classes
        fig = render_method_classes(args.stem, args.base, methods=methods,
                                    mapped=args.mapped, panel_width=args.panel_width)
        fig.savefig(args.out_classes, dpi=args.dpi)
        print("wrote", args.out_classes)
    if args.out_agreement:
        from ..render.surface_methods import render_method_agreement
        stem = str(Path(args.out_agreement).with_suffix(""))
        for cls in ("floodplain", "terrace"):
            out = f"{stem}_{cls}.png"
            fig = render_method_agreement(args.stem, args.base, cls, mapped=args.mapped,
                                          style=MapStyle(figure_size=args.panel_width))
            fig.savefig(out, dpi=args.dpi)
            print("wrote", out)
    if args.out_levels:
        from ..charts.surface_methods import plot_method_levels
        plot_method_levels(args.stem, args.out_levels, methods=methods,
                           log_y=args.levels_log_y)
        print("wrote", args.out_levels)
    if args.out_scores:
        from ..charts.surface_methods import plot_mapped_scores
        plot_mapped_scores(args.stem, args.out_scores)
        print("wrote", args.out_scores)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
