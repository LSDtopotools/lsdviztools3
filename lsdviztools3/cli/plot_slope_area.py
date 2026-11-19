"""CLI: lsdviz3-slope-area - the outputs of an lsdtt-slope-area run.

Reads the run's binned tables (``<prefix>_lambda_<L>m.csv`` for polyfit,
``<prefix>_d8window_<W>m.csv`` for the D8 flow-path gradient) and, when the
run had ``run_low_concavity_detection: true``, its low-concavity outputs,
and writes:

  <out-prefix>_binned.png  binned d50 gradient (d25-d75 band) against drainage
      area, one line per polyfit lambda / D8 window
  <out-prefix>_method_comparison.png  with --compare-data-dir and/or
      --compare-prefix: the binned d50 of this run and the comparison run
  <out-prefix>_theta_distribution.png  with low-concavity outputs: histogram
      of theta for channel pixels with a valid fit, with the +/- threshold
      (low_concavity_max_abs_concavity of the run)
  <out-prefix>_low_concavity_slope_area.png  with the points CSV
      (write_low_concavity_point_table: true): the channel-pixel slope-area
      cloud, low-concavity pixels in a second colour
  <out-prefix>_low_concavity_map.png  with --dem and the segments file: the
      segments on a hillshade, coloured by low_concavity

Everything drawn is read from the files; nothing is recomputed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

#: Used when the run's ingested parameters cannot be read (lsdtt3's default).
DEFAULT_MAX_ABS_CONCAVITY = 0.15


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-slope-area``."""
    p = argparse.ArgumentParser(
        prog="lsdviz3-slope-area",
        description="Plots of an lsdtt-slope-area run: the binned slope-area tables, a comparison "
                    "with a second run, and the low-concavity outputs (theta distribution, "
                    "channel-pixel slope-area cloud, map).")
    p.add_argument("--data-dir", required=True, help="Directory holding the lsdtt-slope-area outputs.")
    p.add_argument("--prefix", required=True, help="The run's write_fname.")
    p.add_argument("--out-dir", required=True, help="Directory for the outputs.")
    p.add_argument("--out-prefix", default=None, help="Prefix of the output names (default --prefix).")
    p.add_argument("--dem", default=None,
                   help="DEM for the low-concavity map's hillshade. Without it no map is drawn.")
    p.add_argument("--compare-data-dir", default=None,
                   help="Directory of a second run for the method comparison (default --data-dir "
                        "when only --compare-prefix is given).")
    p.add_argument("--compare-prefix", default=None,
                   help="write_fname of the second run (default --prefix when only "
                        "--compare-data-dir is given).")
    p.add_argument("--name", default=None,
                   help="Name of this run in the comparison legend (default: the --data-dir name).")
    p.add_argument("--compare-name", default=None,
                   help="Name of the comparison run in the legend (default: the --compare-data-dir "
                        "name).")
    p.add_argument("--comparison-bands", action="store_true",
                   help="Draw the d25-d75 bands on the method comparison too (off by default).")
    p.add_argument("--theta-range", type=float, nargs=2, default=None, metavar=("LOW", "HIGH"),
                   help="x range of the theta histogram (default: 2nd-98th percentile of the valid "
                        "theta, widened to 2 x the threshold either side).")
    p.add_argument("--theta-bins", type=int, default=60, help="Bins of the theta histogram (default 60).")
    p.add_argument("--no-width-by-order", action="store_true",
                   help="Draw every map segment at one width instead of by stream order.")
    p.add_argument("--coordinates", default="local", help="Map coordinates: local/projected or geographic.")
    p.add_argument("--dpi", type=int, default=300, help="DPI of the charts (default 300).")
    return p


def _save(fig, path, dpi, written):
    import matplotlib.pyplot as plt

    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    written.append(path)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    data_dir, out_dir = Path(args.data_dir), Path(args.out_dir)
    stem = args.out_prefix or args.prefix
    if not data_dir.is_dir():
        raise SystemExit(f"{data_dir} is not a directory.")

    from ..io.slope_area import (find_binned_slope_area, load_ingested_params,
                                 load_low_concavity_points, load_low_concavity_rasters,
                                 low_concavity_files, param_float)

    tables = find_binned_slope_area(data_dir, args.prefix)
    lc = low_concavity_files(data_dir, args.prefix)
    has_lc = any(lc[k].exists() for k in ("points", "segments", "flag"))
    if not tables and not has_lc:
        raise SystemExit(f"No lsdtt-slope-area outputs for prefix {args.prefix!r} in {data_dir} "
                         f"(looked for {args.prefix}_lambda_<L>m.csv, {args.prefix}_d8window_<W>m.csv "
                         f"and {args.prefix}_low_concavity_*).")
    out_dir.mkdir(parents=True, exist_ok=True)
    for t in tables:
        print(f"{t['label']}: {len(t['table'])} bins ({t['path'].name})")

    from ..charts.slope_area import (plot_binned_slope_area, plot_low_concavity_slope_area,
                                     plot_slope_area_method_comparison, plot_theta_distribution,
                                     render_low_concavity_map)

    written = []
    if tables:
        fig, _ = plot_binned_slope_area(tables, title=f"{args.prefix}: binned slope-area")
        _save(fig, out_dir / f"{stem}_binned.png", args.dpi, written)
    else:
        print("no binned slope-area tables: no _binned.png")

    if args.compare_data_dir or args.compare_prefix:
        cdir = Path(args.compare_data_dir) if args.compare_data_dir else data_dir
        cprefix = args.compare_prefix or args.prefix
        if not cdir.is_dir():
            raise SystemExit(f"{cdir} is not a directory.")
        ctables = find_binned_slope_area(cdir, cprefix)
        if not ctables:
            raise SystemExit(f"No binned slope-area tables for prefix {cprefix!r} in {cdir}.")
        if not tables:
            raise SystemExit("The method comparison needs binned tables in --data-dir too.")
        name = args.name or data_dir.resolve().name
        cname = args.compare_name or cdir.resolve().name
        if name == cname:
            name, cname = f"{name} ({args.prefix})", f"{cname} ({cprefix})"
        fig, _ = plot_slope_area_method_comparison(tables, ctables, name=name, compare_name=cname,
                                                   title="binned slope-area: method comparison",
                                                   bands=args.comparison_bands)
        _save(fig, out_dir / f"{stem}_method_comparison.png", args.dpi, written)

    if not has_lc:
        print(f"no low-concavity outputs ({args.prefix}_low_concavity_*): no theta, low-concavity "
              "slope-area or map plots (run_low_concavity_detection: true in the run)")
    else:
        params_path = data_dir / f"{args.prefix}_ingested_params.param"
        threshold = None
        if params_path.exists():
            threshold = param_float(load_ingested_params(params_path), "low_concavity_max_abs_concavity")
        if threshold is None:
            threshold = DEFAULT_MAX_ABS_CONCAVITY
            print(f"low_concavity_max_abs_concavity not found in {params_path.name}; "
                  f"drawing lsdtt3's default {threshold:g}")
        else:
            print(f"low_concavity_max_abs_concavity: {threshold:g} (from {params_path.name})")

        points = load_low_concavity_points(lc["points"]) if lc["points"].exists() else None
        if points is not None:
            theta, flag = points["theta"].to_numpy(dtype=float), points["low_concavity"].to_numpy()
            source = lc["points"].name
        elif lc["theta"].exists() and lc["flag"].exists():
            pix = load_low_concavity_rasters(lc["theta"], lc["flag"])
            theta, flag = pix["theta"].to_numpy(dtype=float), pix["low_concavity"].to_numpy()
            source = f"{lc['theta'].name} and {lc['flag'].name}"
        else:
            theta = flag = None
            print("no theta distribution: it needs the points CSV or the theta and flag rasters "
                  "(write_low_concavity_point_table or write_low_concavity_rasters: true)")
        if theta is not None:
            print(f"theta from {source}: {len(theta):,} channel pixels, "
                  f"{int((flag == 1).sum()):,} low concavity, {int((flag == 0).sum()):,} not, "
                  f"{int((flag == -1).sum()):,} no valid fit")
            fig, _ = plot_theta_distribution(theta, flag, threshold, theta_range=args.theta_range,
                                             bins=args.theta_bins)
            _save(fig, out_dir / f"{stem}_theta_distribution.png", args.dpi, written)

        if points is not None:
            fig, _ = plot_low_concavity_slope_area(points)
            _save(fig, out_dir / f"{stem}_low_concavity_slope_area.png", args.dpi, written)
        else:
            print(f"no _low_concavity_slope_area.png: {lc['points'].name} not found "
                  "(write_low_concavity_point_table: true in the run)")

        if args.dem:
            if lc["segments"].exists():
                written.append(render_low_concavity_map(
                    args.dem, lc["segments"], out_dir / f"{stem}_low_concavity_map.png",
                    width_by_order=not args.no_width_by_order, coordinates=args.coordinates,
                    title=f"low concavity: |theta| < {threshold:g}"))
            else:
                print(f"no _low_concavity_map.png: {lc['segments'].name} not found "
                      "(write_low_concavity_segments: true in the run)")
        else:
            print("no _low_concavity_map.png: give --dem for the hillshade")
    for w in written:
        print(w)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
