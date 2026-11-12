"""CLI: lsdviz3-knickpoints - knickpoints of an lsdtt-chi-analysis fit.

Reads the run's knickpoint table (one row per node where two fitted k_sn
segments meet), its suggested thresholds, and the fit segments, and writes:

  <out-prefix>_knickpoint_distribution.png  the signed k_sn changes with the
      suggested thresholds and the chosen one, and the ranked |change|
  <out-prefix>_knickpoint_profile.png  major knickpoints on the stacked chi
      profile (or one basin with --basin)
  <out-prefix>_knickpoints_along_chi.png  where major knickpoints sit along chi
  <out-prefix>_knickpoint_map.png  with --dem: major knickpoints on the map
  <out-prefix>_major_knickpoints.csv  the selected rows and the threshold used

The chi run needs write_chi_knickpoint_table (and write_chi_fit_segment_lines
for the map); both come from a current lsdtt3.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from . import _common
from ._chi_common import add_tributary_selection_args, mn_suffix, tributary_selection_from_args


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-knickpoints``."""
    p = argparse.ArgumentParser(
        prog="lsdviz3-knickpoints",
        description="Knickpoints (where fitted k_sn segments meet) of an lsdtt-chi-analysis run: "
                    "their distribution with suggested thresholds, and the major ones on a map, a "
                    "profile and along chi.")
    p.add_argument("--data-dir", required=True, help="Directory holding the lsdtt-chi-analysis outputs.")
    p.add_argument("--prefix", required=True, help="The chi run's write_fname.")
    p.add_argument("--m-over-n", required=True, help="m/n of the run, 0.45 or the file suffix 0p450.")
    p.add_argument("--out-dir", required=True, help="Directory for the outputs.")
    p.add_argument("--out-prefix", default=None, help="Prefix of the output names (default --prefix).")
    p.add_argument("--dem", default=None, help="Projected DEM for the map. Without it no map is drawn.")
    p.add_argument("--measure", choices=["delta", "log_ratio"], default="delta",
                   help="delta: k_sn upstream minus downstream (default); log_ratio: log10 of their ratio.")
    p.add_argument("--threshold", default="p95",
                   help="Major-knickpoint threshold on |measure|: a value, or a method from the run's "
                        "thresholds table: p90, p95 (default), median_plus_3mad, otsu.")
    p.add_argument("--direction", choices=["both", "steepening", "shallowing"], default="both",
                   help="Which changes count: both (default), steepening upstream only, or shallowing only.")
    p.add_argument("--include-junctions", action="store_true",
                   help="Also count k_sn changes at tributary mouths (left out by default).")
    p.add_argument("--along", choices=["chi", "chi_normalised", "elevation"], default="chi",
                   help="x axis of the along-chi plot (default chi).")
    p.add_argument("--basin", type=int, default=None,
                   help="Draw one basin's profile instead of the stacked plot.")
    p.add_argument("--ksn-vmin", type=float, default=None, help="k_sn colour range lower end (k_sn units).")
    p.add_argument("--ksn-vmax", type=float, default=None, help="k_sn colour range upper end (k_sn units).")
    p.add_argument("--cmap", default="viridis", help="k_sn colour map (default viridis).")
    p.add_argument("--coordinates", default="local", help="Map coordinates: local/projected or geographic.")
    p.add_argument("--dpi", type=int, default=300, help="DPI of the charts (default 300).")
    _common.add_locator_args(p)
    _common.add_class_args(p, flag="--ksn-classes",
                           what="the k_sn network on the map (breaks in k_sn units; schemes work "
                                "in log10 k_sn)")
    add_tributary_selection_args(p)
    return p


def _ksn_class_kwargs(args) -> dict:
    """``--ksn-classes`` and friends as render_knickpoint_map keywords."""
    kw = _common.class_kwargs(args)
    kw["ksn_classes"] = kw.pop("classes")
    return kw


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    for name in ("ksn_vmin", "ksn_vmax"):
        v = getattr(args, name)
        if v is not None and not v > 0:
            raise SystemExit(f"--{name.replace('_', '-')} must be > 0.")
    data_dir, out_dir = Path(args.data_dir), Path(args.out_dir)
    mn = mn_suffix(args.m_over_n)
    stem = args.out_prefix or args.prefix
    base = data_dir / f"{args.prefix}_chi"
    files = {
        "knickpoints": Path(f"{base}_knickpoints_mn_{mn}.csv"),
        "thresholds": Path(f"{base}_knickpoint_thresholds_mn_{mn}.csv"),
        "segments": Path(f"{base}_fit_segments_mn_{mn}.csv"),
        "segment_lines": Path(f"{base}_fit_segment_lines_mn_{mn}.fgb"),
    }
    for key in ("knickpoints", "thresholds", "segments"):
        if not files[key].exists():
            raise SystemExit(f"{files[key]} not found: run lsdtt-chi-analysis (current lsdtt3) with "
                             f"fit_chi_networks, write_fit_tables and write_chi_knickpoint_table at m/n {mn}.")
    out_dir.mkdir(parents=True, exist_ok=True)

    from ..charts.knickpoints import (plot_knickpoint_distribution, plot_knickpoints_along_chi,
                                      plot_knickpoints_on_profile, render_knickpoint_map,
                                      select_major_knickpoints)
    from ..io.chi import load_chi_fit_segments, load_chi_knickpoint_thresholds, load_chi_knickpoints

    removed = tributary_selection_from_args(args, base, mn, out_dir, stem)
    kp = load_chi_knickpoints(files["knickpoints"])
    th = load_chi_knickpoint_thresholds(files["thresholds"])
    try:
        major, value, label = select_major_knickpoints(
            kp, args.threshold, measure=args.measure, direction=args.direction,
            include_junctions=args.include_junctions, thresholds=th, removed=removed)
    except ValueError as exc:
        raise SystemExit(str(exc))
    print(f"threshold {label} on |{args.measure}|: {len(major)} major knickpoints")
    major = major.assign(threshold_measure=args.measure, threshold=value, threshold_label=label,
                         direction=args.direction)
    csv = out_dir / f"{stem}_major_knickpoints.csv"
    major.to_csv(csv, index=False)
    written = [csv]

    import matplotlib.pyplot as plt
    fig, _ = plot_knickpoint_distribution(kp, th, measure=args.measure, chosen=value, chosen_label=label,
                                          include_junctions=args.include_junctions, removed=removed)
    written.append(out_dir / f"{stem}_knickpoint_distribution.png")
    fig.savefig(written[-1], dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    segments = load_chi_fit_segments(files["segments"])
    fig, _ = plot_knickpoints_on_profile(segments, major, basin=args.basin, removed=removed,
                                         show_removed=args.show_removed, vmin=args.ksn_vmin,
                                         vmax=args.ksn_vmax, threshold_label=label)
    written.append(out_dir / (f"{stem}_knickpoint_profile_basin{args.basin}.png" if args.basin is not None
                              else f"{stem}_knickpoint_profile.png"))
    fig.savefig(written[-1], dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    fig, _ = plot_knickpoints_along_chi(kp, major, along=args.along,
                                        include_junctions=args.include_junctions, removed=removed)
    written.append(out_dir / f"{stem}_knickpoints_along_chi.png")
    fig.savefig(written[-1], dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    if args.dem:
        if not files["segment_lines"].exists():
            raise SystemExit(f"{files['segment_lines']} not found: the map needs "
                             "write_chi_fit_segment_lines: true in the chi run.")
        written.append(render_knickpoint_map(
            args.dem, files["segment_lines"], major, out_dir / f"{stem}_knickpoint_map.png",
            removed=removed, coordinates=args.coordinates, cmap=args.cmap,
            ksn_vmin=args.ksn_vmin, ksn_vmax=args.ksn_vmax, threshold_label=label,
            locator=_common.locator_kwargs(args), **_ksn_class_kwargs(args)))
    for w in written:
        print(w)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
