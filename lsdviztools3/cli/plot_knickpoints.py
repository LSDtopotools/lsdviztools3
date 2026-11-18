"""CLI: lsdviz3-knickpoints - knickpoints of an lsdtt-chi-analysis fit.

Reads the run's knickpoint table (one row per node where two fitted k_sn
segments meet), its suggested thresholds, and the fit segments, and writes:

  <out-prefix>_knickpoint_distribution.png  the signed k_sn changes with the
      suggested thresholds and the chosen one, and the ranked |change|
  <out-prefix>_knickpoint_profile.png  major knickpoints on the stacked chi
      profile (or one basin with --basin)
  <out-prefix>_knickpoints_by_elevation.png  the elevations of the major
      knickpoints against all candidates, overall and per basin
  <out-prefix>_knickpoint_order_kde.png  kernel density of the elevation and
      chi of each channel's 1st, 2nd, 3rd, 4th strongest knickpoint, all basins
  <out-prefix>_knickpoints_along_chi.png  where major knickpoints sit along chi
  <out-prefix>_knickpoint_map.png  with --dem: major knickpoints on the map
  <out-prefix>_major_knickpoints.csv  the selected rows and the threshold used
  <out-prefix>_step_long_profiles.png  with steps and the chi data maps: steps
      on the longitudinal profile (elevation against flow distance)
  <out-prefix>_plotted_steps.csv  with steps: the steps drawn, after filtering

Knickzone reaches are drawn on the map when the run wrote them (and on the
profile with --profile-knickzones). Steps (run_step_change_detection) are
drawn on the longitudinal profiles only. Knickpoints
are ranked by lsdtt3 on the ordinary (stage 1) fit only, so the profile and
map draw that fit too; --stepaware-ksn draws the k_sn of the step-aware fit
under the same knickpoints instead. The chi run needs write_chi_knickpoint_table (and
write_chi_fit_segment_lines for the map); both come from a current lsdtt3.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from . import _common
from ._chi_common import (add_ksn_class_basis_arg, add_profile_layout_args, add_tributary_selection_args,
                          ksn_classes_from_args, mn_suffix, parse_basins, parse_gradient_key,
                          tributary_selection_from_args)


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
    p.add_argument("--measure", choices=["split_normalised", "split", "delta", "log_ratio"], default=None,
                   help="split_normalised: the split's chord-weighted k_sn change over the branch chi span, "
                        "lsdtt3's ranking (default; delta for a table without split columns); split: "
                        "the chord-weighted change (m); delta: k_sn upstream minus downstream; log_ratio: "
                        "log10 of their ratio.")
    p.add_argument("--threshold", default=None,
                   help="Major-knickpoint threshold on |measure|: major (lsdtt3's is_major, default with "
                        "split_normalised), a value, or a method from the run's thresholds table: p90, "
                        "p95 (default for the other measures), median_plus_3mad, otsu.")
    p.add_argument("--stepaware-ksn", action="store_true",
                   help="Draw the k_sn of the step-aware (stage 2) fit on the profile and map. The "
                        "knickpoints are the stage 1 ones either way.")
    p.add_argument("--no-knickzones", action="store_true", help="Do not draw knickzones anywhere.")
    p.add_argument("--no-map-knickzones", action="store_true",
                   help="Leave the knickzone reaches off the map (drawn by default).")
    p.add_argument("--profile-knickzones", action="store_true",
                   help="Also draw knickzones (foot to lip) on the chi profile (off by default).")
    p.add_argument("--no-steps", action="store_true", help="Do not draw step changes.")
    p.add_argument("--max-order", type=int, default=4,
                   help="Orders drawn in the knickpoint-order KDE (1st..4th strongest per channel; default 4).")
    p.add_argument("--order-majors-only", action="store_true",
                   help="Order only lsdtt3's major knickpoints in the KDE (default: every boundary).")
    p.add_argument("--min-step-drop", type=float, default=0.0,
                   help="Draw only steps that drop at least this much (m; default 0, all the run found).")
    p.add_argument("--hide-flagged-steps", action="store_true",
                   help="Leave out steps at a basin outlet or tributary mouth, or with uphill pairs "
                        "(drawn grey otherwise).")
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
                           what="k_sn on the map and the profile (breaks in k_sn units; schemes "
                                "work in log10 k_sn)")
    add_ksn_class_basis_arg(p)
    add_profile_layout_args(p)
    add_tributary_selection_args(p)
    return p


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
    from ..io.chi import chi_fit_files
    files = chi_fit_files(base, mn, stage="stepaware" if args.stepaware_ksn else "stage1")
    files["segments"] = files["fit_segments"]
    print(f"k_sn drawn from the {'step-aware (stage 2)' if files['stage'] == 'stepaware' else 'ordinary (stage 1)'}"
          " fit; knickpoints from stage 1")
    for key in ("knickpoints", "thresholds", "segments"):
        if not files[key].exists():
            raise SystemExit(f"{files[key]} not found: run lsdtt-chi-analysis (current lsdtt3) with "
                             f"fit_chi_networks, write_fit_tables and write_chi_knickpoint_table at m/n {mn}.")
    out_dir.mkdir(parents=True, exist_ok=True)

    from ..charts.knickpoints import (plot_knickpoint_distribution, plot_knickpoints_along_chi,
                                      plot_knickpoints_on_profile, render_knickpoint_map,
                                      select_major_knickpoints)
    from ..charts.knickpoints import (knickpoint_selection_order, plot_knickpoint_order_kde,
                                      plot_knickpoints_by_elevation, plot_steps_on_long_profile)
    from ..io.chi import (filter_chi_steps, load_chi_fit_segments, load_chi_knickpoint_thresholds,
                          load_chi_knickpoints, load_chi_knickzones, load_chi_steps, removed_mask)

    removed = tributary_selection_from_args(args, base, mn, out_dir, stem)
    kp = load_chi_knickpoints(files["knickpoints"])
    th = load_chi_knickpoint_thresholds(files["thresholds"])
    if args.measure is None:
        args.measure = "split_normalised" if "split_delta_ksn_normalised" in kp.columns else "delta"
    if args.threshold is None:
        args.threshold = "major" if args.measure == "split_normalised" else "p95"
    # --basin / --basins limit every figure (and the CSVs) to those basins
    basin_sel = [args.basin] if args.basin is not None else parse_basins(args.basins)
    if basin_sel is not None:
        kp = kp[kp["basin_index"].isin(basin_sel)]
    knickzones = None
    if not args.no_knickzones and files["knickzones"].exists():
        knickzones = load_chi_knickzones(files["knickzones"])
        if basin_sel is not None:
            knickzones = knickzones[knickzones["basin_index"].isin(basin_sel)]
        if removed is not None and len(removed):
            knickzones = knickzones[~removed_mask(knickzones, removed)]
        print(f"{len(knickzones)} knickzones")
    steps = None
    if not args.no_steps and files["steps"].exists():
        steps = filter_chi_steps(load_chi_steps(files["steps"]), min_drop=args.min_step_drop,
                                 hide_flagged=args.hide_flagged_steps, removed=removed)
        if basin_sel is not None:
            steps = steps[steps["basin_index"].isin(basin_sel)]
        print(f"{len(steps)} steps drawn ({int(steps['flagged'].sum())} flagged)")
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
    if steps is not None:
        written.append(out_dir / f"{stem}_plotted_steps.csv")
        steps.to_csv(written[-1], index=False)

    import matplotlib.pyplot as plt
    fig, _ = plot_knickpoint_distribution(kp, th, measure=args.measure, chosen=value, chosen_label=label,
                                          include_junctions=args.include_junctions, removed=removed)
    written.append(out_dir / f"{stem}_knickpoint_distribution.png")
    fig.savefig(written[-1], dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    segments = load_chi_fit_segments(files["segments"])
    nodes_path = Path(files["fit_nodes"])

    def load_nodes():
        if not nodes_path.exists():
            raise SystemExit(f"{nodes_path} not found: --ksn-class-basis pixels needs the fit nodes.")
        from ..io.chi import load_chi_fit
        return load_chi_fit(nodes_path, files["segments"])

    ksn_classes = ksn_classes_from_args(args, segments, load_nodes, removed)
    fig, _ = plot_knickpoints_on_profile(segments, major, basin=args.basin, removed=removed,
                                         show_removed=args.show_removed, vmin=args.ksn_vmin,
                                         vmax=args.ksn_vmax, threshold_label=label,
                                         classes=ksn_classes, basins=parse_basins(args.basins),
                                         basins_per_row=args.basins_per_row,
                                         gradient_key=parse_gradient_key(args.ksn_gradient_key),
                                         knickzones=knickzones if args.profile_knickzones else None)
    written.append(out_dir / (f"{stem}_knickpoint_profile_basin{args.basin}.png" if args.basin is not None
                              else f"{stem}_knickpoint_profile.png"))
    fig.savefig(written[-1], dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    fig, _ = plot_knickpoints_by_elevation(kp, major, include_junctions=args.include_junctions,
                                           removed=removed, threshold_label=label)
    written.append(out_dir / f"{stem}_knickpoints_by_elevation.png")
    fig.savefig(written[-1], dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    try:
        ordered = knickpoint_selection_order(kp, measure=args.measure, include_junctions=args.include_junctions,
                                             removed=removed, majors_only=args.order_majors_only)
        fig, _ = plot_knickpoint_order_kde(ordered, max_order=args.max_order, title=stem)
        written.append(out_dir / f"{stem}_knickpoint_order_kde.png")
        fig.savefig(written[-1], dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
    except ValueError as exc:
        print(f"no knickpoint-order KDE: {exc}")

    fig, _ = plot_knickpoints_along_chi(kp, major, along=args.along,
                                        include_junctions=args.include_junctions, removed=removed)
    written.append(out_dir / f"{stem}_knickpoints_along_chi.png")
    fig.savefig(written[-1], dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    if steps is not None and len(steps):
        dm = next((q for q in (Path(f"{base}_data_maps.csv"), Path(f"{base}_data_maps.fgb"))
                   if q.exists()), None)
        if dm is None or not Path(files["fit_nodes"]).exists():
            print("no step long profiles: they need the fit nodes and the chi data maps "
                  "(write_chi_data_maps: true)")
        else:
            from ..io.chi import attach_flow_distance, load_chi_fit_nodes, load_chi_points
            nodes = attach_flow_distance(load_chi_fit_nodes(files["fit_nodes"]), load_chi_points(dm))
            fig, _ = plot_steps_on_long_profile(nodes, steps, basins=basin_sel, removed=removed)
            written.append(out_dir / f"{stem}_step_long_profiles.png")
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
            locator=_common.locator_kwargs(args), ksn_classes=ksn_classes,
            knickzones=None if args.no_map_knickzones else knickzones, basins=basin_sel))
    for w in written:
        print(w)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
