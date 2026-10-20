"""CLI: lsdviz3-channel-bearings — bearing roses and junction angles by group.

Reads the grouped outputs of lsdtt-channel-network-analysis for one run
(``<prefix>_channel_bearing_{rose,group_summary,comparisons}.csv``,
``<prefix>_junction_angle_{group_summary,comparisons}.csv``,
``<prefix>_junction_angles.csv``) and writes, into ``--out-dir``:

* ``<prefix>_rose_<axial|directional>_<count|length>[_<stratum>].png`` — one
  small rose per group (fractions of count / of channel length);
* ``<prefix>_bearing_comparisons_<type>.png`` — lsdtt3's pairwise table;
* ``<prefix>_junction_<measurand>_summary.png`` — boxes from the summary;
* ``<prefix>_junction_<measurand>_distributions.png`` — per-junction
  histograms by order (order grouping only);
* ``<prefix>_junction_<measurand>_comparisons.png``;
* with ``--compare-prefix``: ``<prefix>_vs_<other>_rose_*.png`` and
  ``<prefix>_vs_<other>_junction_<measurand>_summary.png``;
* with ``--map-dem``: ``<prefix>_channel_order_map.png``. With
  ``<prefix>_channel_network_lines.fgb`` (CNA ``write_channel_network_lines``) the
  whole network is drawn, analysed links thick and the rest thin, with
  confluences as circles (filled = analysed; needs
  ``write_junction_network_points``); otherwise the analysed links' bearing
  vectors;
* with ``--map-dem`` and ``--regions``: ``<prefix>_channel_region_map.png``, the
  same network over region polygons (e.g. physiographic provinces).

The same call serves every grouping lsdtt3 writes: ``--group-key
strahler_order`` (default), ``lithology``, ``lithology_interior``, ``basin``,
``region``; with ``--stratum-key`` one rose figure is written per stratum value.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-channel-bearings``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-channel-bearings",
        description="Plot lsdtt3 channel-bearing roses and junction-angle summaries by group.",
    )
    parser.add_argument("--data-dir", required=True, help="Directory with the lsdtt3 outputs.")
    parser.add_argument("--prefix", required=True, help="Run write_fname (e.g. wvbig_order).")
    parser.add_argument("--out-dir", required=True, help="Directory for the PNGs.")
    parser.add_argument("--label", default=None, help="Run label (default: --prefix).")
    parser.add_argument("--compare-prefix", default=None,
                        help="Second run to draw beside this one (e.g. an isotropic null).")
    parser.add_argument("--compare-label", default=None,
                        help="Label for --compare-prefix (default: its prefix).")
    parser.add_argument("--group-key", default="strahler_order",
                        help="group_key(s) to plot, comma-separated (strahler_order, "
                        "lithology, lithology_interior, ...). Non-default keys appear in "
                        "the output file names.")
    parser.add_argument("--stratum-key", default="none",
                        help="stratum_key to plot ('none' = unstratified).")
    parser.add_argument("--measurands", default="alpha",
                        help="Junction measurands, comma-separated (alpha,beta_min,beta_max).")
    parser.add_argument(
        "--plots", default="rose,junction,tables,blocks",
        help="Comma-separated figure kinds: rose, junction, tables, blocks, cdf, pdf "
        "(default: rose,junction,tables,blocks). cdf = empirical CDF of lsdtt3's "
        "per-link bearings (membership file); pdf = lsdtt3's binned fractions.")
    parser.add_argument("--pdf-weight", default="count", choices=["count", "length"],
                        help="PDF weight: fraction of links (count, default) or of "
                        "channel length.")
    parser.add_argument("--map-dem", default=None,
                        help="DEM for a channel map coloured by stream order (the whole "
                        "network from <prefix>_channel_network_lines.fgb if present, "
                        "else <prefix>_channel_link_bearing_vectors.fgb).")
    parser.add_argument("--regions", default=None,
                        help="Region polygons (e.g. the USGS physiographic divisions "
                        "physio.shp) for <prefix>_channel_region_map.png; needs --map-dem.")
    parser.add_argument("--region-field", default="PROVINCE",
                        help="Region name field in --regions (default PROVINCE).")
    parser.add_argument("--region-title", default="Physiographic province",
                        help="Legend header for the regions (default "
                        "'Physiographic province').")
    parser.add_argument("--title", default=None, help="Title prefix.")
    parser.add_argument("--dpi", type=int, default=200, help="Output DPI (default 200).")
    return parser


def _safe(value) -> str:
    """Stratum value as a file-name fragment (spaces/commas/slashes -> '-')."""
    import re

    return re.sub(r"[^A-Za-z0-9+._-]+", "-", str(value)).strip("-")


def _stratum_values(df, group_key, stratum_key):
    from ..io.channel_bearings import GROUP_COLUMNS, NO_STRATUM, select_groups

    if stratum_key == NO_STRATUM:
        return [None]
    sel = select_groups(df, group_key=group_key, stratum_key=stratum_key)
    return list(dict.fromkeys(sel[GROUP_COLUMNS["stratum_value"]].astype(str)))


def main(argv=None) -> int:
    """Entry point: write the bearing / junction figures for one run."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ..charts.channel_bearings import (
        plot_bearing_distributions, plot_comparison_table, plot_junction_angle_distributions,
        plot_junction_angle_summary, plot_rose_comparison,
    )
    from ..io.channel_bearings import (
        AXIAL, DIRECTIONAL, load_channel_geometry_outputs, select_groups,
    )

    args = build_parser().parse_args(argv)
    data = load_channel_geometry_outputs(args.data_dir, args.prefix)
    other = (load_channel_geometry_outputs(args.data_dir, args.compare_prefix)
             if args.compare_prefix else None)
    label = args.label or args.prefix
    olabel = args.compare_label or args.compare_prefix
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sk = args.stratum_key
    written = []

    def _save(fig, name):
        path = out_dir / f"{name}.png"
        fig.savefig(path, dpi=args.dpi)
        plt.close(fig)
        written.append(path)

    def _comparison_strata(df):
        """Stratified runs: tabulate the pooled 'all' stratum if lsdtt3 wrote one."""
        if sk == "none" or df is None or "stratum_value" not in df.columns:
            return [None]
        vals = list(dict.fromkeys(df["stratum_value"].astype(str)))
        return ["all"] if "all" in vals else vals

    measurands = [m.strip() for m in args.measurands.split(",") if m.strip()]
    plots = {p.strip() for p in args.plots.split(",") if p.strip()}
    unknown = plots - {"rose", "junction", "tables", "blocks", "cdf", "pdf"}
    if unknown:
        raise SystemExit(f"Unknown --plots value(s): {sorted(unknown)}")
    for gk in [g.strip() for g in args.group_key.split(",") if g.strip()]:
        # Group key in file names except for the default order grouping, so
        # several groupings (lithology, lithology_interior) can share a prefix.
        gtag = "" if gk == "strahler_order" else f"_{gk}"
        stem = f"{args.prefix}{gtag}"
        if "rose" in plots and data["rose"] is not None:
            for sv in _stratum_values(data["rose"], gk, sk):
                sel = {} if sv is None else {"stratum_key": sk, "stratum_value": sv}
                tag = "" if sv is None else f"_{sk}-{_safe(sv)}"
                for atype in (DIRECTIONAL, AXIAL):
                    for weight in ("fraction_count", "fraction_length"):
                        w = weight.replace("fraction_", "")
                        title = f"{args.title + ' ' if args.title else ''}{label}"
                        if sv is not None:
                            title += f", {sk} {sv}"
                        try:
                            fig, _ = plot_rose_comparison(
                                {"": (data["rose"], data["bearing_summary"])}, group_key=gk,
                                angular_type=atype, weight=weight, ncols=4, title=title, **sel)
                        except ValueError:
                            continue  # no rose rows for this grouping
                        _save(fig, f"{stem}_rose_{atype}_{w}{tag}")
                        if other is not None and other["rose"] is not None:
                            fig, _ = plot_rose_comparison(
                                {label: (data["rose"], data["bearing_summary"]),
                                 olabel: (other["rose"], other["bearing_summary"])},
                                group_key=gk, angular_type=atype, weight=weight,
                                title=args.title, **sel)
                            _save(fig, f"{stem}_vs_{args.compare_prefix}_rose_{atype}_{w}{tag}")
        has_rose = (data["rose"] is not None and not select_groups(
            data["rose"], group_key=gk, stratum_key=sk).empty)
        # Bearing comparison rows carry no group key; only tabulate them for a
        # grouping that actually has bearing rows.
        if "tables" in plots and data["bearing_comparisons"] is not None and has_rose:
            for sv in _comparison_strata(data["bearing_comparisons"]):
                for atype in (DIRECTIONAL, AXIAL):
                    try:
                        fig, _ = plot_comparison_table(
                            data["bearing_comparisons"], kind="bearing", angular_type=atype,
                            group_key=gk, stratum_value=sv,
                            title=f"{label}: {gk} bearing comparisons ({atype}"
                                  + (f", {sk} {sv}" if sv else "")
                                  + "); p-values assume i.i.d. samples, see n_blocks")
                    except ValueError:
                        continue
                    _save(fig, f"{stem}_bearing_comparisons_{atype}")

        for kind in [k for k in ("cdf", "pdf") if k in plots]:
            if not has_rose:
                break
            src = "membership" if kind == "cdf" else "rose"
            if data[src] is None:
                print(f"skip {kind}: no {src} file for {args.prefix}", file=sys.stderr)
                continue
            weight = f"fraction_{args.pdf_weight}"
            wtag = "" if (kind == "cdf" or args.pdf_weight == "count") else "_length"
            sets = [("", {label: data})]
            if other is not None and other[src] is not None:
                sets.append((f"_vs_{args.compare_prefix}", {label: data, olabel: other}))
            for sv in _stratum_values(data["rose"], gk, sk):
                tag = "" if sv is None else f"_{sk}-{_safe(sv)}"
                for atype in (AXIAL, DIRECTIONAL):
                    for vs, runs in sets:
                        title = f"{args.title + ' ' if args.title else ''}" + " vs ".join(runs)
                        try:
                            fig, _ = plot_bearing_distributions(
                                runs, kind=kind, group_key=gk, angular_type=atype,
                                stratum_key=sk, stratum_value=sv, weight=weight, title=title)
                        except ValueError as exc:
                            print(f"skip {stem}{vs} {kind} {atype}{tag}: {exc}", file=sys.stderr)
                            continue
                        _save(fig, f"{stem}{vs}_bearing_{kind}{wtag}_{atype}{tag}")

        js = data["junction_summary"] if "junction" in plots else None
        for m in measurands:
            if js is None:
                break
            try:
                if sk == "none":
                    fig, _ = plot_junction_angle_summary(
                        js, group_key=gk, measurand=m, title=f"{label}: junction {m} by {gk}")
                    _save(fig, f"{stem}_junction_{m}_summary")
                else:
                    # One figure: groups on x, one box per stratum value side by side.
                    strata = _stratum_values(js, gk, sk)
                    runs = {f"{sk.replace('_', ' ')} {v}":
                            js[js["stratum_value"].astype(str) == str(v)] for v in strata}
                    fig, _ = plot_junction_angle_summary(
                        runs, group_key=gk, measurand=m, stratum_key=sk,
                        title=f"{label}: junction {m} by {gk}, per {sk}")
                    _save(fig, f"{stem}_junction_{m}_by_{sk}")
            except ValueError:
                continue  # no rows for this grouping/measurand
            if other is not None and other["junction_summary"] is not None and sk == "none":
                fig, _ = plot_junction_angle_summary(
                    {label: js, olabel: other["junction_summary"]}, group_key=gk, measurand=m,
                    title=f"junction {m} by {gk}: {label} vs {olabel}")
                _save(fig, f"{stem}_vs_{args.compare_prefix}_junction_{m}_summary")
            if data["junctions"] is not None and gk == "strahler_order" and sk == "none":
                try:
                    fig, _ = plot_junction_angle_distributions(
                        data["junctions"], js, measurand=m, title=label)
                    _save(fig, f"{stem}_junction_{m}_distributions")
                except KeyError:
                    pass  # no per-junction column for this measurand
            if data["junction_comparisons"] is not None:
                for sv in _comparison_strata(data["junction_comparisons"]):
                    try:
                        fig, _ = plot_comparison_table(
                            data["junction_comparisons"], kind="junction", measurand=m,
                            group_key=gk, stratum_value=sv,
                            title=f"{label}: junction {m} by {gk}"
                                  + (f" ({sk} {sv})" if sv else "")
                                  + " comparisons; p-values assume i.i.d. samples, see n_blocks")
                    except ValueError:
                        continue
                    _save(fig, f"{stem}_junction_{m}_comparisons")

    if "blocks" in plots and data["block_order_counts"] is not None:
        from ..charts.channel_bearings import plot_block_order_counts

        fig, _ = plot_block_order_counts(data["block_order_counts"],
                                         title=f"{label}: permutation blocks per Strahler order")
        _save(fig, f"{args.prefix}_block_order_counts")

    if args.map_dem:
        from ..io.channel_bearings import load_network_status
        from ..render.base import save_figure
        from ..render.lithology_map import render_channel_order_map
        from ..render.style import MapStyle

        status = load_network_status(args.data_dir, args.prefix)
        vec = data["paths"]["link_vectors"]
        if status is not None or vec.exists():
            style = MapStyle(clean=True, transparent=False, figure_size="14c")
            maps = [("channel_order_map", None)]
            if args.regions:
                maps.append(("channel_region_map", args.regions))
            for name, regions in maps:
                fig = render_channel_order_map(
                    args.map_dem, None if status is not None else vec,
                    network_status=status, regions=regions,
                    region_field=args.region_field, region_header=args.region_title,
                    legend_title="Strahler order", channel_width="0.9p", style=style)
                path = out_dir / f"{args.prefix}_{name}.png"
                save_figure(fig, path, style=style)
                written.append(path)

    for p in written:
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
