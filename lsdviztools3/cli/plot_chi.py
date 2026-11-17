"""CLI: lsdviz3-chi — k_sn map, stacked chi profiles and k_sn distributions from an
lsdtt-chi-analysis run."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..workflows import render_chi_ksn_map, render_chi_profiles, render_chi_removed_map
from . import _common
from ._chi_common import (add_ksn_class_basis_arg, add_profile_layout_args, add_tributary_selection_args,
                          ksn_classes_from_args, mn_suffix, parse_basins, parse_gradient_inset,
                          tributary_selection_from_args)


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-chi``.

    Returns:
        The configured parser (the run's directory, prefix and m/n, the DEM for
        the hillshade, and which figures to draw).
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-chi",
        description=(
            "Plot an lsdtt-chi-analysis run (fit_chi_networks: true): the fitted "
            "chi network coloured by log10 k_sn over a hillshade, a stacked "
            "chi-elevation plot of every basin, and optionally one chi profile "
            "per basin."
        ),
    )
    parser.add_argument("--data-dir", required=True,
                        help="Directory holding the lsdtt-chi-analysis outputs.")
    parser.add_argument("--prefix", required=True, help="The chi run's write_fname.")
    parser.add_argument("--m-over-n", required=True,
                        help="m/n of the run, as a number (0.45) or the file suffix (0p450).")
    parser.add_argument("--dem", default=None,
                        help="Projected DEM for the k_sn map's hillshade. Without it no map is drawn.")
    parser.add_argument("--out-dir", required=True, help="Directory for the PNGs.")
    parser.add_argument("--out-prefix", default=None,
                        help="Prefix of the PNG names (default: --prefix).")
    parser.add_argument("--per-basin", action="store_true",
                        help="Also write one chi-elevation profile per basin.")
    parser.add_argument("--no-stacked", action="store_true", help="Skip the stacked chi plot.")
    parser.add_argument("--cmap", default="viridis", help="Colour map for k_sn (default viridis).")
    parser.add_argument("--ksn-vmin", type=float, default=None,
                        help="Lower end of the k_sn colour scale, in k_sn units (the scale is "
                             "log10 k_sn; default the data's minimum). Lower values take the end colour.")
    parser.add_argument("--ksn-vmax", type=float, default=None,
                        help="Upper end of the k_sn colour scale, in k_sn units (default the data's "
                             "maximum). Higher values take the end colour, with an arrow on the colourbar.")
    add_tributary_selection_args(parser)
    parser.add_argument("--coordinates", default="local",
                        help="Map coordinates: local/projected (default) or geographic.")
    parser.add_argument("--dpi", type=int, default=300, help="DPI of the profile plots (default 300).")
    _common.add_locator_args(parser)
    _common.add_class_args(parser, flag="--ksn-classes",
                           what="k_sn on the map and the chi profiles (breaks in k_sn units; "
                                "schemes work in log10 k_sn)")
    add_ksn_class_basis_arg(parser)
    add_profile_layout_args(parser)
    parser.add_argument("--no-ksn-distribution", action="store_true",
                        help="Skip the k_sn distribution plot (per channel pixel and per segment).")
    return parser


def main(argv=None) -> int:
    """Entry point for ``lsdviz3-chi``.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    for name in ("ksn_vmin", "ksn_vmax"):
        v = getattr(args, name)
        if v is not None and not v > 0:
            raise SystemExit(f"--{name.replace('_', '-')} must be > 0 (k_sn, on a log colour scale).")
    if args.ksn_vmin is not None and args.ksn_vmax is not None and not args.ksn_vmin < args.ksn_vmax:
        raise SystemExit("--ksn-vmin must be below --ksn-vmax.")
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mn = mn_suffix(args.m_over_n)
    stem = args.out_prefix or args.prefix
    base = data_dir / f"{args.prefix}_chi"
    outputs = {
        "fit_nodes": Path(f"{base}_fit_nodes_mn_{mn}.csv"),
        "fit_segments": Path(f"{base}_fit_segments_mn_{mn}.csv"),
        "segment_lines": Path(f"{base}_fit_segment_lines_mn_{mn}.fgb"),
    }
    for key in ("fit_nodes", "fit_segments"):
        if not outputs[key].exists():
            raise SystemExit(f"{outputs[key]} not found: was the chi run made with "
                             f"fit_chi_networks: true and write_fit_tables: true at m/n {mn}?")

    removed = tributary_selection_from_args(args, base, mn, out_dir, stem)

    from ..io.chi import load_chi_fit, load_chi_fit_segments
    segments = load_chi_fit_segments(outputs["fit_segments"])
    nodes = load_chi_fit(outputs["fit_nodes"], outputs["fit_segments"])
    ksn_classes = ksn_classes_from_args(args, segments, nodes, removed)

    written = []
    if not args.no_ksn_distribution:
        import matplotlib.pyplot as plt

        from ..charts.ksn_distribution import plot_ksn_distributions
        fig, _ = plot_ksn_distributions(
            segments, nodes, removed=removed, classes=ksn_classes,
            class_basis=args.ksn_class_basis if ksn_classes is not None and ksn_classes.scheme != "manual" else None,
            cmap=args.cmap, vmin=args.ksn_vmin, vmax=args.ksn_vmax)
        written.append(out_dir / f"{stem}_ksn_distribution.png")
        fig.savefig(written[-1], dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
    if args.dem:
        if not outputs["segment_lines"].exists():
            raise SystemExit(f"{outputs['segment_lines']} not found: the k_sn map needs "
                             f"write_chi_fit_segment_lines: true in the chi run")
        written.append(render_chi_ksn_map(
            args.dem, outputs["segment_lines"], out_dir / f"{stem}_ksn_map.png",
            coordinates=args.coordinates, cmap=args.cmap,
            ksn_vmin=args.ksn_vmin, ksn_vmax=args.ksn_vmax, removed=removed,
            locator=_common.locator_kwargs(args), ksn_classes=ksn_classes))
        if removed is not None:
            written.append(render_chi_removed_map(
                args.dem, outputs["segment_lines"], removed,
                out_dir / f"{stem}_removed_tributaries_map.png",
                coordinates=args.coordinates, cmap=args.cmap,
                ksn_vmin=args.ksn_vmin, ksn_vmax=args.ksn_vmax))
    written += render_chi_profiles(
        outputs, out_dir, stem=f"{stem}_chi", dpi=args.dpi,
        stacked=not args.no_stacked, per_basin=args.per_basin,
        ksn_vmin=args.ksn_vmin, ksn_vmax=args.ksn_vmax,
        removed=removed, show_removed=args.show_removed, ksn_classes=ksn_classes,
        basins=parse_basins(args.basins), basins_per_row=args.basins_per_row,
        gradient_inset=parse_gradient_inset(args.ksn_gradient_inset),
        gradient_inset_loc=args.gradient_inset_loc)
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
