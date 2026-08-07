"""CLI: lsdviz3-cosmo-transient — render transient (step-change) CRN figures."""

from __future__ import annotations

import argparse

from ..workflows import render_cosmo_transient_figures


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-cosmo-transient``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-cosmo-transient",
        description=(
            "Render transient cosmogenic figure families from "
            "lsdtt-cosmo-column-tool step-change sweep CSVs: normalised "
            "concentration vs time, and apparent erosion-rate ratios vs time."
        ),
    )
    parser.add_argument(
        "--sweep-csv",
        required=True,
        nargs="+",
        help="One or more cosmo_transient_sweep_<label>.csv files.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for figure PNGs.",
    )
    parser.add_argument(
        "--value-prefix",
        default="approach_frac_",
        choices=["approach_frac_", "norm_conc_", "norm_full_"],
        help=(
            "Concentration normalisation column prefix for the concentration "
            "figure (default: approach_frac_, which rises 0 -> 1)."
        ),
    )
    parser.add_argument(
        "--families",
        default=["concentration", "erate_ratio"],
        nargs="+",
        choices=["concentration", "erate_ratio"],
        help="Which figure families to render.",
    )
    parser.add_argument(
        "--concentration-xmax",
        type=float,
        default=10000.0,
        help="Upper x limit (yr) for the concentration figure.",
    )
    parser.add_argument(
        "--erate-ratio-xmax",
        type=float,
        default=50000.0,
        help="Upper x limit (yr) for the erosion-rate ratio figure.",
    )
    parser.add_argument(
        "--no-title",
        action="store_true",
        help="Omit the S_t / erosion ratio / final epsilon title.",
    )
    parser.add_argument("--dpi", type=int, default=400, help="Output PNG dpi.")
    return parser


def main(argv=None) -> int:
    """Entry point for transient cosmogenic chart generation."""
    args = build_parser().parse_args(argv)

    out = render_cosmo_transient_figures(
        args.sweep_csv,
        args.outdir,
        value_prefix=args.value_prefix,
        families=tuple(args.families),
        concentration_xlim=(0.0, args.concentration_xmax),
        erate_ratio_xlim=(0.0, args.erate_ratio_xmax),
        show_title=not args.no_title,
        dpi=args.dpi,
    )
    for key, path in out.items():
        print(f"{key}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
