"""CLI: lsdviz3-cosmo-compare — render cosmogenic comparison figures."""

from __future__ import annotations

import argparse

from ..workflows import render_cosmo_comparison_figures


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-cosmo-compare``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-cosmo-compare",
        description=(
            "Render standard lsdtt3 cosmogenic comparison figures from "
            "single-point scan CSVs (+ optional sample-comparison CSV)."
        ),
    )
    parser.add_argument(
        "--scan-csv",
        required=True,
        help="Path to single_point_scan_comparison.csv.",
    )
    parser.add_argument(
        "--sample-csv",
        default=None,
        help="Optional sample-comparison CSV path (e.g. OCTOPUS/riversand summary).",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for figure PNGs.",
    )
    parser.add_argument(
        "--reference",
        default="cronus_v3",
        choices=["cronus_v3", "cronuscalc"],
        help="Reference to plot residuals against in sweep panels.",
    )
    parser.add_argument(
        "--pressure",
        type=float,
        default=None,
        help="Representative pressure slice in hPa (nearest available).",
    )
    parser.add_argument(
        "--erate",
        type=float,
        default=None,
        help="Representative erosion-rate slice in g cm^-2 yr^-1 (nearest available).",
    )
    parser.add_argument(
        "--stem",
        default="cosmo_comparison",
        help="Output filename stem.",
    )
    parser.add_argument("--dpi", type=int, default=220, help="Output PNG dpi.")
    return parser


def main(argv=None) -> int:
    """Entry point for cosmo comparison chart generation."""
    args = build_parser().parse_args(argv)

    out = render_cosmo_comparison_figures(
        scan_comparison_csv=args.scan_csv,
        sample_comparison_csv=args.sample_csv,
        out_dir=args.outdir,
        reference=args.reference,
        representative_pressure_hpa=args.pressure,
        representative_erate_gcm2yr=args.erate,
        stem=args.stem,
        dpi=args.dpi,
    )
    for key, path in out.items():
        print(f"{key}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
