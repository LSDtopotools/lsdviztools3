"""Render lsdtt3 cosmo comparison figures from standard CSV outputs.

Usage example:

    python examples/cosmo_comparison_workflow.py \
        --scan-csv scan_out/single_point_scan_comparison.csv \
        --outdir demo_output/cosmo

Optional sample-level table (OCTOPUS/riversand summary style):

    python examples/cosmo_comparison_workflow.py \
        --scan-csv scan_out/single_point_scan_comparison.csv \
        --sample-csv ../lsdtt3/tests/expected/cosmo/published/dethier2014_riversand_summary.csv \
        --outdir demo_output/cosmo
"""

from __future__ import annotations

import argparse
from pathlib import Path

from lsdviztools3.workflows import render_cosmo_comparison_figures


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scan-csv", required=True, help="single_point_scan_comparison.csv path")
    ap.add_argument("--sample-csv", help="optional sample-comparison CSV path")
    ap.add_argument("--outdir", default="demo_output/cosmo")
    ap.add_argument("--reference", default="cronus_v3", choices=["cronus_v3", "cronuscalc"])
    ap.add_argument("--pressure", type=float, help="representative pressure slice (nearest)")
    ap.add_argument("--erate", type=float, help="representative erosion-rate slice (nearest)")
    args = ap.parse_args()

    outputs = render_cosmo_comparison_figures(
        scan_comparison_csv=Path(args.scan_csv),
        sample_comparison_csv=Path(args.sample_csv) if args.sample_csv else None,
        out_dir=Path(args.outdir),
        reference=args.reference,
        representative_pressure_hpa=args.pressure,
        representative_erate_gcm2yr=args.erate,
    )

    for k, v in outputs.items():
        print(f"{k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
