"""CLI: lsdviz3-cosmo-published — plot residuals from published-summary CSVs."""

from __future__ import annotations

import argparse

from ..workflows import render_cosmo_published_residuals


def build_parser() -> argparse.ArgumentParser:
    """Build parser for ``lsdviz3-cosmo-published``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-cosmo-published",
        description=(
            "Render grouped residual boxplots from lsdtt3 published summary CSVs "
            "(stone/lsdn/riversand)."
        ),
    )
    parser.add_argument(
        "--published-dir",
        required=True,
        help="Directory containing <dataset>_summary.csv style files.",
    )
    parser.add_argument(
        "--datasets",
        required=True,
        help="Comma-separated dataset names (e.g. dethier2014,clapp2002).",
    )
    parser.add_argument(
        "--mode",
        default="auto",
        choices=["auto", "stone", "lsdn", "riversand"],
        help="Summary flavour selection mode.",
    )
    parser.add_argument("--outdir", required=True, help="Output figure directory.")
    parser.add_argument("--stem", default="cosmo_published", help="Output filename stem.")
    parser.add_argument("--dpi", type=int, default=220, help="Output PNG dpi.")
    return parser


def main(argv=None) -> int:
    """Entry point for published-summary residual plotting."""
    args = build_parser().parse_args(argv)
    datasets = [x.strip() for x in args.datasets.split(",") if x.strip()]
    out = render_cosmo_published_residuals(
        published_dir=args.published_dir,
        out_dir=args.outdir,
        datasets=datasets,
        mode=args.mode,
        stem=args.stem,
        dpi=args.dpi,
    )
    for key, val in out.items():
        if isinstance(val, list):
            print(f"{key}:")
            for p in val:
                print(f"  - {p}")
        else:
            print(f"{key}: {val}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
