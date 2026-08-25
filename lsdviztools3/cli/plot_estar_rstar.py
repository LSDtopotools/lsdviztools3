"""CLI: lsdviz3-estar-rstar — E*-R* validation against nonlinear hillslope theory.

Reads a pooled hilltop-metrics table (one or many landscapes) and writes the
E*-R* figure plus the per-group residual summary.

Nothing is recomputed here beyond the theoretical curve, which is closed-form:
E*, R* and the trace statuses all come from lsdtt3's own output.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ..charts.estar_rstar import plot_estar_rstar, summarise


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-estar-rstar``.

    Returns:
        The configured parser.
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-estar-rstar",
        description=(
            "Plot measured hilltop E*-R* against the Roering et al. (2007) "
            "steady nonlinear-hillslope relation, with the residual panel that "
            "says which term is at fault."
        ),
    )
    parser.add_argument("--metrics", required=True,
                        help="Hilltop metrics CSV (pooled across landscapes, "
                             "or a single run).")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for the PNG and the summary CSV.")
    parser.add_argument("--group-col", default="landscape",
                        help="Column to colour by (default: landscape). Use "
                             "U_over_D to check that only the ratio matters.")
    parser.add_argument("--title", default=None, help="Figure title.")
    parser.add_argument("--dpi", type=int, default=160)
    return parser


def main(argv=None) -> int:
    """Entry point for ``lsdviz3-estar-rstar``.

    Returns:
        0 on success; 1 if the table carried no usable rows.
    """
    args = build_parser().parse_args(argv)

    df = pd.read_csv(args.metrics)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    group = args.group_col if args.group_col in df.columns else None

    try:
        fig = plot_estar_rstar(df, group_col=group, title=args.title)
    except ValueError as exc:
        print(f"nothing to plot: {exc}")
        return 1

    png = outdir / "estar_rstar.png"
    fig.savefig(png, dpi=args.dpi, bbox_inches="tight")
    print(f"wrote {png}")

    table = summarise(df, group_col=group)
    csv_path = outdir / "estar_rstar_summary.csv"
    table.to_csv(csv_path, index=False)
    print(f"wrote {csv_path}\n")
    with pd.option_context("display.width", 200,
                           "display.max_columns", 20):
        print(table.to_string(index=False))

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
