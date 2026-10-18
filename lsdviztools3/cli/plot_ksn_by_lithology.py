"""CLI: lsdviz3-ksn-by-lithology — k_sn histograms and box plots by lithologic unit.

Reads the lithology outputs of ``lsdtt-chi-analysis`` for one run
(``<prefix>_chi_ksn_{by_lithology,histogram_by_lithology,lithology_nodes}_mn_<mn>.csv``
and, when present, the reach and unit x basin tables) and writes:

* ``<out-prefix>_ksn_hist_length.png`` — length-weighted histograms;
* ``<out-prefix>_ksn_hist_count.png``  — node-count histograms;
* ``<out-prefix>_ksn_box.png``         — per-unit box plots with sample sizes;
* ``<out-prefix>_reach_ksn_box.png``   — per-reach k_sn (if the reach table exists);
* ``<out-prefix>_ksn_unit_basin.png``  — unit x basin dots (if that table exists).
"""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-ksn-by-lithology``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-ksn-by-lithology",
        description="Plot k_sn by lithologic unit from lsdtt-chi-analysis "
        "lithology CSVs.",
    )
    parser.add_argument("--data-dir", required=True,
                        help="Directory holding the lsdtt3 outputs.")
    parser.add_argument("--prefix", required=True,
                        help="The chi run's write_fname (e.g. wvaf_lithprimary).")
    parser.add_argument("--m-over-n", default=None,
                        help="m/n of the run (float, or suffix like 0p450). "
                        "Default: the only one present.")
    parser.add_argument("--out-dir", required=True, help="Directory for the PNGs.")
    parser.add_argument("--out-prefix", default=None,
                        help="Output file prefix (default: --prefix).")
    parser.add_argument("--order", default="auto",
                        choices=["auto", "median", "reach_median", "code", "name", "length"],
                        help="Unit order. auto (default): node-level figures by the "
                        "summary node median, reach-level figures by the summary "
                        "reach median.")
    parser.add_argument("--layout", default="auto",
                        choices=["auto", "overlay", "facets"],
                        help="Histogram layout (default auto: overlay up to 6 units).")
    parser.add_argument("--raw", action="store_true",
                        help="Histograms in raw length/count rather than "
                        "fraction of each unit's total.")
    parser.add_argument("--box-source", default="auto",
                        choices=["auto", "nodes", "summary"],
                        help="Box plot from per-node values (auto when the node "
                        "table exists) or from the summary quartiles.")
    parser.add_argument(
        "--unit-colours", default=None,
        help="name,colour CSV or JSON {name: colour} fixing unit colours by NAME "
        "(the same file used for lsdviz3-lithology-map).")
    parser.add_argument("--title", default=None, help="Title prefix for every figure.")
    parser.add_argument("--dpi", type=int, default=200, help="Output DPI (default 200).")
    return parser


def main(argv=None) -> int:
    """Entry point: write the k_sn-by-lithology figures."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ..charts.lithology_ksn import (
        plot_ksn_box_by_lithology, plot_ksn_by_lithology_and_basin,
        plot_ksn_histogram_by_lithology, plot_reach_ksn_by_lithology,
    )
    from ..io.lithology import load_ksn_lithology_outputs, read_unit_colours

    args = build_parser().parse_args(argv)
    mn = args.m_over_n
    if mn is not None and "p" not in mn:
        mn = float(mn)
    data = load_ksn_lithology_outputs(args.data_dir, args.prefix, mn)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.out_prefix or args.prefix
    written = []

    def _save(fig, name):
        path = out_dir / f"{stem}_{name}.png"
        fig.savefig(path, dpi=args.dpi)
        plt.close(fig)
        written.append(path)

    s = data["summary"]
    uc = {"unit_colours": read_unit_colours(args.unit_colours)} if args.unit_colours else {}
    node_order = "median" if args.order == "auto" else args.order
    reach_order = "reach_median" if args.order == "auto" else args.order
    if data["histogram"] is not None:
        for weight in ("length", "count"):
            fig, _ = plot_ksn_histogram_by_lithology(
                data["histogram"], s, weight=weight, normalise=not args.raw,
                layout=args.layout, order=node_order, title=args.title, **uc)
            _save(fig, f"ksn_hist_{weight}")
    fig, _ = plot_ksn_box_by_lithology(s, data["nodes"], source=args.box_source,
                                       order=node_order, title=args.title, **uc)
    _save(fig, "ksn_box")
    if data["reaches"] is not None:
        fig, _ = plot_reach_ksn_by_lithology(data["reaches"], s, order=reach_order,
                                             title=args.title, **uc)
        _save(fig, "reach_ksn_box")
    if data["unit_basin"] is not None:
        fig, _ = plot_ksn_by_lithology_and_basin(data["unit_basin"], s,
                                                 order=reach_order, title=args.title, **uc)
        _save(fig, "ksn_unit_basin")
    for p in written:
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
