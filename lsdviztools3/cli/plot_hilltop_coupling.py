"""CLI: lsdviz3-hilltop-coupling — figures for a hilltop-flow-routing run.

Reads one output set (hilltop points, ridgelines, traces, rasters) and writes
the Phase 0/1 inspection figures: the traces over a hillshade, the trace
outcomes, the ordered ridgelines with their two flanking sides, and the
per-divide hilltop profiles.

Nothing is recomputed here: the hillslope lengths, reliefs, outcomes and side
assignments all come from lsdtt3's own output.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..charts.estar_rstar import plot_estar_rstar
from ..charts.hilltop_coupling import (
    plot_length_distributions,
    plot_looped_trace_profiles,
    plot_ridgeline_profiles,
    plot_ridgelines_and_sides,
    plot_trace_outcomes,
    plot_traces_on_hillshade,
)
from ..io.hilltop_coupling import load_hilltop_outputs


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-hilltop-coupling``.

    Returns:
        The configured parser (``--prefix`` output set, ``--outdir``, ``--dpi``).
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-hilltop-coupling",
        description=(
            "Render hilltop / channel coupling outputs: traces over a "
            "hillshade, trace-outcome diagnostics, the ordered ridgelines with "
            "their two flanking hillslopes, and hilltop coverage along each "
            "divide."
        ),
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Directory holding one hilltop-flow-routing output set, or a "
        "<write_path>/<write_fname> stem.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for the figure PNGs.",
    )
    parser.add_argument(
        "--dem",
        default=None,
        help="DEM raster for the map backgrounds. Only needed when the output "
             "set does not sit beside a *dem.tif -- e.g. an ENVI .bil input.",
    )
    parser.add_argument(
        "--channels",
        default=None,
        help="Channel-mask raster to draw under the vectors (optional).",
    )
    parser.add_argument(
        "--max-traces",
        type=int,
        default=None,
        help="Cap the number of traces drawn, for legibility on large runs.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=160,
        help="Raster resolution for the PNGs (default: 160).",
    )
    return parser


def main(argv=None) -> int:
    """Entry point for ``lsdviz3-hilltop-coupling``.

    Args:
        argv: Argument list; ``sys.argv[1:]`` when omitted.

    Returns:
        0 on success, 1 if the output set carried nothing to draw.
    """
    args = build_parser().parse_args(argv)

    outputs = load_hilltop_outputs(args.prefix)
    if args.dem:
        outputs.dem_path = Path(args.dem)
    if args.channels:
        outputs.channels_path = Path(args.channels)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    figures = []

    def add(name, fn):
        # One failing panel must not cost the others: a run may legitimately
        # lack traces, or curvature, or chi.
        try:
            figures.append((name, fn()))
        except Exception as exc:  # noqa: BLE001
            print(f"skipped {name}: {exc}")

    if outputs.dem_path is not None and outputs.any_traces() is not None:
        add("traces_on_hillshade",
            lambda: plot_traces_on_hillshade(outputs, max_traces=args.max_traces))
    if outputs.trace_table() is not None:
        add("trace_outcomes", lambda: plot_trace_outcomes(outputs))
    if outputs.dem_path is not None and outputs.ridgelines is not None:
        add("ridgelines_and_sides", lambda: plot_ridgelines_and_sides(outputs))
    if outputs.any_traces() is not None and outputs.dem_path is not None:
        # Restricted to sources within one hillslope length of a channel:
        # loops further out have no channel to reach and explain nothing.
        add("looped_trace_profiles",
            lambda: plot_looped_trace_profiles(outputs, near_channel_m=75.0))

    if outputs.hilltops is not None:
        add("length_distributions", lambda: plot_length_distributions(outputs))
        add("ridgeline_profiles", lambda: plot_ridgeline_profiles(outputs))
        add("estar_rstar",
            lambda: plot_estar_rstar(outputs.hilltops,
                                     title="Indian Creek hilltop metrics"))

    if not figures:
        print(f"nothing to draw: no recognised outputs under {args.prefix}")
        return 1

    for name, fig in figures:
        path = outdir / f"{name}.png"
        fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
        print(f"wrote {path}")

    if outputs.summary:
        print("\nrun summary:")
        for key, value in outputs.summary.items():
            print(f"  {key:32s} {value}")

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
