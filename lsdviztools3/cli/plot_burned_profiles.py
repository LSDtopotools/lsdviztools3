"""CLI: lsdviz3-burned-profiles - long profiles of burned channels.

Reads ``<write_fname>_burned_profiles.csv`` from one or more lsdtt3
``lsdtt-stream-burning`` runs with ``stream_burn_mode: burn_to_corridor`` and
draws one course (a path and the paths it joins, down to the outlet) per run,
one panel each, on a shared x axis: the input DEM under the channel, the
burned channel, the minima the profile joins, the small pits filled before
sampling, and where the channel was cut or raised.

Choose the course with --name (from the most upstream cell with that river
name; needs vector_name_field in the run), --path-rank (0 = the first,
longest path; the default) or --path-id (one run only).

Example (X = 0 / 60 / 150 m):

    lsdviz3-burned-profiles --csv c_vt2_x0_burned_profiles.csv \\
        --csv c_vt2_x60_burned_profiles.csv --csv c_vt2_x150_burned_profiles.csv \\
        --labels "X = 0 m" "X = 60 m" "X = 150 m" --name "Wistaston Brook" \\
        --out wistaston.png
"""

from __future__ import annotations

import argparse


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-burned-profiles``."""
    p = argparse.ArgumentParser(
        prog="lsdviz3-burned-profiles",
        description="Long profiles of channels burned by lsdtt-stream-burning "
                    "(burn_to_corridor), one panel per run.")
    p.add_argument("--csv", action="append", required=True,
                   help="A <write_fname>_burned_profiles.csv; repeat for one panel per run.")
    p.add_argument("--labels", nargs="+", default=None,
                   help="One panel label per --csv (default the file stems).")
    sel = p.add_mutually_exclusive_group()
    sel.add_argument("--name", default=None, help="River name to follow (vector_name_field).")
    sel.add_argument("--path-rank", type=int, default=None,
                     help="Path rank to follow (0 = the first, longest path; the default).")
    sel.add_argument("--path-id", type=int, default=None, help="Path id to follow (one --csv only).")
    p.add_argument("--x", choices=["distance", "chi"], default="distance",
                   help="x axis: distance downstream (km, default) or chi (m; chi-fit runs).")
    p.add_argument("--mark", default=None,
                   help="x,y,label: a vertical marker at the course cell nearest (x, y), "
                        "in the run's CRS (e.g. 663518,6015882,Forge Valley).")
    p.add_argument("--max-join-labels", type=int, default=6,
                   help="Name at most this many joining rivers (longest first; default 6).")
    p.add_argument("--out", required=True, help="Output PNG.")
    p.add_argument("--title", default=None)
    p.add_argument("--dpi", type=int, default=200)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from ..charts.burned_profiles import plot_burned_profiles
    if args.path_id is not None and len(args.csv) > 1:
        raise SystemExit("--path-id works with one --csv only (path ids differ between runs); "
                         "use --name or --path-rank.")
    if args.labels and len(args.labels) != len(args.csv):
        raise SystemExit("give one --labels entry per --csv.")
    mark = None
    if args.mark:
        parts = args.mark.split(",", 2)
        if len(parts) != 3:
            raise SystemExit("--mark needs x,y,label.")
        mark = (float(parts[0]), float(parts[1]), parts[2])
    plot_burned_profiles(args.csv, args.out, labels=args.labels, path_id=args.path_id,
                         path_rank=args.path_rank, name=args.name, x=args.x, mark=mark,
                         title=args.title, dpi=args.dpi, max_join_labels=args.max_join_labels)
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
