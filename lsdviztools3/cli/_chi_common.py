"""Shared pieces of the chi CLIs (lsdviz3-chi, lsdviz3-knickpoints)."""

from __future__ import annotations

import argparse
from pathlib import Path


def mn_suffix(value: str) -> str:
    """``0.45`` or ``0p450`` -> ``0p450``, the suffix lsdtt-chi-analysis writes."""
    value = str(value).strip()
    if "p" in value:
        return value
    return f"{float(value):.3f}".replace(".", "p")


def add_tributary_selection_args(parser: argparse.ArgumentParser) -> None:
    """The tributary-selection flags (see :func:`tributary_selection_from_args`)."""
    sel = parser.add_argument_group(
        "tributary selection",
        "Remove tributaries (branches of the chi fit, head to junction) using the run's "
        "<prefix>_chi_tributaries_mn_*.csv (write_chi_tributary_table). Main stems are "
        "never removed. Writes <out-prefix>_removed_tributaries.csv.")
    sel.add_argument("--remove-ksn-above", default=None,
                     help="Remove tributaries with any segment k_sn above this: a value, or "
                          "q3_plus_3iqr (the outlier value lsdtt3 suggests).")
    sel.add_argument("--min-tributary-area-m2", type=float, default=None,
                     help="Remove tributaries whose drainage area just before they join is below this (m^2).")
    sel.add_argument("--min-tributary-area-fraction", type=float, default=None,
                     help="Remove tributaries whose area where they join is below this fraction of the basin's.")
    sel.add_argument("--top-tributaries", type=int, default=None,
                     help="Keep only this many tributaries per basin, the largest by area where they join.")
    sel.add_argument("--longest-channels", type=int, default=None,
                     help="Keep only the main stem plus the N-1 longest tributaries per basin (flow "
                          "length head to junction), so 1 keeps just the main stem. Applied after the "
                          "other rules; the kept channels always form one connected network.")
    sel.add_argument("--keep-descendants", action="store_true",
                     help="Keep tributaries that drain into a removed one (default: remove them too).")
    sel.add_argument("--show-removed", choices=["ghost", "hide"], default="ghost",
                     help="Draw removed tributaries thin and grey on the profiles (ghost, default) or not at all.")


def tributary_selection_from_args(args, base: Path, mn: str, out_dir: Path, stem: str):
    """Apply the selection flags. Returns the ``removed`` frame, or ``None``
    when no selection flag was given. Writes ``<stem>_removed_tributaries.csv``."""
    filters = (args.remove_ksn_above, args.min_tributary_area_m2,
               args.min_tributary_area_fraction, args.top_tributaries, args.longest_channels)
    if all(v is None for v in filters):
        return None
    from ..io.chi import load_chi_knickpoint_thresholds, load_chi_tributaries, select_tributaries
    tpath = Path(f"{base}_tributaries_mn_{mn}.csv")
    if not tpath.exists():
        raise SystemExit(f"{tpath} not found: tributary selection needs a chi run with "
                         "write_chi_tributary_table: true (lsdtt3 with tributary tables).")
    ksn_max = args.remove_ksn_above
    thresholds = None
    if ksn_max is not None:
        try:
            ksn_max = float(ksn_max)
        except ValueError:
            thpath = Path(f"{base}_knickpoint_thresholds_mn_{mn}.csv")
            if not thpath.exists():
                raise SystemExit(f"{thpath} not found (needed for --remove-ksn-above {ksn_max}).")
            thresholds = load_chi_knickpoint_thresholds(thpath)
    try:
        _, removed = select_tributaries(
            load_chi_tributaries(tpath), ksn_max=ksn_max,
            min_junction_area_m2=args.min_tributary_area_m2,
            min_area_fraction=args.min_tributary_area_fraction,
            top_n=args.top_tributaries, longest_channels=args.longest_channels,
            remove_descendants=not args.keep_descendants,
            thresholds=thresholds)
    except ValueError as exc:
        raise SystemExit(str(exc))
    out_dir.mkdir(parents=True, exist_ok=True)
    rpath = out_dir / f"{stem}_removed_tributaries.csv"
    removed.to_csv(rpath, index=False)
    print(f"removed {len(removed)} tributaries {removed['reason'].value_counts().to_dict()} -> {rpath}")
    return removed


def add_ksn_class_basis_arg(parser: argparse.ArgumentParser) -> None:
    """``--ksn-class-basis``: which k_sn distribution a class scheme is found in."""
    parser.add_argument(
        "--ksn-class-basis", choices=["segments", "pixels"], default="segments",
        help="Which k_sn distribution a --ksn-classes scheme (natural_breaks, quantile, "
             "equal_interval) is found in: segments (each fitted segment one sample, default) or "
             "pixels (each channel node one sample, so length-weighted). The same classes colour "
             "the map and the chi profiles, and are drawn on the k_sn distribution plot.")


def ksn_classes_from_args(args, segments, nodes, removed):
    """The k_sn :class:`~lsdviztools3.colour_classes.ColourClasses` from
    ``--ksn-classes`` and friends, found on the kept tributaries in the
    ``--ksn-class-basis`` distribution; ``None`` without ``--ksn-classes``.
    ``nodes`` may be a callable that loads them (only called for pixels)."""
    from ..charts.ksn_distribution import ksn_colour_classes
    from . import _common

    kw = _common.class_kwargs(args)
    if kw["classes"] is None:
        return None
    if args.ksn_class_basis == "pixels" and callable(nodes):
        nodes = nodes()
    try:
        cc = ksn_colour_classes(kw["classes"], segments=segments, nodes=nodes,
                                basis=args.ksn_class_basis, n_classes=kw["n_classes"],
                                cmap=args.cmap, class_colors=kw["class_colors"],
                                color_mode=kw["class_color_mode"], vmin=args.ksn_vmin,
                                vmax=args.ksn_vmax, removed=removed)
    except ValueError as exc:
        raise SystemExit(str(exc))
    print(f"k_sn classes ({cc.scheme}): breaks " + ", ".join(f"{b:.4g}" for b in cc.breaks))
    return cc
