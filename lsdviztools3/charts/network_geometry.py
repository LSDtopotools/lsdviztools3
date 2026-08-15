"""Charts for the network-geometry sweep (matplotlib).

The sweep sets an Optimal Channel Network, annealed across temperature, against
stochastic-surface controls (diamond-square, red noise, white noise). The
question these figures answer is which measurable network properties actually
separate an optimised network from an arbitrary one -- and, just as usefully,
which do not.

Three views, deliberately different in kind:

:func:`plot_junction_angle_cdf`
    Empirical CDFs of junction angle. A CDF rather than a histogram because the
    families differ in *where* their distributions sit and how tight they are,
    and overlaid CDFs show both without the bin-width argument a histogram
    invites.

:func:`plot_metric_panels`
    One panel per summary metric, bars grouped and coloured by family. The
    across-family contrast is the claim; within-family variation is the
    yardstick for whether that contrast means anything.

:func:`plot_metric_vs_temperature`
    The OCN temperature series alone, with the control families as horizontal
    reference bands. This is the figure that shows whether a metric responds to
    annealing temperature at all, which is a separate question from whether it
    separates families.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core install stays lean.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np

from ..io.network_geometry import FAMILY_COLOURS, FAMILY_ORDER, family_of

_MPL_MSG = (
    "Network-geometry charts need matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def _ecdf(values: np.ndarray):
    """Empirical CDF as step coordinates.

    Returns ``(x, y)`` with ``y`` running from 1/n to 1, so the curve reaches
    the top of the axis exactly at the largest observation rather than
    stopping short of it.
    """
    x = np.sort(np.asarray(values, dtype=float))
    x = x[np.isfinite(x)]
    if x.size == 0:
        return x, x
    y = np.arange(1, x.size + 1, dtype=float) / float(x.size)
    return x, y


# --------------------------------------------------------------------------
# Junction-angle CDF
# --------------------------------------------------------------------------

def plot_junction_angle_cdf(
    angle_sets: Mapping[str, np.ndarray],
    *,
    labels: Optional[Mapping[str, str]] = None,
    colour_by_family: bool = True,
    reference_angle: Optional[float] = 72.0,
    min_count: int = 30,
    title: Optional[str] = None,
    ax=None,
):
    """Overlaid empirical CDFs of junction angle, one curve per DEM.

    Args:
        angle_sets: Mapping ``stem -> angles`` in degrees, as returned by
            :func:`lsdviztools3.io.network_geometry.load_junction_angle_sets`.
        labels: Optional ``stem -> legend label``. Defaults to the stem.
        colour_by_family: Colour each curve by generation family, so the
            families read as groups rather than as ten unrelated lines. Line
            style then distinguishes members within a family.
        reference_angle: Draw a vertical reference line here. Defaults to 72
            degrees, the optimal-network branching angle of Howard (1971) for
            equal-discharge tributaries; pass ``None`` to omit it.
        min_count: Curves built from fewer junctions than this are still drawn
            but marked in the legend, because a CDF over a handful of junctions
            is a staircase, not a distribution.
        title: Optional title.
        ax: Existing matplotlib ``Axes`` or ``None``.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.6, 5.2))
    else:
        fig = ax.figure

    # Within a family, cycle line styles so overlapping curves stay separable.
    styles = ["-", "--", "-.", ":"]
    seen_in_family: dict = {}

    for stem, values in angle_sets.items():
        x, y = _ecdf(np.asarray(values, dtype=float))
        if x.size == 0:
            continue

        fam = family_of(stem)
        idx = seen_in_family.get(fam, 0)
        seen_in_family[fam] = idx + 1

        colour = FAMILY_COLOURS.get(fam) if colour_by_family else None
        label = (labels or {}).get(stem, stem)
        if x.size < min_count:
            label = f"{label} (n={x.size}, sparse)"

        ax.plot(
            x, y,
            linestyle=styles[idx % len(styles)],
            color=colour,
            linewidth=1.6,
            label=label,
        )

    if reference_angle is not None:
        ax.axvline(
            reference_angle, color="k", linestyle=":", linewidth=1.0, zorder=0,
        )
        ax.text(
            reference_angle, 0.02, f"  {reference_angle:g}$^\\circ$ (Howard 1971)",
            fontsize=8, color="k", rotation=90, va="bottom", ha="left",
        )

    ax.set_xlabel("junction (donor opening) angle  [degrees]")
    ax.set_ylabel("cumulative fraction of junctions")
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower right", framealpha=0.9)
    if title:
        ax.set_title(title)

    fig.tight_layout()
    return fig, ax


def plot_delta_L_cdf(
    pair_sets: Mapping[str, np.ndarray],
    *,
    labels: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    ax=None,
):
    """Overlaid CDFs of lengthwise asymmetry ``delta_L``.

    The published end-members are marked: 0.19 for sub-parallel trunk streams
    and 1.00 for a tortuous network (Goren & Shelef 2024). Those are medians, so
    they are the x-positions at which a curve should cross y = 0.5, not bounds
    on the distribution.
    """
    plt = _import_mpl()

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.6, 5.2))
    else:
        fig = ax.figure

    styles = ["-", "--", "-.", ":"]
    seen_in_family: dict = {}

    for stem, values in pair_sets.items():
        x, y = _ecdf(np.asarray(values, dtype=float))
        if x.size == 0:
            continue
        fam = family_of(stem)
        idx = seen_in_family.get(fam, 0)
        seen_in_family[fam] = idx + 1
        ax.plot(
            x, y,
            linestyle=styles[idx % len(styles)],
            color=FAMILY_COLOURS.get(fam),
            linewidth=1.6,
            label=f"{(labels or {}).get(stem, stem)} (n={x.size})",
        )

    for value, note in ((0.19, "sub-parallel"), (1.00, "tortuous")):
        ax.axvline(value, color="k", linestyle=":", linewidth=1.0, zorder=0)
        ax.text(value, 0.02, f"  {note}", fontsize=8, rotation=90,
                va="bottom", ha="left")

    ax.axhline(0.5, color="k", linewidth=0.6, alpha=0.4, zorder=0)
    ax.set_xlabel(r"lengthwise asymmetry  $\Delta L$")
    ax.set_ylabel("cumulative fraction of pairs")
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower right", framealpha=0.9)
    if title:
        ax.set_title(title)

    fig.tight_layout()
    return fig, ax


def plot_edge_basin_cdf(
    size_sets: Mapping[str, np.ndarray],
    *,
    labels: Optional[Mapping[str, str]] = None,
    area_weighted: bool = True,
    title: Optional[str] = None,
    ax=None,
):
    """CDFs of edge-draining basin size, as a fraction of the tile.

    Args:
        size_sets: Mapping ``stem -> basin_fraction_of_domain`` arrays, from
            :func:`lsdviztools3.io.network_geometry.load_edge_basin_sizes`.
        area_weighted: Weight each basin by its own area, so the y-axis reads
            "fraction of the TILE that drains through basins this size or
            smaller". This is almost always the one you want: unweighted, every
            landscape looks the same because the counts are dominated by
            single-cell boundary basins with negligible area. Set False to see
            that raw count distribution.
        labels: Optional ``stem -> legend label``.
        title: Optional title.
        ax: Existing matplotlib ``Axes`` or ``None``.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.6, 5.2))
    else:
        fig = ax.figure

    styles = ["-", "--", "-.", ":"]
    seen_in_family: dict = {}

    for stem, values in size_sets.items():
        v = np.asarray(values, dtype=float)
        v = np.sort(v[np.isfinite(v)])
        if v.size == 0:
            continue

        if area_weighted:
            total = v.sum()
            y = np.cumsum(v) / total if total > 0 else np.zeros_like(v)
        else:
            y = np.arange(1, v.size + 1, dtype=float) / v.size

        fam = family_of(stem)
        idx = seen_in_family.get(fam, 0)
        seen_in_family[fam] = idx + 1

        ax.plot(v, y, linestyle=styles[idx % len(styles)],
                color=FAMILY_COLOURS.get(fam), linewidth=1.6,
                label=f"{(labels or {}).get(stem, stem)} (n={v.size})")

    ax.set_xscale("log")
    ax.set_xlabel("edge-draining basin size  [fraction of tile]")
    ax.set_ylabel("fraction of tile draining through basins this size or smaller"
                  if area_weighted else "cumulative fraction of basins")
    ax.set_ylim(0.0, 1.0)
    ax.axhline(0.5, color="k", linewidth=0.6, alpha=0.4, zorder=0)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
    if title:
        ax.set_title(title)

    fig.tight_layout()
    return fig, ax


# --------------------------------------------------------------------------
# Summary panels
# --------------------------------------------------------------------------

#: Default panels: ``(column, axis label, reference value or None)``.
#:
#: Hack's exponent h is deliberately NOT here. lsdtt3 writes the two rasters it
#: needs (longest upslope length via ``write_distance_to_farthest_ridge_raster``
#: and drainage area) but does not fit the exponent, and fitting it in the
#: plotting layer would put a topographic computation in lsdviztools3, which is
#: the one thing this package does not do. Add ``("hack_h", ...)`` here once
#: lsdtt3 writes a Hack summary; ``load_hack`` is already waiting for it.
DEFAULT_PANELS = (
    ("median_delta_L", r"lengthwise asymmetry  $\Delta L$", 0.19),
    ("beta", r"drainage-area exponent  $\beta$", 0.45),
    ("median_deg", "median junction angle  [deg]", 72.0),
    ("ar_c_equivalent", "Tokunaga $c$ (from area ratios)", None),
    ("internal_divide_fraction", "fraction of divides that are internal", None),
    ("area_weighted_median_fraction",
     "area-weighted median\nedge-basin size [frac. of tile]", None),
)


def plot_metric_panels(
    table,
    *,
    panels: Sequence = DEFAULT_PANELS,
    labels: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    figsize=None,
):
    """One panel per metric, bars grouped and coloured by generation family.

    Args:
        table: Tidy table from
            :func:`lsdviztools3.io.network_geometry.collect_network_sweep`.
        panels: Sequence of ``(column, ylabel, reference)``. Columns absent
            from the table are skipped with their panel, rather than drawn
            empty.
        labels: Optional ``stem -> tick label``.
        title: Optional suptitle.
        figsize: Defaults to a width that scales with the panel count.

    Returns:
        ``(fig, axes)``.
    """
    plt = _import_mpl()
    import pandas as pd

    df = table.copy()
    if "family" in df.columns:
        # Sort by temperature within family, NOT by stem: a string sort puts
        # "ocn_T001" (T=0.01) before "ocn_T00" (T=0) and silently scrambles the
        # temperature series across the x-axis. Non-OCN rows have NaN
        # temperature, which sorts last within its family and leaves the stem
        # order intact there.
        sort_key = "temperature" if "temperature" in df.columns else "stem"
        df = df.sort_values(["family", sort_key, "stem"], kind="stable",
                            na_position="last")

    usable = [p for p in panels if p[0] in df.columns and df[p[0]].notna().any()]
    if not usable:
        raise ValueError(
            "none of the requested metric columns are present and non-empty; "
            f"table has {list(df.columns)}"
        )

    n = len(usable)
    if figsize is None:
        figsize = (3.6 * n, 5.6)
    fig, axes = plt.subplots(1, n, figsize=figsize)
    if n == 1:
        axes = [axes]

    stems = df["stem"].tolist()
    ticks = [(labels or {}).get(s, s) for s in stems]
    colours = [FAMILY_COLOURS.get(f, "#999999") for f in df["family"]]
    x = np.arange(len(stems))

    for ax, (col, ylabel, reference) in zip(axes, usable):
        values = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        ax.bar(x, values, color=colours)

        if reference is not None:
            ax.axhline(reference, color="crimson", linestyle=":", linewidth=1.2)
            ax.text(
                0.01, reference, f" {reference:g}", color="crimson", fontsize=8,
                transform=ax.get_yaxis_transform(), va="bottom", ha="left",
            )

        ax.set_xticks(x)
        ax.set_xticklabels(ticks, rotation=90, fontsize=7.5)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.3)

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=FAMILY_COLOURS[f])
        for f in FAMILY_ORDER if f in set(df["family"])
    ]
    names = [f for f in FAMILY_ORDER if f in set(df["family"])]
    fig.legend(handles, names, loc="upper center", ncol=len(names), fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 0.99))

    if title:
        fig.suptitle(title, fontweight="bold", y=1.02)

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig, axes


def plot_metric_vs_temperature(
    table,
    *,
    metric: str = "median_delta_L",
    ylabel: Optional[str] = None,
    show_control_bands: bool = True,
    title: Optional[str] = None,
    ax=None,
):
    """The OCN temperature series, against the control families as bands.

    Separating this from :func:`plot_metric_panels` is the point: a metric can
    separate generation families cleanly and still be flat in temperature, and
    conflating the two questions in one figure is how that gets missed. The
    control bands span each stochastic family's min-max, so a temperature curve
    lying inside a band is not distinguishing the OCN from that control.

    Args:
        table: Tidy table from
            :func:`lsdviztools3.io.network_geometry.collect_network_sweep`.
        metric: Column to plot.
        ylabel: Y-axis label; defaults to the metric name.
        show_control_bands: Draw the non-OCN families as shaded bands.
        title: Optional title.
        ax: Existing matplotlib ``Axes`` or ``None``.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    import pandas as pd

    if metric not in table.columns:
        raise KeyError(f"{metric!r} not in table; have {list(table.columns)}")

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
    else:
        fig = ax.figure

    df = table.copy()
    df[metric] = pd.to_numeric(df[metric], errors="coerce")

    ocn = df[df["family"] == "Optimal Channel Network"].dropna(subset=["temperature"])
    ocn = ocn.sort_values("temperature")
    if not ocn.empty:
        ax.plot(
            ocn["temperature"], ocn[metric],
            marker="o", color=FAMILY_COLOURS["Optimal Channel Network"],
            linewidth=1.8, label="Optimal Channel Network",
        )

    if show_control_bands:
        for fam in FAMILY_ORDER:
            if fam == "Optimal Channel Network":
                continue
            sub = df[df["family"] == fam][metric].dropna()
            if sub.empty:
                continue
            lo, hi = float(sub.min()), float(sub.max())
            colour = FAMILY_COLOURS[fam]
            if np.isclose(lo, hi):
                ax.axhline(lo, color=colour, linestyle="--", linewidth=1.4,
                           label=fam)
            else:
                ax.axhspan(lo, hi, color=colour, alpha=0.18, label=fam)

    ax.set_xlabel("OCN annealing temperature scale")
    ax.set_ylabel(ylabel or metric)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, framealpha=0.9)
    if title:
        ax.set_title(title)

    fig.tight_layout()
    return fig, ax
