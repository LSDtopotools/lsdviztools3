"""k_sn by lithologic unit (matplotlib).

Charts of the ``lsdtt-chi-analysis`` lithology outputs (read with
:mod:`lsdviztools3.io.lithology`):

* :func:`plot_ksn_histogram_by_lithology` — the common-bin histogram lsdtt3
  wrote, per unit, by node count or by channel length; overlaid or as small
  multiples. The open-ended outermost bins (``-inf``/``inf`` edges, values
  beyond the pooled 1st-99th percentile) are reported as overflow text, not
  drawn as infinite bars.
* :func:`plot_ksn_box_by_lithology` — per-unit box plots of the per-node k_sn
  values lsdtt3 wrote (or, with no node table, boxes built directly from the
  summary's median/quartiles), annotated with the summary's sample sizes;
  low-sample units are marked.
* :func:`plot_reach_ksn_by_lithology` — box plots of per-reach k_sn (OLS) for
  qualifying reaches, with the summary's reach-median confidence interval.
* :func:`plot_ksn_by_lithology_and_basin` — unit x basin dot plot.

Unit colours come from :func:`chart_unit_colours` (keyed by unit name, shared
with the geologic map); every chart takes ``unit_colours={name: colour}`` to
override them.

Every number printed on a figure (sample sizes, medians, CIs) is read from the
lsdtt3 summary CSV; nothing is recomputed here. Normalised histograms divide a
bin by the unit total that lsdtt3 wrote in the summary.

Node values are each node's **fitted segment** k_sn, so neighbouring nodes on
one segment share a value: they are node-weighted, not independent samples.
Figures say so.

k_sn is in m^(2 theta) with A_0 = 1 m^2, theta = m/n from the files.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from ..io.lithology import (
    ALL_UNITS_CODE, CODE_COLUMN, HISTOGRAM_COLUMNS, NAME_COLUMN, NODE_COLUMNS,
    REACH_COLUMNS, SUMMARY_COLUMNS, UNIT_BASIN_COLUMNS,
)

_MPL_MSG = (
    "Lithology k_sn charts need matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Caption fragment stating what a node value is.
NODE_WEIGHTED_NOTE = (
    "node-weighted: each node carries its fitted segment k$_{sn}$, "
    "so nodes are not independent samples"
)


def _caption(title: Optional[str], text: str, width: int = 95) -> str:
    """Title line (if any) over a wrapped caption."""
    import textwrap

    body = textwrap.fill(text, width=width)
    return f"{title}\n{body}" if title else body


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def ksn_label(m_over_n: Optional[float]) -> str:
    """Axis label: ``k_sn (m^(2 theta))`` with the exponent from m/n."""
    if m_over_n is None or not np.isfinite(m_over_n):
        return "k$_{sn}$ (m$^{2\\theta}$, A$_0$ = 1 m$^2$)"
    return f"k$_{{sn}}$ (m$^{{{2.0 * m_over_n:g}}}$, A$_0$ = 1 m$^2$, $\\theta$ = {m_over_n:g})"


def _mn_from(*tables) -> Optional[float]:
    col = SUMMARY_COLUMNS["m_over_n"]
    for t in tables:
        if t is not None and col in t.columns:
            v = pd.to_numeric(t[col], errors="coerce").dropna()
            if len(v):
                return float(v.iloc[0])
    return None


def unit_rows(summary: pd.DataFrame, *, order: str = "median",
              units: Optional[Sequence[int]] = None,
              include_empty: bool = False) -> pd.DataFrame:
    """The summary's per-unit rows (no ``all``/empty ``unmapped``), ordered.

    Args:
        summary: Per-unit summary table.
        order: ``"median"`` (ascending summary node median),
            ``"reach_median"`` (ascending summary reach median; falls back to
            the node median if absent), ``"code"``, ``"name"`` or ``"length"``
            (descending channel length).
        units: Optional subset of codes to keep.
        include_empty: Keep units with zero nodes.
    """
    s = summary[summary[CODE_COLUMN] != ALL_UNITS_CODE].copy()
    if not include_empty:
        s = s[pd.to_numeric(s[SUMMARY_COLUMNS["n_nodes"]], errors="coerce") > 0]
    if units is not None:
        s = s[s[CODE_COLUMN].isin(list(units))]
    if order == "reach_median" and SUMMARY_COLUMNS["reach_median"] not in s.columns:
        order = "median"
    key = {"median": (SUMMARY_COLUMNS["median"], True),
           "reach_median": (SUMMARY_COLUMNS["reach_median"], True),
           "code": (CODE_COLUMN, True),
           "name": (NAME_COLUMN, True),
           "length": (SUMMARY_COLUMNS["channel_length_m"], False)}
    if order not in key:
        raise ValueError(f"order must be one of {sorted(key)}, got {order!r}.")
    col, asc = key[order]
    s = s.sort_values(col, ascending=asc, kind="stable")
    # Reserved categories (mixed, arms disagree, unmapped) always after the units.
    from ..io.lithology import reserved_category
    from ..render.lithology_colours import RESERVED_ORDER

    rank = [(-1 if reserved_category(n) is None else RESERVED_ORDER.index(reserved_category(n)))
            for n in s[NAME_COLUMN]]
    s = s.assign(_reserved_rank=rank).sort_values("_reserved_rank", kind="stable")
    return s.drop(columns="_reserved_rank").reset_index(drop=True)


def _is_low(row) -> bool:
    col = SUMMARY_COLUMNS["low_sample"]
    try:
        return bool(int(row.get(col, 0)))
    except (TypeError, ValueError):
        return False


def sample_size_text(row, *, reaches: bool = False) -> str:
    """``"n=5484 nodes, 917 seg, 81 basins"`` from one summary row."""
    def _i(key):
        v = row.get(SUMMARY_COLUMNS[key])
        return None if v is None or pd.isna(v) else int(v)

    if reaches:
        q, n = _i("n_reaches_qualifying"), _i("n_reaches")
        parts = [f"{q} of {n} reaches" if q is not None and n is not None
                 else f"{n if n is not None else '?'} reaches"]
    else:
        parts = [f"n={_i('n_nodes')} nodes"]
        if _i("n_segments") is not None:
            parts.append(f"{_i('n_segments')} seg")
    if _i("n_basins") is not None:
        parts.append(f"{_i('n_basins')} basins")
    text = ", ".join(parts)
    if not _is_low(row):
        return text
    # Since the reach-level outputs, lsdtt3's low_sample means too few
    # QUALIFYING reaches (its CI is then NaN); say so when that column exists.
    if SUMMARY_COLUMNS["n_reaches_qualifying"] in row.index:
        return text + " (low sample: few qualifying reaches)"
    return text + " (low sample)"


def _truthy(series: pd.Series) -> np.ndarray:
    """1/1.0/true/yes -> True (lsdtt3 writes 0/1 flags)."""
    num = pd.to_numeric(series, errors="coerce")
    txt = series.astype(str).str.strip().str.lower().isin(("true", "yes"))
    return ((num == 1) | txt).to_numpy()


def buffer_note(summary: Optional[pd.DataFrame]) -> str:
    """``"contact buffer 120 m"`` from the summary (empty if absent)."""
    if summary is None:
        return ""
    v = _mn_from_col(summary, SUMMARY_COLUMNS["contact_buffer_m"])
    return "" if v is None else f"contact buffer {v:g} m"


def _mn_from_col(df: pd.DataFrame, col: str) -> Optional[float]:
    if col not in df.columns:
        return None
    v = pd.to_numeric(df[col], errors="coerce").dropna()
    return float(v.iloc[0]) if len(v) else None


def chart_unit_colours(summary: pd.DataFrame,
                       unit_colours: Optional[Mapping[str, str]] = None) -> dict:
    """``{name: hex}`` for the summary's units, by the shared name-keyed rule.

    Uses :func:`render.lithology_colours.lithology_palette` on the summary's
    unit names, so the colours match the geologic map of the same run.
    """
    from ..render.lithology_colours import lithology_palette

    return lithology_palette(summary[NAME_COLUMN].astype(str), unit_colours)


def _unit_colours(summary, rows, unit_colours=None) -> dict:
    """``{code: hex}`` for the plotted rows (looked up by name)."""
    pal = chart_unit_colours(summary, unit_colours)
    return {c: pal.get(str(n), "#888889")
            for c, n in zip(rows[CODE_COLUMN], rows[NAME_COLUMN].astype(str))}


# --------------------------------------------------------------------------
# Histogram
# --------------------------------------------------------------------------

def _split_bins(h: pd.DataFrame):
    """(inner, below, above) rows of one unit's histogram."""
    lo = h[HISTOGRAM_COLUMNS["bin_lower"]].to_numpy(float)
    hi = h[HISTOGRAM_COLUMNS["bin_upper"]].to_numpy(float)
    inner = np.isfinite(lo) & np.isfinite(hi)
    return h[inner], h[~np.isfinite(lo)], h[~np.isfinite(hi)]


def histogram_overflow(h_unit: pd.DataFrame, weight_col: str) -> dict:
    """Overflow of one unit's histogram: ``{"below": (edge, w), "above": (edge, w)}``.

    ``edge`` is the finite edge of the open bin and ``w`` its weight as written.
    """
    _inner, below, above = _split_bins(h_unit)
    out = {}
    if len(below):
        out["below"] = (float(below[HISTOGRAM_COLUMNS["bin_upper"]].iloc[0]),
                        float(below[weight_col].sum()))
    if len(above):
        out["above"] = (float(above[HISTOGRAM_COLUMNS["bin_lower"]].iloc[0]),
                        float(above[weight_col].sum()))
    return out


def plot_ksn_histogram_by_lithology(
    histogram: pd.DataFrame,
    summary: pd.DataFrame,
    *,
    weight: str = "length",
    normalise: bool = True,
    layout: str = "auto",
    order: str = "median",
    units: Optional[Sequence[int]] = None,
    max_overlay: int = 6,
    ncols: int = 4,
    unit_colours: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
):
    """k_sn histograms by unit from lsdtt3's common-bin histogram CSV.

    Args:
        histogram: ``*_chi_ksn_histogram_by_lithology_mn_*.csv`` table.
        summary: ``*_chi_ksn_by_lithology_mn_*.csv`` table (names, sample
            sizes, unit totals for normalising, ordering).
        weight: ``"length"`` (channel length in the bin) or ``"count"`` (nodes).
        normalise: Divide each bin by the unit's total (length or node count
            from the summary), so units of different extent compare by shape.
        layout: ``"overlay"`` (one axis), ``"facets"`` (one panel per unit) or
            ``"auto"`` (overlay up to ``max_overlay`` units).
        order: Unit order (see :func:`unit_rows`).
        units: Optional subset of unit codes.
        max_overlay: Unit count above which ``auto`` switches to facets.
        ncols: Facet columns.
        title: Figure title.

    Returns:
        ``(fig, axes)`` — ``axes`` is a list.
    """
    plt = _import_mpl()
    if weight not in ("length", "count"):
        raise ValueError(f"weight must be 'length' or 'count', got {weight!r}.")
    wcol = HISTOGRAM_COLUMNS["channel_length_m" if weight == "length" else "n_nodes"]
    tcol = SUMMARY_COLUMNS["channel_length_m" if weight == "length" else "n_nodes"]

    rows = unit_rows(summary, order=order, units=units)
    if rows.empty:
        raise ValueError("No units with nodes in the summary.")
    if layout == "auto":
        layout = "overlay" if len(rows) <= max_overlay else "facets"
    if layout not in ("overlay", "facets"):
        raise ValueError(f"layout must be overlay/facets/auto, got {layout!r}.")
    colours = _unit_colours(summary, rows, unit_colours)

    if layout == "overlay":
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        axes = [ax] * len(rows)
    else:
        nrows = int(np.ceil(len(rows) / ncols))
        fig, grid = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.1 * nrows),
                                 sharex=True, squeeze=False)
        axes = list(grid.ravel())
        for extra in axes[len(rows):]:
            extra.set_visible(False)

    if weight == "length":
        ylab = "fraction of unit's channel length" if normalise else "channel length (km)"
    else:
        ylab = "fraction of unit's nodes" if normalise else "nodes"
    xmin, xmax = np.inf, -np.inf

    for i, (_, row) in enumerate(rows.iterrows()):
        code = row[CODE_COLUMN]
        h = histogram[histogram[CODE_COLUMN] == code].sort_values(
            HISTOGRAM_COLUMNS["bin_index"])
        inner, _below, _above = _split_bins(h)
        if inner.empty:
            continue
        edges = np.append(inner[HISTOGRAM_COLUMNS["bin_lower"]].to_numpy(float),
                          float(inner[HISTOGRAM_COLUMNS["bin_upper"]].iloc[-1]))
        vals = inner[wcol].to_numpy(float)
        total = float(row[tcol]) if normalise else 1.0
        scale = 1.0 / total if (normalise and total > 0) else (1e-3 if weight == "length" else 1.0)
        xmin, xmax = min(xmin, edges[0]), max(xmax, edges[-1])

        ov = histogram_overflow(h, wcol)
        ovtxt = []
        if "below" in ov and ov["below"][1] > 0:
            ovtxt.append(f"<{ov['below'][0]:.3g}: {_fmt_w(ov['below'][1] * scale, normalise)}")
        if "above" in ov and ov["above"][1] > 0:
            ovtxt.append(f">{ov['above'][0]:.3g}: {_fmt_w(ov['above'][1] * scale, normalise)}")

        ax = axes[i]
        low = _is_low(row)
        label = f"{row[NAME_COLUMN]} ({sample_size_text(row)})"
        ax.stairs(vals * scale, edges, color=colours[code], linewidth=1.6,
                  linestyle="--" if low else "-",
                  label=label + (f"; outside: {', '.join(ovtxt)}" if ovtxt and layout == "overlay" else ""))
        if layout == "facets":
            ax.set_title(str(row[NAME_COLUMN]), fontsize=8)
            note = sample_size_text(row)
            if ovtxt:
                note += "\noutside: " + ", ".join(ovtxt)
            ax.text(0.98, 0.95, note, transform=ax.transAxes, ha="right", va="top",
                    fontsize=6, color="0.4" if low else "0.1")
            ax.tick_params(labelsize=7)

    mn = _mn_from(histogram, summary)
    what = "length-weighted" if weight == "length" else "node counts"
    sub = f"{what}; {NODE_WEIGHTED_NOTE}; common bins; values beyond the pooled 1st-99th percentile listed as 'outside'"
    if layout == "overlay":
        ax = axes[0]
        ax.set_xlim(xmin, xmax)
        ax.set_xlabel(ksn_label(mn))
        ax.set_ylabel(ylab)
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(alpha=0.3)
    else:
        for ax in axes[: len(rows)]:
            ax.set_xlim(xmin, xmax)
        fig.supxlabel(ksn_label(mn), fontsize=9)
        fig.supylabel(ylab, fontsize=9)
    fig.suptitle(_caption(title, sub, 110 if layout == "facets" else 95), fontsize=8)
    fig.tight_layout()
    return fig, axes[: len(rows)] if layout == "facets" else [axes[0]]


def _fmt_w(v: float, normalise: bool) -> str:
    return f"{100 * v:.1f}%" if normalise else f"{v:.3g}"


# --------------------------------------------------------------------------
# Box plots
# --------------------------------------------------------------------------

def _summary_box_stats(row, *, whiskers: str = "quartiles") -> dict:
    """matplotlib ``bxp`` stats built only from the lsdtt3 summary row."""
    med = float(row[SUMMARY_COLUMNS["median"]])
    q1 = float(row[SUMMARY_COLUMNS["q25"]])
    q3 = float(row[SUMMARY_COLUMNS["q75"]])
    if whiskers == "minmax":
        lo, hi = float(row[SUMMARY_COLUMNS["min"]]), float(row[SUMMARY_COLUMNS["max"]])
    else:
        lo, hi = q1, q3
    return {"med": med, "q1": q1, "q3": q3, "whislo": lo, "whishi": hi,
            "fliers": [], "label": str(row[NAME_COLUMN])}


def _style_low(artists: dict, flags: Sequence[bool]) -> None:
    for box, low in zip(artists.get("boxes", []), flags):
        if low:
            box.set_hatch("///")
            box.set_alpha(0.55)


def _orient(n: int, orientation: str) -> bool:
    """True if horizontal (units on the y axis)."""
    if orientation == "auto":
        return n > 6
    return orientation == "horizontal"


_LOW_QUALIFYING = " (low sample: few qualifying reaches)"


def _annotate_units(ax, rows, positions, horizontal, texts) -> str:
    """Write the per-unit sample-size text beside each box.

    Returns a note for the caption: vertical tick labels carry only
    "(low sample)", since the long form is wider than a box's slot and
    collides with its neighbours, so the caption explains it once.
    """
    for pos, (_, row), txt in zip(positions, rows.iterrows(), texts):
        colour = "0.45" if _is_low(row) else "0.15"
        if horizontal:
            ax.annotate(txt, xy=(1.0, pos), xycoords=("axes fraction", "data"),
                        xytext=(4, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=6.5, color=colour)
    if not horizontal:
        # Vertical boxes: fold the sample sizes into the tick labels (rotated
        # text above the axes collides with the title).
        shortened = any(_LOW_QUALIFYING in txt for txt in texts)
        labels = [f"{row[NAME_COLUMN]}\n"
                  + txt.replace(_LOW_QUALIFYING, " (low sample)").replace(", ", "\n")
                       .replace(" (low", "\n(low")
                  for (_, row), txt in zip(rows.iterrows(), texts)]
        ax.set_xticks(positions, labels, rotation=0, ha="center", fontsize=7)
        for tick, (_, row) in zip(ax.get_xticklabels(), rows.iterrows()):
            if _is_low(row):
                tick.set_color("0.45")
        if shortened:
            return "grey, hatched = low sample (few qualifying reaches)"
    return ""


def _pooled_line(ax, summary, horizontal, column, label):
    allrow = summary[summary[CODE_COLUMN] == ALL_UNITS_CODE]
    if allrow.empty or column not in allrow.columns:
        return
    v = float(allrow[column].iloc[0])
    if not np.isfinite(v):
        return
    (ax.axvline if horizontal else ax.axhline)(
        v, color="0.3", linestyle=":", linewidth=1.0, label=label, zorder=0)


def plot_ksn_box_by_lithology(
    summary: pd.DataFrame,
    nodes: Optional[pd.DataFrame] = None,
    *,
    source: str = "auto",
    order: str = "median",
    units: Optional[Sequence[int]] = None,
    orientation: str = "auto",
    show_pooled: bool = True,
    unit_colours: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    ax=None,
):
    """Per-unit k_sn box plots, annotated with lsdtt3's sample sizes.

    Args:
        summary: Per-unit summary table (names, n_nodes/n_segments/n_basins,
            low_sample, median/quartiles, pooled ``all`` row).
        nodes: Per-node table. Needed for ``source="nodes"``.
        source: ``"nodes"`` — matplotlib box plots of the per-node values
            lsdtt3 wrote (box = IQR, whiskers 1.5 IQR, no fliers);
            ``"summary"`` — boxes drawn from the summary's median and
            quartiles only (no whiskers); ``"auto"`` = nodes if given.
        order: Unit order (see :func:`unit_rows`).
        units: Optional subset of unit codes.
        orientation: ``"vertical"``, ``"horizontal"`` or ``"auto"`` (horizontal
            above 6 units, so long names stay readable).
        show_pooled: Draw the summary's all-units median as a dotted line.
        title: Axis title.
        ax: Existing axes or ``None``.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    if source == "auto":
        source = "nodes" if nodes is not None else "summary"
    if source not in ("nodes", "summary"):
        raise ValueError(f"source must be nodes/summary/auto, got {source!r}.")
    if source == "nodes" and nodes is None:
        raise ValueError("source='nodes' needs the per-node table.")

    rows = unit_rows(summary, order=order, units=units)
    if rows.empty:
        raise ValueError("No units with nodes in the summary.")
    horizontal = _orient(len(rows), orientation)
    if ax is None:
        size = (7.5, 0.32 * len(rows) + 1.8) if horizontal else (1.0 * len(rows) + 3.0, 5.0)
        fig, ax = plt.subplots(figsize=size)
    else:
        fig = ax.figure

    positions = np.arange(1, len(rows) + 1)
    colours = _unit_colours(summary, rows, unit_colours)
    lows = [_is_low(r) for _, r in rows.iterrows()]
    common = dict(positions=positions, widths=0.6,
                  orientation="horizontal" if horizontal else "vertical",
                  patch_artist=True, showfliers=False,
                  medianprops={"color": "black", "linewidth": 1.4})
    if source == "nodes":
        kcol = NODE_COLUMNS["k_sn"]
        data = [pd.to_numeric(nodes.loc[nodes[CODE_COLUMN] == c, kcol],
                              errors="coerce").dropna().to_numpy()
                for c in rows[CODE_COLUMN]]
        art = ax.boxplot(data, **common)
        how = f"box = IQR, whiskers 1.5 IQR, no outliers; {NODE_WEIGHTED_NOTE}"
    else:
        stats = [_summary_box_stats(r) for _, r in rows.iterrows()]
        art = ax.bxp(stats, **common)
        how = "box = lsdtt3 summary quartiles, line = summary median"
    for box, code in zip(art["boxes"], rows[CODE_COLUMN]):
        box.set_facecolor(colours[code])
        box.set_alpha(0.75)
    _style_low(art, lows)

    names = rows[NAME_COLUMN].astype(str).tolist()
    mn = _mn_from(summary, nodes)
    if horizontal:
        ax.set_yticks(positions, names, fontsize=7)
        ax.set_xlabel(ksn_label(mn))
        ax.grid(axis="x", alpha=0.3)
    else:
        ax.set_xticks(positions, names, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel(ksn_label(mn))
        ax.grid(axis="y", alpha=0.3)
    low_note = _annotate_units(ax, rows, positions, horizontal,
                               [sample_size_text(r) for _, r in rows.iterrows()])
    if low_note:
        how += f"; {low_note}"
    if show_pooled:
        _pooled_line(ax, summary, horizontal, SUMMARY_COLUMNS["median"],
                     "all units median (lsdtt3 summary)")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7, loc="lower right" if horizontal else "upper left")
    ax.set_title(_caption(title, how, 80 if horizontal else 70), fontsize=8)
    fig.tight_layout()
    return fig, ax


def plot_reach_ksn_by_lithology(
    reaches: pd.DataFrame,
    summary: pd.DataFrame,
    *,
    qualifying_only: bool = True,
    order: str = "reach_median",
    units: Optional[Sequence[int]] = None,
    orientation: str = "auto",
    unit_colours: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    ax=None,
):
    """Per-unit box plots of per-reach k_sn (OLS), with the reach-median CI.

    Boxes are drawn from the reach table's ``k_sn_ols`` values (qualifying
    reaches only, by default). Where the summary carries the reach-median and
    its CI columns, the median is drawn as a black diamond with the CI as an
    error bar beside each box; counts come from the summary's
    ``n_reaches``/``n_reaches_qualifying``.
    """
    plt = _import_mpl()
    kcol = REACH_COLUMNS["k_sn"]
    data_tbl = reaches
    if qualifying_only and REACH_COLUMNS["qualifies"] in reaches.columns:
        data_tbl = reaches[_truthy(reaches[REACH_COLUMNS["qualifies"]])]

    rows = unit_rows(summary, order=order, units=units)
    rows = rows[rows[CODE_COLUMN].isin(data_tbl[CODE_COLUMN].unique())].reset_index(drop=True)
    if rows.empty:
        raise ValueError("No units with reaches to plot.")
    horizontal = _orient(len(rows), orientation)
    if ax is None:
        size = (7.5, 0.32 * len(rows) + 1.8) if horizontal else (1.0 * len(rows) + 3.0, 5.0)
        fig, ax = plt.subplots(figsize=size)
    else:
        fig = ax.figure

    positions = np.arange(1, len(rows) + 1)
    data = [pd.to_numeric(data_tbl.loc[data_tbl[CODE_COLUMN] == c, kcol],
                          errors="coerce").dropna().to_numpy()
            for c in rows[CODE_COLUMN]]
    colours = _unit_colours(summary, rows, unit_colours)
    art = ax.boxplot(data, positions=positions, widths=0.5,
                     orientation="horizontal" if horizontal else "vertical",
                     patch_artist=True, showfliers=False,
                     medianprops={"color": "black", "linewidth": 1.2})
    for box, code in zip(art["boxes"], rows[CODE_COLUMN]):
        box.set_facecolor(colours[code])
        box.set_alpha(0.6)
    _style_low(art, [_is_low(r) for _, r in rows.iterrows()])

    med_c, lo_c, hi_c = (SUMMARY_COLUMNS[k] for k in ("reach_median", "reach_ci_low", "reach_ci_high"))
    has_ci = all(c in rows.columns for c in (med_c, lo_c, hi_c))
    if has_ci:
        med = pd.to_numeric(rows[med_c], errors="coerce").to_numpy(float)
        lo = pd.to_numeric(rows[lo_c], errors="coerce").to_numpy(float)
        hi = pd.to_numeric(rows[hi_c], errors="coerce").to_numpy(float)
        ok = np.isfinite(med) & np.isfinite(lo) & np.isfinite(hi)
        off = positions[ok] + 0.38
        err = np.vstack([med[ok] - lo[ok], hi[ok] - med[ok]])
        kw = dict(fmt="D", color="black", markersize=3.5, capsize=2, linewidth=1.0,
                  label="reach median and CI (lsdtt3 summary)")
        if horizontal:
            ax.errorbar(med[ok], off, xerr=err, **kw)
        else:
            ax.errorbar(off, med[ok], yerr=err, **kw)
        ax.legend(fontsize=7, loc="lower right" if horizontal else "upper left")

    names = rows[NAME_COLUMN].astype(str).tolist()
    mn = _mn_from(summary, reaches)
    if horizontal:
        ax.set_yticks(positions, names, fontsize=7)
        ax.set_xlabel("reach " + ksn_label(mn))
        ax.grid(axis="x", alpha=0.3)
    else:
        ax.set_xticks(positions, names, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel("reach " + ksn_label(mn))
        ax.grid(axis="y", alpha=0.3)
    has_counts = SUMMARY_COLUMNS["n_reaches"] in rows.columns
    texts = [sample_size_text(r, reaches=True) if has_counts else
             f"{len(d)} reaches" for (_, r), d in zip(rows.iterrows(), data)]
    low_note = _annotate_units(ax, rows, positions, horizontal, texts)
    which = "qualifying reaches" if qualifying_only else "all reaches"
    note = f"per-reach OLS k$_{{sn}}$, {which}; box = IQR, whiskers 1.5 IQR"
    if has_ci:
        nb = _mn_from_col(summary, SUMMARY_COLUMNS["n_bootstrap"])
        ccol = SUMMARY_COLUMNS["bootstrap_cluster"]
        cl = (summary[ccol].dropna().astype(str) if ccol in summary.columns
              else pd.Series([], dtype=str))
        how = f"{int(nb)} bootstrap replicates" if nb else "bootstrap"
        if len(cl):
            how += f", resampled by {cl.iloc[0]}"
        note += f"; diamond = reach median, bar = its CI (lsdtt3 summary: {how})"
    if buffer_note(summary):
        note += f"; {buffer_note(summary)}"
    if low_note:
        note += f"; {low_note}"
    ax.set_title(_caption(title, note, 80), fontsize=8)
    fig.tight_layout()
    return fig, ax


def plot_ksn_by_lithology_and_basin(
    unit_basin: pd.DataFrame,
    summary: pd.DataFrame,
    *,
    order: str = "reach_median",
    units: Optional[Sequence[int]] = None,
    unit_colours: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    ax=None,
):
    """Unit x basin dot plot of the per-basin REACH k_sn median.

    Each dot is one basin's ``reach_k_sn_median`` in one unit (from the unit x
    basin CSV); dot area scales with that basin-unit's number of qualifying
    reaches; low-sample combinations are hollow. The black bar is the unit's
    reach median from the summary. Basin-units with no reach median are not
    drawn.
    """
    plt = _import_mpl()
    rows = unit_rows(summary, order=order, units=units)
    rows = rows[rows[CODE_COLUMN].isin(unit_basin[CODE_COLUMN].unique())].reset_index(drop=True)
    if rows.empty:
        raise ValueError("No units in the unit x basin table.")
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.5, 0.32 * len(rows) + 1.8))
    else:
        fig = ax.figure
    colours = _unit_colours(summary, rows, unit_colours)
    mcol, lcol = UNIT_BASIN_COLUMNS["median"], UNIT_BASIN_COLUMNS["low_sample"]
    ncol = next((UNIT_BASIN_COLUMNS[k] for k in ("n_reaches_qualifying", "n_nodes")
                 if UNIT_BASIN_COLUMNS[k] in unit_basin.columns), None)
    bar_col = (SUMMARY_COLUMNS["reach_median"] if SUMMARY_COLUMNS["reach_median"] in rows.columns
               else SUMMARY_COLUMNS["median"])
    positions = np.arange(1, len(rows) + 1)
    n_drawn = []
    for pos, (_, row) in zip(positions, rows.iterrows()):
        d = unit_basin[unit_basin[CODE_COLUMN] == row[CODE_COLUMN]]
        x = pd.to_numeric(d[mcol], errors="coerce").to_numpy(float)
        n = (pd.to_numeric(d[ncol], errors="coerce").fillna(0).to_numpy(float)
             if ncol else np.ones_like(x))
        low = _truthy(d[lcol]) if lcol in d.columns else np.zeros_like(x, dtype=bool)
        ok = np.isfinite(x)
        x, n, low = x[ok], n[ok], low[ok]
        n_drawn.append(int(ok.sum()))
        sizes = 6.0 + 6.0 * np.sqrt(n)
        c = colours[row[CODE_COLUMN]]
        ax.scatter(x[~low], np.full((~low).sum(), pos), s=sizes[~low], color=c,
                   alpha=0.75, edgecolors="none")
        ax.scatter(x[low], np.full(low.sum(), pos), s=sizes[low], facecolors="none",
                   edgecolors=c, linewidths=0.8)
        med = pd.to_numeric(pd.Series([row.get(bar_col)]), errors="coerce").iloc[0]
        if np.isfinite(med):
            ax.plot([float(med)] * 2, [pos - 0.35, pos + 0.35], color="black", linewidth=1.6)
    ax.set_yticks(positions, rows[NAME_COLUMN].astype(str).tolist(), fontsize=7)
    ax.set_xlabel("basin reach-median " + ksn_label(_mn_from(summary, unit_basin)))
    ax.grid(axis="x", alpha=0.3)
    _annotate_units(ax, rows, positions, True,
                    [f"{k} basin{'' if k == 1 else 's'} drawn" for k in n_drawn])
    sized = "qualifying reaches" if ncol == UNIT_BASIN_COLUMNS["n_reaches_qualifying"] else "nodes"
    note = (f"dot = one basin's reach k$_{{sn}}$ median (area ~ {sized}; hollow = low "
            "sample); bar = unit reach median (lsdtt3 summary)")
    if buffer_note(summary):
        note += f"; {buffer_note(summary)}"
    ax.set_title(_caption(title, note, 80), fontsize=8)
    fig.tight_layout()
    return fig, ax
