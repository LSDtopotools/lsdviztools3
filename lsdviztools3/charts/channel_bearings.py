"""Channel-bearing roses and junction-angle charts by group (matplotlib).

Draws the grouped outputs of lsdtt-channel-network-analysis read by
:mod:`lsdviztools3.io.channel_bearings` (groups: Strahler order, basin, region,
lithology; optional stratum):

* :func:`plot_binned_rose` — one rose from lsdtt3's **pre-binned** counts or
  lengths (no re-binning), groups overlaid or stacked, with each group's
  circular mean and R-bar from the summary drawn as a radial line.
* :func:`plot_binned_rose_grid` — one small rose per group value.
* :func:`plot_rose_comparison` — rows = runs (e.g. WV vs an isotropic null),
  columns = groups, on a shared radial scale.
* :func:`plot_junction_angle_summary` — per-group junction-angle boxes drawn
  from the summary's quartiles (``bxp``), optionally several runs side by side,
  with the median CI where lsdtt3 wrote one.
* :func:`plot_junction_angle_distributions` — step histograms of the
  per-junction angles lsdtt3 wrote, one per order group.
* :func:`plot_comparison_table` — lsdtt3's pairwise comparison rows as a table.

Statistics shown as numbers (n, n_blocks, means, R-bar, quartiles, CIs,
p-values) are read from lsdtt3's files, never recomputed. p-values keep their
``_iid_approx`` label and are always shown beside ``n_blocks``: the junctions
and links are spatially clustered, so the i.i.d. p-value is an approximation.

North-up, clockwise polar axes, as :func:`charts.rose.rose_diagram`. Axial
roses (0-180) are mirrored to 360.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from ..io.channel_bearings import (
    AXIAL, BEARING_COMPARISON_COLUMNS, BEARING_SUMMARY_COLUMNS, DIRECTIONAL,
    GROUP_COLUMNS, JUNCTION_COLUMNS, JUNCTION_COMPARISON_COLUMNS,
    JUNCTION_SUMMARY_COLUMNS, NO_STRATUM, ROSE_COLUMNS, order_group_mask,
    select_groups, BLOCK_COLUMNS,
)

_MPL_MSG = (
    "Bearing charts need matplotlib. Install the optional 'charts' extra, e.g. "
    "`pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Weight columns of the rose file, by short name.
WEIGHTS = ("count", "length", "fraction_count", "fraction_length")


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def _sort_values(values):
    """Group values in numeric order ("4+" after 3), else as text.

    Reserved lithology categories (``mixed ...``, ``arms disagree ...``,
    ``unmapped``) and the pooled ``all`` come last, in that order.
    """
    from ..io.lithology import reserved_category
    from ..render.lithology_colours import RESERVED_ORDER

    vals = list(dict.fromkeys(str(v) for v in values))

    def key(v):
        cat = reserved_category(v)
        if cat is not None:
            return (2, RESERVED_ORDER.index(cat), 0.0, v)
        if v == "all":
            return (3, 0, 0.0, v)
        try:
            return (0, 0, float(v.rstrip("+")), v)
        except ValueError:
            return (1, 0, 0.0, v.casefold())

    return sorted(vals, key=key)


def is_lithology_key(key: Optional[str]) -> bool:
    """True for lithology groupings (``lithology``, ``lithology_interior``, ...)."""
    return bool(key) and str(key).startswith("lithology")


def group_colours(group_key: str, values: Sequence, *, cmap: str = "viridis",
                  unit_colours: Optional[Mapping[str, str]] = None) -> dict:
    """Colour per group value: the shared unit palette for lithology groupings
    (same colours as the geologic map and k_sn charts; ``mixed`` grey), else
    samples of ``cmap``."""
    values = [str(v) for v in values]
    if is_lithology_key(group_key):
        from ..render.lithology_colours import lithology_palette

        pal = lithology_palette(values, unit_colours)
        return {v: pal.get(v, "#888889") for v in values}
    plt = _import_mpl()
    cm = plt.get_cmap(cmap)
    return {v: cm(i / max(len(values) - 1, 1)) for i, v in enumerate(values)}


def block_note(*tables) -> str:
    """``"blocks: order 3 (auto), min 20 blocks"`` from lsdtt3's appended
    ``block_order, block_order_mode, min_blocks`` columns (empty if absent)."""
    c = BLOCK_COLUMNS
    for t in tables:
        if t is None or c["order"] not in t.columns:
            continue
        rows = t[pd.to_numeric(t[c["order"]], errors="coerce").notna()]
        if rows.empty:
            continue
        r = rows.iloc[0]
        txt = f"blocks: Strahler order {int(float(r[c['order']]))}"
        mode = r.get(c["mode"]) if c["mode"] in t.columns else None
        if mode is not None and not pd.isna(mode) and str(mode).strip():
            txt += f" ({mode})"
        mb = _num(r, c["min_blocks"]) if c["min_blocks"] in t.columns else None
        if mb is not None:
            txt += f", min {int(mb)} blocks"
        return txt
    return ""


def _with_block_note(note: str, *tables) -> str:
    b = block_note(*tables)
    return f"{note}; {b}" if b else note


def _num(row, col) -> Optional[float]:
    if col not in row.index:
        return None
    v = pd.to_numeric(pd.Series([row[col]]), errors="coerce").iloc[0]
    return float(v) if np.isfinite(v) else None


def _fmt_p(p: Optional[float]) -> str:
    if p is None:
        return "nan"
    return f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}"


def _petals(rows: pd.DataFrame, weight: str):
    """(theta_rad, width_rad, height, axial) for one group's bins; axial mirrored."""
    c = ROSE_COLUMNS
    col = {"count": c["count"], "length": c["length"],
           "fraction_count": c["fraction_count"],
           "fraction_length": c["fraction_length"]}[weight]
    types = rows[c["angular_type"]].astype(str).unique()
    if len(types) != 1:
        raise ValueError(f"Rows mix angular types {list(types)}; select one.")
    axial = types[0] == AXIAL
    period = 180.0 if axial else 360.0
    centre = rows[c["bin_centre"]].to_numpy(float)
    lo, hi = rows[c["bin_lower"]].to_numpy(float), rows[c["bin_upper"]].to_numpy(float)
    # Bin 0 is centred on 0 (lower edge -w/2); widths are taken modulo period.
    width = np.mod(hi - lo, period)
    width = np.where(width == 0, period, width)
    h = rows[col].to_numpy(float)
    if axial:
        centre = np.concatenate([centre, centre + 180.0])
        width = np.concatenate([width, width])
        h = np.concatenate([h, h])
    return np.deg2rad(centre), np.deg2rad(width), h, axial


def _setup_polar(ax):
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)


def _mean_line(ax, srow, *, weight: str, rmax: float, colour, axial: bool):
    """Radial line at the summary circular mean, length R-bar x rmax."""
    s = BEARING_SUMMARY_COLUMNS
    use_len = weight in ("length", "fraction_length")
    mean = _num(srow, s["mean_length"] if use_len else s["mean_count"])
    rbar = _num(srow, s["rbar_length"] if use_len else s["rbar_count"])
    if mean is None or rbar is None:
        return
    for a in ([mean, mean + 180.0] if axial else [mean]):
        ax.plot([np.deg2rad(a)] * 2, [0, rbar * rmax], color=colour, linewidth=2.2,
                solid_capstyle="round", zorder=5)


def bearing_stats_text(srow, weight: str = "fraction_length") -> str:
    """``"n=1765 (106 blocks), R=0.05, Rayleigh p_iid_approx=0.01"`` from a summary row."""
    s = BEARING_SUMMARY_COLUMNS
    use_len = weight in ("length", "fraction_length")
    n, nb = _num(srow, s["n"]), _num(srow, s["n_blocks"])
    rbar = _num(srow, s["rbar_length"] if use_len else s["rbar_count"])
    p = _num(srow, s["rayleigh_p"])
    parts = []
    if n is not None:
        parts.append(f"n={int(n)}" + (f" ({int(nb)} blocks)" if nb is not None else ""))
    if rbar is not None:
        parts.append(f"R={rbar:.3f}")
    if s["rayleigh_p"] in srow.index:
        parts.append(f"Rayleigh p_iid_approx={_fmt_p(p)}")
    if s["meets_min"] in srow.index and _num(srow, s["meets_min"]) == 0:
        parts.append("below min group size")
    if _descriptive_only(srow):
        parts.append("descriptive only")
    return ", ".join(parts)


def _descriptive_only(row) -> bool:
    """True if lsdtt3's summary note marks the group as descriptive only
    (e.g. ``mixed`` / ``arms disagree``: never in comparison rows)."""
    note = row.get("note") if hasattr(row, "get") else None
    return note is not None and not pd.isna(note) and "descriptive only" in str(note)


def _group_rows(df, value):
    return df[df[GROUP_COLUMNS["group_value"]].astype(str) == str(value)]


def plot_binned_rose(
    rose: pd.DataFrame,
    summary: Optional[pd.DataFrame] = None,
    *,
    group_key: str = "strahler_order",
    angular_type: str = DIRECTIONAL,
    groups: Optional[Sequence] = None,
    measurand: Optional[str] = None,
    stratum_key: str = NO_STRATUM,
    stratum_value=None,
    weight: str = "fraction_length",
    stacked: bool = False,
    cmap: str = "viridis",
    unit_colours: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    ax=None,
):
    """Rose of lsdtt3's pre-binned bearings, one colour per group value.

    Args:
        rose: ``<wf>_channel_bearing_rose.csv`` table.
        summary: ``<wf>_channel_bearing_group_summary.csv`` table; if given,
            each group's circular mean is drawn as a line of length R-bar
            (fraction of the outer radius) and n / n_blocks go in the legend.
        group_key: Grouping to draw (``strahler_order``, ``lithology``, ...).
        angular_type: ``"directional"`` (0-360) or ``"axial"`` (0-180,
            mirrored). Both are in the same file, so one must be chosen.
        groups: Subset of group values (default all, in order).
        measurand: Optional measurand filter (``link_bearing``).
        stratum_key, stratum_value: Optional stratum (default unstratified).
        weight: ``count``, ``length``, ``fraction_count`` or ``fraction_length``.
        stacked: Stack petals instead of overlaying outlines.
        cmap: Colour map sampled per group (non-lithology groupings).
        unit_colours: ``{unit name: colour}`` overrides for lithology
            groupings, which otherwise use the shared unit palette.
        title: Axis title.
        ax: Existing polar axes or ``None``.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    if weight not in WEIGHTS:
        raise ValueError(f"weight must be one of {WEIGHTS}, got {weight!r}.")
    sel = select_groups(rose, group_key=group_key, measurand=measurand,
                        stratum_key=stratum_key, stratum_value=stratum_value,
                        angular_type=angular_type)
    if sel.empty:
        raise ValueError(f"No {angular_type} rose rows for group_key={group_key!r}.")
    values = _sort_values(sel[GROUP_COLUMNS["group_value"]]) if groups is None \
        else [str(g) for g in groups]
    ssel = None
    if summary is not None:
        ssel = select_groups(summary, group_key=group_key, measurand=measurand,
                             stratum_key=stratum_key, stratum_value=stratum_value,
                             angular_type=angular_type)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6.2, 5.8), subplot_kw={"projection": "polar"})
    else:
        fig = ax.figure
    _setup_polar(ax)
    colours = group_colours(group_key, values, cmap=cmap, unit_colours=unit_colours)
    bottom, rmax, drawn = None, 0.0, []
    for i, v in enumerate(values):
        rows = _group_rows(sel, v).sort_values(ROSE_COLUMNS["bin_centre"])
        if rows.empty:
            continue
        theta, width, h, axial = _petals(rows, weight)
        colour = colours[str(v)]
        label = f"{group_key.replace('_', ' ')} {v}"
        if ssel is not None and not _group_rows(ssel, v).empty:
            label += f": {bearing_stats_text(_group_rows(ssel, v).iloc[0], weight)}"
        if stacked:
            base = np.zeros_like(h) if bottom is None else bottom
            ax.bar(theta, h, width=width, bottom=base, color=colour, edgecolor="white",
                   linewidth=0.3, align="center", label=label)
            bottom = base + h
            rmax = max(rmax, float(np.max(bottom)))
        else:
            ax.bar(theta, h, width=width, bottom=0.0, fill=False, edgecolor=colour,
                   linewidth=1.5, align="center", label=label)
            rmax = max(rmax, float(np.max(h)))
        drawn.append((v, colour, axial))
    if not drawn:
        raise ValueError("None of the requested groups has rose rows.")
    if ssel is not None:
        for v, colour, axial in drawn:
            srows = _group_rows(ssel, v)
            if not srows.empty:
                _mean_line(ax, srows.iloc[0], weight=weight, rmax=rmax, colour=colour,
                           axial=axial)
    ax.tick_params(labelsize=7)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), fontsize=6.5)
    note = f"{angular_type}, {weight.replace('_', ' ')}; azimuth from grid north"
    if ssel is not None:
        note += "; line = circular mean, length R (lsdtt3 summary)"
    note = _with_block_note(note, summary)
    ax.set_title((title + "\n" if title else "") + note, fontsize=8)
    fig.tight_layout()
    return fig, ax


def plot_binned_rose_grid(
    rose: pd.DataFrame,
    summary: Optional[pd.DataFrame] = None,
    *,
    group_key: str = "strahler_order",
    angular_type: str = DIRECTIONAL,
    weight: str = "fraction_length",
    ncols: int = 4,
    color: str = "steelblue",
    title: Optional[str] = None,
    **select_kwargs,
):
    """One small rose per group value, each with its summary mean line."""
    fig, axes = plot_rose_comparison(
        {"": (rose, summary)}, group_key=group_key, angular_type=angular_type,
        weight=weight, color=color, title=title, ncols=ncols, **select_kwargs)
    return fig, axes


def plot_rose_comparison(
    runs: Mapping[str, tuple],
    *,
    group_key: str = "strahler_order",
    angular_type: str = DIRECTIONAL,
    weight: str = "fraction_length",
    color: Optional[str] = None,
    ncols: Optional[int] = None,
    shared_radius: bool = True,
    unit_colours: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    **select_kwargs,
):
    """Small roses: one row per run, one column per group value.

    Args:
        runs: ``{label: (rose_df, summary_df_or_None)}``, e.g. the WV run and
            an isotropic null surface.
        group_key, angular_type, weight: As :func:`plot_binned_rose`.
        color: Petal colour. ``None`` (default): the shared unit palette for
            lithology groupings (``mixed`` grey), else steel blue.
        ncols: With a single run, wrap the groups into this many columns.
        shared_radius: Same radial limit on every rose (default), so petal
            lengths compare across runs; only meaningful for fractions.
        title: Figure title.
        select_kwargs: ``measurand`` / ``stratum_key`` / ``stratum_value``.

    Returns:
        ``(fig, axes)`` with ``axes`` a flat list of the drawn polar axes.
    """
    plt = _import_mpl()
    sels = {}
    for label, (rose, summary) in runs.items():
        r = select_groups(rose, group_key=group_key, angular_type=angular_type, **select_kwargs)
        s = (select_groups(summary, group_key=group_key, angular_type=angular_type,
                           **select_kwargs) if summary is not None else None)
        sels[label] = (r, s)
    values = _sort_values(pd.concat([r[GROUP_COLUMNS["group_value"]] for r, _ in sels.values()]))
    if not values:
        raise ValueError(f"No {angular_type} rose rows for group_key={group_key!r}.")
    single = len(runs) == 1
    if single and ncols:
        nrows, ncol = int(np.ceil(len(values) / ncols)), ncols
    else:
        nrows, ncol = len(runs), len(values)
    fig, grid = plt.subplots(nrows, ncol, figsize=(2.9 * ncol, 3.7 * nrows),
                             subplot_kw={"projection": "polar"}, squeeze=False)
    flat = list(grid.ravel())
    rmax = 0.0
    drawn = []
    if color is None and is_lithology_key(group_key):
        petal = group_colours(group_key, values, unit_colours=unit_colours)
    else:
        petal = {v: color or "steelblue" for v in values}
    for ri, (label, (r, s)) in enumerate(sels.items()):
        for vi, v in enumerate(values):
            ax = flat[vi] if single else grid[ri][vi]
            rows = _group_rows(r, v).sort_values(ROSE_COLUMNS["bin_centre"])
            _setup_polar(ax)
            ax.tick_params(labelsize=6)
            if rows.empty:
                ax.set_title(f"{label} {v}: no data".strip(), fontsize=7)
                continue
            theta, width, h, axial = _petals(rows, weight)
            ax.bar(theta, h, width=width, color=petal[str(v)], edgecolor="white", linewidth=0.3)
            rmax = max(rmax, float(np.max(h)))
            srow = None
            if s is not None and not _group_rows(s, v).empty:
                srow = _group_rows(s, v).iloc[0]
            head = f"{label}: " if label else ""
            head += f"{group_key.replace('_', ' ')} {v}"
            stats = bearing_stats_text(srow, weight) if srow is not None else ""
            ax.set_title(head + ("\n" + stats.replace(", Rayleigh", ",\nRayleigh") if stats else ""),
                         fontsize=6.5)
            drawn.append((ax, srow, axial, float(np.max(h))))
    for ax in flat[len(values):] if single else []:
        ax.set_visible(False)
    for ax, srow, axial, local in drawn:
        lim = rmax if shared_radius else local
        if srow is not None:
            _mean_line(ax, srow, weight=weight, rmax=lim, colour="black", axial=axial)
        # Set the radial limit last: plotting the mean line re-autoscales.
        ax.set_rlim(0.0, lim * 1.05)
        ax.set_autoscale_on(False)
    note = (f"{angular_type}, {weight.replace('_', ' ')}; azimuth from grid north; "
            "line = circular mean, length R (lsdtt3 summary)")
    note = _with_block_note(note, *[s for _r, s in sels.values()])
    fig.suptitle((title + "\n" if title else "") + note, fontsize=8)
    fig.tight_layout(h_pad=2.5)
    return fig, [d[0] for d in drawn]


def plot_junction_angle_summary(
    summary,
    *,
    group_key: str = "strahler_order",
    measurand: str = "alpha",
    stratum_key: str = NO_STRATUM,
    stratum_value=None,
    unit_colours: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    ax=None,
):
    """Junction-angle boxes per group from the summary's quartiles.

    Box = summary q25-q75, line = median, diamond = mean; no whiskers (the
    summary carries none). Where the summary has ``median_ci_low/high_deg`` the
    median CI is drawn as an error bar. ``n`` and ``n_blocks`` label each box;
    groups below lsdtt3's minimum group size are hatched.

    Args:
        summary: One ``<wf>_junction_angle_group_summary.csv`` table, or a
            mapping ``{run label: table}`` to draw runs side by side.
    """
    plt = _import_mpl()
    c = JUNCTION_SUMMARY_COLUMNS
    runs = summary if isinstance(summary, Mapping) else {"": summary}
    sels = {lab: select_groups(df, group_key=group_key, measurand=measurand,
                               stratum_key=stratum_key, stratum_value=stratum_value)
            for lab, df in runs.items()}
    values = _sort_values(pd.concat([s[c["group_value"]] for s in sels.values()]))
    if not values:
        raise ValueError(f"No junction summary rows for {group_key!r}/{measurand!r}.")
    nrun = len(runs)
    width = 0.7 / nrun
    if ax is None:
        fig, ax = plt.subplots(figsize=(min(1.1 * len(values) * max(nrun, 1) + 3.0, 14.0), 4.8))
    else:
        fig = ax.figure
    colours = plt.get_cmap("tab10")
    handles = []
    for ri, (lab, sel) in enumerate(sels.items()):
        stats, pos, rows = [], [], []
        for vi, v in enumerate(values):
            g = _group_rows(sel, v)
            if g.empty:
                continue
            r = g.iloc[0]
            q1, med, q3 = _num(r, c["q25"]), _num(r, c["median"]), _num(r, c["q75"])
            if None in (q1, med, q3):
                continue
            stats.append({"med": med, "q1": q1, "q3": q3, "whislo": q1, "whishi": q3,
                          "fliers": [], "label": v})
            pos.append(vi + 1 + (ri - (nrun - 1) / 2) * width)
            rows.append(r)
        if not stats:
            continue
        col = "lightsteelblue" if nrun == 1 else colours(ri)
        art = ax.bxp(stats, positions=pos, widths=width * 0.85, patch_artist=True,
                     boxprops={"facecolor": col, "alpha": 0.8},
                     medianprops={"color": "black", "linewidth": 1.4})
        if nrun == 1 and is_lithology_key(group_key):
            ucol = group_colours(group_key, values, unit_colours=unit_colours)
            for box, st in zip(art["boxes"], stats):
                box.set_facecolor(ucol[str(st["label"])])
        for box, r in zip(art["boxes"], rows):
            if _num(r, c["meets_min"]) == 0:
                box.set_hatch("///")
        handles.append((art["boxes"][0], lab))
        means = [_num(r, c["mean"]) for r in rows]
        ax.plot(pos, [m if m is not None else np.nan for m in means], "D",
                color="black", markersize=3)
        lo = [_num(r, c["median_ci_low"]) for r in rows]
        hi = [_num(r, c["median_ci_high"]) for r in rows]
        ok = [i for i in range(len(rows)) if lo[i] is not None and hi[i] is not None]
        if ok:
            meds = np.array([stats[i]["med"] for i in ok])
            ax.errorbar([pos[i] + width * 0.3 for i in ok], meds,
                        yerr=np.vstack([meds - np.array([lo[i] for i in ok]),
                                        np.array([hi[i] for i in ok]) - meds]),
                        fmt="none", ecolor="black", capsize=2, linewidth=1.0)
        for p, r in zip(pos, rows):
            n, nb = _num(r, c["n"]), _num(r, c["n_blocks"])
            txt = f"n={int(n)}" if n is not None else ""
            if nb is not None:
                txt += f"\n{int(nb)} blocks"
            if _descriptive_only(r):
                txt += "\ndescriptive"
            ax.annotate(txt, xy=(p, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(0, 3), textcoords="offset points", ha="center",
                        va="bottom", fontsize=6.5)
    ax.set_xticks(np.arange(1, len(values) + 1), values)
    ax.set_xlabel(group_key.replace("_", " "))
    ax.set_ylabel(f"junction angle {measurand} (degrees)")
    ax.grid(axis="y", alpha=0.3)
    if nrun > 1:
        ax.legend([h for h, _ in handles], [lab for _, lab in handles], fontsize=7,
                  loc="lower right")
    note = "box = q25-q75, line = median, diamond = mean (lsdtt3 summary); hatched = below min group size"
    note = _with_block_note(note, *sels.values())
    if any(c["median_ci_low"] in s.columns for s in sels.values()):
        note += "; bar = median CI"
    ax.set_title((title + "\n" if title else "") + note, fontsize=8, pad=26)
    fig.tight_layout()
    return fig, ax


def plot_junction_angle_distributions(
    junctions: pd.DataFrame,
    summary: pd.DataFrame,
    *,
    group_key: str = "strahler_order",
    measurand: str = "alpha",
    bins=None,
    density: bool = True,
    title: Optional[str] = None,
    ax=None,
):
    """Step histograms of per-junction angles, one per order group.

    Groups (including pooled ``"4+"``) are taken from the summary; junctions
    are assigned by ``junction_order``. Legend counts are the summary's n and
    n_blocks. Uses :func:`charts.distributions.plot_distributions`.
    """
    from .distributions import plot_distributions

    if group_key != "strahler_order":
        raise ValueError("Per-junction distributions are grouped by Strahler order only.")
    acol = JUNCTION_COLUMNS.get(measurand)
    if acol is None or acol not in junctions.columns:
        raise KeyError(f"No per-junction column for measurand {measurand!r}.")
    c = JUNCTION_SUMMARY_COLUMNS
    sel = select_groups(summary, group_key=group_key, measurand=measurand)
    series, styles = {}, {}
    for v in _sort_values(sel[c["group_value"]]):
        r = _group_rows(sel, v).iloc[0]
        mask = order_group_mask(junctions[JUNCTION_COLUMNS["order"]], v)
        n, nb = _num(r, c["n"]), _num(r, c["n_blocks"])
        label = f"order {v} (n={int(n) if n is not None else int(mask.sum())}"
        label += f", {int(nb)} blocks)" if nb is not None else ")"
        series[label] = pd.to_numeric(junctions.loc[mask, acol], errors="coerce").to_numpy()
    if bins is None:
        bins = np.arange(0.0, 181.0, 10.0)
    fig, ax = plot_distributions(
        series, bins=bins, density=density, xlabel=f"junction angle {measurand} (degrees)",
        ylabel="density" if density else "junctions", title=None, ax=ax)
    ax.set_title((title + "\n" if title else "")
                 + f"per-junction {measurand} by junction order (lsdtt3 junction table)",
                 fontsize=8)
    fig.tight_layout()
    return fig, ax


def plot_comparison_table(
    comparisons: pd.DataFrame,
    *,
    kind: str = "bearing",
    measurand: Optional[str] = None,
    angular_type: Optional[str] = None,
    group_key: Optional[str] = None,
    stratum_value=None,
    title: Optional[str] = None,
    ax=None,
):
    """lsdtt3's pairwise group comparisons as a table (values as written).

    Columns: groups, n (n_blocks) per group, the distance (Wasserstein), the
    test statistic (Kuiper V / KS D), ``p_iid_approx`` and its Holm-adjusted
    value. ``n_blocks`` sits beside n because the p-values assume independent
    samples.

    ``group_key`` keeps rows whose ``note`` reads ``group_key=<key>`` (how
    lsdtt3 tells e.g. ``lithology`` from ``lithology_interior`` apart in one
    file); rows with no such note are kept for any key. ``stratum_value``
    selects one stratum (e.g. ``"all"``).
    """
    plt = _import_mpl()
    if kind == "bearing":
        c = BEARING_COMPARISON_COLUMNS
        stat, stat_lab = c["kuiper_V"], "Kuiper V"
        p, p_lab, ph = c["kuiper_p"], "kuiper_p_iid_approx", c["kuiper_p_holm"]
    elif kind == "junction":
        c = JUNCTION_COMPARISON_COLUMNS
        stat, stat_lab = c["ks_D"], "KS D"
        p, p_lab, ph = c["ks_p"], "ks_p_iid_approx", c["ks_p_holm"]
    else:
        raise ValueError("kind must be 'bearing' or 'junction'.")
    df = comparisons
    if measurand is not None:
        df = df[df[c["measurand"]].astype(str) == measurand]
    if angular_type is not None and c.get("angular_type", "") in df.columns:
        df = df[df[c["angular_type"]].astype(str) == angular_type]
    if group_key is not None and c["note"] in df.columns:
        note = df[c["note"]].fillna("").astype(str)
        tagged = note.str.contains("group_key=", regex=False)
        df = df[~tagged | note.str.contains(f"group_key={group_key}", regex=False)
                & ~note.str.contains(f"group_key={group_key}_", regex=False)]
    if stratum_value is not None and c["stratum_value"] in df.columns:
        df = df[df[c["stratum_value"]].astype(str) == str(stratum_value)]
    if df.empty:
        raise ValueError("No comparison rows to tabulate.")

    def _f(v, fmt):
        v = pd.to_numeric(pd.Series([v]), errors="coerce").iloc[0]
        return "nan" if not np.isfinite(v) else format(v, fmt)

    header = ["group a", "group b", "n a (blocks)", "n b (blocks)", "Wasserstein (deg)",
              stat_lab, p_lab, "Holm", "permutation_p", "note (lsdtt3)"]

    def _note(r):
        n = r.get(c["note"]) if c["note"] in df.columns else None
        if n is None or pd.isna(n):
            return ""
        # Drop the grouping tag (already in the title); keep lsdtt3's wording.
        parts = [s.strip() for s in str(n).split(";") if not s.strip().startswith("group_key=")]
        return "; ".join(parts)
    cells = []
    for _, r in df.iterrows():
        cells.append([
            str(r[c["group_a"]]), str(r[c["group_b"]]),
            f"{_f(r[c['n_a']], '.0f')} ({_f(r[c['n_blocks_a']], '.0f')})",
            f"{_f(r[c['n_b']], '.0f')} ({_f(r[c['n_blocks_b']], '.0f')})",
            _f(r[c["wasserstein"]], ".3g"), _f(r[stat], ".3g"),
            _fmt_p(_num(r, p)), _fmt_p(_num(r, ph)),
            _fmt_p(_num(r, c["permutation_p"])) if c["permutation_p"] in df.columns else "nan",
            _note(r),
        ])
    if ax is None:
        fig, ax = plt.subplots(figsize=(16.0, 0.3 * len(cells) + 1.2))
    else:
        fig = ax.figure
    ax.axis("off")
    widths = [0.07, 0.05, 0.075, 0.075, 0.07, 0.055, 0.085, 0.05, 0.065, 0.405]
    tab = ax.table(cellText=cells, colLabels=header, loc="center", cellLoc="center",
                   colWidths=widths)
    tab.auto_set_font_size(False)
    tab.set_fontsize(6.5)
    for (row, col), cell in tab.get_celld().items():
        if col == len(header) - 1 and row > 0:
            cell.set_text_props(ha="left")
            cell._loc = "left"
    tab.scale(1.0, 1.25)
    base = title or f"{kind} pairwise comparisons (lsdtt3; p-values assume i.i.d. samples)"
    ax.set_title(base + "\npermutation_p = block permutation (lsdtt3); nan = not defined",
                 fontsize=8)
    fig.tight_layout()
    return fig, ax


def plot_block_order_counts(counts: pd.DataFrame, *, title: Optional[str] = None, ax=None):
    """lsdtt3's ``<wf>_channel_network_block_order_counts.csv`` as a small table.

    One row per Strahler order: the number of basins with records and whether
    lsdtt3 chose that order to define the permutation blocks.
    """
    from ..io.channel_bearings import BLOCK_COUNT_COLUMNS as c

    plt = _import_mpl()
    rows = []
    for _, r in counts.iterrows():
        chosen = str(r.get(c["chosen"], "")).strip().lower() in ("1", "1.0", "true", "yes")
        rows.append([str(r[c["order"]]), f"{int(float(r[c['n_basins']]))}",
                     "chosen" if chosen else ""])
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 0.3 * len(rows) + 1.0))
    else:
        fig = ax.figure
    ax.axis("off")
    tab = ax.table(cellText=rows, colLabels=["Strahler order", "basins with records", "block order"],
                   loc="center", cellLoc="center")
    tab.auto_set_font_size(False)
    tab.set_fontsize(7)
    ax.set_title(title or "permutation blocks: basins per Strahler order (lsdtt3)", fontsize=8)
    fig.tight_layout()
    return fig, ax


# --------------------------------------------------------------------------
# Bearing distributions: empirical CDF (from lsdtt3's per-link records) and
# PDF (lsdtt3's own binned fractions), with lsdtt3's comparison statistics.
# --------------------------------------------------------------------------

#: Caption fragments (kept in one place so wording stays consistent).
ORIGIN_NOTE = ("0° = grid north, clockwise; axial on [0, 180), directional on "
               "[0, 360)")
CDF_NOTE = ("ECDF of lsdtt3's per-link bearings (a display of the data, not a "
            "statistic). A circular ECDF depends on the chosen origin (here north); "
            "Kuiper V does not.")
PDF_NOTE = ("lsdtt3's binned fractions (rose file), drawn as steps; no smoothing. "
            "Bin 0 is centred on 0° and is drawn split at both ends of the axis.")
NAN_NOTE = "nan = not defined (see lsdtt3 note)"


def _period(angular_type: str) -> float:
    return 180.0 if angular_type == AXIAL else 360.0


def membership_groups(membership: pd.DataFrame, *, angular_type: str = DIRECTIONAL,
                      stratum_value=None) -> tuple:
    """Per-group bearing arrays from lsdtt3's membership rows, plus excluded count.

    Uses ``axial_bearing_deg`` (axial) or ``bearing_deg`` (directional) exactly
    as written. ``stratum_value`` (e.g. ``"3"``) keeps rows of that stratum;
    ``None``/``"all"`` keeps all rows. Rows with an empty ``group_value`` are
    excluded records and are left out.

    Returns:
        ``({group: np.ndarray}, n_excluded)``.
    """
    from ..io.channel_bearings import MEMBERSHIP_COLUMNS as M

    col = M["axial_bearing"] if angular_type == AXIAL else M["bearing"]
    df = membership
    if stratum_value is not None and str(stratum_value) != "all":
        df = df[df[M["stratum_value"]].astype(str) == str(stratum_value)]
    excluded = int(df[M["group_value"]].isna().sum())
    df = df[df[M["group_value"]].notna()]
    out = {}
    for g in _sort_values(df[M["group_value"]]):
        v = pd.to_numeric(df.loc[df[M["group_value"]].astype(str) == g, col],
                          errors="coerce").to_numpy(float)
        out[g] = v[np.isfinite(v)]
    return out, excluded


def check_membership_against_summary(groups: Mapping[str, np.ndarray],
                                     summary: Optional[pd.DataFrame], *, group_key: str,
                                     angular_type: str, stratum_key: str = NO_STRATUM,
                                     stratum_value=None) -> None:
    """Raise if the record counts per group differ from lsdtt3's summary ``n``.

    The ECDF is only shown beside lsdtt3's statistics if it is drawn from the
    same records they were computed on.
    """
    if summary is None:
        return
    s = select_groups(summary, group_key=group_key, stratum_key=stratum_key,
                      stratum_value=stratum_value, angular_type=angular_type)
    for g, v in groups.items():
        r = _group_rows(s, g)
        if r.empty:
            raise ValueError(f"group {g!r} is in the membership file but not in the "
                             f"{group_key!r} summary; the membership is for another grouping.")
        n = _num(r.iloc[0], BEARING_SUMMARY_COLUMNS["n"])
        if n is not None and int(n) != len(v):
            raise ValueError(f"group {g!r}: {len(v)} membership records but lsdtt3 summary "
                             f"n={int(n)}; not drawing an ECDF that differs from the "
                             "statistics' records.")


def _comparison_cells(comparisons: Optional[pd.DataFrame], *, angular_type: str,
                      group_key: Optional[str], stratum_value=None) -> list:
    """Rows [a, b, n a (blocks), n b (blocks), V, W1, p_iid, perm_p, note] from lsdtt3."""
    if comparisons is None:
        return []
    c = BEARING_COMPARISON_COLUMNS
    df = comparisons[comparisons[c["angular_type"]].astype(str) == angular_type]
    if c["stratum_value"] in df.columns:
        sv = "all" if stratum_value is None else str(stratum_value)
        df = df[df[c["stratum_value"]].astype(str) == sv]
    if group_key is not None and c["note"] in df.columns:
        note = df[c["note"]].fillna("").astype(str)
        tagged = note.str.contains("group_key=", regex=False)
        df = df[~tagged | note.str.contains(f"group_key={group_key}", regex=False)
                & ~note.str.contains(f"group_key={group_key}_", regex=False)]

    def _i(v):
        v = pd.to_numeric(pd.Series([v]), errors="coerce").iloc[0]
        return "nan" if not np.isfinite(v) else f"{v:.0f}"

    def _g(v):
        v = pd.to_numeric(pd.Series([v]), errors="coerce").iloc[0]
        return "nan" if not np.isfinite(v) else f"{v:.3g}"

    cells = []
    for _, r in df.iterrows():
        n = r.get(c["note"]) if c["note"] in df.columns else None
        n = "" if n is None or pd.isna(n) else "; ".join(
            s.strip() for s in str(n).split(";") if not s.strip().startswith("group_key="))
        cells.append([str(r[c["group_a"]]), str(r[c["group_b"]]),
                      f"{_i(r[c['n_a']])} ({_i(r[c['n_blocks_a']])})",
                      f"{_i(r[c['n_b']])} ({_i(r[c['n_blocks_b']])})",
                      _g(r[c["kuiper_V"]]), _g(r[c["wasserstein"]]),
                      _fmt_p(_num(r, c["kuiper_p"])),
                      _fmt_p(_num(r, c["permutation_p"])) if c["permutation_p"] in df.columns
                      else "nan", n])
    return cells


_TABLE_HEADER = ["group a", "group b", "n a (blocks)", "n b (blocks)", "Kuiper V",
                 "W1 (deg)", "kuiper_p_iid_approx", "permutation_p", "note"]


def _draw_table(ax, cells):
    """Comparison table; each row's lsdtt3 note is referenced as [k] and the
    distinct notes are listed in full beneath the table."""
    import textwrap

    ax.axis("off")
    if not cells:
        ax.text(0.5, 0.9, "no lsdtt3 comparison rows for this selection",
                ha="center", va="top", fontsize=7, transform=ax.transAxes)
        return
    notes = list(dict.fromkeys(r[-1] for r in cells if r[-1]))
    rows = [r[:-1] + ([f"[{notes.index(r[-1]) + 1}]"] if r[-1] else [""]) for r in cells]
    tab = ax.table(cellText=rows, colLabels=_TABLE_HEADER, loc="upper center",
                   cellLoc="center",
                   colWidths=[0.08, 0.08, 0.13, 0.13, 0.08, 0.08, 0.18, 0.13, 0.06])
    tab.auto_set_font_size(False)
    tab.set_fontsize(6)
    tab.scale(1.0, 1.15)
    if notes:
        text = "\n".join(textwrap.fill(f"[{i + 1}] {n}", 120,
                                        subsequent_indent="     ")
                         for i, n in enumerate(notes))
        ax.text(0.0, -0.02, "lsdtt3 notes:\n" + text, transform=ax.transAxes,
                ha="left", va="top", fontsize=5.5)


def plot_bearing_distributions(
    runs: Mapping[str, Mapping],
    *,
    kind: str = "cdf",
    group_key: str = "strahler_order",
    angular_type: str = AXIAL,
    stratum_key: str = NO_STRATUM,
    stratum_value=None,
    weight: str = "fraction_count",
    cmap: str = "viridis",
    unit_colours: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
):
    """Overlaid per-group bearing CDFs or PDFs, one panel per run, with lsdtt3's
    comparison statistics tabulated under each panel.

    Args:
        runs: ``{label: {"membership": df, "rose": df, "bearing_summary": df,
            "bearing_comparisons": df}}`` (as from
            :func:`io.channel_bearings.load_channel_geometry_outputs`). Several
            runs (e.g. WV and an isotropic null) are drawn side by side with
            the same group colours.
        kind: ``"cdf"`` — empirical CDF of the per-link bearings in lsdtt3's
            membership file (``axial_bearing_deg`` / ``bearing_deg`` as
            written; excluded records left out); ``"pdf"`` — lsdtt3's binned
            ``weight`` fractions from the rose file as a step line.
        group_key: Grouping (``strahler_order``, ``lithology``, ...).
        angular_type: ``"axial"`` or ``"directional"``.
        stratum_key, stratum_value: Stratum (``None`` / ``"all"`` = all records).
        weight: PDF weight: ``fraction_count`` (default) or ``fraction_length``.
        cmap: Colour map for non-lithology groups (lithology uses the shared
            unit palette; mixed grey).
        unit_colours: Unit colour overrides.
        title: Figure title.

    Returns:
        ``(fig, axes)`` — the plot axes, one per run.

    Raises:
        ValueError: if a run's membership records do not match lsdtt3's
            summary ``n`` for the selection (the ECDF would not be of the
            records the statistics describe).
    """
    plt = _import_mpl()
    if kind not in ("cdf", "pdf"):
        raise ValueError("kind must be 'cdf' or 'pdf'.")
    if weight not in ("fraction_count", "fraction_length"):
        raise ValueError("weight must be 'fraction_count' or 'fraction_length'.")
    period = _period(angular_type)
    sv_sel = None if stratum_value in (None, "all") else stratum_value

    # Gather per-run group data first so colours are shared across panels.
    prepared = {}
    all_groups = []
    for label, d in runs.items():
        if kind == "cdf":
            if d.get("membership") is None:
                raise ValueError(f"{label}: no membership file for an ECDF.")
            groups, excl = membership_groups(d["membership"], angular_type=angular_type,
                                             stratum_value=sv_sel)
            sk = NO_STRATUM if sv_sel is None and stratum_key == NO_STRATUM else stratum_key
            check_membership_against_summary(
                groups, d.get("bearing_summary"), group_key=group_key,
                angular_type=angular_type, stratum_key=sk,
                stratum_value="all" if sk != NO_STRATUM and sv_sel is None else sv_sel)
            prepared[label] = (groups, excl)
        else:
            sel = select_groups(d["rose"], group_key=group_key, angular_type=angular_type,
                                stratum_key=stratum_key,
                                stratum_value=("all" if stratum_key != NO_STRATUM
                                               and sv_sel is None else sv_sel))
            groups = {g: _group_rows(sel, g).sort_values(ROSE_COLUMNS["bin_centre"])
                      for g in _sort_values(sel[GROUP_COLUMNS["group_value"]])}
            prepared[label] = (groups, None)
        all_groups += list(prepared[label][0])
    order = _sort_values(all_groups)
    colours = group_colours(group_key, order, cmap=cmap, unit_colours=unit_colours)

    n = len(runs)
    fig = plt.figure(figsize=(7.2 * n, 8.0))
    gs = fig.add_gridspec(2, n, height_ratios=[3.0, 1.5], hspace=0.3,
                          top=0.84, bottom=0.12)
    axes = []
    for i, (label, (groups, excl)) in enumerate(prepared.items()):
        ax = fig.add_subplot(gs[0, i], sharey=axes[0] if axes else None)
        for g in order:
            if g not in groups:
                continue
            if kind == "cdf":
                v = np.sort(groups[g])
                if v.size == 0:
                    continue
                y = np.arange(1, v.size + 1) / v.size
                ax.step(np.concatenate([[0.0], v, [period]]),
                        np.concatenate([[0.0], y, [1.0]]), where="post",
                        color=colours[g], linewidth=1.4, label=f"{g} (n={v.size})")
            else:
                rows = groups[g]
                vals = rows[ROSE_COLUMNS[weight]].to_numpy(float)
                centres = rows[ROSE_COLUMNS["bin_centre"]].to_numpy(float)
                w = float(np.median(np.diff(centres))) if len(centres) > 1 else period
                # Bin 0 straddles 0 deg: draw its half-bins at both ends.
                edges = np.concatenate([[0.0], centres[1:] - w / 2, [period - w / 2, period]])
                heights = np.concatenate([vals, [vals[0]]])
                ax.stairs(heights, edges, color=colours[g], linewidth=1.4, label=g,
                          baseline=None)
        ax.set_xlim(0.0, period)
        ax.set_xticks(np.arange(0.0, period + 1, 45.0 if period == 360 else 30.0))
        ax.set_xlabel(f"{angular_type} bearing (degrees; 0 = grid north, clockwise, "
                      f"[0, {period:g}))")
        if kind == "cdf":
            ax.set_ylabel("cumulative fraction of links")
            ax.set_ylim(0.0, 1.0)
        else:
            what = "links" if weight == "fraction_count" else "channel length"
            ax.set_ylabel(f"fraction of group's {what} per bin")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6.5, loc="upper left" if kind == "cdf" else "upper right",
                  title=group_key.replace("_", " "), title_fontsize=6.5)
        head = label or ""
        if excl:
            head += f" ({excl} excluded record(s) with no group left out)"
        ax.set_title(head, fontsize=8)
        axes.append(ax)
        tax = fig.add_subplot(gs[1, i])
        _draw_table(tax, _comparison_cells(runs[label].get("bearing_comparisons"),
                                           angular_type=angular_type, group_key=group_key,
                                           stratum_value=sv_sel))
    strat = "" if sv_sel is None else f"; {stratum_key} {sv_sel}"
    body = (CDF_NOTE if kind == "cdf" else PDF_NOTE) + f" {ORIGIN_NOTE}{strat}. "
    body += f"Table: lsdtt3 comparisons; p-values beside n (n_blocks); {NAN_NOTE}."
    blocks = block_note(*[d.get("bearing_summary") for d in runs.values()])
    if blocks:
        body += f" {blocks}."
    import textwrap

    fig.suptitle((title + "\n" if title else "") + textwrap.fill(body, 105 * n), fontsize=7.5)
    return fig, axes
