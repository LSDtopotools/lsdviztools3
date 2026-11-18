"""Knickpoint plots for lsdtt3 chi fits.

lsdtt-chi-analysis (``write_chi_knickpoint_table``) writes one row per node
where two fitted k_sn segments meet (``*_chi_knickpoints_mn_*.csv``), with
``delta_ksn`` (upstream minus downstream k_sn: positive where the channel
steepens upstream) and ``log_ratio_ksn``, the split that made the boundary
(``split_delta_ksn_metres``: the change weighted by its chords;
``split_delta_ksn_normalised``: that over the branch's chi span, the ranking
measure) and ``is_major``, and a table of suggested thresholds for a "major"
knickpoint (``*_chi_knickpoint_thresholds_mn_*.csv``: p90, p95, median +
3 MAD, Otsu, over the absolute values). It also writes knickzones (foot/lip
pairs) and, with ``run_step_change_detection``, step reaches and a
step-aware (stage 2) fit. All the numbers come from lsdtt3; this module
selects rows and draws them.

``select_major_knickpoints``  rows beyond a threshold (a value, or a method name)
``plot_knickpoint_distribution``  the signed distribution with the suggested
    thresholds, and the ranked absolute values (for choosing by eye)
``plot_knickpoints_on_profile``  major knickpoints on a single-basin or the
    stacked chi profile
``plot_knickpoints_along_chi``  where major knickpoints sit along chi (or
    normalised chi, or elevation), against all candidates
``render_knickpoint_map``  major knickpoints on the hillshade and k_sn network,
    with knickzone reaches and steps
``plot_knickpoints_by_elevation``  the elevations of the major knickpoints
    against all candidates, overall and per basin (a standard plot)
``plot_steps_on_long_profile``  steps on the longitudinal profile (elevation
    against flow distance), where steps are judged
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

MEASURES = {"split_normalised": ("split_delta_ksn_normalised", "abs_split_delta_ksn_normalised",
                                 r"split $\Delta k_{sn}$, chord-weighted / branch $\chi$ span"),
            "split": ("split_delta_ksn_metres", "abs_split_delta_ksn_metres",
                      r"split $\Delta k_{sn}$, chord-weighted (m)"),
            "delta": ("delta_ksn", "abs_delta_ksn", r"$\Delta k_{sn}$ (upstream $-$ downstream)"),
            "log_ratio": ("log_ratio_ksn", "abs_log_ratio_ksn",
                          r"$\log_{10}(k_{sn,up}/k_{sn,down})$")}
#: Measures of the split hierarchy: heavy tailed, so plotted on log axes.
SPLIT_MEASURES = ("split_normalised", "split")
#: ``threshold="major"``: lsdtt3's own ``is_major`` (its
#: knickpoint_threshold_algorithm_choice on |split_delta_ksn_normalised|).
MAJOR = "major"
KNICKZONE_COLOUR = "#eb6834"
STEP_COLOUR = "#1f5fa8"
STEP_MUTED = "0.6"
#: Steps at a basin outlet or tributary mouth (flagged) on the long profile.
STEP_FLAGGED_COLOUR = "#2e9e4f"
#: Knickpoint selection order 1st..4th: one hue, darkest first (sequential).
ORDER_COLOURS = ("#0b3d7a", "#2f6db5", "#6d9fd6", "#a9c6ea")
METHODS = ("p90", "p95", "median_plus_3mad", "otsu")
# Validated together (light surface); every line is labelled in the legend.
METHOD_COLOURS = {"p90": "#2a78d6", "p95": "#1baf7a", "median_plus_3mad": "#9a5bd6", "otsu": "#eb6834"}
INK = "#3d3d3a"
ALL_COLOUR = "0.82"


def _measure(measure: str):
    if measure not in MEASURES:
        raise ValueError(f"measure must be one of {list(MEASURES)}, got {measure!r}")
    return MEASURES[measure]


def resolve_threshold(threshold, measure: str = "delta", thresholds: Optional[pd.DataFrame] = None
                      ) -> Tuple[float, str]:
    """``(value, label)`` for a threshold given as a number or a method name
    (``p90``, ``p95``, ``median_plus_3mad``, ``otsu``) from the lsdtt3
    thresholds table, on the absolute ``measure``. ``"major"`` is lsdtt3's
    ``is_major`` choice and has no single value here (use
    :func:`select_major_knickpoints`)."""
    _, abs_name, _ = _measure(measure)
    if threshold == MAJOR:
        raise ValueError("threshold 'major' is lsdtt3's is_major column, not a value: "
                         "use select_major_knickpoints")
    try:
        v = float(threshold)
        return v, f"{v:g} (user)"
    except (TypeError, ValueError):
        pass
    if thresholds is None:
        raise ValueError(f"threshold {threshold!r} is a method name: give the thresholds table.")
    hit = thresholds[(thresholds["measure"] == abs_name) & (thresholds["method"] == threshold)]
    if hit.empty or not np.isfinite(hit["threshold"].iloc[0]):
        raise ValueError(f"no {abs_name} threshold named {threshold!r} in the thresholds table")
    v = float(hit["threshold"].iloc[0])
    return v, f"{v:.3g} ({threshold})"


def select_major_knickpoints(
    knickpoints: pd.DataFrame,
    threshold,
    *,
    measure: str = "delta",
    direction: str = "both",
    include_junctions: bool = False,
    thresholds: Optional[pd.DataFrame] = None,
    removed: Optional[pd.DataFrame] = None,
) -> Tuple[pd.DataFrame, float, str]:
    """Select the major knickpoints.

    Args:
        knickpoints: The lsdtt3 knickpoint table.
        threshold: A value (in the measure's units), a method name read from
            ``thresholds``, or ``"major"``: the rows lsdtt3 marked ``is_major``
            (only with ``measure="split_normalised"``, the measure it ranks by;
            the threshold value is then the smallest selected |value|).
        measure: ``"split_normalised"`` (``split_delta_ksn_normalised``, the
            lsdtt3 ranking), ``"split"`` (``split_delta_ksn_metres``),
            ``"delta"`` (``delta_ksn``) or ``"log_ratio"``. The split measures
            need a knickpoint table with the split columns.
        direction: ``"both"`` (|value| above), ``"steepening"`` (value above
            the threshold) or ``"shallowing"`` (value below minus it).
        include_junctions: Keep rows at tributary mouths (``at_junction``).
        removed: Tributaries removed from the plots
            (:func:`lsdviztools3.io.chi.select_tributaries`); their knickpoints
            are dropped too.

    Returns:
        ``(rows, value, label)``: the selected rows (with a ``measure_value``
        column), the threshold value and a label for it.
    """
    col, _, _ = _measure(measure)
    if direction not in ("both", "steepening", "shallowing"):
        raise ValueError(f"direction must be both, steepening or shallowing, got {direction!r}")
    if col not in knickpoints.columns:
        raise ValueError(f"the knickpoint table has no {col!r} column: rerun lsdtt-chi-analysis with a "
                         "current lsdtt3, or use measure='delta'")
    k = candidate_knickpoints(knickpoints, include_junctions=include_junctions, removed=removed)
    v = k[col].to_numpy(dtype=float)
    if threshold == MAJOR:
        if measure != "split_normalised":
            raise ValueError("threshold 'major' (lsdtt3's is_major) ranks by split_normalised; "
                             f"give a value or method for measure {measure!r}")
        if "is_major" not in k.columns:
            raise ValueError("the knickpoint table has no is_major column: rerun with a current lsdtt3")
        keep = k["is_major"].to_numpy() != 0
        if direction == "steepening":
            keep &= v > 0
        elif direction == "shallowing":
            keep &= v < 0
        out = k[keep & np.isfinite(v)].copy()
        out["measure_value"] = out[col]
        a = np.abs(out[col].to_numpy(dtype=float))
        value = float(a.min()) if a.size else float("nan")
        return out, value, f"{value:.3g} (lsdtt3 is_major)"
    value, label = resolve_threshold(threshold, measure, thresholds)
    if direction == "both":
        keep = np.abs(v) > value
    elif direction == "steepening":
        keep = v > value
    else:
        keep = v < -value
    out = k[keep & np.isfinite(v)].copy()
    out["measure_value"] = out[col]
    return out, value, label


def candidate_knickpoints(knickpoints, *, include_junctions=False, removed=None) -> pd.DataFrame:
    """The rows the statistics and selections use: without junction rows
    (unless asked) and without removed tributaries."""
    k = knickpoints
    if not include_junctions:
        k = k[k["at_junction"] == 0]
    if removed is not None and len(removed):
        from ..io.chi import removed_mask
        k = k[~removed_mask(k, removed)]
    return k


def plot_knickpoint_distribution(
    knickpoints: pd.DataFrame,
    thresholds: Optional[pd.DataFrame] = None,
    *,
    measure: str = "delta",
    chosen: Optional[float] = None,
    chosen_label: Optional[str] = None,
    include_junctions: bool = False,
    removed: Optional[pd.DataFrame] = None,
    bins: int = 60,
):
    """The distribution of k_sn changes, for choosing a threshold.

    Left: a histogram of the signed measure (steepening right of 0, shallowing
    left), counts on a log axis so the tails show. Right: the absolute values
    ranked largest first, where a break in slope can be picked by eye. Both
    carry the suggested thresholds from lsdtt3 (one colour per method, its n
    above in the legend) and the chosen threshold (black, ±).

    Returns:
        ``(fig, (ax_hist, ax_ranked))``.
    """
    import matplotlib.pyplot as plt

    col, abs_name, xlabel = _measure(measure)
    k = candidate_knickpoints(knickpoints, include_junctions=include_junctions, removed=removed)
    v = k[col].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    fig, (ah, ar) = plt.subplots(1, 2, figsize=(11, 4), gridspec_kw={"width_ratios": [1.4, 1]})
    if v.size == 0:
        raise ValueError("no knickpoint candidates to plot")
    log_axes = measure in SPLIT_MEASURES
    if log_axes:
        # signed log: 0 at the smallest nonzero |value|, so both tails show
        nz = np.abs(v[v != 0])
        lin = float(np.nanmin(nz)) if nz.size else 1.0
        ah.set_xscale("symlog", linthresh=lin)
        top = float(np.abs(v).max()) or 1.0
        pos = np.geomspace(lin, top * 1.0001, max(bins // 2, 2))
        edges = np.concatenate([-pos[::-1], pos])
        ah.hist(v, bins=edges, color=ALL_COLOUR, edgecolor="white", linewidth=0.5)
    else:
        ah.hist(v, bins=bins, color=ALL_COLOUR, edgecolor="white", linewidth=0.5)
    ah.set_yscale("log")
    ah.set_xlabel(xlabel)
    ah.set_ylabel("number of segment boundaries")
    a = np.sort(np.abs(v))[::-1]
    ar.plot(np.arange(1, a.size + 1), a, color=INK, lw=1.2)
    ar.set_xscale("log")
    if log_axes and (a > 0).any():
        ar.set_yscale("log")
    ar.set_xlabel("rank (largest first)")
    ar.set_ylabel(f"|{xlabel}|")
    handles = []
    if thresholds is not None:
        for method in METHODS:
            hit = thresholds[(thresholds["measure"] == abs_name) & (thresholds["method"] == method)]
            if hit.empty or not np.isfinite(hit["threshold"].iloc[0]):
                continue
            t = float(hit["threshold"].iloc[0])
            n = int((np.abs(v) > t).sum())
            c = METHOD_COLOURS[method]
            for x in (-t, t):
                ah.axvline(x, color=c, lw=1.0, ls="--")
            ar.axhline(t, color=c, lw=1.0, ls="--")
            handles.append(plt.Line2D([], [], color=c, ls="--", label=f"{method}: {t:.3g} ({n} above)"))
    if chosen is not None and np.isfinite(chosen):
        n = int((np.abs(v) >= chosen).sum()) if chosen_label and "is_major" in chosen_label \
            else int((np.abs(v) > chosen).sum())
        for x in (-chosen, chosen):
            ah.axvline(x, color="black", lw=1.8)
        ar.axhline(chosen, color="black", lw=1.8)
        handles.append(plt.Line2D([], [], color="black", lw=1.8,
                                  label=f"chosen: {chosen_label or f'{chosen:g}'} ({n} above)"))
    for ax in (ah, ar):
        ax.grid(alpha=0.25, lw=0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    ah.set_title(f"{v.size} segment boundaries" + ("" if include_junctions else " (tributary mouths left out)"),
                 fontsize=9, loc="left")
    if handles:
        ar.legend(handles=handles, fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    return fig, (ah, ar)


def _marker_sizes(values: np.ndarray, lo: float = 20.0, hi: float = 120.0) -> np.ndarray:
    a = np.abs(np.asarray(values, dtype=float))
    if a.size == 0:
        return a
    amin, amax = float(np.nanmin(a)), float(np.nanmax(a))
    if not amax > amin:
        return np.full(a.size, 0.5 * (lo + hi))
    return lo + (hi - lo) * (a - amin) / (amax - amin)


def _draw_markers(ax, x, y, values, zorder=6):
    v = np.asarray(values, dtype=float)
    s = _marker_sizes(v)
    up = v > 0
    h = []
    if up.any():
        ax.scatter(np.asarray(x)[up], np.asarray(y)[up], s=s[up], marker="^", facecolor="white",
                   edgecolor="black", linewidth=0.9, zorder=zorder)
        h.append(("^", "steepening upstream"))
    if (~up).any():
        ax.scatter(np.asarray(x)[~up], np.asarray(y)[~up], s=s[~up], marker="v", facecolor="0.25",
                   edgecolor="black", linewidth=0.9, zorder=zorder)
        h.append(("v", "shallowing upstream"))
    return h


def _marker_legend(ax, kinds, label_extra="", loc="upper left"):
    """Marker legend; ``loc`` a legend location or a dict of legend keywords."""
    import matplotlib.pyplot as plt
    hs = [plt.Line2D([], [], marker=m, ls="none", markerfacecolor=("white" if m == "^" else "0.25"),
                     markeredgecolor="black", markersize=7, label=lab + label_extra) for m, lab in kinds]
    if hs:
        kw = loc if isinstance(loc, dict) else {"loc": loc}
        ax.legend(handles=hs, fontsize=8, frameon=False, **kw)


def _step_sizes(drops, lo: float = 14.0, hi: float = 90.0) -> np.ndarray:
    """Marker areas for steps, scaled by log drop."""
    d = np.log10(np.maximum(np.asarray(drops, dtype=float), 1e-3))
    if d.size == 0:
        return d
    if not d.max() > d.min():
        return np.full(d.size, 0.5 * (lo + hi))
    return lo + (hi - lo) * (d - d.min()) / (d.max() - d.min())


def _draw_knickzones(ax, kz: pd.DataFrame, shift=None, zorder=3):
    """Knickzones on a chi profile: a broad band from foot to lip, solid for a
    split foot, faint where the steep reach runs to the branch base."""
    kinds = set()
    if kz is None or kz.empty:
        return kinds
    sh = np.zeros(len(kz)) if shift is None else np.asarray(shift, dtype=float)
    for (_, z), d in zip(kz.iterrows(), sh):
        base = z["foot_type"] != "split"
        ax.plot([z["foot_chi"] + d, z["lip_chi"] + d], [z["foot_elevation"], z["lip_elevation"]],
                color=KNICKZONE_COLOUR, lw=7, alpha=0.22 if base else 0.45, solid_capstyle="butt",
                zorder=zorder)
        kinds.add("base" if base else "split")
    return kinds


def _draw_steps(ax, st: pd.DataFrame, x_bottom, x_top, zorder=5, markers=True,
                flagged_colour=STEP_MUTED):
    """Steps on a profile: the reach as a heavy line (with ``markers``, a
    square at its top sized by drop); flagged steps in ``flagged_colour``."""
    kinds = set()
    if st is None or st.empty:
        return kinds
    xb, xt = np.asarray(x_bottom, dtype=float), np.asarray(x_top, dtype=float)
    flagged = st["flagged"].to_numpy(dtype=bool) if "flagged" in st else np.zeros(len(st), dtype=bool)
    sizes = _step_sizes(st["drop_m"])
    for f in (True, False):
        sel = flagged == f
        if not sel.any():
            continue
        c = flagged_colour if f else STEP_COLOUR
        zb, zt = st["bottom_elevation"].to_numpy(dtype=float)[sel], st["top_elevation"].to_numpy(dtype=float)[sel]
        for a, b, za, zb_ in zip(xb[sel], xt[sel], zb, zt):
            ax.plot([a, b], [za, zb_], color=c, lw=2.2 if markers else 2.6, zorder=zorder,
                    solid_capstyle="butt")
        if markers:
            ax.scatter(xt[sel], zt, s=sizes[sel], marker="s", facecolor=c, edgecolor="black",
                       linewidth=0.6, zorder=zorder + 0.5)
        kinds.add("flagged" if f else "step")
    return kinds


def _overlay_legend(ax, marker_kinds, kz_kinds, step_kinds, loc="upper left", step_markers=True,
                    flagged_colour=STEP_MUTED):
    import matplotlib.pyplot as plt
    hs = [plt.Line2D([], [], marker=m, ls="none", markerfacecolor=("white" if m == "^" else "0.25"),
                     markeredgecolor="black", markersize=7, label=lab) for m, lab in sorted(marker_kinds)]
    if "split" in kz_kinds:
        hs.append(plt.Line2D([], [], color=KNICKZONE_COLOUR, lw=7, alpha=0.45, label="knickzone (foot to lip)"))
    if "base" in kz_kinds:
        hs.append(plt.Line2D([], [], color=KNICKZONE_COLOUR, lw=7, alpha=0.22,
                             label="knickzone to the branch base"))
    mk = {"marker": "s", "markeredgecolor": "black"} if step_markers else {}
    if "step" in step_kinds:
        hs.append(plt.Line2D([], [], color=STEP_COLOUR, lw=2.2,
                             label="step (square sized by drop)" if step_markers else "step", **mk))
    if "flagged" in step_kinds:
        hs.append(plt.Line2D([], [], color=flagged_colour, lw=2.2,
                             label="step at outlet / mouth / with uphill pairs" if step_markers
                             else "step at basin outlet or tributary mouth", **mk))
    if hs:
        kw = loc if isinstance(loc, dict) else {"loc": loc}
        # framed: on a busy profile the data runs under the legend
        ax.legend(handles=hs, fontsize=8, frameon=True, framealpha=0.85, edgecolor="none", **kw)


def plot_knickpoints_on_profile(
    segments: pd.DataFrame,
    major: pd.DataFrame,
    *,
    basin: Optional[int] = None,
    basins: Optional[Sequence[int]] = None,
    removed: Optional[pd.DataFrame] = None,
    show_removed: str = "ghost",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    threshold_label: str = "",
    classes=None,
    basins_per_row: Optional[int] = None,
    gradient_key=None,
    knickzones: Optional[pd.DataFrame] = None,
    steps: Optional[pd.DataFrame] = None,
):
    """Major knickpoints on the chi profile: one basin (``basin=``) or the
    stacked plot of ``basins`` (default all). Markers: open upward triangles
    where the channel steepens upstream, filled downward where it shallows,
    sized by the size of the change. ``classes`` colours k_sn in classes, as
    :func:`lsdviztools3.charts.chi.plot_chi_profile`; ``basins_per_row`` and
    ``gradient_key`` as :func:`lsdviztools3.charts.chi.plot_chi_stacked`.
    ``knickzones`` (:func:`lsdviztools3.io.chi.load_chi_knickzones`) are drawn
    as bands from foot to lip; ``steps``
    (:func:`lsdviztools3.io.chi.filter_chi_steps`) as heavy reaches with a
    square at the top sized by drop, flagged ones grey.

    Returns:
        ``(fig, ax)``; ``ax`` is an array of row axes for a multi-row stack.
    """
    from .chi import plot_chi_profile, plot_chi_stacked, stack_rows

    legend_loc = "upper left"
    if basin is not None:
        fig, ax = plot_chi_profile(segments, basin=basin, removed=removed, show_removed=show_removed,
                                   vmin=vmin, vmax=vmax, classes=classes, gradient_key=gradient_key)
        m = major[major["basin_index"] == basin]
        kinds = set(_draw_markers(ax, m["chi"], m["elevation"], m["measure_value"]))
        kz_k, st_k = set(), set()
        if knickzones is not None:
            kz_k = _draw_knickzones(ax, knickzones[knickzones["basin_index"] == basin])
        if steps is not None:
            sb = steps[steps["basin_index"] == basin]
            st_k = _draw_steps(ax, sb, sb["bottom_chi"], sb["top_chi"])
        _overlay_legend(ax, kinds, kz_k, st_k, loc=legend_loc)
        title_ax, title, n_drawn = ax, ax.get_title() + "  ", len(m)
    else:
        sel_basins = sorted(segments["basin_index"].unique()) if basins is None else list(basins)
        fig, ax = plot_chi_stacked(segments, basins=sel_basins, removed=removed,
                                   show_removed=show_removed, vmin=vmin, vmax=vmax,
                                   classes=classes, basins_per_row=basins_per_row,
                                   gradient_key=gradient_key)
        axes = list(np.atleast_1d(ax))
        n_drawn, kinds, kz_k, st_k = 0, set(), set(), set()
        for row_ax, shifts in zip(axes, stack_rows(segments, sel_basins, None, basins_per_row)):
            mr = major[major["basin_index"].isin(list(shifts))]
            x = mr["chi"].to_numpy(dtype=float) + mr["basin_index"].map(shifts).to_numpy(dtype=float)
            kinds |= set(_draw_markers(row_ax, x, mr["elevation"], mr["measure_value"]))
            n_drawn += len(mr)
            if knickzones is not None:
                kr = knickzones[knickzones["basin_index"].isin(list(shifts))]
                kz_k |= _draw_knickzones(row_ax, kr, kr["basin_index"].map(shifts))
            if steps is not None:
                sr = steps[steps["basin_index"].isin(list(shifts))]
                d = sr["basin_index"].map(shifts).to_numpy(dtype=float)
                st_k |= _draw_steps(row_ax, sr, sr["bottom_chi"].to_numpy(dtype=float) + d,
                                    sr["top_chi"].to_numpy(dtype=float) + d)
        _overlay_legend(axes[0], kinds, kz_k, st_k, loc=legend_loc)
        title_ax, title = axes[0], ""
    extra = ""
    if knickzones is not None:
        nk = knickzones[knickzones["basin_index"] == basin] if basin is not None else knickzones
        extra += f", {len(nk)} knickzones"
    if steps is not None:
        ns = steps[steps["basin_index"] == basin] if basin is not None else steps
        extra += f", {len(ns)} steps"
    title_ax.set_title(title + f"{n_drawn} major knickpoints"
                       + (f", threshold {threshold_label}" if threshold_label else "") + extra, fontsize=9)
    return fig, ax


def plot_knickpoints_by_elevation(
    knickpoints: pd.DataFrame,
    major: pd.DataFrame,
    *,
    include_junctions: bool = False,
    removed: Optional[pd.DataFrame] = None,
    bins: int = 30,
    threshold_label: str = "",
):
    """The distribution of knickpoints with elevation, elevation upward.

    Left: histograms over elevation of all candidate segment boundaries
    (grey) and of the major knickpoints, shallowing upstream (dark) and
    steepening upstream (open), counts on a log axis so both show. Right:
    each major knickpoint at its elevation in its basin's column (markers as
    on the profile, sized by the change). A plateau rim or a migrating wave
    shows as a cluster at one elevation across basins.

    Returns:
        ``(fig, (ax_hist, ax_strip))``.
    """
    import matplotlib.pyplot as plt

    k = candidate_knickpoints(knickpoints, include_junctions=include_junctions, removed=removed)
    allz = k["elevation"].to_numpy(dtype=float)
    allz = allz[np.isfinite(allz)]
    mz = major["elevation"].to_numpy(dtype=float)
    mv = major["measure_value"].to_numpy(dtype=float)
    ok = np.isfinite(mz)
    if allz.size == 0 and not ok.any():
        raise ValueError("no knickpoints to plot")
    edges = np.histogram_bin_edges(allz if allz.size else mz[ok], bins=bins)
    fig, (ah, asr) = plt.subplots(1, 2, figsize=(10, 6), sharey=True, gridspec_kw={"width_ratios": [1, 1.4]})
    ah.hist(allz, bins=edges, orientation="horizontal", color=ALL_COLOUR, edgecolor="white", linewidth=0.5,
            label=f"all segment boundaries ({allz.size})")
    down, up = ok & (mv < 0), ok & (mv > 0)
    if down.any():
        ah.hist(mz[down], bins=edges, orientation="horizontal", color="0.25", edgecolor="white",
                linewidth=0.5, label=f"major, shallowing upstream ({down.sum()})")
    if up.any():
        ah.hist(mz[up], bins=edges, orientation="horizontal", histtype="step", color="black", linewidth=1.2,
                label=f"major, steepening upstream ({up.sum()})")
    ah.set_xscale("log")
    ah.set_xlabel("count")
    ah.set_ylabel("elevation (m)")
    ah.legend(fontsize=8, frameon=False, loc="lower right")
    basins_sorted = sorted(int(b) for b in major["basin_index"].unique()) if len(major) else []
    pos = {b: i for i, b in enumerate(basins_sorted)}
    x = major["basin_index"].map(pos).to_numpy(dtype=float)
    kinds = _draw_markers(asr, x, mz, mv)
    asr.set_xticks(range(len(basins_sorted)))
    asr.set_xticklabels([str(b) for b in basins_sorted], rotation=90 if len(basins_sorted) > 20 else 0,
                        fontsize=8 if len(basins_sorted) > 20 else 9)
    asr.set_xlim(-0.7, max(len(basins_sorted), 1) - 0.3)
    asr.set_xlabel("basin")
    _marker_legend(asr, kinds)
    asr.set_title(f"{int(ok.sum())} major knickpoints" + (f", threshold {threshold_label}" if threshold_label else ""),
                  fontsize=9, loc="left")
    for ax in (ah, asr):
        ax.grid(alpha=0.25, lw=0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig, (ah, asr)


def plot_steps_on_long_profile(
    nodes: pd.DataFrame,
    steps: pd.DataFrame,
    *,
    basins: Optional[Sequence[int]] = None,
    main_stem_only: bool = False,
    removed: Optional[pd.DataFrame] = None,
    ncols: int = 3,
):
    """Steps on the longitudinal profile, elevation against flow distance,
    one panel per basin: every branch thin grey (main stem black), each step
    reach a heavy blue segment (its length is the step's flow length), green
    where it is at a basin outlet or tributary mouth (or has uphill pairs).

    Args:
        nodes: Chi fit nodes with ``flow_distance``
            (:func:`lsdviztools3.io.chi.attach_flow_distance`).
        steps: Steps to draw (:func:`lsdviztools3.io.chi.filter_chi_steps`).
        basins: Basins to draw (default: those with steps, most steps first,
            at most 12).
        main_stem_only: Draw only the main stem (and its steps).
        removed: Removed tributaries, left off.

    Returns:
        ``(fig, axes)``.
    """
    import matplotlib.pyplot as plt

    from ..io.chi import removed_mask

    if "flow_distance" not in nodes.columns:
        raise ValueError("nodes need flow_distance: see lsdviztools3.io.chi.attach_flow_distance")
    n = nodes
    if removed is not None and len(removed):
        n = n[~removed_mask(n, removed)]
    if basins is None:
        counts = steps["basin_index"].value_counts()
        basins = [int(b) for b in counts.index[:12]] or sorted(n["basin_index"].unique())[:12]
    basins = list(basins)
    if not basins:
        raise ValueError("no basins to draw")
    nc = max(1, min(int(ncols), len(basins)))
    nr = int(np.ceil(len(basins) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.2 * nr), squeeze=False)
    kinds = set()
    for ax, b in zip(axes.ravel(), basins):
        nb = n[n["basin_index"] == b]
        main = nb[nb["is_main_stem"] != 0]
        if not main_stem_only:
            for _, g in nb[nb["is_main_stem"] == 0].groupby("branch_id"):
                g = g.sort_values("flow_distance")
                ax.plot(g["flow_distance"] / 1000.0, g["observed_elevation"], color="0.75", lw=0.6, zorder=1)
        main = main.sort_values("flow_distance")
        ax.plot(main["flow_distance"] / 1000.0, main["observed_elevation"], color=INK, lw=1.0, zorder=2)
        sb = steps[steps["basin_index"] == b]
        if main_stem_only and len(main):
            sb = sb[sb["branch_id"].isin(main["branch_id"].unique())]
        top_fd = sb["bottom_flow_distance_m"].to_numpy(dtype=float) + sb["flow_length_m"].to_numpy(dtype=float)
        kinds |= _draw_steps(ax, sb, sb["bottom_flow_distance_m"].to_numpy(dtype=float) / 1000.0,
                             top_fd / 1000.0, markers=False, flagged_colour=STEP_FLAGGED_COLOUR)
        ax.set_title(f"basin {b}: {len(sb)} steps", fontsize=9, loc="left")
        ax.set_xlabel("flow distance (km)")
        ax.set_ylabel("elevation (m)")
        ax.grid(alpha=0.25, lw=0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in axes.ravel()[len(basins):]:
        ax.set_visible(False)
    _overlay_legend(axes.ravel()[0], set(), set(), kinds, step_markers=False,
                    flagged_colour=STEP_FLAGGED_COLOUR)
    fig.tight_layout()
    return fig, axes


def knickpoint_selection_order(
    knickpoints: pd.DataFrame,
    *,
    measure: str = "split_normalised",
    include_junctions: bool = False,
    removed: Optional[pd.DataFrame] = None,
    majors_only: bool = False,
    basins: Optional[Sequence[int]] = None,
) -> pd.DataFrame:
    """Order the knickpoints of each channel (basin, branch) by strength,
    separately for each direction: ``direction`` is ``"shallowing"``
    (measure < 0, upstream shallower: a lip) or ``"steepening"`` (> 0: a
    foot), and ``selection_order`` 1 is the channel's strongest boundary of
    that direction (largest |measure|), 2 the next, and so on.
    ``majors_only`` orders only rows lsdtt3 marked ``is_major``; ``basins``
    keeps only those basins. Rows without a finite, nonzero measure are
    dropped. Returns the rows with ``direction``, ``selection_order`` and
    ``measure_value``.
    """
    col, _, _ = _measure(measure)
    k = candidate_knickpoints(knickpoints, include_junctions=include_junctions, removed=removed)
    if basins is not None:
        k = k[k["basin_index"].isin(list(basins))]
    if majors_only:
        k = k[k["is_major"] != 0]
    v = k[col].to_numpy(dtype=float)
    k = k[np.isfinite(v) & (v != 0)].copy()
    k["measure_value"] = k[col]
    k["direction"] = np.where(k[col].to_numpy(dtype=float) < 0, "shallowing", "steepening")
    k["_abs"] = np.abs(k[col].to_numpy(dtype=float))
    k = k.sort_values(["basin_index", "branch_id", "direction", "_abs"],
                      ascending=[True, True, True, False], kind="mergesort")
    k["selection_order"] = k.groupby(["basin_index", "branch_id", "direction"]).cumcount() + 1
    return k.drop(columns="_abs")


def _gaussian_kde(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Gaussian kernel density of ``values`` on ``grid`` (Scott's bandwidth)."""
    v = values[np.isfinite(values)]
    if v.size < 2 or not np.std(v) > 0:
        return np.full(grid.size, np.nan)
    h = np.std(v, ddof=1) * v.size ** (-1.0 / 5.0)
    u = (grid[:, None] - v[None, :]) / h
    return np.exp(-0.5 * u * u).sum(axis=1) / (v.size * h * np.sqrt(2.0 * np.pi))


def plot_knickpoint_order_kde(
    ordered: pd.DataFrame,
    *,
    max_order: int = 4,
    title: str = "",
):
    """Kernel density of the elevation and chi of the 1st, 2nd, ...
    strongest knickpoint of each direction on every channel, pooled over the
    basins in ``ordered`` (:func:`knickpoint_selection_order`; pass it
    ``basins=`` for a subset). Top row: shallowing upstream (lips); bottom
    row: steepening upstream (feet). Left: elevation upward, density across;
    right: density against chi. One line per order, darkest for the 1st;
    the legends give n. A plateau rim shows as a 1st-order shallowing peak.

    Returns:
        ``(fig, axes)``: a 2 x 2 array, rows (shallowing, steepening),
        columns (elevation, chi).
    """
    import matplotlib.pyplot as plt

    if not 1 <= int(max_order) <= len(ORDER_COLOURS):
        raise ValueError(f"max_order must be 1 to {len(ORDER_COLOURS)}, got {max_order}")
    if ordered.empty:
        raise ValueError("no knickpoints to plot")
    counts = ordered.groupby(["direction", "selection_order"]).size()
    if not (counts[counts.index.get_level_values(1) <= int(max_order)] >= 2).any():
        raise ValueError("fewer than 2 knickpoints of every order: no density to draw")
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), gridspec_kw={"width_ratios": [1, 1.4]})
    z_all, c_all = ordered["elevation"].to_numpy(dtype=float), ordered["chi"].to_numpy(dtype=float)
    zg = np.linspace(np.nanmin(z_all), np.nanmax(z_all), 400)
    cg = np.linspace(0.0, np.nanmax(c_all), 400)
    names = {1: "1st", 2: "2nd", 3: "3rd", 4: "4th"}
    for row, (direction, what) in enumerate((("shallowing", "shallowing upstream (lips)"),
                                             ("steepening", "steepening upstream (feet)"))):
        az, ac = axes[row]
        d = ordered[ordered["direction"] == direction]
        drawn = 0
        for o in range(1, int(max_order) + 1):
            r = d[d["selection_order"] == o]
            if len(r) < 2:
                continue
            c = ORDER_COLOURS[o - 1]
            lab = f"{names[o]} on its channel (n = {len(r)})"
            az.plot(_gaussian_kde(r["elevation"].to_numpy(dtype=float), zg), zg, color=c, lw=2, label=lab)
            ac.plot(cg, _gaussian_kde(r["chi"].to_numpy(dtype=float), cg), color=c, lw=2, label=lab)
            drawn += 1
        from matplotlib.ticker import MaxNLocator
        az.xaxis.set_major_locator(MaxNLocator(3))   # densities per metre: long tick labels
        az.set_ylim(zg[0], zg[-1])
        ac.set_xlim(cg[0], cg[-1])
        az.set_ylabel("elevation (m)")
        az.set_xlabel("density")
        ac.set_xlabel(r"$\chi$ (m)")
        ac.set_ylabel("density")
        az.set_title(what, fontsize=10, loc="left")
        if drawn:
            ac.legend(fontsize=8, frameon=False)
        else:
            ac.text(0.5, 0.5, "fewer than 2 per order", transform=ac.transAxes, ha="center", color="0.4")
        for ax in (az, ac):
            ax.grid(alpha=0.25, lw=0.5)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
    n_channels = ordered.groupby(["basin_index", "branch_id"]).ngroups
    n_basins = ordered["basin_index"].nunique()
    fig.suptitle((title + ": " if title else "") + f"knickpoints ordered by strength on each channel "
                 f"({n_channels} channels, {n_basins} basins)", fontsize=10, x=0.01, ha="left")
    fig.tight_layout()
    return fig, axes

def plot_knickpoints_along_chi(
    knickpoints: pd.DataFrame,
    major: pd.DataFrame,
    *,
    along: str = "chi",
    include_junctions: bool = False,
    removed: Optional[pd.DataFrame] = None,
    bins: int = 40,
):
    """Where the major knickpoints sit: a histogram of ``along`` (``chi``,
    ``chi_normalised`` or ``elevation``) of the major knickpoints over that
    of all candidates (grey), and below, one row per basin with each major
    knickpoint as a marker. A cluster at one chi value across tributaries or
    basins is what a migrating wave of incision would leave.

    Returns:
        ``(fig, (ax_hist, ax_strip))``.
    """
    import matplotlib.pyplot as plt

    if along not in ("chi", "chi_normalised", "elevation"):
        raise ValueError(f"along must be chi, chi_normalised or elevation, got {along!r}")
    k = candidate_knickpoints(knickpoints, include_junctions=include_junctions, removed=removed)
    allv = k[along].to_numpy(dtype=float)
    allv = allv[np.isfinite(allv)]
    mv = major[along].to_numpy(dtype=float)
    edges = np.histogram_bin_edges(allv if allv.size else mv, bins=bins)
    fig, (ah, asr) = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                                  gridspec_kw={"height_ratios": [1.3, 1]})
    ah.hist(allv, bins=edges, color=ALL_COLOUR, edgecolor="white", linewidth=0.5,
            label=f"all segment boundaries ({allv.size})")
    ah.hist(mv[np.isfinite(mv)], bins=edges, color=INK, edgecolor="white", linewidth=0.5,
            label=f"major knickpoints ({np.isfinite(mv).sum()})")
    ah.set_ylabel("count")
    ah.legend(fontsize=8, frameon=False)
    basins_sorted = sorted(k["basin_index"].unique())
    pos = {b: i for i, b in enumerate(basins_sorted)}
    y = major["basin_index"].map(pos).to_numpy(dtype=float)
    kinds = _draw_markers(asr, major[along], y, major["measure_value"])
    asr.set_yticks(range(len(basins_sorted)))
    asr.set_yticklabels([str(int(b)) for b in basins_sorted])
    asr.set_ylabel("basin")
    asr.set_ylim(-0.7, len(basins_sorted) - 0.3)
    asr.set_xlabel({"chi": r"$\chi$ (m)", "chi_normalised": r"$\chi$ / basin maximum $\chi$",
                    "elevation": "elevation (m)"}[along])
    _marker_legend(asr, kinds)
    for ax in (ah, asr):
        ax.grid(alpha=0.25, lw=0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig, (ah, asr)


def render_knickpoint_map(
    dem_path,
    segment_lines,
    major: pd.DataFrame,
    png_path,
    *,
    removed: Optional[pd.DataFrame] = None,
    coordinates: str = "local",
    cmap: str = "viridis",
    ksn_vmin: Optional[float] = None,
    ksn_vmax: Optional[float] = None,
    threshold_label: str = "",
    high_res: bool = False,
    legend_position: str = "JTR+jTR+o0.2c",
    ksn_classes=None,
    n_classes: int = 5,
    class_colors=None,
    class_color_mode: str = "auto",
    locator: Optional[dict] = None,
    knickzones: Optional[pd.DataFrame] = None,
    steps: Optional[pd.DataFrame] = None,
    basins: Optional[Sequence[int]] = None,
):
    """Major knickpoints on the hillshade and the k_sn network: open upward
    triangles where the channel steepens upstream, filled downward ones where
    it shallows, sized by the size of the change (log-scaled). Removed
    tributaries are left off the network and their knickpoints are not shown.

    ``ksn_classes`` (with ``n_classes``, ``class_colors``, ``class_color_mode``)
    colours the network in k_sn classes, with a class legend outside the right
    edge instead of the colourbar; see
    :func:`lsdviztools3.workflows.render_chi_ksn_map`. ``locator`` adds a
    locator inset, as there. ``knickzones``: the fitted segments between
    each foot and lip drawn as a broad orange band under the network;
    ``steps``: a square at each step's bottom sized by drop (grey if
    flagged). ``basins``: draw only these basins' network, knickpoints,
    knickzones and steps.

    Returns:
        ``Path`` to the written PNG.
    """
    import io as _io

    import geopandas as gpd

    from ..io.chi import removed_mask
    from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from ..io.vector import load_vector
    from ..render.base import grid_to_display, save_figure
    from ..render.channels import render_channels
    from ..render.decorations import add_scalebar
    from ..render.hillshade import render_hillshade
    from ..render.points import render_points
    from ..render.style import MapStyle

    width = MapStyle.MAX_WIDTH
    da = to_coordinates(load_raster(dem_path), coordinates)
    crs = get_crs(da)
    _, extent = grid_to_display(to_gmt_grid(da), crs)
    gdf = load_vector(str(segment_lines))
    if basins is not None:
        keep = list(basins)
        gdf = gdf[gdf["basin_index"].isin(keep)]
        major = major[major["basin_index"].isin(keep)]
        if knickzones is not None:
            knickzones = knickzones[knickzones["basin_index"].isin(keep)]
        if steps is not None:
            steps = steps[steps["basin_index"].isin(keep)]
    if removed is not None and len(removed):
        gdf = gdf[~removed_mask(gdf, removed)]
    fig = render_hillshade(str(dem_path), coordinates=coordinates, shade=True,
                           style=MapStyle(cmap="gray", figure_size=width, transparent=False))
    kz_rows = []
    if knickzones is not None and len(knickzones):
        reach = knickzone_reach_lines(gdf, knickzones)
        if len(reach):
            fig = render_channels(reach, coordinates=coordinates,
                                  style=MapStyle(line_width="4p", line_color="250/140/80",
                                                 figure_size=width), fig=fig)
            kz_rows.append(f"S 0.3c - 0.5c - 4p,250/140/80 0.9c knickzone reach ({len(knickzones)})")
    fig = render_channels(gdf, color_by="k_sn", log_color=True, cmap=cmap, colorbar=True,
                          colorbar_label="log@-10@-k@-sn@-", coordinates=coordinates,
                          style=MapStyle(line_width="0.9p", cmap=cmap, figure_size=width),
                          fig=fig, vmin=ksn_vmin, vmax=ksn_vmax,
                          classes=ksn_classes, n_classes=n_classes, class_colors=class_colors,
                          class_color_mode=class_color_mode, legend_title="k@-sn@-")
    rows = []
    m = major.dropna(subset=["x", "y"]).copy()
    m["abs_value"] = np.abs(m["measure_value"].to_numpy(dtype=float))
    for sign, symbol, fill, label in ((1, "t", "white", "steepening upstream"),
                                      (-1, "i", "gray25", "shallowing upstream")):
        part = m[np.sign(m["measure_value"]) == sign]
        if part.empty:
            continue
        pts = gpd.GeoDataFrame(part[["abs_value"]].copy(),
                               geometry=gpd.points_from_xy(part["x"], part["y"]), crs=gdf.crs)
        fig = render_points(pts, size_by="abs_value", size_range=(5.0, 12.0), symbol=symbol,
                            pen="0.8p,black", coordinates=coordinates,
                            style=MapStyle(line_color=fill, figure_size=width), fig=fig)
        rows.append(f"S 0.3c {symbol} 0.25c {fill} 0.8p,black 0.9c {label} ({len(part)})")
    if steps is not None and len(steps):
        st = steps.dropna(subset=["bottom_x", "bottom_y"]).copy()
        flagged = st["flagged"].to_numpy(dtype=bool) if "flagged" in st else np.zeros(len(st), dtype=bool)
        for f, fill, label in ((True, "gray60", "step, flagged"), (False, "31/95/168", "step")):
            part = st[flagged == f]
            if part.empty:
                continue
            pts = gpd.GeoDataFrame(part[["drop_m"]].copy(),
                                   geometry=gpd.points_from_xy(part["bottom_x"], part["bottom_y"]), crs=gdf.crs)
            fig = render_points(pts, size_by="drop_m", size_range=(3.0, 8.0), symbol="s",
                                pen="0.6p,black", coordinates=coordinates,
                                style=MapStyle(line_color=fill, figure_size=width), fig=fig)
            rows.append(f"S 0.3c s 0.2c {fill} 0.6p,black 0.9c {label} ({len(part)}, sized by drop)")
    rows.extend(kz_rows)
    if rows:
        if threshold_label:
            rows.insert(0, f"L 9p - L major knickpoints, threshold {threshold_label}")
        fig.legend(spec=_io.StringIO("\n".join(rows)), position=legend_position,
                   box="+gwhite@20+p0.5p,black")
    add_scalebar(fig, extent, crs)
    if locator:
        from ..render.basemap import add_locator_inset
        add_locator_inset(fig, dem_path, **locator)
    save_figure(fig, png_path, style=MapStyle(figure_size=width, transparent=False), high_res=high_res)
    return Path(png_path)


def knickzone_reach_lines(segment_lines, knickzones: pd.DataFrame):
    """The fit segment lines of each knickzone's reach: segments of its basin
    and branch whose chi interval lies within the foot-to-lip chi range."""
    keep = np.zeros(len(segment_lines), dtype=bool)
    b = segment_lines["basin_index"].to_numpy()
    br = segment_lines["branch_id"].to_numpy()
    lo = segment_lines["downstream_chi"].to_numpy(dtype=float)
    hi = segment_lines["upstream_chi"].to_numpy(dtype=float)
    eps = 1e-6
    for _, z in knickzones.iterrows():
        a, c = sorted((float(z["foot_chi"]), float(z["lip_chi"])))
        keep |= ((b == z["basin_index"]) & (br == z["branch_id"]) & (lo >= a - eps) & (hi <= c + eps))
    return segment_lines[keep]
