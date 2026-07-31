"""Divide-migration charts (matplotlib).

Draws what ``lsdtt-divide-migration`` already computed. Two views:

* :func:`plot_pair_comparison` — the side-vs-side pair summary: side A and side
  B means for one metric, with Δ = A − B and its uncertainty beside them.
* :func:`plot_head_metric_vs_distance` — per-head metric against distance to the
  nearest divide, split by unit.

Three things this module is careful about, because getting them wrong makes a
figure worse than no figure:

**The sign.** ``delta_* = mean_*_a − mean_*_b`` and side A is always the
lower-numbered unit. Every pair label spells out which unit is A and which is B,
and the Δ panel states the sign rule that implicates A for that metric.

**NaN.** ``relief_m``, ``gradient_to_divide`` and ``divide_elevation_m`` are
legitimately NaN for some heads, and a whole pair-summary row is NaN when a side
has fewer than two valid values. NaN is dropped explicitly and the count of what
was dropped is reported in the legend/caption — it never becomes 0.

**The uncertainty column.** ``std_dev`` / ``std_err`` / ``bootstrap`` are
dispersions and are drawn as error bars; ``ttest`` is a **p-value** and is drawn
as a significance annotation instead. The method is read from the file.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

from ..io.divide_migration import (
    resolve_metric,
    uncertainty_method as _uncertainty_method,
)

_MPL_MSG = (
    "Divide-migration plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Colours for the two sides of a divide. A is the lower-numbered unit.
COLOR_A = "#1b6ca8"
COLOR_B = "#c8531b"


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(_MPL_MSG) from exc
    return plt


def _pair_label(row) -> str:
    """``"1|2\\nA=1  B=2"`` — the sign convention spelled out on every tick."""
    return f"{row['pair_id']}\nA={int(row['unit_a'])}  B={int(row['unit_b'])}"


def aggressor_annotation(row, metric, method, *, alpha: float = 0.05):
    """Describe which side the Δ points at, and whether that can be believed.

    Presentation only: it reads the sign of the Δ the CLI wrote and compares its
    magnitude with the uncertainty column according to what that column means.
    Nothing is recomputed from the head data.

    Args:
        row: One pair-summary row.
        metric: A :class:`~lsdviztools3.io.divide_migration.DivideMetric`.
        method: The run's
            :class:`~lsdviztools3.io.divide_migration.UncertaintyMethod`.
        alpha: p-value threshold used when the method is ``ttest``.

    Returns:
        ``(text, state)`` where ``state`` is ``"strong"`` (Δ clears the
        uncertainty / p < alpha), ``"weak"`` (it does not), ``"untested"`` (the
        column supports no such test, e.g. ``std_dev``), or ``"nodata"``.
    """
    delta = row.get(metric.delta_column, np.nan)
    if not np.isfinite(delta):
        return "no Δ (NaN)", "nodata"

    side_unit = row["unit_a"] if np.sign(delta) == metric.aggressor_delta_sign else row["unit_b"]
    side_name = "A" if side_unit == row["unit_a"] else "B"
    who = f"unit {int(side_unit)} ({side_name})"

    unc = row.get(metric.uncertainty_column, np.nan)
    if method.key == "ttest":
        if not np.isfinite(unc):
            return f"{who}, p = NaN", "untested"
        state = "strong" if unc < alpha else "weak"
        return f"{who}, p = {unc:.3g}", state

    if method.significance_multiplier is None or not np.isfinite(unc):
        return f"{who} (not tested)", "untested"

    state = "strong" if abs(delta) > method.significance_multiplier * unc else "weak"
    return (who if state == "strong" else f"{who} (n.s.)"), state


_STATE_COLORS = {
    "strong": "black",
    "weak": "0.45",
    "untested": "0.45",
    "nodata": "0.6",
}


def plot_pair_comparison(
    summary: pd.DataFrame,
    *,
    metric: str = "chi",
    pairs: Optional[Sequence[str]] = None,
    alpha: float = 0.05,
    figsize: Optional[tuple] = None,
    title: Optional[str] = None,
):
    """Plot side A vs side B means and Δ = A − B for one metric, per divide pair.

    Left panel: the two side means joined by a bar, A in blue and B in orange,
    with each side's valid-value count ``n`` printed beside its marker (that n,
    not ``n_heads_*``, is what the mean came from). Right panel: Δ = A − B with
    its uncertainty, and an annotation naming the side the Δ implicates as the
    aggressor.

    Pairs whose Δ **and** both means are NaN are dropped and counted in the
    caption. A pair with only one valid side keeps its row: the available mean is
    drawn and the empty side is labelled, since "one side has no heads" is a
    result, not a gap.

    Args:
        summary: Pair-summary table from
            :func:`~lsdviztools3.io.divide_migration.load_pair_summary`.
        metric: ``chi`` / ``relief`` / ``gradient`` / ``head_elevation`` (the
            channel-head column names are accepted as aliases).
        pairs: Optional subset/order of ``pair_id`` values.
        alpha: p-value threshold, used only when the run's uncertainty method is
            ``ttest``.
        figsize: Matplotlib figure size; sized from the pair count if ``None``.
        title: Figure title; a default naming the metric is used if ``None``.

    Returns:
        ``(fig, axes)`` with ``axes = (means_ax, delta_ax)``.

    Raises:
        KeyError: if the required summary columns are absent.
        ValueError: if no pair has any plottable value.
    """
    plt = _import_mpl()
    metric = resolve_metric(metric)
    method = _uncertainty_method(summary)

    required = {"pair_id", "unit_a", "unit_b", metric.delta_column,
                metric.mean_a_column, metric.mean_b_column}
    missing = required - set(summary.columns)
    if missing:
        raise KeyError(
            f"Pair summary is missing columns for metric {metric.key!r}: "
            f"{sorted(missing)}"
        )

    work = summary.copy()
    work["pair_id"] = work["pair_id"].astype(str)
    if pairs is not None:
        wanted = [str(p) for p in pairs]
        work = work.set_index("pair_id").reindex(wanted).reset_index()

    values = work[[metric.mean_a_column, metric.mean_b_column, metric.delta_column]]
    keep = np.isfinite(values.to_numpy(dtype=float)).any(axis=1)
    n_dropped = int((~keep).sum())
    work = work[keep]
    if work.empty:
        raise ValueError(
            f"No pair has a finite {metric.key!r} mean or Δ — nothing to plot."
        )

    n_pairs = len(work)
    if figsize is None:
        figsize = (10.0, max(2.6, 0.95 * n_pairs + 1.9))
    fig, (ax_means, ax_delta) = plt.subplots(
        1, 2, figsize=figsize, sharey=True,
        gridspec_kw={"width_ratios": [1.35, 1.0]},
    )

    y = np.arange(n_pairs)[::-1]  # first pair at the top
    labels = [_pair_label(row) for _, row in work.iterrows()]

    # ---------------- left: the two side means ----------------
    for yi, (_, row) in zip(y, work.iterrows()):
        a = float(row[metric.mean_a_column])
        b = float(row[metric.mean_b_column])
        if np.isfinite(a) and np.isfinite(b):
            ax_means.plot([a, b], [yi, yi], color="0.6", linewidth=1.6, zorder=1)
        for value, color, marker, side in (
            (a, COLOR_A, "o", "a"), (b, COLOR_B, "s", "b")
        ):
            if not np.isfinite(value):
                continue
            ax_means.plot(value, yi, marker=marker, color=color, markersize=8,
                          zorder=3, linestyle="none")
            n_col = metric.n_a_column if side == "a" else metric.n_b_column
            n_valid = row.get(n_col, np.nan)
            if np.isfinite(n_valid):
                # A above the marker, B below: the two means can coincide (a
                # stable divide), and overlapping n labels are unreadable.
                offset = (0, 9) if side == "a" else (0, -16)
                ax_means.annotate(
                    f"n={int(n_valid)}", (value, yi), textcoords="offset points",
                    xytext=offset, ha="center", fontsize=7, color=color,
                )
        for value, side_label in ((a, "A"), (b, "B")):
            if not np.isfinite(value):
                ax_means.annotate(
                    f"side {side_label}: no valid values",
                    (0.02, yi), xycoords=("axes fraction", "data"),
                    fontsize=7.5, color="0.4", va="center",
                )

    ax_means.set_yticks(y)
    ax_means.set_yticklabels(labels, fontsize=8)
    ax_means.set_xlabel(metric.axis_label)
    ax_means.grid(True, axis="x", alpha=0.25)
    ax_means.set_ylim(-0.7, n_pairs - 0.3)
    ax_means.plot([], [], marker="o", color=COLOR_A, linestyle="none",
                  label="side A (lower-numbered unit)")
    ax_means.plot([], [], marker="s", color=COLOR_B, linestyle="none",
                  label="side B")
    ax_means.legend(fontsize=7.5, loc="best", framealpha=0.9)

    # ---------------- right: Δ = A − B ----------------
    deltas = work[metric.delta_column].to_numpy(dtype=float)
    uncs = (work[metric.uncertainty_column].to_numpy(dtype=float)
            if metric.uncertainty_column in work.columns
            else np.full(n_pairs, np.nan))

    finite_delta = np.isfinite(deltas)
    draw_bars = method.is_dispersion
    if draw_bars:
        xerr = np.where(np.isfinite(uncs) & finite_delta, uncs, 0.0)
        ax_delta.errorbar(
            np.where(finite_delta, deltas, np.nan), y, xerr=xerr, fmt="D",
            color="black", markersize=5.5, capsize=3.5, linewidth=1.2,
            linestyle="none", zorder=3,
        )
    else:
        # ttest: the column is a p-value, so no error bar is drawn.
        ax_delta.plot(np.where(finite_delta, deltas, np.nan), y, "D",
                      color="black", markersize=5.5, linestyle="none", zorder=3)

    ax_delta.axvline(0.0, color="0.4", linewidth=1.0, linestyle="--")

    # Aggressor annotations, right-aligned inside the panel.
    for yi, (_, row) in zip(y, work.iterrows()):
        text, state = aggressor_annotation(row, metric, method, alpha=alpha)
        ax_delta.annotate(
            text, (0.985, yi), xycoords=("axes fraction", "data"),
            ha="right", va="center", fontsize=7.5,
            color=_STATE_COLORS[state],
            fontweight="bold" if state == "strong" else "normal",
        )

    if draw_bars:
        bar_note = f"error bar: {method.label}"
    else:
        bar_note = f"no error bar — {method.label}"
    ax_delta.set_xlabel(f"Δ {metric.label} = A − B\n{bar_note}", fontsize=9)
    ax_delta.grid(True, axis="x", alpha=0.25)
    # Headroom on the right so the annotations do not collide with the markers.
    lo, hi = ax_delta.get_xlim()
    span = hi - lo if hi > lo else 1.0
    ax_delta.set_xlim(lo - 0.05 * span, hi + 0.75 * span)

    if title is None:
        title = f"Divide pairs: {metric.label.lower()} across the divide"
    fig.suptitle(title, y=0.985, fontsize=11)

    caption = (
        f"Δ = mean(A) − mean(B); A is the lower-numbered unit. "
        f"{metric.aggressor_rule} (and B when the sign reverses); the annotated "
        f"unit is the aggressor, so the divide is expected to migrate away from "
        f"it. Uncertainty method: {method.key} — {method.label}. "
        f"n is the valid-value count per side (NaN heads excluded)."
    )
    if n_dropped:
        caption += f" {n_dropped} pair(s) dropped: all-NaN for this metric."
    if method.key == "std_dev":
        caption += (
            " std_dev is a spread of the head populations, not an interval on "
            "Δ, so no significance is claimed."
        )
    fig.text(0.01, 0.005, caption, fontsize=7.2, color="0.3", ha="left",
             va="bottom", wrap=True)

    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    return fig, (ax_means, ax_delta)


def plot_head_metric_vs_distance(
    heads,
    *,
    metric: str = "chi",
    pair_id: Optional[str] = None,
    group_column: str = "unit_id",
    figsize: tuple = (7.6, 4.8),
    title: Optional[str] = None,
):
    """Scatter each channel head's metric against its distance to the divide.

    This is the evidence behind a pair-summary row: one point per head, split by
    unit so the two sides of a divide can be compared directly. Heads whose
    metric or distance is NaN are dropped and the count is reported per group in
    the legend.

    Args:
        heads: Channel-head table — the metrics CSV, or the point vector
            (a ``GeoDataFrame`` works; only attributes are used).
        metric: ``chi`` / ``relief`` / ``gradient`` / ``head_elevation``.
        pair_id: Restrict to heads whose ``pair_id_nearest`` is this pair.
        group_column: Column to split/colour the heads by (default ``unit_id``).
        figsize: Matplotlib figure size.
        title: Figure title; defaulted from the metric and pair if ``None``.

    Returns:
        ``(fig, ax)``.

    Raises:
        KeyError: if a required column is absent.
        ValueError: if no head survives the NaN filter.
    """
    plt = _import_mpl()
    metric = resolve_metric(metric)

    df = pd.DataFrame(heads).copy()
    needed = {metric.head_column, "distance_to_divide_m", group_column}
    missing = needed - set(df.columns)
    if missing:
        raise KeyError(f"Channel-head table is missing columns: {sorted(missing)}")

    if pair_id is not None:
        if "pair_id_nearest" not in df.columns:
            raise KeyError("Channel-head table has no 'pair_id_nearest' column.")
        df = df[df["pair_id_nearest"].astype(str) == str(pair_id)]

    df[metric.head_column] = pd.to_numeric(df[metric.head_column], errors="coerce")
    df["distance_to_divide_m"] = pd.to_numeric(
        df["distance_to_divide_m"], errors="coerce"
    )

    fig, ax = plt.subplots(figsize=figsize)
    total_dropped = 0
    plotted = 0
    for i, (group, sub) in enumerate(df.groupby(group_column, sort=True)):
        valid = np.isfinite(sub[metric.head_column].to_numpy(dtype=float)) & np.isfinite(
            sub["distance_to_divide_m"].to_numpy(dtype=float)
        )
        dropped = int((~valid).sum())
        total_dropped += dropped
        sub = sub[valid]
        if sub.empty:
            continue
        plotted += len(sub)
        label = f"{group_column}={group} (n={len(sub)}"
        label += f", {dropped} NaN dropped)" if dropped else ")"
        ax.plot(
            sub["distance_to_divide_m"], sub[metric.head_column],
            marker="o", linestyle="none", markersize=4.5, alpha=0.75,
            color=f"C{i % 10}", label=label,
        )

    if plotted == 0:
        raise ValueError(
            f"No channel head has finite {metric.head_column!r} and "
            "distance_to_divide_m — nothing to plot."
        )

    ax.set_xlabel("Distance to nearest divide (m)")
    per_head_label = metric.axis_label.replace("Mean ", "", 1)
    ax.set_ylabel(per_head_label[:1].upper() + per_head_label[1:])
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    if title is None:
        title = f"Channel heads vs distance to divide"
        if pair_id is not None:
            title += f" — pair {pair_id}"
    ax.set_title(title, fontsize=11)
    if total_dropped:
        fig.text(
            0.01, 0.005,
            f"{total_dropped} head(s) with NaN {metric.head_column} or distance "
            "were dropped (nearest divide point on nodata, or zero distance).",
            fontsize=7.2, color="0.3", ha="left", va="bottom",
        )
        fig.tight_layout(rect=(0, 0.05, 1, 1))
    else:
        fig.tight_layout()
    return fig, ax
