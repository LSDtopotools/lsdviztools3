"""Spatial (per-sample) transient cosmogenic scatter charts (matplotlib).

Draws the two scatter families produced by sampling a step-change landscape
along its channel network and inverting each sample's **basin-averaged**
concentration:

* **Family B** — apparent erosion-rate ratio against time since the knickpoint
  passed, coloured and sized by the basin's true mean erosion rate
  (:func:`plot_ratio_vs_steptime`).
* **Family C** — the same ratio against that true mean erosion rate, coloured
  and sized by time since the step (:func:`plot_ratio_vs_mean_erate`).
* **Family D** — Family B with the equivalent 1-D column sweep drawn over it,
  so the routed landscape and the analytical column can be read against each
  other (:func:`plot_ratio_overlay_1d`).

Every point is one detrital sample. Both the concentration behind the apparent
rate and the erosion rate it is compared against are basin averages, because
that is what real cosmogenic sampling measures; the local pixel values in the
joined table are diagnostics and are deliberately not drawn here.

**Visualisation only.** Every quantity drawn is a column that already exists in
the joined CSV -- the apparent rates, their ratios, the upstream mean erosion
rate and its unit variants are all produced upstream by lsdtt3 or by the join.
This module never divides, normalises or inverts anything.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

_MPL_MSG = (
    "Cosmo transient plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Ratio columns that can be drawn on the y axis, with mathtext labels.
RATIO_LABELS = {
    "app_erate_ratio_C14_Be10":
        "apparent erosion $^{14}C$/\napparent erosion $^{10}Be$",
    "app_erate_ratio_C14_Al26":
        "apparent erosion $^{14}C$/\napparent erosion $^{26}Al$",
}

_STEPTIME_YR_COL = "sc_steptime_yr"
_STEPTIME_KYR_COL = "sc_steptime_kyr"
_MEAN_ERATE_COL = "upstream_mean_erate_mm_yr"

#: Time column in a ``lsdtt-cosmo-column-tool`` step-change sweep table. It is
#: the same physical quantity as :data:`_STEPTIME_YR_COL` -- years since the
#: step arrived -- which is what makes the overlay a legitimate comparison.
_SWEEP_TIME_COL = "time_since_step_yr"

_STEPTIME_LABEL = "Time since passage of knickpoint (yrs)"
_MEAN_ERATE_LABEL = "Average upstream erosion (mm/yr)"
_MEAN_ERATE_CBAR = "Average upstream erosion rate (mm/yr)"
_STEPTIME_CBAR = "Time since step change (kyr)"

_AXIS_LABEL_SIZE = 22
_TICK_LABEL_SIZE = 16
_CBAR_LABEL_SIZE = 16

#: Marker areas (pt^2) the size encoding runs between.
_SIZE_MIN = 12.0
_SIZE_MAX = 190.0

#: Default colour map. A rainbow is a poor default in general, but these figures
#: encode a variable that is ALSO shown by marker size, so hue only has to
#: separate neighbouring points rather than carry the reading on its own.
_DEFAULT_CMAP = "jet"

#: Line styles cycled by the 1-D sweep curves in the overlay, in draw order.
#: Neutral greys and dash patterns, never colours: the scatter beneath already
#: owns the colour channel, so a coloured curve would read as a data value.
_SWEEP_STYLES = (
    ("black", "-"),
    ("0.35", "--"),
    ("0.55", "-."),
    ("0.70", ":"),
)

#: Default symlog linear threshold (yr) for the overlay x axis. Below it the
#: axis is linear, so samples at exactly t = 0 -- the anchor population, whose
#: ratio of 1.0 is the strongest check in the experiment -- are drawn rather
#: than silently dropped as a log axis would drop them.
#:
#: 100 rather than a rounder 1000: matplotlib still lays decade ticks inside
#: the linear region, so a wider threshold crushes the 0 and 10^2 labels into
#: each other at the origin.
_OVERLAY_LINTHRESH = 100.0

#: Fraction of the y range added above the data to seat an upper-corner legend
#: clear of the curve peak. Purely axes headroom -- no drawn value moves.
_LEGEND_HEADROOM_FRAC = 0.16


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def load_joined_samples(path: str | Path) -> pd.DataFrame:
    """Load one ``<stem>_joined_samples.csv``.

    Args:
        path: Path to the joined sample table. Leading ``#`` comment lines are
            tolerated.

    Returns:
        The table as-is. Nothing is computed, renamed or dropped.
    """
    return pd.read_csv(path, comment="#")


def _require_columns(df: pd.DataFrame, required: Sequence[str]) -> None:
    missing = set(required) - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")


def _finite_xyc(df: pd.DataFrame, x_col: str, y_col: str, c_col: str):
    """Rows where all three of x, y and the colour/size variable are finite.

    A sample missing any one of them cannot be drawn faithfully -- dropping it
    is honest, whereas plotting it at a default colour or size would invent a
    value. The count of dropped rows is returned so callers can report it.
    """
    x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[y_col], errors="coerce").to_numpy(dtype=float)
    c = pd.to_numeric(df[c_col], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
    return x[keep], y[keep], c[keep], int((~keep).sum())


def _sizes_from(values: np.ndarray) -> np.ndarray:
    """Linear map from the colour variable onto marker area.

    Size is redundant with colour by design: it is what makes the 2022 figure
    family readable in greyscale and at small print size. A degenerate range
    (all values equal) maps to the midpoint rather than dividing by zero.
    """
    lo = float(np.min(values))
    hi = float(np.max(values))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.full(values.shape, 0.5 * (_SIZE_MIN + _SIZE_MAX))
    frac = (values - lo) / (hi - lo)
    return _SIZE_MIN + frac * (_SIZE_MAX - _SIZE_MIN)


def _new_or_adopted_axes(plt, ax, figsize: Tuple[float, float]):
    if ax is None:
        # Constrained layout rather than tight_layout: these figures carry a
        # colourbar with a long label, and tight_layout does not reserve room
        # for it, so the label is clipped at the right edge. Both keep the
        # figure size fixed (unlike bbox_inches="tight"), so the saved pixel
        # size is unaffected.
        fig, ax = plt.subplots(figsize=figsize, layout="constrained")
    else:
        fig = ax.figure
    return fig, ax


def _scatter(
    df: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    c_col: str,
    xlabel: str,
    ylabel: str,
    cbar_label: str,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    show_colorbar: bool,
    figsize: Tuple[float, float],
    ax,
    edgecolor,
    linewidth: float,
):
    plt = _import_mpl()
    _require_columns(df, (x_col, y_col, c_col))

    x, y, c, n_dropped = _finite_xyc(df, x_col, y_col, c_col)
    if x.size == 0:
        raise ValueError(
            f"No rows with finite '{x_col}', '{y_col}' and '{c_col}' to plot."
        )

    fig, ax = _new_or_adopted_axes(plt, ax, figsize)
    sc = ax.scatter(
        x, y, c=c, s=_sizes_from(c), cmap=cmap, vmin=vmin, vmax=vmax,
        edgecolors=edgecolor, linewidths=linewidth,
    )

    ax.set_xlabel(xlabel, fontsize=_AXIS_LABEL_SIZE)
    ax.set_ylabel(ylabel, fontsize=_AXIS_LABEL_SIZE)
    ax.tick_params(labelsize=_TICK_LABEL_SIZE)

    if show_colorbar:
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label(cbar_label, fontsize=_CBAR_LABEL_SIZE)
        cb.ax.tick_params(labelsize=_TICK_LABEL_SIZE)

    # Layout is handled by constrained_layout at figure creation (see
    # _new_or_adopted_axes); an adopted axes keeps whatever its owner set.
    return fig, ax, n_dropped


def plot_ratio_vs_steptime(
    df: pd.DataFrame,
    *,
    ratio_col: str = "app_erate_ratio_C14_Be10",
    figsize: Tuple[float, float] = (10.0, 5.0),
    cmap: str = _DEFAULT_CMAP,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_colorbar: bool = True,
    edgecolor: str = "none",
    linewidth: float = 0.0,
    ax=None,
):
    """Family B: apparent erosion-rate ratio against time since the knickpoint.

    One marker per detrital sample, coloured **and** sized by the basin's true
    mean erosion rate, so a sample's position, hue and area together show how
    far the inferred rate has drifted from the real one and how fast the basin
    is actually eroding.

    Args:
        df: One joined sample table. Must carry ``sc_steptime_yr``,
            ``upstream_mean_erate_mm_yr`` and ``ratio_col`` -- all produced
            upstream, none computed here.
        ratio_col: Which apparent-rate ratio to draw. See :data:`RATIO_LABELS`.
        figsize: Figure size, used only when creating a new figure.
        cmap: Colour map for the erosion-rate encoding.
        vmin: Lower colour limit, or ``None`` to autoscale.
        vmax: Upper colour limit, or ``None`` to autoscale.
        show_colorbar: Draw the colourbar.
        edgecolor: Marker edge colour.
        linewidth: Marker edge width.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Note:
        No gridlines are drawn -- the same deliberate deviation from the house
        ``ax.grid(True, linestyle=":")`` convention as the sweep charts, since
        the published figure family this reproduces has none.

    Returns:
        ``(fig, ax)``. This function never saves -- see
        :func:`lsdviztools3.workflows.render_cosmo_transient_spatial_figures`.
    """
    if ratio_col not in RATIO_LABELS:
        raise KeyError(
            f"Unknown ratio column '{ratio_col}'. "
            f"Expected one of {sorted(RATIO_LABELS)}."
        )
    fig, ax, _ = _scatter(
        df,
        x_col=_STEPTIME_YR_COL, y_col=ratio_col, c_col=_MEAN_ERATE_COL,
        xlabel=_STEPTIME_LABEL, ylabel=RATIO_LABELS[ratio_col],
        cbar_label=_MEAN_ERATE_CBAR, cmap=cmap, vmin=vmin, vmax=vmax,
        show_colorbar=show_colorbar, figsize=figsize, ax=ax,
        edgecolor=edgecolor, linewidth=linewidth,
    )
    return fig, ax


def plot_ratio_vs_mean_erate(
    df: pd.DataFrame,
    *,
    ratio_col: str = "app_erate_ratio_C14_Be10",
    figsize: Tuple[float, float] = (10.0, 5.0),
    cmap: str = _DEFAULT_CMAP,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_colorbar: bool = True,
    edgecolor: str = "none",
    linewidth: float = 0.0,
    ax=None,
):
    """Family C: apparent erosion-rate ratio against true mean erosion rate.

    The same samples as :func:`plot_ratio_vs_steptime` with the roles traded:
    the basin's true mean erosion rate moves to the x axis and time since the
    step becomes the colour and size encoding. This is the view that shows
    which erosion rates the transient bias is worst at.

    Args:
        df: One joined sample table. Must carry ``upstream_mean_erate_mm_yr``,
            ``sc_steptime_kyr`` and ``ratio_col``.
        ratio_col: Which apparent-rate ratio to draw. See :data:`RATIO_LABELS`.
        figsize: Figure size, used only when creating a new figure.
        cmap: Colour map for the time encoding.
        vmin: Lower colour limit, or ``None`` to autoscale.
        vmax: Upper colour limit, or ``None`` to autoscale.
        show_colorbar: Draw the colourbar.
        edgecolor: Marker edge colour.
        linewidth: Marker edge width.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``. This function never saves.
    """
    if ratio_col not in RATIO_LABELS:
        raise KeyError(
            f"Unknown ratio column '{ratio_col}'. "
            f"Expected one of {sorted(RATIO_LABELS)}."
        )
    fig, ax, _ = _scatter(
        df,
        x_col=_MEAN_ERATE_COL, y_col=ratio_col, c_col=_STEPTIME_KYR_COL,
        xlabel=_MEAN_ERATE_LABEL, ylabel=RATIO_LABELS[ratio_col],
        cbar_label=_STEPTIME_CBAR, cmap=cmap, vmin=vmin, vmax=vmax,
        show_colorbar=show_colorbar, figsize=figsize, ax=ax,
        edgecolor=edgecolor, linewidth=linewidth,
    )
    return fig, ax


def plot_scenario_delta(
    df: pd.DataFrame,
    *,
    delta_col: str = "d_app_erate_ratio_C14_Be10",
    scenario_col: str = "scenario",
    scenario_order: Optional[Sequence[str]] = None,
    scenario_labels: Optional[dict] = None,
    xscale: str = "symlog",
    linthresh: float = _OVERLAY_LINTHRESH,
    figsize: Tuple[float, float] = (10.0, 5.0),
    show_legend: bool = True,
    ax=None,
):
    """Family E: paired change in the apparent-rate ratio, one series per scenario.

    Each point is one detrital sample differenced against **itself** in the
    no-landslide control, so the control's own sample-to-sample scatter (which
    spans the whole 1.0-1.55 range of the transient) cancels exactly. That is
    what makes a sub-1 % shift visible at all; an unpaired comparison of two
    distributions would be swamped by the scatter.

    A horizontal line at zero is drawn because it is the decision boundary, not
    decoration: above it the perturbation amplifies the transient signal, below
    it the perturbation damps it.

    Args:
        df: Long-form comparison table -- one row per sample per scenario -- as
            written by the experiment's ``scenario_delta.py``. Must carry
            ``sc_steptime_yr``, ``scenario_col`` and ``delta_col``. The delta is
            read from the file; nothing is subtracted here.
        delta_col: Which paired-difference column to draw.
        scenario_col: Column holding the scenario label.
        scenario_order: Draw order, and hence colour order. ``None`` uses the
            order of first appearance, which keeps the figure stable against a
            re-sorted input file.
        scenario_labels: Optional map from scenario id to legend label.
        xscale: ``"symlog"`` (default) or ``"linear"``. As for
            :func:`plot_ratio_overlay_1d`, a log axis is refused because it
            would silently drop the ``t = 0`` anchor samples -- and here those
            are the samples whose delta MUST be zero, so losing them would hide
            the strongest check in the figure.
        linthresh: Symlog linear-region half-width in years.
        figsize: Figure size, used only when creating a new figure.
        show_legend: Draw the scenario legend.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``. This function never saves.
    """
    plt = _import_mpl()
    if xscale not in ("symlog", "linear"):
        raise ValueError(
            f"xscale must be 'symlog' or 'linear', not {xscale!r}. A log axis "
            "would silently drop the t = 0 anchor samples, which are exactly "
            "the ones whose delta must be zero."
        )
    _require_columns(df, (_STEPTIME_YR_COL, scenario_col, delta_col))

    if scenario_order is None:
        scenario_order = list(pd.unique(df[scenario_col].dropna()))

    fig, ax = _new_or_adopted_axes(plt, ax, figsize)

    # Zero first and behind: it is the decision boundary the eye reads against.
    ax.axhline(0.0, color="0.4", linewidth=1.2, zorder=1)

    for name in scenario_order:
        sub = df[df[scenario_col] == name]
        x = pd.to_numeric(sub[_STEPTIME_YR_COL], errors="coerce").to_numpy(float)
        y = pd.to_numeric(sub[delta_col], errors="coerce").to_numpy(float)
        keep = np.isfinite(x) & np.isfinite(y)
        if not keep.any():
            continue
        label = (scenario_labels or {}).get(name, str(name))
        ax.scatter(x[keep], y[keep], s=26, alpha=0.75, edgecolors="none",
                   label=label, zorder=2)

    if xscale == "symlog":
        ax.set_xscale("symlog", linthresh=linthresh)
    ax.set_xlabel(_STEPTIME_LABEL, fontsize=_AXIS_LABEL_SIZE)
    ax.set_ylabel("change in apparent\nerosion-rate ratio", fontsize=_AXIS_LABEL_SIZE)
    ax.tick_params(labelsize=_TICK_LABEL_SIZE)
    if show_legend and len(scenario_order):
        ax.legend(fontsize=_CBAR_LABEL_SIZE, loc="upper left")
    return fig, ax


def plot_ratio_overlay_1d(
    df: pd.DataFrame,
    sweeps: Sequence[Tuple[str, pd.DataFrame]],
    *,
    ratio_col: str = "app_erate_ratio_C14_Be10",
    xscale: str = "symlog",
    linthresh: float = _OVERLAY_LINTHRESH,
    figsize: Tuple[float, float] = (10.0, 5.0),
    cmap: str = _DEFAULT_CMAP,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_colorbar: bool = True,
    sweep_linewidth: float = 2.0,
    show_legend: bool = True,
    legend_loc: str = "upper left",
    ax=None,
):
    """Family D: the spatial samples with the 1-D column sweep drawn over them.

    Both halves of the step-change experiment predict the same thing -- how the
    apparent erosion-rate ratio evolves after the step -- by wholly different
    routes: one routes an 821 km2 landscape and averages each sample's basin,
    the other advances a single analytical column. Drawing them on one axis is
    the check that neither is wrong, because the x variable is the same
    physical quantity in both (``sc_steptime_yr`` for the samples,
    ``time_since_step_yr`` for the sweep).

    **They are not expected to coincide, and the disagreement is the point.**
    A detrital sample averages a whole basin, and every pixel upstream of the
    sample was reached by the knickpoint *later* than the sample itself was --
    many not at all. So a sample plotted at its own step time carries a mixture
    whose transient is younger than that abscissa suggests, and the scatter
    sits to the right of, and below, the column curve. A scatter falling to the
    *left* of the curve, or an anchor sample away from 1.0, would be the
    finding worth chasing.

    Args:
        df: One joined sample table, as for :func:`plot_ratio_vs_steptime`.
        sweeps: Sequence of ``(label, sweep_df)`` pairs, each a
            ``lsdtt-cosmo-column-tool`` step-change sweep table carrying
            ``time_since_step_yr`` and ``ratio_col``. Drawn in order with
            :data:`_SWEEP_STYLES`. The label is the caller's to write: this
            function cannot know whether a given sweep's erosion rates and
            production scaling actually match the spatial run, so it does not
            imply a match that may not hold.
        ratio_col: Which apparent-rate ratio to draw. The **same** column is
            read from the samples and from every sweep, so the two can never be
            crossed. See :data:`RATIO_LABELS`.
        xscale: ``"symlog"`` (default) or ``"linear"``. Symlog resolves both
            the column's rise over the first few thousand years and the
            samples' full range out to ~500 kyr, and unlike a log axis it can
            draw the ``t = 0`` anchor samples. A plain log scale is refused
            rather than offered, because it would drop those samples silently.
        linthresh: Symlog linear-region half-width in years; ignored when
            ``xscale="linear"``.
        figsize: Figure size, used only when creating a new figure.
        cmap: Colour map for the samples' erosion-rate encoding.
        vmin: Lower colour limit, or ``None`` to autoscale.
        vmax: Upper colour limit, or ``None`` to autoscale.
        show_colorbar: Draw the colourbar for the sample encoding.
        sweep_linewidth: Line width for the sweep curves.
        show_legend: Draw a legend naming the sweep curves.
        legend_loc: Legend position. Defaults to ``"upper left"`` rather than
            ``"best"``: the upper left is the one corner this figure is always
            empty in (the column curve is still near 1.0 at small t, and the
            samples sit far to the right), whereas ``"best"`` scores the dense
            sample cloud as free space and lands on top of it.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Raises:
        KeyError: If ``ratio_col`` is unknown, or a sweep lacks the columns.
        ValueError: If ``xscale`` is not ``"symlog"`` or ``"linear"``, or if a
            sweep contributes no finite points.

    Returns:
        ``(fig, ax)``. This function never saves -- see
        :func:`lsdviztools3.workflows.render_cosmo_transient_spatial_figures`.
    """
    if ratio_col not in RATIO_LABELS:
        raise KeyError(
            f"Unknown ratio column '{ratio_col}'. "
            f"Expected one of {sorted(RATIO_LABELS)}."
        )
    if xscale not in ("symlog", "linear"):
        raise ValueError(
            f"xscale must be 'symlog' or 'linear', not {xscale!r}. A log axis "
            "is deliberately not offered: it would silently drop the t = 0 "
            "anchor samples."
        )

    fig, ax, _ = _scatter(
        df,
        x_col=_STEPTIME_YR_COL, y_col=ratio_col, c_col=_MEAN_ERATE_COL,
        xlabel=_STEPTIME_LABEL, ylabel=RATIO_LABELS[ratio_col],
        cbar_label=_MEAN_ERATE_CBAR, cmap=cmap, vmin=vmin, vmax=vmax,
        show_colorbar=show_colorbar, figsize=figsize, ax=ax,
        edgecolor="none", linewidth=0.0,
    )

    for i, (label, sweep_df) in enumerate(sweeps):
        _require_columns(sweep_df, (_SWEEP_TIME_COL, ratio_col))
        x = pd.to_numeric(sweep_df[_SWEEP_TIME_COL], errors="coerce")
        y = pd.to_numeric(sweep_df[ratio_col], errors="coerce")
        keep = np.isfinite(x) & np.isfinite(y)
        if not keep.any():
            raise ValueError(
                f"Sweep '{label}' has no finite '{_SWEEP_TIME_COL}' / "
                f"'{ratio_col}' pairs to draw."
            )
        x = x[keep].to_numpy(dtype=float)
        y = y[keep].to_numpy(dtype=float)
        order = np.argsort(x, kind="stable")
        colour, style = _SWEEP_STYLES[i % len(_SWEEP_STYLES)]
        # zorder above the default scatter zorder (1) so the curve stays
        # readable where the sample cloud is dense.
        ax.plot(
            x[order], y[order], color=colour, linestyle=style,
            linewidth=sweep_linewidth, label=label, zorder=3,
        )

    if xscale == "symlog":
        ax.set_xscale("symlog", linthresh=linthresh)
    if show_legend and len(sweeps):
        if legend_loc.startswith("upper"):
            # The column curve peaks near the top of the axes and the legend
            # box would sit on that peak -- which is the figure's headline
            # number. Reserving a band above the data moves the box off it
            # without moving any drawn value.
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi + _LEGEND_HEADROOM_FRAC * (hi - lo))
        ax.legend(fontsize=_CBAR_LABEL_SIZE, loc=legend_loc)

    return fig, ax
