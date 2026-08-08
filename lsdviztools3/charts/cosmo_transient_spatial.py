"""Spatial (per-sample) transient cosmogenic scatter charts (matplotlib).

Draws the two scatter families produced by sampling a step-change landscape
along its channel network and inverting each sample's **basin-averaged**
concentration:

* **Family B** — apparent erosion-rate ratio against time since the knickpoint
  passed, coloured and sized by the basin's true mean erosion rate
  (:func:`plot_ratio_vs_steptime`).
* **Family C** — the same ratio against that true mean erosion rate, coloured
  and sized by time since the step (:func:`plot_ratio_vs_mean_erate`).

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
