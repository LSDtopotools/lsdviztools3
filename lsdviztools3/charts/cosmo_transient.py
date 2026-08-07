"""Transient (step-change) cosmogenic time-series charts (matplotlib).

Draws the two figure families produced by the lsdtt3
``lsdtt-cosmo-column-tool`` ``step_change_sweep`` CSVs:

* **Family A** — normalised CRN concentration against time since the step
  change, one series per nuclide (:func:`plot_cosmo_transient_concentration`).
* **Family B** — ratios of apparent erosion rates against time since the step
  change (:func:`plot_cosmo_transient_erate_ratio`).

**Visualisation only.** Every quantity drawn here is a column that already
exists in the CSV: the normalisations, the apparent erosion rates and their
ratios are all computed by lsdtt3. This module never divides, normalises or
inverts anything. If a figure needs a derived quantity, the C++ side should
emit it as a new column.

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

#: Value-column prefixes emitted by lsdtt3 for the Family A chart.
#:
#: * ``approach_frac_`` — ``(ss_old - conc) / (ss_old - ss_new)``, rising 0 -> 1.
#: * ``norm_conc_`` — ``conc / ss_old``, starting at 1.
#: * ``norm_full_`` — ``(conc - ss_new) / (ss_old - ss_new)``, falling 1 -> 0.
CONCENTRATION_PREFIXES = ("approach_frac_", "norm_conc_", "norm_full_")

#: Nuclides in the order they are drawn (so they take tab10 positions 0-2).
NUCLIDES = ("Be10", "Al26", "C14")

#: Mathtext legend labels, matched to :data:`NUCLIDES`.
NUCLIDE_LABELS = (r"$^{10}Be$", r"$^{26}Al$", r"$^{14}C$")

#: Family B series: (column, legend label), in draw order.
ERATE_RATIO_SERIES = (
    ("app_erate_ratio_C14_Be10", r"$^{14}C$ to $^{10}Be$"),
    ("app_erate_ratio_C14_Al26", r"$^{14}C$ to $^{26}Al$"),
)

_TIME_COL = "time_since_step_yr"
_SCALING_COL = "production_scaling"
_RATIO_COL = "erosion_ratio"
_FINAL_ERATE_COL = "eff_erate_new_g_cm2_yr"

_XLABEL = "Time since step change (yrs)"
_AXIS_LABEL_SIZE = 14
_LEGEND_SIZE = 16


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def load_cosmo_transient_table(path: str | Path) -> pd.DataFrame:
    """Load one ``lsdtt-cosmo-column-tool`` step-change sweep CSV.

    Args:
        path: Path to ``cosmo_transient_sweep_<label>.csv``. Leading ``#``
            comment lines are tolerated.

    Returns:
        The table as-is. Nothing is computed, renamed or dropped.
    """
    return pd.read_csv(path, comment="#")


def cosmo_transient_title(*, s_t: float, erosion_ratio: float,
                          final_eps: float) -> str:
    """Build the standard mathtext title shared by both figure families.

    Args:
        s_t: Total production scaling factor (``production_scaling``).
        erosion_ratio: New / old effective erosion rate (``erosion_ratio``).
        final_eps: Post-step effective erosion rate in g cm^-2 yr^-1
            (``eff_erate_new_g_cm2_yr``).

    Returns:
        A raw mathtext string (never usetex).
    """
    return (
        rf"$S_t$ = {s_t:g}, erosion ratio = {erosion_ratio:g}, "
        rf"final $\epsilon$ = {final_eps:g} g cm$^{{-2}}$yr$^{{-1}}$"
    )


def _require_columns(df: pd.DataFrame, required: Sequence[str]) -> None:
    missing = set(required) - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")


def _single_run_label(df: pd.DataFrame) -> None:
    """A title built from mixed run metadata would be a lie, so refuse."""
    if "run_label" not in df.columns:
        return
    labels = pd.unique(df["run_label"].dropna())
    if len(labels) > 1:
        raise ValueError(
            "Expected a single run_label per figure; found "
            f"{sorted(map(str, labels))}. Plot one sweep at a time."
        )


def _constant_scalar(df: pd.DataFrame, column: str) -> float:
    """Read a constant metadata column, refusing genuinely varying values."""
    vals = pd.to_numeric(df[column], errors="coerce")
    vals = vals[np.isfinite(vals)]
    if vals.empty:
        raise ValueError(f"Metadata column '{column}' has no finite values.")
    first = float(vals.iloc[0])
    if not np.allclose(vals.to_numpy(dtype=float), first):
        raise ValueError(
            f"Metadata column '{column}' is not constant across the table; "
            "plot one sweep at a time."
        )
    return first


def _title_from_metadata(df: pd.DataFrame) -> str:
    _require_columns(df, (_SCALING_COL, _RATIO_COL, _FINAL_ERATE_COL))
    return cosmo_transient_title(
        s_t=_constant_scalar(df, _SCALING_COL),
        erosion_ratio=_constant_scalar(df, _RATIO_COL),
        final_eps=_constant_scalar(df, _FINAL_ERATE_COL),
    )


def _finite_xy(df: pd.DataFrame, y_col: str):
    x = pd.to_numeric(df[_TIME_COL], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[y_col], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    order = np.argsort(x[keep], kind="stable")
    return x[keep][order], y[keep][order]


def _new_or_adopted_axes(plt, ax, figsize: Tuple[float, float]):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    return fig, ax


def plot_cosmo_transient_concentration(
    df: pd.DataFrame,
    *,
    value_prefix: str = "approach_frac_",
    xlim: Tuple[float, float] = (0.0, 10000.0),
    show_title: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (5.0, 3.0),
    linewidth: float = 2.0,
    ax=None,
):
    """Family A: normalised CRN concentration against time since step change.

    Three series in the order Be10, Al26, C14, so they take tab10 positions
    0-2 from the default colour cycle (colours are never set explicitly).

    Args:
        df: One step-change sweep table (a single ``run_label``). Must carry
            ``time_since_step_yr`` and ``<value_prefix><nuclide>`` for each of
            Be10, Al26 and C14 — all pre-computed by lsdtt3.
        value_prefix: Which normalisation to draw. Default
            ``"approach_frac_"``: ``(ss_old - conc) / (ss_old - ss_new)``,
            which runs from 0 at the step to ~1, matching the target figures.
            See :data:`CONCENTRATION_PREFIXES` for the alternatives
            (``"norm_conc_"``, ``"norm_full_"``).
        xlim: Hard x limits in years. Default ``(0, 10000)``.
        show_title: Draw the standard metadata title.
        title: Explicit title text, overriding the metadata-built one.
        figsize: Figure size, used only when creating a new figure.
        linewidth: Series line width.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Note:
        No gridlines are drawn. This is a deliberate deviation from the house
        ``ax.grid(True, linestyle=":")`` convention: the published figure
        family this reproduces has none, and all four spines are kept instead.

    Returns:
        ``(fig, ax)``. This function never saves — see
        :func:`lsdviztools3.workflows.render_cosmo_transient_figures`.
    """
    plt = _import_mpl()

    value_cols = [f"{value_prefix}{nuc}" for nuc in NUCLIDES]
    _require_columns(df, [_TIME_COL, *value_cols])
    _single_run_label(df)

    series = []
    for col, label in zip(value_cols, NUCLIDE_LABELS):
        x, y = _finite_xy(df, col)
        if x.size:
            series.append((x, y, label))
    if not series:
        raise ValueError(
            "No finite time / concentration values to plot for prefix "
            f"'{value_prefix}'."
        )

    fig, ax = _new_or_adopted_axes(plt, ax, figsize)
    for x, y, label in series:
        ax.plot(x, y, linewidth=linewidth, label=label)

    ax.set_xlim(*xlim)
    ax.set_xlabel(_XLABEL, fontsize=_AXIS_LABEL_SIZE)
    ax.set_ylabel(
        "Normalised CRN\nconcentration (dimensionless)",
        fontsize=_AXIS_LABEL_SIZE,
    )
    ax.legend(fontsize=_LEGEND_SIZE, loc="best")
    if show_title:
        ax.set_title(title if title is not None else _title_from_metadata(df))
    return fig, ax


def plot_cosmo_transient_erate_ratio(
    df: pd.DataFrame,
    *,
    xlim: Tuple[float, float] = (0.0, 50000.0),
    show_title: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (5.0, 3.0),
    linewidth: float = 2.0,
    ax=None,
):
    """Family B: apparent erosion-rate ratios against time since step change.

    Two series in the order C14/Be10 then C14/Al26, so they take tab10
    positions 0-1 from the default colour cycle. Both ratios are read straight
    from the CSV columns ``app_erate_ratio_C14_Be10`` and
    ``app_erate_ratio_C14_Al26``; nothing is divided here.

    Args:
        df: One step-change sweep table (a single ``run_label``), carrying
            ``time_since_step_yr`` plus the two ratio columns.
        xlim: Hard x limits in years. Default ``(0, 50000)``. Y is autoscaled.
        show_title: Draw the standard metadata title.
        title: Explicit title text, overriding the metadata-built one.
        figsize: Figure size, used only when creating a new figure.
        linewidth: Series line width.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Note:
        No gridlines are drawn — the same deliberate deviation from the house
        ``ax.grid(True, linestyle=":")`` convention as Family A.

    Returns:
        ``(fig, ax)``. This function never saves.
    """
    plt = _import_mpl()

    ratio_cols = [col for col, _ in ERATE_RATIO_SERIES]
    _require_columns(df, [_TIME_COL, *ratio_cols])
    _single_run_label(df)

    series = []
    for col, label in ERATE_RATIO_SERIES:
        x, y = _finite_xy(df, col)
        if x.size:
            series.append((x, y, label))
    if not series:
        raise ValueError("No finite apparent erosion-rate ratios to plot.")

    fig, ax = _new_or_adopted_axes(plt, ax, figsize)
    for x, y, label in series:
        ax.plot(x, y, linewidth=linewidth, label=label)

    ax.set_xlim(*xlim)
    ax.set_xlabel(_XLABEL, fontsize=_AXIS_LABEL_SIZE)
    ax.set_ylabel(
        "Ratio of apparent\nerosion rates (dimensionless)",
        fontsize=_AXIS_LABEL_SIZE,
    )
    ax.legend(fontsize=_LEGEND_SIZE, loc="best")
    if show_title:
        ax.set_title(title if title is not None else _title_from_metadata(df))
    return fig, ax
