"""Distribution comparison plots (matplotlib).

:func:`plot_distributions` overlays normalised step-histograms of several
labelled value series on one axis — for comparing the distribution of a quantity
(junction angles, area ratios, gradients, …) across cases. Optionally takes
``log10`` of the values first (for heavy-tailed quantities like area ratios).

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence, Tuple, Union

import numpy as np

_MPL_MSG = (
    "Distribution plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def _clean_values(values, log_x: bool) -> np.ndarray:
    """Finite values as a float array; if ``log_x`` keep > 0 and take log10."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if log_x:
        v = v[v > 0]
        v = np.log10(v)
    return v


def plot_distributions(
    series_by_label: Mapping[str, Sequence[float]],
    *,
    bins: Union[int, Sequence[float]] = 30,
    value_range: Optional[Tuple[float, float]] = None,
    log_x: bool = False,
    density: bool = True,
    styles: Optional[Mapping[str, dict]] = None,
    xlabel: str = "Value",
    ylabel: str = "Density",
    title: Optional[str] = None,
    linewidth: float = 1.6,
    ax=None,
):
    """Overlay step-histograms of several labelled value series.

    Args:
        series_by_label: Mapping of ``label -> 1-D values``. Empty or all-NaN
            series are skipped. With ``log_x`` each series is reduced to its
            positive values and ``log10``-transformed.
        bins: Bin count (shared edges are computed across all series) or explicit
            bin edges.
        value_range: ``(lo, hi)`` for the shared bin edges when ``bins`` is an
            int (in log10 units if ``log_x``). ``None`` uses the pooled data
            range.
        log_x: Plot ``log10`` of the values (values ``<= 0`` are dropped).
        density: Normalise each histogram to unit area (compare shapes, not
            counts).
        styles: Optional ``label -> matplotlib kwargs`` (e.g.
            ``{"a": {"color": "crimson", "linestyle": "--"}}``) to control each
            series' line. Labels missing here use the default colour cycle.
        xlabel, ylabel, title: Axis labels and title.
        linewidth: Step-histogram line width.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``.

    Raises:
        ValueError: if no series has any plottable values.
    """
    plt = _import_mpl()
    styles = styles or {}

    cleaned = {label: _clean_values(vals, log_x)
               for label, vals in series_by_label.items()}
    cleaned = {label: v for label, v in cleaned.items() if v.size}
    if not cleaned:
        raise ValueError("No plottable values in any series.")

    # Shared bin edges so the overlaid histograms are comparable.
    if np.ndim(bins) == 0:
        if value_range is not None:
            lo, hi = value_range
        else:
            pooled = np.concatenate(list(cleaned.values()))
            lo, hi = float(pooled.min()), float(pooled.max())
            if lo == hi:
                hi = lo + 1.0
        edges = np.linspace(lo, hi, int(bins) + 1)
    else:
        edges = np.asarray(bins, dtype=float)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4.5))
    else:
        fig = ax.figure

    for label, v in cleaned.items():
        ax.hist(v, bins=edges, density=density, histtype="step",
                linewidth=linewidth, label=label, **styles.get(label, {}))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.set_xlim(edges[0], edges[-1])
    ax.legend(fontsize=8)
    return fig, ax
