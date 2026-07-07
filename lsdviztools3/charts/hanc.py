"""Height-above-channel (HANC) profile plotting (matplotlib).

:func:`plot_hanc_distance` draws HANC against distance along the channel for
classified floodplain / terrace pixels, from an ``lsdtt-flat-surface-finder``
per-pixel table (``*_hanc_profile.csv``, written by ``write_hanc_profile_csv``).
Optionally overlays the reach-relative band break levels from the companion
``*_hanc_band_breaks.csv`` table.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Mapping, Optional

import numpy as np
import pandas as pd

_MPL_MSG = (
    "HANC-profile plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Surface-class codes in the lsdtt3 ``_surface_class`` convention.
_DEFAULT_CLASS_LABELS = {1: "floodplain", 2: "terrace"}
_DEFAULT_CLASS_COLORS = {1: "tab:blue", 2: "tab:orange"}


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def plot_hanc_distance(
    df: pd.DataFrame,
    *,
    distance_col: str = "channel_flow_distance_m",
    hanc_col: str = "hanc_m",
    class_col: str = "class",
    source_col: Optional[str] = "source_id",
    source_id: Optional[int] = None,
    class_labels: Mapping[int, str] = _DEFAULT_CLASS_LABELS,
    class_colors: Mapping[int, str] = _DEFAULT_CLASS_COLORS,
    breaks: Optional[pd.DataFrame] = None,
    distance_units: str = "km",
    point_size: float = 2.0,
    alpha: float = 0.25,
    max_points: int = 200_000,
    seed: int = 0,
    log_y: bool = False,
    ax=None,
):
    """Plot HANC vs distance along the channel, coloured by surface class.

    Args:
        df: Per-pixel HANC profile table (``*_hanc_profile.csv``) with at
            least the distance, HANC, and class columns.
        distance_col: Along-channel flow distance column (metres, measured
            from the DEM outlet).
        hanc_col: Height-above-channel column (metres).
        class_col: Integer surface-class column (1 floodplain, 2 terrace).
        source_col: Channel source-id column, used only when ``source_id``
            filters to one channel. ``None`` pools all sources.
        source_id: Plot only pixels attributed to this channel source.
        class_labels: Legend text per class code.
        class_colors: Point colour per class code.
        breaks: Optional band-breaks table (``*_hanc_band_breaks.csv`` with
            ``window_centre_m`` and semicolon-separated ``breaks_m``) to
            overlay the reach-relative break levels.
        distance_units: ``"km"`` (divide metres by 1000) or ``"m"``.
        point_size: Scatter marker size.
        alpha: Scatter marker opacity.
        max_points: Random subsample cap per class (pixel tables can run to
            millions of rows). 0 or negative = no cap.
        seed: Subsample RNG seed (reproducible figures).
        log_y: Logarithmic HANC axis (drops non-positive values).
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``.

    Raises:
        ValueError: if no plottable pixels remain after filtering.
    """
    plt = _import_mpl()

    data = df
    if source_id is not None and source_col is not None:
        data = data[data[source_col] == source_id]

    scale = 1000.0 if distance_units == "km" else 1.0

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.figure

    rng = np.random.default_rng(seed)
    n_drawn = 0
    for code, label in class_labels.items():
        cls = data[data[class_col] == code]
        x = cls[distance_col].to_numpy(dtype=float)
        y = cls[hanc_col].to_numpy(dtype=float)
        keep = np.isfinite(x) & np.isfinite(y)
        if log_y:
            keep &= y > 0
        x, y = x[keep], y[keep]
        if x.size == 0:
            continue
        if 0 < max_points < x.size:
            idx = rng.choice(x.size, size=max_points, replace=False)
            x, y = x[idx], y[idx]
        ax.scatter(
            x / scale, y,
            s=point_size, alpha=alpha, linewidths=0,
            color=class_colors.get(code, "0.5"), label=label,
            rasterized=True,
        )
        n_drawn += x.size

    if n_drawn == 0:
        raise ValueError("No plottable pixels in the HANC profile table.")

    if breaks is not None and not breaks.empty:
        bx, by = [], []
        brk = breaks
        if source_id is not None and "source_id" in brk.columns:
            brk = brk[brk["source_id"] == source_id]
        for _, row in brk.iterrows():
            levels = str(row.get("breaks_m", "") or "")
            for lvl in levels.split(";"):
                if not lvl:
                    continue
                bx.append(float(row["window_centre_m"]) / scale)
                by.append(float(lvl))
        if bx:
            ax.scatter(bx, by, marker="_", s=120, color="0.15",
                       linewidths=1.2, label="band break", zorder=3)

    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel(f"Distance along channel ({distance_units})")
    ax.set_ylabel("Height above channel (m)")
    ax.margins(x=0)
    legend = ax.legend(fontsize=8, loc="best", frameon=False,
                       scatterpoints=1, markerscale=4)
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)
    return fig, ax
