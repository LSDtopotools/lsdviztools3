"""Rose / wind-rose diagrams (matplotlib).

:func:`rose_diagram` draws a polar rose of bearings (degrees), optionally grouped
(e.g. channel bearings by Strahler stream order) as stacked petals or overlaid
outlines. Geographic convention: North-up, clockwise. Bearings may be treated as
``directional`` (full 360 deg, e.g. upstream->downstream channel links) or
``axial`` (orientation only, folded to 180 deg and mirrored, e.g. lineaments).

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence, Union

import numpy as np

_MPL_MSG = (
    "Rose diagrams need matplotlib. Install the optional 'charts' extra, e.g. "
    "`pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def _clean_bearings(values, mode: str) -> np.ndarray:
    """Finite bearings wrapped to [0, 360); folded to [0, 180) if axial."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    v = np.mod(v, 360.0)
    if mode == "axial":
        v = np.mod(v, 180.0)
    return v


def rose_diagram(
    bearings: Optional[Sequence[float]] = None,
    *,
    groups: Optional[Mapping[Union[str, int], Sequence[float]]] = None,
    bin_width_deg: float = 15.0,
    mode: str = "directional",
    density: bool = False,
    stacked: bool = True,
    cmap: str = "viridis",
    colors: Optional[Mapping] = None,
    label_fmt: str = "{}",
    title: Optional[str] = None,
    ax=None,
):
    """Draw a rose diagram of bearings on a North-up, clockwise polar axis.

    Provide either a flat ``bearings`` array (single rose) or ``groups``
    (``{label: bearings}``) to compare/stack sub-populations such as stream
    orders.

    Args:
        bearings: 1-D bearings in degrees. Ignored if ``groups`` is given.
        groups: Mapping ``label -> bearings`` for a grouped rose.
        bin_width_deg: Angular bin width in degrees (default 15).
        mode: ``"directional"`` (full 360 deg) or ``"axial"`` (fold to 180 deg
            and mirror each petal). Channel/outlet bearings are directional.
        density: Normalise each group to fractions (compare shapes) instead of
            raw counts.
        stacked: With groups, stack petals (low order at the base) vs. overlay
            outline-only petals.
        cmap: Colormap sampled for per-group colours.
        colors: Optional explicit ``{label: colour}`` overriding ``cmap``.
        label_fmt: Format string applied to each group key for the legend.
        title: Axis title.
        ax: Existing polar ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)`` with ``ax`` a North-up, clockwise ``PolarAxes``.

    Raises:
        ValueError: if no bearings are supplied / none are finite, or ``mode``
            is invalid.
    """
    plt = _import_mpl()
    if mode not in ("directional", "axial"):
        raise ValueError(f"mode must be 'directional' or 'axial', got {mode!r}.")

    if groups is None:
        if bearings is None:
            raise ValueError("Provide either `bearings` or `groups`.")
        groups = {"": bearings}

    cleaned = {label: _clean_bearings(v, mode) for label, v in groups.items()}
    cleaned = {label: v for label, v in cleaned.items() if v.size}
    if not cleaned:
        raise ValueError("No finite bearings to plot.")

    # Bins on [0, 180) for axial (mirrored at draw time), else [0, 360).
    top = 180.0 if mode == "axial" else 360.0
    edges = np.arange(0.0, top + bin_width_deg, bin_width_deg)
    centres = np.deg2rad(0.5 * (edges[:-1] + edges[1:]))
    width = np.deg2rad(bin_width_deg)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.5, 5.5),
                               subplot_kw={"projection": "polar"})
    else:
        fig = ax.figure
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    labels = list(cleaned.keys())
    cmap_obj = plt.get_cmap(cmap)
    n = len(labels)
    palette = {lab: (colors.get(lab) if colors and lab in colors
                     else cmap_obj(i / max(n - 1, 1)))
               for i, lab in enumerate(labels)}

    bottom = np.zeros_like(centres)
    show_legend = any(str(lab) for lab in labels)
    for lab in labels:
        counts, _ = np.histogram(cleaned[lab], bins=edges)
        heights = counts.astype(float)
        if density and heights.sum() > 0:
            heights = heights / heights.sum()
        theta, h, base = centres, heights, bottom
        if mode == "axial":  # mirror petals to the opposite semicircle
            theta = np.concatenate([centres, centres + np.pi])
            h = np.concatenate([heights, heights])
            base = np.concatenate([bottom, bottom])
        legend_label = label_fmt.format(lab) if str(lab) else None
        if stacked:
            ax.bar(theta, h, width=width, bottom=base, color=palette[lab],
                   edgecolor="white", linewidth=0.3, align="center",
                   label=legend_label)
            bottom = bottom + heights
        else:
            ax.bar(theta, h, width=width, bottom=0.0, fill=False,
                   edgecolor=palette[lab], linewidth=1.6, align="center",
                   label=legend_label)

    if title:
        ax.set_title(title)
    if show_legend:
        ax.legend(loc="upper right", bbox_to_anchor=(1.25, 1.1), fontsize=8)
    return fig, ax
