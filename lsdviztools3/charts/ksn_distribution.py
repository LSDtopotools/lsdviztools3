"""Distributions of k_sn over a chi fit, and k_sn classes drawn from them.

The k_sn of a chi fit can be sampled two ways, which answer different
questions:

* **channel pixels**: every channel node with a k_sn (the fit-nodes table with
  ``k_sn`` attached by :func:`lsdviztools3.io.chi.load_chi_fit`) is one sample,
  so a segment counts in proportion to its length: how much of the channel
  network is steep;
* **segments**: every fitted segment (the fit-segments table) is one sample,
  whatever its length: how many distinct steep reaches the fit found.

``ksn_colour_classes`` finds k_sn classes (natural breaks, quantiles, equal
intervals in log10 k_sn, or breaks you give) from either distribution, for the
k_sn maps and chi profiles. ``plot_ksn_distributions`` draws both
distributions side by side and, given classes, shows them on both.
"""

from __future__ import annotations

import re
from typing import Optional

import numpy as np
import pandas as pd

#: The two ways of sampling k_sn: ``(name, axis title)``.
BASES = {"pixels": "Channel pixels (length-weighted)",
         "segments": "Fitted segments (one sample each)"}
INK = "#3d3d3a"
#: Default k_sn class limits: schemes find their breaks between these.
KSN_CLASS_MIN = 5.0
KSN_CLASS_MAX = 250.0
_RGB = re.compile(r"^(\d{1,3})/(\d{1,3})/(\d{1,3})$")


def mpl_colour(colour):
    """A :class:`~lsdviztools3.colour_classes.ColourClasses` colour (GMT
    ``r/g/b`` 0-255, or a name) as matplotlib takes it."""
    m = _RGB.match(str(colour).strip())
    if m:
        return tuple(int(x) / 255 for x in m.groups())
    return colour


def ksn_values(frame: pd.DataFrame, removed: Optional[pd.DataFrame] = None) -> np.ndarray:
    """The positive, finite ``k_sn`` of a nodes or segments frame, without
    the rows of removed tributaries (:func:`lsdviztools3.io.chi.select_tributaries`)."""
    if removed is not None and len(removed):
        from ..io.chi import removed_mask
        frame = frame[~removed_mask(frame, removed)]
    v = pd.to_numeric(frame["k_sn"], errors="coerce").to_numpy(dtype=float)
    return v[np.isfinite(v) & (v > 0)]


def _check_basis(basis: str) -> str:
    if basis not in BASES:
        raise ValueError(f"basis must be one of {list(BASES)}, got {basis!r}")
    return basis


def ksn_colour_classes(
    classes,
    *,
    segments: Optional[pd.DataFrame] = None,
    nodes: Optional[pd.DataFrame] = None,
    basis: str = "pixels",
    n_classes: int = 5,
    cmap: str = "viridis",
    class_colors=None,
    color_mode: str = "auto",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    removed: Optional[pd.DataFrame] = None,
    class_min: Optional[float] = KSN_CLASS_MIN,
    class_max: Optional[float] = KSN_CLASS_MAX,
    min_end_fraction: float = 0.001,
):
    """k_sn classes for the maps and profiles, found in one of the two
    distributions.

    Args:
        classes: The breaks in k_sn units, or a scheme (``"natural_breaks"``,
            ``"quantile"``, ``"equal_interval"``) worked out in log10 k_sn.
        segments: The fit-segments table (for ``basis="segments"``).
        nodes: The fit nodes with ``k_sn`` (:func:`lsdviztools3.io.chi.load_chi_fit`;
            for ``basis="pixels"``).
        basis: ``"pixels"`` (default: each channel node one sample, so
            length-weighted) or ``"segments"`` (each fitted segment one sample). Only matters
            for a scheme: given breaks are used as they are.
        n_classes, cmap, class_colors, color_mode, vmin, vmax: as
            :func:`lsdviztools3.colour_classes.make_colour_classes`.
        removed: Tributaries left out of the plots; they are not sampled.
        class_min, class_max: The scheme finds its ``n_classes`` classes only
            between these k_sn (default 5 and 250); below ``class_min`` and at
            or above ``class_max`` are end classes, each made only when at
            least ``min_end_fraction`` (default 0.1%) of the samples fall in
            it. ``None`` for no limit. See
            :func:`lsdviztools3.colour_classes.make_colour_classes`.

    Returns:
        A :class:`lsdviztools3.colour_classes.ColourClasses` (``scheme`` is
        suffixed with the basis, e.g. ``natural_breaks (pixels)``).
    """
    from dataclasses import replace

    from ..colour_classes import make_colour_classes

    frame = nodes if _check_basis(basis) == "pixels" else segments
    if frame is None:
        raise ValueError(f"basis={basis!r} needs the {'nodes' if basis == 'pixels' else 'segments'} table")
    cc = make_colour_classes(ksn_values(frame, removed), classes, n_classes=n_classes, cmap=cmap,
                             class_colors=class_colors, vmin=vmin, vmax=vmax, log_scale=True,
                             color_mode=color_mode, class_min=class_min, class_max=class_max,
                             min_end_fraction=min_end_fraction)
    if cc.scheme != "manual":
        cc = replace(cc, scheme=f"{cc.scheme} ({basis})")
    return cc


def _bin_edges(log_values, breaks, bins: int):
    """Even bins in log10 k_sn, split at the class breaks so no bar straddles one."""
    lo, hi = float(np.min(log_values)), float(np.max(log_values))
    if not hi > lo:
        lo, hi = lo - 0.05, hi + 0.05
    edges = np.linspace(lo, hi, bins + 1)
    if breaks is not None and len(breaks):
        b = np.log10(np.asarray(breaks, dtype=float))
        edges = np.unique(np.concatenate([edges, b[(b > lo) & (b < hi)]]))
    return edges


def plot_ksn_distributions(
    segments: pd.DataFrame,
    nodes: pd.DataFrame,
    *,
    removed: Optional[pd.DataFrame] = None,
    classes=None,
    class_basis: Optional[str] = None,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    bins: int = 40,
    title: str = "",
):
    """The k_sn distribution sampled by channel pixel and by fitted segment,
    side by side on a shared log10 k_sn axis.

    Each panel is a histogram in log10 k_sn as a probability density (bins are
    split at class breaks, so they can differ in width), with its median
    marked. Bars are coloured as the maps colour k_sn: by class when
    ``classes`` is given, with the breaks drawn on both panels and a class
    legend, the panel the classes were found in being marked; otherwise along
    ``cmap`` between ``vmin`` and ``vmax`` (k_sn units; default the data's
    range).

    Args:
        segments: Fit-segments table.
        nodes: Fit nodes with ``k_sn`` (:func:`lsdviztools3.io.chi.load_chi_fit`).
        removed: Tributaries left out; they are not sampled.
        classes: A :class:`~lsdviztools3.colour_classes.ColourClasses`
            (e.g. from :func:`ksn_colour_classes`), or ``None``.
        class_basis: ``"pixels"`` or ``"segments"``: which distribution the
            classes came from (marked on that panel); ``None`` for breaks
            given by hand.

    Returns:
        ``(fig, (ax_pixels, ax_segments))``.
    """
    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors
    from matplotlib.patches import Patch

    from ..colour_limits import colour_limits

    samples = {"pixels": ksn_values(nodes, removed), "segments": ksn_values(segments, removed)}
    if any(v.size == 0 for v in samples.values()):
        raise ValueError("no positive k_sn to plot (all tributaries removed?)")
    logs = {k: np.log10(v) for k, v in samples.items()}
    breaks = list(classes.breaks) if classes is not None else None
    edges = _bin_edges(np.concatenate(list(logs.values())), breaks, bins)
    centres = 0.5 * (edges[:-1] + edges[1:])
    if classes is not None:
        idx = classes.classify(10.0 ** centres).astype(int)
        bar_colours = [mpl_colour(classes.colors[i]) for i in idx]
    else:
        lo, hi, _ = colour_limits(np.concatenate(list(logs.values())), vmin, vmax, log_scale=True)
        cm = plt.get_cmap(cmap)
        norm = mcolors.Normalize(vmin=lo, vmax=hi, clip=True)
        bar_colours = [cm(norm(c)) for c in centres]

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True, constrained_layout=True)
    for ax, key in zip(axes, ("pixels", "segments")):
        dens, _ = np.histogram(logs[key], bins=edges, density=True)
        ax.bar(edges[:-1], dens, width=np.diff(edges), align="edge", color=bar_colours,
               edgecolor="white", linewidth=0.4)
        med = float(np.median(logs[key]))
        ax.axvline(med, color=INK, linewidth=1.0)
        ax.annotate(f"median {10 ** med:.3g}", xy=(med, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(3, -3), textcoords="offset points", ha="left", va="top",
                    fontsize=8, color=INK)
        if breaks:
            for b in breaks:
                ax.axvline(np.log10(b), color=INK, linewidth=0.7, linestyle=(0, (3, 2)))
        name = BASES[key] + f"\nn = {samples[key].size:,}"
        if classes is not None and class_basis == key:
            scheme = classes.scheme.split(" (")[0].replace("_", " ")
            name += f" — {scheme} classes found here"
        ax.set_title(name, fontsize=9, loc="left")
        ax.set_xlabel(r"$\log_{10}\,k_{sn}$")
        ax.grid(axis="y", color="0.9", linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel(r"Probability density (per $\log_{10}$ unit)")
    if classes is not None:
        handles = [Patch(facecolor=mpl_colour(c), label=l) for c, l in zip(classes.colors, classes.labels)]
        fig.legend(handles=handles, title=r"$k_{sn}$ class", loc="outside right center",
                   frameon=False, fontsize=8, title_fontsize=9)
    if title:
        fig.suptitle(title, fontsize=10)
    return fig, tuple(axes)
