"""Hydraulic-break analysis and plotting (matplotlib).

The valley cross-section workflow (``lsdtt-valley-centreline`` with
``write_cross_section_features``) floods every cross section stage by stage and
records the *hydraulic breaks*: the water-surface elevations that separate the
flat surfaces (floodplain, inset terraces) from the risers between them. Each
break is written as a local **HAND** value (height above the section's thalweg)
in ``*_valley_cross_section_hydraulic_breaks.csv``.

This module:

* :func:`plot_hydraulic_breaks_downstream` — the break HAND of every section
  against distance downstream, one series per break rank, so floodplain and
  terrace levels can be read as they track down-valley.
* :func:`hydraulic_break_clusters` / :func:`plot_hydraulic_break_clusters` — a
  binned (default 10 sections) 1-D clustering test asking whether the break
  elevations fall into distinct levels within each reach.
* :func:`plot_cross_section_informed_valley_mask` — the new,
  cross-section-informed valley floor (thresholded at the channel's median
  lowest hydraulic break) drawn semi-transparent over a hillshade so it blends
  with the original histogram-derived valley mask where they overlap.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean; the clustering falls back to a pure-numpy
implementation when scikit-learn / scipy are absent.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np
import pandas as pd

_MPL_MSG = (
    "Hydraulic-break plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

_RASTERIO_MSG = (
    "The valley-mask map needs rasterio to read the rasters. Install it e.g. "
    "`pip install rasterio`."
)

#: Fixed, colour-blind-safe series order (Wong palette) for the break ranks.
#: break_index 1 uses the first colour, break_index 2 the second, etc.
_BREAK_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
                 "#D55E00", "#F0E442"]

#: Overlay colours for the cross-section-informed valley mask.
# Original valley floor vs informed floodplain: two saturated, well-separated
# hues drawn semi-transparent so their overlap blends to a clear third colour
# (crimson + blue -> purple).
_ORIGINAL_MASK_COLOR = "#d1495b"   # crimson: original valley floor
_INFORMED_L1_COLOR = "#1f6fb2"     # blue: at/below the median lowest break


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt


def _import_rasterio():
    try:
        import rasterio  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without rio
        raise ImportError(_RASTERIO_MSG) from exc
    return rasterio


def _break_color(break_index: int):
    return _BREAK_COLORS[(int(break_index) - 1) % len(_BREAK_COLORS)]


# --- 1. breaks vs distance downstream ---------------------------------------

def plot_hydraulic_breaks_downstream(
    breaks_df: pd.DataFrame,
    *,
    source_id: int = 0,
    source_col: str = "source_id",
    along_col: str = "along_distance_m",
    hand_col: str = "hand_m",
    break_col: str = "break_index",
    max_break: Optional[int] = None,
    connect: bool = True,
    ax=None,
):
    """Plot the hydraulic-break HAND of every cross section against distance
    downstream, one series per break rank.

    Reveals whether the floodplain (break 1) and successive terrace levels
    (breaks 2, 3, ...) hold a consistent height above the channel as the valley
    runs downstream, or step / merge along it.

    Args:
        breaks_df: Hydraulic-break table
            (``*_valley_cross_section_hydraulic_breaks.csv``) with the source,
            along-distance, ``break_index`` and ``hand_m`` columns. The thalweg
            row (``break_index == 0``, ``hand_m == 0``) is ignored.
        source_id: Channel source to plot.
        source_col, along_col, hand_col, break_col: Column names in
            ``breaks_df``.
        max_break: Highest break rank to draw (``None`` = all present).
        connect: Draw a faint line joining each rank's points down-valley.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``.

    Raises:
        ValueError: if there are no breaks (rank >= 1) for ``source_id``.
    """
    plt = _import_mpl()

    d = breaks_df[breaks_df[source_col] == source_id]
    d = d[d[break_col] >= 1]
    d = d[np.isfinite(d[hand_col].to_numpy())]
    if not len(d):
        raise ValueError(
            f"No hydraulic breaks (break_index >= 1) for "
            f"{source_col}={source_id}.")

    ranks = sorted(int(r) for r in d[break_col].unique())
    if max_break is not None:
        ranks = [r for r in ranks if r <= max_break]

    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 5))
    else:
        fig = ax.figure

    for r in ranks:
        g = d[d[break_col] == r].sort_values(along_col)
        x = g[along_col].to_numpy() / 1000.0     # metres -> km
        y = g[hand_col].to_numpy()
        colour = _break_color(r)
        if connect:
            ax.plot(x, y, "-", lw=0.8, color=colour, alpha=0.35, zorder=2)
        ax.plot(x, y, "o", ms=4.0, color=colour, zorder=3,
                label=f"break {r}")

    ax.set_xlabel("distance downstream (km)")
    ax.set_ylabel("break height above thalweg, HAND (m)")
    ax.grid(True, color="0.9", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    if len(ranks) >= 2:
        ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.8)
    ax.set_title(
        f"channel {source_id}: hydraulic breaks vs distance downstream",
        fontsize=11)
    return fig, ax


# --- 2. binned clustering test ----------------------------------------------

def _kmeans_1d(x: np.ndarray, k: int, iters: int = 100):
    """Deterministic 1-D k-means. Centres are initialised at evenly spaced
    quantiles of the data, so the result is reproducible. Returns
    ``(labels, centres)`` with ``centres`` ascending and ``labels`` aligned to
    the INPUT order of ``x``."""
    xv = np.asarray(x, dtype=float)
    q = (np.arange(k) + 0.5) / k
    centres = np.unique(np.quantile(xv, q))
    # Nudge apart if quantiles collide (many tied values).
    if len(centres) < k:
        lo, hi = float(xv.min()), float(xv.max())
        centres = np.linspace(lo, hi, k) if hi > lo else np.full(k, lo)
    labels = np.zeros(len(xv), dtype=int)
    for _ in range(iters):
        labels = np.argmin(np.abs(xv[:, None] - centres[None, :]), axis=1)
        new = np.array([xv[labels == j].mean() if np.any(labels == j)
                        else centres[j] for j in range(len(centres))])
        if np.allclose(new, centres):
            centres = new
            break
        centres = new
    order = np.argsort(centres)
    remap = np.zeros(len(centres), dtype=int)
    remap[order] = np.arange(len(centres))
    return remap[labels], centres[order]


def _silhouette_1d(x: np.ndarray, labels: np.ndarray) -> float:
    """Mean silhouette coefficient for a 1-D clustering (numpy only). Returns 0
    when fewer than two non-empty clusters exist."""
    xs = np.asarray(x, dtype=float)
    uniq = [u for u in np.unique(labels) if np.any(labels == u)]
    if len(uniq) < 2:
        return 0.0
    sil = np.zeros(len(xs))
    for i in range(len(xs)):
        same = labels == labels[i]
        same[i] = False
        a = (np.abs(xs[same] - xs[i]).mean() if np.any(same) else 0.0)
        b = np.inf
        for u in uniq:
            if u == labels[i]:
                continue
            other = labels == u
            if np.any(other):
                b = min(b, np.abs(xs[other] - xs[i]).mean())
        denom = max(a, b)
        sil[i] = (b - a) / denom if denom > 0 else 0.0
    return float(sil.mean())


def hydraulic_break_clusters(
    breaks_df: pd.DataFrame,
    *,
    source_id: int = 0,
    bin_size: int = 10,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    along_col: str = "along_distance_m",
    hand_col: str = "hand_m",
    break_col: str = "break_index",
    max_clusters: int = 3,
    min_silhouette: float = 0.5,
    min_separation_m: float = 1.0,
) -> pd.DataFrame:
    """Test, in bins of consecutive cross sections, whether the hydraulic-break
    heights fall into distinct clusters (levels).

    The sections of ``source_id`` are ordered by ``station_seq`` and grouped
    into bins of ``bin_size`` consecutive sections. The break HAND values of a
    bin (all ranks >= 1 pooled) are clustered by 1-D k-means for
    ``k = 1..max_clusters``; the ``k`` with the best mean silhouette is chosen,
    but ``k > 1`` is accepted only if that silhouette is at least
    ``min_silhouette`` AND the closest pair of cluster centres is at least
    ``min_separation_m`` apart — so a single diffuse level is not split into
    spurious clusters. Distinct clusters across a bin indicate separate
    floodplain / terrace levels in that reach.

    Uses scikit-learn's ``KMeans``/``silhouette_score`` when importable,
    otherwise a deterministic pure-numpy k-means with the same criterion.

    Args:
        breaks_df: Hydraulic-break table (see
            :func:`plot_hydraulic_breaks_downstream`).
        source_id: Channel source to test.
        bin_size: Number of consecutive sections per bin.
        source_col, station_col, along_col, hand_col, break_col: Column names.
        max_clusters: Largest number of clusters considered per bin.
        min_silhouette: Silhouette floor to accept more than one cluster.
        min_separation_m: Minimum spacing (m) between adjacent cluster centres
            to accept a split.

    Returns:
        One row per bin with columns ``bin_index``, ``bin_centre_along_m``,
        ``along_lo_m``, ``along_hi_m``, ``n_sections``, ``n_values``,
        ``n_clusters``, ``cluster_centres`` (ascending list of HAND values) and
        ``silhouette``.

    Raises:
        ValueError: if there are no breaks (rank >= 1) for ``source_id``.
    """
    d = breaks_df[breaks_df[source_col] == source_id]
    d = d[d[break_col] >= 1]
    d = d[np.isfinite(d[hand_col].to_numpy())]
    if not len(d):
        raise ValueError(
            f"No hydraulic breaks (break_index >= 1) for "
            f"{source_col}={source_id}.")

    stations = sorted(d[station_col].unique())
    rows = []
    for b, lo in enumerate(range(0, len(stations), max(1, bin_size))):
        bin_stations = stations[lo:lo + bin_size]
        sub = d[d[station_col].isin(bin_stations)]
        vals = sub[hand_col].to_numpy(dtype=float)
        along = sub[along_col].to_numpy(dtype=float)
        n_clusters, centres, sil = _cluster_1d(
            vals, max_clusters, min_silhouette, min_separation_m)
        rows.append({
            "bin_index": b,
            "bin_centre_along_m": float(np.mean(along)) if len(along) else np.nan,
            "along_lo_m": float(np.min(along)) if len(along) else np.nan,
            "along_hi_m": float(np.max(along)) if len(along) else np.nan,
            "n_sections": len(bin_stations),
            "n_values": int(len(vals)),
            "n_clusters": int(n_clusters),
            "cluster_centres": [round(float(c), 3) for c in centres],
            "silhouette": float(sil),
        })
    return pd.DataFrame(rows)


def _cluster_1d(vals, max_clusters, min_silhouette, min_separation_m):
    """Choose a cluster count for a 1-D sample. Returns
    ``(n_clusters, centres, silhouette)``."""
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return 0, [], 0.0
    if len(vals) < 2:
        return 1, [float(vals[0])], 0.0

    # Prefer scikit-learn when available; identical criterion either way.
    try:  # pragma: no cover - depends on optional sklearn
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score

        def fit(k):
            km = KMeans(n_clusters=k, n_init=10, random_state=0).fit(
                vals.reshape(-1, 1))
            lab = km.labels_
            cen = km.cluster_centers_.ravel()
            sil = (silhouette_score(vals.reshape(-1, 1), lab)
                   if len(np.unique(lab)) > 1 else 0.0)
            return lab, cen, sil
    except Exception:
        def fit(k):
            lab, cen = _kmeans_1d(vals, k)
            sil = _silhouette_1d(vals, lab)
            return lab, cen, sil

    best_k, best_centres, best_sil = 1, [float(vals.mean())], 0.0
    kmax = min(max_clusters, len(np.unique(vals)))
    for k in range(2, kmax + 1):
        lab, cen, sil = fit(k)
        cen = np.sort(cen)
        # Reject a split that is not well separated (silhouette) or whose
        # closest centres are too near to be distinct levels.
        sep = np.min(np.diff(cen)) if len(cen) > 1 else np.inf
        if sil >= min_silhouette and sep >= min_separation_m and sil > best_sil:
            best_k, best_centres, best_sil = k, list(cen), sil
    return best_k, [float(c) for c in best_centres], float(best_sil)


def plot_hydraulic_break_clusters(
    breaks_df: pd.DataFrame,
    *,
    source_id: int = 0,
    bin_size: int = 10,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    along_col: str = "along_distance_m",
    hand_col: str = "hand_m",
    break_col: str = "break_index",
    clusters_df: Optional[pd.DataFrame] = None,
    ax=None,
    **cluster_kwargs,
):
    """Break-HAND scatter down-valley with the per-bin cluster levels overlaid.

    Draws the same points as :func:`plot_hydraulic_breaks_downstream` (faint,
    coloured by break rank) and, for each bin of ``bin_size`` sections, overlays
    the detected cluster centres as horizontal bars spanning the bin's
    along-distance range — so distinct floodplain / terrace levels stand out and
    reaches where the levels merge are visible.

    Args:
        breaks_df: Hydraulic-break table.
        source_id: Channel source to plot.
        bin_size: Sections per bin (passed to :func:`hydraulic_break_clusters`).
        source_col, station_col, along_col, hand_col, break_col: Column names.
        clusters_df: Pre-computed :func:`hydraulic_break_clusters` output; when
            ``None`` it is computed here (honouring ``cluster_kwargs``).
        ax: Existing matplotlib ``Axes`` or ``None``.
        **cluster_kwargs: Forwarded to :func:`hydraulic_break_clusters`
            (``max_clusters``, ``min_silhouette``, ``min_separation_m``).

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    from matplotlib.lines import Line2D

    fig, ax = plot_hydraulic_breaks_downstream(
        breaks_df, source_id=source_id, source_col=source_col,
        along_col=along_col, hand_col=hand_col, break_col=break_col, ax=ax)
    # Fade the raw points so the cluster levels read on top.
    for ln in ax.get_lines():
        ln.set_alpha(min(ln.get_alpha() or 1.0, 0.5))

    if clusters_df is None:
        clusters_df = hydraulic_break_clusters(
            breaks_df, source_id=source_id, bin_size=bin_size,
            source_col=source_col, station_col=station_col,
            along_col=along_col, hand_col=hand_col, break_col=break_col,
            **cluster_kwargs)

    for _, r in clusters_df.iterrows():
        lo = r["along_lo_m"] / 1000.0
        hi = r["along_hi_m"] / 1000.0
        for c in r["cluster_centres"]:
            ax.plot([lo, hi], [c, c], "-", lw=3.0, color="0.15",
                    solid_capstyle="butt", zorder=5)

    handles = ax.get_legend().legend_handles if ax.get_legend() else []
    handles = list(handles) + [
        Line2D([0], [0], color="0.15", lw=3.0, label="reach cluster level")]
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True,
              framealpha=0.8)
    ax.set_title(
        f"channel {source_id}: hydraulic-break clusters "
        f"(bins of {bin_size} sections)", fontsize=11)
    return fig, ax


# --- 3. cross-section-informed valley mask ----------------------------------

def plot_cross_section_informed_valley_mask(
    informed_mask_path: str,
    original_mask_path: str,
    hillshade_path: str,
    *,
    nodes_df: Optional[pd.DataFrame] = None,
    source_id: Optional[int] = None,
    source_col: str = "source_id",
    x_col: str = "x",
    y_col: str = "y",
    alpha: float = 0.55,
    centreline_color: str = "crimson",
    crop_margin_frac: float = 0.06,
    title: Optional[str] = None,
    ax=None,
):
    """Draw the cross-section-informed valley floor over a hillshade, beside the
    original histogram-derived valley mask.

    The original valley floor and the informed floodplain are each drawn as a
    semi-transparent flat colour at the SAME opacity, so their overlap blends
    into a clear third colour: original-only, informed-only, and both are all
    distinguishable, and the original stays visible everywhere.

    Args:
        informed_mask_path: ``*_cross_section_informed_valley_mask.tif``
            (Byte; 0 = above the break, 1 = at/below the channel's median lowest
            hydraulic break — the floodplain).
        original_mask_path: The original valley mask
            (``*_valley_mask.tif``; > 0 inside the valley).
        hillshade_path: Hillshade GeoTIFF drawn in grey behind.
        nodes_df: Optional valley-centreline nodes to overlay.
        source_id: If given (with ``nodes_df``), only that channel's centreline
            is drawn. ``None`` draws all.
        source_col, x_col, y_col: Column names in ``nodes_df``.
        alpha: Opacity of the mask overlays.
        centreline_color: Colour of the overlaid centreline(s).
        crop_margin_frac: Fraction of the informed-mask extent added as margin.
        title: Figure title; a default is used when None.
        ax: Existing matplotlib ``Axes`` or ``None``.

    Returns:
        ``(fig, ax)``.

    Raises:
        ValueError: if the informed mask has no valley (>0) pixels.
    """
    plt = _import_mpl()
    rasterio = _import_rasterio()
    from matplotlib.colors import to_rgb
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    with rasterio.open(informed_mask_path) as src:
        im = src.read(1).astype(float)
        itr = src.transform
        ib = src.bounds
        ind = src.nodata
    if ind is not None:
        im[im == ind] = 0.0
    informed = im > 0
    if not informed.any():
        raise ValueError(
            f"Informed mask {informed_mask_path} has no valley (>0) pixels.")

    with rasterio.open(original_mask_path) as src:
        om = src.read(1).astype(float)
        ob = src.bounds
        ond = src.nodata
    if ond is not None:
        om[om == ond] = 0.0

    with rasterio.open(hillshade_path) as src:
        hs = src.read(1).astype(float)
        hb = src.bounds
        hnd = src.nodata
    if hnd is not None:
        hs[hs == hnd] = np.nan

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 9))
    else:
        fig = ax.figure

    ax.imshow(hs, cmap="gray", extent=[hb.left, hb.right, hb.bottom, hb.top],
              origin="upper", interpolation="nearest", zorder=0)

    # Original valley floor and the informed floodplain are each drawn as a
    # semi-transparent flat colour, at the SAME alpha, so where they overlap the
    # two blend into a clear third colour (original-only, informed-only, and
    # both are all distinguishable).
    oc = to_rgb(_ORIGINAL_MASK_COLOR)
    orgba = np.zeros((om.shape[0], om.shape[1], 4), dtype=float)
    orgba[om > 0] = (oc[0], oc[1], oc[2], alpha)
    ax.imshow(orgba, extent=[ob.left, ob.right, ob.bottom, ob.top],
              origin="upper", interpolation="nearest", zorder=1)

    ic = to_rgb(_INFORMED_L1_COLOR)
    irgba = np.zeros((im.shape[0], im.shape[1], 4), dtype=float)
    irgba[im == 1] = (ic[0], ic[1], ic[2], alpha)
    ax.imshow(irgba, extent=[ib.left, ib.right, ib.bottom, ib.top],
              origin="upper", interpolation="nearest", zorder=2)

    if nodes_df is not None and len(nodes_df):
        nd = nodes_df
        if source_id is not None and source_col in nd.columns:
            nd = nd[nd[source_col] == source_id]
        if source_col in nd.columns:
            for _, g in nd.groupby(source_col):
                ax.plot(g[x_col].to_numpy(), g[y_col].to_numpy(), "-",
                        lw=1.4, color=centreline_color, zorder=3)
        else:
            ax.plot(nd[x_col].to_numpy(), nd[y_col].to_numpy(), "-",
                    lw=1.4, color=centreline_color, zorder=3)

    # Crop to the informed mask with a margin.
    rows, cols = np.where(informed)
    xs, ys = itr * (cols + 0.5, rows + 0.5)
    xmin, xmax = float(np.min(xs)), float(np.max(xs))
    ymin, ymax = float(np.min(ys)), float(np.max(ys))
    mx = crop_margin_frac * (xmax - xmin + 1e-6)
    my = crop_margin_frac * (ymax - ymin + 1e-6)
    ax.set_xlim(xmin - mx, xmax + mx)
    ax.set_ylim(ymin - my, ymax + my)
    ax.set_aspect("equal")
    ax.set_axis_off()

    # Legend swatches show the actual on-map colours: each layer composited over
    # white, and the two composited together for the overlap ("both").
    def _over(fg, a, bg):
        return tuple(a * f + (1.0 - a) * b for f, b in zip(fg, bg))
    white = (1.0, 1.0, 1.0)
    orig_swatch = _over(oc, alpha, white)
    inf_swatch = _over(ic, alpha, white)
    both_swatch = _over(ic, alpha, _over(oc, alpha, white))
    handles = [
        Patch(facecolor=orig_swatch, label="original valley floor only"),
        Patch(facecolor=inf_swatch, label="informed floodplain only"),
        Patch(facecolor=both_swatch, label="both (overlap)"),
    ]
    if nodes_df is not None and len(nodes_df):
        handles.append(Line2D([0], [0], color=centreline_color, lw=1.4,
                              label="centreline"))
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True,
              framealpha=0.85)
    ax.set_title(
        title if title is not None
        else "cross-section-informed valley floor vs original", fontsize=11)
    return fig, ax
