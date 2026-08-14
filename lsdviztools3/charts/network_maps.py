"""Channel-network map panels for synthetic-network sweeps (matplotlib).

Three panels per case, sharing one layout so a sweep can be read as a strip:

    plot_hillshade          shaded relief
    plot_network            the channel network alone, line width ~ log(area)
    plot_network_on_hillshade   the network drawn over the relief

Each carries a **header above the axes** naming the generation method and its
settings, including the seed for any stochastic method -- a sweep figure is
useless six months later if you cannot tell which panel was which run.

WHY MATPLOTLIB, when the package renders maps with PyGMT
    These are sweep panels, not publication maps, and they need two things PyGMT
    makes awkward: a line width that varies CONTINUOUSLY with drainage area
    (PyGMT would need one plot call per width class), and a multi-line settings
    block outside the plotting area. Keeping all three panels in one engine also
    guarantees they are directly comparable across a sweep, which matters more
    here than GMT's typography. `render/` remains the PyGMT map layer.

NOTHING HERE COMPUTES GEOMORPHOLOGY. Flow directions, drainage area and the
longest-upslope length all come from lsdtt3 rasters. This module turns the D8
direction codes into line segments and draws them.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np

_MPL_MSG = (
    "Network map panels need matplotlib. Install the optional 'charts' extra, "
    "e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

# ArcMap D8 convention, as written by lsdtt-channel-extraction:
# 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE.
# Values are (d_row, d_col) with row increasing SOUTHWARD, matching raster order.
D8_OFFSETS = {
    1: (0, 1),
    2: (1, 1),
    4: (1, 0),
    8: (1, -1),
    16: (0, -1),
    32: (-1, -1),
    64: (-1, 0),
    128: (-1, 1),
}


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def hillshade(z, *, cell_size: float = 1.0, azimuth: float = 315.0,
              altitude: float = 45.0, z_factor: float = 1.0):
    """Horn (1981) hillshade, returned in 0-1.

    Computed here rather than shelled out to GMT so the panels have no GMT
    dependency and the three of them share exactly one relief image.
    """
    z = np.asarray(z, dtype=float) * z_factor
    dzdy, dzdx = np.gradient(z, cell_size, cell_size)
    slope = np.arctan(np.hypot(dzdx, dzdy))
    # atan2 with -dzdx puts aspect in compass terms for a north-up raster.
    aspect = np.arctan2(dzdy, -dzdx)
    az = np.radians(360.0 - azimuth + 90.0)
    alt = np.radians(altitude)
    hs = (np.sin(alt) * np.cos(slope)
          + np.cos(alt) * np.sin(slope) * np.cos(az - aspect))
    return np.clip(hs, 0.0, 1.0)


def network_segments(flowdir, area, *, cell_size: float = 1.0,
                     min_area_cells: float = 0.0):
    """Line segments for every channel link, with a drainage area per segment.

    Args:
        flowdir: ArcMap D8 direction raster (2-D array).
        area: Drainage area raster on the same grid. Units are the caller's --
            only the ORDER matters, because widths use log10.
        cell_size: Used to put segments in map units.
        min_area_cells: Drop links draining less than this, in the same units as
            ``area``. This is what turns "every pixel" into "the channels".

    Returns:
        ``(segments, areas)`` where segments is an (N, 2, 2) array of
        ``[[x0, y0], [x1, y1]]`` in map units with y increasing NORTHWARD, and
        areas is the length-N drainage area of each link.
    """
    fd = np.asarray(flowdir)
    A = np.asarray(area, dtype=float)
    if fd.shape != A.shape:
        raise ValueError(f"flowdir {fd.shape} and area {A.shape} differ in shape")

    nrows, ncols = fd.shape
    keep = np.isfinite(A) & (A >= min_area_cells)

    segs, vals = [], []
    for code, (dr, dc) in D8_OFFSETS.items():
        m = keep & (fd == code)
        if not m.any():
            continue
        r, c = np.nonzero(m)
        r2, c2 = r + dr, c + dc
        # Links that leave the grid are dropped rather than clipped: they are
        # outlets, and a half-length stub at the edge reads as a real channel.
        inside = (r2 >= 0) & (r2 < nrows) & (c2 >= 0) & (c2 < ncols)
        r, c, r2, c2 = r[inside], c[inside], r2[inside], c2[inside]
        if r.size == 0:
            continue
        x0 = c * cell_size
        y0 = (nrows - 1 - r) * cell_size
        x1 = c2 * cell_size
        y1 = (nrows - 1 - r2) * cell_size
        segs.append(np.stack([np.stack([x0, y0], 1),
                              np.stack([x1, y1], 1)], 1))
        vals.append(A[r, c])

    if not segs:
        return np.zeros((0, 2, 2)), np.zeros(0)
    return np.concatenate(segs, 0), np.concatenate(vals, 0)


def _linewidths(areas, *, min_lw: float = 0.15, max_lw: float = 2.6):
    """Line widths from log10(drainage area), scaled into [min_lw, max_lw].

    Log, not linear: drainage area spans orders of magnitude, so a linear width
    makes every channel but the trunk invisible.
    """
    a = np.asarray(areas, dtype=float)
    la = np.log10(np.where(a > 0, a, np.nan))
    lo = np.nanmin(la)
    hi = np.nanmax(la)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.full(a.shape, min_lw)
    t = (la - lo) / (hi - lo)
    return min_lw + t * (max_lw - min_lw)


def _header(fig, ax, title: str, settings: Optional[Mapping[str, object]]):
    """Title plus a settings block ABOVE the axes, never over the map."""
    lines = []
    if settings:
        items = [f"{k} = {v}" for k, v in settings.items()]
        # Wrap to keep the block from running off a narrow figure.
        per_line = 4
        for i in range(0, len(items), per_line):
            lines.append("   ".join(items[i:i + per_line]))
    ax.set_title(title, fontsize=11, fontweight="bold", loc="left", pad=None)
    if lines:
        fig.text(0.01, 0.99, "", va="top")   # anchor
        ax.set_title(title + "\n" + "\n".join(lines),
                     fontsize=9, loc="left", fontfamily="monospace")
        # Restore a bold first line by drawing it separately is fiddly; the
        # monospace block is more legible for a settings dump anyway.
    return ax


def _finish(ax, extent):
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.6)


def _extent(shape, cell_size):
    nrows, ncols = shape
    return [0.0, (ncols - 1) * cell_size, 0.0, (nrows - 1) * cell_size]


def plot_hillshade(z, *, cell_size=1.0, z_factor=1.0, title="", settings=None,
                   figsize=(11.0, 6.0), ax=None):
    """Shaded relief panel.

    ``z_factor`` exists so a sweep can compare surfaces whose relief differs by
    orders of magnitude: a fractal surface with 5 m of relief and an OCN with
    400 m are both legible only if the low one is exaggerated. Always report the
    factor in the header -- an unlabelled exaggeration is a lie about the
    topography.
    """
    plt = _import_mpl()
    fig, ax = (plt.subplots(figsize=figsize) if ax is None else (ax.figure, ax))
    hs = hillshade(z, cell_size=cell_size, z_factor=z_factor)
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1,
              extent=_extent(np.shape(z), cell_size), origin="upper")
    _header(fig, ax, title, settings)
    _finish(ax, _extent(np.shape(z), cell_size))
    fig.tight_layout()
    return fig, ax


def plot_network(flowdir, area, *, cell_size=1.0, min_area_cells=0.0,
                 title="", settings=None, color="#12508c",
                 figsize=(11.0, 6.0), ax=None):
    """Channel network alone, line width scaling with log(drainage area)."""
    plt = _import_mpl()
    from matplotlib.collections import LineCollection

    fig, ax = (plt.subplots(figsize=figsize) if ax is None else (ax.figure, ax))
    segs, areas = network_segments(flowdir, area, cell_size=cell_size,
                                   min_area_cells=min_area_cells)
    ax.set_facecolor("white")
    if segs.shape[0]:
        ax.add_collection(LineCollection(segs, linewidths=_linewidths(areas),
                                         colors=color, capstyle="round"))
    _header(fig, ax, title, settings)
    _finish(ax, _extent(np.shape(flowdir), cell_size))
    fig.tight_layout()
    return fig, ax


def plot_network_on_hillshade(z, flowdir, area, *, cell_size=1.0, z_factor=1.0,
                              min_area_cells=0.0, title="", settings=None,
                              color="#0b3d91", figsize=(11.0, 6.0), ax=None):
    """Channel network drawn over shaded relief. See plot_hillshade re z_factor."""
    plt = _import_mpl()
    from matplotlib.collections import LineCollection

    fig, ax = (plt.subplots(figsize=figsize) if ax is None else (ax.figure, ax))
    ext = _extent(np.shape(z), cell_size)
    ax.imshow(hillshade(z, cell_size=cell_size, z_factor=z_factor),
              cmap="gray", vmin=0, vmax=1, extent=ext, origin="upper")
    segs, areas = network_segments(flowdir, area, cell_size=cell_size,
                                   min_area_cells=min_area_cells)
    if segs.shape[0]:
        ax.add_collection(LineCollection(segs, linewidths=_linewidths(areas),
                                         colors=color, capstyle="round",
                                         alpha=0.95))
    _header(fig, ax, title, settings)
    _finish(ax, ext)
    fig.tight_layout()
    return fig, ax
