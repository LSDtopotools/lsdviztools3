"""Chi-profile plotting (matplotlib).

Two views over ``lsdtt-chi-analysis`` fit output:

* :func:`plot_chi_profile` — a single basin's chi-elevation profile.
* :func:`plot_chi_stacked` — many basins on one axis, each normalised to start
  at chi = 0 and shifted right by a cumulative offset (the "stacked chi plot"
  of lsdviztools), coloured globally by ``k_sn``.

Both accept **either** representation of the fit, so plotting works whether the
analysis was written as CSV or FlatGeobuf:

* a **nodes** table — the CSV ``*_chi_fit_nodes_mn_*.csv`` with per-node
  ``observed_elevation``/``fitted_elevation`` and a ``k_sn`` column attached by
  :func:`lsdviztools3.io.chi.assign_node_ksn`. Drawn as observed points coloured
  by k_sn with the fitted profile overlaid.
* a **segments** table — the CSV ``*_chi_fit_segments_mn_*.csv`` or the fgb
  ``*_chi_fit_segment_lines_mn_*.fgb`` (via
  :func:`lsdviztools3.io.chi.load_chi_segment_lines`). Each fitted segment is
  drawn as a line in (chi, elevation) coloured by its own ``k_sn`` — no observed
  scatter, and no node->segment join needed.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

from ..colour_limits import colour_limits

GHOST_COLOUR = "0.78"   # removed tributaries, drawn under the kept ones


def _split_removed(rows: pd.DataFrame, removed, show_removed: str):
    """(kept rows, removed rows to ghost) for one frame; see plot_chi_profile."""
    if show_removed not in ("hide", "ghost"):
        raise ValueError(f"show_removed must be 'hide' or 'ghost', got {show_removed!r}")
    if removed is None or len(removed) == 0:
        return rows, rows.iloc[0:0]
    from ..io.chi import removed_mask
    m = removed_mask(rows, removed)
    gone = rows[m] if show_removed == "ghost" else rows.iloc[0:0]
    return rows[~m], gone

_MPL_MSG = (
    "Chi plotting needs matplotlib. Install the optional 'charts' extra, "
    "e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
        from matplotlib import colors  # noqa: F401
        from matplotlib.collections import LineCollection  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt, colors, LineCollection


def _is_segments(data: pd.DataFrame) -> bool:
    """A segments table has chi/elevation segment endpoints; nodes do not."""
    return {"downstream_chi", "upstream_chi"}.issubset(data.columns)


def _colour_values(values: np.ndarray, log_scale: bool) -> np.ndarray:
    """Return colour values, optionally log10-scaled (floored like lsdviztools)."""
    v = np.asarray(values, dtype=float)
    if not log_scale:
        return v
    # Mirror lsdviztools: floor tiny/zero k_sn to avoid -inf in log space.
    out = np.full_like(v, np.nan)
    positive = v >= 0.1
    out[positive] = np.log10(v[positive])
    out[(~positive) & np.isfinite(v)] = 0.0
    return out


def _colour_mapping(rows, color_by, log_color, cmap, vmin, vmax, classes):
    """``(values, cmap, norm, extend)`` for colouring ``rows``: continuous
    (log10) k_sn, or the class index when ``classes`` (a
    :class:`~lsdviztools3.colour_classes.ColourClasses`) is given."""
    from matplotlib import colors

    if classes is not None:
        from .ksn_distribution import mpl_colour
        n = classes.n_classes
        return (classes.classify(rows[color_by].to_numpy(dtype=float)),
                colors.ListedColormap([mpl_colour(c) for c in classes.colors]),
                colors.BoundaryNorm(np.arange(n + 1) - 0.5, n), "neither")
    cvals = _colour_values(rows[color_by].to_numpy(), log_color)
    lo, hi, extend = colour_limits(cvals, vmin, vmax, log_scale=log_color)
    return cvals, cmap, colors.Normalize(vmin=lo, vmax=hi), extend


def _add_colorbar(fig, ax, mappable, color_by, log_color, extend, classes, max_height=None):
    """Colourbar beside ``ax`` (an axes or a list of row axes). With
    ``max_height`` (inches) it is no taller than that, at the top of the rows,
    instead of spanning them all."""
    kw = {}
    if max_height is not None:
        axes = list(np.atleast_1d(ax))
        top = max(a.get_position().y1 for a in axes)
        bottom = min(a.get_position().y0 for a in axes)
        span = (top - bottom) * fig.get_size_inches()[1]
        if span > max_height:
            kw = dict(shrink=max_height / span, anchor=(0.0, 1.0), panchor=(1.0, 1.0),
                      fraction=0.15 * min(1.0, 6.0 / fig.get_size_inches()[0]) + 0.02)
    if classes is None:
        fig.colorbar(mappable, ax=ax, label=_colour_label(color_by, log_color), extend=extend, **kw)
        return
    cb = fig.colorbar(mappable, ax=ax, label=r"$k_{sn}$ class" if color_by == "k_sn" else color_by,
                      ticks=np.arange(classes.n_classes), **kw)
    cb.ax.set_yticklabels(list(classes.labels), fontsize=8)
    cb.ax.tick_params(length=0)


def _colour_label(color_by: str, log_scale: bool) -> str:
    if color_by == "k_sn":
        return r"$\log_{10}\,k_{sn}$" if log_scale else r"$k_{sn}$"
    return f"log10 {color_by}" if log_scale else color_by


def _basin_chi_bounds(data: pd.DataFrame, segments: bool):
    """Return (min_chi, max_chi) for a single-basin frame."""
    if segments:
        return float(data["downstream_chi"].min()), float(data["upstream_chi"].max())
    return float(data["chi"].min()), float(data["chi"].max())


def _basin_column(df: pd.DataFrame) -> str:
    """Basin grouping column: ``basin_index`` (fit tables) or ``basin_key`` (data maps)."""
    for name in ("basin_index", "basin_key"):
        if name in df.columns:
            return name
    raise KeyError(
        "Frame has no 'basin_index' or 'basin_key' column; "
        f"columns are {list(df.columns)}."
    )


def _scatter_observed(ax, observed, basin, shift, *, color, size):
    """Overlay observed chi-elevation points for one basin at the given x-shift.

    ``observed`` is a chi data-map frame (e.g. from
    :func:`lsdviztools3.io.chi.load_chi_points`) with ``chi``/``elevation`` and a
    basin column. Coloured neutrally because the data-map points carry no k_sn.
    """
    col = _basin_column(observed)
    pts = observed[observed[col] == basin]
    if pts.empty:
        return
    ax.scatter(
        pts["chi"].to_numpy(dtype=float) + shift, pts["elevation"],
        s=size, c=color, edgecolors="none", zorder=0,
    )


def _segment_lines(data: pd.DataFrame, x_shift: float):
    """Build matplotlib line segments [(chi0,z0),(chi1,z1)] in (chi, elevation)."""
    d_chi = data["downstream_chi"].to_numpy(dtype=float) + x_shift
    u_chi = data["upstream_chi"].to_numpy(dtype=float) + x_shift
    d_z = data["downstream_elevation"].to_numpy(dtype=float)
    u_z = data["upstream_elevation"].to_numpy(dtype=float)
    return [
        [(d_chi[i], d_z[i]), (u_chi[i], u_z[i])] for i in range(len(data))
    ]


def stack_offsets(
    data: pd.DataFrame,
    basins: Optional[Sequence[int]] = None,
    x_offset: Optional[float] = None,
) -> dict:
    """The chi shift of each basin on the stacked plot.

    A basin's x on the stack is ``chi + shift``: each basin starts at zero and
    is placed after the previous one plus ``x_offset`` (default 20% of the
    widest basin's chi span). :func:`plot_chi_stacked` uses exactly these, so
    other layers (e.g. knickpoints) can be placed on the same axis.

    Returns:
        ``{basin_index: shift}`` for the basins present.
    """
    segments = _is_segments(data)
    if basins is None:
        basins = sorted(data["basin_index"].unique())
    bounds = {}
    for b in basins:
        rows = data[data["basin_index"] == b]
        if not rows.empty:
            bounds[b] = _basin_chi_bounds(rows, segments)
    if x_offset is None:
        widths = [hi - lo for lo, hi in bounds.values()]
        x_offset = 0.2 * max(widths) if widths else 0.0
    out, cursor = {}, 0.0
    for b in basins:
        if b not in bounds:
            continue
        lo, hi = bounds[b]
        out[b] = -lo + cursor
        cursor += (hi - lo) + x_offset
    return out


def stack_rows(
    data: pd.DataFrame,
    basins: Optional[Sequence[int]] = None,
    x_offset: Optional[float] = None,
    basins_per_row: Optional[int] = None,
) -> list:
    """The rows of the stacked plot and each basin's chi shift in its row.

    Basins (in the given order) fill rows of at most ``basins_per_row``
    (``None``: one row), e.g. 9 basins at 5 per row make rows of 5 and 4.
    Every row uses the same gap ``x_offset`` (default 20% of the widest basin's
    chi span over all the basins), so the rows share one chi scale.

    Returns:
        A list of ``{basin_index: shift}``, one per row.
    """
    if basins is None:
        basins = sorted(data["basin_index"].unique())
    basins = [b for b in basins if (data["basin_index"] == b).any()]
    if x_offset is None:
        segments = _is_segments(data)
        widths = [np.subtract(*_basin_chi_bounds(data[data["basin_index"] == b], segments)[::-1])
                  for b in basins]
        x_offset = 0.2 * max(widths) if widths else 0.0
    n = max(len(basins), 1) if basins_per_row is None else int(basins_per_row)
    if n < 1:
        raise ValueError(f"basins_per_row must be >= 1, got {basins_per_row}")
    return [stack_offsets(data, basins[i:i + n], x_offset) for i in range(0, len(basins), n)] or [{}]


def reference_ksn(values, n: int = 3) -> list:
    """``n`` k_sn values spanning the data for a gradient inset: the 10th,
    50th and 90th percentiles (for ``n = 3``) rounded to 2 significant figures."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v) & (v > 0)]
    if v.size == 0:
        return []
    out = []
    for q in np.linspace(10, 90, n):
        x = float(np.percentile(v, q))
        x = float(f"{x:.2g}")
        if x > 0 and x not in out:
            out.append(x)
    return out


#: Corners for :func:`add_gradient_inset`; ``"auto"`` picks the emptiest.
INSET_CORNERS = ("upper left", "upper right", "lower left", "lower right")


def _inset_fractions(ax, size, pad, loc):
    """(x0, y0, wf, hf) of an inset of ``size`` inches in corner ``loc``, in axes fraction."""
    bbox = ax.get_position()
    fw, fh = ax.figure.get_size_inches()
    wf = min(0.45, size[0] / (bbox.width * fw))
    hf = min(0.45, size[1] / (bbox.height * fh))
    x0 = pad if "left" in loc else 1 - pad - wf
    y0 = 1 - pad - hf if "upper" in loc else pad
    return x0, y0, wf, hf


def _plotted_points(ax) -> np.ndarray:
    """Points (axes fraction) of what is drawn on ``ax``: scatter offsets,
    line vertices, line-collection segments (sampled along them) and texts."""
    from matplotlib.collections import LineCollection, PathCollection

    pts = []
    for c in ax.collections:
        if isinstance(c, LineCollection):
            for seg in c.get_segments():
                seg = np.asarray(seg, dtype=float)
                if len(seg) >= 2:
                    t = np.linspace(0, 1, 6)[:, None]
                    pts.append(seg[0] + t * (seg[-1] - seg[0]))
        elif isinstance(c, PathCollection):
            pts.append(np.asarray(c.get_offsets(), dtype=float))
    for ln in ax.lines:
        pts.append(np.asarray(ln.get_xydata(), dtype=float))
    for t in ax.texts:
        pts.append(np.atleast_2d(np.asarray(t.get_position(), dtype=float)))
    if not pts:
        return np.empty((0, 2))
    xy = np.concatenate([p.reshape(-1, 2) for p in pts])
    return ax.transLimits.transform(xy[np.all(np.isfinite(xy), axis=1)])


def best_inset_spot(axes, size=(1.3, 1.0), pad: float = 0.02):
    """``(ax, corner)``: the corner, over all ``axes`` (e.g. the rows of a
    stacked plot), whose inset area covers the fewest plotted points (ties go
    to the first row, then the order of :data:`INSET_CORNERS`)."""
    best = None
    for a in axes:
        xy = _plotted_points(a)
        for loc in INSET_CORNERS:
            x0, y0, wf, hf = _inset_fractions(a, size, pad, loc)
            m = 0.03
            n = int(np.sum((xy[:, 0] > x0 - m) & (xy[:, 0] < x0 + wf + m)
                           & (xy[:, 1] > y0 - m) & (xy[:, 1] < y0 + hf + m)))
            if best is None or n < best[0]:
                best = (n, a, loc)
    return best[1], best[2]


def add_gradient_inset(ax, ksn_values, *, loc: str = "auto", size=(1.3, 1.0), pad: float = 0.02):
    """A key of chi-elevation gradients: one line per k_sn value, at the same
    scale as ``ax``.

    In chi-elevation space a channel's slope is its k_sn (dz/dchi, with the
    reference area of the chi run), so the key shows how steep a segment with
    each k_sn looks on this plot. The inset takes ``size`` (inches) in a
    corner (``loc``: one of :data:`INSET_CORNERS`, or ``"auto"`` for the
    corner covering the fewest plotted points), and its axis ranges are the
    same fractions of ``ax``'s ranges, so the scales match. ``ax`` may be a
    list of axes sharing their scales (the rows of a stacked plot): ``"auto"``
    then searches every row's corners, a named corner uses the first row.
    Call it after the limits are final (it freezes them).

    Returns:
        The inset ``Axes``.
    """
    ksn = [float(k) for k in ksn_values if np.isfinite(k) and k > 0]
    if not ksn:
        raise ValueError("give at least one positive k_sn for the gradient inset")
    for a in np.atleast_1d(ax):
        a.set_autoscale_on(False)
    if loc == "auto":
        ax, loc = best_inset_spot(list(np.atleast_1d(ax)), size, pad)
    elif loc not in INSET_CORNERS:
        raise ValueError(f"loc must be 'auto' or one of {INSET_CORNERS}, got {loc!r}")
    else:
        ax = np.atleast_1d(ax)[0]
    x0, y0, wf, hf = _inset_fractions(ax, size, pad, loc)
    (xa, xb), (ya, yb) = ax.get_xlim(), ax.get_ylim()
    w, h = wf * (xb - xa), hf * (yb - ya)
    ins = ax.inset_axes([x0, y0, wf, hf])
    ins.set_xlim(0, w)
    ins.set_ylim(0, h)
    ins.set_xticks([])
    ins.set_yticks([])
    ins.set_facecolor("white")
    for sp in ins.spines.values():
        sp.set_linewidth(0.6)
        sp.set_color("0.4")
    # origin inset a little so the labels fit; lines run to the right/top edge
    ox, oy = 0.08 * w, 0.08 * h
    for k in sorted(ksn):
        run = 0.58 * w
        if oy + k * run > 0.76 * h:
            run = (0.76 * h - oy) / k
        x1, y1 = ox + run, oy + k * run
        ins.plot([ox, x1], [oy, y1], color="#3d3d3a", linewidth=1.2, solid_capstyle="round")
        ins.text(x1 + 0.02 * w, y1, f"{k:g}", fontsize=7, color="#3d3d3a", va="center", ha="left")
    ins.text(0.96 * w, 0.96 * h, r"$k_{sn}$ (dz/d$\chi$)", fontsize=7, color="#3d3d3a",
             ha="right", va="top")
    return ins


def _gradient_inset(ax, gradient_inset, values, loc):
    """Draw the inset for a ``gradient_inset`` argument (False/None, True or
    ``"auto"`` for :func:`reference_ksn` of ``values``, or a list of k_sn)."""
    if gradient_inset is None or gradient_inset is False:
        return None
    if gradient_inset is True or (isinstance(gradient_inset, str) and gradient_inset == "auto"):
        ks = reference_ksn(values)
    else:
        ks = list(gradient_inset)
    return add_gradient_inset(ax, ks, loc=loc) if ks else None


def plot_chi_profile(
    data: pd.DataFrame,
    *,
    basin: Optional[int] = None,
    color_by: str = "k_sn",
    log_color: bool = True,
    cmap: str = "viridis",
    point_size: float = 4.0,
    show_fit: bool = True,
    observed: Optional[pd.DataFrame] = None,
    observed_color: str = "0.6",
    observed_size: float = 2.0,
    ax=None,
    colorbar: bool = True,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    removed: Optional[pd.DataFrame] = None,
    show_removed: str = "ghost",
    classes=None,
    gradient_inset=None,
    gradient_inset_loc: str = "auto",
):
    """Plot a single basin's chi-elevation profile.

    Accepts a **nodes** table (observed points coloured by k_sn + fitted
    overlay) or a **segments** table (fitted segments coloured by their own
    k_sn). See the module docstring.

    Args:
        data: Nodes or segments table for one or more basins.
        basin: ``basin_index`` to plot. If ``None``, the first basin present.
        color_by: Colour column (default ``k_sn``).
        log_color: Colour on a log10 scale (default True; sensible for k_sn).
        cmap: matplotlib colourmap name.
        point_size: Observed-point marker size (nodes table only).
        show_fit: Draw the fitted profile (nodes table only; always drawn for
            a segments table).
        observed: Optional chi data-map frame (e.g. from
            :func:`lsdviztools3.io.chi.load_chi_points` on ``*_chi_data_maps.fgb``)
            whose observed chi-elevation points are overlaid as neutral markers.
            Useful with a segments table to show the data under the fit.
        observed_color: Marker colour for the observed overlay.
        observed_size: Marker size for the observed overlay.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.
        colorbar: Attach a colourbar.
        vmin, vmax: Colour-scale limits in DATA units (k_sn, not log10 k_sn),
            ``None`` for the basin's own range. Values beyond them take the end
            colours and the colourbar shows an extend arrow. Fixed limits make
            separate per-basin plots comparable.
        removed: Tributaries to leave out (the ``removed`` frame of
            :func:`lsdviztools3.io.chi.select_tributaries`, or any frame with
            ``basin_index`` and ``branch_id``). The colour range is taken from
            the kept tributaries only.
        show_removed: ``"ghost"`` (default) draws removed tributaries thin and
            light grey under the kept ones; ``"hide"`` leaves them out.
        classes: Colour k_sn in classes instead of a continuous scale: a
            :class:`~lsdviztools3.colour_classes.ColourClasses` (e.g. from
            :func:`lsdviztools3.charts.ksn_distribution.ksn_colour_classes`,
            so the profiles match a classed map). ``vmin`` / ``vmax`` and
            ``log_color`` are then unused; the colourbar shows the classes.
        gradient_inset: Add a key of chi-elevation gradients at the plot's
            scale (:func:`add_gradient_inset`): ``True`` / ``"auto"`` for
            the 10th, 50th and 90th percentiles of the plotted k_sn, or a list
            of k_sn values; ``None`` (default) for none.
        gradient_inset_loc: Its corner, one of :data:`INSET_CORNERS`, or
            ``"auto"`` (default) for the emptiest corner (of any row).

    Returns:
        ``(fig, ax)``.
    """
    plt, colors, LineCollection = _import_mpl()
    segments = _is_segments(data)

    if basin is None:
        basin = int(data["basin_index"].iloc[0])
    rows = data[data["basin_index"] == basin]
    if rows.empty:
        raise ValueError(f"No data found for basin_index={basin}.")

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = ax.figure

    if observed is not None:
        _scatter_observed(ax, observed, basin, 0.0,
                          color=observed_color, size=observed_size)

    mappable = None
    rows, ghost = _split_removed(rows, removed, show_removed)
    cvals, cmap, norm, extend = _colour_mapping(rows, color_by, log_color, cmap, vmin, vmax, classes)
    if segments and len(ghost):
        ax.add_collection(LineCollection(_segment_lines(ghost, 0.0), colors=GHOST_COLOUR,
                                         linewidths=0.7, zorder=0))
    elif len(ghost):
        ax.scatter(ghost["chi"], ghost["observed_elevation"], c=GHOST_COLOUR,
                   s=point_size * 0.6, edgecolors="none", zorder=0)
    if segments:
        lc = LineCollection(_segment_lines(rows, 0.0), cmap=cmap, norm=norm,
                            linewidths=1.4)
        lc.set_array(cvals)
        ax.add_collection(lc)
        ax.autoscale()
        mappable = lc
    else:
        mappable = ax.scatter(
            rows["chi"], rows["observed_elevation"],
            c=cvals, cmap=cmap, norm=norm, s=point_size, edgecolors="none", zorder=2,
        )
        if show_fit and "fitted_elevation" in rows.columns:
            for _, branch in rows.groupby("branch_id"):
                branch = branch.sort_values("chi")
                ax.plot(branch["chi"], branch["fitted_elevation"],
                        color="0.25", linewidth=0.8, zorder=1)

    ax.set_xlabel(r"$\chi$ (m)")
    ax.set_ylabel("Elevation (m)")
    ax.set_title(f"Basin {basin}")
    if colorbar and mappable is not None:
        _add_colorbar(fig, ax, mappable, color_by, log_color, extend, classes)
    if gradient_inset is not None and gradient_inset is not False:
        ax.autoscale_view()
        _gradient_inset(ax, gradient_inset, rows[color_by].to_numpy(dtype=float), gradient_inset_loc)
    return fig, ax


def plot_chi_stacked(
    data: pd.DataFrame,
    *,
    basins: Optional[Sequence[int]] = None,
    x_offset: Optional[float] = None,
    color_by: str = "k_sn",
    log_color: bool = True,
    cmap: str = "viridis",
    point_size: float = 3.0,
    show_fit: bool = True,
    observed: Optional[pd.DataFrame] = None,
    observed_color: str = "0.6",
    observed_size: float = 1.5,
    label_basins: bool = True,
    ax=None,
    colorbar: bool = True,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    removed: Optional[pd.DataFrame] = None,
    show_removed: str = "ghost",
    classes=None,
    basins_per_row: Optional[int] = None,
    row_height: float = 3.2,
    gradient_inset=None,
    gradient_inset_loc: str = "auto",
):
    """Stacked chi-elevation plot: basins normalised to chi=0 and offset right.

    For each basin the chi axis is shifted to start at zero and then translated
    by a cumulative offset so the profiles sit side-by-side on one axis, as in
    lsdviztools' ``StackedProfilesGradient``. Points/segments are coloured
    globally by ``color_by`` (typically log10 ``k_sn``).

    Accepts a **nodes** or **segments** table (module docstring).

    Args:
        data: Nodes or segments table.
        basins: Ordered ``basin_index`` values to stack. Default: all present.
        x_offset: Horizontal gap (chi units) between successive basins.
            Default: 20% of the widest single-basin chi span.
        color_by: Colour column (default ``k_sn``).
        log_color: Colour on a log10 scale (default True).
        cmap: matplotlib colourmap name.
        point_size: Marker size (nodes table only).
        show_fit: Overlay the fitted profile (nodes table only).
        observed: Optional chi data-map frame (e.g. from
            :func:`lsdviztools3.io.chi.load_chi_points`) whose observed points are
            overlaid as neutral markers, using the same per-basin offset.
        observed_color: Marker colour for the observed overlay.
        observed_size: Marker size for the observed overlay.
        label_basins: Annotate each stack with its ``basin_index``.
        ax: Existing matplotlib ``Axes`` or ``None``.
        colorbar: Attach a colourbar.
        vmin, vmax: Colour-scale limits in DATA units (k_sn, not log10 k_sn),
            ``None`` for the range over the selected basins. Values beyond them
            take the end colours and the colourbar shows an extend arrow.
        removed, show_removed: Tributaries to leave out, and whether to ghost
            them; as :func:`plot_chi_profile`. Basin positions on the stack do
            not change.
        classes: Colour k_sn in classes instead of a continuous scale: a
            :class:`~lsdviztools3.colour_classes.ColourClasses` (e.g. from
            :func:`lsdviztools3.charts.ksn_distribution.ksn_colour_classes`,
            so the profiles match a classed map). ``vmin`` / ``vmax`` and
            ``log_color`` are then unused; the colourbar shows the classes.
        basins_per_row: At most this many basins per row; more wrap onto
            further rows (9 basins at 5 per row: rows of 5 and 4). The rows
            share the chi and elevation scales and one colourbar. ``None``
            (default): one row. Not with ``ax``.
        row_height: Height of each row in inches (default 3.2).
        gradient_inset: Add a key of chi-elevation gradients at the plot's
            scale (:func:`add_gradient_inset`): ``True`` / ``"auto"`` for
            the 10th, 50th and 90th percentiles of the plotted k_sn, or a list
            of k_sn values; ``None`` (default) for none.
        gradient_inset_loc: Its corner, one of :data:`INSET_CORNERS`, or
            ``"auto"`` (default) for the emptiest corner (of any row).
            On a multi-row plot the scales are shared, so it holds for every
            row wherever it goes.

    Returns:
        ``(fig, ax)``; with more than one row, ``ax`` is the array of row axes
        (top first).
    """
    plt, colors, LineCollection = _import_mpl()
    segments = _is_segments(data)

    if basins is None:
        basins = sorted(data["basin_index"].unique())
    if len(basins) == 0:
        raise ValueError("No basins to plot.")

    # Global colour normalisation across all selected basins.
    sel = data[data["basin_index"].isin(basins)]
    sel, _ = _split_removed(sel, removed, show_removed)
    _, cmap, norm, extend = _colour_mapping(sel, color_by, log_color, cmap, vmin, vmax, classes)

    def cvalues(r):
        if classes is not None:
            return classes.classify(r[color_by].to_numpy(dtype=float))
        return _colour_values(r[color_by].to_numpy(), log_color)

    layout = stack_rows(data, basins, x_offset, basins_per_row)
    n_basins_max = max(len(r) for r in layout)

    if ax is None:
        fig, axes = plt.subplots(len(layout), 1, squeeze=False, sharex=True, sharey=True,
                                 figsize=(max(6, 1.6 * n_basins_max), row_height * len(layout) + 0.8))
        axes = list(axes[:, 0])
    elif len(layout) > 1:
        raise ValueError("basins_per_row with more than one row needs ax=None (it makes the figure).")
    else:
        fig, axes = ax.figure, [ax]

    mappable = None
    for row_ax, shifts in zip(axes, layout):
        m = _draw_stack_row(row_ax, data, shifts, segments=segments, removed=removed,
                            show_removed=show_removed, observed=observed,
                            observed_color=observed_color, observed_size=observed_size,
                            cvalues=cvalues, cmap=cmap, norm=norm, point_size=point_size,
                            show_fit=show_fit, label_basins=label_basins,
                            LineCollection=LineCollection)
        mappable = m if m is not None else mappable
        row_ax.set_ylabel("Elevation (m)")
    _share_limits(axes)
    axes[-1].set_xlabel(r"$\chi$ (m, offset per basin)")
    if colorbar and mappable is not None:
        _add_colorbar(fig, axes if len(axes) > 1 else axes[0], mappable, color_by, log_color,
                      extend, classes, max_height=max(row_height, 3.0) if len(axes) > 1 else None)
    if gradient_inset is not None and gradient_inset is not False:
        _gradient_inset(axes, gradient_inset, sel[color_by].to_numpy(dtype=float),
                        gradient_inset_loc)
    return fig, (axes[0] if len(axes) == 1 else np.array(axes))


def _share_limits(axes, margin: float = 0.05):
    """Set every row to the union of the rows' data limits (plus ``margin``),
    explicitly, so the rows share one chi and elevation scale whatever is
    added later."""
    from matplotlib.transforms import Bbox

    boxes = [a.dataLim for a in axes if np.all(np.isfinite(a.dataLim.bounds))
             and a.dataLim.width >= 0 and a.dataLim.height >= 0]
    if not boxes:
        return
    bb = Bbox.union(boxes)
    mx = margin * (bb.width or 1.0)
    my = margin * (bb.height or 1.0)
    for a in axes:
        a.set_xlim(bb.x0 - mx, bb.x1 + mx)
        a.set_ylim(bb.y0 - my, bb.y1 + my)


def _draw_stack_row(ax, data, shifts, *, segments, removed, show_removed, observed,
                    observed_color, observed_size, cvalues, cmap, norm, point_size,
                    show_fit, label_basins, LineCollection):
    """Draw the basins of one stacked row (``{basin: shift}``) on ``ax``;
    returns the colour mappable."""
    mappable = None
    all_lines, all_cvals = [], []
    for b, shift in shifts.items():
        rows = data[data["basin_index"] == b]
        if rows.empty:
            continue
        chi0, chi_max = _basin_chi_bounds(rows, segments)
        width = chi_max - chi0
        cursor = shift + chi0
        rows, ghost = _split_removed(rows, removed, show_removed)
        if segments and len(ghost):
            ax.add_collection(LineCollection(_segment_lines(ghost, shift), colors=GHOST_COLOUR,
                                             linewidths=0.6, zorder=0))
        elif len(ghost):
            ax.scatter(ghost["chi"].to_numpy() + shift, ghost["observed_elevation"],
                       c=GHOST_COLOUR, s=point_size * 0.6, edgecolors="none", zorder=0)
        if observed is not None:
            _scatter_observed(ax, observed, b, shift, color=observed_color, size=observed_size)
        if rows.empty:
            continue
        if segments:
            all_lines.extend(_segment_lines(rows, shift))
            all_cvals.append(cvalues(rows))
            zmin = float(rows["downstream_elevation"].min())
        else:
            mappable = ax.scatter(rows["chi"].to_numpy() + shift, rows["observed_elevation"],
                                  c=cvalues(rows), cmap=cmap, norm=norm, s=point_size,
                                  edgecolors="none", zorder=2)
            if show_fit and "fitted_elevation" in rows.columns:
                for _, branch in rows.groupby("branch_id"):
                    branch = branch.sort_values("chi")
                    ax.plot(branch["chi"].to_numpy() + shift, branch["fitted_elevation"],
                            color="0.25", linewidth=0.6, zorder=1)
            zmin = float(rows["observed_elevation"].min())
        if label_basins:
            ax.text(cursor + 0.05 * (width or 1.0), zmin, str(int(b)),
                    fontsize=8, va="bottom", ha="left")
    if segments and all_lines:
        lc = LineCollection(all_lines, cmap=cmap, norm=norm, linewidths=1.2)
        lc.set_array(np.concatenate(all_cvals))
        ax.add_collection(lc)
        mappable = lc
    return mappable


def plot_ksn_vs_erate(
    data: pd.DataFrame,
    *,
    erate_col: str = "inferred_erate_mm_yr",
    ksn_col: str = "ksn_median",
    erate_err_low_col: str = "erate_err_low_mm_yr",
    erate_err_high_col: str = "erate_err_high_mm_yr",
    ksn_q25_col: str = "ksn_q25",
    ksn_q75_col: str = "ksn_q75",
    basin_label_col: Optional[str] = "basin_index",
    log_axes: bool = True,
    theta: Optional[float] = None,
    ax=None,
):
    """Plot basin-median k_sn versus inferred erosion rate with uncertainties.

    Y-errors use ``ksn_median`` against ``ksn_q25``/``ksn_q75`` (IQR). X-errors
    use the provided erosion-rate uncertainty columns.
    """
    plt, _, _ = _import_mpl()

    cols = [erate_col, ksn_col, ksn_q25_col, ksn_q75_col]
    for col in cols:
        if col not in data.columns:
            raise KeyError(f"Required column missing for k_sn vs erate plot: {col}")

    work = data.copy()
    for col in [erate_col, ksn_col, ksn_q25_col, ksn_q75_col, erate_err_low_col, erate_err_high_col]:
        if col in work.columns:
            work[col] = pd.to_numeric(work[col], errors="coerce")

    work = work[np.isfinite(work[erate_col]) & np.isfinite(work[ksn_col])]
    if work.empty:
        raise ValueError("No finite k_sn/erosion-rate rows available for plotting.")

    yerr_low = (work[ksn_col] - work[ksn_q25_col]).clip(lower=0).to_numpy(dtype=float)
    yerr_high = (work[ksn_q75_col] - work[ksn_col]).clip(lower=0).to_numpy(dtype=float)

    if erate_err_low_col in work.columns:
        xerr_low = work[erate_err_low_col].fillna(0).clip(lower=0).to_numpy(dtype=float)
    else:
        xerr_low = np.zeros(len(work), dtype=float)
    if erate_err_high_col in work.columns:
        xerr_high = work[erate_err_high_col].fillna(0).clip(lower=0).to_numpy(dtype=float)
    else:
        xerr_high = np.zeros(len(work), dtype=float)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.5, 4.0))
    else:
        fig = ax.figure

    ax.errorbar(
        work[erate_col].to_numpy(dtype=float),
        work[ksn_col].to_numpy(dtype=float),
        xerr=np.vstack([xerr_low, xerr_high]),
        yerr=np.vstack([yerr_low, yerr_high]),
        fmt="o",
        ms=4,
        color="0.2",
        ecolor="0.45",
        elinewidth=0.8,
        capsize=2,
    )

    if basin_label_col is not None and basin_label_col in work.columns:
        for _, row in work.iterrows():
            ax.text(row[erate_col], row[ksn_col], str(row[basin_label_col]), fontsize=7, ha="left", va="bottom")

    if log_axes:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel("Inferred erosion rate E (mm/yr)")
    ax.set_ylabel("Median ksn")
    if theta is None:
        ax.set_title("Basin median ksn vs inferred erosion rate")
    else:
        ax.set_title(rf"Basin median ksn vs inferred erosion rate (theta={theta:g})")
    return fig, ax
