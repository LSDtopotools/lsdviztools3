"""Valley cross-section plotting (matplotlib).

:func:`plot_valley_cross_sections` draws a grid of regularly-spaced transects
orthogonal to a valley centreline, from the per-sample table
(``*_valley_cross_section_samples.csv``) written by ``lsdtt-valley-centreline``
with ``write_valley_cross_sections``. Each panel plots elevation against the
signed offset from the centreline, shades the in-valley run, and marks the
detected left / right valley sides. Every transect of a channel shares the same
length (a fixed multiple of the widest valley on that channel), so the panels
are directly comparable.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

_MPL_MSG = (
    "Cross-section plotting needs matplotlib. Install the optional 'charts' "
    "extra, e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


_RASTERIO_MSG = (
    "The cross-section map needs rasterio to read the hillshade. Install it "
    "e.g. `pip install rasterio`."
)


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


def plot_valley_cross_sections(
    df: pd.DataFrame,
    *,
    source_id: int = 0,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    offset_col: str = "offset_m",
    elevation_col: str = "elevation",
    in_valley_col: str = "in_valley",
    along_col: str = "along_distance_m",
    left_col: str = "left_offset_m",
    right_col: str = "right_offset_m",
    max_panels: int = 30,
    n_cols: int = 5,
    elevation_color: str = "0.35",
    valley_color: str = "tab:blue",
    left_color: str = "tab:green",
    right_color: str = "tab:red",
):
    """Grid of valley cross sections orthogonal to the centreline.

    Args:
        df: Per-sample table (``*_valley_cross_section_samples.csv``) with the
            source, station, offset, elevation, in-valley, along-distance and
            left/right-offset columns. Non-finite elevations are dropped from
            the elevation line but kept for the in-valley shading.
        source_id: Channel source to plot.
        source_col, station_col, offset_col, elevation_col, in_valley_col,
        along_col, left_col, right_col: Column names in ``df``.
        max_panels: Cap on the number of transect panels; when a channel has
            more stations they are evenly subsampled (the vector/CSV outputs
            keep every station — this is a visual summary). 0 or negative =
            draw all stations.
        n_cols: Panels per row.
        elevation_color: Colour of the elevation profile line.
        valley_color: Fill colour for the in-valley run.
        left_color, right_color: Colours of the left / right valley-side marks.

    Returns:
        ``(fig, axes)`` — axes as a flat list, one per transect drawn.

    Raises:
        ValueError: if the table holds no stations for ``source_id``.
    """
    plt = _import_mpl()
    from matplotlib.lines import Line2D

    ch = df[df[source_col] == source_id]
    stations = sorted(ch[station_col].unique())
    if not stations:
        raise ValueError(
            f"No cross-section stations for {source_col}={source_id}.")
    if max_panels and max_panels > 0 and len(stations) > max_panels:
        idx = np.linspace(0, len(stations) - 1, max_panels).round().astype(int)
        stations = [stations[i] for i in sorted(set(idx))]

    n = len(stations)
    n_cols = max(1, min(n_cols, n))
    n_rows = (n + n_cols - 1) // n_cols
    fig, axarr = plt.subplots(
        n_rows, n_cols, figsize=(3.1 * n_cols, 2.3 * n_rows), squeeze=False)

    ev = ch[elevation_col].replace([np.inf, -np.inf], np.nan).dropna()
    ymin, ymax = (float(ev.min()), float(ev.max())) if len(ev) else (0.0, 1.0)
    pad = 0.05 * (ymax - ymin + 1e-6)
    xmax = float(ch[offset_col].abs().max())

    axes = []
    for k, sseq in enumerate(stations):
        ax = axarr[k // n_cols][k % n_cols]
        axes.append(ax)
        s = ch[ch[station_col] == sseq].sort_values(offset_col)
        off = s[offset_col].to_numpy(dtype=float)
        elev = s[elevation_col].to_numpy(dtype=float)
        inval = s[in_valley_col].to_numpy(dtype=float)
        lo = float(s[left_col].iloc[0])
        ro = float(s[right_col].iloc[0])
        along = float(s[along_col].iloc[0])

        ax.plot(off, elev, "-", color=elevation_color, lw=1.0)
        ax.fill_between(off, ymin - pad, elev, where=(inval > 0),
                        color=valley_color, alpha=0.25, step="mid")
        if np.isfinite(lo):
            ax.axvline(lo, color=left_color, lw=1.1, ls="--")
        if np.isfinite(ro):
            ax.axvline(ro, color=right_color, lw=1.1, ls="--")
        ax.axvline(0.0, color="k", lw=0.7, ls=":")   # centreline
        ax.set_ylim(ymin - pad, ymax + pad)
        ax.set_xlim(-xmax, xmax)
        ax.set_title(f"{along / 1000:.2f} km  (w={ro - lo:.0f} m)", fontsize=8)
        ax.tick_params(labelsize=7)
        if k % n_cols == 0:
            ax.set_ylabel("elev (m)", fontsize=8)
        if k // n_cols == n_rows - 1:
            ax.set_xlabel("offset (m)", fontsize=8)

    for k in range(n, n_rows * n_cols):
        axarr[k // n_cols][k % n_cols].set_axis_off()

    handles = [
        Line2D([0], [0], color=elevation_color, lw=1.0, label="elevation"),
        Line2D([0], [0], color=left_color, lw=1.1, ls="--", label="left side"),
        Line2D([0], [0], color=right_color, lw=1.1, ls="--", label="right side"),
        Line2D([0], [0], color="k", lw=0.7, ls=":", label="centreline"),
    ]
    fig.legend(handles=handles, loc="upper right", fontsize=8, frameon=False,
               ncol=4, bbox_to_anchor=(0.995, 0.998))
    fig.suptitle(
        f"channel {source_id}: valley cross sections (orthogonal to centreline)",
        fontsize=11, y=0.999)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig, axes


def plot_valley_cross_section_map(
    samples_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
    hillshade_path: str,
    *,
    source_id: int = 0,
    transform=None,
    sides_df: Optional[pd.DataFrame] = None,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    offset_col: str = "offset_m",
    x_col: str = "x",
    y_col: str = "y",
    node_x_col: str = "x",
    node_y_col: str = "y",
    label_every: int = 1,
    centreline_color: str = "crimson",
    transect_color: str = "yellow",
    side_color: str = "deepskyblue",
    crop_margin_frac: float = 0.08,
    ax=None,
):
    """Map of the numbered cross-section transects over a hillshade.

    A companion to :func:`plot_valley_cross_sections`: instead of the transect
    profiles, this shows WHERE each transect sits. It draws the hillshade, the
    valley centreline, and every cross-section transect as a line labelled with
    its station number, so a reader can match a numbered panel in the profile
    grid to its location on the ground.

    Transect endpoints are taken from the per-sample table: for each station
    the samples at the minimum and maximum signed offset are the two ends of
    the line. World coordinates are mapped to raster pixel space with the
    hillshade's affine transform (the plot is drawn in pixel/image space, as
    ``imshow`` requires).

    Args:
        samples_df: Per-sample cross-section table
            (``*_valley_cross_section_samples.csv``) with the source, station,
            offset and world ``x``/``y`` columns.
        nodes_df: Valley-centreline nodes table
            (``*_valley_centreline_nodes.csv``) with source and world
            ``x``/``y`` columns, ordered upstream→downstream.
        hillshade_path: Path to the hillshade GeoTIFF (read with rasterio).
        source_id: Channel source to map.
        transform: Optional affine transform for the world→pixel mapping. When
            ``None`` the hillshade's own transform is used (the normal case);
            pass one only if ``samples_df`` is georeferenced to a different
            grid than the hillshade.
        sides_df: Optional detected valley-side points (loaded from
            ``*_valley_cross_section_sides.fgb``) with ``x``/``y`` columns, to
            dot the left/right valley margins. ``None`` skips them.
        source_col, station_col, offset_col, x_col, y_col: Column names in
            ``samples_df``.
        node_x_col, node_y_col: World-coordinate columns in ``nodes_df``.
        label_every: Label every Nth station (all transect lines are always
            drawn; this only thins the number labels). 1 labels every station.
        centreline_color, transect_color, side_color: Line/point colours.
        crop_margin_frac: Fraction of the data bounding-box size added as a
            margin when cropping the map to the channel. 0 or negative shows
            the full hillshade.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax)``.

    Raises:
        ValueError: if the table holds no stations for ``source_id``.
    """
    plt = _import_mpl()
    rasterio = _import_rasterio()
    from matplotlib.lines import Line2D
    import matplotlib.patheffects as pe

    _label_stroke = [pe.withStroke(linewidth=1.4, foreground="black")]

    ch = samples_df[samples_df[source_col] == source_id]
    stations = sorted(ch[station_col].unique())
    if not stations:
        raise ValueError(
            f"No cross-section stations for {source_col}={source_id}.")

    with rasterio.open(hillshade_path) as src:
        hs = src.read(1).astype(float)
        if src.nodata is not None:
            hs[hs == src.nodata] = np.nan
        tr = transform if transform is not None else src.transform

    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 9))
    else:
        fig = ax.figure
    ax.imshow(hs, cmap="gray", interpolation="nearest")

    all_rows, all_cols = [], []

    # Centreline.
    g = nodes_df[nodes_df[source_col] == source_id]
    if len(g):
        rows, cols = rasterio.transform.rowcol(
            tr, g[node_x_col].to_numpy(), g[node_y_col].to_numpy())
        rows = np.asarray(rows, dtype=float)
        cols = np.asarray(cols, dtype=float)
        ax.plot(cols, rows, "-", lw=1.6, color=centreline_color,
                zorder=3, label="centreline")
        ax.plot(cols[0], rows[0], "o", ms=6, color="black", zorder=4)
        ax.plot(cols[-1], rows[-1], "s", ms=6, color="purple", zorder=4)
        all_rows.append(rows)
        all_cols.append(cols)

    # Transects: each station's two endpoints are its min/max-offset samples.
    for k, sseq in enumerate(stations):
        s = ch[ch[station_col] == sseq].sort_values(offset_col)
        ends = s.iloc[[0, -1]]
        rr, cc = rasterio.transform.rowcol(
            tr, ends[x_col].to_numpy(), ends[y_col].to_numpy())
        rr = np.asarray(rr, dtype=float)
        cc = np.asarray(cc, dtype=float)
        ax.plot(cc, rr, "-", lw=0.9, color=transect_color, alpha=0.9, zorder=2)
        all_rows.append(rr)
        all_cols.append(cc)
        if label_every > 0 and (k % label_every == 0):
            # Label at the first endpoint, nudged along the transect.
            dx, dy = cc[0] - cc[-1], rr[0] - rr[-1]
            norm = float(np.hypot(dx, dy)) or 1.0
            lx = cc[0] + 8.0 * dx / norm
            ly = rr[0] + 8.0 * dy / norm
            ax.text(lx, ly, str(int(sseq)), fontsize=6.5, color=transect_color,
                    ha="center", va="center", zorder=5,
                    path_effects=_label_stroke)

    # Optional valley-side points.
    if sides_df is not None and len(sides_df):
        sd = sides_df
        if source_col in sd.columns:
            sd = sd[sd[source_col] == source_id]
        if len(sd):
            srow, scol = rasterio.transform.rowcol(
                tr, sd[x_col].to_numpy(), sd[y_col].to_numpy())
            ax.plot(np.asarray(scol, dtype=float), np.asarray(srow, dtype=float),
                    ".", ms=3, color=side_color, zorder=3, label="valley side")

    # Crop to the channel with a margin.
    if crop_margin_frac and crop_margin_frac > 0 and all_rows:
        rr = np.concatenate(all_rows)
        cc = np.concatenate(all_cols)
        r0, r1 = np.nanmin(rr), np.nanmax(rr)
        c0, c1 = np.nanmin(cc), np.nanmax(cc)
        mr = crop_margin_frac * max(r1 - r0, 1.0)
        mc = crop_margin_frac * max(c1 - c0, 1.0)
        ax.set_xlim(max(c0 - mc, 0), min(c1 + mc, hs.shape[1]))
        # y inverted (row 0 at top) — imshow default.
        ax.set_ylim(min(r1 + mr, hs.shape[0]), max(r0 - mr, 0))

    handles = [
        Line2D([0], [0], color=centreline_color, lw=1.6, label="centreline"),
        Line2D([0], [0], color=transect_color, lw=1.2,
               label="cross section (numbered)"),
    ]
    if sides_df is not None and len(sides_df):
        handles.append(Line2D([0], [0], color=side_color, lw=0, marker=".",
                              ms=6, label="valley side"))
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True,
              framealpha=0.7)
    ax.set_title(
        f"channel {source_id}: cross-section locations (numbered by station)",
        fontsize=11)
    ax.set_axis_off()
    return fig, ax


# --- Curved (boundary-fitted) cross sections --------------------------------

def _label_height_data_units(ax, fontsize):
    """Height of a two-digit label at ``fontsize``, expressed in the axes' data
    units (metres, when the map is drawn in world coordinates). Needed so the
    adaptive-numbering gap can be stated as a multiple of the label height."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    t = ax.text(float(np.mean(ax.get_xlim())), float(np.mean(ax.get_ylim())),
                "88", fontsize=fontsize)
    h_px = t.get_window_extent(renderer=renderer).height
    t.remove()
    inv = ax.transData.inverted()
    (_, y0) = inv.transform((0.0, 0.0))
    (_, y1) = inv.transform((0.0, h_px))
    return abs(y1 - y0)


def plot_valley_curved_cross_section_map(
    samples_df: pd.DataFrame,
    hillshade_path: str,
    *,
    source_id: int = 0,
    nodes_df: Optional[pd.DataFrame] = None,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    offset_col: str = "offset_m",
    x_col: str = "x",
    y_col: str = "y",
    node_x_col: str = "x",
    node_y_col: str = "y",
    number_labels: bool = True,
    label_gap_factor: float = 2.0,
    label_fontsize: int = 8,
    section_cmap: str = "turbo",
    section_lw: float = 1.4,
    centreline_color: str = "tab:cyan",
    crop_margin_metres: float = 300.0,
    ax=None,
):
    """Map of the CURVED (boundary-fitted) cross sections over a hillshade.

    The companion to :func:`plot_valley_cross_section_map` for the curved
    sections written by ``lsdtt-valley-centreline`` with
    ``curved_cross_sections: true``. Each section is an isocontour of the
    along-valley coordinate, so it is a *polyline* (many samples), not a straight
    chord: this draws the full traced curve of every section from its per-sample
    world coordinates. The map is drawn in world coordinates (metres) via the
    hillshade's ``extent`` so the curves and the adaptive number labels share one
    metric coordinate system.

    Station numbers are placed ADAPTIVELY. Walking downstream, a section is
    labelled only when its label anchor (the consistent left-side endpoint,
    nudged outward) is at least ``label_gap_factor`` times the label height from
    the previously placed label. A sparse network therefore numbers every
    section; a dense one numbers an evenly-spaced subset, with the text never
    closer than ``label_gap_factor`` label-heights — directly honouring "how far
    apart the closest endpoints are".

    Args:
        samples_df: Per-sample cross-section table
            (``*_valley_cross_section_samples.csv``). For a curved section the
            samples trace the curve; ``offset_m`` is the signed arc-length along
            it (negative left, positive right of downstream).
        hillshade_path: Path to the hillshade GeoTIFF (read with rasterio).
        source_id: Channel source to map.
        nodes_df: Optional valley-centreline nodes
            (``*_valley_centreline_nodes.csv``) to draw the centreline. ``None``
            skips it.
        source_col, station_col, offset_col, x_col, y_col: Column names in
            ``samples_df``.
        node_x_col, node_y_col: World-coordinate columns in ``nodes_df``.
        number_labels: Draw adaptive station-number labels.
        label_gap_factor: Minimum spacing between placed labels, as a multiple
            of the label height (2.0 = at least two label-heights apart).
        label_fontsize: Font size of the station-number labels.
        section_cmap: Matplotlib colormap for the sections (coloured
            downstream). Pass a single colour name to draw them all one colour.
        section_lw: Section line width.
        centreline_color: Colour of the centreline (when ``nodes_df`` given).
        crop_margin_metres: Margin (metres) added around the sections' bounding
            box when cropping the map. 0 or negative shows the full hillshade.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.

    Returns:
        ``(fig, ax, n_labels)`` — the figure, axes, and the number of station
        labels actually placed.

    Raises:
        ValueError: if the table holds no stations for ``source_id``.
    """
    plt = _import_mpl()
    rasterio = _import_rasterio()
    from matplotlib.lines import Line2D
    import matplotlib.patheffects as pe

    _stroke = [pe.withStroke(linewidth=2.0, foreground="white")]

    ch = samples_df[samples_df[source_col] == source_id]
    stations = sorted(ch[station_col].unique())
    if not stations:
        raise ValueError(
            f"No cross-section stations for {source_col}={source_id}.")

    with rasterio.open(hillshade_path) as src:
        hs = src.read(1).astype(float)
        if src.nodata is not None:
            hs[hs == src.nodata] = np.nan
        tr = src.transform
        extent = (tr.c, tr.c + tr.a * src.width,
                  tr.f + tr.e * src.height, tr.f)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
    else:
        fig = ax.figure
    ax.imshow(hs, extent=extent, cmap="gray", origin="upper", zorder=0)
    ax.set_aspect("equal")

    # Optional centreline.
    if nodes_df is not None:
        g = nodes_df[nodes_df[source_col] == source_id]
        if len(g):
            ax.plot(g[node_x_col].to_numpy(), g[node_y_col].to_numpy(), "-",
                    color=centreline_color, lw=1.0, alpha=0.85, zorder=2)

    # Per-station geometry: the full ordered polyline, plus the left-side
    # endpoint (min signed offset) and centre (offset nearest 0) for labelling.
    try:
        cmap = plt.get_cmap(section_cmap)
        single_colour = None
    except (ValueError, TypeError):
        cmap, single_colour = None, section_cmap

    label_info = []
    xs_all, ys_all = [], []
    n = len(stations)
    for i, sseq in enumerate(stations):
        s = ch[ch[station_col] == sseq].sort_values(offset_col)
        px = s[x_col].to_numpy(); py = s[y_col].to_numpy()
        if len(px) < 2:
            continue
        colour = single_colour if single_colour is not None \
            else cmap(i / max(1, n - 1))
        ax.plot(px, py, "-", lw=section_lw, color=colour, zorder=3)
        xs_all.append(px); ys_all.append(py)
        # Left endpoint (most negative offset) as the consistent label side;
        # nudge outward from the centre (offset closest to 0).
        c_idx = int(np.argmin(np.abs(s[offset_col].to_numpy())))
        anchor = np.array([px[0], py[0]])
        centre = np.array([px[c_idx], py[c_idx]])
        out = anchor - centre
        nrm = float(np.hypot(*out))
        out = out / nrm if nrm else np.array([-1.0, 0.0])
        label_info.append(dict(seq=int(sseq), anchor=anchor, out=out))

    # Adaptive station numbering.
    n_labels = 0
    if number_labels and label_info:
        H = _label_height_data_units(ax, label_fontsize)
        min_gap = label_gap_factor * H
        last = None
        for it in label_info:
            a = it["anchor"] + it["out"] * (0.75 * H)
            if last is None or float(np.hypot(*(a - last))) >= min_gap:
                ha = "left" if it["out"][0] >= 0 else "right"
                va = "bottom" if it["out"][1] >= 0 else "top"
                ax.text(a[0], a[1], str(it["seq"]), fontsize=label_fontsize,
                        ha=ha, va=va, color="black", zorder=6,
                        path_effects=_stroke)
                last = a
                n_labels += 1

    # Crop to the sections' bounding box with a margin.
    if crop_margin_metres and crop_margin_metres > 0 and xs_all:
        ax_x = np.concatenate(xs_all); ax_y = np.concatenate(ys_all)
        x0, x1 = np.nanmin(ax_x), np.nanmax(ax_x)
        y0, y1 = np.nanmin(ax_y), np.nanmax(ax_y)
        m = crop_margin_metres
        ax.set_xlim(x0 - m, x1 + m)
        ax.set_ylim(y0 - m, y1 + m)

    handles = [
        Line2D([0], [0], color="0.4", lw=section_lw,
               label="curved cross section"),
    ]
    if nodes_df is not None:
        handles.append(Line2D([0], [0], color=centreline_color, lw=1.2,
                              label="centreline"))
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True,
              framealpha=0.7)
    ax.set_title(
        f"channel {source_id}: curved (boundary-fitted) cross sections — "
        f"{len(stations)} sections, {n_labels} labels", fontsize=11)
    ax.tick_params(labelsize=6)
    return fig, ax, n_labels


# --- Feature-annotated cross sections ---------------------------------------

# Marker style per feature type (colour, matplotlib marker).
_FEATURE_STYLE = {
    "channel":        ("tab:blue",   "v"),
    "wall_left":      ("tab:green",  "|"),
    "wall_right":     ("tab:red",    "|"),
    "curvature_peak": ("tab:orange", "o"),
    "hillshade_peak": ("tab:purple", "x"),
}


def plot_valley_cross_section_features(
    samples_df: pd.DataFrame,
    features_df: pd.DataFrame,
    *,
    source_id: int = 0,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    offset_col: str = "offset_m",
    elevation_col: str = "elevation",
    in_valley_col: str = "in_valley",
    along_col: str = "along_distance_m",
    curvature_col: str = "curvature",
    hillshade_col: str = "hillshade",
    feature_type_col: str = "feature_type",
    feature_offset_col: str = "offset_m",
    max_panels: int = 20,
    n_cols: int = 4,
):
    """Individual cross sections annotated with detected features.

    For each station, two stacked panels sharing the offset axis:
      * top    — elevation profile with in-valley shading, and every detected
                 feature marked at its offset: D8 channel crossings, left/right
                 valley walls, maximum-curvature peaks, hillshade peaks;
      * bottom — the maximum-curvature profile (left axis) and multidirectional
                 hillshade profile (right axis) the (c)/(d) peaks come from,
                 with vertical guides at the peak offsets.

    Args:
        samples_df: per-sample table (``*_valley_cross_section_samples.csv``),
            which must include the ``curvature``, ``hillshade`` and
            ``is_channel`` columns written when features are enabled.
        features_df: detected-feature table
            (``*_valley_cross_section_features.csv``) with a ``feature_type``
            and ``offset_m`` per feature.
        source_id: channel source to plot.
        max_panels: cap on stations shown (evenly subsampled if exceeded;
            0/negative draws all). n_cols: station columns.

    Returns:
        ``(fig, axes)`` — ``axes`` a list of ``(elev_ax, feat_ax)`` pairs.

    Raises:
        ValueError: if no stations exist for ``source_id``.
    """
    plt = _import_mpl()
    from matplotlib.lines import Line2D

    ch = samples_df[samples_df[source_col] == source_id]
    stations = sorted(ch[station_col].unique())
    if not stations:
        raise ValueError(
            f"No cross-section stations for {source_col}={source_id}.")
    if max_panels and max_panels > 0 and len(stations) > max_panels:
        idx = np.linspace(0, len(stations) - 1, max_panels).round().astype(int)
        stations = [stations[i] for i in sorted(set(idx))]

    feats = features_df[features_df[source_col] == source_id] \
        if features_df is not None else None

    n = len(stations)
    n_cols = max(1, min(n_cols, n))
    n_block_rows = (n + n_cols - 1) // n_cols
    # 2 grid rows per station row: elevation on top, curv/hillshade below.
    fig, axarr = plt.subplots(
        2 * n_block_rows, n_cols,
        figsize=(3.4 * n_cols, 3.2 * n_block_rows), squeeze=False,
        gridspec_kw={"height_ratios": [2.0, 1.0] * n_block_rows})

    ev = ch[elevation_col].replace([np.inf, -np.inf], np.nan).dropna()
    ymin, ymax = (float(ev.min()), float(ev.max())) if len(ev) else (0.0, 1.0)
    pad = 0.05 * (ymax - ymin + 1e-6)
    xmax = float(ch[offset_col].abs().max())

    axes = []
    for k, sseq in enumerate(stations):
        br, col = k // n_cols, k % n_cols
        eax = axarr[2 * br][col]
        fax = axarr[2 * br + 1][col]
        axes.append((eax, fax))

        s = ch[ch[station_col] == sseq].sort_values(offset_col)
        off = s[offset_col].to_numpy(dtype=float)
        elev = s[elevation_col].to_numpy(dtype=float)
        inval = s[in_valley_col].to_numpy(dtype=float)
        along = float(s[along_col].iloc[0])
        fs = feats[feats[station_col] == sseq] if feats is not None else None

        def elev_at(o):
            return float(np.interp(o, off, elev))

        # --- elevation panel ---
        eax.plot(off, elev, "-", color="0.35", lw=1.1, zorder=2)
        eax.fill_between(off, ymin - pad, elev, where=(inval > 0),
                         color="tab:blue", alpha=0.18, step="mid", zorder=1)
        eax.axvline(0.0, color="k", lw=0.7, ls=":", zorder=1)
        if fs is not None:
            for ftype, (colr, mk) in _FEATURE_STYLE.items():
                fo = fs[fs[feature_type_col] == ftype][feature_offset_col]
                for o in fo.to_numpy(dtype=float):
                    if ftype in ("wall_left", "wall_right"):
                        eax.axvline(o, color=colr, lw=1.1, ls="--", zorder=3)
                    else:
                        eax.plot(o, elev_at(o), mk, color=colr, ms=6,
                                 mew=1.4, zorder=4)
                    if ftype in ("curvature_peak", "hillshade_peak"):
                        fax.axvline(o, color=colr, lw=0.8, ls=":", alpha=0.7)
        eax.set_ylim(ymin - pad, ymax + pad)
        eax.set_xlim(-xmax, xmax)
        eax.set_title(f"station {int(sseq)} — {along / 1000:.2f} km",
                      fontsize=8)
        eax.tick_params(labelsize=6, labelbottom=False)
        if col == 0:
            eax.set_ylabel("elev (m)", fontsize=7)

        # --- curvature / hillshade panel ---
        curv = s[curvature_col].to_numpy(dtype=float) \
            if curvature_col in s else np.full_like(off, np.nan)
        hsv = s[hillshade_col].to_numpy(dtype=float) \
            if hillshade_col in s else np.full_like(off, np.nan)
        fax.plot(off, curv, "-", color="tab:orange", lw=0.9, zorder=2)
        fax.axhline(0.0, color="0.7", lw=0.5)
        fax.set_xlim(-xmax, xmax)
        fax.tick_params(labelsize=6, colors="tab:orange")
        if col == 0:
            fax.set_ylabel("max curv", fontsize=7, color="tab:orange")
        hax = fax.twinx()
        hax.plot(off, hsv, "-", color="tab:purple", lw=0.9, alpha=0.8, zorder=2)
        hax.tick_params(labelsize=6, colors="tab:purple")
        if col == n_cols - 1:
            hax.set_ylabel("hillshade", fontsize=7, color="tab:purple")
        if br == n_block_rows - 1:
            fax.set_xlabel("offset (m)", fontsize=7)

    for k in range(n, n_block_rows * n_cols):
        br, col = k // n_cols, k % n_cols
        axarr[2 * br][col].set_axis_off()
        axarr[2 * br + 1][col].set_axis_off()

    handles = [
        Line2D([0], [0], color="0.35", lw=1.1, label="elevation"),
        Line2D([0], [0], color="tab:green", lw=1.1, ls="--", label="left wall"),
        Line2D([0], [0], color="tab:red", lw=1.1, ls="--", label="right wall"),
        Line2D([0], [0], color="tab:blue", lw=0, marker="v", ms=7,
               label="D8 channel"),
        Line2D([0], [0], color="tab:orange", lw=0, marker="o", ms=6,
               label="max-curvature peak"),
        Line2D([0], [0], color="tab:purple", lw=0, marker="x", ms=6,
               label="hillshade peak"),
    ]
    fig.legend(handles=handles, loc="upper center", fontsize=8, frameon=False,
               ncol=6, bbox_to_anchor=(0.5, 0.998))
    fig.suptitle(
        f"channel {source_id}: cross sections with detected features",
        fontsize=11, y=1.035)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    return fig, axes


_MAP_FEATURE_TYPES = ("channel", "wall_left", "wall_right")


def plot_valley_cross_section_features_map(
    features_df: pd.DataFrame,
    samples_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
    hillshade_path: str,
    *,
    source_id: int = 0,
    boundaries_df: Optional[pd.DataFrame] = None,
    valley_mask_path: Optional[str] = None,
    flat_edges_df: Optional[pd.DataFrame] = None,
    transform=None,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    offset_col: str = "offset_m",
    x_col: str = "x",
    y_col: str = "y",
    feature_type_col: str = "feature_type",
    boundary_type_col: str = "boundary_type",
    step_col: str = "step",
    node_x_col: str = "x",
    node_y_col: str = "y",
    label_every: int = 1,
    hydraulic_label_every: int = 0,
    valley_mask_alpha: float = 0.28,
    centreline_color: str = "crimson",
    transect_color: str = "0.75",
    crop_margin_frac: float = 0.08,
    ax=None,
):
    """Map of cross-section transects with valley-edge and terrace/floodplain
    edge points, over a hillshade.

    ``label_every`` labels every Nth transect with its station number (all
    transect lines are always drawn); 0 disables labels.

    Draws the hillshade, optional valley-mask overlay, the centreline, every
    transect line (from the per-sample endpoints), channel crossings and wall
    boundaries, and hydraulic step-midpoint boundaries. If
    ``hydraulic_label_every`` > 0, hydraulic edges are labelled (station:step)
    every N points.

    Args:
        features_df: detected-feature table with ``feature_type`` and world
            ``x``/``y`` columns (channel and legacy wall points are used).
        boundaries_df: optional boundary table from
            ``*_valley_cross_section_boundaries.fgb`` with
            ``boundary_type``/``step`` and world ``x``/``y`` columns. When
            supplied, walls and hydraulic edges are read from here.
        valley_mask_path: optional path to the original
            ``*_valley_mask.tif``. For source ``k`` this overlays pixels equal
            to ``k+1`` (or non-zero in per-channel masks).
        flat_edges_df: output of :func:`compute_flat_surface_edges` (world
            ``x``/``y`` of the step-midpoint edges); optional.
        samples_df, nodes_df, hillshade_path: as in
            :func:`plot_valley_cross_section_map` (transect endpoints and
            centreline over the hillshade).
        source_id: channel source to map.

    Returns:
        ``(fig, ax)``.
    """
    plt = _import_mpl()
    rasterio = _import_rasterio()
    from matplotlib.lines import Line2D
    import matplotlib.patheffects as pe

    _label_stroke = [pe.withStroke(linewidth=1.4, foreground="black")]

    ch = samples_df[samples_df[source_col] == source_id]
    stations = sorted(ch[station_col].unique())
    if not stations:
        raise ValueError(
            f"No cross-section stations for {source_col}={source_id}.")

    with rasterio.open(hillshade_path) as src:
        hs = src.read(1).astype(float)
        if src.nodata is not None:
            hs[hs == src.nodata] = np.nan
        tr = transform if transform is not None else src.transform

    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 9))
    else:
        fig = ax.figure
    ax.imshow(hs, cmap="gray", interpolation="nearest")

    # Optional valley-mask overlay (original extracted valley raster).
    if valley_mask_path:
        with rasterio.open(valley_mask_path) as src:
            vm = src.read(1).astype(float)
            if src.nodata is not None:
                vm[vm == src.nodata] = np.nan
        vmax = np.nanmax(vm) if np.isfinite(vm).any() else np.nan
        if np.isfinite(vmax) and vmax <= 1.5:
            in_src = vm > 0.5
        else:
            in_src = np.isclose(vm, float(source_id + 1))
        rgba = np.zeros((hs.shape[0], hs.shape[1], 4), dtype=float)
        rgba[..., 1] = 0.70
        rgba[..., 2] = 0.62
        rgba[..., 3] = np.where(in_src, max(0.0, min(1.0, valley_mask_alpha)),
                                0.0)
        ax.imshow(rgba, interpolation="nearest", zorder=1)

    all_rows, all_cols = [], []

    g = nodes_df[nodes_df[source_col] == source_id]
    if len(g):
        rows, cols = rasterio.transform.rowcol(
            tr, g[node_x_col].to_numpy(), g[node_y_col].to_numpy())
        rows = np.asarray(rows, dtype=float); cols = np.asarray(cols, dtype=float)
        ax.plot(cols, rows, "-", lw=1.4, color=centreline_color, zorder=3)
        all_rows.append(rows); all_cols.append(cols)

    for k, sseq in enumerate(stations):
        s = ch[ch[station_col] == sseq].sort_values(offset_col)
        ends = s.iloc[[0, -1]]
        rr, cc = rasterio.transform.rowcol(
            tr, ends[x_col].to_numpy(), ends[y_col].to_numpy())
        rr = np.asarray(rr, dtype=float); cc = np.asarray(cc, dtype=float)
        ax.plot(cc, rr, "-", lw=0.7, color=transect_color, alpha=0.9, zorder=2)
        all_rows.append(rr); all_cols.append(cc)
        if label_every > 0 and (k % label_every == 0):
            dx, dy = cc[0] - cc[-1], rr[0] - rr[-1]
            norm = float(np.hypot(dx, dy)) or 1.0
            lx = cc[0] + 10.0 * dx / norm
            ly = rr[0] + 10.0 * dy / norm
            ax.text(lx, ly, str(int(sseq)), fontsize=7, color="yellow",
                    ha="center", va="center", zorder=6,
                    path_effects=_label_stroke)

    fe = features_df[features_df[source_col] == source_id] \
        if features_df is not None else None
    if fe is not None and len(fe):
        for ftype in _MAP_FEATURE_TYPES:
            colr, mk = _FEATURE_STYLE[ftype]
            sub = fe[fe[feature_type_col] == ftype]
            if not len(sub):
                continue
            rr, cc = rasterio.transform.rowcol(
                tr, sub[x_col].to_numpy(), sub[y_col].to_numpy())
            rr = np.asarray(rr, dtype=float); cc = np.asarray(cc, dtype=float)
            ax.plot(cc, rr, linestyle="none", marker=mk, ms=5, mew=1.2,
                    color=colr, zorder=4, label=ftype.replace("_", " "))
            all_rows.append(rr); all_cols.append(cc)

    # Prefer explicit boundary FGB for walls/hydraulic edges when available.
    bd = boundaries_df[boundaries_df[source_col] == source_id] \
        if boundaries_df is not None else None
    if bd is not None and len(bd):
        style = {
            "wall_left": ("tab:green", "|", 7, 1.4),
            "wall_right": ("tab:red", "|", 7, 1.4),
            "hydraulic_edge": ("magenta", "D", 5, 0.8),
        }
        for btype, (colr, mk, ms, mew) in style.items():
            sub = bd[bd[boundary_type_col] == btype]
            if not len(sub):
                continue
            rr, cc = rasterio.transform.rowcol(
                tr, sub[x_col].to_numpy(), sub[y_col].to_numpy())
            rr = np.asarray(rr, dtype=float)
            cc = np.asarray(cc, dtype=float)
            if btype == "hydraulic_edge":
                ax.plot(cc, rr, linestyle="none", marker=mk, ms=ms, mew=mew,
                        mfc=colr, mec="black", zorder=6,
                        label="hydraulic edge")
                if hydraulic_label_every and hydraulic_label_every > 0:
                    sub = sub.reset_index(drop=True)
                    for i in range(0, len(sub), int(hydraulic_label_every)):
                        st = int(sub.loc[i, station_col]) \
                            if station_col in sub.columns else -1
                        if step_col in sub.columns and pd.notna(sub.loc[i, step_col]):
                            sp = int(sub.loc[i, step_col])
                            txt = f"{st}:{sp}"
                        else:
                            txt = str(st)
                        ax.text(cc[i] + 2.0, rr[i] + 2.0, txt, fontsize=5.5,
                                color="yellow", ha="left", va="bottom",
                                zorder=7, path_effects=_label_stroke)
            else:
                ax.plot(cc, rr, linestyle="none", marker=mk, ms=ms, mew=mew,
                        color=colr, zorder=5,
                        label=btype.replace("_", " "))
            all_rows.append(rr); all_cols.append(cc)

    if flat_edges_df is not None and len(flat_edges_df):
        ed = flat_edges_df[flat_edges_df[source_col] == source_id]
        if len(ed):
            rr, cc = rasterio.transform.rowcol(
                tr, ed[x_col].to_numpy(), ed[y_col].to_numpy())
            rr = np.asarray(rr, dtype=float); cc = np.asarray(cc, dtype=float)
            ax.plot(cc, rr, linestyle="none", marker="D", ms=5, mew=0.7,
                    mfc="magenta", mec="black", zorder=5,
                    label="terrace / floodplain edge")
            all_rows.append(rr); all_cols.append(cc)

    if crop_margin_frac and crop_margin_frac > 0 and all_rows:
        rr = np.concatenate(all_rows); cc = np.concatenate(all_cols)
        r0, r1 = np.nanmin(rr), np.nanmax(rr)
        c0, c1 = np.nanmin(cc), np.nanmax(cc)
        mr = crop_margin_frac * max(r1 - r0, 1.0)
        mc = crop_margin_frac * max(c1 - c0, 1.0)
        ax.set_xlim(max(c0 - mc, 0), min(c1 + mc, hs.shape[1]))
        ax.set_ylim(min(r1 + mr, hs.shape[0]), max(r0 - mr, 0))

    handles = [Line2D([0], [0], color=centreline_color, lw=1.4,
                      label="centreline"),
               Line2D([0], [0], color=transect_color, lw=1.0,
                      label="cross section")]
    if valley_mask_path:
        handles.append(Line2D([0], [0], color="#2cb3a0", lw=3,
                              label="original valley mask"))
    for ftype in _MAP_FEATURE_TYPES:
        colr, mk = _FEATURE_STYLE[ftype]
        handles.append(Line2D([0], [0], color=colr, lw=0, marker=mk, ms=6,
                              label=ftype.replace("_", " ")))
    if bd is not None and len(bd):
        handles.extend([
            Line2D([0], [0], color="tab:green", lw=0, marker="|", ms=8,
                   label="wall left"),
            Line2D([0], [0], color="tab:red", lw=0, marker="|", ms=8,
                   label="wall right"),
            Line2D([0], [0], color="magenta", lw=0, marker="D", ms=6,
                   mec="black", label="hydraulic edge"),
        ])
    if flat_edges_df is not None and len(flat_edges_df):
        handles.append(Line2D([0], [0], color="magenta", lw=0, marker="D",
                              ms=6, mec="black",
                              label="terrace / floodplain edge"))
    # Keep legend unique while preserving insertion order.
    uniq = {}
    for h in handles:
        uniq[h.get_label()] = h
    ax.legend(handles=list(uniq.values()), loc="upper right", fontsize=7, frameon=True,
              framealpha=0.7)
    ax.set_title(f"channel {source_id}: cross-section features", fontsize=11)
    ax.set_axis_off()
    return fig, ax


# --- Stage-rise hydraulic-geometry flat-surface detection -------------------

def cross_section_hydraulics(offset, elevation, stage_step=0.1):
    """Flood a cross-section profile and return stage-hydraulic curves.

    From the profile's lowest point, raise a water-surface stage ``h`` in
    ``stage_step`` (m) increments. At each stage, over the wetted run CONNECTED
    to the thalweg (lowest point), compute wetted area ``A`` (integral of depth
    over width), total flow width ``W``, wetted perimeter ``P`` (bed length
    across the run) and hydraulic depth ``D_h = A / W``.

    Args:
        offset: 1-D array of signed offsets (m), evenly spaced, ascending.
        elevation: 1-D array of elevations (m); non-finite samples are dropped.
        stage_step: stage increment (m).

    Returns:
        dict of numpy arrays ``stage, area, width, perimeter, hyd_depth`` plus
        scalars ``zmin`` and ``ds`` (offset spacing). Empty arrays if the
        profile has < 3 valid samples.
    """
    off = np.asarray(offset, dtype=float)
    z = np.asarray(elevation, dtype=float)
    ok = np.isfinite(z) & np.isfinite(off)
    off, z = off[ok], z[ok]
    order = np.argsort(off)
    off, z = off[order], z[order]
    empty = dict(stage=np.array([]), area=np.array([]), width=np.array([]),
                 perimeter=np.array([]), hyd_depth=np.array([]),
                 zmin=np.nan, ds=np.nan)
    n = len(z)
    if n < 3:
        return empty
    ds = float(np.median(np.diff(off)))
    if not np.isfinite(ds) or ds <= 0:
        return empty
    zmin, zmax = float(z.min()), float(z.max())
    thalweg = int(np.argmin(z))
    stages = np.arange(zmin + stage_step, zmax, stage_step)
    A = np.zeros_like(stages); W = np.zeros_like(stages)
    P = np.zeros_like(stages); Dh = np.zeros_like(stages)
    for i, h in enumerate(stages):
        if z[thalweg] >= h:
            continue
        lo = thalweg
        while lo - 1 >= 0 and z[lo - 1] < h:
            lo -= 1
        hi = thalweg
        while hi + 1 < n and z[hi + 1] < h:
            hi += 1
        seg = z[lo:hi + 1]
        depth = h - seg
        A[i] = float(np.sum(depth)) * ds
        W[i] = len(seg) * ds
        P[i] = (float(np.sum(np.hypot(ds, np.diff(seg))))
                if len(seg) > 1 else ds)
        Dh[i] = A[i] / W[i] if W[i] > 0 else 0.0
    return dict(stage=stages, area=A, width=W, perimeter=P, hyd_depth=Dh,
                zmin=zmin, ds=ds)


def _prominent_peaks(v, min_prom_frac=0.1):
    """Indices of prominent local maxima of a 1-D array (prominence >=
    min_prom_frac of the value range)."""
    v = np.asarray(v, dtype=float)
    n = len(v)
    if n < 3:
        return []
    rng = float(np.nanmax(v) - np.nanmin(v))
    if rng <= 0:
        return []
    thr = min_prom_frac * rng
    out = []
    for i in range(1, n - 1):
        if not (v[i] > v[i - 1] and v[i] >= v[i + 1]):
            continue
        left = v[i]
        for j in range(i - 1, -1, -1):
            if v[j] > v[i]:
                break
            left = min(left, v[j])
        right = v[i]
        for j in range(i + 1, n):
            if v[j] > v[i]:
                break
            right = min(right, v[j])
        if v[i] - max(left, right) >= thr:
            out.append(i)
    return out


def detect_hydraulic_flats(hyd, width_prom_frac=0.12,
                           hyd_depth_flat_frac=1.0):
    """Detect flat-surface stages from stage-hydraulic curves.

    A flat surface is a stage where flow WIDTH increases rapidly (a prominent
    peak of ``dW/dh``) while hydraulic depth stays roughly constant (``|dD_h/dh|``
    at the peak below ``hyd_depth_flat_frac`` x the median ``|dD_h/dh|``).

    Args:
        hyd: output of :func:`cross_section_hydraulics`.
        width_prom_frac: min prominence of a ``dW/dh`` peak (fraction of range).
        hyd_depth_flat_frac: keep a peak only if ``|dD_h/dh|`` there is below
            this multiple of the median ``|dD_h/dh|``. Large (e.g. 1e9) disables
            the constant-depth filter.

    Returns:
        ``(stages, dW_dh)`` — the detected flat-surface stages (elevations, m)
        and the full ``dW/dh`` curve (aligned to ``hyd['stage']``).
    """
    h = hyd["stage"]
    if len(h) < 3:
        return np.array([]), np.array([])
    dW = np.gradient(hyd["width"], h)
    dDh = np.gradient(hyd["hyd_depth"], h)
    med = float(np.median(np.abs(dDh))) or 1e-9
    peaks = _prominent_peaks(dW, width_prom_frac)
    keep = [i for i in peaks
            if abs(dDh[i]) <= hyd_depth_flat_frac * med]
    return h[keep], dW


def combine_width_steps(stage, width, ds, min_step_height_m=2.0,
                        min_step_width_px=3):
    """Combine the steps of a (monotonic) flow-width curve into flat surfaces.

    Candidate boundaries come from PROMINENT peaks of ``dW/dh`` (rapid widening
    of the flooded section), then adjacent candidates are merged unless they
    collectively clear both thresholds since the previous accepted boundary:
    ``>= min_step_width_px`` pixels of added width AND ``>= min_step_height_m``
    of stage rise. This avoids periodic over-splitting on long smooth walls
    where cumulative width growth can otherwise trigger a boundary every few
    metres.

    Args:
        stage: ascending stage elevations (m) from
            :func:`cross_section_hydraulics`.
        width: flow width at each stage (m), monotonic non-decreasing.
        ds: offset spacing (m); ``min_step_width_px`` is measured in these.
        min_step_height_m: minimum vertical rise for a new step boundary.
        min_step_width_px: minimum added width (in pixels) for a new boundary.

    Returns:
        ``(boundaries, midpoints)`` — step-boundary elevations (m, including the
        lowest and highest stage) and the mid-elevation of each step between
        consecutive boundaries.
    """
    stage = np.asarray(stage, dtype=float)
    width = np.asarray(width, dtype=float)
    ok = np.isfinite(stage) & np.isfinite(width)
    stage = stage[ok]
    width = width[ok]
    if len(stage) < 2 or not np.isfinite(ds) or ds <= 0:
        return [], []

    order = np.argsort(stage)
    stage = stage[order]
    width = width[order]

    min_w = min_step_width_px * ds

    # Candidate boundaries = prominent widening events (dW/dh peaks), bounded
    # by the lowest and highest stage.
    if len(stage) >= 3:
        dW = np.gradient(width, stage)
        pk = _prominent_peaks(dW, min_prom_frac=0.10)
    else:
        pk = []
    cand_idx = [0] + [i for i in pk if 0 < i < len(stage) - 1] + [len(stage) - 1]
    cand_idx = sorted(set(cand_idx))

    keep_idx = [cand_idx[0]]
    ref_idx = cand_idx[0]
    for i in cand_idx[1:-1]:
        if ((width[i] - width[ref_idx]) >= min_w
                and (stage[i] - stage[ref_idx]) >= min_step_height_m):
            keep_idx.append(i)
            ref_idx = i
    if keep_idx[-1] != cand_idx[-1]:
        keep_idx.append(cand_idx[-1])

    boundaries = []
    for i in keep_idx:
        z = float(stage[i])
        if not boundaries or (z - boundaries[-1]) > 1e-9:
            boundaries.append(z)

    if len(boundaries) < 2:
        return [], []

    mids = [0.5 * (boundaries[k] + boundaries[k + 1])
            for k in range(len(boundaries) - 1)]
    return boundaries, mids


def _water_line_crossings(offset, elevation, level):
    """Offsets where a profile crosses a horizontal water surface ``level``.

    Returns every left/right bank intersection at that stage (a mid-valley ridge
    poking above the water yields 4 crossings instead of 2), by linear
    interpolation between adjacent samples that straddle ``level``.
    """
    off = np.asarray(offset, dtype=float)
    z = np.asarray(elevation, dtype=float)
    ok = np.isfinite(off) & np.isfinite(z)
    off, z = off[ok], z[ok]
    order = np.argsort(off)
    off, z = off[order], z[order]
    d = z - level
    out = []
    for i in range(len(d) - 1):
        a, b = d[i], d[i + 1]
        if a == 0.0:
            out.append(float(off[i]))
        if (a < 0.0 and b > 0.0) or (a > 0.0 and b < 0.0):
            t = a / (a - b)
            out.append(float(off[i] + t * (off[i + 1] - off[i])))
    if len(d) and d[-1] == 0.0:
        out.append(float(off[-1]))
    return out


def _optimize_step_midpoint(
    offset,
    elevation,
    zmid,
    zlo,
    zhi,
    *,
    n_candidates=81,
    min_crossings=2,
):
    """Nudge a step midpoint within ``(zlo, zhi)`` to reduce crossings.

    The candidate level is searched only inside its own step interval so it
    never crosses adjacent step boundaries.
    """
    if not np.isfinite(zlo) or not np.isfinite(zhi) or zhi <= zlo:
        return float(zmid)

    eps = max(1e-9, 1e-6 * (zhi - zlo))
    lo, hi = zlo + eps, zhi - eps
    if hi <= lo:
        return float(0.5 * (zlo + zhi))

    z0 = float(zmid) if np.isfinite(zmid) else float(0.5 * (lo + hi))
    z0 = float(np.clip(z0, lo, hi))

    n = max(5, int(n_candidates))
    cands = np.linspace(lo, hi, n)
    cands = np.unique(np.concatenate((cands, np.array([z0], dtype=float))))
    counts = np.array([
        len(_water_line_crossings(offset, elevation, float(z)))
        for z in cands
    ], dtype=int)

    # Prefer levels with at least one bank pair when available; then minimise
    # crossings and, on ties, keep the level closest to the original midpoint.
    valid = counts >= int(min_crossings)
    if np.any(valid):
        idx_all = np.where(valid)[0]
        best_local = min(
            range(len(idx_all)),
            key=lambda j: (counts[idx_all[j]], abs(cands[idx_all[j]] - z0)),
        )
        return float(cands[idx_all[best_local]])

    best = min(range(len(cands)), key=lambda i: (counts[i], abs(cands[i] - z0)))
    return float(cands[best])


def _adjust_step_midpoints_for_crossings(offset, elevation, boundaries, mids):
    """Adjust each step midpoint to reduce multi-crossing artifacts."""
    boundaries = np.asarray(boundaries, dtype=float)
    if len(boundaries) < 2:
        return []
    out = []
    for k in range(len(boundaries) - 1):
        zlo, zhi = float(boundaries[k]), float(boundaries[k + 1])
        zmid = float(mids[k]) if k < len(mids) else 0.5 * (zlo + zhi)
        out.append(_optimize_step_midpoint(offset, elevation, zmid, zlo, zhi))
    return out


def compute_flat_surface_edges(
    samples_df: pd.DataFrame,
    *,
    source_id: int = 0,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    offset_col: str = "offset_m",
    elevation_col: str = "elevation",
    x_col: str = "x",
    y_col: str = "y",
    stage_step: float = 0.1,
    min_step_height_m: float = 2.0,
    min_step_width_px: int = 3,
):
    """Terrace / floodplain edges from flow-width step midpoints, per station.

    For each cross section, floods the profile (:func:`cross_section_hydraulics`),
    combines the width steps (:func:`combine_width_steps`), nudges each
    step-mid elevation within its own step bounds to reduce excessive
    water-line crossings, then records every crossing of the profile
    (:func:`_water_line_crossings`) with interpolated world coordinates.

    Returns:
        DataFrame with columns ``source_id, station_seq, step, stage, offset_m,
        x, y`` (one row per crossing).
    """
    ch = samples_df[samples_df[source_col] == source_id]
    cols = ["source_id", "station_seq", "step", "stage", "offset_m", "x", "y"]
    rows = []
    for sseq in sorted(ch[station_col].unique()):
        s = ch[ch[station_col] == sseq].sort_values(offset_col)
        off = s[offset_col].to_numpy(dtype=float)
        elev = s[elevation_col].to_numpy(dtype=float)
        xw = s[x_col].to_numpy(dtype=float)
        yw = s[y_col].to_numpy(dtype=float)
        hyd = cross_section_hydraulics(off, elev, stage_step=stage_step)
        if not len(hyd["stage"]):
            continue
        boundaries, mids = combine_width_steps(
            hyd["stage"], hyd["width"], hyd["ds"],
            min_step_height_m=min_step_height_m,
            min_step_width_px=min_step_width_px)
        mids = _adjust_step_midpoints_for_crossings(off, elev, boundaries, mids)
        okxy = np.isfinite(off) & np.isfinite(xw) & np.isfinite(yw)
        offc, xc, yc = off[okxy], xw[okxy], yw[okxy]
        for step_idx, zmid in enumerate(mids):
            for xoff in _water_line_crossings(off, elev, zmid):
                rows.append(dict(
                    source_id=source_id, station_seq=int(sseq), step=step_idx,
                    stage=float(zmid), offset_m=float(xoff),
                    x=float(np.interp(xoff, offc, xc)),
                    y=float(np.interp(xoff, offc, yc))))
    return pd.DataFrame(rows, columns=cols)


def plot_valley_cross_section_hydraulics(
    samples_df: pd.DataFrame,
    *,
    source_id: int = 0,
    boundaries_df: "pd.DataFrame | None" = None,
    hydraulics_df: "pd.DataFrame | None" = None,
    source_col: str = "source_id",
    station_col: str = "station_seq",
    offset_col: str = "offset_m",
    elevation_col: str = "elevation",
    in_valley_col: str = "in_valley",
    along_col: str = "along_distance_m",
    max_panels: int = 20,
    n_cols: int = 4,
):
    """Flow-width step surfaces & terrace/floodplain edges, per cross section.

    VISUALISATION ONLY — reads the values COMPUTED BY THE C++ ENGINE
    (lsdtt-valley-centreline) and plots them; it does NOT recompute the
    hydraulic geometry. (The Python ``cross_section_hydraulics`` /
    ``combine_width_steps`` / ``detect_hydraulic_flats`` helpers in this module
    were prototypes only; the C++ engine is the single source of truth.)

    Two side-by-side panels per station: (left) the elevation profile with the
    C++ terrace/floodplain edges (magenta diamonds at each ``hydraulic_edge``'s
    offset/stage), (right) the flow-width curve ``W(h)`` ROTATED so
    stage/elevation is the vertical axis. Edge stages are drawn as dashed
    magenta lines on both panels.

    Args:
        samples_df: per-sample cross-section table
            (``_valley_cross_section_samples.csv``).
        source_id: channel source to plot.
        boundaries_df: the C++ boundary points
            (``_valley_cross_section_boundaries.fgb`` read into a DataFrame with
            a ``boundary_type`` column). Rows with ``boundary_type ==
            "hydraulic_edge"`` supply the edge offsets and stages. If None, no
            edges are drawn.
        hydraulics_df: the C++ hydraulic-depth table
            (``_valley_cross_section_hydraulics.csv``, written when
            ``write_cross_section_hydraulics`` is on). Supplies the ``W(h)``
            curve for the right panel. If None, the width panel is left blank.
        max_panels, n_cols: layout (``n_cols`` stations per row).

    Returns:
        ``(fig, results)`` where ``results`` maps station -> list of C++ edge
        stages (m).
    """
    plt = _import_mpl()
    from matplotlib.lines import Line2D

    ch = samples_df[samples_df[source_col] == source_id]
    stations = sorted(ch[station_col].unique())
    if not stations:
        raise ValueError(
            f"No cross-section stations for {source_col}={source_id}.")
    if max_panels and max_panels > 0 and len(stations) > max_panels:
        idx = np.linspace(0, len(stations) - 1, max_panels).round().astype(int)
        stations = [stations[i] for i in sorted(set(idx))]

    # C++ hydraulic edges for this channel (magenta diamonds): offset + stage.
    edges_ch = None
    if boundaries_df is not None and len(boundaries_df) \
            and "boundary_type" in boundaries_df.columns:
        edges_ch = boundaries_df[
            (boundaries_df[source_col] == source_id)
            & (boundaries_df["boundary_type"] == "hydraulic_edge")]
    # C++ width curve W(h) for this channel.
    hyd_ch = None
    if hydraulics_df is not None and len(hydraulics_df):
        hyd_ch = hydraulics_df[hydraulics_df[source_col] == source_id]

    n = len(stations)
    n_cols = max(1, min(n_cols, n))
    n_block_rows = (n + n_cols - 1) // n_cols
    fig, axarr = plt.subplots(
        n_block_rows, 2 * n_cols,
        figsize=(4.2 * n_cols, 3.0 * n_block_rows), squeeze=False,
        gridspec_kw={"width_ratios": [3.0, 1.0] * n_cols})

    xmax = float(ch[offset_col].abs().max())
    results = {}

    for k, sseq in enumerate(stations):
        br, col = k // n_cols, k % n_cols
        eax = axarr[br][2 * col]
        wax = axarr[br][2 * col + 1]

        s = ch[ch[station_col] == sseq].sort_values(offset_col)
        off = s[offset_col].to_numpy(dtype=float)
        elev = s[elevation_col].to_numpy(dtype=float)
        inval = s[in_valley_col].to_numpy(dtype=float)
        along = float(s[along_col].iloc[0])

        # C++-computed edges for this station: (offset, stage) pairs.
        edge_off, edge_stage = [], []
        if edges_ch is not None:
            es = edges_ch[edges_ch[station_col] == sseq]
            edge_off = es[offset_col].to_numpy(dtype=float).tolist()
            edge_stage = es["stage_m"].to_numpy(dtype=float).tolist()
        mids = sorted(set(round(z, 6) for z in edge_stage))
        results[int(sseq)] = list(map(float, mids))

        efin = elev[np.isfinite(elev)]
        ymin, ymax = (float(efin.min()), float(efin.max())) if len(efin) \
            else (0.0, 1.0)
        pad = 0.05 * (ymax - ymin + 1e-6)

        # elevation profile + C++ terrace/floodplain edges
        eax.plot(off, elev, "-", color="0.35", lw=1.1, zorder=2)
        eax.fill_between(off, ymin - pad, elev, where=(inval > 0),
                         color="tab:blue", alpha=0.15, step="mid", zorder=1)
        for zmid in mids:
            eax.axhline(zmid, color="magenta", lw=0.6, ls="--", alpha=0.5,
                        zorder=2)
        if edge_off:
            eax.plot(edge_off, edge_stage, linestyle="none", marker="D", ms=5,
                     mew=0.7, mfc="magenta", mec="black", zorder=5)
        eax.set_ylim(ymin - pad, ymax + pad)
        eax.set_xlim(-xmax, xmax)
        eax.set_title(f"station {int(sseq)} — {along / 1000:.2f} km",
                      fontsize=8)
        eax.tick_params(labelsize=6)
        if col == 0:
            eax.set_ylabel("elev (m)", fontsize=7)
        if br == n_block_rows - 1:
            eax.set_xlabel("offset (m)", fontsize=7)

        # flow width vs stage, rotated (width on x, stage/elevation on y)
        if hyd_ch is not None:
            hs = hyd_ch[hyd_ch[station_col] == sseq].sort_values("stage_m")
            if len(hs):
                wax.plot(hs["width_m"].to_numpy(dtype=float),
                         hs["stage_m"].to_numpy(dtype=float),
                         "-", color="tab:blue", lw=1.0)
        for zmid in mids:
            wax.axhline(zmid, color="magenta", lw=0.8, ls="--", alpha=0.7)
        wax.set_ylim(ymin - pad, ymax + pad)
        wax.tick_params(labelsize=6, labelleft=False)
        wax.set_title("W(h)", fontsize=7)
        if br == n_block_rows - 1:
            wax.set_xlabel("width (m)", fontsize=7)

    for k in range(n, n_block_rows * n_cols):
        br, col = k // n_cols, k % n_cols
        axarr[br][2 * col].set_axis_off()
        axarr[br][2 * col + 1].set_axis_off()

    handles = [
        Line2D([0], [0], color="0.35", lw=1.1, label="elevation"),
        Line2D([0], [0], color="tab:blue", lw=1.0, label="flow width W(h)"),
        Line2D([0], [0], color="magenta", lw=0, marker="D", ms=6,
               mec="black", label="terrace / floodplain edge"),
    ]
    fig.legend(handles=handles, loc="upper center", fontsize=8, frameon=False,
               ncol=3, bbox_to_anchor=(0.5, 0.998))
    fig.suptitle(
        f"channel {source_id}: flow-width step surfaces & terrace edges",
        fontsize=11, y=1.03)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    return fig, results
