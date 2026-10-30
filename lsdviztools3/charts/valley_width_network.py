"""Valley widths for every channel: the map of sections and the width profile.

Reads the outputs of one lsdtt3 ``lsdtt-valley-centreline`` run made with
``run_valley_width_network: true``, by its stem (write_path + write_fname):

  <stem>_valley_cross_sections.fgb      section lines, with segment_id,
                                        parent_segment_id, excluded_reason
  <stem>_valley_widths.csv              one row per station
  <stem>_valley_segments.csv            one row per channel segment
  <stem>_valley_centreline.fgb          one centreline per segment
  <stem>_valley_centreline_nodes.csv    centreline nodes with distance_m
  <stem>_valley_segment_channels.fgb    each segment's extracted D8 channel

Four figures:

``plot_valley_width_map``: each section's width line (left valley side to
right side) on the hillshade, coloured by excluded_reason, with the
centrelines and optionally the valley mask. The kept sections are the widths;
the excluded ones show what the tributary rules removed and why.

``plot_valley_widths_along_channel``: valley width against distance down one
segment's centreline (by default the trunk with the largest drainage area), kept
stations joined, excluded ones as hollow markers, tributary confluences marked,
and drainage area in a second panel below (its own axis, not a second y-scale).

``plot_channel_and_centreline``: each segment's extracted D8 channel and the
centreline its sections hang from, on the hillshade - where the two part (a
meandering channel, a channel hugging one valley wall) the sections follow the
centreline, not the channel.

``plot_centreline_profile``: down one segment's centreline, the elevation of
the DEM pixel under each centreline node against the centreline's baseline (the
reference the relative elevation is measured from), pins marked, and their
difference - the relative elevation along the centreline - in a panel below.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Fixed order: kept, then the two tributary rules; side_not_found is neutral
# grey (no width was measured). Excluded sections are also dashed / hollow, so
# identity does not rest on colour alone.
REASON_COLOURS = {
    "none": "#2a78d6",
    "near_parent_channel": "#eb6834",
    "tributary_crossing": "#1baf7a",
    "side_not_found": "#8a8a85",
    "not_in_valley": "#b4b2a9",
}
REASON_LABELS = {
    "none": "kept",
    "near_parent_channel": "excluded: tributary inside parent valley",
    "tributary_crossing": "excluded: runs up a tributary",
    "side_not_found": "excluded: valley side not found",
    "not_in_valley": "no valley at the station",
}
MASK_COLOUR = "#f2d48a"
CENTRELINE_COLOUR = "#3d3d3a"


def _read(stem, suffix):
    p = Path(f"{stem}{suffix}")
    if not p.exists():
        raise FileNotFoundError(
            f"{p} (run lsdtt-valley-centreline with run_valley_width_network: true).")
    return p


def largest_trunk(widths: pd.DataFrame | None, segments: pd.DataFrame) -> int:
    """The trunk segment (no parent) reaching the largest drainage area.

    Falls back to the trunk with the most channel pixels when there is no widths
    table (a run without cross sections) or no station has a drainage area.
    """
    trunks = segments[segments["parent_segment_id"] < 0]["segment_id"]
    if trunks.empty:
        raise ValueError("no trunk segment in the segments table")
    w = (widths[widths["segment_id"].isin(trunks)] if widths is not None
         else pd.DataFrame({"segment_id": [], "drainage_area_m2": []}))
    area = w.groupby("segment_id")["drainage_area_m2"].max().dropna()
    if not area.empty:
        return int(area.idxmax())
    s = segments[segments["segment_id"].isin(trunks)]
    return int(s.loc[s["n_channel_pixels"].idxmax(), "segment_id"])


def confluence_distances(stem, segment_id: int) -> pd.DataFrame:
    """Where each direct tributary of ``segment_id`` joins it.

    Returns segment_id and distance_m (along the parent's centreline from its
    head) of each child, from the child's end point (its confluence) projected
    onto the nearest parent centreline node.
    """
    segs = pd.read_csv(_read(stem, "_valley_segments.csv"))
    nodes = pd.read_csv(_read(stem, "_valley_centreline_nodes.csv"))
    pn = nodes[nodes["source_id"] == segment_id]
    kids = segs[segs["parent_segment_id"] == segment_id]
    if pn.empty or kids.empty:
        return pd.DataFrame({"segment_id": [], "distance_m": []})
    px, py, pd_ = pn["x"].to_numpy(), pn["y"].to_numpy(), pn["distance_m"].to_numpy()
    rows = []
    for _, k in kids.iterrows():
        i = int(np.argmin((px - k["end_x"]) ** 2 + (py - k["end_y"]) ** 2))
        rows.append((int(k["segment_id"]), float(pd_[i])))
    return pd.DataFrame(rows, columns=["segment_id", "distance_m"])


def plot_valley_widths_along_channel(
    stem,
    out_png,
    segment_id: int | None = None,
    show_excluded: bool = True,
    title: str | None = None,
    figsize=(9.0, 6.0),
    dpi: int = 200,
) -> Path:
    """Valley width down one segment, with drainage area below.

    Args:
        stem: lsdtt-valley-centreline output stem.
        out_png: Output PNG.
        segment_id: Segment to plot; default the trunk with the largest
            drainage area.
        show_excluded: Also plot excluded stations (hollow markers).
        title: Figure title.
        figsize: Figure size (inches).
        dpi: Output resolution.

    Returns:
        The path written.
    """
    widths = pd.read_csv(_read(stem, "_valley_widths.csv"))
    segs = pd.read_csv(_read(stem, "_valley_segments.csv"))
    seg = largest_trunk(widths, segs) if segment_id is None else int(segment_id)
    w = widths[widths["segment_id"] == seg].sort_values("along_distance_m")
    if w.empty:
        raise ValueError(f"segment {seg} has no stations")
    km = w["along_distance_m"] / 1000.0
    kept = w["excluded_reason"] == "none"

    fig, (ax, ax_a) = plt.subplots(
        2, 1, sharex=True, figsize=figsize,
        gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.08})

    ax.plot(km[kept], w.loc[kept, "valley_width_m"], "-o", color=REASON_COLOURS["none"],
            lw=1.5, ms=3.5, label=f"{REASON_LABELS['none']} ({int(kept.sum())})", zorder=3)
    if show_excluded:
        for reason in ("near_parent_channel", "tributary_crossing", "side_not_found"):
            m = w["excluded_reason"] == reason
            if not m.any():
                continue
            ax.plot(km[m], w.loc[m, "valley_width_m"], "o", mfc="none",
                    mec=REASON_COLOURS[reason], mew=1.2, ms=5,
                    label=f"{REASON_LABELS[reason]} ({int(m.sum())})", zorder=2)

    # Stations with no valley (the channel runs outside the mask): ticks on
    # the axis, so the reach is visibly unmeasured rather than zero-width.
    nv = w["excluded_reason"] == "not_in_valley"
    if nv.any():
        ax.plot(km[nv], np.zeros(int(nv.sum())), "|", color=REASON_COLOURS["not_in_valley"],
                ms=7, mew=1.0, clip_on=False,
                label=f"{REASON_LABELS['not_in_valley']} ({int(nv.sum())})", zorder=2)

    conf = confluence_distances(stem, seg)
    for i, d in enumerate(conf["distance_m"] / 1000.0):
        for a in (ax, ax_a):
            a.axvline(d, color="#b4b2a9", lw=0.8, ls="--", zorder=1,
                      label="tributary confluence" if (i == 0 and a is ax) else None)

    ax.set_ylabel("valley width (m)")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", color="#e6e5e0", lw=0.6)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.9)
    ax.set_title(title or f"Valley width along segment {seg}"
                 + (" (largest trunk)" if segment_id is None else ""))

    area = w["drainage_area_m2"] / 1e6
    ax_a.plot(km, area, color="#5f5e5a", lw=1.5)
    if (area > 0).any():
        ax_a.set_yscale("log")
    ax_a.set_ylabel("drainage\narea (km²)")
    ax_a.set_xlabel("distance down the centreline from the channel head (km)")
    ax_a.grid(axis="y", color="#e6e5e0", lw=0.6)
    for a in (ax, ax_a):
        for side in ("top", "right"):
            a.spines[side].set_visible(False)

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png


def plot_valley_width_map(
    stem,
    hillshade,
    out_png,
    valley_mask=None,
    bbox=None,
    segment_ids=None,
    show_excluded: bool = True,
    title: str | None = None,
    figsize=(10.0, 11.0),
    dpi: int = 200,
) -> Path:
    """Section width lines on the hillshade, coloured by excluded_reason.

    Args:
        stem: lsdtt-valley-centreline output stem.
        hillshade: Hillshade raster on the run's grid.
        out_png: Output PNG.
        valley_mask: Optional valley mask raster (cells > 0.5), drawn pale.
        bbox: Optional (xmin, ymin, xmax, ymax); default the sections' extent.
        segment_ids: Optional iterable of segments to draw (default all).
        show_excluded: Also draw excluded sections (dashed).
        title: Figure title.
        figsize: Figure size (inches).
        dpi: Output resolution.

    Returns:
        The path written.
    """
    import geopandas as gpd
    import rasterio
    from matplotlib.collections import LineCollection
    from rasterio.plot import plotting_extent

    with rasterio.open(str(hillshade)) as src:
        hs = src.read(1, masked=True).astype("float64").filled(np.nan)
        extent = plotting_extent(src)

    xs = gpd.read_file(str(_read(stem, "_valley_cross_sections.fgb")))
    cl = gpd.read_file(str(_read(stem, "_valley_centreline.fgb")))
    if segment_ids is not None:
        keep = set(int(s) for s in segment_ids)
        xs = xs[xs["segment_id"].isin(keep)]
        cl = cl[cl["source_id"].isin(keep)]
    if "excluded_reason" not in xs.columns:
        raise ValueError("no excluded_reason column: not a run_valley_width_network run")
    # Stations with no valley have no width line to draw.
    n_no_valley = int((xs["excluded_reason"] == "not_in_valley").sum())
    xs = xs[xs["excluded_reason"] != "not_in_valley"]

    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(hs, cmap="gray", extent=extent, origin="upper", zorder=0)
    if valley_mask is not None:
        with rasterio.open(str(valley_mask)) as src:
            vm = src.read(1, masked=True).filled(0)
            vext = plotting_extent(src)
        # An RGBA image rather than a masked array through a colormap: the
        # downsampling of a large masked array drops whole patches of it.
        rgba = np.zeros(vm.shape + (4,), dtype=np.float32)
        rgba[..., :3] = matplotlib.colors.to_rgb(MASK_COLOUR)
        rgba[..., 3] = np.where(vm > 0.5, 0.6, 0.0)
        ax.imshow(rgba, extent=vext, origin="upper", interpolation="nearest",
                  zorder=1)
    cl.plot(ax=ax, color=CENTRELINE_COLOUR, lw=0.8, zorder=2)

    handles = []
    order = ["side_not_found", "near_parent_channel", "tributary_crossing", "none"]
    for z, reason in enumerate(order):
        if reason != "none" and not show_excluded:
            continue
        sub = xs[xs["excluded_reason"] == reason]
        if sub.empty:
            continue
        segs = [[(r.left_x, r.left_y), (r.right_x, r.right_y)] for r in sub.itertuples()]
        kept = reason == "none"
        ax.add_collection(LineCollection(
            segs, colors=REASON_COLOURS[reason], linewidths=1.1 if kept else 0.9,
            linestyles="solid" if kept else (0, (3, 2)), zorder=3 + z))
        handles.append(plt.Line2D([], [], color=REASON_COLOURS[reason],
                                  lw=1.4 if kept else 1.1, ls="-" if kept else "--",
                                  label=f"{REASON_LABELS[reason]} ({len(sub)})"))
    handles.reverse()
    handles.append(plt.Line2D([], [], color=CENTRELINE_COLOUR, lw=0.8,
                              label="segment centreline"))
    if n_no_valley:
        handles.append(plt.Line2D([], [], color="none",
                                  label=f"(+{n_no_valley} stations with no valley, not drawn)"))
    if valley_mask is not None:
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=MASK_COLOUR, alpha=0.6,
                                     label="valley mask"))
    # Below the map, so it never sits on the sections.
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.07),
              ncol=2, fontsize=8, frameon=False)

    if bbox is not None:
        x0, y0, x1, y1 = bbox
    else:
        x0 = float(np.nanmin(np.r_[xs["left_x"], xs["right_x"]]))
        x1 = float(np.nanmax(np.r_[xs["left_x"], xs["right_x"]]))
        y0 = float(np.nanmin(np.r_[xs["left_y"], xs["right_y"]]))
        y1 = float(np.nanmax(np.r_[xs["left_y"], xs["right_y"]]))
        mx, my = 0.05 * max(x1 - x0, 1.0), 0.05 * max(y1 - y0, 1.0)
        x0, x1, y0, y1 = x0 - mx, x1 + mx, y0 - my, y1 + my
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xlabel("easting (m)")
    ax.set_ylabel("northing (m)")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set_title(title or "Valley cross sections by excluded_reason")

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png


CHANNEL_COLOUR = "#2a78d6"
CENTRELINE_MAP_COLOUR = "#eb6834"


def plot_channel_and_centreline(
    stem,
    hillshade,
    out_png,
    valley_mask=None,
    bbox=None,
    segment_ids=None,
    title: str | None = None,
    figsize=(10.0, 11.0),
    dpi: int = 200,
) -> Path:
    """Each segment's extracted channel and its valley centreline on the hillshade.

    Args:
        stem: lsdtt-valley-centreline output stem.
        hillshade: Hillshade raster on the run's grid.
        out_png: Output PNG.
        valley_mask: Optional valley mask raster (cells > 0.5), drawn pale.
        bbox: Optional (xmin, ymin, xmax, ymax); default the extent of the
            drawn channels.
        segment_ids: Optional iterable of segments to draw (default all).
        title: Figure title.
        figsize: Figure size (inches).
        dpi: Output resolution.

    Returns:
        The path written.
    """
    import geopandas as gpd
    import rasterio
    from rasterio.plot import plotting_extent

    with rasterio.open(str(hillshade)) as src:
        hs = src.read(1, masked=True).astype("float64").filled(np.nan)
        extent = plotting_extent(src)
    ch = gpd.read_file(str(_read(stem, "_valley_segment_channels.fgb")))
    cl = gpd.read_file(str(_read(stem, "_valley_centreline.fgb")))
    if segment_ids is not None:
        keep = set(int(s) for s in segment_ids)
        ch = ch[ch["segment_id"].isin(keep)]
        cl = cl[cl["source_id"].isin(keep)]

    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(hs, cmap="gray", extent=extent, origin="upper", zorder=0)
    handles = []
    if valley_mask is not None:
        with rasterio.open(str(valley_mask)) as src:
            vm = src.read(1, masked=True).filled(0)
            vext = plotting_extent(src)
        rgba = np.zeros(vm.shape + (4,), dtype=np.float32)
        rgba[..., :3] = matplotlib.colors.to_rgb(MASK_COLOUR)
        rgba[..., 3] = np.where(vm > 0.5, 0.6, 0.0)
        ax.imshow(rgba, extent=vext, origin="upper", interpolation="nearest", zorder=1)
    ch.plot(ax=ax, color=CHANNEL_COLOUR, lw=1.0, zorder=2)
    cl.plot(ax=ax, color=CENTRELINE_MAP_COLOUR, lw=1.3, zorder=3)
    handles.append(plt.Line2D([], [], color=CHANNEL_COLOUR, lw=1.2,
                              label=f"extracted channel ({len(ch)} segments)"))
    handles.append(plt.Line2D([], [], color=CENTRELINE_MAP_COLOUR, lw=1.4,
                              label="valley centreline"))
    if valley_mask is not None:
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=MASK_COLOUR, alpha=0.6,
                                     label="valley mask"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.07),
              ncol=3, fontsize=8, frameon=False)

    if bbox is not None:
        x0, y0, x1, y1 = bbox
    else:
        x0, y0, x1, y1 = ch.total_bounds
        mx, my = 0.04 * max(x1 - x0, 1.0), 0.04 * max(y1 - y0, 1.0)
        x0, x1, y0, y1 = x0 - mx, x1 + mx, y0 - my, y1 + my
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xlabel("easting (m)")
    ax.set_ylabel("northing (m)")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set_title(title or "Extracted channel and valley centreline")

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png


def plot_centreline_profile(
    stem,
    out_png,
    segment_id: int | None = None,
    title: str | None = None,
    figsize=(9.0, 6.0),
    dpi: int = 200,
) -> Path:
    """Elevation under the centreline against its baseline, down one segment.

    Works on any lsdtt-valley-centreline run (it reads
    ``<stem>_valley_centreline_nodes.csv``). For a network run the default
    segment is the trunk reaching the largest drainage area; for a run on
    source points it is source 0.

    Args:
        stem: lsdtt-valley-centreline output stem.
        out_png: Output PNG.
        segment_id: Segment / source to plot.
        title: Figure title.
        figsize: Figure size (inches).
        dpi: Output resolution.

    Returns:
        The path written.
    """
    nodes = pd.read_csv(_read(stem, "_valley_centreline_nodes.csv"))
    seg_csv = Path(f"{stem}_valley_segments.csv")
    is_network = seg_csv.exists()
    if segment_id is None:
        if is_network:
            wcsv = Path(f"{stem}_valley_widths.csv")
            widths = pd.read_csv(wcsv) if wcsv.exists() else None
            segment_id = largest_trunk(widths, pd.read_csv(seg_csv))
            default_note = " (largest trunk)"
        else:
            segment_id = 0
            default_note = ""
    else:
        default_note = ""
    n = nodes[nodes["source_id"] == int(segment_id)].sort_values("node_seq")
    if n.empty:
        raise ValueError(f"no centreline nodes for segment {segment_id}")
    km = n["distance_m"] / 1000.0
    z = n["elevation"]
    b = n["baseline_elevation"]
    pins = n["is_pin"] == 1

    fig, (ax, ax_d) = plt.subplots(
        2, 1, sharex=True, figsize=figsize,
        gridspec_kw={"height_ratios": [3, 1.4], "hspace": 0.08})
    ax.plot(km, z, color="#8a8a85", lw=1.0, label="DEM under the centreline", zorder=2)
    ax.plot(km, b, color="#2a78d6", lw=1.8, label="centreline baseline", zorder=3)
    if pins.any():
        ax.plot(km[pins], b[pins], "o", ms=3, mfc="#eb6834", mec="none",
                label=f"centreline node on its channel ({int(pins.sum())})", zorder=4)
    if is_network:
        conf = confluence_distances(stem, int(segment_id))
        for i, d in enumerate(conf["distance_m"] / 1000.0):
            for a in (ax, ax_d):
                a.axvline(d, color="#b4b2a9", lw=0.8, ls="--", zorder=1,
                          label="tributary confluence" if (i == 0 and a is ax) else None)
    ax.set_ylabel("elevation (m)")
    ax.grid(axis="y", color="#e6e5e0", lw=0.6)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
    ax.set_title(title or f"Centreline elevation profile, segment {segment_id}{default_note}")

    ax_d.axhline(0.0, color="#5f5e5a", lw=0.8, zorder=1)
    ax_d.plot(km, z - b, color="#3d3d3a", lw=1.0, zorder=2)
    ax_d.set_ylabel("DEM - baseline\n(m)")
    ax_d.set_xlabel("distance down the centreline from the channel head (km)")
    ax_d.grid(axis="y", color="#e6e5e0", lw=0.6)
    for a in (ax, ax_d):
        for side in ("top", "right"):
            a.spines[side].set_visible(False)

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png
