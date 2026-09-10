"""Where the cross sections are, relative to the centreline and the channel.

The section plots show a profile against an offset from the centreline. They
cannot show whether that section crossed the river, or how far the centreline
sits from it, and both matter when reading a section: a low point away from
offset 0 means one thing if the channel is there and another if the section
missed the channel altogether.

This is the map that answers it. Three layers on the hillshade - the valley
centreline, the channel network, and the numbered sections - plus a marker on
any section that crosses NO channel cell. That last is the case worth catching:
the sections are cast orthogonal to the CENTRELINE, and where the centreline has
wandered from the river (a wide floor, a meander cut off from the medial axis,
the reach above the channel head) a section can run the full width of the valley
without meeting the channel at all. Its low point is then some other low ground,
and its height above nearest drainage is measured to a donor somewhere else
entirely.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_valley_section_map(
    hillshade: str | Path,
    sections_fgb: str | Path,
    out_png: str | Path,
    centreline_fgb: str | Path | None = None,
    channel_fgb: str | Path | None = None,
    samples_csv: str | Path | None = None,
    sides_fgb: str | Path | None = None,
    source_id: int = 0,
    label_every: int = 2,
    title: str | None = None,
    dpi: int = 200,
) -> Path:
    """Map the cross sections against the centreline and the channel.

    Args:
        hillshade: Hillshade raster for the base.
        sections_fgb: `<stem>_valley_cross_sections.fgb`.
        out_png: Output PNG path.
        centreline_fgb: Optional `<stem>_valley_centreline.fgb`.
        channel_fgb: Optional `<stem>_channel_network_lines.fgb`.
        samples_csv: Optional sample table; used to mark sections that cross
            no channel cell.
        sides_fgb: Optional `<stem>_valley_cross_section_sides.fgb`, the
            detected valley margins, dotted on the map.
        source_id: Which channel to draw.
        label_every: Label every Nth station.
        title: Figure title.
        dpi: Output resolution.

    Returns:
        The path written.
    """
    import geopandas as gpd
    import rasterio
    from rasterio.plot import plotting_extent

    with rasterio.open(str(hillshade)) as src:
        hs = src.read(1).astype("float64")
        nod = src.nodata
        if nod is not None:
            hs[hs == nod] = np.nan
        extent = plotting_extent(src)

    sec = gpd.read_file(str(sections_fgb))
    if "source_id" in sec.columns:
        sec = sec[sec["source_id"] == source_id]

    # Sections whose samples never touch a channel cell.
    missing = set()
    if samples_csv is not None:
        df = pd.read_csv(samples_csv)
        df = df[df["source_id"] == source_id]
        for st, s in df.groupby("station_seq"):
            if int(s["is_channel"].max()) == 0:
                missing.add(int(st))

    fig, ax = plt.subplots(figsize=(9, 11))
    ax.imshow(hs, cmap="gray", extent=extent, origin="upper", zorder=0)

    if channel_fgb is not None and Path(channel_fgb).exists():
        gpd.read_file(str(channel_fgb)).plot(
            ax=ax, color="#2171b5", lw=1.6, zorder=2)
    if centreline_fgb is not None and Path(centreline_fgb).exists():
        cl = gpd.read_file(str(centreline_fgb))
        if "source_id" in cl.columns:
            cl = cl[cl["source_id"] == source_id]
        cl.plot(ax=ax, color="#cb181d", lw=1.8, zorder=3)

    # Sections last, so they read on top; the ones with no channel in red.
    for _, row in sec.iterrows():
        st = int(row["station_seq"]) if "station_seq" in sec.columns else -1
        bad = st in missing
        gpd.GeoSeries([row.geometry]).plot(
            ax=ax, color=("#d94801" if bad else "#ffd92f"),
            lw=(1.8 if bad else 1.0), zorder=4)
        if st >= 0 and st % max(1, label_every) == 0:
            c = row.geometry.centroid
            ax.annotate(str(st), xy=(c.x, c.y), fontsize=6,
                        color=("#d94801" if bad else "black"), zorder=5)

    if sides_fgb is not None and Path(sides_fgb).exists():
        sd = gpd.read_file(str(sides_fgb))
        if "source_id" in sd.columns:
            sd = sd[sd["source_id"] == source_id]
        if len(sd):
            ax.plot(sd.geometry.x, sd.geometry.y, ".", ms=2.5,
                    color="#00bfff", zorder=5)

    handles = [
        plt.Line2D([], [], color="#cb181d", lw=1.8, label="valley centreline"),
        plt.Line2D([], [], color="#2171b5", lw=1.6, label="channel network"),
        plt.Line2D([], [], color="#ffd92f", lw=1.2, label="cross section"),
    ]
    if sides_fgb is not None and Path(sides_fgb).exists():
        handles.append(plt.Line2D([], [], color="#00bfff", marker=".",
                                  ls="none", label="valley side"))
    if missing:
        handles.append(plt.Line2D([], [], color="#d94801", lw=1.8,
                                  label="section crossing NO channel cell"))
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.9)

    # Crop to the sections: on a large tile the valley is a small part of the
    # DEM and the full extent renders it unreadable.
    b = sec.total_bounds
    mx = 0.08 * max(b[2] - b[0], 1.0)
    my = 0.08 * max(b[3] - b[1], 1.0)
    ax.set_xlim(b[0] - mx, b[2] + mx)
    ax.set_ylim(b[1] - my, b[3] + my)
    ax.set_aspect("equal")

    ax.set_xlabel("easting (m)")
    ax.set_ylabel("northing (m)")
    ax.set_title(title or
                 f"channel {source_id}: cross sections, centreline and channel")
    fig.tight_layout()
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png
