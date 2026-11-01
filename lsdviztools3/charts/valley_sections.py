"""Valley cross sections: the profile, the chosen valley, and the height field.

One panel per station along a channel. Each shows the ground surface across the
section, the extent picked out as valley floor, and - when a field raster is
supplied - the height above nearest drainage at every sample.

WHAT THE FIGURE IS FOR. A valley width in a table is not checkable. Seeing the
section is: a width that stops short of an obvious break in slope, or runs on
past one, is visible immediately, and so is a valley whose floor the ladder has
cut at the wrong height. This is the figure that says whether the numbers
downstream mean anything.

TWO THINGS IT MARKS THAT ARE EASY TO ASSUME AWAY.

The centreline is at offset 0 by construction, and the LOW POINT of the section
usually is not. The centreline is a medial axis - it runs down the middle of the
valley's planform - while the thalweg runs wherever the water goes, which on a
meandering river is hard against the outside of each bend. The offset between
them is real and is worth seeing; a large one is not an error, but it does mean
"distance from the centreline" and "distance from the river" are different
quantities and should not be used interchangeably.

The field trace (HAND) is drawn against its own axis on the right, because it
shares no units with elevation. Where the valley floor is picked correctly the
field is flat and low across the shaded span and climbs at both edges; a field
that is still climbing inside the shading is the clearest single sign that the
break was set too high.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Band colours. Ordered dark-to-light up the ladder so a tread reads as a step
# in shade as well as in position, and kept away from the green/red the section
# sides use.
_BAND_COLOURS = ["#08519c", "#e6550d", "#6a51a3", "#31a354",
                 "#d94801", "#756bb1", "#78c679"]
_N_BAND_COLOURS = len(_BAND_COLOURS)

# Surface class strips (lsdtt-fluvial-surfaces <stem>_surface_class.tif:
# 1 floodplain, 2 terrace); unclassified cells are left pale grey.
_CLASS_COLOURS = {1: "#2a78d6", 2: "#eb6834"}
_CLASS_NONE = "#e4e3df"
_STRIP_HEIGHT = 0.07   # axes fraction per strip


def _sample_raster(path, xs, ys):
    """Nearest-cell sample of a raster at projected coordinates.

    Args:
        path: Raster path.
        xs: Easting array.
        ys: Northing array.

    Returns:
        Float array of sampled values, NaN off-grid or at nodata.
    """
    # rasterio rather than GDAL's Python bindings: it is what the rest of
    # lsdviztools3 uses and what this environment actually ships.
    import rasterio

    with rasterio.open(str(path)) as src:
        data = src.read(1).astype("float64")
        nodata = src.nodata
        inv = ~src.transform
        cols, rows = inv * (np.asarray(xs, dtype=float),
                            np.asarray(ys, dtype=float))
        cols = np.floor(cols).astype(int)
        rows = np.floor(rows).astype(int)
        h, w = data.shape
        ok = (rows >= 0) & (rows < h) & (cols >= 0) & (cols < w)
        out = np.full(len(cols), np.nan)
        out[ok] = data[rows[ok], cols[ok]]
        if nodata is not None:
            out[out == nodata] = np.nan
    return out


def plot_valley_sections(
    samples_csv: str | Path,
    out_png: str | Path,
    source_id: int = 0,
    field_raster: str | Path | None = None,
    field_label: str = "HAND (m)",
    band_raster: str | Path | None = None,
    class_rasters=None,
    max_panels: int = 30,
    n_cols: int = 5,
    title: str | None = None,
    dpi: int = 150,
) -> Path:
    """Plot valley cross sections for one channel.

    Args:
        samples_csv: A `<stem>_valley_cross_section_samples.csv`.
        out_png: Output PNG path.
        source_id: Which channel to plot.
        field_raster: Optional raster sampled at each point (normally HAND).
        field_label: Axis label for that field.
        band_raster: Optional ladder band-index raster. When given, the field
            points are coloured by band instead of one colour, which is what
            shows whether a band is a terrace tread or a slice of hillslope.
            Band indices are PER WINDOW, so the same colour at two stations far
            apart is not the same height.
        class_rasters: Optional list of (label, path) surface-class rasters
            (1 floodplain, 2 terrace), each drawn as a strip under the ground
            profile so methods can be compared section by section.
        max_panels: Cap on panels, sampled evenly downstream.
        n_cols: Panels per row.
        title: Figure title.
        dpi: Output resolution.

    Returns:
        The path written.
    """
    df = pd.read_csv(samples_csv)
    df = df[df["source_id"] == source_id]
    if df.empty:
        raise ValueError(f"no samples for source_id {source_id}")

    stations = sorted(df["station_seq"].unique())
    if len(stations) > max_panels:
        idx = np.linspace(0, len(stations) - 1, max_panels).astype(int)
        stations = [stations[i] for i in idx]

    bands = None
    if band_raster is not None:
        bands = _sample_raster(band_raster, df["x"].to_numpy(),
                               df["y"].to_numpy())

    field = None
    if field_raster is not None:
        field = _sample_raster(field_raster, df["x"].to_numpy(),
                               df["y"].to_numpy())
        df = df.assign(_field=field)
    if bands is not None:
        df = df.assign(_band=bands)
    class_rasters = list(class_rasters or [])
    for j, (_, path) in enumerate(class_rasters):
        df[f"_class{j}"] = _sample_raster(path, df["x"].to_numpy(),
                                          df["y"].to_numpy())

    n_rows = int(np.ceil(len(stations) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.5 * n_cols, 2.6 * n_rows),
                             squeeze=False)

    for k, st in enumerate(stations):
        ax = axes[k // n_cols][k % n_cols]
        s = df[df["station_seq"] == st].sort_values("offset_m")
        off = s["offset_m"].to_numpy()
        elev = s["elevation"].to_numpy()

        ax.plot(off, elev, "-", color="#555555", lw=1.2, zorder=3)

        if class_rasters and np.isfinite(elev).any():
            # Room under the profile for one strip per method.
            e0, e1 = np.nanmin(elev), np.nanmax(elev)
            frac = _STRIP_HEIGHT * len(class_rasters) + 0.04
            pad = (e1 - e0) * frac / (1.0 - frac - 0.05)
            ax.set_ylim(e0 - pad, e1 + 0.05 * (e1 - e0))
            tr = ax.get_xaxis_transform()
            # Each sample covers half the gap to its neighbours.
            edges = np.r_[off[0], 0.5 * (off[1:] + off[:-1]), off[-1]]
            for j, (label, _) in enumerate(class_rasters):
                y0 = 0.01 + _STRIP_HEIGHT * (len(class_rasters) - 1 - j)
                cv = s[f"_class{j}"].to_numpy()
                for i in range(len(off)):
                    c = _CLASS_COLOURS.get(int(cv[i]) if np.isfinite(cv[i]) else 0,
                                           _CLASS_NONE)
                    ax.fill_between([edges[i], edges[i + 1]], y0,
                                    y0 + _STRIP_HEIGHT * 0.85, color=c, lw=0,
                                    transform=tr, zorder=2)
                ax.text(off[0], y0 + _STRIP_HEIGHT * 0.42, f" {label}",
                        transform=tr, fontsize=5.5, va="center", ha="left",
                        zorder=6)

        # The picked valley, shaded between its two sides.
        lo = s["left_offset_m"].iloc[0]
        hi = s["right_offset_m"].iloc[0]
        width = s["valley_width_m"].iloc[0]
        inside = (off >= lo) & (off <= hi)
        if inside.any():
            ax.fill_between(off[inside], elev[inside], elev[inside].max(),
                            color="#a6cee3", alpha=0.55, zorder=1)
        ax.axvline(lo, ls="--", color="#1a9850", lw=1.2, zorder=4)
        ax.axvline(hi, ls="--", color="#d73027", lw=1.2, zorder=4)
        ax.axvline(0.0, ls=":", color="black", lw=0.9, zorder=4)

        # The low point, which is usually NOT at offset 0 - see the module
        # docstring. Annotated with its offset so the size of the difference is
        # readable rather than merely visible.
        if len(elev) and np.isfinite(elev).any():
            j = int(np.nanargmin(elev))
            ax.plot(off[j], elev[j], "v", color="#08519c", ms=6, zorder=5)
            ax.annotate(f"{off[j]:+.0f} m", xy=(off[j], elev[j]),
                        xytext=(2, -11), textcoords="offset points",
                        fontsize=6.5, color="#08519c")

        if field is not None:
            axf = ax.twinx()
            if bands is not None:
                bv = s["_band"].to_numpy()
                fv = s["_field"].to_numpy()
                for bi in range(_N_BAND_COLOURS):
                    sel = bv == bi
                    if sel.any():
                        axf.plot(off[sel], fv[sel], ".",
                                 color=_BAND_COLOURS[bi], ms=3.0, zorder=2)
                sel = ~np.isin(bv, np.arange(_N_BAND_COLOURS))
                if sel.any():
                    axf.plot(off[sel], fv[sel], ".", color="#bbbbbb",
                             ms=2.0, zorder=2)
            else:
                axf.plot(off, s["_field"].to_numpy(), ".", color="#e6550d",
                         ms=2.5, zorder=2)
            axf.tick_params(axis="y", labelsize=6, colors="#e6550d")
            if k % n_cols == n_cols - 1:
                axf.set_ylabel(field_label, fontsize=7, color="#e6550d")

        ax.set_title(f"{s['along_distance_m'].iloc[0]/1000.0:.2f} km  "
                     f"(w={width:.0f} m)", fontsize=8)
        ax.tick_params(labelsize=7)
        if k % n_cols == 0:
            ax.set_ylabel("elev (m)", fontsize=8)
        if k // n_cols == n_rows - 1:
            ax.set_xlabel("offset (m)", fontsize=8)

    for k in range(len(stations), n_rows * n_cols):
        axes[k // n_cols][k % n_cols].axis("off")

    handles = [
        plt.Line2D([], [], color="#555555", lw=1.2, label="elevation"),
        plt.Line2D([], [], color="#1a9850", ls="--", label="left side"),
        plt.Line2D([], [], color="#d73027", ls="--", label="right side"),
        plt.Line2D([], [], color="black", ls=":", label="centreline"),
        plt.Line2D([], [], color="#08519c", marker="v", ls="none",
                   label="section low point"),
    ]
    if class_rasters:
        handles += [plt.Rectangle((0, 0), 1, 1, fc=_CLASS_COLOURS[1],
                                  label="floodplain"),
                    plt.Rectangle((0, 0), 1, 1, fc=_CLASS_COLOURS[2],
                                  label="terrace"),
                    plt.Rectangle((0, 0), 1, 1, fc=_CLASS_NONE,
                                  label="neither")]
    if field is not None and bands is None:
        handles.append(plt.Line2D([], [], color="#e6550d", marker=".",
                                  ls="none", label=field_label))
    if field is not None and bands is not None:
        for bi in range(_N_BAND_COLOURS):
            handles.append(plt.Line2D([], [], color=_BAND_COLOURS[bi],
                                      marker=".", ls="none",
                                      label=f"band {bi}"))
        handles.append(plt.Line2D([], [], color="#bbbbbb", marker=".",
                                  ls="none", label="unbanded"))
    # Legend below the title rather than beside it: at five or six entries it
    # is wider than the space left over on the title row and overlaps it.
    fig.suptitle(title or
                 f"channel {source_id}: valley cross sections "
                 "(orthogonal to centreline)", y=0.995, fontsize=11)
    fig.legend(handles=handles, loc="upper center",
               bbox_to_anchor=(0.5, 0.972), ncol=len(handles),
               fontsize=8, frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png


def summarise_sections(samples_csv: str | Path,
                       source_id: int = 0) -> pd.DataFrame:
    """Per-station width and centreline-to-thalweg offset.

    Args:
        samples_csv: A `<stem>_valley_cross_section_samples.csv`.
        source_id: Which channel.

    Returns:
        One row per station.
    """
    df = pd.read_csv(samples_csv)
    df = df[df["source_id"] == source_id]
    rows = []
    for st, s in df.groupby("station_seq"):
        s = s.sort_values("offset_m")
        elev = s["elevation"].to_numpy()
        off = s["offset_m"].to_numpy()
        j = int(np.nanargmin(elev)) if np.isfinite(elev).any() else 0
        rows.append({
            "station_seq": st,
            "along_distance_m": s["along_distance_m"].iloc[0],
            "valley_width_m": s["valley_width_m"].iloc[0],
            "left_offset_m": s["left_offset_m"].iloc[0],
            "right_offset_m": s["right_offset_m"].iloc[0],
            "low_point_offset_m": off[j],
            "low_point_elevation_m": elev[j],
        })
    return pd.DataFrame(rows)
