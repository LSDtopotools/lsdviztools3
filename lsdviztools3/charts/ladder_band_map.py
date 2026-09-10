"""Where each ladder band actually lands on the ground.

The per-window ladder table says a window has bands at, say, 2, 50 and 76 m
above the channel. It cannot say whether those bands are terrace treads - broad,
continuous, bounded by risers - or slices cut through a hillslope that happen to
fall between two breaks. Only the map says that.

READ IT AS PATCHES, NOT AS A CONTINUOUS SURFACE. A band is a terrace level if it
draws coherent benches with sharp edges, parallel to the valley and repeating
down it. Speckle scattered across a slope is the same band index meaning nothing
of the sort.

The band index is PER WINDOW. Two patches of the same colour in different
reaches were separated by different breaks and need not sit at the same height;
the along-channel windows are why a level can change colour downstream without
anything physical happening. Compare heights in the window table, not by colour.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

_BAND_COLOURS = ["#08519c", "#e6550d", "#6a51a3", "#31a354",
                 "#fdae6b", "#756bb1", "#78c679", "#cc4c02"]


def plot_ladder_band_map(
    hillshade: str | Path,
    band_raster: str | Path,
    out_png: str | Path,
    channel_fgb: str | Path | None = None,
    centreline_fgb: str | Path | None = None,
    max_band: int | None = None,
    title: str | None = None,
    dpi: int = 200,
) -> Path:
    """Map the ladder band index over a hillshade.

    Args:
        hillshade: Hillshade raster for the base.
        band_raster: `<stem>_band_index.tif`.
        out_png: Output PNG path.
        channel_fgb: Optional channel network lines.
        centreline_fgb: Optional valley centreline.
        max_band: Highest band to colour; above this the cells are drawn in the
            top colour. None uses the raster's own maximum.
        title: Figure title.
        dpi: Output resolution.

    Returns:
        The path written.
    """
    import rasterio
    from rasterio.plot import plotting_extent

    with rasterio.open(str(hillshade)) as src:
        hs = src.read(1).astype("float64")
        if src.nodata is not None:
            hs[hs == src.nodata] = np.nan
        extent = plotting_extent(src)

    with rasterio.open(str(band_raster)) as src:
        bi = src.read(1).astype("float64")
        if src.nodata is not None:
            bi[bi == src.nodata] = np.nan
        bi[bi < 0] = np.nan
        bext = plotting_extent(src)

    top = int(np.nanmax(bi)) if max_band is None else int(max_band)
    top = max(0, min(top, len(_BAND_COLOURS) - 1))
    cols = _BAND_COLOURS[:top + 1]
    cmap = ListedColormap(cols)
    norm = BoundaryNorm(np.arange(-0.5, top + 1.0, 1.0), cmap.N)

    fig, ax = plt.subplots(figsize=(13, 8))
    ax.imshow(hs, cmap="gray", extent=extent, origin="upper", zorder=0)
    im = ax.imshow(np.ma.masked_invalid(np.clip(bi, 0, top)), cmap=cmap,
                   norm=norm, extent=bext, origin="upper", alpha=0.85,
                   interpolation="nearest", zorder=1)

    if channel_fgb is not None and Path(channel_fgb).exists():
        import geopandas as gpd
        gpd.read_file(str(channel_fgb)).plot(ax=ax, color="#2171b5", lw=0.9,
                                             zorder=2)
    if centreline_fgb is not None and Path(centreline_fgb).exists():
        import geopandas as gpd
        gpd.read_file(str(centreline_fgb)).plot(ax=ax, color="black", lw=1.0,
                                                zorder=3)

    cb = fig.colorbar(im, ax=ax, ticks=np.arange(0, top + 1), shrink=0.7,
                      pad=0.02)
    cb.set_label("ladder band (0 = valley floor; per window)", fontsize=9)
    ax.set_xlabel("easting (m)")
    ax.set_ylabel("northing (m)")
    ax.set_aspect("equal")
    ax.set_title(title or "ladder bands over the hillshade")
    fig.tight_layout()
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png
