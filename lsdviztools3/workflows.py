"""Composite lsdtt3 -> lsdviztools3 workflows.

Convenience orchestration that drives lsdtt3 (via :mod:`lsdviztools3.driver`)
and renders the result (via :mod:`lsdviztools3.render`). Kept thin: each step is
a plain function call you could also run yourself.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional, Sequence, Union

import numpy as np
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.ops import unary_union

from .io.raster import load_raster, to_coordinates, to_gmt_grid
from .io.vector import load_vector, to_coordinates as to_vector_coordinates
from .driver import run_lsdtt
from .render.base import grid_to_display, save_figure
from .render.basins import render_basins
from .render.channels import render_channels
from .render.hillshade import render_hillshade
from .render.points import render_points
from .render.style import MapStyle

PathLike = Union[str, Path]


def _overlay_crs(raster: PathLike, coordinates):
    """The CRS a base raster is drawn in, for reprojecting its vector overlays.

    Resolving ``coordinates`` once from the raster (not per vector) keeps every
    overlay on the base map, including over an lsdtt3 oblique frame.
    """
    from .io.raster import resolve_target_crs

    return resolve_target_crs(load_raster(raster), coordinates)


def _valid_data_union(raster: PathLike, coordinates: str):
    """Return a polygon union covering finite cells of a raster."""
    grid = to_gmt_grid(to_coordinates(load_raster(raster), coordinates))
    valid = np.isfinite(grid.values)
    if not valid.any():
        raise ValueError(f"Raster has no finite cells: {raster}")

    transform = grid.rio.transform()
    polys = [
        shape(geom)
        for geom, value in shapes(valid.astype("uint8"), mask=valid, transform=transform)
        if value == 1
    ]
    return unary_union(polys)


def _clip_to_valid(gdf, valid_union):
    """Clip vector geometries to a valid-data polygon union."""
    out = gdf.copy()
    out.geometry = out.geometry.intersection(valid_union)
    out = out[(~out.geometry.isna()) & (~out.geometry.is_empty)]
    return out


def _ensure_complete_basins(
    basins_gdf,
    raster: PathLike,
    *,
    edge_buffer_pixels: int = 3,
    nodata_buffer_pixels: int = 1,
    basin_id_col: Optional[str] = None,
):
    """Raise if any basin polygon touches DEM edge/nodata safety buffers."""
    import rasterio
    from shapely.geometry import box

    with rasterio.open(raster) as src:
        cell = abs(float(src.transform.a))
        bounds = src.bounds
        raster_crs = src.crs

    work = basins_gdf
    if raster_crs is not None and work.crs is not None and str(work.crs) != str(raster_crs):
        work = work.to_crs(raster_crs)

    edge_buf = max(0, int(edge_buffer_pixels)) * cell
    nodata_buf = max(0, int(nodata_buffer_pixels)) * cell

    edge_inner = box(
        bounds.left + edge_buf,
        bounds.bottom + edge_buf,
        bounds.right - edge_buf,
        bounds.top - edge_buf,
    )
    valid_union = _valid_data_union(raster, coordinates="local")
    valid_inner = valid_union.buffer(-nodata_buf) if nodata_buf > 0 else valid_union

    complete_mask = (
        work.geometry.within(edge_inner)
        & work.geometry.within(valid_inner)
    )
    if bool(complete_mask.all()):
        return basins_gdf

    id_col = _infer_basin_id_column(work, basin_id_col)
    if id_col is not None:
        bad_ids = [str(v) for v in work.loc[~complete_mask, id_col].tolist()]
        msg = f"Incomplete basins detected (touch edge/nodata buffer): {', '.join(bad_ids)}"
    else:
        msg = f"Incomplete basins detected: {int((~complete_mask).sum())} polygon(s) touch edge/nodata buffer."
    raise ValueError(msg)


def _render_hillshade_base(base_raster: PathLike, hillshade_raster: Optional[PathLike], *, coordinates: str, style: MapStyle):
    """Render base from an lsdtt3 hillshade when provided, else from DEM."""
    from dataclasses import replace

    base_style = replace(style, cmap="gray")
    if hillshade_raster is not None:
        return render_hillshade(hillshade_raster, coordinates=coordinates, shade=False, style=base_style)
    return render_hillshade(base_raster, coordinates=coordinates, shade=True, style=base_style)


def _display_extent_for_raster(raster: PathLike, coordinates: str = "local"):
    da = to_coordinates(load_raster(raster), coordinates)
    crs = da.rio.crs
    _, extent = grid_to_display(to_gmt_grid(da), crs)
    return extent


def _add_header_lines(fig, lines, extent):
    """Draw multiline header text above the map frame in map coordinates."""
    xmin, xmax, ymin, ymax = extent
    width = xmax - xmin
    height = ymax - ymin
    x = xmin + 0.01 * width
    y = ymax + 0.12 * height
    dy = 0.05 * height
    for line in lines:
        if not str(line).strip():
            continue
        fig.text(
            x=x,
            y=y,
            text=str(line),
            font="14p,Helvetica-Bold,black",
            justify="TL",
            no_clip=True,
        )
        y -= dy


def _downsample_raster_for_display(da, factor: int):
    """Downsample raster data for faster map rendering.

    Uses block means over ``factor x factor`` cells in raster x/y dimensions.
    Intended for quick-look PNGs, not analysis.
    """
    if factor <= 1:
        return da
    x_dim = da.rio.x_dim
    y_dim = da.rio.y_dim
    out = da.coarsen({y_dim: factor, x_dim: factor}, boundary="trim").mean()
    out.attrs = da.attrs.copy()
    return out


def lsdtt3_hillshade_channels_basin(
    dem_path: PathLike,
    workdir: PathLike,
    *,
    stem: str = "dem",
    source_threshold: int = 150,
    lsdtt_dir: Optional[PathLike] = None,
) -> dict:
    """Run hillshade + channel + largest-basin extraction on a DEM with lsdtt3.

    Steps: ``lsdtt-raster-preprocessing`` (reproject to UTM, mask sea at
    ``min_elevation=0``, write hillshade) -> ``lsdtt-channel-extraction``
    (channel network lines) -> ``lsdtt-basin-extraction`` (single largest basin
    as a polygon). The channel/basin steps run on the sea-masked DEM.

    Args:
        dem_path: Input DEM (any GDAL raster; lsdtt3 reprojects to UTM if needed).
        workdir: Directory for lsdtt3 outputs.
        stem: Output filename stem (lsdtt3 ``write_fname``).
        source_threshold: Contributing-pixel threshold (lower => denser network).
        lsdtt_dir: Optional explicit lsdtt3 bin directory.

    Returns:
        Dict with keys ``hillshade``, ``channels``, ``basin`` (Paths).
    """
    dem = Path(dem_path)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    base = {"write_path": str(workdir), "write_fname": stem}

    run_lsdtt(
        "lsdtt-raster-preprocessing",
        dict(base, read_path=str(dem.parent), read_fname=dem.name,
             min_elevation=0.0, write_fill_raster="true", write_hillshade_raster="true"),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )
    processed_dem = f"{stem}_fill.tif"

    run_lsdtt(
        "lsdtt-channel-extraction",
        dict(base, read_path=str(workdir), read_fname=processed_dem,
             raster_is_filled="true",
             threshold_contributing_pixels=source_threshold,
             write_channel_network_lines="true"),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )

    run_lsdtt(
        "lsdtt-basin-extraction",
        dict(base, read_path=str(workdir), read_fname=processed_dem,
             raster_is_filled="true",
             threshold_contributing_pixels=source_threshold,
             select_all_basins="true", filter_all_but_largest_basin="true",
             write_basin_outlines="true", write_basin_outlines_as_polygon="true"),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )

    return {
        "hillshade": workdir / f"{stem}_hillshade.tif",
        "channels": workdir / f"{stem}_channel_network_lines.fgb",
        "basin": workdir / f"{stem}_junction_basin_polygons.fgb",
    }


def lsdtt3_stream_order_basins(
    dem_path: PathLike,
    workdir: PathLike,
    *,
    stem: str = "dem",
    source_threshold: int = 150,
    stream_order: int = 2,
    include_incomplete: bool = False,
    lsdtt_dir: Optional[PathLike] = None,
) -> dict:
    """Extract all junction basins of one Strahler order and write polygons.

    Steps: preprocessing (hillshade + sea mask), channel extraction, then
    ``lsdtt-basin-extraction`` with stream-order filtering.

    Args:
        dem_path: Input DEM.
        workdir: Directory for lsdtt3 outputs.
        stem: Output filename stem.
        source_threshold: Contributing-pixel threshold for channels.
        stream_order: Exact Strahler order to keep (e.g. 2 for 2nd-order basins).
        include_incomplete: Keep basins that touch nodata/edge. Defaults false.
        lsdtt_dir: Optional explicit lsdtt3 bin directory.

    Returns:
        Dict with keys ``hillshade``, ``channels``, ``basin`` and ``basin_order``.
    """
    dem = Path(dem_path)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    base = {"write_path": str(workdir), "write_fname": stem}

    run_lsdtt(
        "lsdtt-raster-preprocessing",
        dict(
            base,
            read_path=str(dem.parent),
            read_fname=dem.name,
            min_elevation=0.0,
            write_fill_raster="true",
            write_hillshade_raster="true",
        ),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )
    processed_dem = f"{stem}_fill.tif"

    run_lsdtt(
        "lsdtt-channel-extraction",
        dict(base, read_path=str(workdir), read_fname=processed_dem,
             raster_is_filled="true",
             threshold_contributing_pixels=source_threshold,
             write_channel_network_lines="true"),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )

    run_lsdtt(
        "lsdtt-basin-extraction",
        dict(
            base,
            read_path=str(workdir), read_fname=processed_dem,
            raster_is_filled="true",
            threshold_contributing_pixels=source_threshold,
            select_all_basins="true",
            filter_basins_by_stream_order="true",
            basins_min_stream_order=stream_order,
            basins_max_stream_order=stream_order,
            filter_incomplete_basins="false" if include_incomplete else "true",
            write_basin_outlines="true",
            write_basin_outlines_as_polygon="true",
        ),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )

    return {
        "hillshade": workdir / f"{stem}_hillshade.tif",
        "channels": workdir / f"{stem}_channel_network_lines.fgb",
        "basin": workdir / f"{stem}_junction_basin_polygons.fgb",
        "basin_order": stream_order,
    }


def _mn_suffix(m_over_n: float) -> str:
    """Filename suffix lsdtt3 uses for an m/n value (e.g. 0.45 -> ``0p450``).

    Mirrors lsdtt-chi-analysis' ``mn_suffix``: three decimals with ``.`` -> ``p``
    and a leading ``-`` -> ``m``.
    """
    s = f"{m_over_n:.3f}"
    return s.replace(".", "p").replace("-", "m")


def lsdtt3_chi_analysis(
    dem_path: PathLike,
    workdir: PathLike,
    *,
    stem: str = "dem",
    source_threshold: int = 300,
    n_largest_basins: int = 5,
    m_over_n: float = 0.45,
    n_pixels_move_outlet_upstream: int = 10,
    lsdtt_dir: Optional[PathLike] = None,
) -> dict:
    """Run chi analysis on a DEM and return the fitted-profile output paths.

    Steps: ``lsdtt-raster-preprocessing`` (reproject to UTM, mask sea at
    ``min_elevation=0``) -> ``lsdtt-chi-analysis`` with the N-largest-complete
    basin selection. Auto-detected coastal outlets are walked
    ``n_pixels_move_outlet_upstream`` pixels up the main stem so island basins
    clear the sea-nodata fringe and pass the completeness test; the ``N`` largest
    complete basins are then fitted at a single ``m_over_n``.

    Args:
        dem_path: Input DEM (any GDAL raster; reprojected to UTM if needed).
        workdir: Directory for lsdtt3 outputs.
        stem: Output filename stem (lsdtt3 ``write_fname``).
        source_threshold: Contributing-pixel threshold for channel heads.
        n_largest_basins: Keep the N largest complete basins (0 = all complete).
        m_over_n: Concavity (m/n) used for the chi fit.
        n_pixels_move_outlet_upstream: Upstream shift for interior outlets.
        lsdtt_dir: Optional explicit lsdtt3 bin directory.

    Returns:
        Dict of Paths: ``fit_nodes``, ``fit_segments``, ``ksn_summary`` (CSVs),
        ``data_maps`` and ``segment_lines`` (FlatGeobuf).
    """
    dem = Path(dem_path)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    base = {"write_path": str(workdir), "write_fname": stem}

    run_lsdtt(
        "lsdtt-raster-preprocessing",
        dict(base, read_path=str(dem.parent), read_fname=dem.name,
             min_elevation=0.0),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )
    sea_removed = f"{stem}_remove_seas.tif"

    run_lsdtt(
        "lsdtt-chi-analysis",
        dict(base, read_path=str(workdir), read_fname=sea_removed,
             threshold_contributing_pixels=source_threshold,
             n_largest_basins=n_largest_basins,
             n_pixels_move_outlet_upstream=n_pixels_move_outlet_upstream,
             write_chi_data_maps="true",
             chi_data_maps_write_vector_format="flatgeobuf",
             fit_chi_networks="true",
             write_fit_tables="true",
             write_chi_fit_segment_lines="true",
             write_vector_format="flatgeobuf",
             mn_min=m_over_n, mn_max=m_over_n, mn_step=0.02),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )

    mn = _mn_suffix(m_over_n)
    return {
        "fit_nodes": workdir / f"{stem}_chi_fit_nodes_mn_{mn}.csv",
        "fit_segments": workdir / f"{stem}_chi_fit_segments_mn_{mn}.csv",
        "ksn_summary": workdir / f"{stem}_chi_ksn_basin_summary_mn_{mn}.csv",
        "data_maps": workdir / f"{stem}_chi_data_maps.fgb",
        "segment_lines": workdir / f"{stem}_chi_fit_segment_lines_mn_{mn}.fgb",
    }


def render_chi_profiles(
    outputs: dict,
    out_dir: PathLike,
    *,
    stem: str = "chi",
    dpi: int = 300,
    stacked: bool = True,
    per_basin: bool = True,
    basins: Optional[Sequence[int]] = None,
    ksn_vmin: Optional[float] = None,
    ksn_vmax: Optional[float] = None,
    removed=None,
    show_removed: str = "ghost",
):
    """Plot fitted chi profiles from :func:`lsdtt3_chi_analysis` outputs.

    Writes a stacked chi plot (all basins on one axis, offset and coloured by
    log10 k_sn) and one chi-elevation profile per basin. Uses the matplotlib
    ``charts`` engine; requires the ``charts`` extra (matplotlib).

    Args:
        outputs: Dict from :func:`lsdtt3_chi_analysis` (needs ``fit_nodes`` and
            ``fit_segments``).
        out_dir: Directory to write PNGs into.
        stem: Output filename stem for the PNGs.
        dpi: Save resolution.
        stacked: Write the stacked chi plot.
        per_basin: Write one profile PNG per basin.
        basins: Optional basin list to render. Defaults to all in the fit table.
        ksn_vmin, ksn_vmax: Colour-scale limits in k_sn units (the scale is
            log10 k_sn), applied to every plot so per-basin profiles share one
            scale; ``None`` for each plot's own range.
        removed, show_removed: Tributaries to leave out (the ``removed`` frame
            of :func:`lsdviztools3.io.chi.select_tributaries`) and whether to
            ghost them (``"ghost"``) or hide them (``"hide"``).

    Returns:
        List of written PNG Paths.
    """
    from .io.chi import load_chi_fit
    from .charts import plot_chi_profile, plot_chi_stacked

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_chi_fit(outputs["fit_nodes"], outputs["fit_segments"])
    written = []

    selected_basins = sorted(nodes["basin_index"].unique()) if basins is None else [int(b) for b in basins]

    if stacked:
        fig, _ = plot_chi_stacked(nodes, basins=selected_basins, vmin=ksn_vmin, vmax=ksn_vmax,
                                  removed=removed, show_removed=show_removed)
        path = out_dir / f"{stem}_stacked.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written.append(path)

    if per_basin:
        for b in selected_basins:
            fig, _ = plot_chi_profile(nodes, basin=int(b), vmin=ksn_vmin, vmax=ksn_vmax,
                                      removed=removed, show_removed=show_removed)
            path = out_dir / f"{stem}_profile_basin{int(b)}.png"
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            written.append(path)

    return written


def render_chi_ksn_map(
    dem_path: PathLike,
    segment_lines: PathLike,
    png_path: PathLike,
    *,
    coordinates: str = "local",
    width: str = MapStyle.MAX_WIDTH,
    cmap: str = "viridis",
    high_res: bool = False,
    ksn_vmin: Optional[float] = None,
    ksn_vmax: Optional[float] = None,
    removed=None,
    ksn_classes=None,
    n_classes: int = 5,
    class_colors=None,
    class_color_mode: str = "auto",
):
    """Render the chi network coloured by log10 k_sn over a hillshade.

    Composes a grey hillshade of the (UTM, sea-masked) DEM, the fitted chi
    segment lines coloured through a CPT by log10 ``k_sn`` with a colourbar
    legend, and a scale bar. ``dem_path`` should be the projected DEM produced by
    preprocessing (e.g. ``<stem>_remove_seas.tif``); ``segment_lines`` is the
    ``*_chi_fit_segment_lines_mn_*.fgb`` from :func:`lsdtt3_chi_analysis`.

    Args:
        dem_path: Projected (UTM) DEM for the hillshade base.
        segment_lines: Fitted chi segment-lines fgb carrying ``k_sn``.
        png_path: Output PNG path.
        coordinates: Coordinate mode (``"local"``: projected km, UTM or an lsdtt3 oblique frame / ``"geographic"``).
        width: Figure width (defaults to the 17 cm max).
        cmap: CPT for the log10 k_sn colour scale.
        high_res: Save at high-res DPI.
        ksn_vmin, ksn_vmax: Colour-scale limits in k_sn units (the scale is
            log10 k_sn); ``None`` for the data's range.
        removed: Tributaries to leave off the map (the ``removed`` frame of
            :func:`lsdviztools3.io.chi.select_tributaries`); the colour range is
            then that of the kept ones. :func:`render_chi_removed_map` shows them.
        ksn_classes: Colour k_sn in classes instead of a continuous scale: the
            breaks in k_sn units (a list, tuple or array), or a scheme
            (``"equal_interval"``, ``"quantile"``, ``"natural_breaks"`` /
            ``"jenks"``) worked out in log10 k_sn. ``ksn_vmin`` / ``ksn_vmax``
            set the ends of the colour gradient. A class legend replaces the
            colourbar. See :mod:`lsdviztools3.colour_classes`.
        n_classes: Number of classes for a scheme (default 5).
        class_colors: One colour per class instead of sampling ``cmap``.
        class_color_mode: ``"gradient"``, ``"distinct"`` or ``"auto"``.

    Returns:
        ``Path`` to the written PNG.
    """
    from .io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from .render.base import grid_to_display, save_figure
    from .render.channels import render_channels
    from .render.decorations import add_scalebar
    from .render.hillshade import render_hillshade

    da = to_coordinates(load_raster(dem_path), coordinates)
    crs = get_crs(da)
    _, extent = grid_to_display(to_gmt_grid(da), crs)

    fig = render_hillshade(
        str(dem_path), coordinates=coordinates, shade=True,
        style=MapStyle(cmap="gray", figure_size=width, transparent=False),
    )
    lines = str(segment_lines)
    if removed is not None and len(removed):
        from .io.chi import removed_mask
        from .io.vector import load_vector
        gdf = load_vector(lines)
        lines = gdf[~removed_mask(gdf, removed)]
    fig = render_channels(
        lines, color_by="k_sn", log_color=True, cmap=cmap,
        colorbar=True, colorbar_label="log@-10@-k@-sn@-",
        coordinates=coordinates,
        style=MapStyle(line_width="1.2p", cmap=cmap, figure_size=width),
        fig=fig, vmin=ksn_vmin, vmax=ksn_vmax,
        classes=ksn_classes, n_classes=n_classes, class_colors=class_colors,
        class_color_mode=class_color_mode, legend_title="k@-sn@-",
    )
    add_scalebar(fig, extent, crs)
    save_figure(fig, png_path, style=MapStyle(figure_size=width, transparent=False),
                high_res=high_res)
    return Path(png_path)


#: Line style (GMT pen) and legend text per removal reason, for
#: :func:`render_chi_removed_map`. Dash patterns as well as colours, so the
#: reasons stay apart over the k_sn colours and in greyscale.
REMOVED_STYLES = {
    "ksn_outlier": ("1.4p,#d6266e,-", "removed: k_sn outlier"),
    "small_area": ("1.0p,black,-", "removed: small drainage area"),
    "small_fraction": ("1.0p,black,.", "removed: small fraction of basin"),
    "rank": ("1.0p,#5a5a5a,-.", "removed: outside the top tributaries"),
    "drains_to_removed": ("0.7p,#5a5a5a,.", "removed: drains into a removed one"),
}


def render_chi_removed_map(
    dem_path: PathLike,
    segment_lines: PathLike,
    removed,
    png_path: PathLike,
    *,
    coordinates: str = "local",
    width: str = MapStyle.MAX_WIDTH,
    cmap: str = "viridis",
    high_res: bool = False,
    ksn_vmin: Optional[float] = None,
    ksn_vmax: Optional[float] = None,
    legend_position: str = "JTR+jTR+o0.2c",
):
    """Map the tributaries removed from the chi plots, and why.

    The kept chi network is drawn coloured by log10 k_sn as in
    :func:`render_chi_ksn_map`. Each removed tributary is drawn in its reason's
    dashed style (:data:`REMOVED_STYLES`), with its junction (where it joins
    its parent) as an open circle, and a legend of the reasons present.

    Args:
        dem_path: Projected DEM for the hillshade.
        segment_lines: ``*_chi_fit_segment_lines_mn_*.fgb``.
        removed: The ``removed`` frame of
            :func:`lsdviztools3.io.chi.select_tributaries` (``basin_index``,
            ``branch_id``, ``reason``, ``junction_x``, ``junction_y``).
        png_path: Output PNG.
        coordinates, width, cmap, high_res, ksn_vmin, ksn_vmax: as
            :func:`render_chi_ksn_map`.
        legend_position: GMT legend position string.

    Returns:
        ``Path`` to the written PNG.
    """
    import io as _io

    import geopandas as gpd

    from .io.chi import removed_mask
    from .io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from .io.vector import load_vector
    from .render.base import grid_to_display, save_figure
    from .render.channels import render_channels
    from .render.decorations import add_scalebar
    from .render.hillshade import render_hillshade
    from .render.points import render_points

    da = to_coordinates(load_raster(dem_path), coordinates)
    crs = get_crs(da)
    _, extent = grid_to_display(to_gmt_grid(da), crs)
    gdf = load_vector(str(segment_lines))
    gone = removed_mask(gdf, removed)

    fig = render_hillshade(
        str(dem_path), coordinates=coordinates, shade=True,
        style=MapStyle(cmap="gray", figure_size=width, transparent=False),
    )
    fig = render_channels(
        gdf[~gone], color_by="k_sn", log_color=True, cmap=cmap,
        colorbar=True, colorbar_label="log@-10@-k@-sn@-", coordinates=coordinates,
        style=MapStyle(line_width="1.0p", cmap=cmap, figure_size=width),
        fig=fig, vmin=ksn_vmin, vmax=ksn_vmax,
    )
    rows = []
    if removed is not None and len(removed):
        keyed = gdf[gone].merge(removed[["basin_index", "branch_id", "reason"]],
                                on=["basin_index", "branch_id"], how="left")
        for reason, (pen, label) in REMOVED_STYLES.items():
            part = keyed[keyed["reason"] == reason]
            if part.empty:
                continue
            width_p, colour = pen.split(",")[0], ",".join(pen.split(",")[1:])
            fig = render_channels(part, coordinates=coordinates,
                                  style=MapStyle(line_width=width_p, line_color=colour,
                                                 figure_size=width), fig=fig)
            rows.append(f"S 0.3c - 0.6c - {pen} 0.9c {label} ({part[['basin_index', 'branch_id']].drop_duplicates().shape[0]})")
        j = removed.dropna(subset=["junction_x", "junction_y"])
        if len(j):
            pts = gpd.GeoDataFrame(j[["reason"]].copy(),
                                   geometry=gpd.points_from_xy(j["junction_x"], j["junction_y"]),
                                   crs=gdf.crs)
            fig = render_points(pts, coordinates=coordinates, symbol="c", pen="0.8p,black",
                                style=MapStyle(point_size="0.14c", line_color="white",
                                               figure_size=width), fig=fig)
            rows.append("S 0.3c c 0.14c white 0.8p,black 0.9c junction of a removed tributary")
    if rows:
        fig.legend(spec=_io.StringIO("\n".join(rows)), position=legend_position,
                   box="+gwhite@20+p0.5p,black")
    add_scalebar(fig, extent, crs)
    save_figure(fig, png_path, style=MapStyle(figure_size=width, transparent=False),
                high_res=high_res)
    return Path(png_path)


def render_chi_ksn_relief3d(
    dem_path: PathLike,
    segment_lines: PathLike,
    png_path: PathLike,
    *,
    coordinates: str = "local",
    azimuth: float = 135.0,
    elevation: float = 30.0,
    vertical_exaggeration: float = 2.0,
    cmap: str = "viridis",
    width: str = MapStyle.MAX_WIDTH,
    high_res: bool = False,
    ksn_vmin: Optional[float] = None,
    ksn_vmax: Optional[float] = None,
):
    """Drape the chi network (log10 k_sn) over a 3D perspective relief.

    Draws a grey hillshaded 3D surface of the DEM (``render_relief3d``) and
    overlays the fitted chi segments on the surface as 3D lines coloured by
    log10 ``k_sn`` (their z is sampled from the DEM), with a colourbar. A
    rayshader-style "draped network" view.

    Args:
        dem_path: Projected (UTM) DEM for the relief surface.
        segment_lines: ``*_chi_fit_segment_lines_mn_*.fgb`` carrying ``k_sn``.
        png_path: Output PNG path.
        coordinates: Coordinate mode (``"local"``: projected km, UTM or an lsdtt3 oblique frame recommended).
        azimuth: Viewing azimuth in degrees.
        elevation: Viewing elevation angle in degrees.
        vertical_exaggeration: True vertical:horizontal scale ratio.
        cmap: CPT for the log10 k_sn colour scale.
        width: Figure width.
        high_res: Save at high-res DPI.
        ksn_vmin, ksn_vmax: Colour-scale limits in k_sn units (the scale is
            log10 k_sn); ``None`` for the data's range.

    Returns:
        ``Path`` to the written PNG.
    """
    import numpy as np
    from lsdviztools3._gmt import require_pygmt
    pygmt = require_pygmt()
    import xarray as xr

    from .io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
    from .render.base import (
        auto_projection, gdf_to_display, grid_to_display, is_projected, save_figure,
    )
    from .render.channels import _iter_line_coords
    from .render.decorations import add_colourbar
    from .render.relief3d import _zsize_from_ve, render_relief3d

    # Shared 3D framing (must match between the surface and the draped lines).
    da = to_coordinates(load_raster(dem_path), coordinates)
    crs = get_crs(da)
    grid, extent = grid_to_display(to_gmt_grid(da), crs)
    zmin = float(np.nanmin(grid.values))
    zmax = float(np.nanmax(grid.values))
    region = [extent[0], extent[1], extent[2], extent[3], zmin, zmax]
    projection = auto_projection(extent, crs, width)
    perspective = [azimuth, elevation]
    if is_projected(crs):
        zsize = _zsize_from_ve(vertical_exaggeration, extent, zmax - zmin, width, crs)
    else:
        zsize = "4c"

    # Grey hillshaded relief surface.
    fig = render_relief3d(
        dem_path, coordinates=coordinates, azimuth=azimuth, elevation=elevation,
        zsize=zsize, shade=True, cmap="gray",
        style=MapStyle(cmap="gray", figure_size=width, transparent=False),
    )

    # Chi network in display (km) coords, with k_sn -> log10 colour scale.
    gdf = gdf_to_display(to_vector_coordinates(load_vector(segment_lines), crs), crs)
    ksn = gdf["k_sn"].to_numpy(dtype=float)
    z = np.full_like(ksn, np.nan)
    z[ksn > 0] = np.log10(ksn[ksn > 0])
    finite = z[np.isfinite(z)]
    if finite.size:
        from .colour_limits import colour_limits
        lo, hi, extend = colour_limits(z, ksn_vmin, ksn_vmax, log_scale=True)
        pygmt.makecpt(cmap=cmap, series=[lo, hi], background=True)

        for geom, zval in zip(gdf.geometry, z):
            if not np.isfinite(zval):
                continue
            for xs, ys in _iter_line_coords(geom):
                # Sample DEM elevation at each line vertex (pointwise interp),
                # lifting the segment slightly so it sits on the surface.
                elev = grid.interp(
                    x=xr.DataArray(xs, dims="pts"),
                    y=xr.DataArray(ys, dims="pts"),
                ).values
                good = np.isfinite(elev)
                if good.sum() < 2:
                    continue
                fig.plot3d(
                    x=np.asarray(xs)[good], y=np.asarray(ys)[good],
                    z=elev[good] + 0.01 * (zmax - zmin),
                    zvalue=zval, cmap=True, pen="1.5p+z",
                    region=region, projection=projection,
                    perspective=perspective, zsize=zsize,
                )
        add_colourbar(fig, label="log@-10@-k@-sn@-", extend=extend)

    save_figure(fig, png_path, style=MapStyle(figure_size=width, transparent=False),
                high_res=high_res)
    return Path(png_path)


def lsdtt3_swath_profile(
    dem_path: PathLike,
    workdir: PathLike,
    *,
    stem: str = "dem",
    source_threshold: int = 300,
    half_width_metres: float = 1000.0,
    bin_width_metres: float = 500.0,
    lsdtt_dir: Optional[PathLike] = None,
) -> dict:
    """Run a swath profile along the largest basin's main channel.

    Steps: ``lsdtt-raster-preprocessing`` (UTM + sea masked) ->
    ``lsdtt-chi-analysis`` (largest complete basin only, writing the chi data-map
    channel lines) -> pick the main stem (longest line) as the baseline ->
    ``lsdtt-swath-profile`` swathing the DEM elevation along that baseline.

    Args:
        dem_path: Input DEM (any GDAL raster; reprojected to UTM if needed).
        workdir: Directory for lsdtt3 outputs.
        stem: Output filename stem.
        source_threshold: Contributing-pixel threshold for channel heads.
        half_width_metres: Swath half-width (0 = whole valid extent).
        bin_width_metres: Along-axis bin width.
        lsdtt_dir: Optional explicit lsdtt3 bin directory.

    Returns:
        Dict of Paths: ``swath_csv``, ``baseline``, ``along_axis``,
        ``perp_distance``, ``signed_perp``.
    """
    from .io.swath import extract_main_channel_baseline

    dem = Path(dem_path)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    base = {"write_path": str(workdir), "write_fname": stem}

    run_lsdtt(
        "lsdtt-raster-preprocessing",
        dict(base, read_path=str(dem.parent), read_fname=dem.name,
             min_elevation=0.0),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )
    sea_removed = f"{stem}_remove_seas.tif"

    # Largest complete basin's channel network as primary-path lines.
    run_lsdtt(
        "lsdtt-chi-analysis",
        dict(base, read_path=str(workdir), read_fname=sea_removed,
             threshold_contributing_pixels=source_threshold,
             n_largest_basins=1,
             write_chi_data_maps="true",
             chi_data_maps_write_vector_format="flatgeobuf",
             write_chi_data_map_lines="true",
             write_vector_format="flatgeobuf",
             fit_chi_networks="false"),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )

    baseline = extract_main_channel_baseline(
        workdir / f"{stem}_chi_data_map_lines.fgb",
        workdir / f"{stem}_baseline.fgb",
        basin=0,
    )

    run_lsdtt(
        "lsdtt-swath-profile",
        dict(base, read_path=str(workdir), read_fname=sea_removed,
             baseline_vector_fname=str(baseline),
             half_width_metres=half_width_metres,
             bin_width_metres=bin_width_metres,
             write_distance_rasters="true"),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )

    return {
        "swath_csv": workdir / f"{stem}_swath_profile.csv",
        "baseline": Path(baseline),
        "along_axis": workdir / f"{stem}_swath_along_axis.tif",
        "perp_distance": workdir / f"{stem}_swath_perp_distance.tif",
        "signed_perp": workdir / f"{stem}_swath_signed_perp.tif",
    }


def render_swath_profile(
    outputs: dict,
    png_path: PathLike,
    *,
    value_label: str = "Elevation (m)",
    dpi: int = 300,
):
    """Plot the swath envelope from :func:`lsdtt3_swath_profile` outputs.

    Writes a distance-vs-elevation swath plot (min/max envelope, IQR band,
    median + mean). Uses the matplotlib ``charts`` engine.

    Args:
        outputs: Dict from :func:`lsdtt3_swath_profile` (needs ``swath_csv``).
        png_path: Output PNG path.
        value_label: Y-axis label for the swathed quantity.
        dpi: Save resolution.

    Returns:
        ``Path`` to the written PNG.
    """
    from .io.swath import load_swath
    from .charts import plot_swath

    df = load_swath(outputs["swath_csv"])
    fig, _ = plot_swath(df, value_label=value_label)
    png_path = Path(png_path)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    return png_path


def render_hand_profile(
    hand_profile_csv: PathLike,
    png_path: PathLike,
    *,
    band_breaks_csv: Optional[PathLike] = None,
    dpi: int = 300,
    **kwargs,
):
    """Plot HAND vs distance along the channel for classified pixels.

    Reads the per-pixel ``*_hand_profile.csv`` written by
    ``lsdtt-flat-surface-finder`` (``write_hand_profile_csv``) and draws the
    floodplain / terrace HAND populations along the channel; optionally
    overlays the reach-relative band break levels from
    ``*_hand_band_breaks.csv``. Uses the matplotlib ``charts`` engine.

    Args:
        hand_profile_csv: Per-pixel HAND profile table.
        png_path: Output PNG path.
        band_breaks_csv: Optional band-breaks table to overlay.
        dpi: Save resolution.
        **kwargs: Passed through to
            :func:`lsdviztools3.charts.plot_hand_distance` (e.g. ``source_id``,
            ``log_y``, ``max_points``).

    Returns:
        ``Path`` to the written PNG.
    """
    from .io.table import load_table
    from .charts import plot_hand_distance

    df = load_table(hand_profile_csv)
    breaks = load_table(band_breaks_csv) if band_breaks_csv else None
    fig, _ = plot_hand_distance(df, breaks=breaks, **kwargs)
    png_path = Path(png_path)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    return png_path


def render_hand_window_populations(
    window_hist_csv: PathLike,
    png_path: PathLike,
    *,
    breaks_csv: Optional[PathLike] = None,
    source_id: int = 0,
    dpi: int = 300,
    **kwargs,
):
    """Plot per-window HAND populations with the selected valley breaks.

    Reads the ``*_valley_hand_window_hist.csv`` (and optionally the
    ``*_valley_breaks.csv``) written by ``lsdtt-valley-centreline`` and draws
    one panel per along-channel window, so the reach-adaptive threshold
    selection is inspectable. Uses the matplotlib ``charts`` engine.

    Args:
        window_hist_csv: Per-window HAND histogram table.
        png_path: Output PNG path.
        breaks_csv: Optional per-window break table to overlay.
        source_id: Channel to plot.
        dpi: Save resolution.
        **kwargs: Passed through to
            :func:`lsdviztools3.charts.plot_hand_window_populations`.

    Returns:
        ``Path`` to the written PNG.
    """
    from .io.table import load_table
    from .charts import plot_hand_window_populations

    hist_df = load_table(window_hist_csv)
    breaks_df = load_table(breaks_csv) if breaks_csv else None
    fig, _ = plot_hand_window_populations(
        hist_df, breaks_df=breaks_df, source_id=source_id, **kwargs)
    png_path = Path(png_path)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    return png_path


def render_hillshade_channels_basin(
    outputs: dict,
    png_path: PathLike,
    *,
    coordinates: str = "local",
    width: str = MapStyle.SINGLE_COLUMN,
    high_res: bool = False,
):
    """Render hillshade base + channel overlay + semi-transparent largest basin.

    Saves at publication DPI (300, or 600 if ``high_res``) on an opaque white
    background to avoid transparency artifacts in downstream viewers/printers.
    ``width`` defaults to single-column (8.5 cm).
    """
    valid_union = _valid_data_union(outputs["hillshade"], coordinates)

    overlay_crs = _overlay_crs(outputs["hillshade"], coordinates)
    channels = _clip_to_valid(
        to_vector_coordinates(load_vector(outputs["channels"]), overlay_crs),
        valid_union,
    )
    basins = _clip_to_valid(
        to_vector_coordinates(load_vector(outputs["basin"]), overlay_crs),
        valid_union,
    )
    if not basins.empty:
        basins = basins.copy()
        basins.geometry = basins.geometry.buffer(0)
        basins = basins[(~basins.geometry.isna()) & (~basins.geometry.is_empty)]

    fig = render_hillshade(
        str(outputs["hillshade"]),
        coordinates=coordinates,
        shade=False,
        style=MapStyle(cmap="gray", figure_size=width, transparent=False),
    )

    if not channels.empty:
        fig = render_channels(
            channels,
            coordinates=coordinates,
            style=MapStyle(line_width="0.8p", line_color="dodgerblue4", figure_size=width),
            fig=fig,
        )

    if not basins.empty:
        fig = render_basins(
            basins,
            coordinates=coordinates,
            fill="blue",
            transparency=70,
            style=MapStyle(line_width="1.2p", line_color="navy", figure_size=width),
            fig=fig,
        )
    save_figure(
        fig,
        png_path,
        style=MapStyle(figure_size=width, transparent=False),
        high_res=high_res,
    )
    return fig


def render_hillshade_stream_order_basins(
    outputs: dict,
    png_path: PathLike,
    *,
    coordinates: str = "local",
    width: str = MapStyle.SINGLE_COLUMN,
    high_res: bool = False,
):
    """Render hillshade + channels + all selected-order basin polygons.

    Works with outputs from :func:`lsdtt3_stream_order_basins`.
    """
    order = outputs.get("basin_order", "selected")
    valid_union = _valid_data_union(outputs["hillshade"], coordinates)

    overlay_crs = _overlay_crs(outputs["hillshade"], coordinates)
    channels = _clip_to_valid(
        to_vector_coordinates(load_vector(outputs["channels"]), overlay_crs),
        valid_union,
    )
    basins = _clip_to_valid(
        to_vector_coordinates(load_vector(outputs["basin"]), overlay_crs),
        valid_union,
    )
    if not basins.empty:
        basins = basins.copy()
        basins.geometry = basins.geometry.buffer(0)
        basins = basins[(~basins.geometry.isna()) & (~basins.geometry.is_empty)]

    fig = render_hillshade(
        str(outputs["hillshade"]),
        coordinates=coordinates,
        shade=False,
        style=MapStyle(cmap="gray", figure_size=width, transparent=False),
    )

    if not channels.empty:
        fig = render_channels(
            channels,
            coordinates=coordinates,
            style=MapStyle(line_width="0.6p", line_color="deepskyblue4", figure_size=width),
            fig=fig,
        )

    if not basins.empty:
        fig = render_basins(
            basins,
            coordinates=coordinates,
            fill="orange",
            transparency=60,
            style=MapStyle(line_width="0.8p", line_color="darkred", figure_size=width),
            fig=fig,
        )

    from .render.decorations import add_colourbar
    add_colourbar(fig, label=f"Hillshade; order {order} basins", placement="right")
    save_figure(
        fig,
        png_path,
        style=MapStyle(figure_size=width, transparent=False),
        high_res=high_res,
    )
    return fig


def _infer_layer_kind(gdf) -> str:
    """Guess a layer's geometry kind ("polygon"/"line"/"point") from a gdf."""
    types = set(gdf.geom_type.dropna().unique())
    if any("Polygon" in t for t in types):
        return "polygon"
    if any("Line" in t for t in types):
        return "line"
    return "point"


def hillshade_with_vectors(
    base_raster: PathLike,
    layers: list,
    savepath: Optional[PathLike] = None,
    *,
    coordinates: str = "local",
    shade: bool = True,
    base_cmap: Optional[str] = None,
    width: str = MapStyle.DOUBLE_COLUMN,
    region=None,
    clip_to_valid: bool = False,
    legend: bool = True,
    legend_position: str = "JTR+jTR+o0.2c",
    high_res: bool = False,
    display_downsample: int = 1,
    fig=None,
):
    """Render a hillshade and overlay any number of vector layers.

    A general, dataset-agnostic map builder: give it a base raster and a list of
    vector layers (lines, points, and/or polygons, in any mix) and it draws each
    in its own style with an optional legend. Reusable for any dataset — channel
    networks, faults, basins, sample sites, roads, etc.

    Args:
        base_raster: path to a DEM (``shade=True``; shaded relief) or a
            pre-computed hillshade (``shade=False``).
        layers: list of layer dicts. Recognised keys::

                path        vector file path (or a GeoDataFrame)         [required]
                layer       sub-layer name for multi-layer files (gpkg)  [optional]
                kind        "auto" (default) | "line" | "point" | "polygon"
                label       legend label (omit to skip in legend)
                color       line/outline/point colour          (default "black")
                width       line / outline width               (default "1p")
                size        point symbol size                  (default "0.15c")
                symbol      GMT point symbol code (c,s,t,...)   (default "c")
                fill        polygon/point fill colour          (default = color)
                transparency  polygon fill transparency %      (default 60)

            Layers are drawn in list order (first = bottom).
        savepath: write a PNG here (``None`` = don't save, just return the fig).
        coordinates: ``"local"`` (default: the base raster's own CRS when
            projected, UTM only for a geographic raster) / ``"utm"`` /
            ``"geographic"`` / explicit CRS. Base and all layers are put in
            this CRS so they align.
        shade: treat ``base_raster`` as elevation (True) or a ready hillshade.
        base_cmap: colour map for the base (default ``"gray"`` for a clean
            backdrop; pass e.g. ``"geo"`` for a coloured elevation base).
        region: ``[xmin, xmax, ymin, ymax]`` in the raster's NATIVE CRS to zoom;
            crops the base raster and everything follows.
        clip_to_valid: clip overlays to the base raster's finite-data footprint
            (useful when a vector tile is larger than the raster).
        legend: draw a legend from the layers that have a ``label``.
        legend_position: GMT legend position string.
        display_downsample: integer factor for plot-time raster downsampling
            before rendering (``1`` = full resolution, ``2`` = half resolution
            in x/y, etc.). Use this to speed up large-area QC maps.
        fig: existing ``pygmt.Figure`` to compose onto (``None`` = new).

    Returns:
        The ``pygmt.Figure``.
    """
    if display_downsample < 1:
        raise ValueError("display_downsample must be >= 1")

    # Optional prep: zoom and/or display downsampling. We write a temp raster so
    # hillshade render and valid-footprint clipping use the same grid.
    if region is not None or display_downsample > 1:
        import tempfile

        da = load_raster(str(base_raster))
        if region is not None:
            xmin, xmax, ymin, ymax = region
            da = da.rio.clip_box(xmin, ymin, xmax, ymax)
        if display_downsample > 1:
            da = _downsample_raster_for_display(da, int(display_downsample))

        tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
        da.rio.to_raster(tmp.name)
        base_raster = tmp.name

    valid_union = (
        _valid_data_union(base_raster, coordinates) if clip_to_valid else None
    )

    fig = render_hillshade(
        str(base_raster),
        coordinates=coordinates,
        shade=shade,
        style=MapStyle(cmap=base_cmap or "gray", figure_size=width,
                       transparent=False),
        fig=fig,
    )

    legend_rows = []
    overlay_crs = _overlay_crs(base_raster, coordinates)
    for layer in layers:
        src = layer["path"]
        gdf = load_vector(str(src), layer.get("layer")) if not hasattr(
            src, "geometry") else src
        gdf = to_vector_coordinates(gdf, overlay_crs)
        if valid_union is not None:
            gdf = _clip_to_valid(gdf, valid_union)
        if gdf.empty:
            continue

        kind = layer.get("kind", "auto")
        if kind == "auto":
            kind = _infer_layer_kind(gdf)

        color = layer.get("color", "black")
        lw = layer.get("width", "1p")
        size = layer.get("size", "0.15c")
        symbol = layer.get("symbol", "c")
        fill = layer.get("fill", color)
        transp = layer.get("transparency", 60)

        if kind == "polygon":
            fig = render_basins(
                gdf, coordinates=coordinates, fill=fill, transparency=transp,
                style=MapStyle(line_width=lw, line_color=color,
                               figure_size=width), fig=fig)
            legend_rows.append(
                (layer.get("label"), f"S 0.3c s 0.35c {fill} {lw},{color}"))
        elif kind == "point":
            # line_color carries the symbol fill; pen (color) draws the outline
            # and is what outline-only symbols (x, +, ...) render with.
            fig = render_points(
                gdf, coordinates=coordinates, symbol=symbol,
                pen=f"{lw},{color}",
                style=MapStyle(point_size=size, line_color=fill, cmap=base_cmap
                               or "gray", figure_size=width), fig=fig)
            legend_rows.append(
                (layer.get("label"), f"S 0.3c {symbol} {size} {fill} {lw},{color}"))
        else:  # line
            fig = render_channels(
                gdf, coordinates=coordinates,
                style=MapStyle(line_width=lw, line_color=color,
                               figure_size=width), fig=fig)
            legend_rows.append((layer.get("label"), f"S 0.3c - 0.6c - {lw},{color}"))

    rows = [f"{sym} 0.9c {lbl}" for lbl, sym in legend_rows if lbl]
    if legend and rows:
        import io as _io
        fig.legend(spec=_io.StringIO("\n".join(rows)), position=legend_position,
                   box="+gwhite@20+p0.5p,black")

    if savepath is not None:
        save_figure(fig, savepath,
                    style=MapStyle(figure_size=width, transparent=False),
                    high_res=high_res)
    return fig


def compare_stream_networks(
    base_raster: PathLike,
    networks: list,
    png_path: PathLike,
    *,
    clip_to_valid: bool = True,
    **kwargs,
):
    """Stream-burning preset over :func:`hillshade_with_vectors`.

    Thin wrapper that overlays a set of stream networks on a hillshade to see
    where the flow-routed network departs from the conditioning vectors. Each
    entry of ``networks`` is a layer dict (see ``hillshade_with_vectors``);
    ``clip_to_valid`` defaults True because vector tiles usually exceed the DEM.
    """
    return hillshade_with_vectors(
        base_raster, networks, png_path,
        clip_to_valid=clip_to_valid, **kwargs)


def plot_segment_classification(
    base_raster: PathLike,
    segments,
    png_path: Optional[PathLike] = None,
    *,
    status_field: str = "status",
    followed_value=0,
    burned_gaps=None,
    flow_traces=None,
    extra_layers=None,
    followed_color="dodgerblue4",
    unresolved_color="red",
    gap_color="orange",
    trace_color="green",
    clip_to_valid: bool = True,
    **kwargs,
):
    """Map the trace-divergence burn segment-by-segment over a hillshade.

    The single most useful QC plot for ``lsdtt-stream-burning``'s
    ``snap_to_accumulation`` mode. It colours each fixed-length vector segment
    by how the burn treated it, so you can see at a glance where the LiDAR
    already carried the channel and where the burn could not fix the flow:

    - **followed** (``status == followed_value``): the D8 trace stayed within
      tolerance of the vector; nothing was burned.
    - **unresolved** (otherwise): the trace still diverged after all passes —
      the cases to investigate.

    Args:
        base_raster: DEM (shaded relief) or hillshade for the backdrop.
        segments: path to (or GeoDataFrame of) the ``_segments.fgb`` written by
            ``write_burn_diagnostics``; must carry ``status_field``.
        png_path: output PNG (``None`` = return the figure only).
        status_field: attribute holding the per-segment status code.
        followed_value: status value that means "followed".
        burned_gaps: optional ``_burned_gaps.fgb`` (the carved trenches) to draw.
        flow_traces: optional ``_flow_traces.fgb`` (sampled D8 traces) to draw.
        extra_layers: optional list of additional ``hillshade_with_vectors``
            layer dicts drawn on top (e.g. the source vector network).
        followed_color, unresolved_color, gap_color, trace_color: line colours.
        clip_to_valid: clip overlays to the raster footprint (default True).
        **kwargs: forwarded to :func:`hillshade_with_vectors` (e.g. ``region=``).

    Returns:
        The ``pygmt.Figure``.
    """
    gdf = load_vector(str(segments)) if not hasattr(segments, "geometry") \
        else segments
    if status_field not in gdf.columns:
        raise KeyError(
            f"segments layer has no '{status_field}' field; columns: "
            f"{list(gdf.columns)}")

    followed = gdf[gdf[status_field] == followed_value]
    unresolved = gdf[gdf[status_field] != followed_value]

    layers = []
    if flow_traces is not None:
        layers.append({"path": flow_traces, "kind": "line",
                       "color": trace_color, "width": "0.6p",
                       "label": "D8 flow trace"})
    if not followed.empty:
        layers.append({"path": followed, "kind": "line",
                       "color": followed_color, "width": "1.2p",
                       "label": "followed (LiDAR carries it)"})
    if burned_gaps is not None:
        layers.append({"path": burned_gaps, "kind": "line",
                       "color": gap_color, "width": "2p",
                       "label": "burned trench"})
    if not unresolved.empty:
        layers.append({"path": unresolved, "kind": "line",
                       "color": unresolved_color, "width": "2p",
                       "label": "unresolved (still diverging)"})
    if extra_layers:
        layers.extend(extra_layers)

    return hillshade_with_vectors(
        base_raster, layers, png_path,
        clip_to_valid=clip_to_valid, **kwargs)


def plot_algoviz(
    base_raster: PathLike,
    d8_traces,
    source_path,
    png_path: Optional[PathLike] = None,
    *,
    burns=None,
    pass_filter: Optional[int] = 0,
    outcome_field: str = "outcome",
    zoom_pad_m: float = 150.0,
    label_iterations: bool = True,
    clip_to_valid: bool = False,
    **kwargs,
):
    """Replay the snap_to_accumulation source-path burn for one channel.

    Visualises the per-iteration behaviour of ``burn_vectormap_source_paths``
    (the ``_algoviz_jXX_*`` layers): the vector source path, the D8 flow traces
    used to test divergence (coloured by outcome), and the segments that were
    burned. This is the diagnostic for "why does this channel fail downstream":
    you can see segment by segment where the D8 follows the vector, where it
    strays (and to what), and where a burn was (or was not) applied.

    Args:
        base_raster: DEM (shaded relief) or hillshade for the backdrop.
        d8_traces: path to (or GeoDataFrame of) ``_algoviz_jXX_d8.fgb``. Must
            carry ``pass``, ``iteration`` and ``outcome`` fields.
        source_path: path to ``_algoviz_jXX_source.fgb`` (the vector source line).
        png_path: output PNG (``None`` = return the figure only).
        burns: optional ``_algoviz_jXX_burns.fgb`` (carved segments).
        pass_filter: keep only this re-flow ``pass`` (``None`` = all passes). The
            burn re-traces every pass, so a single pass is usually clearest.
        outcome_field: trace attribute; 0 followed, 1 burned, 2 skip-conflict.
        zoom_pad_m: pad (m) around the source path bbox for the map region
            (``region=`` in kwargs overrides this).
        label_iterations: annotate each D8 trace start with its iteration index.
        clip_to_valid: clip overlays to the raster footprint (default False).
        **kwargs: forwarded to :func:`hillshade_with_vectors`.

    Returns:
        The ``pygmt.Figure``.
    """
    traces = load_vector(str(d8_traces)) if not hasattr(d8_traces, "geometry") \
        else d8_traces
    src = load_vector(str(source_path)) if not hasattr(source_path, "geometry") \
        else source_path
    if pass_filter is not None and "pass" in traces.columns:
        traces = traces[traces["pass"] == pass_filter]

    # Auto-zoom to the source path unless the caller fixed a region.
    if "region" not in kwargs or kwargs["region"] is None:
        minx, miny, maxx, maxy = src.total_bounds
        kwargs["region"] = [minx - zoom_pad_m, maxx + zoom_pad_m,
                            miny - zoom_pad_m, maxy + zoom_pad_m]

    foll = traces[traces[outcome_field] == 0]
    burned = traces[traces[outcome_field] == 1]
    conflict = traces[traces[outcome_field] == 2]

    layers = [
        {"path": src, "kind": "line", "color": "cyan", "width": "3p",
         "label": "VML source path"},
    ]
    if not foll.empty:
        layers.append({"path": foll, "kind": "line", "color": "green",
                       "width": "1p", "label": "D8 trace: follows"})
    if not conflict.empty:
        layers.append({"path": conflict, "kind": "line", "color": "red",
                       "width": "1.2p", "label": "D8 trace: diverged, skipped"})
    if not burned.empty:
        layers.append({"path": burned, "kind": "line", "color": "orange",
                       "width": "1.2p", "label": "D8 trace: diverged -> burned"})
    if burns is not None:
        bg = load_vector(str(burns)) if not hasattr(burns, "geometry") else burns
        if not bg.empty:
            layers.append({"path": bg, "kind": "line", "color": "magenta",
                           "width": "3p", "label": "burned segment"})

    fig = hillshade_with_vectors(
        base_raster, layers, None, clip_to_valid=clip_to_valid, **kwargs)

    # Annotate each trace's start point with its iteration number.
    if label_iterations and "iteration" in traces.columns and not traces.empty:
        coords = kwargs.get("coordinates", "local")
        tr = to_vector_coordinates(traces, _overlay_crs(base_raster, coords))
        for _, row in tr.iterrows():
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue
            pt = geom.coords[0]
            fig.text(x=pt[0], y=pt[1], text=str(int(row["iteration"])),
                     font="6p,Helvetica-Bold,black", fill="white@30",
                     offset="0.1c/0.1c")

    if png_path is not None:
        save_figure(fig, png_path,
                    style=MapStyle(figure_size=MapStyle.DOUBLE_COLUMN,
                                   transparent=False),
                    high_res=kwargs.get("high_res", True))
    return fig


def sample_dem_along_line(dem_path, line, step=None):
    """Sample a DEM's elevation at regular spacing along a polyline.

    Args:
        dem_path: path to a single-band raster (DEM).
        line: a shapely LineString (in the DEM's CRS).
        step: spacing in map units (default = the DEM cell size).

    Returns:
        ``(dist, elev)`` numpy arrays: cumulative distance (m) and elevation,
        with DEM nodata converted to ``np.nan``.
    """
    import numpy as np
    import rasterio
    with rasterio.open(str(dem_path)) as ds:
        if step is None:
            step = abs(ds.transform.a)
        nod = ds.nodata
        L = line.length
        dist = np.arange(0.0, L + step, step)
        pts = [(line.interpolate(d).x, line.interpolate(d).y) for d in dist]
        elev = np.array([v[0] for v in ds.sample(pts)], dtype=float)
    if nod is not None:
        elev[elev == nod] = np.nan
    return dist, elev


def trace_d8_flowpath(dem_path, start_xy, max_steps=200000):
    """Trace a D8 steepest-descent path from a start point over a DEM.

    A lightweight follower for QC long profiles: from the start cell it steps to
    the lowest of the 8 neighbours while one is strictly lower, accumulating map
    distance (diagonals * sqrt(2)). Stops at a pit/flat (no lower neighbour), the
    grid edge, or nodata. Best run on a hydrologically conditioned (burned/filled)
    DEM, where the descent is well defined.

    Args:
        dem_path: path to a single-band DEM.
        start_xy: ``(x, y)`` start point in the DEM CRS.
        max_steps: safety cap on the number of cells walked.

    Returns:
        ``(line, dist, elev)``: a shapely LineString of the path and the
        cumulative-distance / elevation arrays along it.
    """
    import numpy as np
    import rasterio
    from shapely.geometry import LineString
    with rasterio.open(str(dem_path)) as ds:
        Z = ds.read(1).astype(float)
        nod = ds.nodata
        if nod is not None:
            Z[Z == nod] = np.nan
        tr = ds.transform
        cs_x, cs_y = abs(tr.a), abs(tr.e)
        nrow, ncol = Z.shape
        r, c = ds.index(start_xy[0], start_xy[1])
        if not (0 <= r < nrow and 0 <= c < ncol):
            raise ValueError("start point outside the DEM")
        xs, ys, zs = [], [], []
        seen = set()
        nbrs = [(-1, -1), (-1, 0), (-1, 1), (0, -1),
                (0, 1), (1, -1), (1, 0), (1, 1)]
        for _ in range(max_steps):
            if (r, c) in seen:
                break
            seen.add((r, c))
            x, y = ds.xy(r, c)
            xs.append(x); ys.append(y); zs.append(Z[r, c])
            # Prefer the steepest strictly-lower neighbour; if none (a flat or a
            # gentle burned trench), step to the lowest neighbour that is <=
            # current and unvisited so the trace crosses flats instead of stopping.
            best = None; best_z = Z[r, c]
            flat = None; flat_z = np.inf
            for dr, dc in nbrs:
                rr, cc = r + dr, c + dc
                if not (0 <= rr < nrow and 0 <= cc < ncol):
                    continue
                z = Z[rr, cc]
                if np.isnan(z):
                    continue
                if z < best_z:
                    best_z = z; best = (rr, cc)
                if z <= Z[r, c] + 1e-6 and (rr, cc) not in seen and z < flat_z:
                    flat_z = z; flat = (rr, cc)
            nxt = best if best is not None else flat
            if nxt is None:
                break
            r, c = nxt
    xs = np.array(xs); ys = np.array(ys); zs = np.array(zs)
    if len(xs) < 2:
        return None, np.array([0.0]), zs
    seg = np.hypot(np.diff(xs), np.diff(ys))
    dist = np.concatenate([[0.0], np.cumsum(seg)])
    return LineString(np.column_stack([xs, ys])), dist, zs


def plot_long_profile(series, png_path=None, *, title=None,
                      xlabel="Distance downstream (m)",
                      ylabel="Elevation (m)", figsize=(11, 5),
                      annotations=None, dpi=200):
    """Plot one or more elevation-vs-distance long profiles on shared axes.

    Args:
        series: list of dicts, each ``{"dist": arr, "elev": arr, "label": str,
            "color": str?, "ls": str?, "lw": float?}``.
        png_path: output PNG (``None`` = just return the matplotlib fig/ax).
        title, xlabel, ylabel, figsize, dpi: usual plot controls.
        annotations: optional list of ``{"x": float, "text": str, "color": str?}``
            drawn as labelled vertical lines (e.g. a road crossing location).

    Returns:
        ``(fig, ax)``.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=figsize)
    for s in series:
        ax.plot(s["dist"], s["elev"], label=s.get("label"),
                color=s.get("color"), linestyle=s.get("ls", "-"),
                lw=s.get("lw", 1.6))
    if annotations:
        for a in annotations:
            ax.axvline(a["x"], color=a.get("color", "0.4"), ls=":", lw=1.0)
            ax.text(a["x"], ax.get_ylim()[1], a["text"], rotation=90,
                    va="top", ha="right", fontsize=8, color=a.get("color", "0.4"))
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if png_path is not None:
        fig.savefig(str(png_path), dpi=dpi)
        plt.close(fig)
    return fig, ax


def render_cosmo_comparison_figures(
    scan_comparison_csv: PathLike,
    out_dir: PathLike,
    *,
    sample_comparison_csv: Optional[PathLike] = None,
    reference: str = "cronus_v3",
    representative_pressure_hpa: Optional[float] = None,
    representative_erate_gcm2yr: Optional[float] = None,
    stem: str = "cosmo_comparison",
    dpi: int = 220,
) -> dict:
    """Render standard cosmo comparison figures from lsdtt3 workflow CSVs.

    Args:
        scan_comparison_csv: Path to ``single_point_scan_comparison.csv``.
        out_dir: Output directory for figure PNGs.
        sample_comparison_csv: Optional path to a real-sample comparison table
            (for example OCTOPUS/riversand summary-style CSVs).
        reference: ``cronus_v3`` or ``cronuscalc`` for sweep residual panels.
        representative_pressure_hpa: Optional nearest-pressure slice for the
            sweep panel figure.
        representative_erate_gcm2yr: Optional nearest-erate slice for the sweep
            panel figure.
        stem: Filename stem for outputs.
        dpi: PNG output dpi.

    Returns:
        Dict with written figure paths. Always includes ``sweep_panels`` and,
        when ``sample_comparison_csv`` is provided, ``sample_residuals``.
    """
    from .charts import (
        load_cosmo_table,
        plot_cosmo_sample_residuals,
        plot_cosmo_sweep_panels,
    )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scan_df = load_cosmo_table(scan_comparison_csv)
    fig_sweep, _ = plot_cosmo_sweep_panels(
        scan_df,
        reference=reference,
        representative_pressure_hpa=representative_pressure_hpa,
        representative_erate_gcm2yr=representative_erate_gcm2yr,
    )
    sweep_png = out_dir / f"{stem}_sweep_{reference}.png"
    fig_sweep.savefig(str(sweep_png), dpi=dpi)

    out = {"sweep_panels": sweep_png}

    if sample_comparison_csv is not None:
        sample_df = load_cosmo_table(sample_comparison_csv)
        fig_sample, _ = plot_cosmo_sample_residuals(sample_df)
        sample_png = out_dir / f"{stem}_samples.png"
        fig_sample.savefig(str(sample_png), dpi=dpi)
        out["sample_residuals"] = sample_png

    return out


def render_cosmo_transient_figures(
    sweep_csvs: Union[PathLike, Sequence[PathLike]],
    out_dir: PathLike,
    *,
    value_prefix: str = "approach_frac_",
    families: Sequence[str] = ("concentration", "erate_ratio"),
    concentration_xlim: tuple = (0.0, 10000.0),
    erate_ratio_xlim: tuple = (0.0, 50000.0),
    show_title: bool = True,
    dpi: int = 400,
) -> dict:
    """Render the transient (step-change) cosmogenic figure families.

    One PNG per family per sweep CSV. At the default ``figsize=(5.0, 3.0)``
    and ``dpi=400`` these come out at 2000x1200 px, matching the published
    figure family.

    Args:
        sweep_csvs: One path, or a sequence of paths, to
            ``cosmo_transient_sweep_<label>.csv`` files. Each is plotted
            separately (one ``run_label`` per figure).
        out_dir: Output directory for figure PNGs (created if absent).
        value_prefix: Family A normalisation column prefix; see
            :data:`lsdviztools3.charts.cosmo_transient.CONCENTRATION_PREFIXES`.
        families: Which families to render: ``"concentration"`` (A) and/or
            ``"erate_ratio"`` (B).
        concentration_xlim: Hard x limits for Family A, in years.
        erate_ratio_xlim: Hard x limits for Family B, in years.
        show_title: Draw the standard metadata title on each figure.
        dpi: PNG output dpi (400 gives the 2000x1200 px target).

    Returns:
        Dict of ``"<stem>_<family>" -> Path`` for every figure written.
    """
    from .charts import (
        load_cosmo_transient_table,
        plot_cosmo_transient_concentration,
        plot_cosmo_transient_erate_ratio,
    )

    known = {"concentration", "erate_ratio"}
    unknown = set(families) - known
    if unknown:
        raise ValueError(f"Unknown figure families: {sorted(unknown)}")
    if not families:
        raise ValueError("No figure families requested.")

    if isinstance(sweep_csvs, (str, Path)):
        sweep_csvs = [sweep_csvs]

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written: dict = {}
    for csv_path in sweep_csvs:
        csv_path = Path(csv_path)
        df = load_cosmo_transient_table(csv_path)
        stem = csv_path.stem

        if "concentration" in families:
            fig, _ = plot_cosmo_transient_concentration(
                df,
                value_prefix=value_prefix,
                xlim=concentration_xlim,
                show_title=show_title,
            )
            png = out_dir / f"{stem}_norm_conc.png"
            fig.savefig(str(png), dpi=dpi)
            written[f"{stem}_concentration"] = png

        if "erate_ratio" in families:
            fig, _ = plot_cosmo_transient_erate_ratio(
                df,
                xlim=erate_ratio_xlim,
                show_title=show_title,
            )
            png = out_dir / f"{stem}_app_erate_ratio.png"
            fig.savefig(str(png), dpi=dpi)
            written[f"{stem}_erate_ratio"] = png

    return written


def render_cosmo_transient_spatial_figures(
    joined_csvs: Union[PathLike, Sequence[PathLike]],
    out_dir: PathLike,
    *,
    ratio_cols: Sequence[str] = ("app_erate_ratio_C14_Be10",),
    families: Sequence[str] = ("vs_steptime", "vs_mean_erate"),
    overlay_sweeps: Optional[Mapping[str, Sequence[tuple]]] = None,
    cmap: str = "jet",
    dpi: int = 320,
) -> dict:
    """Render the per-sample transient cosmogenic scatter families.

    One PNG per family per ratio column per joined CSV. At the default
    ``figsize=(10.0, 5.0)`` and ``dpi=320`` these come out at 3200x1600 px,
    matching the published figure family.

    Args:
        joined_csvs: One path, or a sequence of paths, to
            ``<stem>_joined_samples.csv`` files. Each is plotted separately.
        out_dir: Output directory for figure PNGs (created if absent).
        ratio_cols: Which apparent-rate ratio columns to draw; see
            :data:`lsdviztools3.charts.cosmo_transient_spatial.RATIO_LABELS`.
        families: Which families to render: ``"vs_steptime"`` (B),
            ``"vs_mean_erate"`` (C) and/or ``"overlay_1d"`` (D).
        overlay_sweeps: Required when ``"overlay_1d"`` is in ``families``.
            Maps a joined CSV's **stem** to a sequence of
            ``(curve_label, sweep_csv_path)`` pairs to draw over that run's
            samples. Keyed per stem rather than shared, because the sweep that
            matches one run's erosion rates does not match another's; a stem
            with no entry raises rather than silently rendering a bare scatter
            that looks like an overlay.
        cmap: Colour map for the colour/size encoding.
        dpi: PNG output dpi (320 gives the 3200x1600 px target).

    Returns:
        Dict of ``"<stem>_<ratio>_<family>" -> Path`` for every figure written.
    """
    from .charts import (
        load_cosmo_transient_table,
        load_joined_samples,
        plot_ratio_overlay_1d,
        plot_ratio_vs_mean_erate,
        plot_ratio_vs_steptime,
    )

    known = {"vs_steptime", "vs_mean_erate", "overlay_1d"}
    unknown = set(families) - known
    if unknown:
        raise ValueError(f"Unknown figure families: {sorted(unknown)}")
    if not families:
        raise ValueError("No figure families requested.")
    if not ratio_cols:
        raise ValueError("No ratio columns requested.")
    if "overlay_1d" in families and not overlay_sweeps:
        raise ValueError(
            "families includes 'overlay_1d' but overlay_sweeps is empty; "
            "there would be no 1-D curve to overlay."
        )

    if isinstance(joined_csvs, (str, Path)):
        joined_csvs = [joined_csvs]

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # The ratio column name is long and repeats the whole prefix; the figure
    # filenames use just the nuclide pair so they stay readable.
    def _short(col: str) -> str:
        return col.replace("app_erate_ratio_", "")

    plotters = {
        "vs_steptime": plot_ratio_vs_steptime,
        "vs_mean_erate": plot_ratio_vs_mean_erate,
    }

    written: dict = {}
    for csv_path in joined_csvs:
        csv_path = Path(csv_path)
        df = load_joined_samples(csv_path)
        stem = csv_path.stem

        # Loaded once per joined CSV rather than once per ratio column: the
        # sweep tables are large and carry every ratio already.
        sweeps: list = []
        if "overlay_1d" in families:
            entries = (overlay_sweeps or {}).get(stem)
            if not entries:
                raise KeyError(
                    f"No overlay_sweeps entry for '{stem}'. Every joined CSV "
                    "needs its own matched sweep; give one, or drop "
                    "'overlay_1d' from families."
                )
            sweeps = [
                (label, load_cosmo_transient_table(Path(p)))
                for label, p in entries
            ]

        for ratio_col in ratio_cols:
            for family in families:
                if family == "overlay_1d":
                    fig, _ = plot_ratio_overlay_1d(
                        df, sweeps, ratio_col=ratio_col, cmap=cmap,
                    )
                else:
                    fig, _ = plotters[family](
                        df, ratio_col=ratio_col, cmap=cmap,
                    )
                png = out_dir / f"{stem}_{_short(ratio_col)}_{family}.png"
                fig.savefig(str(png), dpi=dpi)
                written[f"{stem}_{_short(ratio_col)}_{family}"] = png

    return written


def render_cosmo_published_residuals(
    published_dir: PathLike,
    out_dir: PathLike,
    *,
    datasets: list[str],
    mode: str = "auto",
    stem: str = "cosmo_published",
    dpi: int = 220,
) -> dict:
    """Render residual boxplots from lsdtt3 published-summary CSVs.

    Args:
        published_dir: Directory containing `<dataset>_summary.csv` style files.
        out_dir: Output directory for PNG figures.
        datasets: Dataset names (e.g. `dethier2014`, `clapp2002`).
        mode: Which summary flavour to use: `stone`, `lsdn`, `riversand`, or
            `auto` (prefer riversand, then lsdn, then stone per dataset).
        stem: Output filename stem.
        dpi: PNG output dpi.

    Returns:
        Dict with written paths and selected source CSVs.
    """
    from .charts import load_cosmo_table, plot_cosmo_sample_residuals

    if not datasets:
        raise ValueError("datasets must contain at least one dataset name")
    mode = mode.lower()
    if mode not in {"auto", "stone", "lsdn", "riversand"}:
        raise ValueError("mode must be one of: auto, stone, lsdn, riversand")

    published_dir = Path(published_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _candidates(ds: str) -> list[Path]:
        if mode == "stone":
            return [published_dir / f"{ds}_summary.csv"]
        if mode == "lsdn":
            return [published_dir / f"{ds}_lsdn_summary.csv"]
        if mode == "riversand":
            return [published_dir / f"{ds}_riversand_summary.csv"]
        return [
            published_dir / f"{ds}_riversand_summary.csv",
            published_dir / f"{ds}_lsdn_summary.csv",
            published_dir / f"{ds}_summary.csv",
        ]

    frames = []
    used = []
    for ds in datasets:
        chosen = next((p for p in _candidates(ds) if p.is_file()), None)
        if chosen is None:
            raise FileNotFoundError(
                f"No published summary found for dataset '{ds}' in mode '{mode}' under {published_dir}"
            )
        df = load_cosmo_table(chosen)
        if "dataset" not in df.columns:
            df = df.copy()
            df["dataset"] = ds
        if "scaling_scheme" not in df.columns:
            if chosen.name.endswith("_lsdn_summary.csv"):
                df = df.copy()
                df["scaling_scheme"] = "lsdn"
            elif chosen.name.endswith("_riversand_summary.csv"):
                df = df.copy()
                df["scaling_scheme"] = "riversand"
            else:
                df = df.copy()
                df["scaling_scheme"] = "stone"
        frames.append(df)
        used.append(chosen)

    if len(frames) == 1:
        all_df = frames[0]
    else:
        import pandas as pd
        all_df = pd.concat(frames, ignore_index=True)

    fig, _ = plot_cosmo_sample_residuals(all_df)
    out_png = out_dir / f"{stem}_residuals.png"
    fig.savefig(str(out_png), dpi=dpi)

    return {
        "sample_residuals": out_png,
        "source_csvs": used,
    }


def _pick_column(df, requested, candidates, label):
    if requested is not None:
        if requested not in df.columns:
            raise KeyError(f"{label} column not found: {requested}")
        return requested
    lower_map = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    raise KeyError(
        f"Could not infer {label} column. Available columns: {list(df.columns)}"
    )


def _infer_basin_id_column(gdf, requested=None):
    if requested is not None:
        if requested not in gdf.columns:
            raise KeyError(f"Basin ID column not found: {requested}")
        return requested

    candidates = [
        "basin_index",
        "basin_key",
        "junction",
        "id",
        "outlet_junction",
    ]
    lower_map = {str(c).lower(): c for c in gdf.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def join_basins_with_cosmo_erates(
    basins: PathLike,
    cosmo_csv: PathLike,
    *,
    basin_id_col: Optional[str] = None,
    erate_col: Optional[str] = None,
    lon_col: Optional[str] = None,
    lat_col: Optional[str] = None,
    sample_col: Optional[str] = None,
    out_erate_col: str = "inferred_erate_mm_yr",
    predicate: str = "intersects",
    erate_input_units: str = "auto",
    assign_unmatched_to_nearest: bool = True,
):
    """Spatially join cosmo sample erosion rates onto basin polygons.

    Aggregates per-basin means and propagated erosion-rate uncertainty estimates:
    ``erate_err_low_mm_yr`` and ``erate_err_high_mm_yr``.

    Unit handling:
    - ``erate_input_units='auto'`` treats lsdtt-cosmo-tool's legacy
      ``apparent_e_*_mm_yr`` columns as mm/kyr and converts to mm/yr.
    - ``'mm_per_kyr'`` always divides by 1000.
    - ``'mm_per_yr'`` uses values unchanged.
    """
    import geopandas as gpd
    import pandas as pd

    def _quadrature_mean(series):
        values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return np.nan
        return float(np.sqrt(np.sum(finite * finite)) / finite.size)

    basins_gdf = basins.copy() if hasattr(basins, "geometry") else load_vector(basins)
    if basins_gdf.empty:
        raise ValueError("Cannot join erosion rates: basin polygons are empty.")

    cosmo_df = pd.read_csv(cosmo_csv)
    if cosmo_df.empty:
        raise ValueError("Cannot join erosion rates: cosmo CSV has no rows.")

    lon_name = _pick_column(
        cosmo_df,
        lon_col,
        ["longitude", "lon", "x", "outlet_longitude"],
        "longitude",
    )
    lat_name = _pick_column(
        cosmo_df,
        lat_col,
        ["latitude", "lat", "y", "outlet_latitude"],
        "latitude",
    )
    erate_name = _pick_column(
        cosmo_df,
        erate_col,
        [
            "apparent_e_mm_yr",
            "inferred_erate_mm_yr",
            "erosion_rate_mm_yr",
            "erate_mm_yr",
            "erate",
        ],
        "erosion rate",
    )

    unit_mode = erate_input_units.lower().strip()
    if unit_mode not in {"auto", "mm_per_yr", "mm_per_kyr"}:
        raise ValueError("erate_input_units must be one of: auto, mm_per_yr, mm_per_kyr")
    if unit_mode == "mm_per_kyr":
        erate_scale = 1.0 / 1000.0
    elif unit_mode == "mm_per_yr":
        erate_scale = 1.0
    else:
        erate_scale = 1.0 / 1000.0 if erate_name in {
            "apparent_e_mm_yr",
            "apparent_e_minus_mm_yr",
            "apparent_e_plus_mm_yr",
        } else 1.0

    sample_name_col = None
    if sample_col is not None and sample_col in cosmo_df.columns:
        sample_name_col = sample_col
    elif sample_col is None:
        for cand in ["sample_name", "sample", "name", "id"]:
            if cand in cosmo_df.columns:
                sample_name_col = cand
                break

    unc_sym_col = next(
        (
            c
            for c in ["apparent_e_unc_ext_mm_yr", "apparent_e_unc_int_mm_yr", "erate_unc_mm_yr"]
            if c in cosmo_df.columns
        ),
        None,
    )
    minus_col = next((c for c in ["apparent_e_minus_mm_yr", "erate_minus_mm_yr"] if c in cosmo_df.columns), None)
    plus_col = next((c for c in ["apparent_e_plus_mm_yr", "erate_plus_mm_yr"] if c in cosmo_df.columns), None)

    keep = [lon_name, lat_name, erate_name]
    for extra in [sample_name_col, unc_sym_col, minus_col, plus_col]:
        if extra is not None and extra not in keep:
            keep.append(extra)
    pts = cosmo_df[keep].copy()

    if sample_name_col is not None:
        pts[sample_name_col] = cosmo_df[sample_name_col]

    pts[lon_name] = pd.to_numeric(pts[lon_name], errors="coerce")
    pts[lat_name] = pd.to_numeric(pts[lat_name], errors="coerce")
    pts[erate_name] = pd.to_numeric(pts[erate_name], errors="coerce")
    pts[erate_name] = pts[erate_name] * erate_scale

    if unc_sym_col is not None:
        unc = pd.to_numeric(pts[unc_sym_col], errors="coerce").abs()
        if unit_mode == "mm_per_kyr" or (unit_mode == "auto" and unc_sym_col.startswith("apparent_e_unc_")):
            unc = unc / 1000.0
        pts["erate_err_low_mm_yr"] = unc
        pts["erate_err_high_mm_yr"] = unc
    elif minus_col is not None and plus_col is not None:
        minus = pd.to_numeric(pts[minus_col], errors="coerce")
        plus = pd.to_numeric(pts[plus_col], errors="coerce")
        if erate_scale != 1.0:
            minus = minus * erate_scale
            plus = plus * erate_scale
        central = pts[erate_name]
        pts["erate_err_low_mm_yr"] = (central - minus).clip(lower=0)
        pts["erate_err_high_mm_yr"] = (plus - central).clip(lower=0)
    else:
        pts["erate_err_low_mm_yr"] = np.nan
        pts["erate_err_high_mm_yr"] = np.nan

    pts = pts.dropna(subset=[lon_name, lat_name, erate_name])
    if pts.empty:
        raise ValueError("No finite lon/lat/erosion-rate rows available for joining.")

    points_gdf = gpd.GeoDataFrame(
        pts,
        geometry=gpd.points_from_xy(pts[lon_name], pts[lat_name]),
        crs="EPSG:4326",
    )

    if basins_gdf.crs is not None:
        points_gdf = points_gdf.to_crs(basins_gdf.crs)

    inferred_col = _infer_basin_id_column(basins_gdf, basin_id_col)
    join_key = inferred_col if inferred_col else "_basin_join_id"
    basins_join = basins_gdf.copy()
    if inferred_col is None:
        basins_join[join_key] = basins_join.index.astype(int)

    right = basins_join[[join_key, basins_join.geometry.name]]
    joined = gpd.sjoin(points_gdf, right, how="left", predicate=predicate)

    if assign_unmatched_to_nearest and joined[join_key].isna().any():
        unmatched = joined[joined[join_key].isna()].copy()
        if not unmatched.empty:
            left = unmatched
            polys = basins_join[[join_key, basins_join.geometry.name]].copy()
            if left.crs is not None and str(left.crs).lower().endswith("4326"):
                left = left.to_crs("EPSG:3857")
                polys = polys.to_crs("EPSG:3857")
            assigned = []
            for geom in left.geometry:
                d = polys.geometry.distance(geom)
                nearest_idx = d.idxmin()
                assigned.append(polys.loc[nearest_idx, join_key])
            joined.loc[unmatched.index, join_key] = assigned

    joined = joined.dropna(subset=[join_key])
    if joined.empty:
        raise ValueError("No cosmo samples intersect basin polygons.")

    grouped = (
        joined.groupby(join_key, dropna=False)
        .agg(
            **{
                out_erate_col: (erate_name, "mean"),
                "n_cosmo_samples": (erate_name, "count"),
                "erate_err_low_mm_yr": ("erate_err_low_mm_yr", _quadrature_mean),
                "erate_err_high_mm_yr": ("erate_err_high_mm_yr", _quadrature_mean),
            }
        )
        .reset_index()
    )
    if sample_name_col is not None:
        sample_group = (
            joined.groupby(join_key, dropna=False)[sample_name_col]
            .apply(lambda s: ",".join(sorted({str(v) for v in s if pd.notna(v)})))
            .reset_index(name="sample_names")
        )
        grouped = grouped.merge(sample_group, on=join_key, how="left")

    out = basins_join.merge(grouped, on=join_key, how="left")
    out["n_cosmo_samples"] = out["n_cosmo_samples"].fillna(0).astype(int)
    if inferred_col is None:
        out = out.drop(columns=[join_key])
    return out


def render_cosmo_basin_erate_map(
    base_raster: PathLike,
    basin_polygons: PathLike,
    cosmo_csv: PathLike,
    out_png: PathLike,
    *,
    hillshade_raster: Optional[PathLike] = None,
    joined_basins_out: Optional[PathLike] = None,
    basin_id_col: Optional[str] = None,
    erate_col: Optional[str] = None,
    lon_col: Optional[str] = None,
    lat_col: Optional[str] = None,
    sample_col: Optional[str] = None,
    out_erate_col: str = "inferred_erate_mm_yr",
    scaling_scheme: str = "stone",
    muon_scheme: str = "BraucherBorchers",
    atmosphere: str = "ERA40",
    plot_title: str = "",
    coordinates: str = "local",
    transparency: float = 65.0,
    fill: str = "orange",
    colorbar: bool = True,
    colorbar_label: str = "Inferred erosion rate (mm/yr)",
    style: Optional[MapStyle] = None,
    erate_input_units: str = "auto",
    require_complete_basins: bool = True,
    complete_edge_buffer_pixels: int = 3,
    complete_nodata_buffer_pixels: int = 1,
):
    """Render the standard hillshade + basin-erate overlay product.

    Returns:
        Dict with ``out_png``, optional ``joined_basins``, and ``n_joined_basins``.
    """
    style = style or MapStyle()
    if style.cmap == "geo":
        style = MapStyle(**{**style.__dict__, "cmap": "viridis"})
    basins_gdf = load_vector(basin_polygons)
    if require_complete_basins:
        _ensure_complete_basins(
            basins_gdf,
            base_raster,
            edge_buffer_pixels=complete_edge_buffer_pixels,
            nodata_buffer_pixels=complete_nodata_buffer_pixels,
            basin_id_col=basin_id_col,
        )

    joined = join_basins_with_cosmo_erates(
        basins_gdf,
        cosmo_csv,
        basin_id_col=basin_id_col,
        erate_col=erate_col,
        lon_col=lon_col,
        lat_col=lat_col,
        sample_col=sample_col,
        out_erate_col=out_erate_col,
        erate_input_units=erate_input_units,
    )

    fig = _render_hillshade_base(
        base_raster,
        hillshade_raster,
        coordinates=coordinates,
        style=style,
    )
    extent = _display_extent_for_raster(base_raster, coordinates=coordinates)
    fig = render_basins(
        joined,
        coordinates=coordinates,
        fill=fill,
        transparency=transparency,
        color_by=out_erate_col,
        colorbar=colorbar,
        colorbar_label=colorbar_label,
        cmap=style.cmap,
        style=style,
        fig=fig,
    )

    title_prefix = plot_title.strip()
    header = []
    if title_prefix:
        header.append(title_prefix)
    header.extend([
        f"scaling: {scaling_scheme}",
        f"atmosphere: {atmosphere}",
        f"muon: {muon_scheme}",
    ])
    _add_header_lines(
        fig,
        header,
        extent,
    )

    save_figure(fig, out_png, style=style)

    out = {
        "out_png": Path(out_png),
        "n_joined_basins": int((joined["n_cosmo_samples"] > 0).sum()),
    }
    if joined_basins_out is not None:
        joined_path = Path(joined_basins_out)
        joined_path.parent.mkdir(parents=True, exist_ok=True)
        joined.to_file(joined_path)
        out["joined_basins"] = joined_path
    return out


def summarize_basin_ksn(
    fit_nodes_csv: PathLike,
    fit_segments_csv: PathLike,
    *,
    basins: Optional[Sequence[int]] = None,
):
    """Summarize per-basin k_sn statistics from chi fit tables.

    Uses per-node ``k_sn`` values (segment value assigned to each chi node) and
    reports basin-wise median, 25th and 75th percentiles.
    """
    import pandas as pd

    from .io.chi import load_chi_fit

    nodes = load_chi_fit(fit_nodes_csv, fit_segments_csv)
    if basins is not None:
        basin_set = {int(b) for b in basins}
        nodes = nodes[nodes["basin_index"].isin(basin_set)]

    nodes = nodes[np.isfinite(pd.to_numeric(nodes["k_sn"], errors="coerce"))]
    if nodes.empty:
        raise ValueError("No finite k_sn values were found in chi fit nodes.")

    grouped = (
        nodes.groupby("basin_index", dropna=False)["k_sn"]
        .agg(
            ksn_median="median",
            ksn_q25=lambda s: float(np.nanpercentile(s, 25)),
            ksn_q75=lambda s: float(np.nanpercentile(s, 75)),
            n_ksn_nodes="count",
        )
        .reset_index()
        .sort_values("basin_index")
    )
    grouped["ksn_err_low"] = grouped["ksn_median"] - grouped["ksn_q25"]
    grouped["ksn_err_high"] = grouped["ksn_q75"] - grouped["ksn_median"]
    return grouped


def render_chi_basin_ksn_map(
    base_raster: PathLike,
    basin_polygons: PathLike,
    fit_nodes_csv: PathLike,
    fit_segments_csv: PathLike,
    out_png: PathLike,
    *,
    hillshade_raster: Optional[PathLike] = None,
    joined_basins_out: Optional[PathLike] = None,
    basin_id_col: Optional[str] = None,
    coordinates: str = "local",
    transparency: float = 65.0,
    colorbar: bool = True,
    colorbar_label: str = "Median k@-sn@-",
    plot_title: str = "",
    scaling_scheme: str = "stone",
    atmosphere: str = "ERA40",
    muon_scheme: str = "BraucherBorchers",
    style: Optional[MapStyle] = None,
    theta: Optional[float] = None,
    require_complete_basins: bool = True,
    complete_edge_buffer_pixels: int = 3,
    complete_nodata_buffer_pixels: int = 1,
):
    """Render hillshade + basins colored by basin-median k_sn."""
    import geopandas as gpd
    import pandas as pd
    import rasterio
    from rasterio.transform import xy as raster_xy
    from .io.table import load_table

    style = style or MapStyle()
    if style.cmap == "geo":
        style = MapStyle(**{**style.__dict__, "cmap": "viridis"})
    ksn = summarize_basin_ksn(fit_nodes_csv, fit_segments_csv)

    basins_gdf = load_vector(basin_polygons)
    if require_complete_basins:
        _ensure_complete_basins(
            basins_gdf,
            base_raster,
            edge_buffer_pixels=complete_edge_buffer_pixels,
            nodata_buffer_pixels=complete_nodata_buffer_pixels,
            basin_id_col=basin_id_col,
        )

    inferred_col = _infer_basin_id_column(basins_gdf, basin_id_col)
    basins_join = basins_gdf.copy()

    if inferred_col is not None:
        join_key = inferred_col
        ksn_df = ksn.rename(columns={"basin_index": join_key})
        joined = basins_join.merge(ksn_df, on=join_key, how="left")
    else:
        fit_nodes = load_table(fit_nodes_csv)
        outlet = (
            fit_nodes[["basin_index", "outlet_flat"]]
            .dropna()
            .drop_duplicates(subset=["basin_index"])
        )
        if outlet.empty:
            raise ValueError("Cannot map k_sn to basin polygons: fit_nodes has no basin_index/outlet_flat rows.")

        with rasterio.open(base_raster) as src:
            n_cols = int(src.width)
            transform = src.transform
            raster_crs = src.crs

        rows = (outlet["outlet_flat"].astype(int) // n_cols).to_numpy()
        cols = (outlet["outlet_flat"].astype(int) % n_cols).to_numpy()
        xs = []
        ys = []
        for r, c in zip(rows, cols):
            x, y = raster_xy(transform, int(r), int(c), offset="center")
            xs.append(float(x))
            ys.append(float(y))

        outlets_gdf = gpd.GeoDataFrame(
            outlet[["basin_index"]].copy(),
            geometry=gpd.points_from_xy(xs, ys),
            crs=raster_crs,
        )
        if basins_join.crs is not None and outlets_gdf.crs is not None and basins_join.crs != outlets_gdf.crs:
            outlets_gdf = outlets_gdf.to_crs(basins_join.crs)

        poly_idx_col = "_poly_idx"
        basins_join[poly_idx_col] = basins_join.index.astype(int)
        matched = gpd.sjoin(
            outlets_gdf,
            basins_join[[poly_idx_col, basins_join.geometry.name]],
            how="left",
            predicate="within",
        ).dropna(subset=[poly_idx_col])

        if matched.empty:
            raise ValueError("Could not map chi basin outlets into supplied basin polygons.")

        basin_map = matched[[poly_idx_col, "basin_index"]].drop_duplicates(subset=[poly_idx_col])
        joined = basins_join.merge(basin_map, on=poly_idx_col, how="left")
        joined = joined.merge(ksn, on="basin_index", how="left")
        joined = joined.drop(columns=[poly_idx_col])

    fig = _render_hillshade_base(
        base_raster,
        hillshade_raster,
        coordinates=coordinates,
        style=style,
    )
    extent = _display_extent_for_raster(base_raster, coordinates=coordinates)
    fig = render_basins(
        joined,
        coordinates=coordinates,
        fill="orange",
        transparency=transparency,
        color_by="ksn_median",
        colorbar=colorbar,
        colorbar_label=colorbar_label,
        cmap=style.cmap,
        style=style,
        fig=fig,
    )

    title_text = f"@~q@~={theta:g}" if theta is not None else ""
    if plot_title.strip():
        title_text = plot_title.strip()
    _add_header_lines(
        fig,
        [title_text],
        extent,
    )

    save_figure(fig, out_png, style=style)

    out = {
        "out_png": Path(out_png),
        "n_joined_basins": int(pd.to_numeric(joined["n_ksn_nodes"], errors="coerce").fillna(0).gt(0).sum()),
        "ksn_summary": ksn,
        "joined_basins": joined,
    }
    if joined_basins_out is not None:
        joined_path = Path(joined_basins_out)
        joined_path.parent.mkdir(parents=True, exist_ok=True)
        joined.to_file(joined_path)
        out["joined_basins_path"] = joined_path
    return out


def render_chi_and_cosmo_basin_workflow(
    outputs: dict,
    out_dir: PathLike,
    *,
    stem: str = "chi",
    dpi: int = 300,
    basin_polygons: Optional[PathLike] = None,
    base_raster: Optional[PathLike] = None,
    basin_id_col: Optional[str] = None,
    cosmo_csv: Optional[PathLike] = None,
    erate_col: Optional[str] = None,
    lon_col: Optional[str] = None,
    lat_col: Optional[str] = None,
    sample_col: Optional[str] = None,
    theta: Optional[float] = None,
    scaling_scheme: str = "stone",
    atmosphere: str = "ERA40",
    muon_scheme: str = "BraucherBorchers",
    erate_input_units: str = "auto",
    hillshade_raster: Optional[PathLike] = None,
    require_complete_basins: bool = True,
    complete_edge_buffer_pixels: int = 3,
    complete_nodata_buffer_pixels: int = 1,
    style: Optional[MapStyle] = None,
):
    """Render per-basin chi PNGs, basin k_sn summary, and optional k_sn-vs-cosmo.

    If basin polygons and a base raster are provided, also writes a basin-median
    k_sn map. If cosmo CSV is also provided, writes a k_sn-vs-erate comparison
    using k_sn IQR (25-75th) and cosmo-reported erosion-rate uncertainties.
    """
    import pandas as pd

    from .charts import plot_chi_profile, plot_chi_stacked, plot_ksn_vs_erate
    from .io.chi import load_chi_fit

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_chi_fit(outputs["fit_nodes"], outputs["fit_segments"])
    basins = sorted({int(b) for b in nodes["basin_index"].dropna().unique()})

    written = {
        "chi_profiles": [],
    }

    fig_stacked, _ = plot_chi_stacked(nodes, basins=basins)
    stacked_png = out_dir / f"{stem}_stacked.png"
    fig_stacked.savefig(stacked_png, dpi=dpi, bbox_inches="tight")
    written["chi_stacked"] = stacked_png

    for basin in basins:
        fig, _ = plot_chi_profile(nodes, basin=basin)
        path = out_dir / f"{stem}_profile_basin{basin}.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written["chi_profiles"].append(path)

    ksn_summary = summarize_basin_ksn(outputs["fit_nodes"], outputs["fit_segments"], basins=basins)
    ksn_csv = out_dir / f"{stem}_basin_ksn_summary.csv"
    ksn_summary.to_csv(ksn_csv, index=False)
    written["ksn_summary_csv"] = ksn_csv

    if basin_polygons is not None and base_raster is not None:
        ksn_map_png = out_dir / f"{stem}_basin_ksn.png"
        joined_fgb = out_dir / f"{stem}_basin_ksn.fgb"
        ksn_map = render_chi_basin_ksn_map(
            base_raster,
            basin_polygons,
            outputs["fit_nodes"],
            outputs["fit_segments"],
            ksn_map_png,
            hillshade_raster=hillshade_raster,
            joined_basins_out=joined_fgb,
            basin_id_col=basin_id_col,
            plot_title="",
            scaling_scheme=scaling_scheme,
            atmosphere=atmosphere,
            muon_scheme=muon_scheme,
            theta=theta,
            require_complete_basins=require_complete_basins,
            complete_edge_buffer_pixels=complete_edge_buffer_pixels,
            complete_nodata_buffer_pixels=complete_nodata_buffer_pixels,
            style=style,
        )
        written["ksn_map_png"] = ksn_map_png
        written["ksn_map_joined_basins"] = joined_fgb

        if cosmo_csv is not None:
            cosmo_joined = join_basins_with_cosmo_erates(
                basin_polygons,
                cosmo_csv,
                basin_id_col=basin_id_col,
                erate_col=erate_col,
                lon_col=lon_col,
                lat_col=lat_col,
                sample_col=sample_col,
                erate_input_units=erate_input_units,
            )

            if basin_id_col is None:
                inferred_id = _infer_basin_id_column(cosmo_joined)
                if inferred_id is None:
                    cosmo_joined = cosmo_joined.copy()
                    cosmo_joined["_basin_join_id"] = cosmo_joined.index.astype(int)
                    basin_join_key = "_basin_join_id"
                else:
                    basin_join_key = inferred_id
            else:
                basin_join_key = basin_id_col

            if basin_join_key not in ksn_map["joined_basins"].columns:
                ksn_join = ksn_map["joined_basins"].copy()
                ksn_join[basin_join_key] = ksn_join.index.astype(int)
            else:
                ksn_join = ksn_map["joined_basins"]

            compare = ksn_join.merge(
                pd.DataFrame(cosmo_joined.drop(columns=cosmo_joined.geometry.name)),
                on=basin_join_key,
                how="inner",
                suffixes=("_ksn", "_cosmo"),
            )
            compare = compare[
                np.isfinite(pd.to_numeric(compare["ksn_median"], errors="coerce"))
                & np.isfinite(pd.to_numeric(compare["inferred_erate_mm_yr"], errors="coerce"))
                & (pd.to_numeric(compare["n_cosmo_samples"], errors="coerce") > 0)
            ]

            if not compare.empty:
                fig_cmp, _ = plot_ksn_vs_erate(compare, theta=theta)
                cmp_png = out_dir / f"{stem}_ksn_vs_erate.png"
                fig_cmp.savefig(cmp_png, dpi=dpi, bbox_inches="tight")
                compare_csv = out_dir / f"{stem}_ksn_vs_erate.csv"
                compare.to_csv(compare_csv, index=False)
                written["ksn_vs_erate_png"] = cmp_png
                written["ksn_vs_erate_csv"] = compare_csv

    return written

def render_divide_migration_figures(
    prefix: PathLike,
    out_dir: PathLike,
    *,
    base_raster: Optional[PathLike] = None,
    metric: str = "chi",
    coordinates: str = "local",
    base_is_hillshade: bool = False,
    cmap: str = "roma",
    point_size: Optional[str] = None,
    divide_color: str = "black",
    alpha: float = 0.05,
    stem: str = "divide_migration",
    style: Optional[MapStyle] = None,
    high_res: bool = False,
    dpi: int = 220,
) -> dict:
    """Render the standard figure set for an ``lsdtt-divide-migration`` run.

    Reads an output set written by the CLI (nothing is recomputed) and writes:

    * ``<stem>_map_<metric>.png`` — grey hillshade + divide polylines + channel
      heads coloured by ``metric`` (only when ``base_raster`` is given);
    * ``<stem>_pairs_<metric>.png`` — side A vs side B means with Δ = A − B and
      its uncertainty, per divide pair (needs the pair-summary CSV);
    * ``<stem>_heads_<metric>.png`` — per-head metric vs distance to divide.

    Args:
        prefix: The ``<write_path>/<write_fname>`` stem, or a directory holding
            one output set (see ``io.divide_migration.find_outputs``).
        out_dir: Output directory for the PNGs.
        base_raster: DEM (or hillshade) for the map layer. Without it the map is
            skipped and only the charts are written.
        metric: ``chi`` / ``relief`` / ``gradient`` / ``head_elevation``.
        coordinates: Coordinate mode for the map layers.
        base_is_hillshade: ``base_raster`` is already a hillshade.
        cmap: Colour map for the channel-head metric.
        point_size: GMT symbol size for the heads (e.g. ``"0.12c"``).
        divide_color: Pen colour for the divide polylines.
        alpha: p-value threshold, used only when the run's uncertainty method is
            ``ttest``.
        stem: Filename stem for the outputs.
        style: Map styling options.
        high_res: Save the map at the high-res DPI.
        dpi: PNG dpi for the matplotlib charts.

    Returns:
        Dict of written figure paths, plus ``n_heads_nan_dropped`` when the map
        was drawn.
    """
    from .charts.divide_migration import (
        plot_head_metric_vs_distance,
        plot_pair_comparison,
    )
    from .io.divide_migration import (
        find_outputs, load_channel_heads, load_pair_summary, resolve_metric,
    )
    from .render.divide_migration import divide_migration_layers

    metric_key = resolve_metric(metric).key
    outputs = find_outputs(prefix)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = {}

    if base_raster is not None:
        if not outputs["divides"] or not outputs["channel_heads"]:
            raise FileNotFoundError(
                "The map needs both <stem>_divides and <stem>_channel_heads "
                f"vectors; found divides={outputs['divides']}, "
                f"heads={outputs['channel_heads']}."
            )
        fig, n_dropped = divide_migration_layers(
            base_raster,
            outputs["divides"],
            outputs["channel_heads"],
            metric=metric_key,
            coordinates=coordinates,
            base_is_hillshade=base_is_hillshade,
            cmap=cmap,
            point_size=point_size,
            divide_color=divide_color,
            style=style,
        )
        map_png = out_dir / f"{stem}_map_{metric_key}.png"
        save_figure(fig, str(map_png), style=style, high_res=high_res)
        written["map"] = map_png
        written["n_heads_nan_dropped"] = n_dropped

    if outputs["pair_summary"]:
        summary = load_pair_summary(outputs["pair_summary"])
        fig_pairs, _ = plot_pair_comparison(summary, metric=metric_key, alpha=alpha)
        pairs_png = out_dir / f"{stem}_pairs_{metric_key}.png"
        fig_pairs.savefig(str(pairs_png), dpi=dpi)
        written["pair_comparison"] = pairs_png

    heads_path = outputs["metrics_csv"] or outputs["channel_heads"]
    if heads_path:
        heads = load_channel_heads(heads_path)
        fig_heads, _ = plot_head_metric_vs_distance(heads, metric=metric_key)
        heads_png = out_dir / f"{stem}_heads_{metric_key}.png"
        fig_heads.savefig(str(heads_png), dpi=dpi)
        written["head_scatter"] = heads_png

    if not written:
        raise FileNotFoundError(
            f"Nothing to render for divide-migration prefix {prefix}."
        )
    return written
