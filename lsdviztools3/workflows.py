"""Composite lsdtt3 -> lsdviztools3 workflows.

Convenience orchestration that drives lsdtt3 (via :mod:`lsdviztools3.driver`)
and renders the result (via :mod:`lsdviztools3.render`). Kept thin: each step is
a plain function call you could also run yourself.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.ops import unary_union

from .io.raster import load_raster, to_coordinates, to_gmt_grid
from .io.vector import load_vector, to_coordinates as to_vector_coordinates
from .driver import run_lsdtt
from .render.base import save_figure
from .render.basins import render_basins
from .render.channels import render_channels
from .render.hillshade import render_hillshade
from .render.points import render_points
from .render.style import MapStyle

PathLike = Union[str, Path]


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
    ``minimum_elevation=0``, write hillshade) -> ``lsdtt-channel-extraction``
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
             minimum_elevation=0.0, write_hillshade_raster="true"),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )
    sea_removed = f"{stem}_remove_seas.tif"

    run_lsdtt(
        "lsdtt-channel-extraction",
        dict(base, read_path=str(workdir), read_fname=sea_removed,
             source_pixels_threshold=source_threshold,
             write_channel_network_lines="true"),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )

    run_lsdtt(
        "lsdtt-basin-extraction",
        dict(base, read_path=str(workdir), read_fname=sea_removed,
             source_pixels_threshold=source_threshold,
             select_all_basins="true", filter_all_but_largest_basin="true",
             write_basin_outlines="true", write_basin_outlines_as_polygon="true"),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )

    return {
        "hillshade": workdir / f"{stem}_hillshade.tif",
        "channels": workdir / f"{stem}_channel_network_lines.fgb",
        "basin": workdir / f"{stem}_junction_basin_polygons.fgb",
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
    ``minimum_elevation=0``) -> ``lsdtt-chi-analysis`` with the N-largest-complete
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
             minimum_elevation=0.0),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )
    sea_removed = f"{stem}_remove_seas.tif"

    run_lsdtt(
        "lsdtt-chi-analysis",
        dict(base, read_path=str(workdir), read_fname=sea_removed,
             source_pixels_threshold=source_threshold,
             n_largest_basins=n_largest_basins,
             n_pixels_move_outlet_upstream=n_pixels_move_outlet_upstream,
             print_chi_data_maps="true",
             chi_data_maps_output_format="flatgeobuf",
             fit_chi_networks="true",
             write_fit_tables="true",
             write_chi_fit_segment_lines="true",
             chi_fit_segment_lines_output_format="flatgeobuf",
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

    Returns:
        List of written PNG Paths.
    """
    from .io.chi import load_chi_fit
    from .charts import plot_chi_profile, plot_chi_stacked

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_chi_fit(outputs["fit_nodes"], outputs["fit_segments"])
    written = []

    if stacked:
        fig, _ = plot_chi_stacked(nodes)
        path = out_dir / f"{stem}_stacked.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written.append(path)

    if per_basin:
        for b in sorted(nodes["basin_index"].unique()):
            fig, _ = plot_chi_profile(nodes, basin=int(b))
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
        coordinates: Coordinate mode (``"local"`` UTM-km / ``"geographic"``).
        width: Figure width (defaults to the 17 cm max).
        cmap: CPT for the log10 k_sn colour scale.
        high_res: Save at high-res DPI.

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
    fig = render_channels(
        str(segment_lines), color_by="k_sn", log_color=True, cmap=cmap,
        colorbar=True, colorbar_label="log@-10@-k@-sn@-",
        coordinates=coordinates,
        style=MapStyle(line_width="1.2p", cmap=cmap, figure_size=width),
        fig=fig,
    )
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
        coordinates: Coordinate mode (``"local"`` UTM-km recommended).
        azimuth: Viewing azimuth in degrees.
        elevation: Viewing elevation angle in degrees.
        vertical_exaggeration: True vertical:horizontal scale ratio.
        cmap: CPT for the log10 k_sn colour scale.
        width: Figure width.
        high_res: Save at high-res DPI.

    Returns:
        ``Path`` to the written PNG.
    """
    import numpy as np
    import pygmt
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
    gdf = gdf_to_display(to_vector_coordinates(load_vector(segment_lines), coordinates),
                         crs)
    ksn = gdf["k_sn"].to_numpy(dtype=float)
    z = np.full_like(ksn, np.nan)
    z[ksn > 0] = np.log10(ksn[ksn > 0])
    finite = z[np.isfinite(z)]
    if finite.size:
        pygmt.makecpt(cmap=cmap, series=[float(finite.min()), float(finite.max())])

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
        add_colourbar(fig, label="log@-10@-k@-sn@-")

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
             minimum_elevation=0.0),
        data_path=str(workdir), lsdtt_dir=lsdtt_dir,
    )
    sea_removed = f"{stem}_remove_seas.tif"

    # Largest complete basin's channel network as primary-path lines.
    run_lsdtt(
        "lsdtt-chi-analysis",
        dict(base, read_path=str(workdir), read_fname=sea_removed,
             source_pixels_threshold=source_threshold,
             n_largest_basins=1,
             print_chi_data_maps="true",
             chi_data_maps_output_format="flatgeobuf",
             write_chi_data_map_lines="true",
             chi_data_map_lines_output_format="flatgeobuf",
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

    channels = _clip_to_valid(
        to_vector_coordinates(load_vector(outputs["channels"]), coordinates),
        valid_union,
    )
    basins = _clip_to_valid(
        to_vector_coordinates(load_vector(outputs["basin"]), coordinates),
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
        coordinates: ``"local"`` (native UTM, default) / ``"geographic"`` /
            explicit CRS. Base and all layers are put in this CRS so they align.
        shade: treat ``base_raster`` as elevation (True) or a ready hillshade.
        base_cmap: colour map for the base (default ``"gray"`` for a clean
            backdrop; pass e.g. ``"geo"`` for a coloured elevation base).
        region: ``[xmin, xmax, ymin, ymax]`` in the raster's NATIVE CRS to zoom;
            crops the base raster and everything follows.
        clip_to_valid: clip overlays to the base raster's finite-data footprint
            (useful when a vector tile is larger than the raster).
        legend: draw a legend from the layers that have a ``label``.
        legend_position: GMT legend position string.
        fig: existing ``pygmt.Figure`` to compose onto (``None`` = new).

    Returns:
        The ``pygmt.Figure``.
    """
    # Optional zoom: crop the base raster (region in its NATIVE CRS) to a temp
    # GeoTIFF so the hillshade region and the valid-data clip follow the zoom.
    if region is not None:
        import tempfile
        xmin, xmax, ymin, ymax = region
        da = load_raster(str(base_raster)).rio.clip_box(xmin, ymin, xmax, ymax)
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
    for layer in layers:
        src = layer["path"]
        gdf = load_vector(str(src), layer.get("layer")) if not hasattr(
            src, "geometry") else src
        gdf = to_vector_coordinates(gdf, coordinates)
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
            fig = render_points(
                gdf, coordinates=coordinates,
                style=MapStyle(point_size=size, line_color=color, cmap=base_cmap
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
