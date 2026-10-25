"""Map of lsdtt3's erodibility-calibration sample basins (PyGMT).

:func:`render_calibration_basin_map` draws the basin polygons lsdtt3 wrote
(``*_erodibility_calibration_basins.fgb``) over a gray hillshade of the run's
DEM, coloured by the final sample status (from the samples CSV) or by the
dominant lithologic unit (shared name -> colour palette), with the snapped
outlets, a legend, a credit line, the minimal scale bar and a north arrow.
"""

from __future__ import annotations

from typing import Mapping, Optional

import numpy as np
import pandas as pd

from .style import MapStyle

#: Fixed colours per calibration status (lsdtt3's status vocabulary).
STATUS_COLOURS = {
    "used": "#1b9e77",
    "too_few_reaches": "#d95f02",
    "missing_uncertainty": "#e6ab02",
    "area_mismatch": "#7570b3",
    "nested_excluded": "#66a61e",
    "incomplete_basin": "#e7298a",
    "bad_rate": "#a6761d",
    "snap_failed": "#666667",
    "outside_dem": "#999998",
}
#: Status display order (used first).
STATUS_ORDER = tuple(STATUS_COLOURS)


def _legend_label(text: str) -> str:
    return str(text).replace("@", "@@")


def render_calibration_basin_map(
    samples: pd.DataFrame,
    basins,
    *,
    dem=None,
    colour_by: str = "status",
    unit_names=None,
    unit_colours: Optional[Mapping[str, str]] = None,
    credit: Optional[str] = None,
    display_downsample: int = 1,
    transparency: float = 25.0,
    credit_anchor: str = "TL",
    style: Optional[MapStyle] = None,
    fig=None,
):
    """Sample basins + outlets over a hillshade (or a plain frame without ``dem``).

    Args:
        samples: Samples table (final ``status``, snapped outlet lat/lon,
            ``dominant_lith_name``).
        basins: Basin polygons (GeoDataFrame or path), DEM CRS.
        dem: DEM for the hillshade; ``None`` draws a plain framed background.
        colour_by: ``"status"`` (final status from the samples file) or
            ``"unit"`` (dominant lithologic unit, shared palette).
        unit_names: All unit names of the run, so unit colours match the
            run's other figures.
        unit_colours: ``{name: colour}`` overrides.
        credit: Credit line(s) (``"\\n"`` separates lines).
        display_downsample: Integer block-mean factor applied to the DEM for
            display only (as ``workflows.hillshade_with_vectors``).
        transparency: Polygon fill transparency (percent).
        credit_anchor: Corner for the credit (default top-left, clear of the
            scale bar and arrow).
        style: Map style.
        fig: Existing figure.

    Returns:
        The ``pygmt.Figure``.
    """
    import io as _io
    import tempfile
    from dataclasses import replace

    import geopandas as gpd
    import pygmt

    from ..io.erodibility_calibration import BASIN_COLUMNS, SAMPLE_COLUMNS
    from ..io.raster import get_crs, load_raster, to_gmt_grid
    from ..io.vector import load_vector
    from .base import (add_crs_annotations, auto_projection, display_frame, gdf_to_display,
                       grid_to_display, new_figure)
    from .categorical import neutral_shading
    from .decorations import add_bar_scalebar, add_map_credit, add_north_arrow
    from .hillshade import render_hillshade
    from .lithology_colours import lithology_palette

    if colour_by not in ("status", "unit"):
        raise ValueError("colour_by must be 'status' or 'unit'.")
    style = style or MapStyle()
    c = SAMPLE_COLUMNS
    g = basins if isinstance(basins, gpd.GeoDataFrame) else load_vector(basins)
    g = g.copy()
    g[BASIN_COLUMNS["sample_id"]] = g[BASIN_COLUMNS["sample_id"]].astype(str)
    s = samples.copy()
    s[c["sample_id"]] = s[c["sample_id"]].astype(str)
    # Final status and dominant unit come from the samples file (the polygon
    # file's status is the pre-k_sn one).
    g = g.drop(columns=[BASIN_COLUMNS["status"]], errors="ignore").merge(
        s[[c["sample_id"], c["status"], c["dominant_lith_name"]]],
        left_on=BASIN_COLUMNS["sample_id"], right_on=c["sample_id"], how="left")

    if dem is not None:
        da = dem if hasattr(dem, "rio") else load_raster(dem)
        if display_downsample > 1:
            from ..workflows import _downsample_raster_for_display

            da = _downsample_raster_for_display(da, int(display_downsample))
            tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
            da.rio.to_raster(tmp.name)
            da = load_raster(tmp.name)
        crs = get_crs(da)
        with neutral_shading():
            fig = render_hillshade(da, shade=True, style=replace(style, cmap="lightgray"),
                                   fig=fig, decorate=False)
        _, extent = grid_to_display(to_gmt_grid(da), crs)
        projection = style.projection or auto_projection(extent, crs, style.figure_size)
    else:
        crs = g.crs
        gd = gdf_to_display(g, crs)
        w, so, e, n = gd.total_bounds
        pad = 0.08 * max(e - w, n - so)
        extent = [w - pad, e + pad, so - pad, n + pad]
        projection = style.projection or auto_projection(extent, crs, style.figure_size)
        fig = fig or new_figure(style)
        fig.basemap(region=extent, projection=projection, frame=display_frame(crs, style))
        if not style.clean:
            add_crs_annotations(fig, extent, crs, style)

    gd = gdf_to_display(g.to_crs(crs), crs)

    # Category -> colour.
    if colour_by == "status":
        present = set(s[c["status"]].astype(str))  # incl. samples with no polygon
        cats = [st for st in STATUS_ORDER if st in present]
        cats += sorted(present - set(cats))
        colours = {k: STATUS_COLOURS.get(k, "#555556") for k in cats}
        key = gd[c["status"]].astype(str)
        skey = s[c["status"]].astype(str)
        title = "sample status (lsdtt3)"
    else:
        names = gd[c["dominant_lith_name"]].dropna().astype(str)
        pal = lithology_palette(list(unit_names) if unit_names is not None else names,
                                unit_colours)
        for nm in names:
            if nm not in pal:
                pal = lithology_palette(list(pal) + [nm], unit_colours)
        key = gd[c["dominant_lith_name"]].fillna("no dominant unit").astype(str)
        skey = s[c["dominant_lith_name"]].fillna("no dominant unit").astype(str)
        cats = [n for n in pal if n in set(skey)] + (
            ["no dominant unit"] if "no dominant unit" in set(skey) else [])
        colours = {k: pal.get(k, "#555556") for k in cats}
        title = "dominant unit (lsdtt3)"

    # Largest basins first so nested/smaller ones stay visible on top.
    gd = gd.assign(_k=key.to_numpy(), _a=gd.geometry.area).sort_values("_a", ascending=False)
    for _, r in gd.iterrows():
        fig.plot(data=gpd.GeoDataFrame(geometry=[r.geometry], crs=gd.crs),
                 fill=colours[r["_k"]], pen="0.5p,black", transparency=transparency)

    # Snapped outlets (EPSG:4326 in the samples file) -> map CRS, display units.
    lat, lon = c["snapped_lat"], c["snapped_lon"]
    ok = s[lat].notna() & s[lon].notna()
    if ok.any():
        pts = gpd.GeoDataFrame(s[ok], geometry=gpd.points_from_xy(s.loc[ok, lon], s.loc[ok, lat]),
                               crs="EPSG:4326").to_crs(crs)
        pts = gdf_to_display(pts, crs)
        pkey = skey[ok].to_numpy()
        for k in cats:  # outlet filled with its category colour
            sel = pkey == k
            if sel.any():
                fig.plot(x=pts.geometry.x.to_numpy()[sel], y=pts.geometry.y.to_numpy()[sel],
                         style="t0.22c", fill=colours[k], pen="0.6p,black")

    rows = [f"H 8p,Helvetica-Bold {_legend_label(title)}"]
    rows += [f"S 0.15c s 0.28c {colours[k]} 0.25p,black 0.4c {_legend_label(k)}" for k in cats]
    rows.append("S 0.15c t 0.22c white 0.6p,black 0.4c snapped outlet (fill = category; "
                "no polygon = no complete basin)")
    with pygmt.config(FONT_ANNOT_PRIMARY="7p"):
        fig.legend(spec=_io.StringIO("\n".join(rows) + "\n"), position="JMR+jML+o0.3c/0c",
                   box="+gwhite+p0.5p,black")
    add_bar_scalebar(fig, extent, crs, label_position="below", label_halo="white")
    add_north_arrow(fig, extent, crs, projection=style.projection, label_halo="white")
    add_map_credit(fig, extent, credit or "", anchor=credit_anchor)
    return fig
