"""Floodplain polygon rendering: outlines over a hillshade.

Built for comparing floodplain boundaries from lsdtt3 ``lsdtt-fluvial-surfaces``
(``write_floodplain_polygons``): the first-pass and second-pass polygons
(``floodplain_second_pass_algorithm_choice qq_upper_departure``), or the
floodplains of different methods. Each polygon set is drawn twice - a faint fill
so the hillshade stays readable underneath, then a strong outline so the
boundary is visible without hiding
the ground at it. An optional coloured drape (HAND, say) goes between the two.
"""

from __future__ import annotations

import io as _io
from dataclasses import replace
from typing import Optional, Sequence

from .channels import _prepare_overlay
from .hillshade import render_hillshade
from .style import MapStyle

# Distinct on a gray hillshade and to most colour-vision deficiencies.
DEFAULT_COLOURS = ("#1f6fd1", "#e8890c", "#2ca02c", "#c51b7d")


def render_floodplain_polygons(
    base_raster,
    polygons: Sequence,
    *,
    labels: Optional[Sequence[str]] = None,
    colours: Optional[Sequence[str]] = None,
    coordinates: str = "local",
    base_is_hillshade: bool = True,
    fill_transparency: float = 85.0,
    outline_transparency: float = 40.0,
    outline_width: str = "0.35p",
    drape_raster=None,
    drape_cmap: str = "viridis",
    drape_vmin: Optional[float] = None,
    drape_vmax: Optional[float] = None,
    drape_transparency: float = 30.0,
    drape_label: Optional[str] = None,
    legend: bool = True,
    style: Optional[MapStyle] = None,
    fig=None,
):
    """Draw polygon sets over a gray hillshade (optionally a coloured drape):
    faint fill, thin semi-transparent outline.

    Args:
        base_raster: Hillshade (``base_is_hillshade=True``, e.g. lsdtt3
            ``_hillshade.tif``) or elevation to illuminate - ``DataArray`` or path.
        polygons: Polygon sets, drawn in order (first underneath): each a
            ``GeoDataFrame`` or a path to a vector file.
        labels: Legend label per set (defaults to "set 1", "set 2", ...).
        colours: GMT colour per set (defaults to ``DEFAULT_COLOURS``).
        coordinates: ``"local"`` / ``"geographic"`` / explicit CRS, for the base
            and every polygon set alike.
        base_is_hillshade: Draw the base flat (``True``) or illuminate it.
        fill_transparency: Percent transparency of the fills (0 opaque, 100
            invisible). High by default so the hillshade shows through.
        outline_transparency: Percent transparency of the outlines. Partly
            transparent by default so a boundary close to another (the
            floodplain against a fan) stays readable.
        outline_width: GMT pen width of the outlines.
        drape_raster: Optional raster coloured over the hillshade under the
            polygons (e.g. lsdtt3 ``_reference_field.tif``, HAND) - ``DataArray``
            or path, on the base's grid. Values outside ``drape_vmin`` ..
            ``drape_vmax`` are clamped to the end colours.
        drape_cmap / drape_vmin / drape_vmax / drape_transparency: its colour
            map, range and transparency (percent).
        drape_label: Colourbar label for the drape (no colourbar when None).
        legend: Add a legend box of the sets.
        style: Styling options (figure size, frame); the base uses a gray CPT.
        fig: Existing ``pygmt.Figure`` to draw onto, or ``None``.

    Returns:
        The ``pygmt.Figure``.
    """
    import pygmt

    style = style or MapStyle()
    n = len(polygons)
    if n == 0:
        raise ValueError("render_floodplain_polygons needs at least one polygon set.")
    labels = list(labels) if labels is not None else [f"set {i + 1}" for i in range(n)]
    colours = list(colours) if colours is not None else list(DEFAULT_COLOURS)
    if len(labels) != n:
        raise ValueError(f"{len(labels)} labels for {n} polygon sets.")
    if len(colours) < n:
        raise ValueError(f"{len(colours)} colours for {n} polygon sets.")

    if drape_raster is None:
        fig = render_hillshade(
            base_raster,
            coordinates=coordinates,
            shade=not base_is_hillshade,
            style=replace(style, cmap="gray"),
            fig=fig,
            decorate=False,
        )
    else:
        from .channels import render_drape
        fig = render_drape(
            base_raster, drape_raster,
            coordinates=coordinates,
            base_is_hillshade=base_is_hillshade,
            transparency=drape_transparency,
            cmap=drape_cmap, vmin=drape_vmin, vmax=drape_vmax,
            out_of_bounds="clamp",
            colorbar=drape_label is not None,
            colorbar_label=drape_label,
            style=style,
            fig=fig,
        )

    for data, colour in zip(polygons, colours):
        gdf, _ = _prepare_overlay(data, coordinates, style)
        if gdf.empty:
            continue
        fig.plot(data=gdf, fill=colour, transparency=fill_transparency)
        fig.plot(data=gdf, pen=f"{outline_width},{colour}",
                 transparency=outline_transparency)

    if legend:
        rows = [f"S 0.15c s 0.28c {c}@{fill_transparency:g} {outline_width},{c} 0.4c {lab}"
                for lab, c in zip(labels, colours)]
        with pygmt.config(FONT_ANNOT_PRIMARY="8p"):
            fig.legend(spec=_io.StringIO("\n".join(rows) + "\n"),
                       position="JTR+jTR+o0.2c", box="+gwhite@20+p0.5p,black")
    return fig
