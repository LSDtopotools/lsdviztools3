"""STI-specific map rendering helpers."""

from __future__ import annotations

from typing import Optional

from .channels import render_drape
from .style import MapStyle


def render_sti_p95_drape(
    hillshade_raster,
    sti_p95_raster,
    *,
    coordinates: str = "local",
    transparency: float = 35.0,
    cmap: str = "YlOrRd",
    colorbar: bool = False,
    style: Optional[MapStyle] = None,
    fig=None,
):
    """Drape STI-P95 values over a hillshade.

    Intended for lsdtt3 outputs such as ``*_sti_p95.tif`` or
    ``*_sti_proxy_p95.tif``. Lower values map to yellow and higher values map to
    red by default, and the colourbar is off by default.

    Args:
        hillshade_raster: Precomputed hillshade raster (``*_hillshade.tif``), as
            an xarray ``DataArray`` or a path.
        sti_p95_raster: STI threshold raster (``*_sti_p95.tif`` or
            ``*_sti_proxy_p95.tif``), as an xarray ``DataArray`` or a path.
        coordinates: Coordinate mode (``local`` by default) or explicit CRS.
        transparency: Drape transparency in percent (0 opaque, 100 invisible).
        cmap: Colour map used for STI values (default ``YlOrRd``).
        colorbar: Add a colourbar. Default ``False``.
        style: Optional ``MapStyle`` overrides.
        fig: Existing ``pygmt.Figure`` for composition, or ``None``.

    Returns:
        The ``pygmt.Figure``.
    """
    return render_drape(
        base_raster=hillshade_raster,
        drape_raster=sti_p95_raster,
        coordinates=coordinates,
        base_is_hillshade=True,
        log_color=False,
        transparency=transparency,
        cmap=cmap,
        colorbar=colorbar,
        style=style,
        fig=fig,
    )
