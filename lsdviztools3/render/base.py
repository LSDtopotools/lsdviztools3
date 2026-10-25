"""Figure lifecycle and display helpers for the PyGMT render layer.

Figure construction (with publication fonts), automatic projection inference,
display-unit handling (km for projected CRS, with a ``CRS is EPSG:XXXXX`` title),
and DPI/transparency-aware saving. A raster in an lsdtt3 oblique frame is a
projected raster too: it is drawn in frame km with along-/across-feature axis
labels, an ``Oblique Mercator, azimuth ...`` title and, unless switched off in
``MapStyle``, the reference easting/northing grid (``graticule``).
"""

from __future__ import annotations

from typing import Optional

from ..frames import frame_azimuth, is_oblique_frame
from .style import MapStyle

#: Metres per kilometre — projected axes are shown in km to avoid 6-7 digit labels.
KM = 1000.0


def new_figure(style: Optional[MapStyle] = None):
    """Create a new ``pygmt.Figure`` with publication font settings applied.

    Args:
        style: Styling options. Defaults applied if ``None``.

    Returns:
        A new ``pygmt.Figure``.
    """
    import pygmt

    style = style or MapStyle()
    fig = pygmt.Figure()
    pygmt.config(
        FONT_ANNOT_PRIMARY=style.font_annotation,
        FONT_LABEL=style.font_label,
        FONT_TITLE=style.font_title,
        COLOR_NAN="white",
    )
    return fig


def is_projected(crs) -> bool:
    """True if the CRS is a projected (metric) CRS."""
    return crs is not None and bool(getattr(crs, "is_projected", False))


def _is_geographic(crs) -> bool:
    """True if the CRS is geographic (lat/long) or unknown."""
    if crs is None:
        return True
    try:
        return bool(crs.is_geographic)
    except AttributeError:
        return True


def auto_projection(extent, crs, width: str = "12c") -> str:
    """Infer a GMT projection string from extent and CRS.

    Geographic CRS -> Mercator (``M<width>``). Projected CRS -> a Cartesian
    linear projection with equal x/y scaling (``X<width>/0``).
    """
    if _is_geographic(crs):
        return f"M{width}"
    return f"X{width}/0"


def auto_frame(extent):
    """Fallback frame: GMT-chosen tick intervals with a boxed frame."""
    return "af"


def grid_to_display(grid, crs):
    """Return a grid and extent in display units.

    Projected grids are rescaled to km (coordinates / 1000); geographic grids
    are left in degrees.

    Returns:
        ``(grid, [west, east, south, north])`` in display units.
    """
    if is_projected(crs):
        grid = grid.assign_coords(x=grid.x / KM, y=grid.y / KM)
    extent = [
        float(grid.x.min()),
        float(grid.x.max()),
        float(grid.y.min()),
        float(grid.y.max()),
    ]
    return grid, extent


def gdf_to_display(gdf, crs):
    """Scale a GeoDataFrame to display units (km for projected, else unchanged).

    Scaling is about the origin (0, 0) so it matches ``grid_to_display`` exactly,
    keeping vector overlays aligned with the base raster.
    """
    if is_projected(crs):
        scaled = gdf.geometry.scale(1 / KM, 1 / KM, origin=(0, 0))
        return gdf.set_geometry(scaled)
    return gdf


def display_frame(crs, style: MapStyle):
    """Build a publication frame: labelled axes, and a CRS title if projected.

    Projected: axes in km with ``Easting (km)`` / ``Northing (km)`` labels;
    an oblique frame: ``Along-feature (km)`` / ``Across-feature (km)`` (the
    feature runs along +x). Geographic: degree axes labelled Longitude /
    Latitude.
    """
    if style.frame is not None:
        return style.frame
    if getattr(style, "clean", False):
        # Plain border only: no annotations, ticks, or axis labels.
        return ["wsen"]
    if is_oblique_frame(crs):
        return [
            "WSne",
            "xaf+lAlong-feature (km)",
            "yaf+lAcross-feature (km)",
        ]
    if is_projected(crs):
        return [
            "WSne",
            "xaf+lEasting (km)",
            "yaf+lNorthing (km)",
        ]
    return ["WSne", "xaf+lLongitude", "yaf+lLatitude"]


def crs_label_text(crs) -> str:
    """The CRS label: ``Oblique Mercator, azimuth 30.0@.`` for a frame (``@.``
    is GMT's degree sign), ``CRS is EPSG:XXXXX``, or ``CRS is (projected)``."""
    azimuth = frame_azimuth(crs)
    if azimuth is not None:
        return f"Oblique Mercator, azimuth {azimuth:.1f}@."
    epsg = crs.to_epsg()
    return f"CRS is EPSG:{epsg}" if epsg else "CRS is (projected)"


def add_crs_label(fig, extent, crs):
    """Add a right-aligned CRS label just above the map frame."""
    if not is_projected(crs):
        return fig
    text = crs_label_text(crs)
    fig.text(
        x=extent[1],
        y=extent[3],
        text=text,
        justify="BR",
        offset="0c/0.12c",
        no_clip=True,
        font="9p,Helvetica,black",
    )
    return fig


def add_crs_annotations(fig, extent, crs, style: Optional[MapStyle] = None):
    """The CRS label and, over an oblique frame, the reference graticule."""
    from .graticule import add_reference_graticule

    add_reference_graticule(fig, extent, crs, style)
    return add_crs_label(fig, extent, crs)


def set_grid_gtype(grid, crs):
    """Tell GMT whether a grid is geographic (1) or Cartesian/projected (0)."""
    try:
        grid.gmt.gtype = 0 if is_projected(crs) else 1
    except (AttributeError, ValueError):
        pass
    return grid


def save_figure(fig, path, *, style: Optional[MapStyle] = None, high_res: bool = False):
    """Save a figure at publication DPI with a transparent background.

    Args:
        fig: The ``pygmt.Figure``.
        path: Output path (e.g. PNG).
        style: Styling options (DPI and transparency come from here).
        high_res: Use ``style.high_res_dpi`` (600) instead of ``style.dpi`` (300).
    """
    style = style or MapStyle()
    dpi = style.high_res_dpi if high_res else style.dpi
    fig.savefig(path, dpi=dpi, transparent=style.transparent)
    return fig


def finalize(fig, savepath: Optional[str] = None, *, style: Optional[MapStyle] = None):
    """Finalise a figure: save to ``savepath`` if given (at publication DPI)."""
    if savepath:
        save_figure(fig, savepath, style=style)
    return fig
