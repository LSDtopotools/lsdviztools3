"""Derived-grid preparation (hillshade / gradient).

This is still IO-level data prep: it produces ``xarray.DataArray`` grids ready
for rendering. NOTE: ``make_hillshade`` may call ``pygmt.grdgradient`` to
compute the illumination grid — this is the one permitted PyGMT touch-point in
the core, and it is data prep rather than rendering.

This module is a deferred extension point in v1.
"""

from __future__ import annotations

import xarray as xr

_DEFERRED_MSG = (
    "lsdviztools3.io.grids is deferred in v1 and not implemented yet."
)


def make_hillshade(
    da: xr.DataArray,
    *,
    azimuth: float = 315.0,
    altitude: float = 45.0,
    normalize: bool = True,
) -> xr.DataArray:
    """Compute a hillshade/illumination grid from an elevation grid.

    Uses gradient-based illumination (may delegate to ``pygmt.grdgradient``).

    Args:
        da: Elevation grid as an ``xarray.DataArray``.
        azimuth: Light source azimuth in degrees (default 315).
        altitude: Light source altitude in degrees (default 45).
        normalize: Whether to normalise the gradient output.

    Returns:
        Hillshade grid as an ``xarray.DataArray``.
    """
    raise NotImplementedError(_DEFERRED_MSG)
