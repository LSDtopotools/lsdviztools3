"""Network-coverage measurement for lsdtt3 stream/burn outputs.

Quantifies how well a derived network (e.g. a D8 ``*_connected_network_lines``)
follows a set of reference reaches (e.g. ``*_accepted_clean``): each reach is
sampled at a fixed spacing and a sample counts as covered when it lies within a
tolerance of the network linework. This is the QC metric the stream-burning work
is judged on; it was previously copy-pasted across ad-hoc scripts. No PyGMT here.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
from shapely.ops import unary_union

from .vector import load_vector


def _as_gdf(data, layer=None) -> gpd.GeoDataFrame:
    """Return a ``GeoDataFrame`` from either a path or an existing frame."""
    if isinstance(data, gpd.GeoDataFrame):
        return data
    return load_vector(data, layer=layer)


def _coverage_counts(geom, network_union, sample_step: float, tol: float):
    """Count (covered, total) samples along one line geometry.

    Samples ``geom`` every ``sample_step`` m (always including the start); a
    sample is covered when within ``tol`` of ``network_union``.
    """
    length = geom.length
    total = 0
    covered = 0
    for s in np.arange(0.0, length, sample_step):
        total += 1
        if geom.interpolate(s).distance(network_union) <= tol:
            covered += 1
    if total == 0:  # degenerate / zero-length reach: sample the single point
        total = 1
        covered = int(geom.interpolate(0.0).distance(network_union) <= tol)
    return covered, total


def network_coverage(reaches, network, sample_step: float = 2.0, tol: float = 3.0) -> dict:
    """Measure how well ``network`` covers the reference ``reaches``.

    Each reach is sampled every ``sample_step`` m; a sample is "covered" when it
    lies within ``tol`` of the (unioned) ``network`` linework. The network union
    is built once.

    Args:
        reaches: Reference reaches — a file path (``.fgb`` / ``.gpkg`` / ``.shp``)
            or a ``GeoDataFrame``.
        network: Derived network — a file path or a ``GeoDataFrame``.
        sample_step: Sampling spacing along each reach, in CRS units (m).
        tol: Coverage tolerance, in CRS units (m).

    Returns:
        ``dict`` with ``covered`` and ``total`` sample counts, ``fraction``
        (0..1), and ``length_m`` / ``n_features`` describing the network.
    """
    reaches_gdf = _as_gdf(reaches)
    network_gdf = _as_gdf(network)
    network_union = unary_union(network_gdf.geometry.values)

    covered = 0
    total = 0
    for geom in reaches_gdf.geometry:
        c, t = _coverage_counts(geom, network_union, sample_step, tol)
        covered += c
        total += t

    return {
        "covered": covered,
        "total": total,
        "fraction": (covered / total) if total else 0.0,
        "length_m": float(network_gdf.geometry.length.sum()),
        "n_features": int(len(network_gdf)),
    }


def per_reach_coverage(
    reaches,
    network,
    *,
    sample_step: float = 2.0,
    tol: float = 3.0,
    attrs=("strahler", "edge_id"),
) -> gpd.GeoDataFrame:
    """Per-reach coverage of ``reaches`` by ``network``.

    Same sampling kernel as :func:`network_coverage`, reported per reach.

    Args:
        reaches: Reference reaches — a file path or a ``GeoDataFrame``.
        network: Derived network — a file path or a ``GeoDataFrame``.
        sample_step: Sampling spacing along each reach, in CRS units (m).
        tol: Coverage tolerance, in CRS units (m).
        attrs: Attribute columns to carry through if present (others ignored).

    Returns:
        A copy of ``reaches`` with an added float ``coverage`` column (0..100 %),
        retaining its geometry and whichever of ``attrs`` exist.
    """
    reaches_gdf = _as_gdf(reaches)
    network_gdf = _as_gdf(network)
    network_union = unary_union(network_gdf.geometry.values)

    coverage = np.empty(len(reaches_gdf), dtype=float)
    for pos, geom in enumerate(reaches_gdf.geometry):
        c, t = _coverage_counts(geom, network_union, sample_step, tol)
        coverage[pos] = 100.0 * c / t

    keep = [a for a in attrs if a in reaches_gdf.columns]
    out = reaches_gdf[keep + [reaches_gdf.geometry.name]].copy()
    out["coverage"] = coverage
    return out
