"""Chi-analysis IO for lsdtt3 outputs.

Readers for the ``lsdtt-chi-analysis`` products used by the chi profile plots:

* ``*_chi_fit_nodes_mn_*.csv``    — per-node fitted data (chi, observed and
  fitted elevation, ``basin_index``, ``branch_id``, …).
* ``*_chi_fit_segments_mn_*.csv`` — per-segment piecewise-linear fits, carrying
  ``k_sn`` (dz/dchi) plus the chi/elevation endpoints.
* ``*_chi_data_maps.fgb`` / ``.csv`` — per-node chi point network (chi,
  elevation, flow_distance, drainage_area, ``basin_key``, ``source_key``).

The fit tables are self-consistent: both group by ``basin_index`` and
``branch_id``, and ``k_sn`` (in the segments) is attached to nodes by locating
the segment whose chi interval contains each node. No PyGMT imports here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .table import load_table

#: Columns expected in a chi fit-nodes table.
FIT_NODE_COLUMNS = (
    "basin_index", "outlet_flat", "flat", "receiver_flat", "chi",
    "observed_elevation", "fitted_elevation", "residual_elevation",
    "abs_residual_elevation", "branch_id", "is_main_stem", "is_source",
    "is_junction", "is_split_point",
)

#: Columns expected in a chi fit-segments table.
FIT_SEGMENT_COLUMNS = (
    "basin_index", "outlet_flat", "segment_id", "parent_segment_id",
    "branch_id", "downstream_flat", "upstream_flat", "downstream_chi",
    "upstream_chi", "downstream_elevation", "upstream_elevation", "k_sn",
    "max_abs_residual", "n_nodes",
)


def load_chi_fit_nodes(path) -> pd.DataFrame:
    """Load a ``*_chi_fit_nodes_mn_*.csv`` table.

    Args:
        path: Path to the fit-nodes CSV.

    Returns:
        The per-node fitted data as a ``pandas.DataFrame``.
    """
    return load_table(path)


def load_chi_fit_segments(path) -> pd.DataFrame:
    """Load a ``*_chi_fit_segments_mn_*.csv`` table.

    Args:
        path: Path to the fit-segments CSV.

    Returns:
        The per-segment fitted data as a ``pandas.DataFrame``.
    """
    return load_table(path)


def assign_node_ksn(nodes: pd.DataFrame, segments: pd.DataFrame) -> pd.DataFrame:
    """Attach a per-node ``k_sn`` from the segment that contains each node.

    For each ``(basin_index, branch_id)`` group, every node is matched to the
    fitted segment whose ``[downstream_chi, upstream_chi]`` interval contains
    the node's chi. Nodes that fall just outside any interval (e.g. boundary
    rounding) are assigned the nearest segment in that branch.

    Args:
        nodes: Fit-nodes table (must have ``basin_index``, ``branch_id``, ``chi``).
        segments: Fit-segments table (must have the same grouping keys plus
            ``downstream_chi``, ``upstream_chi``, ``k_sn``).

    Returns:
        A copy of ``nodes`` with an added ``k_sn`` column (``NaN`` where the
        branch has no segments).
    """
    out = nodes.copy()
    out["k_sn"] = np.nan
    if segments.empty or nodes.empty:
        return out

    seg_by_group = {
        key: grp for key, grp in segments.groupby(["basin_index", "branch_id"])
    }

    for key, node_grp in out.groupby(["basin_index", "branch_id"]):
        seg_grp = seg_by_group.get(key)
        if seg_grp is None or seg_grp.empty:
            continue

        lo = seg_grp["downstream_chi"].to_numpy(dtype=float)
        hi = seg_grp["upstream_chi"].to_numpy(dtype=float)
        ksn = seg_grp["k_sn"].to_numpy(dtype=float)
        # Guard against endpoints stored in either order.
        lo, hi = np.minimum(lo, hi), np.maximum(lo, hi)

        chi = node_grp["chi"].to_numpy(dtype=float)
        # Interval containment: segment s contains chi if lo[s] <= chi <= hi[s].
        contained = (chi[:, None] >= lo[None, :]) & (chi[:, None] <= hi[None, :])
        seg_idx = np.where(
            contained.any(axis=1),
            contained.argmax(axis=1),
            # Fallback: nearest segment by interval midpoint.
            np.abs(chi[:, None] - ((lo + hi) / 2.0)[None, :]).argmin(axis=1),
        )
        out.loc[node_grp.index, "k_sn"] = ksn[seg_idx]

    return out


def load_chi_fit(nodes_path, segments_path) -> pd.DataFrame:
    """Load fit nodes and attach segment ``k_sn`` in one step.

    Convenience wrapper around :func:`load_chi_fit_nodes`,
    :func:`load_chi_fit_segments`, and :func:`assign_node_ksn`.

    Args:
        nodes_path: Path to the ``*_chi_fit_nodes_mn_*.csv``.
        segments_path: Path to the ``*_chi_fit_segments_mn_*.csv``.

    Returns:
        Per-node fitted data with a ``k_sn`` column.
    """
    nodes = load_chi_fit_nodes(nodes_path)
    segments = load_chi_fit_segments(segments_path)
    return assign_node_ksn(nodes, segments)


def load_chi_segment_lines(path) -> pd.DataFrame:
    """Load fitted chi segments from a ``*_chi_fit_segment_lines_mn_*.fgb``.

    The segment-lines vector encodes the fitted chi-elevation profile in its
    attributes (``downstream_chi``/``upstream_chi`` and the matching elevations,
    plus ``k_sn`` per segment) — the geometry itself is in map space, so only
    the attribute table is returned. The columns overlap the CSV
    ``*_chi_fit_segments_mn_*.csv`` table, so the chi plotters accept either.

    Args:
        path: Path to the segment-lines ``.fgb`` (or ``.gpkg``).

    Returns:
        The segment attribute table as a ``pandas.DataFrame`` (geometry dropped).
    """
    from .vector import load_vector

    gdf = load_vector(path)
    return pd.DataFrame(gdf.drop(columns=gdf.geometry.name))


def load_chi_points(path) -> pd.DataFrame:
    """Load a chi data-map point network (``.fgb``/``.gpkg`` or ``.csv``).

    The per-node chi network carries ``chi``, ``elevation``, ``flow_distance``,
    ``drainage_area``, ``source_key`` and ``basin_key``. For vector formats the
    geometry's x/y are exposed as ``x``/``y`` columns and the geometry column is
    dropped (this stays a plain ``DataFrame`` with no PyGMT dependency).

    Args:
        path: Path to the chi data-map file.

    Returns:
        The per-node chi network as a ``pandas.DataFrame``.
    """
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return load_table(path)

    # Vector formats: read with geopandas, then flatten geometry to x/y.
    from .vector import load_vector

    gdf = load_vector(path)
    df = pd.DataFrame(gdf.drop(columns=gdf.geometry.name))
    pts = gdf.geometry.representative_point()
    df["x"] = pts.x.to_numpy()
    df["y"] = pts.y.to_numpy()
    return df
