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
    "is_junction", "is_split_point", "source_key",
)

#: Columns expected in a chi fit-segments table.
FIT_SEGMENT_COLUMNS = (
    "basin_index", "outlet_flat", "segment_id", "parent_segment_id",
    "branch_id", "downstream_flat", "upstream_flat", "downstream_chi",
    "upstream_chi", "downstream_elevation", "upstream_elevation", "k_sn",
    "max_abs_residual", "n_nodes", "source_key",
)

#: Columns of ``*_chi_tributaries_mn_*.csv`` (one row per tributary = fit branch).
TRIBUTARY_COLUMNS = (
    "basin_index", "basin_key", "source_key", "branch_id", "parent_branch_id",
    "branch_order", "head_flat", "junction_flat", "head_x", "head_y",
    "junction_x", "junction_y", "junction_drainage_area_m2",
    "basin_drainage_area_m2", "area_fraction", "area_rank_in_basin", "length_m",
    "chi_span", "n_segments", "ksn_min", "ksn_median", "ksn_max",
    "ksn_max_segment_id", "ksn_max_chi", "n_descendants",
)

#: Columns of ``*_chi_knickpoints_mn_*.csv`` (one row per k_sn segment boundary).
KNICKPOINT_COLUMNS = (
    "basin_index", "basin_key", "source_key", "branch_id", "upstream_segment_id",
    "downstream_segment_id", "flat", "x", "y", "chi", "elevation",
    "flow_distance_m", "drainage_area_m2", "ksn_upstream", "ksn_downstream",
    "delta_ksn", "log_ratio_ksn", "at_junction", "chi_from_branch_head",
    "chi_normalised",
)

#: Columns of ``*_chi_knickpoint_thresholds_mn_*.csv``.
THRESHOLD_COLUMNS = ("measure", "method", "threshold", "n_above", "n_total")


def _load_with_columns(path, columns, what):
    df = load_table(path)
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: not an lsdtt-chi-analysis {what} table "
                         f"(missing {missing}); rerun with a current lsdtt3.")
    return df


def load_chi_tributaries(path) -> pd.DataFrame:
    """Load ``*_chi_tributaries_mn_*.csv`` (``write_chi_tributary_table``).

    One row per tributary (a branch of the chi fit, head to junction), keyed
    by ``source_key``: drainage area where it joins, its fraction of the
    basin, rank by area, parent, order, k_sn summary, number of descendants.
    """
    return _load_with_columns(path, TRIBUTARY_COLUMNS, "tributary")


def load_chi_knickpoints(path) -> pd.DataFrame:
    """Load ``*_chi_knickpoints_mn_*.csv`` (``write_chi_knickpoint_table``).

    One row per node where two fitted k_sn segments meet, with ``delta_ksn``
    (upstream minus downstream), ``log_ratio_ksn`` and ``at_junction`` (1 at a
    tributary mouth).
    """
    return _load_with_columns(path, KNICKPOINT_COLUMNS, "knickpoint")


def load_chi_knickpoint_thresholds(path) -> pd.DataFrame:
    """Load ``*_chi_knickpoint_thresholds_mn_*.csv``: suggested thresholds
    (p90, p95, median_plus_3mad, otsu for abs_delta_ksn and abs_log_ratio_ksn;
    q3_plus_3iqr for ksn) computed by lsdtt3."""
    return _load_with_columns(path, THRESHOLD_COLUMNS, "knickpoint-threshold")


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


#: Removal reasons, in the order they are tested (the first that applies is
#: recorded).
REMOVAL_REASONS = ("ksn_outlier", "small_area", "small_fraction", "rank", "drains_to_removed")


def select_tributaries(
    tributaries: pd.DataFrame,
    *,
    ksn_max=None,
    min_junction_area_m2=None,
    min_area_fraction=None,
    top_n=None,
    remove_descendants: bool = True,
    thresholds: pd.DataFrame = None,
):
    """Choose tributaries to remove from chi plots (selection only, no fitting).

    Works on the lsdtt3 tributary table (:func:`load_chi_tributaries`). A
    tributary is a branch of the chi fit, from its head to where it joins its
    parent; main stems (``branch_id`` 0) are never removed. Rules, each
    ``None`` to skip, applied to every other tributary:

    * ``ksn_max``: remove a tributary with any segment steeper than this
      (``ksn_max`` column), e.g. a link routed over a dam onto a steep face.
      A number (k_sn), or the name of a suggestion in the lsdtt3 thresholds
      table, ``"q3_plus_3iqr"`` (needs ``thresholds``);
    * ``min_junction_area_m2``: remove one whose drainage area just before it
      joins its parent is below this;
    * ``min_area_fraction``: remove one whose area at the junction is below
      this fraction of its basin's;
    * ``top_n``: keep only the ``top_n`` largest tributaries of each basin by
      area at the junction (``area_rank_in_basin``);
    * ``remove_descendants``: also remove every tributary that drains into a
      removed one, at any depth (reason ``drains_to_removed``).

    Returns:
        ``(kept, removed)``: the kept rows, and the removed rows with a
        ``reason`` column (the first rule that removed it, in
        :data:`REMOVAL_REASONS` order).
    """
    t = tributaries.copy()
    if isinstance(ksn_max, str):
        if thresholds is None:
            raise ValueError(f"ksn_max={ksn_max!r} needs the thresholds table "
                             "(load_chi_knickpoint_thresholds).")
        name = ksn_max.removeprefix("ksn_")
        hit = thresholds[(thresholds["measure"] == "ksn") & (thresholds["method"] == name)]
        if hit.empty or not np.isfinite(hit["threshold"].iloc[0]):
            raise ValueError(f"no k_sn threshold named {ksn_max!r} in the thresholds table")
        ksn_max = float(hit["threshold"].iloc[0])
    trib = t["branch_id"] != 0
    reason = pd.Series("", index=t.index, dtype=object)

    def mark(mask, why):
        m = mask & trib & (reason == "")
        reason[m] = why

    if ksn_max is not None:
        mark(t["ksn_max"] > float(ksn_max), "ksn_outlier")
    if min_junction_area_m2 is not None:
        mark(t["junction_drainage_area_m2"] < float(min_junction_area_m2), "small_area")
    if min_area_fraction is not None:
        mark(t["area_fraction"] < float(min_area_fraction), "small_fraction")
    if top_n is not None:
        mark(t["area_rank_in_basin"] > int(top_n), "rank")
    if remove_descendants:
        parent = {(b, br): p for b, br, p in zip(t["basin_index"], t["branch_id"], t["parent_branch_id"])}
        removed_keys = {(b, br) for b, br, r in zip(t["basin_index"], t["branch_id"], reason) if r}
        for idx, b, br in zip(t.index, t["basin_index"], t["branch_id"]):
            if reason[idx] or br == 0:
                continue
            p = parent.get((b, br), -1)
            for _ in range(len(t) + 1):
                if p is None or p < 0:
                    break
                if (b, p) in removed_keys:
                    reason[idx] = "drains_to_removed"
                    break
                p = parent.get((b, p), -1)
    t["reason"] = reason
    return t[reason == ""].drop(columns="reason"), t[reason != ""].copy()


def removed_mask(frame: pd.DataFrame, removed: pd.DataFrame) -> np.ndarray:
    """True for rows of a chi fit nodes / segments frame (or segment-line
    attributes) that belong to a removed tributary, matched on
    ``(basin_index, branch_id)``."""
    if removed is None or len(removed) == 0:
        return np.zeros(len(frame), dtype=bool)
    keys = set(zip(removed["basin_index"].astype(int), removed["branch_id"].astype(int)))
    return np.fromiter(((int(b), int(br)) in keys
                        for b, br in zip(frame["basin_index"], frame["branch_id"])),
                       dtype=bool, count=len(frame))
