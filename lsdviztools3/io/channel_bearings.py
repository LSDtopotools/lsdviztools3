"""Channel-bearing and junction-angle summaries from lsdtt3 (by order / lithology).

Readers for the grouped network-geometry outputs of lsdtt3:

* ``<wf>_channel_bearing_rose.csv``          — pre-binned rose counts/lengths per
  group (and optional stratum);
* ``<wf>_channel_bearing_group_summary.csv`` — circular statistics per group;
* ``<wf>_junction_angle_group_summary.csv``  — linear statistics of junction
  angles per group;
* ``<wf>_junction_angles.csv``               — per-junction angles (existing CNA
  output).

* ``<wf>_channel_bearing_comparisons.csv`` / ``<wf>_junction_angle_comparisons.csv``
  — pairwise group comparisons;
* ``<wf>_channel_bearing_group_membership.csv`` — per-link group assignment.

Every column and file name used anywhere in the package is in the
``*_COLUMNS`` / ``SUFFIXES`` constants here — edit them, nothing else, if
lsdtt3 renames a column. Tables are read by header; extra columns are ignored
and appended optional columns (e.g. the junction median CI) may be absent.

Value conventions: ``measurand`` is ``link_bearing`` for bearings and
``alpha``/``beta_min``/``beta_max`` for junctions; high Strahler orders pool as
``"4+"``; ``nan`` is written as text; p-value columns end ``_iid_approx`` and
must be shown with that label and with ``n_blocks`` beside them.

Bearings are compass azimuths: 0 = grid north, clockwise, degrees. Nothing here
computes a statistic; circular means and R-bar come from the summary file.
No PyGMT imports.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from .table import load_table

#: Grouping keys lsdtt3 writes in ``group_key`` / ``stratum_key``.
GROUP_KEYS = ("strahler_order", "basin", "region", "lithology")
#: ``stratum_key`` value meaning "not stratified".
NO_STRATUM = "none"

#: Columns shared by the grouped tables.
GROUP_COLUMNS = {
    "measurand": "measurand",
    "group_key": "group_key",
    "group_value": "group_value",
    "stratum_key": "stratum_key",
    "stratum_value": "stratum_value",
}

#: ``<wf>_channel_bearing_rose.csv``.
ROSE_COLUMNS = {
    **GROUP_COLUMNS,
    "angular_type": "angular_type",
    "bin_centre": "bin_centre_deg",
    "bin_lower": "bin_lower_deg",
    "bin_upper": "bin_upper_deg",
    "count": "count",
    "fraction_count": "fraction_count",
    "length": "length_m",
    "fraction_length": "fraction_length",
}

#: ``<wf>_channel_bearing_group_summary.csv``.
BEARING_SUMMARY_COLUMNS = {
    **GROUP_COLUMNS,
    "angular_type": "angular_type",
    "period": "period_deg",
    "n": "n",
    "n_blocks": "n_blocks",
    "total_length": "total_length_m",
    "mean_count": "circular_mean_deg_count",
    "rbar_count": "mean_resultant_length_count",
    "mean_length": "circular_mean_deg_length",
    "rbar_length": "mean_resultant_length_length",
    "n_eff_length": "n_eff_length",
    "rayleigh_p": "rayleigh_p_iid_approx",
    "fraction_d8": "fraction_on_d8_multiples",
    "meets_min": "meets_min_group_size",
}

#: ``<wf>_junction_angle_group_summary.csv`` (linear stats per measurand).
#: The three ``median_ci_*``/``n_bootstrap`` columns are appended later and
#: are optional.
JUNCTION_SUMMARY_COLUMNS = {
    **GROUP_COLUMNS,
    "n": "n",
    "n_blocks": "n_blocks",
    "mean": "mean_deg",
    "median": "median_deg",
    "sd": "sd_deg",
    "q25": "q25_deg",
    "q75": "q75_deg",
    "iqr": "iqr_deg",
    "meets_min": "meets_min_group_size",
    "median_ci_low": "median_ci_low_deg",
    "median_ci_high": "median_ci_high_deg",
    "n_bootstrap": "n_bootstrap",
}

#: Junction-angle measurands in the summary.
JUNCTION_MEASURANDS = ("alpha", "beta_min", "beta_max")

#: Per-junction table (existing lsdtt-channel-network-analysis output):
#: measurand -> angle column, and the order column the summary groups by
#: (``junction_order``; its counts match the summary's order groups).
JUNCTION_COLUMNS = {
    "alpha": "alpha_deg",
    "beta_1": "beta_1_deg",
    "beta_2": "beta_2_deg",
    "order": "junction_order",
}

#: ``<wf>_channel_bearing_comparisons.csv``.
BEARING_COMPARISON_COLUMNS = {
    "measurand": "measurand",
    "stratum_key": "stratum_key",
    "stratum_value": "stratum_value",
    "angular_type": "angular_type",
    "mode": "comparison_mode",
    "group_a": "group_a",
    "group_b": "group_b",
    "n_a": "n_a",
    "n_b": "n_b",
    "n_blocks_a": "n_blocks_a",
    "n_blocks_b": "n_blocks_b",
    "wasserstein": "circular_wasserstein_deg",
    "rotation_invariant": "rotation_invariant",
    "optimal_rotation": "optimal_rotation_deg",
    "kuiper_V": "kuiper_V",
    "kuiper_p": "kuiper_p_iid_approx",
    "kuiper_p_holm": "kuiper_p_holm",
    "n_family": "n_comparisons_in_family",
    "permutation_p": "permutation_p",
    "note": "note",
}

#: ``<wf>_junction_angle_comparisons.csv``.
JUNCTION_COMPARISON_COLUMNS = {
    "measurand": "measurand",
    "stratum_key": "stratum_key",
    "stratum_value": "stratum_value",
    "mode": "comparison_mode",
    "group_a": "group_a",
    "group_b": "group_b",
    "n_a": "n_a",
    "n_b": "n_b",
    "n_blocks_a": "n_blocks_a",
    "n_blocks_b": "n_blocks_b",
    "ks_D": "ks_D",
    "ks_p": "ks_p_iid_approx",
    "ks_p_holm": "ks_p_holm",
    "n_family": "n_comparisons_in_family",
    "wasserstein": "linear_wasserstein_deg",
    "permutation_p": "permutation_p",
    "note": "note",
}

#: ``<wf>_channel_bearing_group_membership.csv``. An empty ``group_value``
#: means the record was excluded from grouping.
MEMBERSHIP_COLUMNS = {
    "upstream_junction": "upstream_junction",
    "downstream_junction": "downstream_junction",
    "stream_order": "stream_order",
    "length": "total_length_m",
    "bearing": "bearing_deg",
    "axial_bearing": "axial_bearing_deg",
    "group_value": "group_value",
    "stratum_value": "stratum_value",
    "basin_id": "basin_id",
    "block_id": "block_id",
    "region_value": "region_value",
    "region_purity": "region_purity",
}

#: Angular types.
DIRECTIONAL, AXIAL = "directional", "axial"

SUFFIXES = {
    "rose": "_channel_bearing_rose.csv",
    "bearing_summary": "_channel_bearing_group_summary.csv",
    "bearing_comparisons": "_channel_bearing_comparisons.csv",
    "membership": "_channel_bearing_group_membership.csv",
    "junction_summary": "_junction_angle_group_summary.csv",
    "junction_comparisons": "_junction_angle_comparisons.csv",
    "junctions": "_junction_angles.csv",
    "link_vectors": "_channel_link_bearing_vectors.fgb",
    "strahler_channels": "_strahler_channels.csv",
}


def channel_geometry_paths(directory, prefix: str) -> dict:
    """``{kind: Path}`` for the grouped network-geometry outputs of ``prefix``."""
    directory = Path(directory)
    return {k: directory / f"{prefix}{s}" for k, s in SUFFIXES.items()}


def _require(df: pd.DataFrame, cols, what: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"{what} is missing column(s) {missing}; has {list(df.columns)}")


def load_bearing_rose(path) -> pd.DataFrame:
    """Load ``<wf>_channel_bearing_rose.csv``."""
    df = load_table(path)
    c = ROSE_COLUMNS
    _require(df, (c["group_key"], c["group_value"], c["angular_type"], c["bin_centre"],
                  c["bin_lower"], c["bin_upper"], c["count"]), Path(path).name)
    return df


def load_bearing_group_summary(path) -> pd.DataFrame:
    """Load ``<wf>_channel_bearing_group_summary.csv``."""
    df = load_table(path)
    c = BEARING_SUMMARY_COLUMNS
    _require(df, (c["group_key"], c["group_value"], c["angular_type"]), Path(path).name)
    return df


def load_junction_angle_group_summary(path) -> pd.DataFrame:
    """Load ``<wf>_junction_angle_group_summary.csv``."""
    df = load_table(path)
    c = JUNCTION_SUMMARY_COLUMNS
    _require(df, (c["measurand"], c["group_key"], c["group_value"], c["median"],
                  c["q25"], c["q75"]), Path(path).name)
    return df


def load_bearing_comparisons(path) -> pd.DataFrame:
    """Load ``<wf>_channel_bearing_comparisons.csv``."""
    df = load_table(path)
    c = BEARING_COMPARISON_COLUMNS
    _require(df, (c["group_a"], c["group_b"], c["n_blocks_a"], c["n_blocks_b"]),
             Path(path).name)
    return df


def load_junction_angle_comparisons(path) -> pd.DataFrame:
    """Load ``<wf>_junction_angle_comparisons.csv``."""
    df = load_table(path)
    c = JUNCTION_COMPARISON_COLUMNS
    _require(df, (c["group_a"], c["group_b"], c["n_blocks_a"], c["n_blocks_b"]),
             Path(path).name)
    return df


def load_junction_angles_table(path) -> pd.DataFrame:
    """Load the per-junction ``<wf>_junction_angles.csv`` (CNA output)."""
    df = load_table(path)
    _require(df, (JUNCTION_COLUMNS["alpha"], JUNCTION_COLUMNS["order"]), Path(path).name)
    return df


def order_group_mask(orders: pd.Series, group_value) -> pd.Series:
    """Rows whose order belongs to an lsdtt3 order group (``"3"`` or pooled ``"4+"``)."""
    o = pd.to_numeric(orders, errors="coerce")
    g = str(group_value).strip()
    if g.endswith("+"):
        return o >= float(g[:-1])
    return o == float(g)


def load_channel_geometry_outputs(directory, prefix: str) -> dict:
    """Load whichever grouped bearing/junction tables exist for ``prefix``.

    Returns ``{kind: DataFrame or None, ..., "paths": {...}}`` with kinds
    ``rose``, ``bearing_summary``, ``bearing_comparisons``, ``junction_summary``,
    ``junction_comparisons``, ``junctions``.
    """
    paths = channel_geometry_paths(directory, prefix)
    loaders = {
        "rose": load_bearing_rose,
        "bearing_summary": load_bearing_group_summary,
        "bearing_comparisons": load_bearing_comparisons,
        "junction_summary": load_junction_angle_group_summary,
        "junction_comparisons": load_junction_angle_comparisons,
        "junctions": load_junction_angles_table,
    }
    out = {k: (f(paths[k]) if paths[k].exists() else None) for k, f in loaders.items()}
    out["paths"] = paths
    return out


def select_groups(df: pd.DataFrame, *, group_key: str, measurand: Optional[str] = None,
                  stratum_key: str = NO_STRATUM, stratum_value=None,
                  angular_type: Optional[str] = None) -> pd.DataFrame:
    """Rows of a grouped table for one grouping (and optional stratum/measurand).

    ``stratum_key="none"`` keeps unstratified rows (a missing/blank stratum
    column counts as unstratified).
    """
    c = GROUP_COLUMNS
    out = df[df[c["group_key"]].astype(str) == str(group_key)]
    if measurand is not None and c["measurand"] in out.columns:
        out = out[out[c["measurand"]].astype(str) == str(measurand)]
    if c["stratum_key"] in out.columns:
        sk = out[c["stratum_key"]].fillna(NO_STRATUM).astype(str)
        out = out[sk == str(stratum_key)]
        if stratum_value is not None and c["stratum_value"] in out.columns:
            out = out[out[c["stratum_value"]].astype(str) == str(stratum_value)]
    if angular_type is not None and "angular_type" in out.columns:
        out = out[out["angular_type"].astype(str) == angular_type]
    return out
