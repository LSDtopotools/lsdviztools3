"""Swath-profile IO for lsdtt3 outputs.

Reader for the ``lsdtt-swath-profile`` per-bin statistics CSV
(``*_swath_profile.csv``) and a helper to pull a single main-stem baseline
polyline out of a channel/chi line layer for use as the swath axis. No PyGMT
imports here.

The swath CSV columns are::

    along_axis_centre, n_pixels, mean, p0, p5, p10, ... p95, p100

i.e. the bin-centre distance along the baseline (m), the pixel count, the mean
value, and a 21-level percentile CDF (p0 = min, p50 = median, p100 = max).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .table import load_table

#: The 21 fixed percentile levels lsdtt3 writes per swath bin.
SWATH_PERCENTILE_LEVELS = tuple(range(0, 101, 5))

#: Percentile column names in order (``p0`` .. ``p100``).
SWATH_PERCENTILE_COLUMNS = tuple(f"p{lvl}" for lvl in SWATH_PERCENTILE_LEVELS)


def load_swath(path) -> pd.DataFrame:
    """Load a ``*_swath_profile.csv`` per-bin statistics table.

    Empty bins (``n_pixels == 0``) carry ``NaN`` percentiles/mean.

    Args:
        path: Path to the swath profile CSV.

    Returns:
        The per-bin statistics as a ``pandas.DataFrame``.
    """
    return load_table(path)


def swath_percentile_columns(df: pd.DataFrame) -> list:
    """Return the percentile column names present in ``df``, in ascending order."""
    return [c for c in SWATH_PERCENTILE_COLUMNS if c in df.columns]


def extract_main_channel_baseline(
    lines_path,
    out_path,
    *,
    basin: int = 0,
    basin_col: str = "basin_key",
):
    """Write the longest channel line of one basin as a single-feature baseline.

    Reads a channel/chi line layer (e.g. ``*_chi_data_map_lines.fgb``, one
    primary-donor longest-path line per source), keeps the features for the
    requested basin, and writes the single longest line — the main stem — to
    ``out_path`` as a one-feature FlatGeobuf in the source CRS. That line is the
    swath baseline fed to ``lsdtt-swath-profile``.

    Args:
        lines_path: Channel/chi line layer (``.fgb``/``.gpkg``/``.shp``).
        out_path: Output path for the single-line baseline (suffix picks driver).
        basin: Basin id to select (default 0 — the largest, with n_largest=1).
        basin_col: Attribute naming the basin id (``basin_key`` for chi lines).

    Returns:
        ``Path`` to the written baseline file.
    """
    from .vector import load_vector

    gdf = load_vector(lines_path)
    if basin_col in gdf.columns:
        sel = gdf[gdf[basin_col] == basin]
        if sel.empty:
            sel = gdf
    else:
        sel = gdf
    if sel.empty:
        raise ValueError(f"No channel lines found in {lines_path}.")

    # Longest line = main stem (headwater source to outlet).
    longest = sel.loc[sel.geometry.length.idxmax()]
    baseline = sel.loc[[longest.name]].copy()

    out_path = Path(out_path)
    baseline.to_file(out_path)
    return out_path
