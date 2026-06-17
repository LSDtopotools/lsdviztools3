"""Tabular IO for lsdtt3 point/table CSVs.

Reads CSV data into ``pandas.DataFrame``. lsdtt3 point CSVs use EPSG:4326
(lat/long) coordinates; chi fit tables (nodes/segments) are plain attribute
CSVs. No PyGMT imports here.
"""

from __future__ import annotations

import pandas as pd

#: Common spellings for latitude/longitude columns in lsdtt3 CSVs.
_LAT_NAMES = ("latitude", "lat", "y")
_LON_NAMES = ("longitude", "long", "lon", "lng", "x")


def load_table(path) -> pd.DataFrame:
    """Load a CSV table as a ``pandas.DataFrame``.

    lsdtt3 point CSVs use EPSG:4326 (lat/long) coordinates; attribute tables
    (e.g. chi fit nodes/segments) carry no geometry.

    Args:
        path: Path to the CSV file.

    Returns:
        The table as a ``pandas.DataFrame``.
    """
    return pd.read_csv(path)


def latlon_columns(df: pd.DataFrame) -> tuple:
    """Detect the latitude and longitude columns in a DataFrame.

    Matching is case-insensitive against common lsdtt3 spellings.

    Args:
        df: Input table.

    Returns:
        ``(lat_col, lon_col)`` column-name tuple.

    Raises:
        KeyError: if either column cannot be identified.
    """
    lower = {c.lower(): c for c in df.columns}
    lat = next((lower[n] for n in _LAT_NAMES if n in lower), None)
    lon = next((lower[n] for n in _LON_NAMES if n in lower), None)
    if lat is None or lon is None:
        raise KeyError(
            "Could not find latitude/longitude columns in table; "
            f"columns are {list(df.columns)}."
        )
    return lat, lon
