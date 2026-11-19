"""Slope-area IO for ``lsdtt-slope-area`` outputs.

Readers only; nothing here computes a gradient, a concavity or a bin. Every
number comes from the files lsdtt3 wrote. No PyGMT imports here.

Files of one run (``<prefix>`` is the run's ``write_fname``):

* binned tables, one row per non-empty log-area bin, columns
  ``drainage_area,d05_slope,d25_slope,d50_slope,d75_slope,d95_slope,number_of_pixels``
  (``drainage_area`` is the geometric bin centre, m^2):

  - ``<prefix>_lambda_<L>m.csv`` -- polyfit gradient at lambda ``L`` m
    (``7_5`` stands for 7.5);
  - ``<prefix>_d8window_<W>m.csv`` -- D8 flow-path gradient over a ``W`` m
    window;

* ``<prefix>_low_concavity_points.csv`` -- one row per channel pixel
  (``flat,x,y,drainage_area_m2,slope,theta,r2,n_pixels,log10_area_range,
  n_junctions_in_window,low_concavity``; ``low_concavity`` is 1, 0, or -1 for
  no valid fit);
* ``<prefix>_low_concavity_segments.fgb`` -- runs of channel pixels with one
  ``low_concavity`` value;
* ``<prefix>_low_concavity_theta.tif`` / ``_low_concavity_flag.tif`` -- the
  same per-pixel values as rasters;
* ``<prefix>_ingested_params.param`` -- the ``key: value`` parameters of the
  run.

:func:`load_slope_area` moved here from :mod:`lsdviztools3.io.ridgeline`,
which re-exports it unchanged.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

#: Columns of the ``lsdtt-slope-area`` binned table used here.
SLOPE_AREA_COLUMNS = (
    "drainage_area", "d05_slope", "d25_slope", "d50_slope", "d75_slope", "d95_slope",
)

#: Columns of ``<prefix>_low_concavity_points.csv``.
LOW_CONCAVITY_POINT_COLUMNS = (
    "flat", "x", "y", "drainage_area_m2", "slope", "theta", "r2", "n_pixels",
    "log10_area_range", "n_junctions_in_window", "low_concavity",
)

#: Attributes of ``<prefix>_low_concavity_segments.fgb``.
LOW_CONCAVITY_SEGMENT_COLUMNS = (
    "segment_id", "link_id", "low_concavity", "n_pixels", "length_m", "median_theta",
    "median_slope", "min_area_m2", "max_area_m2", "stream_order", "upstream_flat",
    "downstream_flat",
)

#: ``slope_algorithm_choice`` names, as the ``method`` of a binned table.
POLYFIT = "polyfit"
D8_FLOW_PATH = "d8_flow_path"

_BINNED_NAME = {
    POLYFIT: re.compile(r"_lambda_([0-9]+(?:_[0-9]+)?)m\.csv$"),
    D8_FLOW_PATH: re.compile(r"_d8window_([0-9]+(?:_[0-9]+)?)m\.csv$"),
}


# -----------------------------------------------------------------------------
# Binned tables
# -----------------------------------------------------------------------------
def load_slope_area(path) -> pd.DataFrame:
    """Load the ``lsdtt-slope-area`` binned table, sorted by drainage area.

    Raises:
        ValueError: on missing or non-finite columns.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Slope-area CSV does not exist: {path}")
    df = pd.read_csv(path, encoding="utf-8-sig")
    missing = [c for c in SLOPE_AREA_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Slope-area CSV {path} is missing column(s): {', '.join(missing)}")
    if df.empty:
        raise ValueError(f"Slope-area CSV contains no rows: {path}")
    for column in SLOPE_AREA_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")
        if not np.all(np.isfinite(df[column].to_numpy(dtype=float))):
            raise ValueError(f"Non-finite value in {column!r} in {path}")
    return df.sort_values("drainage_area", kind="stable").reset_index(drop=True)


def scale_label(method: str, scale_m: float) -> str:
    """Legend label of a binned table: ``polyfit, lambda 90 m`` or
    ``D8 flow path, window 150 m``."""
    if method == POLYFIT:
        return f"polyfit, λ {scale_m:g} m"
    if method == D8_FLOW_PATH:
        return f"D8 flow path, window {scale_m:g} m"
    raise ValueError(f"unknown slope-area method {method!r}")


def find_binned_slope_area(data_dir, prefix: str, method: Optional[str] = None) -> list:
    """The binned slope-area tables of one run, labelled.

    Args:
        data_dir: Directory holding the ``lsdtt-slope-area`` outputs.
        prefix: The run's ``write_fname``.
        method: ``"polyfit"`` or ``"d8_flow_path"`` to keep one kind; None
            for both.

    Returns:
        ``[{"method", "scale_m", "label", "path", "table"}, ...]``, polyfit
        before D8, each sorted by scale. ``table`` is
        :func:`load_slope_area` of ``path``. Empty when the run wrote none.
    """
    data_dir = Path(data_dir)
    methods = (POLYFIT, D8_FLOW_PATH) if method is None else (method,)
    found = []
    for m in methods:
        if m not in _BINNED_NAME:
            raise ValueError(f"unknown slope-area method {m!r}")
        pattern = _BINNED_NAME[m]
        rows = []
        for path in data_dir.iterdir():
            if not path.name.startswith(f"{prefix}_") or not path.is_file():
                continue
            rest = path.name[len(prefix):]
            match = pattern.fullmatch(rest)
            if match is None:
                continue
            scale = float(match.group(1).replace("_", "."))
            rows.append({"method": m, "scale_m": scale, "label": scale_label(m, scale),
                         "path": path})
        for row in sorted(rows, key=lambda r: r["scale_m"]):
            row["table"] = load_slope_area(row["path"])
            found.append(row)
    return found


# -----------------------------------------------------------------------------
# Low concavity
# -----------------------------------------------------------------------------
def low_concavity_files(data_dir, prefix: str) -> dict:
    """Paths of a run's low-concavity outputs (they may not exist):
    ``points``, ``segments``, ``theta``, ``flag``."""
    base = Path(data_dir) / f"{prefix}_low_concavity"
    return {"points": Path(f"{base}_points.csv"), "segments": Path(f"{base}_segments.fgb"),
            "theta": Path(f"{base}_theta.tif"), "flag": Path(f"{base}_flag.tif")}


def load_low_concavity_points(path) -> pd.DataFrame:
    """Load ``<prefix>_low_concavity_points.csv``.

    Numeric columns are parsed as floats (``nan`` where lsdtt3 wrote no fit);
    ``low_concavity`` is an integer (1, 0, -1).

    Raises:
        ValueError: on missing columns.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Low-concavity points CSV does not exist: {path}")
    df = pd.read_csv(path, encoding="utf-8-sig")
    missing = [c for c in LOW_CONCAVITY_POINT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Low-concavity points CSV {path} is missing column(s): "
                         f"{', '.join(missing)}")
    for column in LOW_CONCAVITY_POINT_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df["low_concavity"] = df["low_concavity"].fillna(-1).astype(int)
    return df


def load_low_concavity_segments(path):
    """Load ``<prefix>_low_concavity_segments.fgb`` as a ``GeoDataFrame``.

    Raises:
        ValueError: when the ``low_concavity`` attribute is missing.
    """
    from .vector import load_vector

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Low-concavity segments file does not exist: {path}")
    gdf = load_vector(str(path))
    if "low_concavity" not in gdf.columns:
        raise ValueError(f"{path} has no 'low_concavity' attribute")
    gdf["low_concavity"] = pd.to_numeric(gdf["low_concavity"], errors="coerce").fillna(-1).astype(int)
    return gdf


def load_low_concavity_rasters(theta_path, flag_path) -> pd.DataFrame:
    """The channel pixels of the theta and flag rasters as a table.

    Pixels that are nodata in the flag raster (off the channel network) are
    left out. Theta is NaN where it is nodata.

    Returns:
        ``DataFrame`` with ``theta`` and ``low_concavity`` (1, 0, -1).
    """
    import rasterio

    def read(path):
        with rasterio.open(path) as src:
            band = src.read(1).astype(float)
            nodata = src.nodata
        bad = ~np.isfinite(band)
        if nodata is not None:
            bad |= band == nodata
        band[bad] = np.nan
        return band

    flag = read(flag_path)
    theta = read(theta_path)
    if flag.shape != theta.shape:
        raise ValueError("the theta and flag rasters are not on the same grid")
    on = np.isfinite(flag)
    return pd.DataFrame({"theta": theta[on], "low_concavity": flag[on].astype(int)})


# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
def load_ingested_params(path) -> dict:
    """The ``key: value`` pairs of an ``*_ingested_params.param`` file.

    Comment (``#``) and blank lines are skipped; values are stripped strings
    (empty for ``key:`` with no value).
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Parameter file does not exist: {path}")
    params = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        text = line.strip()
        if not text or text.startswith("#") or ":" not in text:
            continue
        key, value = text.split(":", 1)
        params[key.strip()] = value.strip()
    return params


def param_float(params: dict, key: str, default: Optional[float] = None) -> Optional[float]:
    """``params[key]`` as a float, or ``default`` when absent or not a number."""
    try:
        return float(params[key])
    except (KeyError, ValueError):
        return default
