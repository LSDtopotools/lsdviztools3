"""Ridgeline IO for ``lsdtt-ridgeline-analysis`` outputs.

Readers for the ridgeline profile CSVs (``*_ridgelines.csv`` and
``*_with_polyfit.csv``), the triple-divide junction (TDJ) CSVs and the
``lsdtt-slope-area`` table, plus the profile selection rules the ridgeline
charts share (one ridgeline, optionally reversed, optionally clipped to a
distance range). No PyGMT imports here.

Profile CSV columns (one row per ridgeline node)::

    ridgeline_number (basin_junction_number for single-basin runs), row, col,
    latitude, longitude, ridgeline_type, distance, elevation, fitted_elevation,
    detrended_elevation, smoothed_detrended_elevation, smoothed_detrended_slope,
    [moving_window_fitted_elevation_<N>px, moving_window_detrended_elevation_<N>px]
    [<metric>_lambda_<L>m ... in *_with_polyfit.csv; a bare <metric> column when
     one polyfit scale was run]

Distances are metres along the ridge crest. ``row``/``col`` index the divide
grid, whose nodes are at DEM pixel corners.

Derived statistics
------------------
The zero-crossing spacing of maximal curvature and the slope-area peak
(:func:`curvature_spacing_quantiles`, :func:`periodic_median_crossing_spacing`,
:func:`slope_area_peak`, :func:`spacing_summary`) are ported unchanged from
R. Vervoordeldonk's ``ridgeline_plotting`` package. They are analysis, not
display, and are due to move into ``lsdtt-ridgeline-analysis`` /
``lsdtt-slope-area`` as output columns; the charts will then only read them.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

#: Columns that can identify the ridgeline of a profile row, in order of
#: preference (interfluvial / range runs write the first, single-basin runs the
#: second).
RIDGELINE_ID_COLUMNS = ("ridgeline_number", "basin_junction_number")

#: Quantiles of the zero-crossing spacing (name, percentile).
SPACING_QUANTILES = (("d05", 5), ("d25", 25), ("d50", 50), ("d75", 75), ("d95", 95))

#: Columns of the ``lsdtt-slope-area`` binned table used here.
SLOPE_AREA_COLUMNS = (
    "drainage_area", "d05_slope", "d25_slope", "d50_slope", "d75_slope", "d95_slope",
)

_LAMBDA_COLUMN = r"{metric}(?:_lambda_([0-9]+(?:\.[0-9]+)?)m)?"
_MOVING_WINDOW = re.compile(r"moving_window_detrended_elevation_([0-9]+)px$")


# -----------------------------------------------------------------------------
# Profile CSVs
# -----------------------------------------------------------------------------
def ridgeline_id_column(columns: Iterable[str]) -> str:
    """The column naming the ridgeline of each row.

    Args:
        columns: Column names of a ridgeline profile CSV.

    Returns:
        ``"ridgeline_number"`` or ``"basin_junction_number"``.

    Raises:
        ValueError: if neither is present.
    """
    columns = set(columns)
    for name in RIDGELINE_ID_COLUMNS:
        if name in columns:
            return name
    raise ValueError(
        "Profile CSV has no ridgeline identifier column "
        f"({' or '.join(RIDGELINE_ID_COLUMNS)})")


def load_ridgeline_profile(path, ridgeline_number: int = 0,
                           required: Sequence[str] = ("distance",)) -> pd.DataFrame:
    """Load the rows of one ridgeline from a ridgeline profile CSV.

    Rows are sorted by ``distance``.

    Args:
        path: ``*_ridgelines.csv`` or ``*_with_polyfit.csv`` from
            ``lsdtt-ridgeline-analysis``.
        ridgeline_number: The ridgeline to keep (``ridgeline_number``, or
            ``basin_junction_number`` for single-basin runs).
        required: Columns that must be present and numeric (finite) for the
            selected ridgeline.

    Returns:
        The selected rows as a ``pandas.DataFrame``.

    Raises:
        ValueError: on missing columns, no rows for the ridgeline, or
            non-numeric / non-finite values in ``required`` columns.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Profile CSV does not exist: {path}")
    df = pd.read_csv(path, encoding="utf-8-sig")
    id_column = ridgeline_id_column(df.columns)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Profile CSV {path} is missing required column(s): {', '.join(missing)}")
    selected = df[df[id_column].astype(str) == str(int(ridgeline_number))].copy()
    if selected.empty:
        raise ValueError(f"Profile CSV {path} contains no rows for ridgeline {ridgeline_number}")
    for column in required:
        values = pd.to_numeric(selected[column], errors="coerce")
        if not np.all(np.isfinite(values.to_numpy(dtype=float))):
            raise ValueError(f"Profile CSV {path} column {column!r} has non-numeric or "
                             "non-finite values for this ridgeline")
        selected[column] = values
    return selected.sort_values("distance", kind="stable").reset_index(drop=True)


def lambda_columns(columns: Iterable[str], metric: str = "maximum_curvature") -> list:
    """The polyfit scales present for one metric.

    Args:
        columns: Column names of a ``*_with_polyfit.csv``.
        metric: Polyfit metric name, e.g. ``"maximum_curvature"``,
            ``"laplacian"``.

    Returns:
        ``[(lambda_m, column), ...]`` sorted by lambda. A bare ``<metric>``
        column (one polyfit scale) is reported with lambda 0.
    """
    pattern = re.compile(_LAMBDA_COLUMN.format(metric=re.escape(metric)) + "$")
    found = []
    for column in columns:
        match = pattern.fullmatch(column)
        if match:
            found.append((float(match.group(1)) if match.group(1) else 0.0, column))
    return sorted(found)


def lambda_column(metric: str, lambda_m: float) -> str:
    """Column name of a polyfit metric at one scale, as lsdtt3 writes it."""
    return f"{metric}_lambda_{lambda_m:g}m"


def moving_window_sizes(columns: Iterable[str]) -> list:
    """Pixel window sizes of the ``moving_window_detrended_elevation_<N>px`` columns."""
    return sorted(int(m.group(1)) for c in columns if (m := _MOVING_WINDOW.fullmatch(c)))


def moving_window_column(window_pixels: int) -> str:
    """Column name of the moving-window detrended elevation for one window."""
    return f"moving_window_detrended_elevation_{int(window_pixels)}px"


def display_distance(distance, reverse: bool = False) -> np.ndarray:
    """Distances along the displayed profile.

    Unreversed, distances are kept as written. Reversed, they run from 0 at the
    former far end: ``max - (d - min) - min``.

    Args:
        distance: Profile distances (m).
        reverse: Reverse the profile direction.

    Returns:
        The display distances (not re-sorted).
    """
    d = np.asarray(distance, dtype=float)
    if not reverse:
        return d.copy()
    lo, hi = float(np.min(d)), float(np.max(d))
    return hi - (d - lo) - lo


def select_profile(profile: pd.DataFrame, reverse: bool = False,
                   distance_range: Optional[tuple] = None) -> pd.DataFrame:
    """Apply the shared reverse / distance-range rules to one ridgeline.

    Adds a ``display_distance`` column, sorts by it, and keeps rows within the
    inclusive ``distance_range`` (in display distance).

    Args:
        profile: Rows of one ridgeline (:func:`load_ridgeline_profile`).
        reverse: Reverse the profile direction.
        distance_range: ``(minimum, maximum)`` display distance, or None.

    Returns:
        The selected rows, sorted by ``display_distance``.

    Raises:
        ValueError: if the range is inverted or leaves fewer than two rows.
    """
    out = profile.copy()
    out["display_distance"] = display_distance(out["distance"].to_numpy(), reverse)
    out = out.sort_values("display_distance", kind="stable").reset_index(drop=True)
    if distance_range is not None:
        lo, hi = (float(v) for v in distance_range)
        if not (math.isfinite(lo) and math.isfinite(hi)) or lo > hi:
            raise ValueError("The distance range must be finite with minimum <= maximum")
        out = out[(out["display_distance"] >= lo) & (out["display_distance"] <= hi)]
        out = out.reset_index(drop=True)
    if len(out) < 2:
        raise ValueError("The selected profile must contain at least two samples")
    return out


def node_spacing(distance, tolerance: float = 0.05) -> float:
    """Median distance between consecutive ridgeline nodes (the grid spacing).

    Args:
        distance: Increasing profile distances (m).
        tolerance: Relative difference allowed between the median and mean
            node spacing.

    Returns:
        The median spacing in metres.

    Raises:
        ValueError: if distances do not increase strictly, or the spacing is
            not uniform within ``tolerance``.
    """
    steps = np.diff(np.asarray(distance, dtype=float))
    if steps.size == 0 or np.any(steps <= 0):
        raise ValueError("Profile distances must increase strictly")
    median = float(np.median(steps))
    if not math.isclose(median, float(np.mean(steps)), rel_tol=tolerance):
        raise ValueError(f"Profile node spacing varies by more than {tolerance:.0%}; "
                         "cannot derive one grid spacing")
    return median


def divide_grid_spacing(profile: pd.DataFrame) -> float:
    """Grid spacing from profile distance per divide-grid node step (row/col)."""
    d = profile["distance"].to_numpy(dtype=float)
    rows = profile["row"].to_numpy(dtype=float)
    cols = profile["col"].to_numpy(dtype=float)
    pixel = np.hypot(np.diff(rows), np.diff(cols))
    dist = np.diff(d)
    ok = np.isfinite(pixel) & (pixel > 0) & np.isfinite(dist) & (dist > 0)
    if not np.any(ok):
        raise ValueError("Cannot infer the grid spacing from consecutive ridgeline nodes")
    dx = float(np.median(dist[ok] / pixel[ok]))
    if not (math.isfinite(dx) and dx > 0):
        raise ValueError("The inferred grid spacing must be positive and finite")
    return dx


# -----------------------------------------------------------------------------
# Triple divide junctions
# -----------------------------------------------------------------------------
_TDJ_THRESHOLD = {
    "native": re.compile(r"area_([0-9]+(?:\.[0-9]+)?)m2", re.IGNORECASE),
    "topological": re.compile(r"triple_divides_([0-9]+(?:\.[0-9]+)?)m", re.IGNORECASE),
}


def tdj_threshold_from_filename(path, kind: str) -> float:
    """The threshold in a TDJ file name.

    Native files are ``*_native_triple_divides_area_<A>m2.csv`` (A in m²);
    topological files ``*_topological_distance_triple_divides_<D>m.csv``.

    Args:
        path: TDJ CSV path.
        kind: ``"native"`` or ``"topological"``.

    Returns:
        The threshold (m² or m).

    Raises:
        ValueError: if the name holds no single threshold.
    """
    if kind not in _TDJ_THRESHOLD:
        raise ValueError("kind must be 'native' or 'topological'")
    matches = _TDJ_THRESHOLD[kind].findall(Path(path).name)
    if len(matches) != 1:
        expected = "area_<threshold>m2" if kind == "native" else "triple_divides_<threshold>m"
        raise ValueError(f"Cannot parse a {kind} TDJ threshold from {Path(path).name!r}; "
                         f"expected a single {expected}")
    return float(matches[0])


def load_tdj(path, kind: str, ridgeline_number: int) -> dict:
    """Load one TDJ CSV, keeping the junctions on one ridgeline.

    Args:
        path: TDJ CSV from ``lsdtt-ridgeline-analysis``.
        kind: ``"native"`` or ``"topological"``.
        ridgeline_number: The ridgeline to keep.

    Returns:
        ``{"kind", "threshold", "path", "rows"}`` with ``rows`` a DataFrame
        holding a numeric ``distance`` column (``distance_along_ridgeline_m``,
        or ``distance``) plus ``latitude`` / ``longitude``.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"TDJ CSV does not exist: {path}")
    df = pd.read_csv(path, encoding="utf-8-sig")
    needed = {"ridgeline_number", "latitude", "longitude"}
    if not needed.issubset(df.columns):
        raise ValueError(f"TDJ CSV {path} is missing column(s): "
                         f"{', '.join(sorted(needed - set(df.columns)))}")
    dist_col = ("distance_along_ridgeline_m" if "distance_along_ridgeline_m" in df.columns
                else "distance" if "distance" in df.columns else None)
    if dist_col is None:
        raise ValueError(f"TDJ CSV {path} must contain distance_along_ridgeline_m or distance")
    rows = df[df["ridgeline_number"].astype(str) == str(int(ridgeline_number))].copy()
    rows["distance"] = pd.to_numeric(rows[dist_col], errors="coerce")
    for column in ("distance", "latitude", "longitude"):
        rows[column] = pd.to_numeric(rows[column], errors="coerce")
        if rows[column].isna().any():
            raise ValueError(f"TDJ CSV {path} has non-numeric {column} values")
    return {"kind": kind, "threshold": tdj_threshold_from_filename(path, kind),
            "path": path, "rows": rows.reset_index(drop=True)}


def tdj_rows_in_range(tdj: dict, profile_distance, reverse: bool = False,
                      distance_range: Optional[tuple] = None) -> pd.DataFrame:
    """The junctions of a TDJ set that fall on the displayed profile.

    Uses the same reverse rule as :func:`display_distance`, anchored on the
    full profile's distance range (``profile_distance`` is the unclipped
    ridgeline), and keeps junctions within ``distance_range`` (default: the
    whole profile).

    Returns:
        The kept rows with a ``display_distance`` column.
    """
    d = np.asarray(profile_distance, dtype=float)
    lo, hi = float(np.min(d)), float(np.max(d))
    rows = tdj["rows"].copy()
    pos = rows["distance"].to_numpy(dtype=float)
    if reverse:
        pos = hi - (pos - lo) - lo
    rows["display_distance"] = pos
    shown = display_distance(d, reverse)
    rng = distance_range if distance_range is not None else (float(shown.min()),
                                                             float(shown.max()))
    return rows[(pos >= rng[0]) & (pos <= rng[1])].reset_index(drop=True)


def tdj_display_positions(tdj: dict, profile_distance, reverse: bool = False,
                          distance_range: Optional[tuple] = None) -> np.ndarray:
    """Display distances of a TDJ set along a ridgeline profile
    (:func:`tdj_rows_in_range`)."""
    return tdj_rows_in_range(tdj, profile_distance, reverse,
                             distance_range)["display_distance"].to_numpy(dtype=float)


# -----------------------------------------------------------------------------
# Slope-area
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


# -----------------------------------------------------------------------------
# Derived statistics (ported unchanged; due to move into lsdtt3)
# -----------------------------------------------------------------------------
def curvature_spacing_quantiles(profile: pd.DataFrame, reverse: bool = False,
                                distance_range: Optional[tuple] = None,
                                absolute_curvature: bool = False) -> dict:
    """Spacing of maximal-curvature sign changes at every polyfit scale.

    For each ``maximum_curvature_lambda_<L>m`` column: the distances between
    successive sign changes (a run of zeros counts as one), or, with
    ``absolute_curvature``, the absolute curvature values; summarised by the
    5/25/50/75/95th percentiles. Scales with fewer than two crossings are
    skipped. The grid spacing is :func:`node_spacing` of the whole (reversed)
    profile, taken before the distance range is applied.

    Args:
        profile: Rows of one ridgeline (:func:`load_ridgeline_profile`).
        reverse: Reverse the profile direction.
        distance_range: ``(minimum, maximum)`` display distance, or None.
        absolute_curvature: Summarise |curvature| instead of spacing.

    Returns:
        ``{"polyfit_lambda", "d05", "d25", "d50", "d75", "d95"}`` arrays, and
        ``"grid_resolution"`` (m).
    """
    scales = [(lam, col) for lam, col in lambda_columns(profile.columns) if lam > 0]
    if not scales:
        raise ValueError("No maximum_curvature_lambda_<N>m columns found")
    d = profile["distance"].to_numpy(dtype=float)
    finite = np.isfinite(d)
    rows = profile[finite].copy()
    d = d[finite]
    order = np.argsort(d, kind="stable")
    rows, d = rows.iloc[order], d[order]
    if d.size < 2:
        raise ValueError("At least two valid distance samples are needed")
    lo, hi = float(d.min()), float(d.max())
    disp = lo + hi - d if reverse else d.copy()
    order = np.argsort(disp, kind="stable")
    rows, disp = rows.iloc[order], disp[order]
    grid_resolution = node_spacing(disp)
    if distance_range is not None:
        keep = (disp >= distance_range[0]) & (disp <= distance_range[1])
        rows, disp = rows[keep], disp[keep]
    if disp.size < 2:
        raise ValueError("The selected distance range must contain at least two samples")

    out = {"polyfit_lambda": [], **{name: [] for name, _ in SPACING_QUANTILES}}
    for lam, column in scales:
        values = pd.to_numeric(rows[column], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(values)
        values = values[valid]
        if absolute_curvature:
            values = np.abs(values)
        else:
            signs = np.sign(values)
            nonzero = np.flatnonzero(signs != 0)
            if nonzero.size < 2:
                continue
            signs = signs[nonzero]
            crossing = nonzero[1:][signs[1:] != signs[:-1]]
            if crossing.size < 2:
                continue
            values = np.diff(disp[valid][crossing])
        if values.size == 0:
            continue
        out["polyfit_lambda"].append(lam)
        for name, pct in SPACING_QUANTILES:
            out[name].append(float(np.percentile(values, pct)))
    result = {k: np.asarray(v, dtype=float) for k, v in out.items()}
    result["grid_resolution"] = grid_resolution
    return result


def periodic_median_crossing_spacing(distance, curvature) -> float:
    """Median spacing of curvature zero crossings on a periodic profile.

    Crossings are interpolated linearly between nodes, and the profile wraps
    (its last node joins its first), as for a main drainage divide across a
    periodic model domain. Nodata (-9999) and non-finite values are dropped.

    Returns:
        The median gap in metres, or NaN with fewer than two crossings.
    """
    d = np.asarray(distance, dtype=float)
    c = np.asarray(curvature, dtype=float)
    ok = np.isfinite(d) & np.isfinite(c) & (c != -9999)
    d, c = d[ok], c[ok]
    if d.size < 2:
        return float("nan")
    order = np.argsort(d)
    d, c = d[order], c[order]
    period = d[-1] - d[0]
    if period <= 0:
        return float("nan")
    crossings = []
    for i in range(d.size):
        j = (i + 1) % d.size
        x0, x1 = d[i], d[j] + (period if j == 0 else 0.0)
        y0, y1 = c[i], c[j]
        if y0 == 0:
            crossings.append(x0)
        elif y0 * y1 < 0 and y1 != y0:
            crossings.append(x0 - y0 * (x1 - x0) / (y1 - y0))
    if len(crossings) < 2:
        return float("nan")
    crossings = np.sort(np.asarray(crossings) % period)
    gaps = np.diff(np.r_[crossings, crossings[0] + period])
    return float(np.median(gaps))


def slope_area_peak(slope_area: pd.DataFrame) -> float:
    """Drainage area of the bin with the largest median slope (positive areas only)."""
    area = slope_area["drainage_area"].to_numpy(dtype=float)
    positive = area > 0
    if not np.any(positive):
        raise ValueError("The slope-area table has no positive drainage areas")
    median = slope_area["d50_slope"].to_numpy(dtype=float)[positive]
    return float(area[positive][int(np.nanargmax(median))])


def characteristic_length(area_peak: float, grid_resolution: float) -> float:
    """``L = sqrt((A_peak / dx)^2)``, as defined in the original plotting code."""
    return math.sqrt((area_peak / grid_resolution) ** 2)


def spacing_summary(name: str, spacing: dict, slope_area: pd.DataFrame) -> dict:
    """One site's characteristic-length and spacing summary row.

    Returns:
        ``site``, ``A_max_median_slope``, ``dx``, ``L``, ``2xL``,
        ``closest_lambda_to_2xL``, ``median_zero_crossing_spacing_at_closest_lambda``
        and ``median_zero_spacing_divided_by_L`` (the last three empty when no
        scale had a spacing).
    """
    area_peak = slope_area_peak(slope_area)
    dx = float(spacing["grid_resolution"])
    length = characteristic_length(area_peak, dx)
    lambdas = spacing["polyfit_lambda"]
    row = {"site": name, "A_max_median_slope": area_peak, "dx": dx, "L": length,
           "2xL": 2 * length, "closest_lambda_to_2xL": "",
           "median_zero_crossing_spacing_at_closest_lambda": "",
           "median_zero_spacing_divided_by_L": ""}
    if lambdas.size:
        i = int(np.argmin(np.abs(lambdas - 2 * length)))
        row["closest_lambda_to_2xL"] = float(lambdas[i])
        row["median_zero_crossing_spacing_at_closest_lambda"] = float(spacing["d50"][i])
        row["median_zero_spacing_divided_by_L"] = float(spacing["d50"][i]) / length
    return row
