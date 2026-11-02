"""IO for ``lsdtt-divide-migration`` outputs.

Readers plus the small amount of *metadata* needed to label a figure honestly:
which summary columns belong to which metric, which sign of ``delta`` points at
which side, and what each ``uncertainty_method`` actually reports. Nothing here
computes a divide-migration statistic — that is lsdtt3's job (DESIGN §1); this
module only reads what the CLI already wrote.

An output set shares one ``<write_fname>`` stem:

* ``<stem>_divides.fgb`` / ``.gpkg``          — divide polylines,
* ``<stem>_channel_heads.fgb`` / ``.gpkg``    — channel-head points,
* ``<stem>_divide_migration_metrics.csv``     — the same head metrics as a table,
* ``<stem>_divide_pair_summary.csv``          — one side-vs-side row per pair.

No PyGMT imports here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

from .table import load_table
from .vector import load_vector

#: Vector extensions written by ``write_vector_format`` (in preference order).
_VECTOR_SUFFIXES = (".fgb", ".gpkg", ".shp")


# --------------------------------------------------------------------------
# Metric metadata
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class DivideMetric:
    """A single divide-migration metric and how it reads across a divide.

    Attributes:
        key: Canonical short name (``chi``, ``relief``, ``gradient``,
            ``head_elevation``) — the stem used by the pair-summary columns.
        head_column: Column name in the channel-head table/vector.
        label: Human-readable axis label (no units).
        units: Unit string for axis labels, or ``""``.
        aggressor_delta_sign: Sign of ``delta_<key>`` (= A − B) that makes
            **side A** the aggressor: ``-1`` for head elevation and chi (A's
            heads sit lower / at lower chi), ``+1`` for gradient and relief.
    """

    key: str
    head_column: str
    label: str
    units: str
    aggressor_delta_sign: int

    @property
    def mean_a_column(self) -> str:
        return f"mean_{self.key}_a"

    @property
    def mean_b_column(self) -> str:
        return f"mean_{self.key}_b"

    @property
    def delta_column(self) -> str:
        return f"delta_{self.key}"

    @property
    def uncertainty_column(self) -> str:
        return f"uncertainty_{self.key}"

    @property
    def n_a_column(self) -> str:
        return f"n_{self.key}_a"

    @property
    def n_b_column(self) -> str:
        return f"n_{self.key}_b"

    @property
    def axis_label(self) -> str:
        return f"{self.label} ({self.units})" if self.units else self.label

    @property
    def aggressor_rule(self) -> str:
        """One-line statement of which delta sign implicates side A."""
        sign = "< 0" if self.aggressor_delta_sign < 0 else "> 0"
        return f"A is the aggressor when Δ {sign}"


#: The four metrics the tool reports, keyed by their pair-summary stem.
DIVIDE_METRICS = {
    "chi": DivideMetric("chi", "chi", "Mean channel-head χ", "", -1),
    "relief": DivideMetric(
        "relief", "relief_m", "Mean relief to divide", "m", +1
    ),
    "gradient": DivideMetric(
        "gradient", "gradient_to_divide", "Mean gradient to divide", "m/m", +1
    ),
    "head_elevation": DivideMetric(
        "head_elevation", "head_elevation_m", "Mean channel-head elevation",
        "m", -1
    ),
}

#: Accept the channel-head column names as metric aliases too, so a user can
#: say ``--metric relief_m`` (what they see in the CSV) or ``--metric relief``
#: (what they see in the pair summary).
_METRIC_ALIASES = {
    "chi": "chi",
    "relief": "relief",
    "relief_m": "relief",
    "gradient": "gradient",
    "gradient_to_divide": "gradient",
    "head_elevation": "head_elevation",
    "head_elevation_m": "head_elevation",
    "elevation": "head_elevation",
}


def resolve_metric(metric) -> DivideMetric:
    """Resolve a metric name (canonical key or head-column alias).

    Args:
        metric: ``"chi"``, ``"relief"`` / ``"relief_m"``, ``"gradient"`` /
            ``"gradient_to_divide"``, ``"head_elevation"`` /
            ``"head_elevation_m"``, or a :class:`DivideMetric`.

    Returns:
        The :class:`DivideMetric`.

    Raises:
        KeyError: if the name is not a known metric.
    """
    if isinstance(metric, DivideMetric):
        return metric
    key = _METRIC_ALIASES.get(str(metric).strip().lower())
    if key is None:
        raise KeyError(
            f"Unknown divide-migration metric {metric!r}; choose from "
            f"{sorted(set(_METRIC_ALIASES))}."
        )
    return DIVIDE_METRICS[key]


#: Every metric name a CLI may accept, for argparse ``choices``.
METRIC_CHOICES = tuple(sorted(_METRIC_ALIASES))


# --------------------------------------------------------------------------
# Uncertainty-method semantics
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class UncertaintyMethod:
    """What an ``uncertainty_*`` column means for a given ``uncertainty_method``.

    Attributes:
        key: The value written in the ``uncertainty_method`` column.
        label: Axis/legend label describing the number.
        is_dispersion: ``True`` when the column is a spread that can be drawn as
            an error bar; ``False`` for ``ttest``, where it is a **p-value** and
            drawing it as a bar would be plainly wrong.
        significance_multiplier: ``|delta| > multiplier x uncertainty`` marks a
            separation of the two side means. ``None`` when the column supports
            no such test (``std_dev`` is a spread of the head populations, not
            of the delta, and does not shrink with n).
    """

    key: str
    label: str
    is_dispersion: bool
    significance_multiplier: Optional[float]


#: The four methods ``pair_uncertainty_method`` can produce.
UNCERTAINTY_METHODS = {
    "std_dev": UncertaintyMethod(
        "std_dev",
        "pooled SD of the head populations (not an interval on Δ)",
        True,
        None,
    ),
    "std_err": UncertaintyMethod(
        "std_err", "standard error of Δ (±1 s.e.)", True, 2.0
    ),
    "bootstrap": UncertaintyMethod(
        "bootstrap", "bootstrap 95% CI half-width of Δ", True, 1.0
    ),
    "ttest": UncertaintyMethod(
        "ttest", "Welch t-test p-value (not a dispersion)", False, None
    ),
}

#: Used when a file carries an unrecognised or missing ``uncertainty_method``.
UNKNOWN_UNCERTAINTY = UncertaintyMethod(
    "unknown", "uncertainty (method not recorded)", True, None
)


def uncertainty_method(summary: pd.DataFrame) -> UncertaintyMethod:
    """Read the ``uncertainty_method`` recorded in a pair-summary table.

    The method is a property of the run, so the first non-null value wins. An
    unrecognised or absent method falls back to :data:`UNKNOWN_UNCERTAINTY`,
    which claims nothing about significance.

    Args:
        summary: Pair-summary table (see :func:`load_pair_summary`).

    Returns:
        The :class:`UncertaintyMethod`.
    """
    if "uncertainty_method" not in summary.columns:
        return UNKNOWN_UNCERTAINTY
    values = summary["uncertainty_method"].dropna()
    if values.empty:
        return UNKNOWN_UNCERTAINTY
    return UNCERTAINTY_METHODS.get(
        str(values.iloc[0]).strip().lower(), UNKNOWN_UNCERTAINTY
    )


# --------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------

def load_divides(path, layer=None):
    """Load the divide polylines (``<stem>_divides.fgb`` / ``.gpkg``).

    Args:
        path: Path to the divide vector file.
        layer: Optional layer name (GeoPackage).

    Returns:
        ``GeoDataFrame`` with ``pair_id``, ``unit_a``, ``unit_b``, ``junc_a``,
        ``junc_b`` (``-1`` on merged features), ``segment_id``, ``n_fragments``,
        ``n_vertices`` and ``length_m``.
    """
    return load_vector(path, layer=layer)


def load_channel_heads(path, *, layer=None, crs=None):
    """Load channel heads from the point vector or the metrics CSV.

    A ``.csv`` carries ``head_x`` / ``head_y`` in the DEM's CRS but records no
    CRS of its own, so it is returned as a plain ``DataFrame`` unless ``crs`` is
    supplied, in which case the points are rebuilt as a ``GeoDataFrame``.

    Args:
        path: ``<stem>_channel_heads.fgb`` / ``.gpkg``, or
            ``<stem>_divide_migration_metrics.csv``.
        layer: Optional layer name (GeoPackage).
        crs: CRS to attach when reading a CSV (e.g. the DEM's CRS).

    Returns:
        ``GeoDataFrame`` (vector input, or CSV + ``crs``) else ``DataFrame``.
    """
    path = Path(path)
    if path.suffix.lower() != ".csv":
        return load_vector(path, layer=layer)

    df = load_table(path)
    if crs is None:
        return df
    import geopandas as gpd

    return gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["head_x"], df["head_y"]),
        crs=crs,
    )


def load_pair_summary(path) -> pd.DataFrame:
    """Load ``<stem>_divide_pair_summary.csv`` (one row per divide pair).

    Numeric columns are coerced with ``errors="coerce"`` so an empty or ``nan``
    field becomes ``NaN`` rather than an object-dtype string — a whole row is
    legitimately NaN when a side has fewer than two valid values.

    Args:
        path: Path to the pair-summary CSV.

    Returns:
        The table as a ``DataFrame`` with ``pair_id`` as ``str``.
    """
    df = pd.read_csv(path)
    for column in df.columns:
        if column in ("pair_id", "uncertainty_method"):
            continue
        df[column] = pd.to_numeric(df[column], errors="coerce")
    if "pair_id" in df.columns:
        df["pair_id"] = df["pair_id"].astype(str)
    return df


def find_outputs(prefix) -> dict:
    """Locate a divide-migration output set from a stem or a directory.

    Args:
        prefix: Either the ``<write_path>/<write_fname>`` stem (e.g.
            ``/data/divides/sc``) or a directory holding exactly one output set,
            whose stem is inferred from ``*_divide_pair_summary.csv`` /
            ``*_divide_migration_metrics.csv`` / ``*_divides.*``.

    Returns:
        Dict with the keys ``stem``, ``divides``, ``channel_heads``,
        ``metrics_csv`` and ``pair_summary``. A value is ``None`` when that
        output was not written (the switches are independent).

    Raises:
        FileNotFoundError: if no divide-migration output can be found.
    """
    prefix = Path(prefix)
    if prefix.is_dir():
        stem = _infer_stem(prefix)
    else:
        stem = prefix

    out = {
        "stem": stem,
        "divides": _first_existing_vector(stem, "_divides"),
        "channel_heads": _first_existing_vector(stem, "_channel_heads"),
        "metrics_csv": _existing(stem.with_name(stem.name + "_divide_migration_metrics.csv")),
        "pair_summary": _existing(stem.with_name(stem.name + "_divide_pair_summary.csv")),
    }
    if not any(out[k] for k in ("divides", "channel_heads", "metrics_csv", "pair_summary")):
        raise FileNotFoundError(
            f"No lsdtt-divide-migration outputs found for prefix {prefix}. "
            "Expected files named <stem>_divides.fgb, <stem>_channel_heads.fgb, "
            "<stem>_divide_migration_metrics.csv or <stem>_divide_pair_summary.csv."
        )
    return out


def _existing(path: Path):
    return path if path.exists() else None


def _first_existing_vector(stem: Path, tag: str):
    for suffix in _VECTOR_SUFFIXES:
        candidate = stem.with_name(stem.name + tag + suffix)
        if candidate.exists():
            return candidate
    return None


def _infer_stem(directory: Path) -> Path:
    """Infer the output stem from a directory holding one output set."""
    markers = (
        ("_divide_pair_summary.csv", "*_divide_pair_summary.csv"),
        ("_divide_migration_metrics.csv", "*_divide_migration_metrics.csv"),
    )
    for tag, pattern in markers:
        matches = sorted(directory.glob(pattern))
        if matches:
            return directory / matches[0].name[: -len(tag)]
    for suffix in _VECTOR_SUFFIXES:
        matches = sorted(directory.glob(f"*_divides{suffix}"))
        if matches:
            return directory / matches[0].name[: -len(f"_divides{suffix}")]
    raise FileNotFoundError(
        f"Could not infer a divide-migration output stem in {directory}."
    )
