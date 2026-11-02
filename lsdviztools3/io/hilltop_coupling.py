"""IO for hilltop / channel coupling outputs.

Readers for the hilltop-flow-routing output set: the hilltop points, the ordered
ridgelines, the traces and their summary, and the rasters that go under them.
Nothing here computes a hillslope quantity — that is lsdtt3's job (DESIGN §1);
this module only reads what the tool already wrote.

An output set shares one directory (or one ``<write_fname>`` stem):

* ``hilltop_points``    (.fgb / .gpkg / .csv)  — one row per hilltop pixel,
* ``ridgelines``        (.fgb / .gpkg / .csv)  — the ordered divides,
* ``hillslope_traces``  (.fgb / .gpkg)         — thinned traced flow paths,
* ``trace_summary``     (.csv)                 — dev driver only,
* ``trace_paths``       (.csv)                 — dev driver only,
* ``run_summary``       (.csv)                 — counters, name/value pairs,
* ``*_dem`` / ``*_channels`` / ``*_ridge_side`` (.tif) — rasters.

TWO SHAPES OF THE SAME THING.  ``lsdtt-hillslope-channel-coupling`` writes
FlatGeobuf with real geometry: ridgelines and traces are LineString features and
hilltops are Points carrying their attributes.  The Phase 0/1 development driver
wrote the same information as long-form CSV tables (one row per vertex, with
``x``/``y`` columns).  Both are accepted and normalised here, because a chart
should not care which produced the file.  ``.fgb`` is tried first.

No PyGMT imports here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

#: Vector extensions written by ``write_vector_format``, in preference order.
_VECTOR_SUFFIXES = (".fgb", ".gpkg", ".shp", ".csv")

#: Trace outcome codes, as emitted by ``lsdtt::streamline_outcome``.
OUTCOME_NAME = {
    0: "reached channel",
    1: "hit edge",
    2: "loop detected",
    3: "no aspect",
    4: "not traced (DEM border)",
    5: "stalled on flat (no channel)",
}

#: Stable colours so the same outcome reads the same way in every figure.
OUTCOME_COLOUR = {
    0: "#2a7f3f",
    1: "#d62728",
    2: "#9467bd",
    3: "#ff7f0e",
    4: "#7f7f7f",
    # Distinct from the loop purple: a stall is a statement about the channel
    # network's coverage, not about the routing, and the figures should not
    # invite reading them as the same failure.
    5: "#17becf",
}

#: ``ridge_side`` is left/right of travel along the ordered ridgeline, which
#: partitions a divide's hilltop pixels into the two opposing hillslopes.
SIDE_NAME = {0: "left of travel", 1: "right of travel"}
SIDE_COLOUR = {0: "#1f77b4", 1: "#e6550d"}


def _find(base: Path, name: str) -> Optional[Path]:
    """Locate one member of an output set by stem, trying each vector suffix."""
    if base.is_dir():
        for suffix in _VECTOR_SUFFIXES:
            candidate = base / f"{name}{suffix}"
            if candidate.exists():
                return candidate
        return None

    for suffix in _VECTOR_SUFFIXES:
        candidate = base.with_name(f"{base.name}_{name}{suffix}")
        if candidate.exists():
            return candidate
    return None


def _load_table(path: Path) -> pd.DataFrame:
    """Read one table, whatever container it arrived in.

    Vector files come back as a GeoDataFrame so the geometry survives; CSVs come
    back as a plain DataFrame. Callers that need vertices should go through
    :func:`line_vertices`, which handles both.
    """
    if path.suffix == ".csv":
        return pd.read_csv(path)

    import geopandas as gpd  # imported lazily: CSV-only users need no GDAL

    return gpd.read_file(path)


def has_geometry(frame) -> bool:
    """True if ``frame`` is a GeoDataFrame carrying usable geometry."""
    return frame is not None and hasattr(frame, "geometry") and \
        "geometry" in getattr(frame, "columns", [])


def line_vertices(frame, id_col: str):
    """Yield ``(feature_id, xs, ys, attributes)`` for each line in ``frame``.

    Normalises the two shapes described in the module docstring:

    * a GeoDataFrame of LineStrings — vertices come from the geometry;
    * a long-form table with ``x``/``y`` columns — vertices come from the rows,
      ordered by ``node_index`` or ``vertex`` when present.

    Yields nothing for an empty or missing frame, so callers can loop without
    a guard.
    """
    if frame is None or len(frame) == 0:
        return

    if has_geometry(frame):
        for _, row in frame.iterrows():
            geom = row["geometry"]
            if geom is None or geom.is_empty:
                continue
            # A MultiLineString is emitted as its constituent parts rather than
            # joined, because joining would draw a false segment between them.
            parts = geom.geoms if geom.geom_type.startswith("Multi") else [geom]
            for part in parts:
                xs, ys = part.xy
                yield row.get(id_col, 0), list(xs), list(ys), row
        return

    order_col = next((c for c in ("node_index", "vertex") if c in frame.columns), None)
    for key, group in frame.groupby(id_col):
        ordered = group.sort_values(order_col) if order_col else group
        yield key, ordered["x"].tolist(), ordered["y"].tolist(), ordered.iloc[0]


@dataclass
class HilltopCouplingOutputs:
    """One hilltop-flow-routing output set, loaded.

    Attributes:
        hilltops: One row per hilltop pixel, carrying the ordering attributes
            (``ridgeline_id``, ``ridgeline_node_index``,
            ``distance_along_ridgeline_m``, ``ridge_side``).
        ridgelines: The ordered divide nodes, long-form, one row per node.
        traces: One row per trace (``length_m``, ``relief_m``, ``outcome``).
        paths: Thinned trace geometry as a long-form CSV (dev driver).
        trace_lines: Thinned trace geometry as LineStrings (product CLI).
        summary: The run counters as a name -> value mapping.
        dem_path, channels_path, ridge_side_path: Raster paths, when present.
    """

    hilltops: Optional[pd.DataFrame] = None
    ridgelines: Optional[pd.DataFrame] = None
    traces: Optional[pd.DataFrame] = None
    paths: Optional[pd.DataFrame] = None
    trace_lines: Optional[pd.DataFrame] = None
    summary: dict = field(default_factory=dict)

    dem_path: Optional[Path] = None
    channels_path: Optional[Path] = None
    ridge_side_path: Optional[Path] = None

    def any_traces(self):
        """Whichever trace geometry is present, product or dev-driver."""
        if self.trace_lines is not None and len(self.trace_lines):
            return self.trace_lines
        return self.paths

    def trace_table(self) -> Optional[pd.DataFrame]:
        """The per-trace attribute table.

        The product CLI has no separate trace summary -- the outcome and the
        metrics live on the hilltop points -- so fall back to those.
        """
        if self.traces is not None and len(self.traces):
            return self.traces
        return self.hilltops

    def outcome_counts(self) -> pd.Series:
        """Traces per outcome code, ascending by code.

        Returns:
            A Series indexed by outcome code; empty if no traces were loaded.
        """
        table = self.trace_table()
        if table is None:
            return pd.Series(dtype=int)
        col = "outcome" if "outcome" in table.columns else "trace_status"
        if col not in table.columns:
            return pd.Series(dtype=int)
        return table[col].value_counts().sort_index()

    def successful_traces(self) -> pd.DataFrame:
        """The traces that reached a channel.

        Every other outcome leaves ``length_m`` meaningless, so any distribution
        of hillslope length must be taken over this subset and the discarded
        count reported alongside it.
        """
        table = self.trace_table()
        if table is None:
            return pd.DataFrame()
        col = "outcome" if "outcome" in table.columns else "trace_status"
        if col not in table.columns:
            return pd.DataFrame()
        return table[table[col] == 0]

    def border_dropped(self) -> int:
        """Hilltop pixels dropped before tracing because they sit on the border.

        The tracer enumerates sources strictly inside the outer ring, so a
        hilltop pixel on the DEM edge produces no trace and no record at all.
        These pixels appear in ``hilltops`` and in no trace figure, so the count
        belongs on any figure that shows trace outcomes.
        """
        return int(self.summary.get("n_dropped_at_border", 0))


def load_hilltop_outputs(prefix: str | Path) -> HilltopCouplingOutputs:
    """Load one hilltop-flow-routing output set.

    Args:
        prefix: A directory holding one output set, or a ``<write_path>/
            <write_fname>`` stem.

    Returns:
        The loaded set. Members that are absent are left as ``None`` rather than
        raising, so a partial run can still be inspected.
    """
    base = Path(prefix)
    out = HilltopCouplingOutputs()

    for attr, stem in (
        ("hilltops", "hilltop_points"),
        ("ridgelines", "ridgelines"),
        ("traces", "trace_summary"),
        ("paths", "trace_paths"),
        ("trace_lines", "hillslope_traces"),
    ):
        found = _find(base, stem)
        if found is not None:
            setattr(out, attr, _load_table(found))

    summary_path = _find(base, "run_summary")
    if summary_path is not None:
        table = pd.read_csv(summary_path)
        if {"quantity", "value"}.issubset(table.columns):
            out.summary = dict(zip(table["quantity"], table["value"]))

    search_dir = base if base.is_dir() else base.parent

    # Prefer a VECTOR channel network: it is what the extraction step already
    # produced, and drawing it needs no rasterisation onto the analysis grid.
    for pattern in ("*channel_network_lines.fgb", "*channel_network_lines.gpkg",
                    "*channel_network_points.fgb", "*channel_mask.tif"):
        found = sorted(search_dir.glob(pattern))
        if found:
            out.channels_path = found[0]
            break

    for attr, needle in (
        ("dem_path", "dem"),
        ("ridge_side_path", "ridge_side"),
    ):
        for candidate in sorted(search_dir.glob(f"*{needle}.tif")):
            setattr(out, attr, candidate)
            break

    return out
