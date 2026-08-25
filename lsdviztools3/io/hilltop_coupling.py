"""IO for hilltop / channel coupling outputs.

Readers for the hilltop-flow-routing output set: the hilltop points, the ordered
ridgelines, the traces and their summary, and the rasters that go under them.
Nothing here computes a hillslope quantity — that is lsdtt3's job (DESIGN §1);
this module only reads what the tool already wrote.

An output set shares one directory (or one ``<write_fname>`` stem):

* ``hilltop_points`` (.fgb / .gpkg / .csv)  — one row per hilltop pixel,
* ``ridgelines``     (.fgb / .gpkg / .csv)  — the ordered divides themselves,
* ``trace_summary``  (.csv)                 — one row per trace,
* ``trace_paths``    (.fgb / .gpkg / .csv)  — thinned trace geometry,
* ``run_summary``    (.csv)                 — the counters, as name/value pairs,
* ``*_dem`` / ``*_channels`` / ``*_ridge_side`` (.tif) — the rasters.

The product CLI writes FlatGeobuf by default; the Phase 0/1 development driver
writes the same tables as CSV.  Both are accepted, `.fgb` first.

No PyGMT imports here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

#: Vector extensions written by ``vector_output_format``, in preference order.
_VECTOR_SUFFIXES = (".fgb", ".gpkg", ".shp", ".csv")

#: Trace outcome codes, as emitted by ``lsdtt::streamline_outcome``.
OUTCOME_NAME = {
    0: "reached channel",
    1: "hit edge",
    2: "loop detected",
    3: "no aspect",
}

#: Stable colours so the same outcome reads the same way in every figure.
OUTCOME_COLOUR = {
    0: "#2a7f3f",
    1: "#d62728",
    2: "#9467bd",
    3: "#ff7f0e",
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
    """Read one table, whatever container it arrived in."""
    if path.suffix == ".csv":
        return pd.read_csv(path)

    import geopandas as gpd  # imported lazily: CSV-only users need no GDAL

    return gpd.read_file(path)


@dataclass
class HilltopCouplingOutputs:
    """One hilltop-flow-routing output set, loaded.

    Attributes:
        hilltops: One row per hilltop pixel, carrying the ordering attributes
            (``ridgeline_id``, ``ridgeline_node_index``,
            ``distance_along_ridgeline_m``, ``ridge_side``).
        ridgelines: The ordered divide nodes, long-form, one row per node.
        traces: One row per trace (``length_m``, ``relief_m``, ``outcome``).
        paths: Thinned trace geometry, long-form, one row per vertex.
        summary: The run counters as a name -> value mapping.
        dem_path, channels_path, ridge_side_path: Raster paths, when present.
    """

    hilltops: Optional[pd.DataFrame] = None
    ridgelines: Optional[pd.DataFrame] = None
    traces: Optional[pd.DataFrame] = None
    paths: Optional[pd.DataFrame] = None
    summary: dict = field(default_factory=dict)

    dem_path: Optional[Path] = None
    channels_path: Optional[Path] = None
    ridge_side_path: Optional[Path] = None

    def outcome_counts(self) -> pd.Series:
        """Traces per outcome code, ascending by code.

        Returns:
            A Series indexed by outcome code; empty if no traces were loaded.
        """
        if self.traces is None or "outcome" not in self.traces:
            return pd.Series(dtype=int)
        return self.traces["outcome"].value_counts().sort_index()

    def successful_traces(self) -> pd.DataFrame:
        """The traces that reached a channel.

        Every other outcome leaves ``length_m`` meaningless, so any distribution
        of hillslope length must be taken over this subset and the discarded
        count reported alongside it.
        """
        if self.traces is None:
            return pd.DataFrame()
        return self.traces[self.traces["outcome"] == 0]

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
    for attr, needle in (
        ("dem_path", "dem"),
        ("channels_path", "channels"),
        ("ridge_side_path", "ridge_side"),
    ):
        for candidate in sorted(search_dir.glob(f"*{needle}.tif")):
            setattr(out, attr, candidate)
            break

    return out
