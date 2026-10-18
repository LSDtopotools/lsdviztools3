"""Lithology IO for lsdtt3 outputs: unit names, attribution, k_sn-by-unit tables.

Readers for the products of ``lsdtt-fetch-lithology``,
``lsdtt-lithology-erodibility`` and the lithology option of
``lsdtt-chi-analysis`` (``lithology_fname`` + ``lithology_field_name``):

* ``<prefix>_lithology_code.tif``        — run-local integer unit code raster.
* ``<prefix>_lithology_units.csv``       — ``code,name,map_values,...``.
* ``<w>_chi_ksn_by_lithology_mn_<mn>.csv``           — per-unit summary.
* ``<w>_chi_ksn_histogram_by_lithology_mn_<mn>.csv`` — common-bin histogram.
* ``<w>_chi_ksn_lithology_nodes_mn_<mn>.csv``        — per-node segment k_sn.
* ``<w>_chi_lithology_reaches_mn_<mn>.csv``          — per-reach fits.
* ``<w>_chi_ksn_by_lithology_and_basin_mn_<mn>.csv`` — per unit x basin.
* ``<w>_chi_lithology_segments_mn_<mn>.csv``         — per-segment majority unit.
* ``<prefix>_ATTRIBUTION.txt``           — geologic-map citation (fetch-lithology).

Every table is read **by header**: lsdtt3 appends columns over time, so nothing
here relies on column order, and the column names the charts need are kept in
the ``*_COLUMNS`` constants below (the one place to edit if lsdtt3 renames
them). Optional tables that are absent (older runs) return ``None``.

Nothing here computes a statistic. Names, sample sizes, medians, quantiles and
confidence intervals all come from the CSVs lsdtt3 wrote. No PyGMT imports.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, Optional

import pandas as pd

from .table import load_table

# --------------------------------------------------------------------------
# Column names (the single place to change if lsdtt3 renames a column)
# --------------------------------------------------------------------------

#: Code -> name columns, in every lsdtt3 lithology table that carries them.
CODE_COLUMN = "lith_code"
NAME_COLUMN = "lith_name"

#: Code -> name columns of a standalone legend / units / crosswalk CSV.
LEGEND_CODE_COLUMN = "code"
LEGEND_NAME_COLUMN = "name"
#: Optional colour column in a legend / units CSV (any GMT colour string).
LEGEND_COLOUR_COLUMN = "color"

#: Reserved codes written by lsdtt-chi-analysis.
UNMAPPED_CODE = -1
UNMAPPED_NAME = "unmapped"
ALL_UNITS_CODE = -2

#: Per-unit summary (``*_chi_ksn_by_lithology_mn_*.csv``).
SUMMARY_COLUMNS = {
    "n_nodes": "n_nodes",
    "n_segments": "n_segments",
    "n_basins": "n_basins",
    "channel_length_m": "channel_length_m",
    "median": "k_sn_median",
    "q25": "k_sn_q25",
    "q75": "k_sn_q75",
    "mean": "k_sn_mean",
    "min": "k_sn_min",
    "max": "k_sn_max",
    "low_sample": "low_sample",
    "m_over_n": "m_over_n",
    "A_0": "A_0",
    # Contact buffer / reach-level columns (appended by lsdtt3; optional so
    # first-cut files still load). low_sample = fewer than lsdtt3's minimum
    # number of QUALIFYING reaches; the reach CI is NaN for those units.
    "contact_buffer_m": "contact_buffer_metres",
    "min_reach_pixels": "min_reach_pixels",
    "n_excluded_near_contact": "n_nodes_excluded_near_contact",
    "n_reaches": "n_reaches",
    "n_reaches_qualifying": "n_reaches_qualifying",
    "median_drainage_area": "median_drainage_area_m2",
    "reach_median": "reach_k_sn_median",
    "reach_p25": "reach_k_sn_p25",
    "reach_p75": "reach_k_sn_p75",
    "reach_lw_median": "reach_k_sn_length_weighted_median",
    "reach_ci_low": "reach_k_sn_median_ci_low",
    "reach_ci_high": "reach_k_sn_median_ci_high",
    "n_bootstrap": "n_bootstrap",
    "bootstrap_cluster": "bootstrap_cluster",
    "node_lw_median": "node_k_sn_length_weighted_median",
    "node_lw_p25": "node_k_sn_length_weighted_p25",
    "node_lw_p75": "node_k_sn_length_weighted_p75",
}

#: Common-bin histogram (``*_chi_ksn_histogram_by_lithology_mn_*.csv``).
HISTOGRAM_COLUMNS = {
    "bin_index": "bin_index",
    "bin_lower": "bin_lower",
    "bin_upper": "bin_upper",
    "n_nodes": "n_nodes",
    "channel_length_m": "channel_length_m",
}

#: Per-node table (``*_chi_ksn_lithology_nodes_mn_*.csv``).
NODE_COLUMNS = {
    "basin": "basin_index",
    "x": "x",
    "y": "y",
    "k_sn": "k_sn",
    "segment_id": "segment_id",
    "distance_to_contact_m": "distance_to_contact_m",
    "latitude": "latitude",
    "longitude": "longitude",
    "drainage_area_m2": "drainage_area_m2",
    "branch_id": "branch_id",
    "near_contact": "near_contact",
    "reach_id": "reach_id",
}

#: Per-reach table (``*_chi_lithology_reaches_mn_*.csv``). Only rows with
#: ``qualifies == 1`` belong in reach-level distributions.
REACH_COLUMNS = {
    "reach_id": "reach_id",
    "basin": "basin_index",
    "branch_id": "branch_id",
    "n_nodes": "n_nodes",
    "k_sn": "k_sn_ols",
    "k_sn_se": "k_sn_ols_se",
    "r_squared": "r_squared",
    "k_sn_chi_weighted": "k_sn_fit_chi_weighted",
    "length_m": "length_metres",
    "upstream_end": "upstream_end",
    "downstream_end": "downstream_end",
    "qualifies": "qualifies",
}

#: Per unit x basin table (``*_chi_ksn_by_lithology_and_basin_mn_*.csv``).
#: The unit x basin view plots the REACH median, not the node-level median.
UNIT_BASIN_COLUMNS = {
    "basin": "basin_index",
    "median": "reach_k_sn_median",
    "node_median": "k_sn_median",
    "n_nodes": "n_nodes",
    "n_reaches": "n_reaches",
    "n_reaches_qualifying": "n_reaches_qualifying",
    "low_sample": "low_sample",
}

#: Per-segment lithology table (``*_chi_lithology_segments_mn_*.csv``).
SEGMENT_COLUMNS = {
    "basin": "basin_index",
    "segment_id": "segment_id",
    "k_sn": "k_sn",
    "length_m": "length_metres",
    "majority_fraction": "lith_majority_fraction",
    "n_units": "n_lith_units",
}

#: Units table (``<w>_lithology_units.csv``): code,name,map_values,...
UNITS_COLUMNS = {
    "code": "code",
    "name": "name",
    "map_values": "map_values",
    "n_pixels": "n_pixels",
    "area_km2": "area_km2",
    "n_pixels_in_analysed_basins": "n_pixels_in_analysed_basins",
}

#: File-name suffixes (after ``<prefix>``; ``{mn}`` is the m/n suffix).
SUFFIXES = {
    "summary": "_chi_ksn_by_lithology_mn_{mn}.csv",
    "histogram": "_chi_ksn_histogram_by_lithology_mn_{mn}.csv",
    "nodes": "_chi_ksn_lithology_nodes_mn_{mn}.csv",
    "reaches": "_chi_lithology_reaches_mn_{mn}.csv",
    "unit_basin": "_chi_ksn_by_lithology_and_basin_mn_{mn}.csv",
    "segments": "_chi_lithology_segments_mn_{mn}.csv",
    "units": "_lithology_units.csv",
    "code_raster": "_lithology_code.tif",
    "legend": "_lithology_legend.csv",
    "attribution": "_ATTRIBUTION.txt",
    "ingested_params": "_ingested_params.param",
}


def mn_suffix(m_over_n: float) -> str:
    """Filename suffix lsdtt3 uses for an m/n value (``0.45 -> "0p450"``).

    Mirrors lsdtt-chi-analysis' ``mn_suffix``: three decimals with ``.`` -> ``p``
    and a leading ``-`` -> ``m``.
    """
    return f"{m_over_n:.3f}".replace(".", "p").replace("-", "m")


def find_mn_suffixes(directory, prefix: str) -> list:
    """List the m/n suffixes present for ``prefix`` (from the summary CSV names)."""
    pattern = re.compile(
        re.escape(prefix) + r"_chi_ksn_by_lithology_mn_([0-9pm]+)\.csv$")
    out = []
    for p in sorted(Path(directory).glob(f"{prefix}_chi_ksn_by_lithology_mn_*.csv")):
        m = pattern.match(p.name)
        if m:
            out.append(m.group(1))
    return out


def lithology_paths(directory, prefix: str, mn: Optional[str] = None) -> Dict[str, Path]:
    """Expected paths of the lithology outputs for one lsdtt3 run.

    Args:
        directory: Directory holding the outputs.
        prefix: The run's ``write_fname`` (e.g. ``"wvaf_lithprimary"``).
        mn: m/n suffix (e.g. ``"0p450"``) or a float m/n. ``None`` picks the
            only suffix present; raises if there are several.

    Returns:
        ``{kind: Path}`` for every kind in :data:`SUFFIXES` (files may not exist).
    """
    directory = Path(directory)
    if mn is None:
        found = find_mn_suffixes(directory, prefix)
        if len(found) > 1:
            raise ValueError(
                f"Several m/n values for {prefix!r} in {directory}: {found}; "
                "pass one explicitly.")
        mn = found[0] if found else ""
    elif not isinstance(mn, str):
        mn = mn_suffix(float(mn))
    return {k: directory / f"{prefix}{s.format(mn=mn)}" for k, s in SUFFIXES.items()}


def _read_optional(path) -> Optional[pd.DataFrame]:
    """Read a CSV if it exists, else ``None``."""
    if path is None or not Path(path).exists():
        return None
    return load_table(path)


def _require(df: pd.DataFrame, columns: Iterable[str], what: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"{what} is missing column(s) {missing}; has {list(df.columns)}")


# --------------------------------------------------------------------------
# Table readers
# --------------------------------------------------------------------------

def load_ksn_lithology_summary(path) -> pd.DataFrame:
    """Load ``*_chi_ksn_by_lithology_mn_*.csv`` (one row per unit, plus -1/-2)."""
    df = load_table(path)
    _require(df, (CODE_COLUMN, NAME_COLUMN, SUMMARY_COLUMNS["n_nodes"],
                  SUMMARY_COLUMNS["median"]), Path(path).name)
    return df


def load_ksn_lithology_histogram(path) -> pd.DataFrame:
    """Load ``*_chi_ksn_histogram_by_lithology_mn_*.csv``.

    The open edges of the outermost bins are written as ``-inf``/``inf``; they
    are parsed to float infinities here (charts treat them as overflow).
    """
    df = load_table(path)
    _require(df, (CODE_COLUMN, NAME_COLUMN, *HISTOGRAM_COLUMNS.values()),
             Path(path).name)
    for c in (HISTOGRAM_COLUMNS["bin_lower"], HISTOGRAM_COLUMNS["bin_upper"]):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_ksn_lithology_nodes(path) -> pd.DataFrame:
    """Load ``*_chi_ksn_lithology_nodes_mn_*.csv`` (per-node segment k_sn)."""
    df = load_table(path)
    _require(df, (CODE_COLUMN, NAME_COLUMN, NODE_COLUMNS["k_sn"]), Path(path).name)
    return df


def load_lithology_reaches(path) -> Optional[pd.DataFrame]:
    """Load ``*_chi_lithology_reaches_mn_*.csv`` if present."""
    df = _read_optional(path)
    if df is not None:
        _require(df, (CODE_COLUMN, NAME_COLUMN, REACH_COLUMNS["k_sn"]), Path(path).name)
    return df


def load_ksn_lithology_by_basin(path) -> Optional[pd.DataFrame]:
    """Load ``*_chi_ksn_by_lithology_and_basin_mn_*.csv`` if present."""
    df = _read_optional(path)
    if df is not None:
        _require(df, (CODE_COLUMN, NAME_COLUMN, UNIT_BASIN_COLUMNS["basin"],
                      UNIT_BASIN_COLUMNS["median"]), Path(path).name)
    return df


def load_lithology_segments(path) -> Optional[pd.DataFrame]:
    """Load ``*_chi_lithology_segments_mn_*.csv`` if present."""
    df = _read_optional(path)
    if df is not None:
        _require(df, (CODE_COLUMN, NAME_COLUMN, SEGMENT_COLUMNS["k_sn"]), Path(path).name)
    return df


def load_ksn_lithology_outputs(directory, prefix: str, mn=None) -> dict:
    """Load every lithology k_sn table for one run that exists.

    Returns:
        ``{"summary", "histogram", "nodes", "reaches", "unit_basin", "segments",
        "m_over_n", "A_0", "contact_buffer_m", "paths"}``. Missing tables are
        ``None``. ``m_over_n``, ``A_0`` and ``contact_buffer_m`` are read from
        the summary file.
    """
    paths = lithology_paths(directory, prefix, mn)
    if not paths["summary"].exists():
        raise FileNotFoundError(f"No lithology summary at {paths['summary']}")
    summary = load_ksn_lithology_summary(paths["summary"])
    out = {
        "summary": summary,
        "histogram": (load_ksn_lithology_histogram(paths["histogram"])
                      if paths["histogram"].exists() else None),
        "nodes": (load_ksn_lithology_nodes(paths["nodes"])
                  if paths["nodes"].exists() else None),
        "reaches": load_lithology_reaches(paths["reaches"]),
        "unit_basin": load_ksn_lithology_by_basin(paths["unit_basin"]),
        "segments": load_lithology_segments(paths["segments"]),
        "paths": paths,
    }
    out["m_over_n"] = _first_value(summary, SUMMARY_COLUMNS["m_over_n"])
    out["A_0"] = _first_value(summary, SUMMARY_COLUMNS["A_0"])
    out["contact_buffer_m"] = _first_value(summary, SUMMARY_COLUMNS["contact_buffer_m"])
    return out


def _first_value(df: pd.DataFrame, column: str) -> Optional[float]:
    if column not in df.columns or df.empty:
        return None
    v = pd.to_numeric(df[column], errors="coerce").dropna()
    return float(v.iloc[0]) if len(v) else None


# --------------------------------------------------------------------------
# Unit names
# --------------------------------------------------------------------------

def _pairs(df: pd.DataFrame, code_col: str, name_col: str) -> Dict[int, str]:
    sub = df[[code_col, name_col]].dropna()
    out: Dict[int, str] = {}
    for code, name in zip(sub[code_col], sub[name_col]):
        try:
            c = int(code)
        except (TypeError, ValueError):
            continue
        out.setdefault(c, str(name))
    return out


def read_code_names(path) -> Dict[int, str]:
    """Read a ``code -> name`` mapping from any lsdtt3 lithology CSV.

    Accepts a units/legend/crosswalk CSV (``code,name`` columns) or a chi
    lithology table (``lith_code,lith_name`` columns), whichever header is
    present. Rows without a usable integer code are skipped.
    """
    df = load_table(path)
    if LEGEND_CODE_COLUMN in df.columns and LEGEND_NAME_COLUMN in df.columns:
        return _pairs(df, LEGEND_CODE_COLUMN, LEGEND_NAME_COLUMN)
    if CODE_COLUMN in df.columns and NAME_COLUMN in df.columns:
        return _pairs(df, CODE_COLUMN, NAME_COLUMN)
    raise KeyError(
        f"{Path(path).name}: no code/name columns (expected "
        f"{LEGEND_CODE_COLUMN},{LEGEND_NAME_COLUMN} or {CODE_COLUMN},{NAME_COLUMN})")


def read_code_colours(path) -> Dict[int, str]:
    """Optional ``code -> colour`` from a units/legend CSV with a colour column."""
    df = load_table(path)
    if LEGEND_COLOUR_COLUMN not in df.columns or LEGEND_CODE_COLUMN not in df.columns:
        return {}
    return _pairs(df, LEGEND_CODE_COLUMN, LEGEND_COLOUR_COLUMN)


def name_sources(directory, prefix: str, *, crosswalk=None, legend=None) -> list:
    """Candidate code->name files for ``prefix``, in priority order.

    1. ``<prefix>_lithology_units.csv`` (``code,name,...``);
    2. the chi k_sn CSVs (``lith_code,lith_name``), summary first;
    3. an explicit crosswalk CSV (``code,name``);
    4. ``lsdtt-lithology-erodibility``'s ``<prefix>_lithology_legend.csv`` (or
       an explicit ``legend`` path).

    Only files that exist are returned.
    """
    directory = Path(directory)
    cands = [directory / f"{prefix}{SUFFIXES['units']}"]
    for kind in ("summary", "histogram", "nodes"):
        pat = SUFFIXES[kind].format(mn="*")
        cands += sorted(directory.glob(f"{prefix}{pat}"))
    if crosswalk is not None:
        cands.append(Path(crosswalk))
    cands.append(Path(legend) if legend is not None
                 else directory / f"{prefix}{SUFFIXES['legend']}")
    return [p for p in cands if p.exists()]


def resolve_code_names(directory=None, prefix: Optional[str] = None, *,
                       names_csv=None, crosswalk=None, legend=None) -> Dict[int, str]:
    """Resolve lithology code -> unit name, by the documented priority order.

    An explicit ``names_csv`` wins outright. Otherwise the first existing file
    from :func:`name_sources` that yields a mapping is used. Code ``-1`` is
    always ``"unmapped"``.

    Raises:
        FileNotFoundError: if no source is found — a legend of bare integer
            codes is never produced silently.
    """
    sources = []
    if names_csv is not None:
        sources.append(Path(names_csv))
    if directory is not None and prefix is not None:
        sources += name_sources(directory, prefix, crosswalk=crosswalk, legend=legend)
    elif crosswalk is not None or legend is not None:
        sources += [Path(p) for p in (crosswalk, legend) if p is not None]

    for src in sources:
        try:
            mapping = read_code_names(src)
        except KeyError:
            continue
        if mapping:
            mapping[UNMAPPED_CODE] = UNMAPPED_NAME
            mapping.pop(ALL_UNITS_CODE, None)
            return mapping
    raise FileNotFoundError(
        f"No lithology code->name source found for prefix={prefix!r} in "
        f"{directory}; pass a units/legend/crosswalk CSV.")


def read_unit_colours(path) -> Dict[str, str]:
    """Read a user ``name -> colour`` mapping (CSV or JSON).

    CSV: columns ``name`` and ``colour`` (or ``color``). JSON: an object
    ``{"name": "colour", ...}``. Colours are any matplotlib name, ``#rrggbb``
    or GMT ``r/g/b``. Keyed by unit NAME because codes are run-local.
    """
    import json

    path = Path(path)
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError(f"{path.name}: expected a JSON object name -> colour")
        return {str(k): str(v) for k, v in data.items()}
    df = load_table(path)
    ccol = next((c for c in ("colour", "color") if c in df.columns), None)
    if LEGEND_NAME_COLUMN not in df.columns or ccol is None:
        raise KeyError(f"{path.name}: needs 'name' and 'colour' (or 'color') columns")
    sub = df[[LEGEND_NAME_COLUMN, ccol]].dropna()
    return {str(n): str(c) for n, c in zip(sub[LEGEND_NAME_COLUMN], sub[ccol])}


# --------------------------------------------------------------------------
# Parameter files and attribution
# --------------------------------------------------------------------------

def read_param_file(path) -> Dict[str, str]:
    """Parse an lsdtt3 ``key: value`` parameter file (comments ``#`` ignored)."""
    out: Dict[str, str] = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        key, _, value = line.partition(":")
        out[key.strip()] = value.strip()
    return out


#: Human-readable names of DEM products, keyed by lsdtt-fetch-raster settings.
DEM_PRODUCT_NAMES = {
    "cop30_aws": "Copernicus GLO-30",
    "cop-dem-glo-30": "Copernicus GLO-30",
    "cop-dem-glo-90": "Copernicus GLO-90",
    "COP30": "Copernicus GLO-30",
    "COP90": "Copernicus GLO-90",
    "nasadem": "NASADEM",
    "NASADEM": "NASADEM",
    "SRTMGL1": "SRTM GL1",
    "SRTMGL3": "SRTM GL3",
    "SRTMGL1_E": "SRTM GL1 (ellipsoidal)",
    "AW3D30": "ALOS AW3D30",
    "AW3D30_E": "ALOS AW3D30 (ellipsoidal)",
    "EU_DTM": "EU DTM",
    "GEDI_L3": "GEDI L3",
}


def dem_source_from_params(path) -> Optional[str]:
    """DEM product name recorded in an ``lsdtt-fetch-raster`` ingested-params file.

    Reads ``dem_source`` and, for providers that host several products, the
    product key (``pc_collection`` / ``ot_dem_type``). Returns ``None`` if the
    file records no ``dem_source`` (e.g. it is not a fetch run).
    """
    params = read_param_file(path)
    source = params.get("dem_source")
    if not source:
        return None
    if source == "planetary_computer":
        key = params.get("pc_collection", "")
    elif source == "opentopography":
        key = params.get("ot_dem_type", "")
    else:
        key = source
    return DEM_PRODUCT_NAMES.get(key, key or source)


def discover_dem_params(dem_path) -> Optional[Path]:
    """``<stem>_ingested_params.param`` beside a fetched ``<stem>_DEM.tif``, if any."""
    dem_path = Path(dem_path)
    stem = re.sub(r"_DEM$", "", dem_path.stem, flags=re.IGNORECASE)
    cand = dem_path.with_name(f"{stem}{SUFFIXES['ingested_params']}")
    return cand if cand.exists() else None


def discover_attribution(directory, prefix: Optional[str] = None,
                         params=None) -> Optional[Path]:
    """Find the ``lsdtt-fetch-lithology`` ``*_ATTRIBUTION.txt`` for a run.

    Looks at ``<prefix>_ATTRIBUTION.txt`` first, then follows the chi run's
    ``lithology_fname`` (from ``params`` or ``<prefix>_ingested_params.param``):
    ``<stem>_macrostrat_units.fgb`` -> ``<stem>_ATTRIBUTION.txt``.
    """
    directory = Path(directory)
    if prefix:
        direct = directory / f"{prefix}{SUFFIXES['attribution']}"
        if direct.exists():
            return direct
        if params is None:
            cand = directory / f"{prefix}{SUFFIXES['ingested_params']}"
            params = cand if cand.exists() else None
    if params is None:
        return None
    lith = read_param_file(params).get("lithology_fname")
    if not lith:
        return None
    stem = re.sub(r"(_macrostrat_units)?\.[^.]+$", "", Path(lith).name)
    cand = directory / f"{stem}{SUFFIXES['attribution']}"
    return cand if cand.exists() else None


#: Publisher abbreviations used when condensing a source map's series line.
PUBLISHER_ABBREVIATIONS = {
    "U.S. Geological Survey": "USGS",
    "US Geological Survey": "USGS",
    "British Geological Survey": "BGS",
    "Geological Survey of Canada": "GSC",
    "Geoscience Australia": "GA",
}


def parse_attribution(text: str) -> dict:
    """Parse the text of an ``lsdtt-fetch-lithology`` ``_ATTRIBUTION.txt``.

    Returns:
        ``{"provider", "licence", "sources": [{"name", "authors", "year",
        "title", "series", "doi"}], "text"}``. Fields not found are ``None``.
    """
    out = {"provider": None, "licence": None, "sources": [], "text": text}
    first = text.strip().splitlines()[0] if text.strip() else ""
    m = re.search(r"from\s+([^(,]+?)\s*\(", first)
    if m:
        out["provider"] = m.group(1).strip()
    # Licence runs to the end of its line, minus a sentence-final full stop
    # (so "CC BY 4.0." -> "CC BY 4.0").
    m = re.search(r"licen[sc]ed\s+(.+?)\.?\s*$", text, flags=re.MULTILINE)
    if m:
        out["licence"] = m.group(1).strip()

    blocks = re.split(r"^\s*-\s*source_id\s+\S+?:\s*", text, flags=re.MULTILINE)[1:]
    for block in blocks:
        lines = [ln.strip() for ln in block.strip().splitlines() if ln.strip()]
        if not lines:
            continue
        src = {"name": lines[0], "authors": None, "year": None,
               "title": None, "series": None, "doi": None}
        rest = []
        for ln in lines[1:]:
            if ln.lower().startswith("doi"):
                src["doi"] = ln.split(":", 1)[-1].strip()
            elif ln.lower().startswith("http"):
                continue
            else:
                rest.append(ln)
        if rest:
            ym = re.search(r"\((\d{4})\)", rest[0])
            if ym:
                src["year"] = ym.group(1)
                src["authors"] = rest[0][:ym.start()].strip()
            else:
                src["authors"] = rest[0]
        if len(rest) > 1:
            src["title"] = rest[1]
        if len(rest) > 2:
            src["series"] = rest[2]
        out["sources"].append(src)
    return out


def read_attribution(path) -> dict:
    """Read and parse an ``_ATTRIBUTION.txt`` (see :func:`parse_attribution`)."""
    return parse_attribution(Path(path).read_text(encoding="utf-8", errors="replace"))


def _short_authors(authors: Optional[str]) -> Optional[str]:
    if not authors:
        return None
    first = authors.split(",")[0].strip()
    first = first.split(" and ")[0].strip()
    many = ("," in authors) or (" and " in authors) or ("&" in authors)
    return f"{first} et al." if many else first


def _short_source(src: dict) -> str:
    """``"USGS SGMC, Horton et al. 2017"`` from one parsed source block."""
    name = src.get("name") or ""
    acr = re.search(r"\(([A-Z][A-Za-z0-9]{1,10})\)", src.get("title") or "")
    label = acr.group(1) if acr else name
    series = src.get("series") or ""
    for full, abbr in PUBLISHER_ABBREVIATIONS.items():
        if full in series or full in name:
            label = f"{abbr} {label}"
            break
    who = _short_authors(src.get("authors"))
    tail = " ".join(v for v in (who, src.get("year")) if v)
    return f"{label}, {tail}" if tail else label


def condense_attribution(attribution: dict, *, dem_source: Optional[str] = None,
                         max_sources: int = 3) -> str:
    """Condense a parsed attribution into a one-line map credit.

    e.g. ``"Geology: Macrostrat (CC BY 4.0); USGS SGMC, Horton et al. 2017.
    DEM: Copernicus GLO-30"``. Beyond ``max_sources`` source maps the line ends
    ``"+ N more"`` (the full list stays in the ATTRIBUTION file).
    """
    parts = []
    provider = attribution.get("provider")
    if provider:
        lic = attribution.get("licence")
        parts.append(f"{provider} ({lic})" if lic else provider)
    sources = attribution.get("sources") or []
    parts += [_short_source(s) for s in sources[:max_sources]]
    if len(sources) > max_sources:
        parts.append(f"+ {len(sources) - max_sources} more")
    line = ("Geology: " + "; ".join(parts)) if parts else ""
    if dem_source:
        line = f"{line}. DEM: {dem_source}" if line else f"DEM: {dem_source}"
    return line


def lithology_map_credit(*, attribution=None, credit: Optional[str] = None,
                         dem_source: Optional[str] = None,
                         dem_params=None) -> str:
    """Build the map credit line from whichever inputs are supplied.

    ``credit`` (an explicit string) wins outright. Otherwise the attribution
    file (path or parsed dict) is condensed, and the DEM product is taken from
    ``dem_source`` or, failing that, from the fetch run's ``dem_params`` file.
    """
    if credit:
        return credit
    if dem_source is None and dem_params is not None:
        dem_source = dem_source_from_params(dem_params)
    parsed = {}
    if attribution is not None:
        parsed = attribution if isinstance(attribution, dict) else read_attribution(attribution)
    return condense_attribution(parsed, dem_source=dem_source)
