"""Erodibility-calibration IO for lsdtt3 outputs (K from basin erosion rates).

Readers for the outputs of ``lsdtt-chi-analysis`` with ``erosion_rate_fname``
(see the lsdtt3 chi reference, "Erodibility calibration from basin erosion
rates":
https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_chi_analysis.html#erodibility-calibration-from-basin-erosion-rates):

* ``<w>_erodibility_calibration_samples_mn_<mn>.csv``   — one row per sample;
* ``<w>_erodibility_calibration_basin_K_mn_<mn>.csv``   — per sample and n;
* ``<w>_erodibility_calibration_fit_mn_<mn>.csv``       — fixed-n rows and the
  fitted (York) row;
* ``<w>_erodibility_calibration_by_lithology_mn_<mn>.csv`` — per unit;
* ``<w>_erodibility_calibration_basin_lithology.csv``   — unit fractions per basin;
* ``<w>_erodibility_calibration_basins.fgb``            — basin polygons.

lsdtt3's model is E = K A0^m k_sn^n with m = (m/n) n, K in m^(1 - 2m) yr^-1.
Column names live in the ``*_COLUMNS`` constants; tables are read by header.
Nothing here computes a statistic. No PyGMT imports.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from .lithology import mn_suffix, read_param_file
from .table import load_table

#: Samples table.
SAMPLE_COLUMNS = {
    "sample_id": "sample_id",
    "status": "status",
    "status_detail": "status_detail",
    "snapped_lat": "snapped_latitude",
    "snapped_lon": "snapped_longitude",
    "nesting_group": "nesting_group",
    "E": "erosion_rate_m_per_yr",
    "E_sigma": "erosion_rate_uncertainty_m_per_yr",
    "m_over_n": "m_over_n",
    "A_0": "A_0",
    "k_sn_statistic": "k_sn_statistic",
    "k_sn": "k_sn",
    "k_sn_ci_low": "k_sn_ci_low",
    "k_sn_ci_high": "k_sn_ci_high",
    "dominant_lith_name": "dominant_lith_name",
    "dominant_area_fraction": "dominant_area_fraction",
}

#: Fit table (one fixed_n row per n, one fitted_n row).
FIT_COLUMNS = {
    "fit_type": "fit_type",
    "method": "method",
    "status": "status",
    "N_basins": "N_basins",
    "N_groups": "N_groups",
    "m_over_n": "m_over_n",
    "A_0": "A_0",
    "n": "n",
    "n_ci_low": "n_ci_low",
    "n_ci_high": "n_ci_high",
    "m": "m",
    "K": "K",
    "K_ci_low": "K_ci_low",
    "K_ci_high": "K_ci_high",
    "corr_n_lnK": "corr_n_lnK",
    "mswd": "mswd",
    "n_bootstrap": "n_bootstrap",
    "bootstrap_cluster": "bootstrap_cluster",
    "K_units": "K_units",
}

#: Basin polygons (fgb): the status there is the pre-k_sn one ("kept").
BASIN_COLUMNS = {"sample_id": "sample_id", "status": "status",
                 "nesting_group": "nesting_group"}

#: The status that enters the fits.
USED_STATUS = "used"

SUFFIXES = {
    "samples": "_erodibility_calibration_samples_mn_{mn}.csv",
    "basin_K": "_erodibility_calibration_basin_K_mn_{mn}.csv",
    "fit": "_erodibility_calibration_fit_mn_{mn}.csv",
    "by_lithology": "_erodibility_calibration_by_lithology_mn_{mn}.csv",
    "basin_lithology": "_erodibility_calibration_basin_lithology.csv",
    "basins": "_erodibility_calibration_basins.fgb",
    "basins_gpkg": "_erodibility_calibration_basins.gpkg",
}


def calibration_paths(directory, prefix: str, mn) -> Dict[str, Path]:
    """Expected output paths of one calibration run (files may not exist)."""
    directory = Path(directory)
    if not isinstance(mn, str):
        mn = mn_suffix(float(mn))
    return {k: directory / f"{prefix}{s.format(mn=mn)}" for k, s in SUFFIXES.items()}


def find_calibration_mn(directory, prefix: str) -> list:
    """m/n suffixes present for ``prefix`` (from the samples file names)."""
    pat = re.compile(re.escape(prefix) + r"_erodibility_calibration_samples_mn_([0-9pm]+)\.csv$")
    return [m.group(1) for p in sorted(Path(directory).glob(
        f"{prefix}_erodibility_calibration_samples_mn_*.csv")) if (m := pat.match(p.name))]


def _require(df, cols, what):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"{what} is missing column(s) {missing}; has {list(df.columns)}")


def load_calibration_samples(path) -> pd.DataFrame:
    """Load ``*_erodibility_calibration_samples_mn_*.csv``."""
    df = load_table(path)
    c = SAMPLE_COLUMNS
    _require(df, (c["sample_id"], c["status"], c["E"], c["k_sn"]), Path(path).name)
    df[c["sample_id"]] = df[c["sample_id"]].astype(str)
    return df


def load_calibration_fit(path) -> pd.DataFrame:
    """Load ``*_erodibility_calibration_fit_mn_*.csv``."""
    df = load_table(path)
    c = FIT_COLUMNS
    _require(df, (c["fit_type"], c["status"], c["n"], c["K"]), Path(path).name)
    return df


def load_calibration_outputs(directory, prefix: str, mn=None) -> dict:
    """Load whichever calibration outputs exist for ``prefix`` at one m/n.

    Returns ``{"samples", "fit", "basin_K", "by_lithology", "basin_lithology",
    "basins", "paths", "mn"}`` (missing tables ``None``; ``basins`` a
    GeoDataFrame).
    """
    if mn is None:
        found = find_calibration_mn(directory, prefix)
        if len(found) != 1:
            raise ValueError(f"Expected one m/n for {prefix!r}, found {found}; pass one.")
        mn = found[0]
    paths = calibration_paths(directory, prefix, mn)
    if not paths["samples"].exists():
        raise FileNotFoundError(f"No calibration samples at {paths['samples']}")
    out = {"samples": load_calibration_samples(paths["samples"]),
           "fit": load_calibration_fit(paths["fit"]) if paths["fit"].exists() else None,
           "paths": paths, "mn": mn}
    for k in ("basin_K", "by_lithology", "basin_lithology"):
        out[k] = load_table(paths[k]) if paths[k].exists() else None
    poly = paths["basins"] if paths["basins"].exists() else (
        paths["basins_gpkg"] if paths["basins_gpkg"].exists() else None)
    if poly is not None:
        from .vector import load_vector

        g = load_vector(poly)
        g[BASIN_COLUMNS["sample_id"]] = g[BASIN_COLUMNS["sample_id"]].astype(str)
        out["basins"] = g
    else:
        out["basins"] = None
    return out


# --------------------------------------------------------------------------
# Credits
# --------------------------------------------------------------------------

def discover_erosion_rate_attribution(directory, prefix: str) -> Optional[Path]:
    """``<stem>_ATTRIBUTION.txt`` of the erosion-rate file named in the chi run's
    ingested params (``erosion_rate_fname: <stem>_erosion_rates.csv``)."""
    directory = Path(directory)
    params = directory / f"{prefix}_ingested_params.param"
    if not params.exists():
        return None
    fname = read_param_file(params).get("erosion_rate_fname")
    if not fname:
        return None
    stem = re.sub(r"(_erosion_rates)?\.[^.]+$", "", Path(fname).name)
    cand = directory / f"{stem}_ATTRIBUTION.txt"
    return cand if cand.exists() else None


def erosion_rate_credit(path) -> str:
    """Short credit from an ``lsdtt-fetch-erosion-rates`` ATTRIBUTION file.

    e.g. ``"Erosion rates: OCTOPUS (CC BY 4.0); Duxburry 2009; Portenga 2019"``
    — the provider and licence from the header, and the original studies
    (author, year) from the "Original studies" block.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    m = re.search(r"from\s+([^(,\n]+?)\s*\(", text)
    provider = m.group(1).strip() if m else None
    m = re.search(r"licen[sc]ed\s+(.+?)\.?\s*$", text, flags=re.MULTILINE)
    licence = m.group(1).strip() if m else None
    parts = []
    if provider:
        parts.append(f"{provider} ({licence})" if licence else provider)
    block = re.search(r"Original studies[^\n]*\n((?:[ \t]+\S[^\n]*\n?)+)", text)
    if block:
        for line in block.group(1).splitlines():
            f = [s.strip() for s in line.split(",")]
            if len(f) >= 3 and f[2].isdigit():
                parts.append(f"{f[1]} {f[2]}")
    return ("Erosion rates: " + "; ".join(parts)) if parts else ""
