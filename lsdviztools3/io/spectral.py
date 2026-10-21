"""IO for ``lsdtt-spectral-analysis`` 2D significant-frequency output.

``<stem>_significant_frequencies_2d.csv`` (``write_significant_frequencies_2d:
true``, the default) lists every cell of the 2D power spectrum whose normalised
PSD exceeds ``-ln(1 - significant_frequencies_2d_confidence)``. Columns:

* ``bearing_rad``, ``bearing_deg`` -- wave-vector direction, clockwise from
  north. The crests of that wave run perpendicular to it.
* ``frequency_cycles_per_m``, ``wavelength_m``
* ``normalised_psd`` -- observed / background PSD

Four ``# key,value`` lines above the header record the confidence, the
threshold, the background seed and the number of background iterations.

The legacy LSDTopoTools ``<stem>_significant_freq.csv`` (columns
``bearing_rad, bearing_deg, frequency, wavelength, normalised_power``, 99.99%
threshold, no metadata) is read too, and renamed to the lsdtt3 columns, so old
and new runs plot with the same code.

The spectrum of a real surface is point-symmetric, so every significant
frequency appears twice, at bearings theta and theta + 180. A file holding only
the header is a valid result -- nothing passed the threshold -- and loads as an
empty table rather than raising.

Nothing here recomputes significance; this module only reads what the tool wrote.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

#: lsdtt3 column names, in file order.
COLUMNS = ["bearing_rad", "bearing_deg", "frequency_cycles_per_m",
           "wavelength_m", "normalised_psd"]

#: Legacy ``_significant_freq.csv`` column -> lsdtt3 column.
_LEGACY_RENAME = {"frequency": "frequency_cycles_per_m",
                  "wavelength": "wavelength_m",
                  "normalised_power": "normalised_psd"}

#: Threshold hard-coded in the legacy writer: 0.5 * 18.42185 (chi^2, 99.99%).
LEGACY_THRESHOLD = 9.210925

_SUFFIX = "_significant_frequencies_2d.csv"
_LEGACY_SUFFIX = "_significant_freq.csv"


@dataclass
class SignificantFrequencies2D:
    """One run's significant 2D spectrum cells.

    Attributes:
        table: one row per significant cell, columns :data:`COLUMNS`. Empty
            (with those columns) when nothing passed the threshold.
        metadata: the ``# key,value`` lines, values as floats where they parse.
            For legacy files: ``{"normalised_psd_threshold": 9.210925,
            "confidence": 0.9999}``.
        stem: output stem the file came from, for labelling.
        legacy: True when read from a legacy ``_significant_freq.csv``.
    """

    table: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=COLUMNS, dtype=float))
    metadata: dict = field(default_factory=dict)
    stem: str = ""
    legacy: bool = False

    @property
    def threshold(self) -> float:
        return float(self.metadata.get("normalised_psd_threshold", float("nan")))

    @property
    def n_cells(self) -> int:
        return int(len(self.table))


def _resolve(prefix: str | Path) -> Path:
    """Map a stem, a CSV path or a directory to the CSV to read."""
    p = Path(prefix)
    if p.is_file():
        return p
    if p.is_dir():
        found = sorted(p.glob("*" + _SUFFIX)) or sorted(p.glob("*" + _LEGACY_SUFFIX))
        if len(found) != 1:
            raise FileNotFoundError(
                f"expected exactly one *{_SUFFIX} (or legacy *{_LEGACY_SUFFIX}) "
                f"in {p}, found {len(found)}")
        return found[0]
    for suffix in (_SUFFIX, _LEGACY_SUFFIX):
        candidate = Path(f"{p}{suffix}")
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"no {_SUFFIX} or {_LEGACY_SUFFIX} for stem {p}")


def _parse_value(text: str):
    try:
        return float(text)
    except ValueError:
        return text


def load_significant_frequencies_2d(prefix: str | Path) -> SignificantFrequencies2D:
    """Load one run's significant 2D frequencies (lsdtt3 or legacy format).

    Args:
        prefix: ``<write_path>/<write_fname>`` stem, the CSV itself, or a
            directory holding exactly one such CSV.

    Returns:
        A :class:`SignificantFrequencies2D`.

    Raises:
        FileNotFoundError: if no matching file is found.
        ValueError: if the file lacks the expected columns.
    """
    path = _resolve(prefix)

    metadata = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            key, _, value = line[1:].strip().partition(",")
            if key:
                metadata[key.strip()] = _parse_value(value.strip())

    table = pd.read_csv(path, comment="#", skip_blank_lines=True)
    legacy = "normalised_power" in table.columns
    if legacy:
        table = table.rename(columns=_LEGACY_RENAME)
        metadata.setdefault("normalised_psd_threshold", LEGACY_THRESHOLD)
        metadata.setdefault("confidence", 0.9999)

    missing = [c for c in COLUMNS if c not in table.columns]
    if missing:
        raise ValueError(f"{path} is missing columns {missing}")
    table = table[COLUMNS].astype(float).reset_index(drop=True)

    name = path.name
    for suffix in (_SUFFIX, _LEGACY_SUFFIX):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return SignificantFrequencies2D(table=table, metadata=metadata,
                                    stem=name, legacy=legacy)


def load_significant_frequencies_2d_sweep(paths) -> pd.DataFrame:
    """Stack several runs (e.g. image tiles) into one long table.

    Args:
        paths: iterable of stems, CSV paths or directories.

    Returns:
        The rows of every run with a leading ``stem`` column. Runs with no
        significant cells contribute no rows; use
        :func:`summarise_significant_frequencies_2d` to keep one row per run.
    """
    frames = []
    for p in paths:
        sig = load_significant_frequencies_2d(p)
        if sig.n_cells:
            frames.append(sig.table.assign(stem=sig.stem))
    if not frames:
        return pd.DataFrame(columns=["stem"] + COLUMNS)
    out = pd.concat(frames, ignore_index=True)
    return out[["stem"] + COLUMNS]


def summarise_significant_frequencies_2d(paths) -> pd.DataFrame:
    """One row per run: how many cells passed, and the strongest one.

    Args:
        paths: iterable of stems, CSV paths or directories.

    Returns:
        Columns ``stem, n_cells, threshold, peak_wavelength_m,
        peak_bearing_deg, peak_normalised_psd``. Peak columns are NaN for a run
        with no significant cells. ``n_cells`` counts both members of each
        symmetric pair, as written.
    """
    rows = []
    for p in paths:
        sig = load_significant_frequencies_2d(p)
        row = {"stem": sig.stem, "n_cells": sig.n_cells,
               "threshold": sig.threshold,
               "peak_wavelength_m": float("nan"),
               "peak_bearing_deg": float("nan"),
               "peak_normalised_psd": float("nan")}
        if sig.n_cells:
            peak = sig.table.loc[sig.table["normalised_psd"].idxmax()]
            row.update(peak_wavelength_m=float(peak["wavelength_m"]),
                       peak_bearing_deg=float(peak["bearing_deg"]),
                       peak_normalised_psd=float(peak["normalised_psd"]))
        rows.append(row)
    return pd.DataFrame(rows)
