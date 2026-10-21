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

:func:`load_spectral_analysis` gathers everything one run wrote -- the radial
spectrum, the rollover summary, the 2D observed and background PSD rasters, the
significant 2D cells and the input raster -- for
:func:`~lsdviztools3.charts.spectral.plot_spectral_tile_summary`.

The 2D PSD rasters are georeferenced in cycles/m with the zero-frequency cell
centred on (0, 0), and are north-up as written: row 0 holds the most northward
frequency, so a raster's own bounds place every cell correctly with no flipping.

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

#: Threshold lsdtt-spectral-analysis applies to the binned radial spectrum
#: (_significant_frequencies.csv): 0.5 * chi^2(2 DOF, 95%).
BINNED_THRESHOLD = 2.9955

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


@dataclass
class SpectralAnalysisRun:
    """Everything one ``lsdtt-spectral-analysis`` run wrote, as far as found.

    Members whose file is absent are None, so a run with some outputs switched
    off still loads.

    Attributes:
        stem: the output stem.
        radial: ``_radial_spectrum.csv`` (``write_radial_spectrum: true``).
        rollover: the ``#`` lines of ``_significant_frequencies.csv``
            (``rollover_frequency_cycles_per_m``, ``rollover_wavelength_m``,
            ``rollover_beta``, ``sub_rollover_beta`` and, from late 2026,
            ``normalised_psd_threshold`` and ``n_significant_frequencies``),
            as floats.
        binned_significant: the rows of ``_significant_frequencies.csv``.
        psd_2d, background_psd_2d: the 2D rasters as ``xarray.DataArray``
            (``write_2d_spectrum_rasters: true``), in cycles/m, north-up.
        significant_2d: ``_significant_frequencies_2d.csv``.
        input_raster: the analysed raster, found from
            ``_ingested_params.param`` unless given explicitly.
    """

    stem: str = ""
    radial: pd.DataFrame | None = None
    rollover: dict = field(default_factory=dict)
    binned_significant: pd.DataFrame | None = None
    psd_2d: object | None = None
    background_psd_2d: object | None = None
    significant_2d: SignificantFrequencies2D | None = None
    input_raster: object | None = None

    @property
    def normalised_psd_2d(self):
        """psd_2d / background_psd_2d, or None if either raster is missing."""
        if self.psd_2d is None or self.background_psd_2d is None:
            return None
        return self.psd_2d / self.background_psd_2d


def _read_params(path: Path) -> dict:
    out = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        key, _, value = line.partition(":")
        out[key.strip()] = value.strip()
    return out


def _input_from_params(base: Path) -> Path | None:
    params = Path(f"{base}_ingested_params.param")
    if not params.is_file():
        return None
    p = _read_params(params)
    fname = p.get("read_fname", "")
    if not fname:
        return None
    candidate = Path(p.get("read_path", "")) / fname
    if candidate.is_file():
        return candidate
    # Relative read_path: resolve against the directory holding the outputs.
    candidate = base.parent / p.get("read_path", "") / fname
    return candidate if candidate.is_file() else None


def load_spectral_analysis(prefix: str | Path,
                           input_path: str | Path | None = None) -> SpectralAnalysisRun:
    """Load one ``lsdtt-spectral-analysis`` output set.

    Args:
        prefix: ``<write_path>/<write_fname>`` stem.
        input_path: the analysed raster. If None it is read from
            ``<stem>_ingested_params.param`` (``read_path`` / ``read_fname``);
            if that cannot be found ``input_raster`` stays None.

    Returns:
        A :class:`SpectralAnalysisRun`.

    Raises:
        FileNotFoundError: if none of the run's outputs exist.
    """
    from .raster import load_raster

    base = Path(prefix)
    run = SpectralAnalysisRun(stem=base.name)
    found = False

    radial = Path(f"{base}_radial_spectrum.csv")
    if radial.is_file():
        run.radial = pd.read_csv(radial)
        found = True

    binned = Path(f"{base}_significant_frequencies.csv")
    if binned.is_file():
        with open(binned, encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                key, _, value = line[1:].strip().partition(",")
                try:
                    run.rollover[key.strip()] = float(value)
                except ValueError:
                    run.rollover[key.strip()] = float("inf") if value.strip() == "inf" else value
        run.binned_significant = pd.read_csv(binned, comment="#")
        found = True

    for attr, suffix in (("psd_2d", "_psd_2d"), ("background_psd_2d", "_background_psd_2d")):
        for ext in (".tif", ".bil", ".png"):
            path = Path(f"{base}{suffix}{ext}")
            if path.is_file():
                setattr(run, attr, load_raster(path))
                found = True
                break

    for suffix in (_SUFFIX, _LEGACY_SUFFIX):
        if Path(f"{base}{suffix}").is_file():
            run.significant_2d = load_significant_frequencies_2d(base)
            found = True
            break

    if not found:
        raise FileNotFoundError(f"no lsdtt-spectral-analysis outputs for stem {base}")

    src = Path(input_path) if input_path is not None else _input_from_params(base)
    if src is not None and src.is_file():
        run.input_raster = load_raster(src)
    return run
