"""IO for ``lsdtt-spectral-analysis`` valley-spacing output.

Two files per run, both written by ``compute_valley_spacing: true``:

* ``<stem>_valley_spacing.csv``          -- the scalar measures, name/value/units
* ``<stem>_valley_spacing_spectrum.csv`` -- frequency, wavelength, power

Nothing here computes a spacing. The measures come from lsdtt3
(``lsdtt_valley_spacing_cpu``); this module only reads what the tool wrote, so a
chart cannot silently disagree with the number the CLI printed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

#: Bonetti et al. (2020, PNAS 117:1375) critical channelization index.
CI_CRITICAL_LOW = 32.0
CI_CRITICAL_HIGH = 37.0


@dataclass
class ValleySpacing:
    """One landscape's valley-spacing measurement.

    Attributes:
        measures: name -> value, as written by the CLI.
        spectrum: the averaged power spectrum (frequency, wavelength, power).
        stem: the output stem the files came from, for labelling.
    """

    measures: dict = field(default_factory=dict)
    spectrum: Optional[pd.DataFrame] = None
    stem: str = ""

    def _get(self, key, default=float("nan")):
        v = self.measures.get(key, default)
        try:
            return float(v)
        except (TypeError, ValueError):
            return default

    @property
    def hillslope_length(self):
        """Distance-transform hillslope length (p90) -- the primary measure.

        Robust to branching, unlike the row-transect measures: oblique
        tributaries do not inflate it and it does not mix first-order hillslopes
        with inter-basin divides.

        The p90 rather than the mean: cells sit roughly uniformly at 0..L from
        the valley, so the mean is L/2 and would understate the length twofold by
        construction. Falls back to 2 x mean for output written before the p90
        was recorded.
        """
        v = self._get("hillslope_length")
        if v == v and v > 0:          # not NaN
            return v
        return 2.0 * self._get("hillslope_length_mean")

    @property
    def hillslope_length_median(self):
        return self._get("hillslope_length_median")

    # The three measures, named as the papers name them.
    @property
    def convex_wavelength(self):
        """Robust measure: 2x the mean convex run. Survives a small domain."""
        return self._get("convex_wavelength")

    @property
    def spectral_wavelength(self):
        """Perron et al. (2008): power-weighted mean of the spectral peaks."""
        return self._get("spectral_wavelength_power_weighted")

    @property
    def spectral_peak_wavelength(self):
        """The single strongest peak -- shown only to contrast with the above."""
        return self._get("spectral_wavelength_peak")

    @property
    def alpha(self):
        return self._get("spectral_exponent_alpha")

    @property
    def alpha_predicted(self):
        return self._get("spectral_exponent_predicted")

    @property
    def channelization_index(self):
        return self._get("channelization_index")

    @property
    def valleys_in_domain(self):
        return self._get("valleys_in_domain")

    @property
    def under_sampled(self) -> bool:
        """True when too few valleys fit for the SPECTRAL measure to mean much."""
        return bool(self._get("spectrum_under_sampled", 0.0))


def load_valley_spacing(prefix: str | Path) -> ValleySpacing:
    """Load one valley-spacing output set.

    Args:
        prefix: ``<write_path>/<write_fname>`` stem, or a directory holding
            exactly one such set.

    Returns:
        The loaded measures and spectrum. Missing members are left empty rather
        than raising, so a partial run can still be inspected.
    """
    base = Path(prefix)
    if base.is_dir():
        found = sorted(base.glob("*_valley_spacing.csv"))
        if not found:
            return ValleySpacing()
        base = Path(str(found[0])[: -len("_valley_spacing.csv")])

    out = ValleySpacing(stem=base.name)

    measures_path = Path(f"{base}_valley_spacing.csv")
    if measures_path.exists():
        table = pd.read_csv(measures_path)
        if {"quantity", "value"}.issubset(table.columns):
            out.measures = dict(zip(table["quantity"], table["value"]))

    spectrum_path = Path(f"{base}_valley_spacing_spectrum.csv")
    if spectrum_path.exists():
        out.spectrum = pd.read_csv(spectrum_path)

    return out


def load_valley_spacing_sweep(paths) -> pd.DataFrame:
    """Load several runs into one frame, one row per landscape.

    Args:
        paths: iterable of stems or directories.

    Returns:
        A DataFrame with the scalar measures as columns, sorted by
        channelization index so it plots directly against CI.
    """
    rows = []
    for p in paths:
        vs = load_valley_spacing(p)
        if not vs.measures:
            continue
        rows.append({
            "stem": vs.stem,
            "hillslope_length": vs.hillslope_length,
            "hillslope_length_median": vs.hillslope_length_median,
            "channelization_index": vs.channelization_index,
            "convex_wavelength": vs.convex_wavelength,
            "spectral_wavelength": vs.spectral_wavelength,
            "spectral_peak_wavelength": vs.spectral_peak_wavelength,
            "alpha": vs.alpha,
            "alpha_predicted": vs.alpha_predicted,
            "valleys_in_domain": vs.valleys_in_domain,
            "under_sampled": vs.under_sampled,
        })
    frame = pd.DataFrame(rows)
    if len(frame):
        frame = frame.sort_values("channelization_index").reset_index(drop=True)
    return frame
