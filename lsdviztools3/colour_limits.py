"""Colour-scale limits given in data units, for log or linear colour scales.

The chi plots colour k_sn on a log10 scale. A user sets the range in k_sn
units (``vmin=10, vmax=500``); this converts it to the colour scale the plot
uses and says which ends of the data fall outside it, so the colourbar can
show extend arrows. Used by the matplotlib chi charts and by the GMT channel
renderer, so both treat the limits the same way.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


def colour_limits(
    colour_values: np.ndarray,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    *,
    log_scale: bool = True,
) -> Tuple[float, float, str]:
    """Resolve the colour-scale range.

    Args:
        colour_values: The values as coloured (already log10-transformed when
            ``log_scale``); NaNs are ignored.
        vmin, vmax: Limits in DATA units (e.g. k_sn), or ``None`` for the
            data's own minimum / maximum. With ``log_scale`` they must be > 0.
        log_scale: Whether ``colour_values`` are log10 of the data.

    Returns:
        ``(lo, hi, extend)``: the range on the colour scale, and the
        matplotlib-style extend flag (``"neither"``, ``"min"``, ``"max"``,
        ``"both"``) for data beyond it.

    Raises:
        ValueError: a limit <= 0 on a log scale, or ``vmin >= vmax``.
    """
    for name, v in (("vmin", vmin), ("vmax", vmax)):
        if v is not None and log_scale and not v > 0:
            raise ValueError(f"{name} must be > 0 on a log colour scale (got {v}); "
                             "give it in data units, e.g. k_sn.")
    if vmin is not None and vmax is not None and not vmin < vmax:
        raise ValueError(f"vmin ({vmin}) must be below vmax ({vmax}).")
    c = np.asarray(colour_values, dtype=float)
    finite = c[np.isfinite(c)]
    dlo, dhi = (float(finite.min()), float(finite.max())) if finite.size else (0.0, 1.0)
    to_scale = (lambda v: float(np.log10(v))) if log_scale else float
    lo = to_scale(vmin) if vmin is not None else dlo
    hi = to_scale(vmax) if vmax is not None else dhi
    if lo == hi:
        hi = lo + 1.0
    below = finite.size > 0 and dlo < lo
    above = finite.size > 0 and dhi > hi
    extend = "both" if (below and above) else "min" if below else "max" if above else "neither"
    return lo, hi, extend


def gmt_extend(extend: str) -> str:
    """The GMT colourbar ``-D`` modifier for an extend flag ("" for none):
    ``+eb`` (low end, background), ``+ef`` (high end, foreground), ``+e``."""
    return {"min": "+eb", "max": "+ef", "both": "+e"}.get(extend, "")
