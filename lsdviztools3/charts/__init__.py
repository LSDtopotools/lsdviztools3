"""Statistics (x-y) plotting.

Chi-profile plotting is implemented here (:func:`plot_chi_profile`,
:func:`plot_chi_stacked`) on a matplotlib engine (DESIGN decision #1), imported
lazily so the core maps install stays lean. The remaining x-y plots (channel
longitudinal profiles, slope-area) are still deferred and raise
``NotImplementedError``.
"""

from __future__ import annotations

from .chi import plot_chi_profile, plot_chi_stacked
from .distributions import plot_distributions
from .swath import plot_swath

__all__ = [
    "plot_chi_profile", "plot_chi_stacked", "plot_swath",
    "plot_distributions",
    "plot_profile", "plot_slope_area",
]

_DEFERRED_MSG = (
    "This lsdviztools3.charts x-y plot is deferred to a later release and is "
    "not implemented yet."
)


def plot_profile(*args, **kwargs):
    """Longitudinal channel profile. Deferred to a later release."""
    raise NotImplementedError(_DEFERRED_MSG)


def plot_slope_area(*args, **kwargs):
    """Slope-area plot. Deferred to a later release."""
    raise NotImplementedError(_DEFERRED_MSG)
