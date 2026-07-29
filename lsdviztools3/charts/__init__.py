"""Statistics (x-y) plotting.

Chi-profile plotting is implemented here (:func:`plot_chi_profile`,
:func:`plot_chi_stacked`) on a matplotlib engine (DESIGN decision #1), imported
lazily so the core maps install stays lean. The remaining x-y plots (channel
longitudinal profiles, slope-area) are still deferred and raise
``NotImplementedError``.
"""

from __future__ import annotations

from .chi import plot_chi_profile, plot_chi_stacked, plot_ksn_vs_erate
from .circular import (
    circular_wasserstein, kuiper_two_sample, plot_circular_comparison,
)
from .cross_sections import (
    plot_valley_cross_section_map, plot_valley_cross_sections,
    plot_valley_curved_cross_section_map,
    plot_valley_cross_section_features,
    plot_valley_cross_section_features_map,
    plot_valley_cross_section_hydraulics,
    cross_section_hydraulics, detect_hydraulic_flats,
    combine_width_steps, compute_flat_surface_edges,
)
from .cosmo_comparison import (
    load_cosmo_table, plot_cosmo_sweep_panels, plot_cosmo_sample_residuals,
)
from .distributions import plot_distributions
from .hand import (
    plot_hand_distance, plot_hand_window_populations, plot_valley_break_metric,
    plot_valley_break_fraction, plot_valley_floor_map, plot_valley_reach_map,
)
from .hydraulic_breaks import (
    plot_hydraulic_breaks_downstream, hydraulic_break_clusters,
    plot_hydraulic_break_clusters, plot_cross_section_informed_valley_mask,
)
from .rose import rose_diagram
from .swath import plot_swath

__all__ = [
    "plot_chi_profile", "plot_chi_stacked", "plot_ksn_vs_erate", "plot_swath",
    "plot_hand_distance", "plot_hand_window_populations",
    "plot_valley_break_metric", "plot_valley_break_fraction",
    "plot_valley_cross_sections",
    "plot_valley_cross_section_map",
    "plot_valley_curved_cross_section_map",
    "plot_valley_cross_section_features",
    "plot_valley_cross_section_features_map",
    "plot_valley_cross_section_hydraulics",
    "cross_section_hydraulics", "detect_hydraulic_flats",
    "combine_width_steps", "compute_flat_surface_edges",
    "plot_valley_reach_map",
    "plot_valley_floor_map",
    "plot_hydraulic_breaks_downstream", "hydraulic_break_clusters",
    "plot_hydraulic_break_clusters",
    "plot_cross_section_informed_valley_mask",
    "plot_distributions", "rose_diagram",
    "circular_wasserstein", "kuiper_two_sample", "plot_circular_comparison",
    "load_cosmo_table", "plot_cosmo_sweep_panels", "plot_cosmo_sample_residuals",
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
