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
from .cosmo_transient import (
    cosmo_transient_title, load_cosmo_transient_table,
    plot_cosmo_transient_concentration, plot_cosmo_transient_erate_ratio,
    plot_nuclide_removal_response,
)
from .cosmo_transient_spatial import (
    load_joined_samples, plot_ratio_overlay_1d, plot_ratio_vs_mean_erate,
    plot_scenario_delta,
    plot_ratio_vs_steptime,
)
from .distributions import plot_distributions
from .divide_migration import (
    aggressor_annotation, plot_head_metric_vs_distance, plot_pair_comparison,
)
from .hand import (
    plot_hand_distance, plot_hand_window_populations, plot_valley_break_metric,
    plot_valley_break_fraction, plot_valley_floor_map, plot_valley_reach_map,
)
from .hydraulic_breaks import (
    plot_hydraulic_breaks_downstream, hydraulic_break_clusters,
    plot_hydraulic_break_clusters, plot_cross_section_informed_valley_mask,
)
from .landslides import (
    load_landslide_inventory,
    plot_landslide_area_volume_larsen,
    plot_sample_effect_apparent_vs_true,
    plot_sample_effect_ratio,
)
from .network_geometry import (
    plot_delta_L_cdf, plot_edge_basin_cdf, plot_junction_angle_cdf,
    plot_metric_panels, plot_metric_vs_temperature,
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
    "load_landslide_inventory",
    "plot_landslide_area_volume_larsen",
    "plot_sample_effect_apparent_vs_true",
    "plot_sample_effect_ratio",
    "circular_wasserstein", "kuiper_two_sample", "plot_circular_comparison",
    "load_cosmo_table", "plot_cosmo_sweep_panels", "plot_cosmo_sample_residuals",
    "load_cosmo_transient_table", "cosmo_transient_title",
    "plot_cosmo_transient_concentration", "plot_cosmo_transient_erate_ratio",
    "plot_nuclide_removal_response",
    "load_joined_samples", "plot_ratio_vs_steptime", "plot_ratio_vs_mean_erate",
    "plot_ratio_overlay_1d", "plot_scenario_delta",
    "plot_pair_comparison", "plot_head_metric_vs_distance",
    "aggressor_annotation",
    "plot_junction_angle_cdf", "plot_delta_L_cdf", "plot_edge_basin_cdf",
    "plot_metric_panels", "plot_metric_vs_temperature",
    "plot_profile", "plot_slope_area", "plot_gradient_curvature",
    "plot_valley_spectrum", "plot_spacing_vs_ci", "plot_landscape_panel",
]

_DEFERRED_MSG = (
    "This lsdviztools3.charts x-y plot is deferred to a later release and is "
    "not implemented yet."
)


def plot_profile(*args, **kwargs):
    """Longitudinal channel profile. Deferred to a later release."""
    raise NotImplementedError(_DEFERRED_MSG)


def plot_slope_area(*args, **kwargs):
    """Slope-area "boomerang" -- see :mod:`lsdviztools3.charts.slope_area`."""
    from .slope_area import plot_slope_area as _impl
    return _impl(*args, **kwargs)

def plot_valley_spectrum(*args, **kwargs):
    """One landscape's spectrum with both spacing measures marked.

    See :mod:`lsdviztools3.charts.valley_spacing`.
    """
    from .valley_spacing import plot_valley_spectrum as _impl
    return _impl(*args, **kwargs)

def plot_spacing_vs_ci(*args, **kwargs):
    """Valley spacing against the Bonetti channelization index, for a sweep.

    See :mod:`lsdviztools3.charts.valley_spacing`.
    """
    from .valley_spacing import plot_spacing_vs_ci as _impl
    return _impl(*args, **kwargs)

def plot_landscape_panel(*args, **kwargs):
    """Hillshade a set of landscapes side by side, for a visual sniff test.

    See :mod:`lsdviztools3.charts.valley_spacing`.
    """
    from .valley_spacing import plot_landscape_panel as _impl
    return _impl(*args, **kwargs)

def plot_gradient_curvature(*args, **kwargs):
    """Gradient-curvature diagnostic after Roering (2008).

    See :mod:`lsdviztools3.charts.gradient_curvature`.
    """
    from .gradient_curvature import plot_gradient_curvature as _impl
    return _impl(*args, **kwargs)
