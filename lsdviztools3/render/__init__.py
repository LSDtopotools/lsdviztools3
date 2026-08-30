"""Thin PyGMT render layer.

Each ``render_*`` function consumes core (xarray/geopandas/pandas) objects and
returns a ``pygmt.Figure``. Functions accept an optional ``fig=`` so layers can
be composed onto a single figure.
"""
from .hillslope_comparison_animation import (
    discover_frame_pairs, render_hillslope_comparison_frame,
    render_hillslope_comparison_animation, animate_hillslope_comparison,
)

__all__ = [
    "discover_frame_pairs", "render_hillslope_comparison_frame",
    "render_hillslope_comparison_animation", "animate_hillslope_comparison",
]
