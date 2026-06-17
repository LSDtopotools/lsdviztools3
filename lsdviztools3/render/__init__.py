"""Thin PyGMT render layer.

Each ``render_*`` function consumes core (xarray/geopandas/pandas) objects and
returns a ``pygmt.Figure``. Functions accept an optional ``fig=`` so layers can
be composed onto a single figure.
"""
