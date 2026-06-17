"""Backend-agnostic IO/data core.

Uses rioxarray/xarray, geopandas, and pandas. This subpackage must **not**
import PyGMT (the one nuance: ``grids.make_hillshade`` may call
``pygmt.grdgradient`` for data prep — see its docstring).
"""
