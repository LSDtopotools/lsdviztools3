"""Readers for the outputs of ``lsdtt-slope-stability``.

The CLI writes rasters and CSV tables named ``<write_fname>_<suffix>``. This
module finds them in a run directory and reads them; nothing is recomputed. The
stability classes, factors of safety, probabilities, depths and catalogue
statistics are all lsdtt3's.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

#: Output suffixes, by key. A key whose file is absent is simply not drawn.
OUTPUTS: Dict[str, str] = {
    # SINMAP
    "sinmap_class": "_sinmap_stability_class.tif",
    "sinmap_class_areas": "_sinmap_class_areas.csv",
    "sinmap_curves": "_sinmap_slope_area_curves.csv",
    "dinf_slope": "_dinf_slope.tif",
    "specific_catchment_area": "_specific_catchment_area.tif",
    # Rupture planes (scaLr)
    "rupture_depth": "_rupture_depth.tif",
    "rupture_safety_factor": "_rupture_safety_factor.tif",
    "rupture_failure_probability": "_rupture_failure_probability.tif",
    "landslides": "_landslide_catalogue.csv",
    "landslides_amalgamated": "_landslide_catalogue_amalgamated.csv",
    # Slip surfaces
    "slip_safety_factor": "_slip_safety_factor.tif",
    "slip_evaluated": "_slip_evaluated.tif",
    "slip_failure_probability": "_slip_failure_probability.tif",
    "slip_max_depth": "_slip_max_depth.tif",
    "slip_body_depth": "_slip_body_depth.tif",
    "slip_cover_safety_factor": "_slip_cover_safety_factor.tif",
    "slip_cover_failure_probability": "_slip_cover_failure_probability.tif",
    "slip_bodies": "_slip_body_catalogue.csv",
    "slip_bodies_amalgamated": "_slip_body_catalogue_amalgamated.csv",
    # Groundwater
    "water_table_depth": "_water_table_depth.tif",
}

#: Columns a landslide or slip-body catalogue must have to be plotted.
CATALOGUE_COLUMNS = ("area_m2", "volume_m3", "mean_depth_m", "max_depth_m", "length_m", "width_m",
                     "length_to_width")

#: SINMAP stability classes (manual Table 4), in order from stable to defended.
SINMAP_CLASSES = {
    1: "stable",
    2: "moderately stable",
    3: "quasi-stable",
    4: "lower threshold",
    5: "upper threshold",
    6: "defended",
}


def slope_stability_files(data_dir, prefix: str) -> Dict[str, Path]:
    """Path of every known output of the run ``prefix`` in ``data_dir``, present or not."""
    d = Path(data_dir)
    return {key: d / f"{prefix}{suffix}" for key, suffix in OUTPUTS.items()}


def present(files: Dict[str, Path]) -> Dict[str, Path]:
    """The subset of ``files`` that exist."""
    return {k: p for k, p in files.items() if p.exists()}


def read_band(path) -> np.ndarray:
    """Read the first band of a raster, with nodata (and lsdtt3's -9999) as NaN."""
    try:
        import rasterio
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Reading rasters needs rasterio.") from exc
    with rasterio.open(path) as src:
        band = src.read(1).astype(float)
        nodata = src.nodata
    if nodata is not None:
        band = np.where(band == nodata, np.nan, band)
    return np.where(band <= -9990.0, np.nan, band)


def load_catalogue(path) -> pd.DataFrame:
    """A landslide or slip-body catalogue CSV, checked for the plotted columns."""
    df = pd.read_csv(path)
    missing = [c for c in CATALOGUE_COLUMNS if c not in df.columns]
    if missing:
        raise KeyError(f"{Path(path).name}: missing columns {missing}")
    return df


def load_sinmap_curves(path) -> pd.DataFrame:
    """The SINMAP slope - specific catchment area class boundaries."""
    df = pd.read_csv(path)
    missing = {"region_id", "curve", "slope", "specific_catchment_area_m"} - set(df.columns)
    if missing:
        raise KeyError(f"{Path(path).name}: missing columns {sorted(missing)}")
    return df


def load_sinmap_class_areas(path) -> pd.DataFrame:
    """The area in each SINMAP stability class, by region."""
    df = pd.read_csv(path)
    missing = {"region_id", "stability_class", "area_m2", "fraction_of_region"} - set(df.columns)
    if missing:
        raise KeyError(f"{Path(path).name}: missing columns {sorted(missing)}")
    return df
