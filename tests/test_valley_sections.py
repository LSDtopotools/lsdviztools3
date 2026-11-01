"""Valley cross-section panels, with surface-class strips."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin

from lsdviztools3.charts.valley_sections import plot_valley_sections
from lsdviztools3.cli import plot_valley_sections as cli

CRS = "EPSG:32632"


def _write(tmp_path):
    """Two east-west sections across a 100 x 100 m (1 m) grid; the class
    raster is floodplain in the middle and terrace to the east."""
    rows = []
    for st, y in enumerate((70.5, 30.5)):
        for off in np.arange(-40.0, 41.0, 2.0):
            rows.append({"source_id": 0, "station_seq": st, "x": 50.0 + off, "y": y,
                         "offset_m": off, "elevation": 10.0 + 0.01 * off ** 2,
                         "left_offset_m": -20.0, "right_offset_m": 20.0,
                         "valley_width_m": 40.0, "along_distance_m": 100.0 * st})
    pd.DataFrame(rows).to_csv(tmp_path / "s.csv", index=False)
    cls = np.zeros((100, 100), dtype="float32")
    cls[:, 35:65] = 1
    cls[:, 70:85] = 2
    with rasterio.open(tmp_path / "cls.tif", "w", driver="GTiff", height=100, width=100,
                       count=1, dtype="float32", crs=CRS,
                       transform=from_origin(0.0, 100.0, 1.0, 1.0)) as dst:
        dst.write(cls, 1)


def test_class_strips_write(tmp_path):
    _write(tmp_path)
    out = plot_valley_sections(tmp_path / "s.csv", tmp_path / "o.png",
                               class_rasters=[("A", tmp_path / "cls.tif"),
                                              ("B", tmp_path / "cls.tif")], dpi=60)
    assert out.stat().st_size > 0


def test_cli_classes_need_label_and_path(tmp_path):
    _write(tmp_path)
    assert cli.main(["--samples", str(tmp_path / "s.csv"), "--out",
                     str(tmp_path / "o.png"), "--classes",
                     f"Clubb={tmp_path / 'cls.tif'}", "--dpi", "60"]) == 0
    with pytest.raises(SystemExit, match="LABEL=PATH"):
        cli.main(["--samples", str(tmp_path / "s.csv"), "--out",
                  str(tmp_path / "o.png"), "--classes", "nolabel"])
