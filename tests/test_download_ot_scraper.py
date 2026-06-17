"""Unit tests for lsdviztools3.download.ot_scraper input validation."""

import pytest

from lsdviztools3.download.ot_scraper import download_dem


def test_download_dem_invalid_bbox_order_raises(tmp_path):
    with pytest.raises(ValueError, match="west must be < east"):
        download_dem((1, 0, 0, 1), out_file=tmp_path / "x.tif")


def test_download_dem_invalid_bbox_range_raises(tmp_path):
    with pytest.raises(ValueError, match=r"expected \[-90, 90\]"):
        download_dem((-6, -95, -5, -94), out_file=tmp_path / "x.tif")


def test_download_dem_unknown_source_raises(tmp_path):
    with pytest.raises(ValueError, match="Unknown DEM source"):
        download_dem((-6.2, 56.8, -6.1, 56.9), source="unknown", out_file=tmp_path / "x.tif")
