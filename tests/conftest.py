"""Shared pytest fixtures.

Provides a small DEM downloaded once from the keyless Microsoft Planetary
Computer (Copernicus GLO-30, Isle of Eigg) and cached under ``tests/_cache``.
No data is committed. Tests that need the network are marked ``network`` and
skip cleanly when the download is unavailable (e.g. offline).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lsdviztools3.regions import DEFAULT_TEST_REGION, REGIONS

CACHE_DIR = Path(__file__).parent / "_cache"


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "network: test requires network access (DEM download)."
    )
    config.addinivalue_line(
        "markers", "slow: slow test (large DEM / multi-step lsdtt3 analysis)."
    )


@pytest.fixture(scope="session")
def small_dem() -> Path:
    """Path to a small Eigg GLO-30 DEM, downloaded once and cached.

    Skips the test if the download fails (offline / STAC / GDAL issues).
    """
    from lsdviztools3.download.ot_scraper import download_dem

    region = REGIONS[DEFAULT_TEST_REGION]
    CACHE_DIR.mkdir(exist_ok=True)
    dem_path = CACHE_DIR / f"{DEFAULT_TEST_REGION}_glo30.tif"
    if not dem_path.exists():
        try:
            download_dem(
                region.bbox, source="planetary_computer", out_file=dem_path
            )
        except Exception as exc:  # noqa: BLE001 - network/STAC/GDAL -> skip
            pytest.skip(f"Could not download test DEM: {exc}")
    return dem_path


@pytest.fixture(scope="session")
def guadalcanal_dem() -> Path:
    """Path to a Guadalcanal (Solomon Is.) GLO-30 DEM, downloaded once and cached.

    An island with dense dendritic drainage and complete (sea-bounded) basins —
    the chi-analysis test landscape. Larger than Eigg; spans several GLO-30
    tiles that are merged on download. Skips if the download is unavailable.
    """
    from lsdviztools3.download.ot_scraper import download_dem

    region = REGIONS["guadalcanal"]
    CACHE_DIR.mkdir(exist_ok=True)
    dem_path = CACHE_DIR / "guadalcanal_glo30.tif"
    if not dem_path.exists():
        try:
            download_dem(
                region.bbox, source="planetary_computer", out_file=dem_path
            )
        except Exception as exc:  # noqa: BLE001 - network/STAC/GDAL -> skip
            pytest.skip(f"Could not download Guadalcanal DEM: {exc}")
    return dem_path
