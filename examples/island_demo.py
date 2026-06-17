"""End-to-end demo: download a small island DEM and render shaded relief to PNG.

No data is committed — DEMs are downloaded on demand from the keyless Microsoft
Planetary Computer (Copernicus GLO-30). Runs headless (e.g. WSL2).

    python examples/island_demo.py                    # Eigg, local (UTM) coords
    python examples/island_demo.py rum geographic     # Rum, EPSG:4326
    python examples/island_demo.py viti_levu local
    python examples/island_demo.py guadalcanal

Second argument chooses coordinates: "local" (native UTM, default) or
"geographic" (EPSG:4326).
"""

from __future__ import annotations

import sys
from pathlib import Path

from lsdviztools3.download.ot_scraper import download_dem
from lsdviztools3.io.raster import load_raster
from lsdviztools3.regions import DEFAULT_TEST_REGION, REGIONS
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.render.style import MapStyle


def main(argv=None) -> int:
    """Download a region's DEM and save a shaded-relief PNG."""
    argv = argv if argv is not None else sys.argv[1:]
    key = argv[0] if argv else DEFAULT_TEST_REGION
    coordinates = argv[1] if len(argv) > 1 else "local"
    if key not in REGIONS:
        raise SystemExit(
            f"Unknown region '{key}'. Choose from: {', '.join(REGIONS)}"
        )
    region = REGIONS[key]

    out_dir = Path("demo_output")
    out_dir.mkdir(exist_ok=True)
    dem_path = out_dir / f"{key}_glo30.tif"
    png_path = out_dir / f"{key}_hillshade_{coordinates}.png"

    print(f"Downloading {region.name} DEM (GLO-30) -> {dem_path}")
    download_dem(region.bbox, source="planetary_computer", out_file=dem_path)

    style = MapStyle(cmap="geo", figure_size=MapStyle.SINGLE_COLUMN, transparent=False)
    dem = load_raster(dem_path)
    dem = dem.where(dem > 0)
    fig = render_hillshade(dem, coordinates=coordinates, style=style)
    save_figure(fig, png_path, style=style)
    print(f"Wrote {png_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
