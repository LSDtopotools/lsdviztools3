"""End-to-end lsdtt3 -> lsdviztools3 workflow.

Downloads a small island DEM (keyless Planetary Computer), then drives lsdtt3 to
write a hillshade (sea masked), extract the channel network, and extract the
single largest basin as a polygon — then renders a composite: grey hillshade +
channel overlay + semi-transparent largest-basin polygon. Headless (e.g. WSL2).

Requires lsdtt3 on ``LSDTT3_DIR`` / config / PATH (see lsdviztools3.driver.config).

    export LSDTT3_DIR=~/.local/lsdtt3/bin
    python examples/lsdtt3_workflow.py            # Isle of Eigg
    python examples/lsdtt3_workflow.py rum
"""

from __future__ import annotations

import sys
from pathlib import Path

from lsdviztools3.download.ot_scraper import download_dem
from lsdviztools3.regions import DEFAULT_TEST_REGION, REGIONS
from lsdviztools3.workflows import (
    lsdtt3_hillshade_channels_basin,
    render_hillshade_channels_basin,
)


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    key = argv[0] if argv else DEFAULT_TEST_REGION
    if key not in REGIONS:
        raise SystemExit(f"Unknown region '{key}'. Choose from: {', '.join(REGIONS)}")

    region = REGIONS[key]
    workdir = Path("demo_output") / f"lsdtt3_{key}"
    workdir.mkdir(parents=True, exist_ok=True)

    dem = workdir / f"{key}_glo30.tif"
    if not dem.exists():
        print(f"Downloading {region.name} DEM (GLO-30) -> {dem}")
        download_dem(region.bbox, source="planetary_computer", out_file=dem)

    outputs = lsdtt3_hillshade_channels_basin(dem, workdir, stem=key)
    png = Path("demo_output") / f"{key}_lsdtt3_composite.png"
    render_hillshade_channels_basin(outputs, png)
    print(f"Wrote {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
