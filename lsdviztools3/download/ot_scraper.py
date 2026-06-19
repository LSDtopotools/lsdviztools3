"""Multi-source DEM download.

Mirrors lsdtt's multi-source design: keyless Microsoft Planetary Computer
(STAC) and the OpenTopography REST API (requires a key resolved via
``apikey.resolve_ot_api_key``).

Planetary Computer is the default because it needs no credentials. For a small
bounding box the intersecting Copernicus GLO-30 tile is read remotely
(GDAL ``/vsicurl``) and clipped to the bbox, so only a small GeoTIFF is written.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Sequence, Union

from .apikey import mask_key_in_url, resolve_ot_api_key

logger = logging.getLogger(__name__)

PathLike = Union[str, Path]
BBox = Sequence[float]

PC_CATALOG_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
OT_GLOBAL_URL = "https://portal.opentopography.org/API/globaldem"


def _validate_bbox(bbox: BBox) -> tuple[float, float, float, float]:
    """Validate and normalize ``bbox`` as (west, south, east, north)."""
    if len(bbox) != 4:
        raise ValueError(
            f"bbox must contain exactly 4 values [west, south, east, north], got {bbox!r}"
        )
    west, south, east, north = (float(v) for v in bbox)
    if west >= east:
        raise ValueError(f"Invalid bbox: west must be < east (got west={west}, east={east}).")
    if south >= north:
        raise ValueError(
            f"Invalid bbox: south must be < north (got south={south}, north={north})."
        )
    if not (-180.0 <= west <= 180.0 and -180.0 <= east <= 180.0):
        raise ValueError(
            f"Invalid bbox longitudes: expected [-180, 180], got west={west}, east={east}."
        )
    if not (-90.0 <= south <= 90.0 and -90.0 <= north <= 90.0):
        raise ValueError(
            f"Invalid bbox latitudes: expected [-90, 90], got south={south}, north={north}."
        )
    return west, south, east, north


def download_dem(
    bbox: BBox,
    *,
    source: str = "planetary_computer",
    dataset: Optional[str] = None,
    out_file: PathLike,
    directory: Optional[PathLike] = None,
    api_key_file: Optional[str] = None,
) -> Path:
    """Download a DEM for a bounding box.

    Args:
        bbox: ``[west, south, east, north]`` in WGS84 decimal degrees.
        source: ``"planetary_computer"`` (keyless, default) or
            ``"opentopography"``.
        dataset: Dataset/collection name. Defaults to ``"cop-dem-glo-30"`` for
            Planetary Computer and ``"COP30"`` for OpenTopography.
        out_file: Output GeoTIFF path. Relative when ``directory`` is set.
        directory: Target directory for the download (created if absent). A
            relative ``out_file`` is written inside it.
        api_key_file: Optional OT key file (ignored for Planetary Computer).
            The key itself is resolved securely via ``apikey``.

    Returns:
        Path to the downloaded file.
    """
    from ..paths import resolve_in_dir

    out_file = resolve_in_dir(directory, out_file)
    west, south, east, north = _validate_bbox(bbox)
    if source == "planetary_computer":
        return _download_planetary_computer(
            (west, south, east, north),
            collection=dataset or "cop-dem-glo-30",
            out_file=out_file,
        )
    if source == "opentopography":
        return _download_opentopography(
            (west, south, east, north),
            dataset=dataset or "COP30",
            out_file=out_file,
            api_key_file=api_key_file,
        )
    raise ValueError(
        f"Unknown DEM source: {source!r} "
        "(expected 'planetary_computer' or 'opentopography')."
    )


def _download_planetary_computer(
    bbox: BBox, *, collection: str, out_file: PathLike
) -> Path:
    """Read the GLO-30 tile(s) intersecting ``bbox`` and clip to the bbox."""
    import planetary_computer
    import rioxarray  # noqa: F401  (registers the .rio accessor)
    from pystac_client import Client
    from rioxarray.merge import merge_arrays

    west, south, east, north = bbox
    catalog = Client.open(PC_CATALOG_URL)
    search = catalog.search(
        collections=[collection], bbox=[west, south, east, north]
    )
    items = list(search.items())
    if not items:
        raise RuntimeError(
            f"No '{collection}' items found for bbox {list(bbox)}."
        )

    clipped = []
    for raw_item in items:
        item = planetary_computer.sign(raw_item)
        href = item.assets["data"].href
        logger.info("Reading PC tile %s", item.id)
        da = rioxarray.open_rasterio(href, masked=True)
        da = da.rio.clip_box(minx=west, miny=south, maxx=east, maxy=north)
        clipped.append(da)

    merged = clipped[0] if len(clipped) == 1 else merge_arrays(clipped)

    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.rio.to_raster(out_path, driver="GTiff", compress="deflate")
    return out_path


def _download_opentopography(
    bbox: BBox, *, dataset: str, out_file: PathLike, api_key_file: Optional[str]
) -> Path:
    """Download a DEM from the OpenTopography REST API (requires a key)."""
    import requests

    west, south, east, north = bbox
    api_key = resolve_ot_api_key(api_key_file)
    params = {
        "demtype": dataset,
        "west": west,
        "south": south,
        "east": east,
        "north": north,
        "outputFormat": "GTiff",
        "API_Key": api_key,
    }

    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(OT_GLOBAL_URL, params=params, stream=True, timeout=600) as resp:
        logger.info("OpenTopography request: %s", mask_key_in_url(resp.url, api_key))
        resp.raise_for_status()
        with out_path.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
    return out_path
