"""CLI: lsdviz3-grabdem — download a DEM (Planetary Computer / OpenTopography).

The API key is never accepted as a CLI argument; pass ``--ot-api-key-file`` or
set ``OT_API_KEY`` (see ``download.apikey``).
"""

from __future__ import annotations

import argparse

from . import _common
from ..download.ot_scraper import download_dem


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-grabdem``.

    Returns:
        The configured parser (bbox/coords, ``--source``, ``--dataset``,
        ``--out``, ``--ot-api-key-file``).
    """
    parser = argparse.ArgumentParser(
        prog="lsdviz3-grabdem",
        description="Download a DEM from Planetary Computer or OpenTopography.",
    )
    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        metavar=("WEST", "SOUTH", "EAST", "NORTH"),
        required=True,
        help="Bounding box in WGS84 decimal degrees.",
    )
    parser.add_argument(
        "--source",
        choices=("planetary_computer", "opentopography"),
        default="planetary_computer",
        help="DEM source (default: planetary_computer).",
    )
    parser.add_argument(
        "--dataset",
        default=None,
        help="Dataset/collection name (source-dependent default if omitted).",
    )
    parser.add_argument("--out", required=True, help="Output GeoTIFF path.")
    parser.add_argument(
        "--ot-api-key-file",
        default=None,
        help="OpenTopography key file (ignored for planetary_computer).",
    )
    _common.add_target_dir_arg(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: download a DEM to the requested output file.

    Args:
        argv: Optional argument list.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    out = download_dem(
        args.bbox,
        source=args.source,
        dataset=args.dataset,
        out_file=args.out,
        directory=args.directory,
        api_key_file=args.ot_api_key_file,
    )
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
