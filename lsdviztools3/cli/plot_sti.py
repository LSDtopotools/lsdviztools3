"""CLI: lsdviz3-sti — drape STI tifs over hillshade to PNG."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from . import _common
from ..render.channels import render_drape
from ..render.sti import render_sti_p95_drape


def _non_p95_partner(sti_path: Path) -> Optional[Path]:
    """Return the non-P95 STI partner path for a P95 raster, if applicable."""
    name = sti_path.name
    if name.endswith("_sti_p95.tif"):
        return sti_path.with_name(name.replace("_sti_p95.tif", "_sti.tif"))
    if name.endswith("_sti_proxy_p95.tif"):
        return sti_path.with_name(name.replace("_sti_proxy_p95.tif", "_sti_proxy.tif"))
    return None


def _auto_base_hillshade(sti_dir: Path) -> Path:
    candidates = sorted(sti_dir.glob("*_hillshade.tif"))
    if not candidates:
        raise FileNotFoundError(
            f"No hillshade found in {sti_dir}. Pass --base explicitly."
        )
    return candidates[0]


def _discover_sti_rasters(sti_dir: Path) -> list[Path]:
    excluded_suffixes = (
        "_sti_series_codes.tif",
        "_sti_landuse_codes.tif",
        "_sti_valid_mask_ksd.tif",
        "_sti_valid_mask_spr.tif",
    )
    rasters = [
        p for p in sorted(sti_dir.glob("*_sti*.tif"))
        if ".tmp_cog" not in p.name
        and not p.name.endswith(excluded_suffixes)
    ]
    return rasters


def _style_for_sti_raster(path: Path, primary_colorbar: bool):
    name = path.name
    if name.endswith("_sti_twi.tif"):
        return {"cmap": "YlOrRd", "colorbar": True, "label": "STI TWI"}
    if name.endswith("_sti_term_ksd.tif"):
        return {"cmap": "YlOrRd", "colorbar": True, "label": "STI term (ksd)"}
    if name.endswith("_sti_term_spr.tif"):
        return {"cmap": "YlOrRd", "colorbar": True, "label": "STI term (spr)"}
    if name.endswith("_sti_valid_mask_ksd.tif") or name.endswith("_sti_valid_mask_spr.tif"):
        return {"cmap": "gray", "colorbar": True, "label": "STI valid mask"}
    if name.endswith("_sti_series_codes.tif"):
        return {"cmap": "turbo", "colorbar": True, "label": "STI series codes"}
    if name.endswith("_sti_landuse_codes.tif"):
        return {"cmap": "turbo", "colorbar": True, "label": "STI land-use codes"}
    if name.endswith("_sti_proxy_p95.tif"):
        return {"cmap": "YlOrRd", "colorbar": primary_colorbar, "label": "STI proxy p95"}
    if name.endswith("_sti_p95.tif"):
        return {"cmap": "YlOrRd", "colorbar": primary_colorbar, "label": "STI p95"}
    if name.endswith("_sti_proxy.tif"):
        return {"cmap": "YlOrRd", "colorbar": True, "label": "STI proxy"}
    if name.endswith("_sti.tif"):
        return {"cmap": "YlOrRd", "colorbar": True, "label": "STI"}
    return {"cmap": "YlOrRd", "colorbar": primary_colorbar, "label": "STI"}


def _render_intermediate(
    *,
    base_raster: str,
    drape_raster: str,
    out_path: Path,
    coordinates: str,
    transparency: float,
    cmap: str,
    colorbar: bool,
    colorbar_label: Optional[str],
    style,
    high_res: bool,
    directory,
):
    fig = render_drape(
        base_raster,
        drape_raster,
        coordinates=coordinates,
        base_is_hillshade=True,
        log_color=False,
        transparency=transparency,
        cmap=cmap,
        colorbar=colorbar,
        colorbar_label=colorbar_label,
        style=style,
    )
    _common.save_figure(fig, str(out_path), style=style, high_res=high_res,
                        directory=directory)


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-sti``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-sti",
        description=(
            "Drape an STI threshold raster (for example *_sti_p95.tif or "
            "*_sti_proxy_p95.tif) over a hillshade and render to PNG."
        ),
    )
    parser.add_argument(
        "--base",
        dest="base_raster",
        required=False,
        help=(
            "Base hillshade raster path (for example *_hillshade.tif). "
            "Required for single-raster mode; optional with --sti-dir "
            "(auto-detected if omitted)."
        ),
    )
    parser.add_argument(
        "--sti",
        dest="sti_raster",
        required=False,
        help="STI P95 raster path (for example *_sti_p95.tif or *_sti_proxy_p95.tif).",
    )
    parser.add_argument("--out", required=False, help="Output PNG path (single-raster mode).")
    parser.add_argument(
        "--sti-dir",
        default=None,
        help=(
            "Directory containing lsdtt3 STI outputs. When set, the CLI finds "
            "all *_sti*.tif rasters it can and plots each one."
        ),
    )
    parser.add_argument(
        "--transparency",
        type=float,
        default=35.0,
        help="STI drape transparency, 0 (opaque) to 100 (invisible); default 35.",
    )
    parser.add_argument(
        "--colorbar",
        action="store_true",
        help="Add a colourbar (default is no colourbar).",
    )
    _common.add_style_args(parser)
    parser.set_defaults(cmap="YlOrRd")
    _common.add_target_dir_arg(parser)
    return parser


def main(argv=None) -> int:
    """Entry point: render STI P95 drape over hillshade and save PNG."""
    args = build_parser().parse_args(argv)
    style = _common.style_from_args(args)

    if args.sti_dir:
        sti_dir = Path(args.sti_dir)
        if not sti_dir.is_dir():
            raise FileNotFoundError(f"STI directory not found: {sti_dir}")

        base_raster = Path(args.base_raster) if args.base_raster else _auto_base_hillshade(sti_dir)
        rasters = _discover_sti_rasters(sti_dir)
        if not rasters:
            raise FileNotFoundError(f"No STI rasters found in: {sti_dir}")

        out_dir = Path(args.directory) if args.directory else (sti_dir / "pngs")
        out_dir.mkdir(parents=True, exist_ok=True)

        for raster in rasters:
            cfg = _style_for_sti_raster(raster, primary_colorbar=args.colorbar)
            out_png = out_dir / f"{raster.stem}_drape.png"
            _render_intermediate(
                base_raster=str(base_raster),
                drape_raster=str(raster),
                out_path=out_png,
                coordinates=args.coordinates,
                transparency=args.transparency,
                cmap=cfg["cmap"],
                colorbar=cfg["colorbar"],
                colorbar_label=cfg["label"],
                style=style,
                high_res=args.high_res,
                directory=None,
            )
            print(f"Wrote STI plot: {out_png}")
        return 0

    if not args.base_raster or not args.sti_raster or not args.out:
        raise ValueError(
            "Single-raster mode requires --base, --sti, and --out. "
            "Or use --sti-dir for directory mode."
        )

    fig = render_sti_p95_drape(
        args.base_raster,
        args.sti_raster,
        coordinates=args.coordinates,
        transparency=args.transparency,
        cmap=args.cmap,
        colorbar=args.colorbar,
        style=style,
    )
    _common.save_figure(fig, args.out, style=style, high_res=args.high_res,
                        directory=args.directory)

    # If the primary input is a P95 raster, also plot the matching full STI
    # raster when present.
    sti_path = Path(args.sti_raster)
    partner = _non_p95_partner(sti_path)
    if partner is not None and partner.is_file():
        out_main = Path(args.out)
        out_base = out_main.with_suffix("")
        out_partner = Path(f"{out_base}_full.png")
        fig_partner = render_drape(
            args.base_raster,
            str(partner),
            coordinates=args.coordinates,
            base_is_hillshade=True,
            log_color=False,
            transparency=args.transparency,
            cmap=args.cmap,
            colorbar=True,
            colorbar_label="STI",
            style=style,
        )
        _common.save_figure(fig_partner, str(out_partner), style=style,
                            high_res=args.high_res, directory=args.directory)
        print(f"Wrote non-P95 STI plot: {out_partner}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
