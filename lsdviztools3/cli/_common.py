"""Shared argparse helpers for the CLIs.

Centralises common flags (style options, output path) and the build-and-save
pattern so each CLI module stays thin.
"""

from __future__ import annotations

import argparse

from ..render.base import save_figure as _save_figure
from ..render.style import MapStyle

def add_style_args(parser: argparse.ArgumentParser) -> None:
    """Add shared styling flags (cmap, projection, figure size, ...) to a parser.

    Args:
        parser: The argparse parser to extend.
    """
    defaults = MapStyle()
    parser.add_argument(
        "--coordinates",
        default="local",
        help=(
            "Coordinate mode: local/projected (default), geographic, or an "
            "explicit CRS such as EPSG:3857."
        ),
    )
    parser.add_argument(
        "--cmap",
        default=defaults.cmap,
        help=f"Colour map for raster rendering (default: {defaults.cmap}).",
    )
    parser.add_argument(
        "--projection",
        default=None,
        help="Optional explicit GMT projection string (auto by default).",
    )
    parser.add_argument(
        "--figure-size",
        default=defaults.figure_size,
        help=f"Figure width (GMT size string, default: {defaults.figure_size}).",
    )
    parser.add_argument(
        "--line-width",
        default=defaults.line_width,
        help=f"Line width for vector overlays (default: {defaults.line_width}).",
    )
    parser.add_argument(
        "--line-color",
        default=defaults.line_color,
        help=f"Line colour for vector overlays (default: {defaults.line_color}).",
    )
    parser.add_argument(
        "--font-annotation",
        default=defaults.font_annotation,
        help=f"Tick-label font (default: {defaults.font_annotation}).",
    )
    parser.add_argument(
        "--font-label",
        default=defaults.font_label,
        help=f"Axis-label font (default: {defaults.font_label}).",
    )
    parser.add_argument(
        "--font-title",
        default=defaults.font_title,
        help=f"Title font (default: {defaults.font_title}).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=defaults.dpi,
        help=f"Output DPI for normal mode (default: {defaults.dpi}).",
    )
    parser.add_argument(
        "--high-res",
        action="store_true",
        help=f"Use high-res output DPI ({defaults.high_res_dpi}).",
    )
    parser.add_argument(
        "--opaque",
        action="store_true",
        help="Write an opaque background instead of transparent output.",
    )


def style_from_args(args: argparse.Namespace) -> MapStyle:
    """Build a ``MapStyle`` from parsed CLI arguments.

    Args:
        args: Parsed argparse namespace.

    Returns:
        A populated ``MapStyle``.
    """
    return MapStyle(
        projection=getattr(args, "projection", None),
        cmap=getattr(args, "cmap", MapStyle().cmap),
        line_width=getattr(args, "line_width", MapStyle().line_width),
        line_color=getattr(args, "line_color", MapStyle().line_color),
        figure_size=getattr(args, "figure_size", MapStyle().figure_size),
        font_annotation=getattr(args, "font_annotation", MapStyle().font_annotation),
        font_label=getattr(args, "font_label", MapStyle().font_label),
        font_title=getattr(args, "font_title", MapStyle().font_title),
        dpi=getattr(args, "dpi", MapStyle().dpi),
        high_res_dpi=MapStyle().high_res_dpi,
        transparent=not getattr(args, "opaque", False),
    )


def save_figure(fig, out_path: str, *, style: MapStyle, high_res: bool = False) -> None:
    """Save a figure to PNG (headless).

    Args:
        fig: The ``pygmt.Figure``.
        out_path: Output PNG path.
        style: Figure styling/output configuration.
        high_res: Use high-res output DPI.
    """
    _save_figure(fig, out_path, style=style, high_res=high_res)
