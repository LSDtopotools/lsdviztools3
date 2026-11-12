"""Shared argparse helpers for the CLIs.

Centralises common flags (style options, output path) and the build-and-save
pattern so each CLI module stays thin.
"""

from __future__ import annotations

import argparse

from ..render.base import save_figure as _save_figure
from ..render.style import MapStyle

def add_target_dir_arg(parser: argparse.ArgumentParser) -> None:
    """Add the shared ``--dir`` target-directory flag to a parser.

    When set, outputs (and lsdtt3 read/write) are sent into this directory, which
    is created if it does not exist. A relative ``--out`` is resolved inside it.
    """
    parser.add_argument(
        "--dir",
        dest="directory",
        default=None,
        help="Target directory for outputs (created if absent). A relative "
        "--out is written inside it.",
    )


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
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Clean map: no frame ticks/labels, no CRS title, no colourbar — "
        "just the map with a small solid-rectangle scale bar.",
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
        clean=getattr(args, "clean", False),
    )


def save_figure(
    fig,
    out_path: str,
    *,
    style: MapStyle,
    high_res: bool = False,
    directory=None,
) -> None:
    """Save a figure to PNG (headless).

    Args:
        fig: The ``pygmt.Figure``.
        out_path: Output PNG path.
        style: Figure styling/output configuration.
        high_res: Use high-res output DPI.
        directory: Optional target directory (created if absent). A relative
            ``out_path`` is written inside it.
    """
    from ..paths import resolve_in_dir

    out_path = resolve_in_dir(directory, out_path)
    _save_figure(fig, out_path, style=style, high_res=high_res)


def add_class_args(parser: argparse.ArgumentParser, *, flag: str = "--classes",
                   what: str = "the values") -> None:
    """Add the colour-class options shared by the map commands.

    ``flag`` is the main option (``--classes``, or ``--ksn-classes`` for the
    k_sn maps); the others are ``--n-classes``, ``--class-colors`` and
    ``--class-color-mode``. Read them back with :func:`class_kwargs`.
    """
    parser.add_argument(
        flag, default=None, dest="classes",
        help=(f"Colour {what} in classes, with a class legend instead of a colourbar: "
              "comma-separated breaks (e.g. '50,100,200,400': 4 breaks make 5 classes), or a "
              "scheme that finds the breaks in the data: equal_interval, quantile, "
              "natural_breaks (or jenks). The colour limits (vmin/vmax) set the ends of the "
              "colour gradient and the range the schemes use."))
    parser.add_argument(
        "--n-classes", type=int, default=5,
        help="Number of classes for a scheme (default 5). Ignored with breaks.")
    parser.add_argument(
        "--class-colors", default=None,
        help="Comma-separated colours, one per class, instead of sampling the colour map "
             "(names, hex such as '#1b9e77', or r/g/b).")
    parser.add_argument(
        "--class-color-mode", default="auto", choices=["auto", "gradient", "distinct"],
        help="gradient: class colours along the colour map at each class's values; distinct: "
             "a qualitative map's colours in order (tab10, Set1, Dark2, ...); auto (default): "
             "distinct for matplotlib's qualitative maps, else gradient.")


def class_kwargs(args: argparse.Namespace) -> dict:
    """The ``classes``, ``n_classes``, ``class_colors``, ``class_color_mode``
    keyword arguments from :func:`add_class_args` options.

    Raises:
        SystemExit: a ``--classes`` value that is neither a scheme nor numbers.
    """
    from ..colour_classes import parse_classes_arg

    try:
        classes = parse_classes_arg(args.classes)
    except ValueError as err:
        raise SystemExit(str(err)) from None
    colours = None
    if args.class_colors:
        colours = [c.strip() for c in args.class_colors.split(",") if c.strip()]
    return dict(classes=classes, n_classes=args.n_classes, class_colors=colours,
                class_color_mode=args.class_color_mode)

