"""Maps comparing the floodplain / terrace methods of one lsdtt3 run.

lsdtt3 ``lsdtt-fluvial-surfaces`` with ``fluvial_surface_algorithm_choice`` a
list (``ladder,clubb2017,...``) writes, under ``<write_fname>``:

  <stem>_method_summary.csv                      one row per method, list order
  <stem>_<method>_surface_class.tif              0 none, 1 floodplain, 2 terrace
  <stem>_method_agreement_floodplain.tif         bit k = method k put floodplain
  <stem>_method_agreement_terrace.tif            ... terrace
  <stem>_mapped_surface_class.tif                (mapped_surfaces_fname) 1 / 2

``render_method_classes`` draws each method's classes in its own panel, side
by side, over the same hillshade; ``render_method_agreement`` draws one class's
agreement raster with a legend naming the methods behind each colour. Mapped
polygons, when given, are outlined on every panel. Visualisation only: every
value is read from lsdtt3's outputs.
"""

from __future__ import annotations

import io as _io
from dataclasses import replace
from pathlib import Path
from typing import Optional, Sequence

from .channels import _prepare_overlay, render_drape
from .style import MapStyle

CLASS_COLOURS = {1: "#1f6fd1", 2: "#e8890c"}   # floodplain, terrace
# Up to 7 combinations (3 methods); Okabe-Ito, then extras for 4 methods.
AGREEMENT_COLOURS = ("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7",
                     "#56B4E9", "#F0E442", "#882255", "#44AA99", "#117733",
                     "#999999", "#332288", "#DDCC77", "#AA4499", "#000000")


def read_method_list(stem) -> list:
    """Methods of a run, in list order (bit order of the agreement rasters).

    ``<stem>_method_summary.csv`` for a several-method run; a single-method run
    has no summary, and its one method is read from ``<stem>_method_record.csv``.
    """
    import pandas as pd

    summary = Path(f"{stem}_method_summary.csv")
    if summary.exists():
        return list(pd.read_csv(summary)["method"])
    record = Path(f"{stem}_method_record.csv")
    if record.exists():
        return [str(pd.read_csv(record)["method"].iloc[0])]
    raise FileNotFoundError(f"No {summary} or {record}.")


def method_stem(stem, method: str, methods: Sequence[str]) -> str:
    """Output stem of one method: ``<stem>_<method>`` in a list run."""
    return f"{stem}_{method}" if len(methods) > 1 else str(stem)


def _outline(fig, mapped, coordinates, style, pen="0.6p,black,4_2"):
    gdf, _ = _prepare_overlay(mapped, coordinates, style)
    if not gdf.empty:
        fig.plot(data=gdf, pen=pen)


def _panel_title(fig, text):
    fig.text(text=text, position="TL", offset="0.15c/-0.15c", justify="TL",
             font="9p,Helvetica-Bold,black", fill="white@20")


def render_method_classes(
    stem,
    base_raster,
    *,
    methods: Optional[Sequence[str]] = None,
    mapped=None,
    coordinates: str = "local",
    panel_width: str = "10c",
    gap_cm: float = 2.2,
    transparency: float = 25.0,
    style: Optional[MapStyle] = None,
):
    """One panel per method: floodplain and terrace over the hillshade.

    Args:
        stem: lsdtt3 output stem (read_path + write_fname).
        base_raster: Hillshade (e.g. lsdtt3 ``_hillshade.tif``), path or
            ``DataArray``, on the run's grid.
        methods: Methods to draw; default all, in list order.
        mapped: Optional mapped-surface polygons (path or ``GeoDataFrame``),
            outlined black dashed on every panel.
        coordinates: ``"local"`` / ``"geographic"`` / explicit CRS.
        panel_width: GMT width of each panel.
        gap_cm: Horizontal gap between panels (room for the next panel's
            northing labels).
        transparency: Percent transparency of the class colours.
        style: Styling; ``figure_size`` is replaced by ``panel_width``.

    Returns:
        The ``pygmt.Figure``.
    """
    import pygmt

    all_methods = read_method_list(stem)
    methods = list(methods) if methods is not None else all_methods
    style = replace(style or MapStyle(), figure_size=panel_width)
    width = float(panel_width.rstrip("c"))
    fig = pygmt.Figure()
    for k, m in enumerate(methods):
        cls = Path(f"{method_stem(stem, m, all_methods)}_surface_class.tif")
        if not cls.exists():
            raise FileNotFoundError(f"{cls} (write_surface_class_raster: true).")
        if k > 0:
            fig.shift_origin(xshift=f"{width + gap_cm}c")
        render_drape(base_raster, cls, coordinates=coordinates,
                     transparency=transparency,
                     categorical_gt_breaks=[0.5, 1.5],
                     categorical_colors=[CLASS_COLOURS[1], CLASS_COLOURS[2]],
                     style=style, fig=fig)
        if mapped is not None:
            _outline(fig, mapped, coordinates, style)
        _panel_title(fig, m)
    rows = [f"S 0.15c s 0.28c {CLASS_COLOURS[1]} - 0.4c floodplain",
            f"S 0.15c s 0.28c {CLASS_COLOURS[2]} - 0.4c terrace"]
    if mapped is not None:
        rows.append("S 0.15c - 0.4c - 0.6p,black,4_2 0.4c mapped")
    with pygmt.config(FONT_ANNOT_PRIMARY="8p"):
        fig.legend(spec=_io.StringIO("\n".join(rows) + "\n"),
                   position="JMR+jML+o0.3c/0c", box="+gwhite+p0.5p,black")
    return fig


def agreement_combinations(methods: Sequence[str]) -> list:
    """(value, label) for every non-zero bit pattern of ``methods``."""
    out = []
    for v in range(1, 2 ** len(methods)):
        names = [m for k, m in enumerate(methods) if v & (1 << k)]
        out.append((v, " + ".join(names) if len(names) < len(methods) else "all"))
    return out


def render_method_agreement(
    stem,
    base_raster,
    surface_class: str = "floodplain",
    *,
    mapped=None,
    coordinates: str = "local",
    transparency: float = 20.0,
    style: Optional[MapStyle] = None,
):
    """Which methods put ``surface_class`` (floodplain / terrace) on each cell.

    Reads ``<stem>_method_agreement_<surface_class>.tif``; each bit pattern is
    one colour, named in the legend by the methods that set it (bit k =
    method k in list order, from the summary CSV). Up to 4 methods.
    """
    import pygmt

    if surface_class not in ("floodplain", "terrace"):
        raise ValueError("surface_class is floodplain or terrace.")
    methods = read_method_list(stem)
    if len(methods) > 4:
        raise ValueError(f"{len(methods)} methods: the legend holds at most 4.")
    raster = Path(f"{stem}_method_agreement_{surface_class}.tif")
    if not raster.exists():
        raise FileNotFoundError(f"{raster} (a several-method run writes it).")
    combos = agreement_combinations(methods)
    colours = list(AGREEMENT_COLOURS[:len(combos)])
    style = style or MapStyle()
    fig = render_drape(base_raster, raster, coordinates=coordinates,
                       transparency=transparency,
                       categorical_gt_breaks=[v - 0.5 for v, _ in combos],
                       categorical_colors=colours, style=style)
    if mapped is not None:
        _outline(fig, mapped, coordinates, style)
    _panel_title(fig, f"{surface_class}: agreement")
    rows = [f"S 0.15c s 0.28c {c} - 0.4c {lab}" for (_, lab), c in zip(combos, colours)]
    if mapped is not None:
        rows.append("S 0.15c - 0.4c - 0.6p,black,4_2 0.4c mapped")
    with pygmt.config(FONT_ANNOT_PRIMARY="8p"):
        fig.legend(spec=_io.StringIO("\n".join(rows) + "\n"),
                   position="JMR+jML+o0.3c/0c", box="+gwhite+p0.5p,black")
    return fig
