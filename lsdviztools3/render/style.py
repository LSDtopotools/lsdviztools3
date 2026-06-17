"""Map styling configuration.

``MapStyle`` collects the styling knobs shared across render functions:
projection, frame, colour map, line/colour widths, fonts, output resolution,
figure size, and transparency.

Publication defaults follow common journal conventions:
- page width <= 17 cm; single-column ~8.5 cm (see ``SINGLE_COLUMN`` /
  ``DOUBLE_COLUMN``); the default is a comfortable 12 cm.
- fonts no smaller than 8 pt: tick labels 9 pt, axis labels 11 pt.
- DPI 300 (low-res) / 600 (high-res).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar, Optional


@dataclass
class MapStyle:
    """Styling options for a PyGMT map.

    Attributes:
        projection: GMT projection string. ``None`` triggers auto-projection.
        frame: GMT frame override. ``None`` triggers the auto publication frame.
        cmap: Colour map name or CPT path.
        line_width / line_color: default line styling.
        point_size: default point/symbol size.
        scalebar: scalebar options or ``None``.
        figure_size: figure width (GMT size string). Use ``SINGLE_COLUMN`` /
            ``DOUBLE_COLUMN`` for journal column widths.
        font_annotation: tick-label font (>= 8 pt; default 9 pt).
        font_label: axis-label font (default 11 pt).
        font_title: title font, e.g. the CRS label (default 11 pt).
        dpi / high_res_dpi: output resolution for low-/high-res saves.
        transparent: write a transparent background (with nodata transparent).
    """

    SINGLE_COLUMN: ClassVar[str] = "8.5c"
    DOUBLE_COLUMN: ClassVar[str] = "17c"
    MAX_WIDTH: ClassVar[str] = "17c"

    projection: Optional[str] = None
    frame: Optional[Any] = None
    cmap: str = "geo"
    line_width: str = "1p"
    line_color: str = "black"
    point_size: str = "0.2c"
    scalebar: Optional[dict] = field(default=None)
    figure_size: str = "12c"
    font_annotation: str = "9p"
    font_label: str = "11p"
    font_title: str = "11p"
    dpi: int = 300
    high_res_dpi: int = 600
    transparent: bool = True
