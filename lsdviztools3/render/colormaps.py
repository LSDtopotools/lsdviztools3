"""Colour-map resolution for the render layer.

GMT's ``makecpt`` accepts GMT's own master CPTs (``geo``, ``roma``, ``batlow``,
``viridis``, ``seis``, ...) and ``.cpt`` files, but **not** arbitrary matplotlib
colour maps (``cividis``, ``Spectral``, ``RdBu``, ``twilight``, ...). This module
bridges that gap: :func:`resolve_cmap` turns a matplotlib colour-map name into a
GMT CPT file on the fly, passes GMT names/files/colour-lists straight through,
and expands a single named colour (e.g. ``"black"``) into a solid CPT — handy for
threshold masks.
"""

from __future__ import annotations

import os
import tempfile
from functools import lru_cache
from typing import Optional


def _matplotlib_to_cpt(
    name: str,
    n: int = 256,
    *,
    avoid_pure_gray_bug: bool = True,
) -> str:
    """Write a continuous GMT CPT (z in 0..1) sampled from a matplotlib cmap.

    The file is cached per name; ``makecpt(cmap=<file>, series=[vmin, vmax])``
    rescales it to the data range.
    """
    import numpy as np
    from matplotlib import colormaps

    cmap = colormaps[name]
    rgb = (np.asarray([cmap(i / (n - 1))[:3] for i in range(n)]) * 255).round().astype(int)

    # Break pure-grey rows (r == g == b) only when we need to guard GMT's
    # nan_transparent behaviour. A grey CPT (e.g. "gray", "Greys", "binary")
    # combined with a white NaN slot can make GMT paint NaN cells solid red.
    # A +-1/255 nudge on blue is visually tiny but not perfectly neutral, so
    # callers can now disable this when transparency is not in play.
    if avoid_pure_gray_bug:
        grey = (rgb[:, 0] == rgb[:, 1]) & (rgb[:, 1] == rgb[:, 2])
        rgb[grey, 2] = np.where(rgb[grey, 2] >= 255, 254, rgb[grey, 2] + 1)

    z = np.linspace(0.0, 1.0, n)
    lines = ["# COLOR_MODEL = RGB", f"# Generated from matplotlib colormap {name!r}"]
    for i in range(n - 1):
        r0, g0, b0 = rgb[i]
        r1, g1, b1 = rgb[i + 1]
        lines.append(f"{z[i]:.6f} {r0}/{g0}/{b0} {z[i + 1]:.6f} {r1}/{g1}/{b1}")
    # Background/foreground/NaN slots (B/F clamp to the end colours).
    r0, g0, b0 = rgb[0]
    r1, g1, b1 = rgb[-1]
    lines += [f"B {r0}/{g0}/{b0}", f"F {r1}/{g1}/{b1}", "N 255/255/255"]

    fd, path = tempfile.mkstemp(prefix=f"lsdviz3_{name}_", suffix=".cpt")
    with os.fdopen(fd, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def _solid_cpt(color: str) -> str:
    """Write a GMT CPT (z 0..1) that reads as one solid colour throughout.

    Used for threshold masks (e.g. paint everything above a percentile black).
    The table is a single segment ramping from the requested colour to a *barely*
    perturbed, faintly non-grey version of it. Two GMT quirks force this:

    - a perfectly flat (zero colour-range) table makes ``nan_transparent`` paint
      NaN cells solid red, or segfaults GMT outright; and
    - a pure-grey table (``r == g == b``) with a white NaN slot triggers the same
      red-NaN bug.

    The +2/+1/+3 perturbation (toward mid-grey so it never clamps) is visually
    indistinguishable from the requested colour but keeps the table non-flat and
    non-grey, so NaN stays transparent.
    """
    from matplotlib.colors import to_rgb

    r, g, b = (int(round(c * 255)) for c in to_rgb(color))

    def _nudge(c: int, d: int) -> int:
        return c - d if c > 127 else c + d

    c0 = f"{r}/{g}/{b}"
    c1 = f"{_nudge(r, 2)}/{_nudge(g, 1)}/{_nudge(b, 3)}"
    lines = [
        "# COLOR_MODEL = RGB",
        f"# Near-solid colour {color!r}",
        f"0.0 {c0} 1.0 {c1}",
        f"B {c0}",
        f"F {c1}",
        "N 255/255/255",
    ]

    fd, path = tempfile.mkstemp(prefix=f"lsdviz3_solid_{color}_", suffix=".cpt")
    with os.fdopen(fd, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


@lru_cache(maxsize=None)
def resolve_cmap(
    cmap: Optional[str],
    *,
    avoid_pure_gray_bug: bool = True,
) -> Optional[str]:
    """Resolve a colour-map spec to something ``makecpt(cmap=...)`` accepts.

    Resolution order:

    - ``None``, a comma-list (``"red,blue"``), or an existing file path: passed
      through unchanged.
    - A matplotlib colour-map name (e.g. ``"cividis"``, ``"Spectral"``): converted
      to a GMT CPT file (matplotlib must be installed — the ``charts`` extra).
    - A single named colour (e.g. ``"black"``, ``"red"``): expanded to a solid
      one-colour CPT — useful for threshold masks.
    - Anything else (``"geo"``, ``"roma"``, ``"seis"``, a GMT ``.cpt`` name):
      passed through for GMT to resolve.

    Args:
        cmap: Colour-map name / colour / path, or ``None``.
        avoid_pure_gray_bug: If ``True`` (default), generated matplotlib grey
            CPTs are nudged very slightly off pure-grey to avoid GMT's
            ``nan_transparent`` red-NaN artefact. Set ``False`` for perfectly
            neutral greys when transparency is disabled.

    Returns:
        A GMT-acceptable colour-map spec (name or CPT file path), or ``None``.
    """
    if not cmap or "," in cmap or os.path.exists(cmap):
        return cmap
    try:
        from matplotlib import colormaps
        from matplotlib.colors import get_named_colors_mapping
    except ImportError:
        return cmap  # matplotlib unavailable: let GMT try the name as-is.

    if cmap in colormaps:
        return _matplotlib_to_cpt(cmap, avoid_pure_gray_bug=avoid_pure_gray_bug)
    if cmap in get_named_colors_mapping():
        return _solid_cpt(cmap)
    return cmap
