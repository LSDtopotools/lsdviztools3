"""One unit -> colour assignment, keyed by unit NAME, shared by maps and charts.

Lithology codes are local to an lsdtt3 run (the same unit can be code 3 in one
run and 30 in another), so colours are keyed by the unit's name. Both the
geologic map (:mod:`render.lithology_map`) and the k_sn charts
(:mod:`charts.lithology_ksn`) call :func:`lithology_palette` with the same
list of names, so a unit has the same colour in every figure of a run.

Rule: the names (minus ``unmapped``) are sorted case-insensitively and take
the qualitative palette in that order; ``unmapped`` is always light grey. A
user mapping (``{name: colour}``, or a ``name,colour`` CSV / JSON file read by
:func:`lsdviztools3.io.lithology.read_unit_colours`) overrides any entry.

Colours are returned as ``#rrggbb``, which both GMT and matplotlib accept.
No PyGMT import (charts use this too).
"""

from __future__ import annotations

from typing import Iterable, Mapping, Optional

#: Name of the "no unit" class written by lsdtt3 (code -1).
UNMAPPED_NAME = "unmapped"
#: Its colour. Not pure grey (GMT red-NaN bug with grey CPT rows).
UNMAPPED_HEX = "#c8c8ca"
#: Names that are pooled summaries, never units.
NON_UNIT_NAMES = {"all"}


def qualitative_hex(n: int) -> list:
    """``n`` distinct colours: tab20 dark half, tab20 light half, tab20b, tab20c.

    Pure greys are nudged off grey (GMT red-NaN bug); colours cycle after 60.
    """
    from matplotlib import colormaps

    rgb = []
    for name in ("tab20", "tab20b", "tab20c"):
        cmap = colormaps[name]
        rgb += [tuple(int(round(255 * c)) for c in cmap(i)[:3]) for i in range(cmap.N)]
    order = list(range(0, 20, 2)) + list(range(1, 20, 2)) + list(range(20, len(rgb)))
    rgb = [rgb[i] for i in order]
    out = []
    for i in range(n):
        r, g, b = rgb[i % len(rgb)]
        if r == g == b:
            b = b + 1 if b < 255 else b - 1
        out.append(f"#{r:02x}{g:02x}{b:02x}")
    return out


def _to_hex(colour: str) -> str:
    """Any matplotlib colour or GMT ``r/g/b`` string -> ``#rrggbb``."""
    c = str(colour).strip()
    if "/" in c:
        r, g, b = (int(float(v)) for v in c.split("/")[:3])
        return f"#{r:02x}{g:02x}{b:02x}"
    from matplotlib.colors import to_hex

    return to_hex(c)


def lithology_palette(names: Iterable[str],
                      overrides: Optional[Mapping[str, str]] = None) -> dict:
    """``{name: "#rrggbb"}`` for every unit name, by the shared rule.

    Args:
        names: Unit names of the run (e.g. the values of
            :func:`io.lithology.resolve_code_names`, or a summary's
            ``lith_name`` column). Duplicates and the pooled ``all`` row are
            ignored.
        overrides: Optional ``{name: colour}`` taking precedence.

    Returns:
        ``{name: hex}`` including ``unmapped`` if it was in ``names``.
    """
    names = [str(n) for n in dict.fromkeys(names) if str(n) not in NON_UNIT_NAMES]
    units = sorted((n for n in names if n != UNMAPPED_NAME), key=str.casefold)
    pal = dict(zip(units, qualitative_hex(len(units))))
    if UNMAPPED_NAME in names:
        pal[UNMAPPED_NAME] = UNMAPPED_HEX
    for name, colour in (overrides or {}).items():
        if str(name) in pal:
            pal[str(name)] = _to_hex(colour)
    return pal


def code_colours(code_names: Mapping[int, str],
                 overrides: Optional[Mapping[str, str]] = None) -> dict:
    """``{code: hex}`` for a run's code->name mapping, via :func:`lithology_palette`."""
    pal = lithology_palette(code_names.values(), overrides)
    return {int(c): pal[str(n)] for c, n in code_names.items() if str(n) in pal}
