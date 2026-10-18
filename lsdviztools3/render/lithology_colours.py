"""One unit -> colour assignment, keyed by unit NAME, shared by maps and charts.

Lithology codes are local to an lsdtt3 run (the same unit can be code 3 in one
run and 30 in another), so colours are keyed by the unit's name. The geologic
map (:mod:`render.lithology_map`), the k_sn charts (:mod:`charts.lithology_ksn`)
and the bearing / junction charts (:mod:`charts.channel_bearings`) all call
:func:`lithology_palette`, so a unit has the same colour in every figure of a
run.

Rule: real unit names are sorted case-insensitively and take the qualitative
palette in that order. Reserved categories have fixed neutral colours and come
after the units, in this order:

* ``mixed (<N% one unit)`` (code -3; matched by the prefix ``mixed``) — mid grey;
* ``arms disagree ...`` (code -4; prefix ``arms disagree``) — dark grey;
* ``unmapped`` (code -1) — light grey.

A user mapping (``{name: colour}``, or a ``name,colour`` CSV / JSON file read by
:func:`lsdviztools3.io.lithology.read_unit_colours`) overrides any entry.

Colours are returned as ``#rrggbb``, which both GMT and matplotlib accept.
No PyGMT import (charts use this too).
"""

from __future__ import annotations

from typing import Iterable, Mapping, Optional

from ..io.lithology import reserved_category

#: Name of the "no unit" class written by lsdtt3 (code -1).
UNMAPPED_NAME = "unmapped"
#: Reserved category colours. None is pure grey (GMT red-NaN bug with grey
#: CPT rows), and each is distinguishable from the others.
RESERVED_HEX = {
    "mixed": "#8c8c8e",          # mid grey
    "arms_disagree": "#5a5a5c",  # dark grey
    "unmapped": "#c8c8ca",       # light grey
}
UNMAPPED_HEX = RESERVED_HEX["unmapped"]
#: Display order of the reserved categories (after the real units).
RESERVED_ORDER = ("mixed", "arms_disagree", "unmapped")
#: Names that are pooled summaries, never units.
NON_UNIT_NAMES = {"all"}


def qualitative_hex(n: int) -> list:
    """``n`` distinct colours: tab20 dark half, tab20 light half, tab20b, tab20c.

    Pure greys are nudged off grey (GMT red-NaN bug), and the two greys of
    tab20 are skipped so units are never confused with the reserved greys.
    Colours cycle after the palette is exhausted.
    """
    from matplotlib import colormaps

    rgb = []
    for name in ("tab20", "tab20b", "tab20c"):
        cmap = colormaps[name]
        rgb += [tuple(int(round(255 * c)) for c in cmap(i)[:3]) for i in range(cmap.N)]
    order = list(range(0, 20, 2)) + list(range(1, 20, 2)) + list(range(20, len(rgb)))
    rgb = [rgb[i] for i in order]
    rgb = [c for c in rgb if not (max(c) - min(c) <= 3)]  # drop greys
    out = []
    for i in range(n):
        r, g, b = rgb[i % len(rgb)]
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


def order_unit_names(names: Iterable[str]) -> list:
    """Real units sorted case-insensitively, then mixed, arms-disagree, unmapped.

    ``all`` (the pooled row) is dropped.
    """
    names = [str(n) for n in dict.fromkeys(names) if str(n) not in NON_UNIT_NAMES]
    units = sorted((n for n in names if reserved_category(n) is None), key=str.casefold)
    reserved = sorted((n for n in names if reserved_category(n) is not None),
                      key=lambda n: (RESERVED_ORDER.index(reserved_category(n)), n))
    return units + reserved


def lithology_palette(names: Iterable[str],
                      overrides: Optional[Mapping[str, str]] = None) -> dict:
    """``{name: "#rrggbb"}`` for every unit name, by the shared rule.

    Args:
        names: Unit names of the run (e.g. the values of
            :func:`io.lithology.resolve_code_names`, or a summary's
            ``lith_name`` column, or a grouped bearing table's lithology
            ``group_value`` column). Duplicates and ``all`` are ignored.
        overrides: Optional ``{name: colour}`` taking precedence.

    Returns:
        ``{name: hex}`` in :func:`order_unit_names` order.
    """
    ordered = order_unit_names(names)
    units = [n for n in ordered if reserved_category(n) is None]
    pal = dict(zip(units, qualitative_hex(len(units))))
    for n in ordered:
        cat = reserved_category(n)
        if cat is not None:
            pal[n] = RESERVED_HEX[cat]
    for name, colour in (overrides or {}).items():
        if str(name) in pal:
            pal[str(name)] = _to_hex(colour)
    return {n: pal[n] for n in ordered}


def code_colours(code_names: Mapping[int, str],
                 overrides: Optional[Mapping[str, str]] = None) -> dict:
    """``{code: hex}`` for a run's code->name mapping, via :func:`lithology_palette`."""
    pal = lithology_palette(code_names.values(), overrides)
    return {int(c): pal[str(n)] for c, n in code_names.items() if str(n) in pal}
