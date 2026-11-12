"""Colour data in classes: breaks, a colour per class, and legend labels.

The GIS way of colouring a map by value in a few steps instead of a continuous
scale. Rasters (:func:`lsdviztools3.render.channels.render_drape`), lines
(:func:`~lsdviztools3.render.channels.render_channels`, so the k_sn maps of
``lsdviz3-chi`` and ``lsdviz3-knickpoints``) and points
(:func:`~lsdviztools3.render.points.render_points`) all take the same
``classes``, ``n_classes`` and ``class_colors`` arguments and pass them here.

``classes`` is either

- the **breaks** themselves: a list, tuple or NumPy array of increasing values.
  ``n`` breaks make ``n + 1`` classes: below the first break, between each pair,
  and above the last; e.g. ``[50, 100, 200]`` gives
  ``below 50``, ``50 to 100``, ``100 to 200``, ``200 and above``;
- or a **scheme** that finds ``n_classes`` classes in the data:

  - ``"equal_interval"``: equal steps from the minimum to the maximum;
  - ``"quantile"``: the same number of values in every class;
  - ``"natural_breaks"`` (or ``"jenks"``): the Fisher-Jenks optimal breaks,
    which make the classes as internally similar as possible (the least sum
    of squared deviations from the class means).

A class runs from its lower break up to, but not including, its upper break, so
a value equal to a break goes in the class above it. Values that are not finite
(nodata) get no class.

With ``log_scale=True`` the schemes work on log10 of the (positive) values: equal
intervals become equal ratios, the natural scale for k_sn. Breaks you give are
used as they are, in data units.

Colours come from a colour map in one of two ways (``color_mode``):

- ``"gradient"``: each class takes the colour the continuous scale from
  ``vmin`` to ``vmax`` gives at the centre of the class, so unevenly spaced
  breaks give unevenly spaced colours, and a classed map matches the continuous
  map with the same limits. The end classes are centred between the limit and
  their break. Without limits, the data's minimum and maximum are the limits.
  With ``log_scale`` the positions are measured in log10.
- ``"distinct"``: the colours of a qualitative map (``tab10``, ``Set1``,
  ``Dark2``, ``Paired``, ...) in order, one per class, for classes that should
  look different rather than ordered. ``"auto"`` (the default) chooses this for
  matplotlib's qualitative maps and ``"gradient"`` for everything else.

Or give ``class_colors``: one colour per class.

``vmin`` and ``vmax`` also bound the schemes: ``equal_interval`` runs from
``vmin`` to ``vmax``, and ``quantile`` and ``natural_breaks`` use only the
values between them. Values beyond the limits still fall in the end classes.

This module only prepares the colours; it does not draw. It needs GMT only to
sample a GMT colour-map name.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, Union

import numpy as np

#: The data-driven schemes, and their accepted names.
SCHEMES = {
    "equal_interval": "equal_interval",
    "equal": "equal_interval",
    "quantile": "quantile",
    "quantiles": "quantile",
    "natural_breaks": "natural_breaks",
    "natural": "natural_breaks",
    "jenks": "natural_breaks",
}

#: Natural breaks are computed on at most this many values (an even sample of
#: the sorted data, which always includes the minimum and maximum). The exact
#: algorithm costs n_values^2 per class; GIS packages sample the same way.
NATURAL_BREAKS_MAX_VALUES = 2000

_RGB_RE = re.compile(r"^\d{1,3}(\.\d*)?/\d{1,3}(\.\d*)?/\d{1,3}(\.\d*)?$")


@dataclass(frozen=True)
class ColourClasses:
    """Classes for colouring: their breaks, colours and legend labels.

    Attributes:
        breaks: The ``n_classes - 1`` breaks between classes, increasing, in
            data units.
        colors: One GMT ``r/g/b`` colour per class, lowest class first.
        labels: One legend label per class.
        scheme: ``"manual"`` for breaks given by the user, else the scheme name.
        limits: ``(vmin, vmax)`` the gradient colours were placed between, in
            data units, or ``None`` for distinct or user-given colours.
    """

    breaks: Tuple[float, ...]
    colors: Tuple[str, ...]
    labels: Tuple[str, ...]
    scheme: str = "manual"
    limits: Optional[Tuple[float, float]] = None

    @property
    def n_classes(self) -> int:
        """Number of classes (one more than the number of breaks)."""
        return len(self.breaks) + 1

    def classify(self, values) -> np.ndarray:
        """Class index (0 = lowest class) of each value; NaN where not finite.

        Args:
            values: Array-like of numbers, any shape.

        Returns:
            A float array of the same shape holding ``0 .. n_classes - 1``, or
            NaN for values that are not finite.
        """
        v = np.asarray(values, dtype=float)
        out = np.digitize(v, np.asarray(self.breaks, dtype=float), right=False).astype(float)
        out[~np.isfinite(v)] = np.nan
        return out

    def legend_entries(self):
        """``[(class_index, label, colour)]``, for ``add_categorical_legend``."""
        return [(i, lbl, col) for i, (lbl, col) in enumerate(zip(self.labels, self.colors))]


ClassesSpec = Union[str, Sequence[float], np.ndarray, ColourClasses]


def _check_limits(vmin, vmax, log_scale: bool):
    for name, v in (("vmin", vmin), ("vmax", vmax)):
        if v is not None and not np.isfinite(float(v)):
            raise ValueError(f"{name} must be a finite number, got {v}.")
        if v is not None and log_scale and not float(v) > 0:
            raise ValueError(f"{name} must be > 0 on a log scale, got {v}.")
    if vmin is not None and vmax is not None and not float(vmin) < float(vmax):
        raise ValueError(f"vmin ({vmin}) must be below vmax ({vmax}).")


def _finite(values, log_scale: bool, vmin=None, vmax=None) -> np.ndarray:
    """Finite values within the limits, as log10 when ``log_scale``."""
    v = np.asarray(values, dtype=float).ravel()
    v = v[np.isfinite(v)]
    if vmin is not None:
        v = v[v >= float(vmin)]
    if vmax is not None:
        v = v[v <= float(vmax)]
    if log_scale:
        v = v[v > 0]
        v = np.log10(v)
    return v


def natural_breaks(values, n_classes: int, *, max_values: int = NATURAL_BREAKS_MAX_VALUES) -> np.ndarray:
    """Fisher-Jenks natural breaks of ``values`` into ``n_classes`` classes.

    The exact optimum by dynamic programming on the sorted values: the
    partition into ``n_classes`` runs of consecutive values with the least
    total sum of squared deviations from the class means. With more than
    ``max_values`` values an even sample of the sorted data is used.

    Args:
        values: Finite numbers.
        n_classes: Number of classes, at least 2.
        max_values: Sample size above which the data are sampled.

    Returns:
        The ``n_classes - 1`` breaks: each is the smallest value of the class
        above it, so :meth:`ColourClasses.classify` reproduces the partition.
        Fewer when the data have fewer distinct values than classes.
    """
    x = np.sort(np.asarray(values, dtype=float).ravel())
    x = x[np.isfinite(x)]
    if x.size > max_values:
        x = x[np.unique(np.linspace(0, x.size - 1, max_values).round().astype(int))]
    n = x.size
    k = int(n_classes)
    if np.unique(x).size <= k:
        u = np.unique(x)
        return u[1:]
    # Cumulative sums give the sum of squared deviations of any run x[i..j].
    s1 = np.concatenate([[0.0], np.cumsum(x)])
    s2 = np.concatenate([[0.0], np.cumsum(x * x)])
    i = np.arange(n)[:, None]       # run start
    j = np.arange(n)[None, :]       # run end (inclusive)
    count = (j - i + 1).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        ssd = (s2[j + 1] - s2[i]) - (s1[j + 1] - s1[i]) ** 2 / count
    ssd[i > j] = np.inf
    # cost[c][j]: least cost of splitting x[0..j] into c + 1 classes.
    cost = ssd[0].copy()
    starts = []                     # starts[c-1][j]: start of the last class
    for _c in range(1, k):
        total = np.full((n, n), np.inf)
        total[1:, :] = cost[:-1, None] + ssd[1:, :]   # last class starts at i >= 1
        best = np.argmin(total, axis=0)
        cost = total[best, np.arange(n)]
        starts.append(best)
    # Walk back from the last value to find where each class starts.
    bounds = []
    end = n - 1
    for best in reversed(starts):
        s = int(best[end])
        bounds.append(s)
        end = s - 1
    return np.array([x[s] for s in reversed(bounds)])


def class_breaks(values, scheme: str, n_classes: int = 5, *, log_scale: bool = False,
                 vmin: Optional[float] = None, vmax: Optional[float] = None) -> np.ndarray:
    """Breaks of a data-driven scheme, in data units.

    Args:
        values: The data (any shape; non-finite values are ignored).
        scheme: ``"equal_interval"``, ``"quantile"`` or ``"natural_breaks"``
            (``"jenks"``); see the module docstring.
        n_classes: Number of classes, at least 2.
        log_scale: Compute the breaks on log10 of the positive values.
        vmin, vmax: Limits in data units. ``equal_interval`` runs between them
            (default: the data's minimum and maximum); the other schemes use
            only the values between them.

    Returns:
        Increasing breaks above the bottom of the range; fewer than
        ``n_classes - 1`` when the data have too few distinct values (repeated
        breaks are merged).

    Raises:
        ValueError: an unknown scheme, ``n_classes < 2``, bad limits, or no
            usable values.
    """
    key = SCHEMES.get(str(scheme).strip().lower())
    if key is None:
        raise ValueError(f"Unknown classification scheme {scheme!r}: use one of "
                         "equal_interval, quantile, natural_breaks (jenks), or give the breaks.")
    if int(n_classes) < 2:
        raise ValueError(f"n_classes must be at least 2, got {n_classes}.")
    _check_limits(vmin, vmax, log_scale)
    v = _finite(values, log_scale, vmin, vmax)
    if v.size == 0:
        raise ValueError("No finite values to classify"
                         + (" within vmin/vmax" if vmin is not None or vmax is not None else "")
                         + (" (log scale: no positive values)." if log_scale else "."))
    to = (lambda x: np.log10(float(x))) if log_scale else float
    lo = to(vmin) if vmin is not None else float(v.min())
    hi = to(vmax) if vmax is not None else float(v.max())
    if key == "equal_interval":
        b = np.linspace(lo, hi, int(n_classes) + 1)[1:-1]
    elif key == "quantile":
        # Data values, not interpolated ones, as a GIS quantile scheme gives.
        b = np.quantile(v, np.linspace(0, 1, int(n_classes) + 1)[1:-1], method="higher")
    else:
        b = natural_breaks(v, int(n_classes))
    b = np.unique(b)
    b = b[b > lo]                   # a break at the bottom of the range would make an empty class
    return 10.0 ** b if log_scale else b


def _plain(v: float, digits: int) -> str:
    """``v`` to ``digits`` significant digits, written out in full (``1010``, not
    ``1.01e+03``) unless it is very large or very small."""
    if v == 0 or not (1e-4 <= abs(v) < 1e7):
        return f"{v:.{digits}g}"
    rounded = float(f"{v:.{digits}g}")
    text = f"{rounded:.10f}".rstrip("0").rstrip(".")
    return "0" if text in ("", "-0") else text


def _format_values(values: Sequence[float]) -> list:
    """Shortest formats (3 to 8 significant digits) that keep the values distinct."""
    for digits in range(3, 9):
        txt = [_plain(float(v), digits) for v in values]
        if len(set(txt)) == len(txt):
            return txt
    return [repr(float(v)) for v in values]


def class_labels(breaks: Sequence[float], label_format=None) -> list:
    """Legend labels for the classes made by ``breaks``.

    ``below a``, ``a to b``, ..., ``z and above`` (plain ASCII, so GMT legends
    show them as written; "to" rather than a dash, which reads as a minus sign
    next to negative numbers).

    Args:
        breaks: Increasing breaks.
        label_format: Optional format for the numbers: a format spec such as
            ``".0f"``, or a function from number to text. Default: the fewest
            significant digits (3 or more) that keep the breaks distinct.
    """
    b = list(breaks)
    if not b:
        return ["all values"]
    if label_format is None:
        t = _format_values(b)
    elif callable(label_format):
        t = [str(label_format(v)) for v in b]
    else:
        t = [format(v, label_format) for v in b]
    labels = [f"below {t[0]}"]
    labels += [f"{lo} to {hi}" for lo, hi in zip(t[:-1], t[1:])]
    labels.append(f"{t[-1]} and above")
    return labels


def _rgb_text(rgb) -> str:
    """``r/g/b`` (0-255) from a 0-1 RGB triple, nudged off pure grey.

    A pure-grey colour in a CPT with a NaN slot trips GMT's red-NaN bug (see
    ``render.colormaps``), so ``r == g == b`` moves the blue channel by one.
    """
    r, g, b = (int(round(255 * float(c))) for c in rgb[:3])
    if r == g == b:
        b = b + 1 if b < 255 else b - 1
    return f"{r}/{g}/{b}"


#: matplotlib's qualitative colour maps: their colours are used in order, one
#: per class, with ``color_mode="auto"``.
QUALITATIVE_CMAPS = frozenset({
    "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3",
    "tab10", "tab20", "tab20b", "tab20c",
})


def _gmt_cpt_colours(cmap: str, positions) -> list:
    """Colours of a GMT colour map (by name) at positions in [0, 1]."""
    import os
    import tempfile

    from ._gmt import require_pygmt
    pygmt = require_pygmt()

    n_slices = 256
    fd, path = tempfile.mkstemp(suffix=".cpt")
    os.close(fd)
    try:
        pygmt.makecpt(cmap=cmap, series=[0, n_slices, 1], output=path, color_model="r")
        slices = []
        for line in open(path):
            parts = line.split()
            if not parts or parts[0][0] in "#BFN":
                continue
            slices.append(parts[1])
    finally:
        os.remove(path)
    if not slices:
        raise ValueError(f"Could not read colours from the GMT colour map {cmap!r}.")
    out = []
    for t in positions:
        c = slices[min(len(slices) - 1, int(round(float(t) * (len(slices) - 1))))]
        if _RGB_RE.match(c):
            r, g, b = (float(x) for x in c.split("/"))
            out.append(_rgb_text((r / 255, g / 255, b / 255)))
        else:
            out.append(c)
    return out


def _matplotlib_cmap(cmap: str):
    try:
        from matplotlib import colormaps
        return colormaps[cmap]
    except (KeyError, ImportError):
        return None


def gradient_colours(cmap: str, positions) -> list:
    """Colours of a colour map at positions along it (0 = start, 1 = end).

    Args:
        cmap: A matplotlib colour-map name (``"viridis"``, ``"magma_r"``, ...)
            or a GMT one (``"geo"``, ``"roma"``, ...).
        positions: Numbers in [0, 1] (values outside are clamped).

    Returns:
        GMT ``r/g/b`` strings.
    """
    t = np.clip(np.asarray(positions, dtype=float), 0.0, 1.0)
    cm = _matplotlib_cmap(cmap)
    if cm is None:
        return _gmt_cpt_colours(cmap, t)
    return [_rgb_text(cm(float(x))) for x in t]


def distinct_colours(cmap: str, n: int) -> list:
    """The first ``n`` colours of a qualitative (listed) matplotlib colour map.

    Raises:
        ValueError: the map is not a listed one, or has fewer than ``n`` colours.
    """
    cm = _matplotlib_cmap(cmap)
    listed = getattr(cm, "colors", None) if cm is not None else None
    # Continuous maps such as viridis are also stored as 256 listed colours.
    if listed is None or len(listed) > 32:
        raise ValueError(f"color_mode='distinct' needs a qualitative matplotlib colour map "
                         f"such as tab10, Set1 or Dark2; {cmap!r} is not one.")
    if len(listed) < n:
        raise ValueError(f"{cmap!r} has {len(listed)} colours but there are {n} classes: "
                         "use a map with more colours (tab20 has 20), fewer classes, or class_colors.")
    return [_rgb_text(cm(i)) for i in range(n)]


def class_positions(breaks: Sequence[float], vmin: float, vmax: float, *,
                    log_scale: bool = False) -> np.ndarray:
    """Where each class's colour sits along a gradient from ``vmin`` to ``vmax``.

    The centre of each class, as a fraction of the range: ``(centre - vmin) /
    (vmax - vmin)``, in log10 when ``log_scale``. The lowest class runs from
    ``vmin`` to the first break and the highest from the last break to
    ``vmax``. Clamped to [0, 1], so breaks beyond the limits share an end
    colour.

    Args:
        breaks: Increasing breaks, in data units.
        vmin, vmax: The ends of the gradient, in data units.
        log_scale: Measure positions in log10 (breaks and limits must be > 0).

    Returns:
        ``len(breaks) + 1`` positions, lowest class first.
    """
    _check_limits(vmin, vmax, log_scale)
    b = np.asarray(breaks, dtype=float)
    if log_scale:
        if np.any(b <= 0):
            raise ValueError("Breaks must be > 0 on a log scale.")
        b, lo, hi = np.log10(b), np.log10(float(vmin)), np.log10(float(vmax))
    else:
        lo, hi = float(vmin), float(vmax)
    edges = np.concatenate([[lo], b, [hi]])
    centres = 0.5 * (edges[:-1] + edges[1:])
    return np.clip((centres - lo) / (hi - lo), 0.0, 1.0)


def sample_colours(cmap: str, n: int) -> list:
    """``n`` colours evenly spaced along a colour map, lowest first."""
    return gradient_colours(cmap, [0.5] if n == 1 else np.linspace(0.0, 1.0, n))


def _colour_text(colour) -> str:
    """One colour given by the user, as GMT accepts it."""
    c = str(colour).strip()
    if _RGB_RE.match(c):
        r, g, b = (float(x) for x in c.split("/"))
        return _rgb_text((r / 255, g / 255, b / 255))
    try:
        from matplotlib.colors import to_rgb
        return _rgb_text(to_rgb(c))
    except (ValueError, ImportError):
        return c                     # leave GMT-only names (e.g. "darkseagreen4") to GMT


def make_colour_classes(
    values=None,
    classes: ClassesSpec = "natural_breaks",
    *,
    n_classes: int = 5,
    cmap: str = "viridis",
    class_colors: Optional[Sequence[str]] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    log_scale: bool = False,
    color_mode: str = "auto",
    label_format=None,
) -> ColourClasses:
    """Breaks, colours and labels for colouring ``values`` in classes.

    Args:
        values: The data to classify (any shape; non-finite values are
            ignored). Needed for the schemes, and for the gradient limits when
            ``vmin`` / ``vmax`` are not given.
        classes: The breaks (list, tuple or array of increasing numbers), or a
            scheme (``"equal_interval"``, ``"quantile"``, ``"natural_breaks"``
            / ``"jenks"``), or a :class:`ColourClasses` (returned unchanged).
        n_classes: Number of classes for a scheme (at least 2). Ignored when
            the breaks are given: they set the number.
        cmap: Colour map for the class colours (matplotlib or GMT name).
        class_colors: One colour per class instead of ``cmap`` (matplotlib
            colour names or hex, GMT names, or ``r/g/b``).
        vmin, vmax: Limits in data units: the ends of the colour gradient, and
            the range the schemes work in. Default: the data's minimum and
            maximum (for breaks given without data, the first and last break).
        log_scale: Schemes and gradient positions work in log10 of the values
            (k_sn, drainage area, ...).
        color_mode: ``"gradient"``, ``"distinct"`` or ``"auto"``; see the
            module docstring.
        label_format: Number format for the legend labels; see
            :func:`class_labels`.

    Returns:
        A :class:`ColourClasses`.

    Raises:
        ValueError: breaks that are not finite and strictly increasing, an
            unknown scheme or colour mode, too few classes, bad limits, no
            usable data for a scheme, a ``class_colors`` list of the wrong
            length, or too many classes for a distinct colour map.
    """
    if isinstance(classes, ColourClasses):
        return classes
    if color_mode not in ("auto", "gradient", "distinct"):
        raise ValueError(f"color_mode must be 'auto', 'gradient' or 'distinct', got {color_mode!r}.")
    _check_limits(vmin, vmax, log_scale)
    if isinstance(classes, str):
        if values is None:
            raise ValueError(f"The {classes!r} scheme needs the data values.")
        breaks = class_breaks(values, classes, n_classes, log_scale=log_scale, vmin=vmin, vmax=vmax)
        scheme = SCHEMES[classes.strip().lower()]
        if breaks.size == 0:
            raise ValueError(f"The data have a single value in range: {classes!r} cannot make classes.")
    else:
        breaks = np.asarray(list(classes) if not isinstance(classes, np.ndarray) else classes,
                            dtype=float).ravel()
        scheme = "manual"
        if breaks.size == 0:
            raise ValueError("Give at least one break.")
        if not np.all(np.isfinite(breaks)):
            raise ValueError(f"Breaks must be finite numbers, got {breaks.tolist()}.")
        if breaks.size > 1 and not np.all(np.diff(breaks) > 0):
            raise ValueError(f"Breaks must be strictly increasing, got {breaks.tolist()}.")
        if log_scale and np.any(breaks <= 0):
            raise ValueError(f"Breaks must be > 0 on a log scale, got {breaks.tolist()}.")
    n = breaks.size + 1
    limits = None
    if class_colors is not None:
        colours = [_colour_text(c) for c in class_colors]
        if len(colours) != n:
            raise ValueError(f"class_colors has {len(colours)} colours but there are {n} classes "
                             f"({breaks.size} breaks).")
    else:
        mode = color_mode
        if mode == "auto":
            mode = "distinct" if cmap in QUALITATIVE_CMAPS else "gradient"
        if mode == "distinct":
            colours = distinct_colours(cmap, n)
        else:
            lo, hi = vmin, vmax
            data = _finite(values, False) if values is not None else np.array([])
            if log_scale:
                data = data[data > 0]
            if lo is None:
                lo = float(data.min()) if data.size else float(breaks[0])
            if hi is None:
                hi = float(data.max()) if data.size else float(breaks[-1])
            if not lo < hi:
                # One break and no data, or data all at one value: give the
                # range some width so the two classes take the two ends.
                lo, hi = (lo / 10, hi * 10) if log_scale else (lo - 1.0, hi + 1.0)
            limits = (float(lo), float(hi))
            colours = gradient_colours(cmap, class_positions(breaks, lo, hi, log_scale=log_scale))
    return ColourClasses(breaks=tuple(float(b) for b in breaks), colors=tuple(colours),
                         labels=tuple(class_labels(breaks, label_format)), scheme=scheme,
                         limits=limits)


def parse_classes_arg(text: Optional[str]):
    """Read a ``--classes`` command-line value: a scheme name, or breaks.

    ``"natural_breaks"`` -> the scheme; ``"50,100,200"`` -> ``[50.0, 100.0, 200.0]``;
    ``None`` or ``""`` -> ``None``.

    Raises:
        ValueError: neither a known scheme nor a comma-separated list of numbers.
    """
    if text is None or not str(text).strip():
        return None
    t = str(text).strip()
    if t.lower() in SCHEMES:
        return t.lower()
    try:
        return [float(x) for x in t.split(",") if x.strip()]
    except ValueError:
        raise ValueError(f"{text!r} is neither a classification scheme (equal_interval, quantile, "
                         "natural_breaks / jenks) nor a comma-separated list of breaks.") from None
