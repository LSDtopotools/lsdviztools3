"""Shared pieces of the ridgeline CLIs: driver files and output paths.

The ridgeline plots were first written as driver-file programs. The
``lsdviz3-ridgeline-*`` commands take ordinary flags, and ``--driver`` reads
the same ``key: value`` files: each key fills the flag of the same name
(``write_path`` -> ``--write-path``). Flags given on the command line win over
the driver. Relative paths resolve from the working directory, as in lsdtt3.
"""

from __future__ import annotations

import argparse
import ast
import re
from pathlib import Path
from typing import Any, Iterable, Optional


def _unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        try:
            return str(ast.literal_eval(value))
        except (SyntaxError, ValueError):
            return value[1:-1]
    return value


def _remove_comment(value: str) -> str:
    """Remove an unquoted ``#`` comment, keeping hashes inside quotes."""
    quote: Optional[str] = None
    escaped = False
    for index, char in enumerate(value):
        if escaped:
            escaped = False
        elif char == "\\" and quote:
            escaped = True
        elif quote:
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
        elif char == "#":
            return value[:index].rstrip()
    return value.strip()


def read_driver_file(path) -> dict[str, Any]:
    """Read ``key: value`` lines from a ridgeline plotting driver file.

    Blank lines and ``#`` comments are skipped; ``true``/``false``, integers
    and floats are converted, quoted values stay strings, and ``\\n`` in
    ``y_axis_label`` becomes a newline.

    Raises:
        ValueError: on a line without ``:``, an empty key or a repeated key.
    """
    path = Path(path).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"Driver file does not exist: {path}")
    values: dict[str, Any] = {}
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            raise ValueError(f"Invalid driver entry at {path}:{number}; expected key: value")
        key, raw_value = (part.strip() for part in line.split(":", 1))
        if not key:
            raise ValueError(f"Missing key at {path}:{number}")
        if key in values:
            raise ValueError(f"Duplicate driver key {key!r} at {path}:{number}")
        raw_value = _remove_comment(raw_value)
        if key == "y_axis_label":
            raw_value = raw_value.replace(r"\n", "\n")
        value: Any = _unquote(raw_value)
        if value.lower() in ("true", "false"):
            value = value.lower() == "true"
        elif re.fullmatch(r"[+-]?\d+", value):
            value = int(value)
        elif re.fullmatch(r"[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[eE][+-]?\d+)?|[+-]?\d+[eE][+-]?\d+",
                          value):
            value = float(value)
        values[key] = value
    return values


def add_driver_arg(parser: argparse.ArgumentParser) -> None:
    """Add ``--driver``."""
    parser.add_argument("--driver", default=None,
                        help="A key: value driver file; each key sets the flag of the same "
                             "name (write_path -> --write-path). Command-line flags win.")


def parse_with_driver(parser: argparse.ArgumentParser, argv=None, *,
                      extra_key_patterns: Iterable[str] = ()) -> tuple:
    """Parse arguments, taking defaults from ``--driver`` when given.

    Driver values are converted with each flag's ``type`` (lists given as
    comma-separated text are kept as text for the CLI to split).

    Args:
        parser: Parser with :func:`add_driver_arg`.
        argv: Argument list (None: ``sys.argv``).
        extra_key_patterns: Regexes for driver keys that are not flags (e.g.
            ``site_\\d+_.*``); they are returned instead of rejected.

    Returns:
        ``(args, extra)``: the namespace and the matched extra driver keys.
    """
    first, _ = parser.parse_known_args(argv)
    extra: dict[str, Any] = {}
    if first.driver:
        driver = read_driver_file(first.driver)
        actions = {a.dest: a for a in parser._actions}
        defaults = {}
        for key, value in driver.items():
            if any(re.fullmatch(p, key) for p in extra_key_patterns):
                extra[key] = value
                continue
            action = actions.get(key)
            if action is None or key in ("help", "driver"):
                parser.error(f"unknown driver key {key!r} in {first.driver}")
            if isinstance(action, argparse._StoreTrueAction) or action.type is _boolean:
                if not isinstance(value, bool):
                    parser.error(f"driver key {key!r} must be true or false")
            elif action.type is not None and not isinstance(value, bool):
                try:
                    value = action.type(str(value) if action.type is not str else value)
                except (TypeError, ValueError) as exc:
                    parser.error(f"driver key {key!r}: {exc}")
            defaults[key] = value
        parser.set_defaults(**defaults)
    return parser.parse_args(argv), extra


def _boolean(text: str) -> bool:
    lowered = str(text).strip().lower()
    if lowered in ("true", "1"):
        return True
    if lowered in ("false", "0"):
        return False
    raise argparse.ArgumentTypeError("expected true or false")


#: ``type=`` for boolean flags written ``--flag true|false`` (driver style).
boolean = _boolean


def number_list(text) -> list:
    """Split ``"25,49,99"`` (or a single number) into floats."""
    if isinstance(text, (list, tuple)):
        return [float(v) for v in text]
    parts = [p.strip() for p in str(text).split(",") if p.strip()]
    try:
        return [float(p) for p in parts]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected comma-separated numbers: {text!r}") from exc


def add_profile_selection_args(parser: argparse.ArgumentParser) -> None:
    """Add ``--ridgeline-number``, ``--reverse-ridgeline`` and the distance range."""
    parser.add_argument("--ridgeline-number", dest="ridgeline_number", type=int, default=0,
                        help="Ridgeline to plot (ridgeline_number, or basin_junction_number "
                             "in single-basin output; default 0).")
    parser.add_argument("--reverse-ridgeline", dest="reverse_ridgeline", type=boolean,
                        default=False, help="Reverse the distance direction (true/false).")
    parser.add_argument("--use-distance-range", dest="use_distance_range", type=boolean,
                        default=False, help="Clip to --minimum-distance..--maximum-distance.")
    parser.add_argument("--minimum-distance", dest="minimum_distance", type=float, default=0.0,
                        help="Inclusive start of the displayed distance range (m).")
    parser.add_argument("--maximum-distance", dest="maximum_distance", type=float, default=0.0,
                        help="Inclusive end of the displayed distance range (m).")


def distance_range(args) -> Optional[tuple]:
    """The ``(min, max)`` range from the selection flags, or None."""
    if not args.use_distance_range:
        return None
    if args.minimum_distance > args.maximum_distance:
        raise ValueError("minimum_distance must not exceed maximum_distance")
    return (args.minimum_distance, args.maximum_distance)


def add_output_args(parser: argparse.ArgumentParser, default_name: Optional[str]) -> None:
    """Add ``--write-path`` and ``--write-name``."""
    parser.add_argument("--write-path", dest="write_path", default=None,
                        help="Output directory; png/, svg/ and pdf/ are created in it.")
    parser.add_argument("--write-name", dest="write_name", default=default_name,
                        help="File stem of the outputs"
                             + (f" (default: {default_name})." if default_name else "."))


def require(parser: argparse.ArgumentParser, args, *names: str) -> None:
    """Fail like argparse when a flag that may come from the driver is missing."""
    missing = [n for n in names if getattr(args, n, None) in (None, "")]
    if missing:
        parser.error("missing required value(s): "
                     + ", ".join(f"--{n.replace('_', '-')}" for n in missing))
