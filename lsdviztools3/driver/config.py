"""Configuration for locating the lsdtt3 command-line executables.

The lsdtt3 tools (``lsdtt-channel-extraction``, ``lsdtt-basin-extraction``, …)
are separate executables that are often *not* on ``PATH`` (a common install
location is ``~/.local/lsdtt3/bin``). This module resolves the directory that
holds them, mirroring the API-key resolution pattern used in
:mod:`lsdviztools3.download.apikey`.

Resolution order for the lsdtt3 ``bin`` directory (first found wins):

1. An explicit ``lsdtt_dir`` argument.
2. The environment variable ``LSDTT3_DIR``.
3. The ``[lsdtt] dir`` key in a config file (see :func:`config_locations`).
4. The system ``PATH`` (``shutil.which``) as a last resort.
"""

from __future__ import annotations

import os
import re
import shutil
import tomllib
from pathlib import Path
from typing import Optional, Union

ENV_LSDTT_DIR = "LSDTT3_DIR"
ENV_CONFIG = "LSDVIZTOOLS3_CONFIG"

PathLike = Union[str, Path]

#: lsdtt3 documentation and installation guide, quoted in messages and help.
LSDTT3_REPOSITORY_URL = "https://github.com/LSDtopotools/lsdtt3"
LSDTT3_DOCS_URL = "https://lsdtopotools.github.io/lsdtt3/"
LSDTT3_INSTALL_URL = "https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/install_lsdtt3.html"

#: lsdtt3 command-line tools (the lsdtt_cli() targets of lsdtt3 0.9). Used only
#: for friendlier errors/validation; a tool missing here still runs.
KNOWN_TOOLS = (
    "lsdtt-advection-diffusion",
    "lsdtt-bankfull-width",
    "lsdtt-basin-extraction",
    "lsdtt-basin-stats",
    "lsdtt-basin-surface-metrics",
    "lsdtt-channel-extraction",
    "lsdtt-channel-network-analysis",
    "lsdtt-chi-analysis",
    "lsdtt-cosmo-column-tool",
    "lsdtt-cosmo-landslides",
    "lsdtt-cosmo-prediction",
    "lsdtt-cosmo-tool",
    "lsdtt-cross-sections",
    "lsdtt-distributed-hydrology",
    "lsdtt-divide-migration",
    "lsdtt-drainage-analysis",
    "lsdtt-feature-finder",
    "lsdtt-fetch-erosion-rates",
    "lsdtt-fetch-lithology",
    "lsdtt-fetch-raster",
    "lsdtt-fetch-river-lines",
    "lsdtt-flat-surface-finder",
    "lsdtt-fluvial-surfaces",
    "lsdtt-graphflood",
    "lsdtt-hillslope-channel-coupling",
    "lsdtt-internal-drainage",
    "lsdtt-lem",
    "lsdtt-linear-stats",
    "lsdtt-lithology-erodibility",
    "lsdtt-marsh-topography",
    "lsdtt-mountain-front",
    "lsdtt-outlet-matching",
    "lsdtt-raster-calculator",
    "lsdtt-raster-preprocessing",
    "lsdtt-ridgeline-analysis",
    "lsdtt-sediment-routing",
    "lsdtt-slope-area",
    "lsdtt-spectral-analysis",
    "lsdtt-stream-burning",
    "lsdtt-surface-metrics",
    "lsdtt-swath-profile",
    "lsdtt-synthetic-networks",
    "lsdtt-synthetic-surfaces",
    "lsdtt-valley-centreline",
    "lsdtt-water-centreline",
    "lsdtt-water-line-repair",
)


def config_locations() -> list[Path]:
    """Return candidate config-file paths, highest precedence first."""
    locs: list[Path] = []
    env = os.environ.get(ENV_CONFIG)
    if env:
        locs.append(Path(env).expanduser())
    locs.append(Path.cwd() / "lsdviztools3.toml")
    locs.append(Path.home() / ".config" / "lsdviztools3" / "config.toml")
    return locs


def _read_config_lsdtt_dir() -> Optional[Path]:
    """Return the ``[lsdtt] dir`` value from the first readable config file."""
    for path in config_locations():
        if not path.is_file():
            continue
        try:
            with path.open("rb") as fh:
                data = tomllib.load(fh)
        except (OSError, tomllib.TOMLDecodeError):
            continue
        value = data.get("lsdtt", {}).get("dir")
        if value:
            return Path(value).expanduser()
    return None


def resolve_lsdtt_dir(lsdtt_dir: Optional[PathLike] = None) -> Optional[Path]:
    """Resolve the directory holding the lsdtt-* executables.

    Returns ``None`` when nothing is configured, in which case callers should
    fall back to a bare ``PATH`` lookup.
    """
    if lsdtt_dir:
        return Path(lsdtt_dir).expanduser()
    env = os.environ.get(ENV_LSDTT_DIR)
    if env:
        return Path(env).expanduser()
    return _read_config_lsdtt_dir()


def resolve_lsdtt_tool(tool: str, lsdtt_dir: Optional[PathLike] = None) -> str:
    """Resolve an lsdtt3 tool name to a runnable executable path.

    Args:
        tool: tool name, e.g. ``"lsdtt-channel-extraction"``.
        lsdtt_dir: optional explicit bin directory (highest precedence).

    Returns:
        An absolute path (if a directory is configured) or the name as found on
        ``PATH``.

    Raises:
        FileNotFoundError: if the tool cannot be located, with guidance on how
            to configure the lsdtt3 directory.
    """
    directory = resolve_lsdtt_dir(lsdtt_dir)
    if directory is not None:
        candidate = directory / tool
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
        raise FileNotFoundError(_not_found_message(tool, directory))
    found = shutil.which(tool)
    if found:
        return found
    raise FileNotFoundError(_not_found_message(tool, None))


def _not_found_message(tool: str, directory: Optional[Path]) -> str:
    where = f"in configured directory '{directory}'" if directory else "on PATH"
    locs = ", ".join(str(p) for p in config_locations())
    return (
        f"Could not find lsdtt3 tool '{tool}' {where}.\n"
        f"Point lsdviztools3 at your lsdtt3 'bin' directory by any of:\n"
        f"  - export {ENV_LSDTT_DIR}=/path/to/lsdtt3/bin\n"
        f"  - pass lsdtt_dir=... to run_lsdtt()\n"
        f"  - add to a config file:\n"
        f"        [lsdtt]\n        dir = \"/path/to/lsdtt3/bin\"\n"
        f"    (searched, in order: {locs})\n"
        f"  - or run: lsdviztools3.driver.config.write_user_config('/path/to/lsdtt3/bin')\n"
        f"lsdtt3 not installed yet? See {LSDTT3_INSTALL_URL}"
    )


def write_user_config(lsdtt_dir: PathLike) -> Path:
    """Persist the lsdtt3 bin directory to the user config file.

    Writes or updates ``[lsdtt] dir`` in ``~/.config/lsdviztools3/config.toml``
    while preserving unrelated settings in the file.

    Returns:
        The path written.
    """
    target = Path.home() / ".config" / "lsdviztools3" / "config.toml"
    target.parent.mkdir(parents=True, exist_ok=True)
    resolved = Path(lsdtt_dir).expanduser().resolve()
    resolved_text = str(resolved).replace("\\", "\\\\").replace('"', '\\"')

    if not target.exists():
        target.write_text(f'[lsdtt]\ndir = "{resolved_text}"\n')
        return target

    text = target.read_text()
    updated = _upsert_lsdtt_dir(text, resolved_text)
    target.write_text(updated)
    return target


def _upsert_lsdtt_dir(text: str, resolved_text: str) -> str:
    """Upsert ``[lsdtt] dir = ...`` in TOML text while preserving other content."""
    lines = text.splitlines(keepends=True)

    lsdtt_start = None
    for idx, line in enumerate(lines):
        if line.strip() == "[lsdtt]":
            lsdtt_start = idx
            break

    if lsdtt_start is None:
        body = text if text.endswith("\n") or not text else f"{text}\n"
        sep = "" if (not body or body.endswith("\n\n")) else "\n"
        return f'{body}{sep}[lsdtt]\ndir = "{resolved_text}"\n'

    lsdtt_end = len(lines)
    for idx in range(lsdtt_start + 1, len(lines)):
        stripped = lines[idx].strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            lsdtt_end = idx
            break

    dir_line = f'dir = "{resolved_text}"\n'
    dir_pattern = re.compile(r"^\s*dir\s*=")
    for idx in range(lsdtt_start + 1, lsdtt_end):
        if dir_pattern.match(lines[idx]):
            lines[idx] = dir_line
            return "".join(lines)

    lines.insert(lsdtt_end, dir_line)
    return "".join(lines)
