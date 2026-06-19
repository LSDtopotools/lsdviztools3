"""Target-directory helpers.

A small, shared way to send all of a run's artifacts — downloaded DEMs, lsdtt3
outputs, and rendered figures — into a single user-chosen directory, creating it
(and any parents) if it does not exist.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

PathLike = Union[str, Path]


def ensure_dir(directory: PathLike) -> Path:
    """Return ``directory`` as a ``Path``, creating it (and parents) if absent.

    Args:
        directory: Directory path.

    Returns:
        The directory as a ``Path`` (now guaranteed to exist).
    """
    path = Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_in_dir(target_dir: Optional[PathLike], path: PathLike) -> Path:
    """Resolve an output ``path`` inside ``target_dir``, creating dirs as needed.

    - ``target_dir is None``: ``path`` is used as given (its parent is created).
    - ``target_dir`` set and ``path`` relative: returns ``target_dir / path``.
    - ``target_dir`` set and ``path`` absolute: the absolute path wins (it is not
      relocated), but ``target_dir`` is still created.

    Either way the returned path's parent directory exists on return.

    Args:
        target_dir: Base directory for outputs, or ``None``.
        path: Output path (filename, relative path, or absolute path).

    Returns:
        The resolved output ``Path``.
    """
    resolved = Path(path)
    if target_dir is not None:
        base = ensure_dir(target_dir)
        if not resolved.is_absolute():
            resolved = base / resolved
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved
