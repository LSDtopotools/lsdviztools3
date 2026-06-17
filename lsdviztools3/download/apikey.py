"""Secure API-key resolution and URL masking.

Mirrors the ``download_dem_tiles.py`` pattern. The OpenTopography key is
resolved in strict order (first found wins): environment variable, then key
file. The key is **never** accepted as a direct CLI argument, and any URL that
embeds it is masked before logging.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

ENV_OT_API_KEY = "OT_API_KEY"


def resolve_ot_api_key(key_file: Optional[str] = None) -> str:
    """Resolve the OpenTopography API key.

    Resolution order (first found wins):
      1. Environment variable ``OT_API_KEY``.
      2. ``key_file``: first non-empty, non-comment (``#``) line.

    The key is intentionally never read from a CLI argument.

    Args:
        key_file: Optional path to a key file.

    Returns:
        The resolved API key.

    Raises:
        ValueError / FileNotFoundError: if no key can be resolved.
    """
    env = os.environ.get(ENV_OT_API_KEY)
    if env and env.strip():
        return env.strip()

    if key_file:
        path = Path(key_file).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"OpenTopography key file not found: {path}")
        for line in path.read_text().splitlines():
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                return stripped
        raise ValueError(f"No API key found in '{path}'")

    raise ValueError(
        "No OpenTopography API key found. Set the OT_API_KEY environment "
        "variable or pass a key file via api_key_file. (Planetary Computer is "
        "keyless — use source='planetary_computer' to avoid needing a key.)"
    )


def mask_key_in_url(url: str, key: str) -> str:
    """Replace occurrences of ``key`` in ``url`` with a masked placeholder.

    Args:
        url: A URL that may embed the API key.
        key: The secret to mask.

    Returns:
        The URL with the key replaced by ``***``.
    """
    if not key:
        return url
    return url.replace(key, "***")
