"""Command-line entry points.

Each CLI builds a ``pygmt.Figure`` via the render layer and calls
``.savefig(png)`` — suitable for headless WSL2 / server use.
"""

import sys as _sys

from .._gmt import GMTNotAvailableError as _GMTNotAvailableError

_previous_excepthook = _sys.excepthook


def _gmt_excepthook(exc_type, exc, tb):
    """Print a missing GMT as one clear message, not a traceback.

    Every lsdviz3-* command imports this package first, so this covers them
    all. Any other exception goes to the previous hook unchanged.
    """
    if issubclass(exc_type, _GMTNotAvailableError):
        print(f"error: {exc}", file=_sys.stderr)
        return
    _previous_excepthook(exc_type, exc, tb)


_sys.excepthook = _gmt_excepthook
