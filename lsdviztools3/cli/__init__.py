"""Command-line entry points.

Each CLI builds a ``pygmt.Figure`` via the render layer and calls
``.savefig(png)`` — suitable for headless WSL2 / server use.
"""

import os as _os
import sys as _sys

# The commands only write files, so matplotlib's non-interactive Agg backend is
# always right. Set it before anything imports matplotlib: a notebook kernel
# exports MPLBACKEND for its own inline plots (Colab:
# module://matplotlib_inline.backend_inline), a command run from it with "!"
# inherits it, and matplotlib then fails on import with "is not a valid value
# for backend" when that module is not in the command's environment. Every
# lsdviz3-* command imports this package first, so this covers them all.
_os.environ["MPLBACKEND"] = "Agg"

from .._gmt import GMTNotAvailableError as _GMTNotAvailableError  # noqa: E402

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
