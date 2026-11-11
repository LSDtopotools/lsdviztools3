"""Import PyGMT with a clear message when GMT itself is missing.

PyGMT is a Python wrapper around the GMT C library. pip installs PyGMT (it is a
dependency of lsdviztools3) but cannot install GMT, so a plain ``pip install``
gives a PyGMT that fails on import with ``GMTCLibNotFoundError``. Every map
renderer imports PyGMT through :func:`require_pygmt`, which turns that, and a
missing PyGMT, into :class:`GMTNotAvailableError` saying how to install GMT.

The import happens inside each renderer, not when lsdviztools3 is imported, so
the matplotlib charts keep working without GMT.
"""

from __future__ import annotations

INSTALL_HELP = (
    "lsdviztools3 draws maps with PyGMT, which needs the GMT library. pip cannot "
    "install GMT: install it from conda-forge (for example "
    "`conda install -c conda-forge gmt pygmt`, or use the pixi environment in "
    "the lsdviztools3 repository), Homebrew (`brew install gmt`) or your system "
    "package manager (GMT 6.5 or later). The matplotlib charts work without it. "
    "See https://lsdtopotools.github.io/lsdviztools3/ for installation."
)


class GMTNotAvailableError(ImportError):
    """PyGMT, or the GMT library it wraps, cannot be loaded.

    A subclass of ImportError, so ``pytest.importorskip`` and code that already
    catches ImportError treat a missing GMT like a missing package.
    """


def require_pygmt():
    """Return the ``pygmt`` module, or raise :class:`GMTNotAvailableError`."""
    try:
        import pygmt
    except ImportError as err:  # PyGMT itself is not installed
        raise GMTNotAvailableError(f"PyGMT is not installed. {INSTALL_HELP}") from err
    except Exception as err:  # noqa: BLE001 - PyGMT raises its own types on load
        # GMTCLibNotFoundError / GMTCLibError: PyGMT is there but cannot load
        # libgmt (missing, or too old). Match on the name: pygmt.exceptions
        # cannot be imported when importing pygmt itself has just failed.
        if type(err).__name__.startswith("GMTCLib"):
            raise GMTNotAvailableError(
                f"PyGMT is installed but cannot load the GMT library ({err}). {INSTALL_HELP}"
            ) from err
        raise
    return pygmt
