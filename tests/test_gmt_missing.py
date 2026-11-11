"""A missing GMT library gives one clear error, and the charts still import."""

import builtins
import subprocess
import sys

import pytest

from lsdviztools3 import _gmt


def _fake_import(exc):
    real_import = builtins.__import__

    def fake(name, *args, **kwargs):
        if name == "pygmt":
            raise exc
        return real_import(name, *args, **kwargs)

    return fake


class GMTCLibNotFoundError(Exception):
    """Stand-in with the name of PyGMT's own exception."""


def test_missing_gmt_library_is_explained(monkeypatch):
    monkeypatch.delitem(sys.modules, "pygmt", raising=False)
    monkeypatch.setattr(builtins, "__import__",
                        _fake_import(GMTCLibNotFoundError("libgmt not found")))
    with pytest.raises(_gmt.GMTNotAvailableError, match="pip cannot install GMT"):
        _gmt.require_pygmt()


def test_missing_pygmt_is_an_import_error(monkeypatch):
    monkeypatch.delitem(sys.modules, "pygmt", raising=False)
    monkeypatch.setattr(builtins, "__import__",
                        _fake_import(ImportError("No module named 'pygmt'")))
    with pytest.raises(ImportError, match="PyGMT is not installed"):
        _gmt.require_pygmt()


def test_other_errors_are_not_hidden(monkeypatch):
    monkeypatch.delitem(sys.modules, "pygmt", raising=False)
    monkeypatch.setattr(builtins, "__import__", _fake_import(RuntimeError("boom")))
    with pytest.raises(RuntimeError, match="boom"):
        _gmt.require_pygmt()


def test_cli_prints_message_not_traceback():
    code = (
        "import lsdviztools3.cli\n"
        "from lsdviztools3._gmt import GMTNotAvailableError\n"
        "raise GMTNotAvailableError('no GMT here')\n"
    )
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert run.returncode == 1
    assert "error: no GMT here" in run.stderr
    assert "Traceback" not in run.stderr
