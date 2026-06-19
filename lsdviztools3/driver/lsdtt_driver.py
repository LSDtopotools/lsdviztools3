"""Subprocess wrapper for running lsdtt3 tools.

The lsdtt3 executable is located via :mod:`lsdviztools3.driver.config`
(explicit ``lsdtt_dir`` > ``LSDTT3_DIR`` env var > config file > ``PATH``).
lsdtt3 tools are invoked as ``<exe> <path> <parameter_file>`` (the ``<path>``
defaults to the parameter file's own directory).
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Sequence, Union

from .config import resolve_lsdtt_tool

ParamSpec = Union[dict, str, Path]
PathLike = Union[str, Path]

# GDAL/PROJ environment variables a conda/pixi env sets that point at its own
# (incompatible) data and plugins. lsdtt3 is an external native binary with its
# own GDAL, so these must not leak into it.
_CONDA_GDAL_VARS = ("GDAL_DRIVER_PATH", "GDAL_DATA", "PROJ_LIB", "PROJ_DATA")


def _external_env(clean_gdal_env: bool) -> dict:
    """Environment for the lsdtt3 subprocess, optionally dropping conda GDAL vars."""
    env = os.environ.copy()
    if clean_gdal_env:
        for var in _CONDA_GDAL_VARS:
            env.pop(var, None)
    return env


def _write_parameter_file(parameters: dict, directory: Path) -> Path:
    """Serialise a parameter dict to a temporary lsdtt ``.driver`` file.

    lsdtt parameter files are simple ``key: value`` lines.
    """
    directory.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=".driver", dir=directory, delete=False
    )
    with handle as fh:
        for key, value in parameters.items():
            fh.write(f"{key}: {value}\n")
    return Path(handle.name)


def run_lsdtt(
    tool: str,
    parameters: ParamSpec,
    *,
    directory: Optional[PathLike] = None,
    data_path: Optional[PathLike] = None,
    lsdtt_dir: Optional[PathLike] = None,
    extra_args: Optional[Sequence[str]] = None,
    cwd: Optional[PathLike] = None,
    check: bool = True,
    clean_gdal_env: bool = True,
) -> "subprocess.CompletedProcess":
    """Run an lsdtt3 command-line tool with a parameter file or dict.

    Args:
        tool: lsdtt3 tool name, e.g. ``"lsdtt-channel-extraction"``.
        parameters: path to an existing parameter file, or a dict that is
            serialised to a temporary ``.driver`` file.
        directory: target directory for the computation (created if absent). When
            set it is the default for both ``data_path`` and ``cwd``, so lsdtt3
            reads/writes there.
        data_path: the ``<path>`` argument lsdtt3 expects before the parameter
            file. Defaults to ``directory``, else the parameter file's own dir.
        lsdtt_dir: directory containing the lsdtt-* executables; overrides the
            ``LSDTT3_DIR`` env var and config file for this call.
        extra_args: extra CLI flags appended to the invocation.
        cwd: working directory for the subprocess. Defaults to ``directory``.
        check: if True, raise ``CalledProcessError`` on a non-zero exit.
        clean_gdal_env: drop conda/pixi GDAL/PROJ env vars so the external
            lsdtt3 binary uses its own GDAL (avoids plugin/ABI clashes).

    Returns:
        The ``subprocess.CompletedProcess`` (stdout/stderr captured as text).
    """
    if directory is not None:
        from ..paths import ensure_dir

        directory = ensure_dir(directory)
        if data_path is None:
            data_path = directory
        if cwd is None:
            cwd = directory

    exe = resolve_lsdtt_tool(tool, lsdtt_dir)

    tmp_to_clean: Optional[Path] = None
    if isinstance(parameters, dict):
        out_dir = Path(data_path) if data_path else Path(cwd or ".")
        param_file = _write_parameter_file(parameters, out_dir)
        tmp_to_clean = param_file
    else:
        param_file = Path(parameters)
        if not param_file.is_file():
            raise FileNotFoundError(f"Parameter file not found: {param_file}")

    path_arg = str(data_path) if data_path else str(param_file.parent)
    argv = [exe, path_arg, param_file.name]
    if extra_args:
        argv.extend(extra_args)

    try:
        return subprocess.run(
            argv,
            cwd=cwd,
            check=check,
            capture_output=True,
            text=True,
            env=_external_env(clean_gdal_env),
        )
    finally:
        if tmp_to_clean is not None:
            tmp_to_clean.unlink(missing_ok=True)
