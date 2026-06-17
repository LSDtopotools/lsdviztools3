"""Tests for lsdviztools3.driver — lsdtt3 location config and invocation."""

import os
import stat
import subprocess

import pytest

from lsdviztools3.driver import config, lsdtt_driver


def _make_fake_tool(directory, name="lsdtt-channel-extraction"):
    """Create a dummy executable in ``directory`` and return its path."""
    exe = directory / name
    exe.write_text("#!/bin/sh\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return exe


def test_resolve_dir_explicit_arg_wins(tmp_path, monkeypatch):
    monkeypatch.setenv(config.ENV_LSDTT_DIR, "/should/not/win")
    assert config.resolve_lsdtt_dir(tmp_path) == tmp_path


def test_resolve_dir_from_env(tmp_path, monkeypatch):
    monkeypatch.setenv(config.ENV_LSDTT_DIR, str(tmp_path))
    assert config.resolve_lsdtt_dir() == tmp_path


def test_resolve_dir_from_config_file(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    cfg = tmp_path / "lsdviztools3.toml"
    cfg.write_text(f'[lsdtt]\ndir = "{bin_dir}"\n')
    monkeypatch.delenv(config.ENV_LSDTT_DIR, raising=False)
    monkeypatch.setenv(config.ENV_CONFIG, str(cfg))
    assert config.resolve_lsdtt_dir() == bin_dir


def test_resolve_tool_finds_executable(tmp_path, monkeypatch):
    monkeypatch.delenv(config.ENV_LSDTT_DIR, raising=False)
    exe = _make_fake_tool(tmp_path)
    assert config.resolve_lsdtt_tool("lsdtt-channel-extraction", tmp_path) == str(exe)


def test_resolve_tool_missing_raises_with_guidance(tmp_path, monkeypatch):
    monkeypatch.delenv(config.ENV_LSDTT_DIR, raising=False)
    with pytest.raises(FileNotFoundError) as exc:
        config.resolve_lsdtt_tool("lsdtt-channel-extraction", tmp_path)
    assert config.ENV_LSDTT_DIR in str(exc.value)


def test_run_lsdtt_builds_argv(tmp_path, monkeypatch):
    exe = _make_fake_tool(tmp_path)
    params = tmp_path / "job.driver"
    params.write_text("write_hillshade: true\n")

    captured = {}

    def fake_run(argv, **kwargs):
        captured["argv"] = argv
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(subprocess, "run", fake_run)
    lsdtt_driver.run_lsdtt("lsdtt-channel-extraction", params, lsdtt_dir=tmp_path)

    assert captured["argv"] == [str(exe), str(tmp_path), "job.driver"]


def test_run_lsdtt_serialises_dict(tmp_path, monkeypatch):
    _make_fake_tool(tmp_path)
    written = {}

    def fake_run(argv, **kwargs):
        # The temp .driver file still exists at call time; capture its contents.
        param_path = tmp_path / argv[2]
        written["contents"] = param_path.read_text()
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(subprocess, "run", fake_run)
    lsdtt_driver.run_lsdtt(
        "lsdtt-channel-extraction",
        {"write_hillshade": "true", "remove_seas": "true"},
        data_path=tmp_path,
        lsdtt_dir=tmp_path,
    )

    assert "write_hillshade: true" in written["contents"]
    assert "remove_seas: true" in written["contents"]


def test_write_user_config(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    path = config.write_user_config("/opt/lsdtt3/bin")
    assert path.read_text().strip() == '[lsdtt]\ndir = "/opt/lsdtt3/bin"'


def test_write_user_config_preserves_unrelated_keys(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    cfg_dir = tmp_path / ".config" / "lsdviztools3"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    cfg = cfg_dir / "config.toml"
    cfg.write_text('[download]\nsource = "planetary_computer"\n\n[lsdtt]\ndir = "/old/path"\n')

    config.write_user_config("/opt/lsdtt3/bin")
    text = cfg.read_text()
    assert '[download]\nsource = "planetary_computer"' in text
    assert 'dir = "/opt/lsdtt3/bin"' in text
