"""Tests for the target-directory helpers and CLI --dir plumbing."""

from __future__ import annotations

from pathlib import Path

from lsdviztools3.paths import ensure_dir, resolve_in_dir


def test_ensure_dir_creates_missing(tmp_path):
    target = tmp_path / "a" / "b" / "c"
    assert not target.exists()
    out = ensure_dir(target)
    assert out == target and target.is_dir()


def test_resolve_in_dir_joins_relative_and_creates(tmp_path):
    base = tmp_path / "run"
    out = resolve_in_dir(base, "fig.png")
    assert out == base / "fig.png"
    assert base.is_dir()  # created on resolve


def test_resolve_in_dir_none_keeps_path_but_makes_parent(tmp_path):
    p = tmp_path / "sub" / "fig.png"
    out = resolve_in_dir(None, p)
    assert out == p and p.parent.is_dir()


def test_resolve_in_dir_absolute_wins_but_dir_made(tmp_path):
    base = tmp_path / "run"
    absolute = tmp_path / "elsewhere" / "fig.png"
    out = resolve_in_dir(base, absolute)
    assert out == absolute
    assert base.is_dir()          # target dir still created
    assert absolute.parent.is_dir()


def test_grab_dem_cli_passes_directory(monkeypatch):
    from lsdviztools3.cli import grab_dem

    calls = {}

    def fake_download_dem(bbox, **kwargs):
        calls["kwargs"] = kwargs
        return "out.tif"

    monkeypatch.setattr(grab_dem, "download_dem", fake_download_dem)
    code = grab_dem.main(
        ["--bbox", "-6.2", "56.8", "-6.1", "56.9", "--out", "dem.tif",
         "--dir", "/tmp/lsdviz3_target"]
    )
    assert code == 0
    assert calls["kwargs"]["directory"] == "/tmp/lsdviz3_target"


def test_run_lsdtt_cli_passes_directory(monkeypatch):
    import subprocess

    from lsdviztools3.cli import run_lsdtt

    captured = {}

    def fake_run(*args, **kwargs):
        captured.update(kwargs)
        return subprocess.CompletedProcess(["tool"], 0, "", "")

    monkeypatch.setattr(run_lsdtt, "run_lsdtt", fake_run)
    code = run_lsdtt.main(
        ["lsdtt-channel-extraction", "--params", "job.driver",
         "--dir", "/tmp/lsdviz3_target"]
    )
    assert code == 0
    assert captured["directory"] == "/tmp/lsdviz3_target"
