"""Smoke tests for CLI entrypoints."""

from __future__ import annotations

import subprocess

import pytest

from lsdviztools3.cli import (
    grab_dem,
    plot_basemap,
    plot_basins,
    plot_channels,
    plot_hillshade,
    run_lsdtt,
)


@pytest.mark.parametrize(
    "build_parser",
    [
        plot_hillshade.build_parser,
        plot_channels.build_parser,
        plot_basins.build_parser,
        plot_basemap.build_parser,
        grab_dem.build_parser,
        run_lsdtt.build_parser,
    ],
)
def test_cli_help_exits_cleanly(build_parser):
    with pytest.raises(SystemExit) as exc:
        build_parser().parse_args(["--help"])
    assert exc.value.code == 0


def test_plot_hillshade_main_smoke(monkeypatch):
    calls = {}

    def fake_render(path, **kwargs):
        calls["render"] = (path, kwargs)
        return object()

    def fake_save(fig, out_path, *, style, high_res=False):
        calls["save"] = (fig, out_path, style, high_res)

    monkeypatch.setattr(plot_hillshade, "render_hillshade", fake_render)
    monkeypatch.setattr(plot_hillshade._common, "save_figure", fake_save)

    code = plot_hillshade.main(["--in", "dem.tif", "--out", "hill.png", "--opaque"])
    assert code == 0
    assert calls["render"][0] == "dem.tif"
    assert calls["save"][1] == "hill.png"


def test_plot_channels_main_smoke(monkeypatch):
    calls = {}

    class DummyFig:
        pass

    def fake_hillshade(path, **kwargs):
        calls["hillshade"] = (path, kwargs)
        return DummyFig()

    def fake_channels(path, **kwargs):
        calls["channels"] = (path, kwargs)
        return kwargs["fig"]

    def fake_save(fig, out_path, *, style, high_res=False):
        calls["save"] = (fig, out_path, style, high_res)

    monkeypatch.setattr(plot_channels, "render_hillshade", fake_hillshade)
    monkeypatch.setattr(plot_channels, "render_channels", fake_channels)
    monkeypatch.setattr(plot_channels._common, "save_figure", fake_save)

    code = plot_channels.main([
        "--base",
        "dem.tif",
        "--lines",
        "lines.fgb",
        "--out",
        "channels.png",
        "--color-by",
        "stream_order",
    ])
    assert code == 0
    assert calls["channels"][0] == "lines.fgb"
    assert calls["save"][1] == "channels.png"


def test_plot_basins_main_smoke(monkeypatch):
    calls = {}

    class DummyFig:
        pass

    def fake_hillshade(path, **kwargs):
        calls["hillshade"] = (path, kwargs)
        return DummyFig()

    def fake_basins(path, **kwargs):
        calls["basins"] = (path, kwargs)
        return kwargs["fig"]

    def fake_save(fig, out_path, *, style, high_res=False):
        calls["save"] = (fig, out_path, style, high_res)

    monkeypatch.setattr(plot_basins, "render_hillshade", fake_hillshade)
    monkeypatch.setattr(plot_basins, "render_basins", fake_basins)
    monkeypatch.setattr(plot_basins._common, "save_figure", fake_save)

    code = plot_basins.main([
        "--base",
        "dem.tif",
        "--basins",
        "basins.fgb",
        "--out",
        "basins.png",
        "--label-col",
        "id",
    ])
    assert code == 0
    assert calls["basins"][0] == "basins.fgb"
    assert calls["save"][1] == "basins.png"


def test_plot_basemap_main_smoke(monkeypatch):
    calls = {}

    def fake_basemap(path, **kwargs):
        calls["basemap"] = (path, kwargs)
        return object()

    def fake_save(fig, out_path, *, style, high_res=False):
        calls["save"] = (fig, out_path, style, high_res)

    monkeypatch.setattr(plot_basemap, "render_orthographic_basemap", fake_basemap)
    monkeypatch.setattr(plot_basemap._common, "save_figure", fake_save)

    code = plot_basemap.main(
        [
            "--in",
            "dem.tif",
            "--out",
            "basemap.png",
            "--mode",
            "aerial",
            "--elevation-km",
            "50",
            "--look-bearing",
            "0",
            "--look-vertical-angle",
            "0",
        ]
    )
    assert code == 0
    assert calls["basemap"][0] == "dem.tif"
    assert calls["basemap"][1]["mode"] == "aerial"
    assert calls["save"][1] == "basemap.png"


def test_grab_dem_main_smoke(monkeypatch, capsys):
    calls = {}

    def fake_download_dem(bbox, **kwargs):
        calls["download"] = (bbox, kwargs)
        return "out.tif"

    monkeypatch.setattr(grab_dem, "download_dem", fake_download_dem)

    code = grab_dem.main(
        ["--bbox", "-6.2", "56.8", "-6.1", "56.9", "--out", "out.tif"]
    )
    assert code == 0
    assert calls["download"][0] == [-6.2, 56.8, -6.1, 56.9]
    assert "out.tif" in capsys.readouterr().out


def test_run_lsdtt_set_dir_smoke(monkeypatch, capsys):
    monkeypatch.setattr(run_lsdtt, "write_user_config", lambda _: "cfg.toml")
    code = run_lsdtt.main(["--set-lsdtt-dir", "/opt/lsdtt3/bin"])
    assert code == 0
    assert "cfg.toml" in capsys.readouterr().out


def test_run_lsdtt_main_smoke(monkeypatch, capsys):
    result = subprocess.CompletedProcess(["tool"], 0, "ok\n", "")
    monkeypatch.setattr(run_lsdtt, "run_lsdtt", lambda *a, **k: result)

    code = run_lsdtt.main(["lsdtt-channel-extraction", "--params", "job.driver"])
    assert code == 0
    assert "ok" in capsys.readouterr().out
