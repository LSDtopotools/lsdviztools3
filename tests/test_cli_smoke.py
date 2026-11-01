"""Smoke tests for CLI entrypoints."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from lsdviztools3.cli import (
    grab_dem,
    plot_basemap,
    plot_basins,
    plot_channels,
    plot_chi_cosmo,
    plot_cosmo_basins,
    plot_cosmo_compare,
    plot_cosmo_published,
    plot_divide_migration,
    plot_drape,
    plot_channel_bearings,
    plot_erodibility_calibration,
    plot_hillshade,
    plot_ksn_by_lithology,
    plot_lithology_map,
    plot_sti,
    plot_perimeter_sources,
    plot_valley_widths,
    run_lsdtt,
)


@pytest.mark.parametrize(
    "build_parser",
    [
        plot_hillshade.build_parser,
        plot_drape.build_parser,
        plot_sti.build_parser,
        plot_channels.build_parser,
        plot_basins.build_parser,
        plot_divide_migration.build_parser,
        plot_basemap.build_parser,
        plot_chi_cosmo.build_parser,
        plot_cosmo_basins.build_parser,
        plot_cosmo_compare.build_parser,
        plot_cosmo_published.build_parser,
        grab_dem.build_parser,
        run_lsdtt.build_parser,
        plot_lithology_map.build_parser,
        plot_ksn_by_lithology.build_parser,
        plot_channel_bearings.build_parser,
        plot_erodibility_calibration.build_parser,
        plot_valley_widths.build_parser,
        plot_perimeter_sources.build_parser,
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

    def fake_save(fig, out_path, *, style, high_res=False, directory=None):
        calls["save"] = (fig, out_path, style, high_res)

    monkeypatch.setattr(plot_hillshade, "render_hillshade", fake_render)
    monkeypatch.setattr(plot_hillshade._common, "save_figure", fake_save)

    code = plot_hillshade.main(["--in", "dem.tif", "--out", "hill.png", "--opaque"])
    assert code == 0
    assert calls["render"][0] == "dem.tif"
    assert calls["save"][1] == "hill.png"


def test_plot_drape_main_smoke(monkeypatch):
    calls = {}

    def fake_drape(base, drape, **kwargs):
        calls["drape"] = (base, drape, kwargs)
        return object()

    def fake_save(fig, out_path, *, style, high_res=False, directory=None):
        calls["save"] = (fig, out_path, style, high_res)

    monkeypatch.setattr(plot_drape, "render_drape", fake_drape)
    monkeypatch.setattr(plot_drape._common, "save_figure", fake_save)

    code = plot_drape.main([
        "--base", "hs.tif",
        "--drape", "curv.tif",
        "--out", "drape.png",
        "--cmap", "roma",
        "--log",
        "--transparency", "40",
        "--colorbar",
    ])
    assert code == 0
    assert calls["drape"][0] == "hs.tif"
    assert calls["drape"][1] == "curv.tif"
    assert calls["drape"][2]["log_color"] is True
    assert calls["drape"][2]["transparency"] == 40.0
    assert calls["drape"][2]["cmap"] == "roma"
    assert calls["save"][1] == "drape.png"


def test_plot_sti_main_smoke(monkeypatch):
    calls = {}

    def fake_sti(base, sti_path, **kwargs):
        calls["sti"] = (base, sti_path, kwargs)
        return object()

    def fake_save(fig, out_path, *, style, high_res=False, directory=None):
        calls["save"] = (fig, out_path, style, high_res)

    monkeypatch.setattr(plot_sti, "render_sti_p95_drape", fake_sti)
    monkeypatch.setattr(plot_sti._common, "save_figure", fake_save)

    code = plot_sti.main([
        "--base", "hs.tif",
        "--sti", "sti_p95.tif",
        "--out", "sti.png",
        "--transparency", "30",
    ])
    assert code == 0
    assert calls["sti"][0] == "hs.tif"
    assert calls["sti"][1] == "sti_p95.tif"
    assert calls["sti"][2]["transparency"] == 30.0
    assert calls["sti"][2]["cmap"] == "YlOrRd"
    assert calls["save"][1] == "sti.png"


def test_plot_sti_writes_non_p95_partner_when_present(monkeypatch, tmp_path):
    calls = {"drape": [], "save": []}

    def fake_sti(base, sti_path, **kwargs):
        return object()

    def fake_drape(base, drape, **kwargs):
        calls["drape"].append((base, drape, kwargs))
        return object()

    def fake_save(fig, out_path, *, style, high_res=False, directory=None):
        calls["save"].append(out_path)

    monkeypatch.setattr(plot_sti, "render_sti_p95_drape", fake_sti)
    monkeypatch.setattr(plot_sti, "render_drape", fake_drape)
    monkeypatch.setattr(plot_sti._common, "save_figure", fake_save)

    sti_p95 = tmp_path / "Kildavie_sti_p95.tif"
    sti_full = tmp_path / "Kildavie_sti.tif"
    sti_p95.write_text("x")
    sti_full.write_text("x")

    code = plot_sti.main([
        "--base", "hs.tif",
        "--sti", str(sti_p95),
        "--out", "sti.png",
    ])
    assert code == 0
    assert len(calls["drape"]) == 1
    assert calls["drape"][0][1].endswith("Kildavie_sti.tif")
    assert calls["save"][-1].endswith("sti_full.png")


def test_plot_sti_directory_mode_plots_all_found(monkeypatch, tmp_path):
    calls = []

    def fake_intermediate(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr(plot_sti, "_render_intermediate", fake_intermediate)

    (tmp_path / "Kildavie_sti.tif").write_text("x")
    (tmp_path / "Kildavie_sti_p95.tif").write_text("x")
    (tmp_path / "Kildavie_sti_twi.tif").write_text("x")
    (tmp_path / "Kildavie_sti_series_codes.tif").write_text("x")
    (tmp_path / "Kildavie_sti_landuse_codes.tif").write_text("x")
    (tmp_path / "Kildavie_sti_valid_mask_ksd.tif").write_text("x")
    (tmp_path / "Kildavie_sti_valid_mask_spr.tif").write_text("x")
    (tmp_path / "Kildavie_sti.tif.tmp_cog.tif").write_text("x")

    code = plot_sti.main([
        "--sti-dir", str(tmp_path),
        "--base", "hs.tif",
    ])
    assert code == 0
    assert len(calls) == 3
    assert any(c["drape_raster"].endswith("Kildavie_sti.tif") for c in calls)
    assert any(c["drape_raster"].endswith("Kildavie_sti_p95.tif") for c in calls)
    assert any(c["drape_raster"].endswith("Kildavie_sti_twi.tif") for c in calls)
    assert all(not c["drape_raster"].endswith("_sti_series_codes.tif") for c in calls)
    assert all(not c["drape_raster"].endswith("_sti_landuse_codes.tif") for c in calls)
    assert all(not c["drape_raster"].endswith("_sti_valid_mask_ksd.tif") for c in calls)
    assert all(not c["drape_raster"].endswith("_sti_valid_mask_spr.tif") for c in calls)

    by_name = {Path(c["drape_raster"]).name: c for c in calls}
    assert by_name["Kildavie_sti.tif"]["colorbar"] is True
    assert by_name["Kildavie_sti_twi.tif"]["cmap"] == "YlOrRd"


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

    def fake_save(fig, out_path, *, style, high_res=False, directory=None):
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
        """Records the title text the CLI writes onto the figure."""

        def __init__(self):
            self.texts = []

        def text(self, **kwargs):
            self.texts.append(kwargs)

    def fake_hillshade(path, **kwargs):
        calls["hillshade"] = (path, kwargs)
        return DummyFig()

    def fake_basins(path, **kwargs):
        calls["basins"] = (path, kwargs)
        return kwargs["fig"]

    def fake_save(fig, out_path, *, style, high_res=False, directory=None):
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
    # The CLI writes one title line carrying the cosmo-scheme metadata.
    texts = calls["save"][0].texts
    assert len(texts) == 1 and "scaling=" in texts[0]["text"]


def test_plot_basemap_main_smoke(monkeypatch):
    calls = {}

    def fake_basemap(path, **kwargs):
        calls["basemap"] = (path, kwargs)
        return object()

    def fake_save(fig, out_path, *, style, high_res=False, directory=None):
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


def test_plot_cosmo_compare_main_smoke(monkeypatch, capsys):
    def fake_render(**kwargs):
        assert kwargs["scan_comparison_csv"] == "scan.csv"
        assert kwargs["sample_comparison_csv"] == "samples.csv"
        return {
            "sweep_panels": "out/sweep.png",
            "sample_residuals": "out/samples.png",
        }

    monkeypatch.setattr(plot_cosmo_compare, "render_cosmo_comparison_figures", fake_render)
    code = plot_cosmo_compare.main(
        [
            "--scan-csv",
            "scan.csv",
            "--sample-csv",
            "samples.csv",
            "--outdir",
            "out",
            "--reference",
            "cronuscalc",
            "--pressure",
            "900",
            "--erate",
            "0.01",
            "--stem",
            "cmp",
        ]
    )
    assert code == 0
    out = capsys.readouterr().out
    assert "sweep_panels" in out and "sample_residuals" in out


def test_plot_cosmo_basins_main_smoke(monkeypatch, capsys):
    def fake_render(**kwargs):
        assert kwargs["base_raster"] == "dem.tif"
        assert kwargs["basin_polygons"] == "basins.fgb"
        assert kwargs["cosmo_csv"] == "erate.csv"
        return {
            "out_png": "out/basins.png",
            "joined_basins": "out/joined.fgb",
            "n_joined_basins": 4,
        }

    monkeypatch.setattr(plot_cosmo_basins, "render_cosmo_basin_erate_map", fake_render)
    code = plot_cosmo_basins.main(
        [
            "--base",
            "dem.tif",
            "--basins",
            "basins.fgb",
            "--cosmo-csv",
            "erate.csv",
            "--out",
            "out/basins.png",
            "--joined-basins-out",
            "out/joined.fgb",
            "--colorbar",
        ]
    )
    assert code == 0
    out = capsys.readouterr().out
    assert "n_joined_basins" in out and "joined_basins" in out


def test_plot_chi_cosmo_main_smoke(monkeypatch, capsys):
    def fake_render(outputs, out_dir, **kwargs):
        assert outputs["fit_nodes"] == "fit_nodes.csv"
        assert outputs["fit_segments"] == "fit_segments.csv"
        assert str(out_dir).endswith("out")
        assert kwargs["base_raster"] == "dem.tif"
        assert kwargs["basin_polygons"] == "basins.fgb"
        assert kwargs["cosmo_csv"] == "erate.csv"
        return {
            "chi_stacked": "out/chi_stacked.png",
            "chi_profiles": ["out/chi_profile_basin0.png"],
            "ksn_summary_csv": "out/chi_basin_ksn_summary.csv",
        }

    monkeypatch.setattr(plot_chi_cosmo, "render_chi_and_cosmo_basin_workflow", fake_render)
    code = plot_chi_cosmo.main(
        [
            "--fit-nodes-csv",
            "fit_nodes.csv",
            "--fit-segments-csv",
            "fit_segments.csv",
            "--outdir",
            "out",
            "--base",
            "dem.tif",
            "--basins",
            "basins.fgb",
            "--cosmo-csv",
            "erate.csv",
        ]
    )
    assert code == 0
    out = capsys.readouterr().out
    assert "chi_stacked" in out and "ksn_summary_csv" in out


def test_plot_cosmo_published_main_smoke(monkeypatch, capsys):
    def fake_render(**kwargs):
        assert kwargs["published_dir"] == "pub"
        assert kwargs["datasets"] == ["dethier2014", "clapp2002"]
        assert kwargs["mode"] == "auto"
        return {
            "sample_residuals": "out/cosmo_published_residuals.png",
            "source_csvs": ["pub/dethier2014_summary.csv", "pub/clapp2002_summary.csv"],
        }

    monkeypatch.setattr(plot_cosmo_published, "render_cosmo_published_residuals", fake_render)
    code = plot_cosmo_published.main(
        [
            "--published-dir",
            "pub",
            "--datasets",
            "dethier2014,clapp2002",
            "--mode",
            "auto",
            "--outdir",
            "out",
        ]
    )
    assert code == 0
    out = capsys.readouterr().out
    assert "sample_residuals" in out and "source_csvs" in out
