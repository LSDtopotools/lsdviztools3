"""Tests for the erodibility-calibration readers, E-k_sn chart, basin map and CLI.

Synthetic files follow lsdtt-chi-analysis' calibration headers exactly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lsdviztools3.io import erodibility_calibration as EC

mpl = pytest.importorskip("matplotlib")
mpl.use("Agg")

SAMPLE_HEADER = (
    "sample_id,input_row,input_latitude,input_longitude,snapped_latitude,snapped_longitude,"
    "snapping_method,snap_distance_m,outlet_flat,status,status_detail,reported_area_km2,"
    "dem_area_km2,area_ratio,area_check,nesting_group,nested_in,n_samples_inside,"
    "duplicate_outlet_of,erosion_rate_input,erosion_rate_uncertainty_input,"
    "erosion_rate_input_units,erosion_rate_m_per_yr,erosion_rate_uncertainty_m_per_yr,"
    "m_over_n,A_0,k_sn_statistic,k_sn,k_sn_ci_low,k_sn_ci_high,n_reaches,"
    "n_reaches_qualifying,k_sn_reach_iqr_over_median,k_sn_chi_z_slope,chi_z_r_squared,"
    "k_sn_segment_median,k_sn_chi_weighted_mean,dominant_lith_code,dominant_lith_name,"
    "dominant_area_fraction,dominant_channel_length_fraction,unmapped_area_fraction"
).split(",")

FIT_HEADER = (
    "fit_type,method,status,N_basins,N_groups,ksn_range_ratio,m_over_n,A_0,n,n_ci_low,"
    "n_ci_high,n_se_york,n_se_york_mswd_scaled,m,K,K_ci_low,K_ci_high,ln_K_se,"
    "ln_K_se_mswd_scaled,corr_n_lnK,E_at_median_ksn_m_per_yr,mswd,ln_K_dispersion,"
    "corr_lnK_lnk,n_bootstrap,bootstrap_cluster,K_units"
).split(",")


def _samples(a0=1.0) -> pd.DataFrame:
    rows = []
    specs = [("A1", "used", 20.0, 5e-6, "shale"), ("A2", "used", 40.0, 1e-5, "sandstone"),
             ("A3", "used", 60.0, 1.5e-5, "shale"), ("A4", "used", 30.0, 8e-6, "limestone"),
             ("B1", "too_few_reaches", 50.0, 4e-6, "shale"),
             ("C1", "incomplete_basin", np.nan, 6e-6, np.nan),
             ("C2", "area_mismatch", np.nan, 7e-6, "dolostone")]
    for i, (sid, st, k, e, unit) in enumerate(specs):
        r = {h: np.nan for h in SAMPLE_HEADER}
        r.update(sample_id=sid, input_row=i, status=st, snapped_latitude=39.0 + 0.002 * i,
                 snapped_longitude=-79.0 + 0.002 * i, erosion_rate_m_per_yr=e,
                 erosion_rate_uncertainty_m_per_yr=0.2 * e, m_over_n=0.45, A_0=a0,
                 k_sn_statistic="reach_median", k_sn=k, k_sn_ci_low=0.8 * k,
                 k_sn_ci_high=1.2 * k, dominant_lith_name=unit, nesting_group=i)
        rows.append(r)
    return pd.DataFrame(rows, columns=SAMPLE_HEADER)


def _fit(a0=1.0, fitted_status="ok") -> pd.DataFrame:
    rows = []
    for kind, n, K, st in (("fixed_n", 1.0, 2.5e-7, "ok"), ("fixed_n", 2.0, 7e-9, "ok"),
                           ("fitted_n", 0.9, 3e-7, fitted_status)):
        r = {h: np.nan for h in FIT_HEADER}
        r.update(fit_type=kind, method="york" if kind == "fitted_n" else "weighted_log_mean",
                 status=st, N_basins=4, N_groups=4, m_over_n=0.45, A_0=a0, n=n,
                 m=0.45 * n, K=K, K_ci_low=0.5 * K, K_ci_high=2 * K, n_bootstrap=1000,
                 bootstrap_cluster="nesting_group", K_units=f"m^({1 - 0.9 * n:g}) yr^-1")
        if kind == "fitted_n":
            r.update(n_ci_low=0.5, n_ci_high=1.4, mswd=2.5, corr_n_lnK=-0.98)
        rows.append(r)
    return pd.DataFrame(rows, columns=FIT_HEADER)


@pytest.fixture
def run_dir(tmp_path):
    _samples().to_csv(tmp_path / "cal_erodibility_calibration_samples_mn_0p450.csv",
                      index=False)
    _fit().to_csv(tmp_path / "cal_erodibility_calibration_fit_mn_0p450.csv", index=False)
    return tmp_path


def test_load_outputs(run_dir):
    out = EC.load_calibration_outputs(run_dir, "cal")
    assert out["mn"] == "0p450" and len(out["samples"]) == 7
    assert out["basins"] is None and out["by_lithology"] is None


def test_fit_line_is_lsdtt3_formula():
    from lsdviztools3.charts.erodibility_calibration import fit_line

    row = _fit(a0=2.0).iloc[2]
    k = np.array([10.0, 100.0])
    expected = row["K"] * 2.0 ** row["m"] * k ** row["n"]
    assert np.allclose(fit_line(row, k), expected)
    assert np.allclose(fit_line(_fit().iloc[0], k), 2.5e-7 * k)


def test_erosion_rate_vs_ksn_plot():
    from matplotlib.colors import to_hex

    from lsdviztools3.charts.erodibility_calibration import NO_BAND_NOTE, plot_erosion_rate_vs_ksn
    from lsdviztools3.render.lithology_colours import lithology_palette

    names = ["shale", "sandstone", "limestone", "dolostone"]
    fig, ax = plot_erosion_rate_vs_ksn(_samples(), _fit(), unit_names=names)
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"
    labels = ax.get_legend_handles_labels()[1]
    fitted = [lbl for lbl in labels if lbl.startswith("fitted n")]
    assert fitted == ["fitted n = 0.9 [0.5, 1.4], K = 3e-07 [1.5e-07, 6e-07] "
                      "m^(0.19) yr^-1, MSWD 2.5"]
    assert any(lbl.startswith("fixed n = 1,") for lbl in labels)
    assert any(lbl.startswith("fixed n = 2,") for lbl in labels)
    assert "not used (status as labelled)" in labels
    # No band: nothing filled between lines.
    assert not any(c.__class__.__name__ == "FillBetweenPolyCollection" or
                   c.__class__.__name__ == "PolyCollection" for c in ax.collections)
    title = ax.get_title(loc="left").replace("\n", " ")
    assert NO_BAND_NOTE.split(":")[0] in title
    assert "C1 (incomplete_basin)" in title and "C2 (area_mismatch)" in title
    assert "corr(n, ln K) -0.98" in title
    # Points coloured by the shared unit palette.
    pal = lithology_palette(names)
    shale = [h for h, lbl in zip(*ax.get_legend_handles_labels()) if lbl == "shale"][0]
    assert to_hex(shale.lines[0].get_color()) == pal["shale"]
    assert any(t.get_text() == "too_few_reaches" for t in ax.texts)
    mpl.pyplot.close(fig)


def test_fit_not_ok_is_not_drawn():
    from lsdviztools3.charts.erodibility_calibration import plot_erosion_rate_vs_ksn

    fig, ax = plot_erosion_rate_vs_ksn(_samples(), _fit(fitted_status="not_fitted_too_few_groups"))
    labels = ax.get_legend_handles_labels()[1]
    assert not any(lbl.startswith("fitted n") for lbl in labels)
    assert "status not_fitted_too_few_groups, not drawn" in ax.get_title(loc="left").replace("\n", " ")
    mpl.pyplot.close(fig)


def test_erosion_rate_credit(tmp_path):
    p = tmp_path / "oct_ATTRIBUTION.txt"
    p.write_text(
        "Basin-averaged cosmogenic erosion rates from OCTOPUS (https://octopusdata.org),\n"
        "provided by the University of Wollongong and licensed CC BY 4.0.\n\n"
        "Original studies (STUDYID, AUTHOR, PUBYEAR, REFDOI, n samples):\n"
        "  S046, Duxburry, 2009, (no DOI; REFID phdthesis), 11\n"
        "  S186, Portenga, 2019, 10.1130/b31840.1, 17\n\nCite the original studies.\n")
    assert EC.erosion_rate_credit(p) == ("Erosion rates: OCTOPUS (CC BY 4.0); Duxburry 2009; "
                                         "Portenga 2019")
    (tmp_path / "x_ingested_params.param").write_text("erosion_rate_fname: oct_erosion_rates.csv\n")
    assert EC.discover_erosion_rate_attribution(tmp_path, "x") == p


def _basins_and_dem(tmp_path):
    import geopandas as gpd
    import rioxarray  # noqa: F401
    import xarray as xr
    from shapely.geometry import Polygon

    s = _samples()
    pts = gpd.GeoSeries(gpd.points_from_xy(s.snapped_longitude, s.snapped_latitude),
                        crs="EPSG:4326").to_crs("EPSG:32617")
    polys = []
    for sid, p in zip(s.sample_id, pts):
        if sid == "C1":  # incomplete: no polygon
            continue
        polys.append((sid, Polygon([(p.x, p.y), (p.x + 300, p.y), (p.x + 300, p.y + 300),
                                    (p.x, p.y + 300)])))
    g = gpd.GeoDataFrame({"sample_id": [a for a, _ in polys], "status": "kept",
                          "nesting_group": 0}, geometry=[b for _, b in polys], crs="EPSG:32617")
    x0, y0 = pts.x.min() - 1000, pts.y.max() + 1500
    x = x0 + 30.0 * np.arange(120) + 15
    y = y0 - 30.0 * np.arange(120) - 15
    dem = xr.DataArray(np.add.outer(np.arange(120), np.arange(120)).astype("float32"),
                       dims=("y", "x"), coords={"y": y, "x": x}).rio.write_crs("EPSG:32617")
    dem.rio.to_raster(tmp_path / "cal_DEM.tif")
    g.to_file(tmp_path / "cal_erodibility_calibration_basins.fgb", driver="FlatGeobuf")
    return g


@pytest.mark.parametrize("colour_by", ["status", "unit"])
@pytest.mark.parametrize("with_dem", [True, False])
def test_basin_map(run_dir, monkeypatch, colour_by, with_dem):
    import pygmt

    from lsdviztools3.render import decorations
    from lsdviztools3.render.erodibility_calibration_map import render_calibration_basin_map
    from lsdviztools3.render.style import MapStyle

    g = _basins_and_dem(run_dir)
    legends, calls = [], {}
    orig_leg = pygmt.Figure.legend

    def _leg(self, spec=None, **k):
        legends.append(spec.getvalue())
        return orig_leg(self, spec=spec, **k)

    monkeypatch.setattr(pygmt.Figure, "legend", _leg)
    for name in ("add_bar_scalebar", "add_north_arrow", "add_map_credit"):
        orig = getattr(decorations, name)
        monkeypatch.setattr(decorations, name,
                            lambda *a, _o=orig, _n=name, **k: (calls.__setitem__(_n, k), _o(*a, **k))[1])
    fig = render_calibration_basin_map(
        _samples(), g, dem=run_dir / "cal_DEM.tif" if with_dem else None, colour_by=colour_by,
        credit="Geology: X\nErosion rates: OCTOPUS (CC BY 4.0)\nDEM: Y",
        style=MapStyle(clean=True, transparent=False))
    fig.savefig(run_dir / f"m_{colour_by}_{with_dem}.png", dpi=50)
    spec = legends[-1]
    if colour_by == "status":
        assert "used" in spec and "incomplete_basin" in spec  # outlet-only status listed
    else:
        assert "shale" in spec and "no dominant unit" in spec
    assert calls["add_bar_scalebar"]["label_position"] == "below"
    assert "add_north_arrow" in calls and calls["add_map_credit"]["anchor"] == "TL"


def test_map_credit_multiline():
    from lsdviztools3.render.decorations import add_map_credit

    class F:
        def __init__(self):
            self.t = []

        def text(self, **k):
            self.t.append(k)

    f = F()
    add_map_credit(f, [0, 10, 0, 10], "line one\nline two", anchor="TL")
    assert [k["text"] for k in f.t] == ["line one", "line two"]
    assert "offset" not in f.t[0] and f.t[1]["offset"].startswith("0c/-")


def test_cli(run_dir):
    from lsdviztools3.cli import plot_erodibility_calibration as cli

    _basins_and_dem(run_dir)
    out = run_dir / "figs"
    assert cli.main(["--data-dir", str(run_dir), "--prefix", "cal", "--out-dir", str(out),
                     "--dem", str(run_dir / "cal_DEM.tif"), "--dpi", "50",
                     "--credit", "test credit"]) == 0
    names = sorted(p.name for p in out.glob("*.png"))
    assert names == ["cal_basin_map_status.png", "cal_basin_map_unit.png",
                     "cal_erosion_rate_vs_ksn.png"]
