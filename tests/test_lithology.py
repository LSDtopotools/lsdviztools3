"""Tests for the lithology IO, k_sn-by-unit charts and the lithology unit map.

All inputs are small synthetic files written in the shape lsdtt3 writes them
(header-keyed CSVs, an int code raster, an ATTRIBUTION.txt); nothing here needs
the network or lsdtt3.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lsdviztools3.io import lithology as L

ATTRIBUTION_TEXT = """Geological map data from Macrostrat (https://macrostrat.org), licensed CC BY 4.0.
Cite Macrostrat: Peters, S. E., Husson, J. M., & Czaplewski, J. (2018). Macrostrat: A
platform for geological data integration and deep-time Earth crust research.
Geochemistry, Geophysics, Geosystems, 19, 1393-1409. doi:10.1029/2018GC007467

and each source map used:

- source_id 133: State Geologic Map Compilation
  Horton, J.D., C.A. San Juan, and D.B. Stoeser (2017)
  The State Geologic Map Compilation (SGMC) geodatabase of the conterminous United States
  U.S. Geological Survey Data Series 1052
  doi: 10.3133/ds1052
  https://pubs.er.usgs.gov/publication/ds1052
"""

UNITS = {1: "fine alluvium", 2: "limestone", 3: "sandstone"}


def _summary(extra_cols=None) -> pd.DataFrame:
    rows = []
    for code, name, n, med, low in [(1, "fine alluvium", 12, 5.0, 1),
                                    (2, "limestone", 300, 45.0, 0),
                                    (3, "sandstone", 500, 58.0, 0)]:
        rows.append({"lith_code": code, "lith_name": name, "n_nodes": n,
                     "n_segments": n // 6, "n_basins": 4, "channel_length_m": 30.0 * n,
                     "k_sn_median": med, "k_sn_q25": med * 0.6, "k_sn_q75": med * 1.4,
                     "k_sn_mean": med * 1.1, "k_sn_min": 0.0, "k_sn_max": med * 3,
                     "k_sn_length_weighted_mean": med, "low_sample": low,
                     "m_over_n": 0.45, "A_0": 1})
    rows.append({"lith_code": -1, "lith_name": "unmapped", "n_nodes": 0, "n_segments": 0,
                 "n_basins": 0, "channel_length_m": 0, "k_sn_median": np.nan,
                 "k_sn_q25": np.nan, "k_sn_q75": np.nan, "k_sn_mean": np.nan,
                 "k_sn_min": np.nan, "k_sn_max": np.nan, "k_sn_length_weighted_mean": np.nan,
                 "low_sample": 1, "m_over_n": 0.45, "A_0": 1})
    rows.append({"lith_code": -2, "lith_name": "all", "n_nodes": 812, "n_segments": 135,
                 "n_basins": 9, "channel_length_m": 24360.0, "k_sn_median": 50.0,
                 "k_sn_q25": 20.0, "k_sn_q75": 80.0, "k_sn_mean": 52.0, "k_sn_min": 0.0,
                 "k_sn_max": 174.0, "k_sn_length_weighted_mean": 51.0, "low_sample": 0,
                 "m_over_n": 0.45, "A_0": 1})
    df = pd.DataFrame(rows)
    for k, v in (extra_cols or {}).items():
        df[k] = v
    return df


def _histogram() -> pd.DataFrame:
    edges = np.linspace(0.0, 100.0, 6)
    rows = []
    for code, name in UNITS.items():
        rows.append([code, name, -1, -np.inf, edges[0], 2, 60.0])
        for i in range(5):
            rows.append([code, name, i, edges[i], edges[i + 1], 10 + i, 300.0 + 30 * i])
        rows.append([code, name, 5, edges[-1], np.inf, 3 if code != 1 else 0,
                     90.0 if code != 1 else 0.0])
    df = pd.DataFrame(rows, columns=["lith_code", "lith_name", "bin_index", "bin_lower",
                                     "bin_upper", "n_nodes", "channel_length_m"])
    df["m_over_n"], df["A_0"] = 0.45, 1
    return df


def _nodes() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    rows = []
    for code, name in UNITS.items():
        for i in range(40):
            rows.append({"basin_index": i % 3, "node_flat": i, "row": 0, "col": i,
                         "x": 1000.0 + 30 * i, "y": 2000.0 + 30 * code, "chi": 0.1 * i,
                         "elevation": 100 + i, "segment_id": i // 5,
                         "k_sn": float(rng.uniform(1, 90)),
                         "flow_length_to_receiver_m": 30.0, "lith_code": code,
                         "lith_name": name, "distance_to_contact_m": 60.0})
    return pd.DataFrame(rows)


@pytest.fixture
def run_dir(tmp_path) -> Path:
    """A directory laid out like an lsdtt3 lithology chi run (prefix ``run``)."""
    _summary().to_csv(tmp_path / "run_chi_ksn_by_lithology_mn_0p450.csv", index=False)
    _histogram().to_csv(tmp_path / "run_chi_ksn_histogram_by_lithology_mn_0p450.csv",
                        index=False)
    _nodes().to_csv(tmp_path / "run_chi_ksn_lithology_nodes_mn_0p450.csv", index=False)
    (tmp_path / "fetch_ATTRIBUTION.txt").write_text(ATTRIBUTION_TEXT)
    (tmp_path / "run_ingested_params.param").write_text(
        "# lsdtt3 parameter file\nread_fname: fetch_DEM.tif\n"
        "lithology_fname: fetch_macrostrat_units.fgb\nm_over_n: 0.45\n")
    (tmp_path / "fetch_ingested_params.param").write_text(
        "read_fname: dem\nwrite_fname: fetch\ndem_source: cop30_aws\n")
    return tmp_path


# --------------------------------------------------------------------------
# IO: names
# --------------------------------------------------------------------------

def test_names_from_chi_csv(run_dir):
    names = L.resolve_code_names(run_dir, "run")
    assert names[2] == "limestone"
    assert names[-1] == "unmapped"
    assert -2 not in names  # the pooled "all" row is not a unit


def test_names_units_csv_takes_priority(run_dir):
    pd.DataFrame({"code": [1, 2, 3], "name": ["A", "B", "C"], "area_km2": [1, 2, 3]}).to_csv(
        run_dir / "run_lithology_units.csv", index=False)
    names = L.resolve_code_names(run_dir, "run")
    assert names[2] == "B"


def test_names_explicit_csv_wins(run_dir, tmp_path):
    p = tmp_path / "explicit.csv"
    pd.DataFrame({"lith_code": [2], "lith_name": ["dolostone"]}).to_csv(p, index=False)
    assert L.resolve_code_names(run_dir, "run", names_csv=p)[2] == "dolostone"


def test_names_fall_back_to_crosswalk_then_legend(tmp_path):
    cw = tmp_path / "cw.csv"
    pd.DataFrame({"value": ["shale"], "code": [7], "name": ["shale"]}).to_csv(cw, index=False)
    assert L.resolve_code_names(tmp_path, "nothing", crosswalk=cw)[7] == "shale"
    pd.DataFrame({"code": [9], "name": ["chert"], "K": [1e-5]}).to_csv(
        tmp_path / "nothing_lithology_legend.csv", index=False)
    assert L.resolve_code_names(tmp_path, "nothing")[9] == "chert"


def test_names_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        L.resolve_code_names(tmp_path, "absent")


def test_names_read_by_header_not_position(tmp_path):
    p = tmp_path / "x.csv"
    pd.DataFrame({"extra": [0], "lith_name": ["marl"], "new_col": [1],
                  "lith_code": [4]}).to_csv(p, index=False)
    assert L.read_code_names(p) == {4: "marl"}


# --------------------------------------------------------------------------
# IO: attribution / credit
# --------------------------------------------------------------------------

def test_parse_attribution():
    a = L.parse_attribution(ATTRIBUTION_TEXT)
    assert a["provider"] == "Macrostrat"
    assert a["licence"] == "CC BY 4.0"
    assert len(a["sources"]) == 1
    src = a["sources"][0]
    assert src["year"] == "2017" and src["doi"] == "10.3133/ds1052"
    assert src["authors"].startswith("Horton")


def test_condensed_credit_line():
    line = L.condense_attribution(L.parse_attribution(ATTRIBUTION_TEXT),
                                  dem_source="Copernicus GLO-30")
    assert line == ("Geology: Macrostrat (CC BY 4.0); USGS SGMC, Horton et al. 2017. "
                    "DEM: Copernicus GLO-30")


def test_credit_discovery_from_run(run_dir):
    attr = L.discover_attribution(run_dir, "run")
    assert attr is not None and attr.name == "fetch_ATTRIBUTION.txt"
    dem_params = L.discover_dem_params(run_dir / "fetch_DEM.tif")
    assert L.dem_source_from_params(dem_params) == "Copernicus GLO-30"
    line = L.lithology_map_credit(attribution=attr, dem_params=dem_params)
    assert "Macrostrat (CC BY 4.0)" in line and line.endswith("DEM: Copernicus GLO-30")


def test_explicit_credit_overrides():
    assert L.lithology_map_credit(attribution={"provider": "X"}, credit="Mine") == "Mine"


def test_dem_source_providers(tmp_path):
    p = tmp_path / "p.param"
    p.write_text("dem_source: opentopography\not_dem_type: SRTMGL1\n")
    assert L.dem_source_from_params(p) == "SRTM GL1"
    p.write_text("dem_source: planetary_computer\npc_collection: cop-dem-glo-90\n")
    assert L.dem_source_from_params(p) == "Copernicus GLO-90"
    p.write_text("m_over_n: 0.45\n")
    assert L.dem_source_from_params(p) is None


# --------------------------------------------------------------------------
# IO: tables
# --------------------------------------------------------------------------

def test_load_outputs_and_optional_tables(run_dir):
    out = L.load_ksn_lithology_outputs(run_dir, "run")
    assert out["m_over_n"] == pytest.approx(0.45)
    assert out["reaches"] is None and out["unit_basin"] is None
    h = out["histogram"]
    assert np.isneginf(h["bin_lower"]).sum() == 3 and np.isposinf(h["bin_upper"]).sum() == 3


def test_mn_suffix_and_discovery(run_dir):
    assert L.mn_suffix(0.45) == "0p450"
    assert L.find_mn_suffixes(run_dir, "run") == ["0p450"]
    assert L.lithology_paths(run_dir, "run", 0.45)["nodes"].exists()


def test_missing_required_column_raises(tmp_path):
    p = tmp_path / "bad_chi_ksn_by_lithology_mn_0p450.csv"
    pd.DataFrame({"lith_code": [1], "lith_name": ["a"]}).to_csv(p, index=False)
    with pytest.raises(KeyError, match="n_nodes"):
        L.load_ksn_lithology_summary(p)


# --------------------------------------------------------------------------
# Charts
# --------------------------------------------------------------------------

mpl = pytest.importorskip("matplotlib")
mpl.use("Agg")


def test_histogram_overflow_not_drawn_as_bars(run_dir):
    from lsdviztools3.charts.lithology_ksn import (
        histogram_overflow, plot_ksn_histogram_by_lithology,
    )

    h, s = _histogram(), _summary()
    ov = histogram_overflow(h[h.lith_code == 2], "n_nodes")
    assert ov["below"] == (0.0, 2.0) and ov["above"] == (100.0, 3.0)
    for weight in ("length", "count"):
        fig, axes = plot_ksn_histogram_by_lithology(h, s, weight=weight, layout="overlay")
        lo, hi = axes[0].get_xlim()
        assert np.isfinite(lo) and np.isfinite(hi) and hi == pytest.approx(100.0)
        labels = axes[0].get_legend_handles_labels()[1]
        assert any("outside" in lbl and ">100" in lbl for lbl in labels)
        assert any("low sample" in lbl for lbl in labels)
        mpl.pyplot.close(fig)


def test_histogram_normalised_uses_summary_totals():
    from lsdviztools3.charts.lithology_ksn import plot_ksn_histogram_by_lithology

    fig, axes = plot_ksn_histogram_by_lithology(_histogram(), _summary(), weight="count",
                                                layout="overlay", order="code")
    # Unit 2 (limestone) is drawn second; its first inner bin is 10 of 300 nodes.
    step = [p for p in axes[0].patches][1]
    assert step.get_data().values[0] == pytest.approx(10 / 300)
    mpl.pyplot.close(fig)


def test_histogram_facets():
    from lsdviztools3.charts.lithology_ksn import plot_ksn_histogram_by_lithology

    fig, axes = plot_ksn_histogram_by_lithology(_histogram(), _summary(), layout="facets",
                                                ncols=2)
    assert len(axes) == 3
    mpl.pyplot.close(fig)


def _all_text(ax):
    """Annotation texts plus tick labels (sample sizes live in either)."""
    return " | ".join([t.get_text().replace("\n", ", ") for t in ax.texts]
                      + [t.get_text().replace("\n", ", ") for t in ax.get_xticklabels()]
                      + [t.get_text() for t in ax.get_yticklabels()])


@pytest.mark.parametrize("source", ["nodes", "summary"])
def test_box_plot_annotates_sample_sizes(source):
    from lsdviztools3.charts.lithology_ksn import plot_ksn_box_by_lithology

    fig, ax = plot_ksn_box_by_lithology(_summary(), _nodes(), source=source)
    text = _all_text(ax)
    assert "n=500 nodes, 83 seg, 4 basins" in text
    assert "(low sample)" in text
    assert "node-weighted" in ax.get_title() or source == "summary"
    ticks = [t.get_text().split("\n")[0] for t in ax.get_xticklabels()]
    assert ticks == ["fine alluvium", "limestone", "sandstone"]  # by summary median
    mpl.pyplot.close(fig)


def test_box_plot_summary_uses_lsdtt3_quartiles():
    from lsdviztools3.charts.lithology_ksn import plot_ksn_box_by_lithology

    fig, ax = plot_ksn_box_by_lithology(_summary(), source="summary", order="code")
    # Median lines are horizontal segments at the summary median (5, 45, 58);
    # no whiskers beyond the summary quartiles are drawn.
    flat_y = {round(float(ln.get_ydata()[0]), 6) for ln in ax.lines
              if len(set(np.round(ln.get_ydata(), 6))) == 1}
    assert {5.0, 45.0, 58.0} <= flat_y
    all_y = np.concatenate([ln.get_ydata() for ln in ax.lines])
    assert all_y.min() >= 3.0 and all_y.max() <= 58.0 * 1.4 + 1e-9
    mpl.pyplot.close(fig)


def test_reach_box_and_ci_when_present():
    from lsdviztools3.charts.lithology_ksn import plot_reach_ksn_by_lithology

    reaches = pd.DataFrame({
        "reach_id": range(9), "basin_index": [0] * 9, "lith_code": [1, 2, 2, 2, 3, 3, 3, 3, 1],
        "lith_name": ["fine alluvium", "limestone", "limestone", "limestone", "sandstone",
                      "sandstone", "sandstone", "sandstone", "fine alluvium"],
        "k_sn_ols": [4, 40, 50, 45, 60, 55, 70, 65, 6.0], "k_sn_ols_se": 1.0,
        "qualifies": [1, 1, 1, 1, 1, 1, 1, 0, 0], "m_over_n": 0.45})
    s = _summary({"reach_k_sn_median": [4, 45, 60, np.nan, 50],
                  "reach_k_sn_median_ci_low": [3, 40, 55, np.nan, 45],
                  "reach_k_sn_median_ci_high": [5, 50, 65, np.nan, 55],
                  "n_reaches": [2, 3, 4, 0, 9], "n_reaches_qualifying": [1, 3, 3, 0, 7]})
    fig, ax = plot_reach_ksn_by_lithology(reaches, s)
    assert "3 of 4 reaches" in _all_text(ax)
    assert "CI" in ax.get_legend_handles_labels()[1][0]
    mpl.pyplot.close(fig)


def test_unit_basin_dot_plot_uses_reach_median():
    from lsdviztools3.charts.lithology_ksn import plot_ksn_by_lithology_and_basin

    ub = pd.DataFrame({"lith_code": [2, 2, 3, 3],
                       "lith_name": ["limestone", "limestone", "sandstone", "sandstone"],
                       "basin_index": [0, 1, 0, 1], "n_nodes": [100, 5, 200, 10],
                       "k_sn_median": [1.0, 1.0, 1.0, 1.0],
                       "n_reaches": [8, 1, 9, 2], "n_reaches_qualifying": [6, 1, 7, 0],
                       "reach_k_sn_median": [40.0, 50.0, 60.0, np.nan],
                       "low_sample": [0, 1, 0, 1]})
    s = _summary({"reach_k_sn_median": [4, 45, 60, np.nan, 50]})
    fig, ax = plot_ksn_by_lithology_and_basin(ub, s)
    assert [t.get_text() for t in ax.get_yticklabels()] == ["limestone", "sandstone"]
    xs = sorted(float(x) for coll in ax.collections for x, _y in coll.get_offsets())
    assert xs == [40.0, 50.0, 60.0]  # reach medians; NaN basin not drawn; node medians unused
    assert "2 basins drawn" in _all_text(ax) and "1 basin drawn" in _all_text(ax)
    mpl.pyplot.close(fig)


def test_reach_ci_from_summary_not_recomputed():
    """The error bar spans exactly the summary's ci_low..ci_high."""
    from lsdviztools3.charts.lithology_ksn import plot_reach_ksn_by_lithology

    reaches = pd.DataFrame({"reach_id": range(4), "lith_code": [2, 2, 3, 3],
                            "lith_name": ["limestone", "limestone", "sandstone", "sandstone"],
                            "k_sn_ols": [40, 50, 60, 70.0], "qualifies": [1.0, 1.0, 1.0, 1.0]})
    s = _summary({"reach_k_sn_median": [np.nan, 45, 65, np.nan, 55],
                  "reach_k_sn_median_ci_low": [np.nan, 11, 22, np.nan, 33],
                  "reach_k_sn_median_ci_high": [np.nan, 99, 88, np.nan, 77],
                  "n_reaches": [0, 2, 2, 0, 4], "n_reaches_qualifying": [0, 2, 2, 0, 4],
                  "n_bootstrap": 1000, "bootstrap_cluster": "reach",
                  "contact_buffer_metres": 120})
    fig, ax = plot_reach_ksn_by_lithology(reaches, s, orientation="vertical")
    segs = [c for c in ax.collections if c.__class__.__name__ == "LineCollection"]
    spans = sorted(tuple(sorted(np.asarray(seg)[:, 1])) for c in segs for seg in c.get_segments())
    assert (11.0, 99.0) in spans and (22.0, 88.0) in spans
    assert "1000 bootstrap replicates" in ax.get_title()
    assert "contact buffer 120 m" in ax.get_title()
    mpl.pyplot.close(fig)


def test_real_schema_optional_columns(run_dir):
    out = L.load_ksn_lithology_outputs(run_dir, "run")
    assert out["segments"] is None and out["contact_buffer_m"] is None


def test_ksn_cli_writes_figures(run_dir, tmp_path):
    from lsdviztools3.cli import plot_ksn_by_lithology

    out = tmp_path / "figs"
    assert plot_ksn_by_lithology.main(["--data-dir", str(run_dir), "--prefix", "run",
                                       "--out-dir", str(out), "--dpi", "60"]) == 0
    names = sorted(p.name for p in out.glob("*.png"))
    assert names == ["run_ksn_box.png", "run_ksn_hist_count.png", "run_ksn_hist_length.png"]


# --------------------------------------------------------------------------
# Lithology map (PyGMT)
# --------------------------------------------------------------------------

def _code_raster(tmp_path, *, with_unmapped=True) -> Path:
    import rioxarray  # noqa: F401
    import xarray as xr

    ny, nx = 40, 50
    codes = np.ones((ny, nx), dtype="float32")
    codes[:, 15:30] = 2
    codes[:, 30:] = 3
    if with_unmapped:
        codes[:5, :5] = -1
    codes[-3:, -3:] = -9999
    x = 600000.0 + 30.0 * np.arange(nx) + 15
    y = 4300000.0 - 30.0 * np.arange(ny) - 15
    da = xr.DataArray(codes, dims=("y", "x"), coords={"y": y, "x": x})
    da = da.rio.write_crs("EPSG:32617").rio.write_nodata(-9999.0)
    p = tmp_path / "run_lithology_code.tif"
    da.rio.to_raster(p)
    dem = xr.DataArray((np.add.outer(np.arange(ny), np.arange(nx)) * 3.0).astype("float32"),
                       dims=("y", "x"), coords={"y": y, "x": x}).rio.write_crs("EPSG:32617")
    dem.rio.to_raster(tmp_path / "fetch_DEM.tif")
    return p


def test_categorical_entries_use_names_and_refuse_bare_codes():
    from lsdviztools3.render.categorical import categorical_legend_entries

    entries = categorical_legend_entries([3, -1, 1], {1: "a", 3: ("c", "red")})
    assert [(c, lbl) for c, lbl, _ in entries] == [(1, "a"), (3, "c"), (-1, "unmapped")]
    assert entries[1][2] == "red"
    with pytest.raises(KeyError, match="No name"):
        categorical_legend_entries([1, 2], {1: "a"})


@pytest.mark.parametrize("hillshade", [True, False])
@pytest.mark.parametrize("channels", [None, "csv", "lines_by_order"])
def test_render_lithology_map_combinations(run_dir, hillshade, channels, monkeypatch):
    """Every hillshade x channel combination draws legend, credit, bar and arrow."""
    import geopandas as gpd
    from shapely.geometry import LineString

    from lsdviztools3.render import decorations, lithology_map
    from lsdviztools3.render.style import MapStyle

    code = _code_raster(run_dir)
    calls = {}
    for name in ("add_bar_scalebar", "add_north_arrow", "add_map_credit"):
        orig = getattr(decorations, name)

        def _spy(*a, _orig=orig, _name=name, **k):
            calls[_name] = k
            return _orig(*a, **k)

        monkeypatch.setattr(decorations, name, _spy)
    legend_calls = []
    from lsdviztools3.render import categorical
    orig_leg = categorical.add_categorical_legend

    def _leg(fig, entries, **k):
        legend_calls.append((entries, k))
        return orig_leg(fig, entries, **k)

    monkeypatch.setattr(categorical, "add_categorical_legend", _leg)

    chan = None
    if channels == "csv":
        chan = pd.DataFrame({"x": [600100.0, 600800.0], "y": [4299500.0, 4299200.0]})
    elif channels == "lines_by_order":
        chan = gpd.GeoDataFrame(
            {"strahler_order": [1, 2]},
            geometry=[LineString([(600100, 4299500), (600500, 4299300)]),
                      LineString([(600500, 4299300), (601200, 4299000)])],
            crs="EPSG:32617")
    names = {1: "fine alluvium", 2: "limestone", 3: "sandstone", -1: "unmapped"}
    fig = lithology_map.render_lithology_map(
        code, dem=run_dir / "fetch_DEM.tif" if hillshade else None, names=names,
        hillshade=hillshade, channels=chan, channel_order_column="auto",
        credit="Geology: test", style=MapStyle(clean=True, transparent=False))
    out = run_dir / f"map_{hillshade}_{channels}.png"
    fig.savefig(out, dpi=60)
    assert out.stat().st_size > 0
    assert calls["add_bar_scalebar"]["label_position"] == "below"
    assert "add_north_arrow" in calls and "add_map_credit" in calls
    entries, kwargs = legend_calls[-1]
    assert [lbl for _c, lbl, _col in entries] == ["fine alluvium", "limestone",
                                                  "sandstone", "unmapped"]
    extra = " ".join(kwargs.get("extra_rows", []))
    if channels == "csv":
        assert "Channels" in extra
    elif channels == "lines_by_order":
        assert "Order 1" in extra and "Order 2" in extra
    else:
        assert extra == ""


def test_lithology_map_cli(run_dir):
    from lsdviztools3.cli import plot_lithology_map

    code = _code_raster(run_dir)
    out = run_dir / "cli_map.png"
    args = ["--codes", str(code), "--dem", str(run_dir / "fetch_DEM.tif"),
            "--out", str(out), "--clean", "--opaque", "--dpi", "60",
            "--channels", str(run_dir / "run_chi_ksn_lithology_nodes_mn_0p450.csv")]
    assert plot_lithology_map.main(args) == 0
    assert out.stat().st_size > 0
    resolved = plot_lithology_map.resolve_inputs(plot_lithology_map.build_parser().parse_args(args))
    assert resolved["prefix"] == "run"
    assert resolved["credit"] == ("Geology: Macrostrat (CC BY 4.0); USGS SGMC, Horton et al. "
                                  "2017. DEM: Copernicus GLO-30")


def test_lithology_map_cli_needs_dem_for_hillshade(run_dir):
    from lsdviztools3.cli import plot_lithology_map

    code = _code_raster(run_dir)
    with pytest.raises(SystemExit):
        plot_lithology_map.main(["--codes", str(code), "--out", str(run_dir / "x.png")])


@pytest.mark.parametrize("relief,expected", [("illumination", "lightgray"), ("elevation", "gray")])
def test_lithology_map_relief_base(run_dir, monkeypatch, relief, expected):
    """illumination = uniform base colour (brightness from slope only)."""
    from lsdviztools3.render import categorical, lithology_map

    seen = {}
    orig = categorical.render_categorical_drape

    def _spy(*a, **k):
        seen["base_cmap"] = k.get("base_cmap")
        return orig(*a, **k)

    monkeypatch.setattr(categorical, "render_categorical_drape", _spy)
    code = _code_raster(run_dir)
    from lsdviztools3.render.style import MapStyle

    lithology_map.render_lithology_map(code, dem=run_dir / "fetch_DEM.tif",
                                       names={1: "a", 2: "b", 3: "c", -1: "unmapped"},
                                       relief=relief, style=MapStyle(clean=True))
    assert seen["base_cmap"] == expected


def test_channel_order_map(run_dir):
    import geopandas as gpd
    from shapely.geometry import LineString

    from lsdviztools3.render.lithology_map import render_channel_order_map
    from lsdviztools3.render.style import MapStyle

    _code_raster(run_dir)
    chan = gpd.GeoDataFrame(
        {"stream_order": [1, 2]},
        geometry=[LineString([(600100, 4299500), (600500, 4299300)]),
                  LineString([(600500, 4299300), (601200, 4299000)])], crs="EPSG:32617")
    fig = render_channel_order_map(run_dir / "fetch_DEM.tif", chan,
                                   style=MapStyle(clean=True, transparent=False))
    out = run_dir / "order_map.png"
    fig.savefig(out, dpi=50)
    assert out.stat().st_size > 0


# --------------------------------------------------------------------------
# Shared name-keyed unit colours (maps and charts)
# --------------------------------------------------------------------------

def test_palette_keyed_by_name_not_order_or_code():
    from lsdviztools3.render.lithology_colours import code_colours, lithology_palette

    a = lithology_palette(["shale", "sandstone", "unmapped", "all"])
    b = lithology_palette(["unmapped", "sandstone", "shale"])
    assert a == b and "all" not in a
    assert a["unmapped"] == "#c8c8ca"
    # Same names under different run-local codes -> same colours.
    assert code_colours({1: "shale", 2: "sandstone"}) == {1: a["shale"], 2: a["sandstone"]}
    assert code_colours({30: "shale", 10: "sandstone"})[30] == a["shale"]
    assert lithology_palette(["shale"], {"shale": "red"})["shale"] == "#ff0000"
    assert lithology_palette(["shale"], {"shale": "0/0/255"})["shale"] == "#0000ff"


def test_read_unit_colours_csv_and_json(tmp_path):
    p = tmp_path / "c.csv"
    pd.DataFrame({"name": ["shale"], "colour": ["#123456"]}).to_csv(p, index=False)
    assert L.read_unit_colours(p) == {"shale": "#123456"}
    j = tmp_path / "c.json"
    j.write_text('{"Tonoloway, Wills Creek": "red"}')
    assert L.read_unit_colours(j) == {"Tonoloway, Wills Creek": "red"}


@pytest.mark.parametrize("overrides", [None, {"limestone": "#00ff00"}])
def test_map_legend_and_chart_colours_agree(run_dir, monkeypatch, overrides):
    """Same unit name -> same colour on the map legend and in the k_sn charts."""
    from lsdviztools3.charts.lithology_ksn import chart_unit_colours, plot_ksn_box_by_lithology
    from lsdviztools3.render import categorical, lithology_map
    from lsdviztools3.render.style import MapStyle

    legend = {}
    orig = categorical.add_categorical_legend

    def _leg(fig, entries, **k):
        legend.update({lbl: col for _c, lbl, col in entries})
        return orig(fig, entries, **k)

    monkeypatch.setattr(categorical, "add_categorical_legend", _leg)
    code = _code_raster(run_dir)
    names = L.resolve_code_names(run_dir, "run")  # from the chi summary
    lithology_map.render_lithology_map(code, names=names, hillshade=False,
                                       unit_colours=overrides, style=MapStyle(clean=True))
    summary = _summary()
    chart = chart_unit_colours(summary, overrides)
    for name in ("fine alluvium", "limestone", "sandstone", "unmapped"):
        assert legend[name] == chart[name], name
    if overrides:
        assert legend["limestone"] == "#00ff00"
    # The box faces actually use those colours.
    fig, ax = plot_ksn_box_by_lithology(summary, _nodes(), unit_colours=overrides, order="code")
    from matplotlib.colors import to_hex

    faces = [to_hex(p.get_facecolor(), keep_alpha=False) for p in ax.patches[:3]]
    assert faces == [chart["fine alluvium"], chart["limestone"], chart["sandstone"]]
    mpl.pyplot.close(fig)


# --------------------------------------------------------------------------
# k_sn-coloured channels
# --------------------------------------------------------------------------

def test_ksn_colour_limits_percentile_and_log():
    from lsdviztools3.render.lithology_map import ksn_colour_limits, ksn_colourbar_label

    v = np.arange(0.0, 101.0)
    z, lo, hi = ksn_colour_limits(v, (2, 98))
    assert (lo, hi) == (2.0, 98.0) and np.array_equal(z, v)
    z, lo, hi = ksn_colour_limits(np.array([-5.0, 0.0, 10.0, 100.0]), (0, 100), log_scale=True)
    assert np.isnan(z[0]) and np.isnan(z[1]) and (lo, hi) == (1.0, 2.0)
    lab = ksn_colourbar_label(0.45, 1.0)
    assert "m@+0.9@+" in lab and "m/n = 0.45" in lab and "A@-0@- = 1" in lab


@pytest.mark.parametrize("hillshade", [True, False])
def test_lithology_map_ksn_channels(run_dir, monkeypatch, hillshade):
    from lsdviztools3.render import categorical, decorations, lithology_map
    from lsdviztools3.render.style import MapStyle

    bars, legend, drape = [], [], {}
    orig_bar, orig_leg = decorations.add_colourbar, categorical.add_categorical_legend
    orig_drape = categorical.render_categorical_drape
    monkeypatch.setattr(decorations, "add_colourbar",
                        lambda fig, **k: (bars.append(k), orig_bar(fig, **k))[1])
    monkeypatch.setattr(categorical, "add_categorical_legend",
                        lambda fig, e, **k: (legend.append((e, k)), orig_leg(fig, e, **k))[1])

    def _d(*a, **k):
        drape.update(k)
        return orig_drape(*a, **k)

    monkeypatch.setattr(categorical, "render_categorical_drape", _d)
    code = _code_raster(run_dir)
    nodes = pd.DataFrame({"x": 600100.0 + 30 * np.arange(20), "y": 4299500.0,
                          "k_sn": np.linspace(1, 200, 20), "m_over_n": 0.45, "A_0": 1})
    fig = lithology_map.render_lithology_map(
        code, dem=run_dir / "fetch_DEM.tif" if hillshade else None,
        names={1: "a", 2: "b", 3: "c", -1: "unmapped"}, hillshade=hillshade,
        channels=nodes, channel_color_by="k_sn", style=MapStyle(clean=True))
    fig.savefig(run_dir / f"ksn_{hillshade}.png", dpi=50)
    assert "m/n = 0.45" in bars[0]["label"] and "2-98th pct" in bars[0]["label"]
    entries, kwargs = legend[-1]
    assert all(col.endswith("@65") for _c, _l, col in entries)  # muted swatches
    assert "k_sn" in " ".join(kwargs["extra_rows"])
    if hillshade:
        assert drape["transparency"] == 65.0
    else:
        assert drape["flat_transparency"] == 65.0


def test_lithology_map_cli_ksn_and_unit_colours(run_dir):
    from lsdviztools3.cli import plot_lithology_map

    code = _code_raster(run_dir)
    uc = run_dir / "colours.csv"
    pd.DataFrame({"name": ["limestone"], "colour": ["#00ff00"]}).to_csv(uc, index=False)
    out = run_dir / "ksn_cli.png"
    args = ["--codes", str(code), "--no-hillshade", "--out", str(out), "--clean", "--dpi", "50",
            "--channels", str(run_dir / "run_chi_ksn_lithology_nodes_mn_0p450.csv"),
            "--channel-color-by", "k_sn", "--ksn-log", "--unit-colours", str(uc)]
    assert plot_lithology_map.main(args) == 0 and out.stat().st_size > 0
    r = plot_lithology_map.resolve_inputs(plot_lithology_map.build_parser().parse_args(args))
    assert r["ksn_m_over_n"] == pytest.approx(0.45) and r["unit_colours"] == {"limestone": "#00ff00"}


def test_palette_reserved_mixed_and_order():
    from lsdviztools3.render.lithology_colours import (
        RESERVED_HEX, lithology_palette, order_unit_names,
    )

    names = ["unmapped", "mixed (<70% one unit)", "shale", "arms disagree", "Basalt", "all"]
    assert order_unit_names(names) == ["Basalt", "shale", "mixed (<70% one unit)",
                                       "arms disagree", "unmapped"]
    pal = lithology_palette(names)
    assert pal["mixed (<70% one unit)"] == RESERVED_HEX["mixed"]
    assert len({RESERVED_HEX["mixed"], RESERVED_HEX["unmapped"], RESERVED_HEX["arms_disagree"],
                pal["shale"], pal["Basalt"]}) == 5
    # Adding a reserved category does not shift unit colours.
    assert lithology_palette(["shale", "Basalt"])["shale"] == pal["shale"]


def test_categorical_entries_reserved_codes_last():
    from lsdviztools3.render.categorical import categorical_legend_entries
    from lsdviztools3.render.lithology_colours import RESERVED_HEX

    e = categorical_legend_entries([-1, -3, 30, 10, -4],
                                   {10: "a", 30: "b", -3: "mixed (<80% one unit)",
                                    -4: "arms disagree", -1: "unmapped"})
    assert [c for c, _l, _col in e] == [10, 30, -3, -4, -1]
    assert e[2][2] == RESERVED_HEX["mixed"] and e[3][2] == RESERVED_HEX["arms_disagree"]


def test_ksn_charts_order_mixed_last():
    from lsdviztools3.charts.lithology_ksn import chart_unit_colours, unit_rows
    from lsdviztools3.render.lithology_colours import RESERVED_HEX

    s = _summary()
    extra = s.iloc[[1]].copy()
    extra["lith_code"], extra["lith_name"], extra["k_sn_median"] = -3, "mixed (<80% one unit)", 1.0
    s = pd.concat([s, extra], ignore_index=True)
    rows = unit_rows(s, order="median")
    assert rows["lith_name"].iloc[-1] == "mixed (<80% one unit)"  # lowest median, still last
    assert chart_unit_colours(s)["mixed (<80% one unit)"] == RESERVED_HEX["mixed"]
