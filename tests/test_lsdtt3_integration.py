"""Integration tests: drive lsdtt3 and render its outputs.

Skipped entirely unless the lsdtt3 executables are visible (set ``LSDTT3_DIR``
or configure a `[lsdtt] dir`). The workflow test also needs the network (DEM
download) and GMT/PyGMT.
"""

from pathlib import Path

import geopandas as gpd
import pytest

from lsdviztools3.driver import config

REQUIRED_TOOLS = [
    "lsdtt-raster-preprocessing",
    "lsdtt-channel-extraction",
    "lsdtt-basin-extraction",
]


def _lsdtt_available() -> bool:
    try:
        for tool in REQUIRED_TOOLS:
            config.resolve_lsdtt_tool(tool)
        return True
    except FileNotFoundError:
        return False


pytestmark = pytest.mark.skipif(
    not _lsdtt_available(),
    reason="lsdtt3 executables not found (set LSDTT3_DIR or [lsdtt] dir in config).",
)


def test_lsdtt3_clis_visible():
    """The lsdtt3 tools the workflow needs resolve to real executables."""
    for tool in REQUIRED_TOOLS:
        resolved = config.resolve_lsdtt_tool(tool)
        assert Path(resolved).exists()


@pytest.mark.network
def test_lsdtt3_workflow_renders(small_dem, tmp_path):
    """Full pipeline: hillshade + channels + largest basin -> composite PNG."""
    from lsdviztools3.workflows import (
        lsdtt3_hillshade_channels_basin,
        render_hillshade_channels_basin,
    )

    outputs = lsdtt3_hillshade_channels_basin(small_dem, tmp_path, stem="eigg")
    assert outputs["hillshade"].exists()
    assert outputs["channels"].exists()
    assert outputs["basin"].exists()

    # filter_all_but_largest_basin -> exactly one polygon.
    basin = gpd.read_file(outputs["basin"])
    assert len(basin) == 1
    assert basin.geometry.iloc[0].geom_type in ("Polygon", "MultiPolygon")

    png = tmp_path / "composite.png"
    render_hillshade_channels_basin(outputs, png)
    assert png.exists() and png.stat().st_size > 0


def _chi_available() -> bool:
    try:
        config.resolve_lsdtt_tool("lsdtt-chi-analysis")
        return True
    except FileNotFoundError:
        return False


def _run_and_check_chi(dem, workdir, *, stem, source_threshold, n_largest):
    """Drive chi analysis, then load + plot; assert the products are sane."""
    import matplotlib
    matplotlib.use("Agg")

    from lsdviztools3.io.chi import (
        FIT_NODE_COLUMNS, FIT_SEGMENT_COLUMNS,
        load_chi_fit, load_chi_fit_nodes, load_chi_fit_segments,
        load_chi_points, load_chi_segment_lines,
    )
    from lsdviztools3.charts import plot_chi_stacked
    from lsdviztools3.workflows import lsdtt3_chi_analysis, render_chi_profiles

    outputs = lsdtt3_chi_analysis(
        dem, workdir, stem=stem,
        source_threshold=source_threshold, n_largest_basins=n_largest,
        m_over_n=0.45,
    )
    assert outputs["fit_nodes"].exists()
    assert outputs["fit_segments"].exists()
    assert outputs["segment_lines"].exists()

    # Column contract (docs/chi_analysis.md, "Inputs: which lsdtt3 call produces which file"): every declared column is
    # present, in the declared order, in what lsdtt3 actually wrote. Extra
    # trailing columns are tolerated; a rename, removal or reorder on the C++
    # side fails here so the loaders and the docs table get updated with it.
    def _assert_declared_columns(path, declared, loader):
        actual = list(loader(path).columns)
        assert set(declared).issubset(actual), (
            f"{path.name} is missing {sorted(set(declared) - set(actual))}"
        )
        assert [c for c in actual if c in declared] == list(declared), (
            f"{path.name} reorders the declared columns: {actual}"
        )

    _assert_declared_columns(
        outputs["fit_nodes"], FIT_NODE_COLUMNS, load_chi_fit_nodes)
    _assert_declared_columns(
        outputs["fit_segments"], FIT_SEGMENT_COLUMNS, load_chi_fit_segments)
    # k_sn is a segments-only column until assign_node_ksn grafts it on.
    assert "k_sn" not in load_chi_fit_nodes(outputs["fit_nodes"]).columns

    nodes = load_chi_fit(outputs["fit_nodes"], outputs["fit_segments"])
    assert not nodes.empty
    for col in ("chi", "observed_elevation", "fitted_elevation",
                "basin_index", "branch_id", "k_sn"):
        assert col in nodes.columns

    n_basins = nodes["basin_index"].nunique()
    assert 1 <= n_basins <= n_largest
    # k_sn was attached to every node from its containing segment.
    assert nodes["k_sn"].notna().all()

    pngs = render_chi_profiles(outputs, workdir / "plots", stem=stem)
    # One stacked plot + one profile per basin.
    assert len(pngs) == 1 + n_basins
    for p in pngs:
        assert p.exists() and p.stat().st_size > 0

    # The plotters are also compatible with FlatGeobuf output: the fitted
    # segment-lines fgb carries chi/elevation endpoints + k_sn, so it can drive
    # the same stacked plot without the CSV node/segment tables.
    segs = load_chi_segment_lines(outputs["segment_lines"])
    assert not segs.empty
    for col in ("basin_index", "branch_id", "k_sn",
                "downstream_chi", "upstream_chi",
                "downstream_elevation", "upstream_elevation"):
        assert col in segs.columns
    # Observed points from the data-map fgb overlay the fitted segments (same
    # per-basin offset). Basin keys align with the segment basin_index.
    observed = load_chi_points(outputs["data_maps"])
    assert {"chi", "elevation", "basin_key"}.issubset(observed.columns)
    fig, _ = plot_chi_stacked(segs, observed=observed)
    fgb_png = workdir / "plots" / f"{stem}_stacked_fgb.png"
    fig.savefig(fgb_png, dpi=120, bbox_inches="tight")
    assert fgb_png.exists() and fgb_png.stat().st_size > 0
    return nodes


@pytest.mark.skipif(not _chi_available(),
                    reason="lsdtt-chi-analysis not found.")
@pytest.mark.network
def test_chi_profiles_eigg(small_dem, tmp_path):
    """Fast path: chi fit + stacked/individual plots on the small Eigg island."""
    _run_and_check_chi(
        small_dem, tmp_path, stem="eigg",
        source_threshold=100, n_largest=3,
    )


@pytest.mark.skipif(not _chi_available(),
                    reason="lsdtt-chi-analysis not found.")
@pytest.mark.network
@pytest.mark.slow
def test_chi_profiles_guadalcanal(guadalcanal_dem, tmp_path):
    """Chi fit + plots on Guadalcanal: the N-largest complete-basin landscape."""
    _run_and_check_chi(
        guadalcanal_dem, tmp_path, stem="guadalcanal",
        source_threshold=1000, n_largest=5,
    )


@pytest.mark.skipif(not _chi_available(),
                    reason="lsdtt-chi-analysis not found.")
@pytest.mark.network
def test_chi_ksn_map_eigg(small_dem, tmp_path):
    """Chi k_sn network coloured by log10 k_sn over a hillshade -> PNG."""
    from lsdviztools3.workflows import lsdtt3_chi_analysis, render_chi_ksn_map

    outputs = lsdtt3_chi_analysis(
        small_dem, tmp_path, stem="eigg",
        source_threshold=100, n_largest_basins=3, m_over_n=0.45,
    )
    dem_utm = tmp_path / "eigg_remove_seas.tif"
    assert dem_utm.exists() and outputs["segment_lines"].exists()

    png = render_chi_ksn_map(dem_utm, outputs["segment_lines"],
                             tmp_path / "eigg_ksn_map.png")
    assert png.exists() and png.stat().st_size > 0


def _swath_available() -> bool:
    try:
        config.resolve_lsdtt_tool("lsdtt-swath-profile")
        config.resolve_lsdtt_tool("lsdtt-chi-analysis")
        return True
    except FileNotFoundError:
        return False


def _run_and_check_swath(dem, workdir, *, stem, source_threshold,
                         half_width, bin_width):
    """Drive a swath profile along the main channel, then load + plot it."""
    import matplotlib
    matplotlib.use("Agg")

    from lsdviztools3.io.swath import load_swath, swath_percentile_columns
    from lsdviztools3.workflows import lsdtt3_swath_profile, render_swath_profile

    outputs = lsdtt3_swath_profile(
        dem, workdir, stem=stem, source_threshold=source_threshold,
        half_width_metres=half_width, bin_width_metres=bin_width,
    )
    assert outputs["swath_csv"].exists()
    assert outputs["baseline"].exists()
    assert outputs["along_axis"].exists()

    df = load_swath(outputs["swath_csv"])
    assert not df.empty
    for col in ("along_axis_centre", "n_pixels", "mean", "p0", "p50", "p100"):
        assert col in df.columns
    assert len(swath_percentile_columns(df)) == 21
    # Along-axis distance is monotonic and the envelope is ordered (p0 <= p100).
    assert df["along_axis_centre"].is_monotonic_increasing
    populated = df[df["n_pixels"] > 0]
    assert (populated["p0"] <= populated["p100"]).all()

    png = render_swath_profile(outputs, workdir / f"{stem}_swath.png")
    assert png.exists() and png.stat().st_size > 0


@pytest.mark.skipif(not _swath_available(),
                    reason="lsdtt-swath-profile / lsdtt-chi-analysis not found.")
@pytest.mark.network
def test_swath_profile_eigg(small_dem, tmp_path):
    """Fast path: swath along the Eigg main channel -> envelope plot."""
    _run_and_check_swath(
        small_dem, tmp_path, stem="eigg",
        source_threshold=100, half_width=500, bin_width=100,
    )


@pytest.mark.skipif(not _swath_available(),
                    reason="lsdtt-swath-profile / lsdtt-chi-analysis not found.")
@pytest.mark.network
@pytest.mark.slow
def test_swath_profile_guadalcanal(guadalcanal_dem, tmp_path):
    """Swath along Guadalcanal's main channel."""
    _run_and_check_swath(
        guadalcanal_dem, tmp_path, stem="guadalcanal",
        source_threshold=1000, half_width=2000, bin_width=500,
    )
