"""Loader and polar chart for lsdtt-spectral-analysis 2D significant frequencies."""

from __future__ import annotations

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

NEW_CSV = """\
# confidence,0.9999
# normalised_psd_threshold,9.21034
# background_seed,1000
# n_background_iterations,20
bearing_rad,bearing_deg,frequency_cycles_per_m,wavelength_m,normalised_psd
0.785398,45,0.5,2,12.5
3.92699,225,0.5,2,12.5
1.5708,90,0.25,4,30
4.71239,270,0.25,4,30
"""

# Exactly what lsdtt-spectral-analysis writes when nothing is significant.
HEADER_ONLY_CSV = """\
# confidence,0.9999
# normalised_psd_threshold,9.21034
# background_seed,1000
# n_background_iterations,20
# n_significant_frequencies,0
bearing_rad,bearing_deg,frequency_cycles_per_m,wavelength_m,normalised_psd
# No significant frequencies were found in this raster: no cell of the 2D spectrum exceeded the normalised PSD threshold of 9.21034 (confidence 0.9999). This is a valid result, not an error.
"""

# The same before n_significant_frequencies and the message line were added.
HEADER_ONLY_CSV_V1 = """\
# confidence,0.9999
# normalised_psd_threshold,9.21034
# background_seed,1000
# n_background_iterations,20
bearing_rad,bearing_deg,frequency_cycles_per_m,wavelength_m,normalised_psd
"""

# The legacy writer ended the file with a blank line.
LEGACY_CSV = """\
bearing_rad,bearing_deg,frequency,wavelength,normalised_power
0.244979,14.0362,0.901929,1.10873,10.0774
3.38657,194.036,0.901929,1.10873,10.0774

"""

LEGACY_EMPTY_CSV = "bearing_rad,bearing_deg,frequency,wavelength,normalised_power\n\n"


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


# --- loader -----------------------------------------------------------------

def test_load_new_format_reads_metadata_and_columns(tmp_path):
    from lsdviztools3.io.spectral import COLUMNS, load_significant_frequencies_2d

    _write(tmp_path, "tile_a_significant_frequencies_2d.csv", NEW_CSV)
    sig = load_significant_frequencies_2d(tmp_path / "tile_a")
    assert list(sig.table.columns) == COLUMNS
    assert sig.n_cells == 4
    assert sig.stem == "tile_a"
    assert not sig.legacy
    assert sig.threshold == pytest.approx(9.21034)
    assert sig.metadata["background_seed"] == 1000
    assert sig.metadata["n_background_iterations"] == 20


def test_header_only_is_an_empty_result_not_an_error(tmp_path):
    from lsdviztools3.io.spectral import COLUMNS, load_significant_frequencies_2d

    for text in (HEADER_ONLY_CSV, HEADER_ONLY_CSV_V1):
        p = _write(tmp_path, "t_significant_frequencies_2d.csv", text)
        sig = load_significant_frequencies_2d(p)
        assert sig.n_cells == 0
        assert list(sig.table.columns) == COLUMNS
        assert sig.threshold == pytest.approx(9.21034)
    assert sig.metadata == {"confidence": 0.9999, "normalised_psd_threshold": 9.21034,
                            "background_seed": 1000.0, "n_background_iterations": 20.0}
    p = _write(tmp_path, "t_significant_frequencies_2d.csv", HEADER_ONLY_CSV)
    assert load_significant_frequencies_2d(p).metadata["n_significant_frequencies"] == 0


def test_legacy_format_is_renamed_and_carries_legacy_threshold(tmp_path):
    from lsdviztools3.io.spectral import (COLUMNS, LEGACY_THRESHOLD,
                                          load_significant_frequencies_2d)

    _write(tmp_path, "old_significant_freq.csv", LEGACY_CSV)
    sig = load_significant_frequencies_2d(tmp_path / "old")
    assert sig.legacy
    assert list(sig.table.columns) == COLUMNS
    assert sig.n_cells == 2                      # trailing blank line ignored
    assert sig.threshold == pytest.approx(LEGACY_THRESHOLD)
    assert sig.table["wavelength_m"].iloc[0] == pytest.approx(1.10873)
    assert sig.stem == "old"

    _write(tmp_path, "empty_significant_freq.csv", LEGACY_EMPTY_CSV)
    assert load_significant_frequencies_2d(tmp_path / "empty").n_cells == 0


def test_directory_resolution_and_missing_file(tmp_path):
    from lsdviztools3.io.spectral import load_significant_frequencies_2d

    _write(tmp_path, "only_significant_frequencies_2d.csv", NEW_CSV)
    assert load_significant_frequencies_2d(tmp_path).stem == "only"
    with pytest.raises(FileNotFoundError):
        load_significant_frequencies_2d(tmp_path / "nope")


def test_sweep_and_summary_keep_empty_runs_only_in_summary(tmp_path):
    from lsdviztools3.io.spectral import (load_significant_frequencies_2d_sweep,
                                          summarise_significant_frequencies_2d)

    a = _write(tmp_path, "a_significant_frequencies_2d.csv", NEW_CSV)
    b = _write(tmp_path, "b_significant_frequencies_2d.csv", HEADER_ONLY_CSV)
    c = _write(tmp_path, "c_significant_freq.csv", LEGACY_CSV)

    long = load_significant_frequencies_2d_sweep([a, b, c])
    assert list(long["stem"].unique()) == ["a", "c"]
    assert len(long) == 6

    summary = summarise_significant_frequencies_2d([a, b, c]).set_index("stem")
    assert summary.loc["a", "n_cells"] == 4
    assert summary.loc["a", "peak_wavelength_m"] == pytest.approx(4.0)
    assert summary.loc["a", "peak_normalised_psd"] == pytest.approx(30.0)
    assert summary.loc["b", "n_cells"] == 0
    assert np.isnan(summary.loc["b", "peak_wavelength_m"])
    assert summary.loc["c", "peak_bearing_deg"] == pytest.approx(14.0362)


# --- polar chart ------------------------------------------------------------

def _offsets(ax):
    pts = ax.collections[0].get_offsets()
    return np.asarray(pts[:, 0]), np.asarray(pts[:, 1])


def test_polar_is_north_up_clockwise_with_points_at_bearing_and_frequency(tmp_path):
    from lsdviztools3.charts import plot_significant_frequencies_polar
    from lsdviztools3.io.spectral import load_significant_frequencies_2d

    sig = load_significant_frequencies_2d(
        _write(tmp_path, "t_significant_frequencies_2d.csv", NEW_CSV))
    fig, ax = plot_significant_frequencies_polar(sig)
    assert ax.name == "polar"
    assert ax.get_theta_direction() == -1
    assert ax.get_theta_offset() == pytest.approx(np.pi / 2)

    theta, r = _offsets(ax)
    got = sorted(zip(np.round(np.rad2deg(theta), 3), np.round(r, 6)))
    want = sorted([(45.0, 0.5), (225.0, 0.5), (90.0, 0.25), (270.0, 0.25)])
    assert got == want


def test_polar_crest_orientation_adds_90_and_half_keeps_one_of_each_pair(tmp_path):
    from lsdviztools3.charts import plot_significant_frequencies_polar
    from lsdviztools3.io.spectral import load_significant_frequencies_2d

    sig = load_significant_frequencies_2d(
        _write(tmp_path, "t_significant_frequencies_2d.csv", NEW_CSV))

    _, ax = plot_significant_frequencies_polar(sig, orientation="crest")
    theta, _ = _offsets(ax)
    assert sorted(np.round(np.rad2deg(theta), 3)) == [0.0, 135.0, 180.0, 315.0]

    _, ax = plot_significant_frequencies_polar(sig, half=True)
    theta, _ = _offsets(ax)
    assert sorted(np.round(np.rad2deg(theta), 3)) == [45.0, 90.0]


def test_polar_wavelength_radius_and_wavelength_ring_labels(tmp_path):
    from lsdviztools3.charts import plot_significant_frequencies_polar
    from lsdviztools3.io.spectral import load_significant_frequencies_2d

    sig = load_significant_frequencies_2d(
        _write(tmp_path, "t_significant_frequencies_2d.csv", NEW_CSV))

    _, ax = plot_significant_frequencies_polar(sig, radius="wavelength")
    _, r = _offsets(ax)
    assert sorted(np.round(r, 6)) == [2.0, 2.0, 4.0, 4.0]

    _, ax = plot_significant_frequencies_polar(sig, r_max=2.0)  # Nyquist at 0.25 m
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert labels and all(t.endswith(" m") for t in labels)
    assert ax.get_ylim()[1] == pytest.approx(2.0)


def test_polar_empty_run_draws_a_note_not_an_error(tmp_path):
    from lsdviztools3.charts import plot_significant_frequencies_polar
    from lsdviztools3.io.spectral import load_significant_frequencies_2d

    sig = load_significant_frequencies_2d(
        _write(tmp_path, "t_significant_frequencies_2d.csv", HEADER_ONLY_CSV))
    fig, ax = plot_significant_frequencies_polar(sig)
    assert not ax.collections
    assert any("no cells above threshold" in t.get_text() for t in ax.texts)


def test_polar_accepts_a_sweep_slice_and_rejects_bad_options(tmp_path):
    from lsdviztools3.charts import plot_significant_frequencies_polar
    from lsdviztools3.io.spectral import load_significant_frequencies_2d_sweep

    a = _write(tmp_path, "a_significant_frequencies_2d.csv", NEW_CSV)
    long = load_significant_frequencies_2d_sweep([a])
    _, ax = plot_significant_frequencies_polar(long[long["stem"] == "a"])
    assert len(ax.collections[0].get_offsets()) == 4

    with pytest.raises(ValueError, match="orientation"):
        plot_significant_frequencies_polar(long, orientation="ridge")
    with pytest.raises(ValueError, match="radius"):
        plot_significant_frequencies_polar(long, radius="period")
    with pytest.raises(ValueError, match="missing columns"):
        plot_significant_frequencies_polar(long.drop(columns="normalised_psd"))


# --- whole-run loader and the 2 x 2 tile summary ----------------------------

N_FREQ = 32          # padded spectrum size
DF = 1.0 / 8.0       # cycles/m per cell (8 m padded tile at 0.25 m)
PEAK_EAST, PEAK_NORTH = 5, 3   # peak offset from zero frequency, in cells


def _write_freq_raster(path, data):
    """Write as lsdtt-spectral-analysis does: DC centre at (0, 0), north-up."""
    import rasterio
    from rasterio.transform import from_origin

    half = N_FREQ // 2
    x_min = -(half + 0.5) * DF
    y_max = (half + 0.5) * DF
    with rasterio.open(path, "w", driver="GTiff", height=N_FREQ, width=N_FREQ,
                       count=1, dtype="float32", nodata=-9999.0,
                       transform=from_origin(x_min, y_max, DF, DF)) as dst:
        dst.write(data.astype("float32"), 1)


def _write_run(tmp_path, *, radial=True, rasters=True, input_tif=True):
    """A synthetic output set with one significant cell pair at a known spot."""
    import rasterio
    from rasterio.transform import from_origin

    out = tmp_path / "out"
    out.mkdir()
    stem = out / "syn"
    half = N_FREQ // 2
    if rasters:
        psd = np.ones((N_FREQ, N_FREQ))
        # Row 0 is the most NORTHWARD frequency, so north +3 is row half - 3.
        psd[half - PEAK_NORTH, half + PEAK_EAST] = 50.0
        psd[half + PEAK_NORTH, half - PEAK_EAST] = 50.0
        _write_freq_raster(f"{stem}_psd_2d.tif", psd)
        _write_freq_raster(f"{stem}_background_psd_2d.tif", np.ones((N_FREQ, N_FREQ)))
    if radial:
        f = np.array([0.125, 0.25, 0.5, 1.0])
        pd_ = __import__("pandas")
        pd_.DataFrame({"frequency_cycles_per_m": f, "wavelength_m": 1 / f,
                       "psd": [1.0, 0.5, 0.2, 0.05], "background_psd": [1.0, 0.4, 0.2, 0.05],
                       "normalised_psd": [1.0, 3.5, 1.0, 1.0],
                       "beta_by_bin": [2.0] * 4, "r_sq_by_bin": [0.9] * 4}
                      ).to_csv(f"{stem}_radial_spectrum.csv", index=False)
    (out / "syn_significant_frequencies.csv").write_text(
        "# rollover_frequency_cycles_per_m,0.25\n# rollover_wavelength_m,4\n"
        "# rollover_beta,3.1\n# sub_rollover_beta,1.2\n"
        "frequency_cycles_per_m,wavelength_m,normalised_psd\n0.25,4,3.5\n")
    fe, fn = PEAK_EAST * DF, PEAK_NORTH * DF
    f = float(np.hypot(fe, fn))
    b = float(np.degrees(np.arctan2(fe, fn)))
    (out / "syn_significant_frequencies_2d.csv").write_text(
        HEADER_ONLY_CSV_V1
        + f"{np.radians(b)},{b},{f},{1 / f},50\n{np.radians(b + 180)},{b + 180},{f},{1 / f},50\n")
    src = tmp_path / "input" / "syn.tif"
    if input_tif:
        src.parent.mkdir()
        with rasterio.open(src, "w", driver="GTiff", height=16, width=16, count=1,
                           dtype="float32",
                           transform=from_origin(300000, 800004, 0.25, 0.25)) as dst:
            dst.write(np.arange(256, dtype="float32").reshape(16, 16), 1)
    (out / "syn_ingested_params.param").write_text(
        f"# lsdtt3 parameter file - ingested parameters\nread_path: {src.parent}/\n"
        "read_fname: syn.tif\nwrite_fname: syn\n")
    return stem


def test_load_spectral_analysis_finds_every_output_and_the_input(tmp_path):
    from lsdviztools3.io.spectral import load_spectral_analysis

    run = load_spectral_analysis(_write_run(tmp_path))
    assert run.stem == "syn"
    assert len(run.radial) == 4
    assert run.rollover["rollover_frequency_cycles_per_m"] == pytest.approx(0.25)
    assert run.rollover["rollover_beta"] == pytest.approx(3.1)
    assert len(run.binned_significant) == 1
    assert run.psd_2d.shape == (N_FREQ, N_FREQ)
    assert run.significant_2d.n_cells == 2
    assert run.input_raster is not None and run.input_raster.shape == (16, 16)
    norm = run.normalised_psd_2d
    assert float(norm.max()) == pytest.approx(50.0)


def test_load_spectral_analysis_partial_run_and_missing_run(tmp_path):
    from lsdviztools3.io.spectral import load_spectral_analysis

    run = load_spectral_analysis(_write_run(tmp_path, radial=False, rasters=False,
                                            input_tif=False))
    assert run.radial is None and run.psd_2d is None and run.input_raster is None
    assert run.normalised_psd_2d is None
    assert run.significant_2d.n_cells == 2
    with pytest.raises(FileNotFoundError):
        load_spectral_analysis(tmp_path / "nothing_here")


def test_normalised_spectrum_2d_puts_a_north_east_peak_north_east(tmp_path):
    """The peak written at north +3, east +5 cells must be drawn at (+5, +3) * df.

    Checks both the image (via its extent and orientation) and the cyan
    outline, so the two cannot disagree with each other or with the bearing.
    """
    from lsdviztools3.charts import plot_normalised_spectrum_2d
    from lsdviztools3.io.spectral import load_spectral_analysis

    run = load_spectral_analysis(_write_run(tmp_path))
    fig, ax = plot_normalised_spectrum_2d(run)
    im = ax.images[0]
    arr = np.asarray(im.get_array(), dtype=float)
    left, right, bottom, top = im.get_extent()
    assert im.origin == "upper"
    dx = (right - left) / arr.shape[1]
    dy = (top - bottom) / arr.shape[0]
    candidates = np.argwhere(np.nan_to_num(arr) == np.nanmax(arr))
    drawn = {(round(left + (c + 0.5) * dx, 6), round(top - (r + 0.5) * dy, 6))
             for r, c in candidates}
    want_ne = (round(PEAK_EAST * DF, 6), round(PEAK_NORTH * DF, 6))
    assert want_ne in drawn
    assert (-want_ne[0], -want_ne[1]) in drawn

    # Cyan outline around the same two cells.
    verts = np.concatenate([p.vertices for c in ax.collections
                            for p in (c.get_paths() if hasattr(c, "get_paths") else [])])
    centres = [(round(v[0] / DF), round(v[1] / DF)) for v in verts]
    assert any(abs(x - PEAK_EAST) <= 1 and abs(y - PEAK_NORTH) <= 1 for x, y in centres)
    assert not any(abs(x - PEAK_EAST) <= 1 and abs(y + PEAK_NORTH) <= 1 for x, y in centres)


def test_radial_spectrum_counts_bins_and_marks_rollover(tmp_path):
    from lsdviztools3.charts import plot_radial_spectrum
    from lsdviztools3.io.spectral import load_spectral_analysis

    run = load_spectral_analysis(_write_run(tmp_path))
    fig, ax = plot_radial_spectrum(run)
    assert "1 bins > 2.9955" in ax.get_title()
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any(t.startswith("rollover 4") for t in labels)

    (tmp_path / "b").mkdir()
    bare = load_spectral_analysis(_write_run(tmp_path / "b", radial=False))
    with pytest.raises(ValueError, match="radial"):
        plot_radial_spectrum(bare)


def test_tile_summary_is_2x2_and_tolerates_missing_outputs(tmp_path):
    from lsdviztools3.charts import plot_spectral_tile_summary

    fig = plot_spectral_tile_summary(_write_run(tmp_path))
    panels = [a for a in fig.axes
              if a.get_subplotspec() is not None and a.get_label() != "<colorbar>"]
    grid = [a.get_subplotspec().get_geometry() for a in panels]
    # The radial panel's twin axis shares cell 2.
    assert {g[2] for g in grid if g[:2] == (2, 2)} == {0, 1, 2, 3}
    polar = [a for a in fig.axes if a.name == "polar"][0]
    nyquist = (N_FREQ // 2 + 0.5) * DF
    assert polar.get_ylim()[1] == pytest.approx(np.sqrt(2) * nyquist)
    assert fig._suptitle.get_text() == "syn"

    (tmp_path / "p").mkdir()
    partial = plot_spectral_tile_summary(
        _write_run(tmp_path / "p", radial=False, rasters=False, input_tif=False))
    notes = [t.get_text() for a in partial.axes for t in a.texts]
    assert any("input raster not found" in t for t in notes)
    assert any("no 2D spectrum rasters" in t for t in notes)
    assert any("no radial spectrum" in t for t in notes)
