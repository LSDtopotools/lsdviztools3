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

HEADER_ONLY_CSV = """\
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

    p = _write(tmp_path, "t_significant_frequencies_2d.csv", HEADER_ONLY_CSV)
    sig = load_significant_frequencies_2d(p)
    assert sig.n_cells == 0
    assert list(sig.table.columns) == COLUMNS
    assert sig.threshold == pytest.approx(9.21034)


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
