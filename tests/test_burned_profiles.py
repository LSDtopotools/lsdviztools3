"""Burned-channel long profiles (lsdtt3 burn_to_corridor _burned_profiles.csv)."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from lsdviztools3.charts.burned_profiles import (
    PIT_COLOUR,
    course_rows,
    draw_burned_profile,
    plot_burned_profiles,
    read_burned_profiles,
    select_course,
)
from lsdviztools3.cli import plot_burned_profiles as cli


def _write(tmp_path, chi=True, stem="r"):
    """Path 0 (rank 0): 10 cells east along y = 0, x = 0..90, falling 1 m per
    cell, a two-cell pit (pit_id 7) at cells 2-3. Path 1 (rank 1): 5 cells down
    from the north into path 0's cell 4 (x = 40), its last cell being that
    cell. The name "Brook" covers path 1 and path 0 downstream of the join."""
    rows = []
    for i in range(10):
        z = 100.0 - i
        pit = i in (2, 3)
        rows.append(dict(path_id=0, path_rank=0, reach_id=0, name="Brook" if i >= 4 else "Main",
                         point_index=i, distance_m=10.0 * i, flow_distance_m=10.0 * (9 - i),
                         x=10.0 * i, y=0.0, z_original=z - (0.5 if pit else 0.0), z_cleaned=z,
                         z_profile=z, z_final=z, pit_filled=int(pit), pit_id=7 if pit else -1,
                         is_record_low=int(not pit), is_rejected_low=0,
                         drainage_area_m2=1000.0 * (i + 1) if chi else "", chi=float(i) if chi else "",
                         joins_path_id=-1, join_distance_m="", starts_on_path_id=-1))
    for i in range(5):
        z = 110.0 - 2.5 * i if i < 4 else 96.0
        rows.append(dict(path_id=1, path_rank=1, reach_id=1, name="Brook", point_index=i,
                         distance_m=10.0 * i, flow_distance_m=50.0 + 10.0 * (4 - i),
                         x=40.0, y=10.0 * (4 - i), z_original=z, z_cleaned=z, z_profile=z, z_final=z,
                         pit_filled=0, pit_id=-1, is_record_low=1, is_rejected_low=0,
                         drainage_area_m2=500.0 * (i + 1) if chi else "", chi=0.5 * i if chi else "",
                         joins_path_id=0, join_distance_m=40.0, starts_on_path_id=-1))
    p = tmp_path / f"{stem}_burned_profiles.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def test_course_follows_the_join(tmp_path):
    df = read_burned_profiles(_write(tmp_path))
    c = course_rows(df, 1)
    # path 1 without its last cell (the join cell), then path 0 from the join
    assert list(c.path_id) == [1, 1, 1, 1] + [0] * 6
    assert list(c.point_index[c.path_id == 0]) == [4, 5, 6, 7, 8, 9]
    assert np.all(np.diff(c.course_distance_m) > 0)


def test_course_chi_is_continuous_and_monotone(tmp_path):
    df = read_burned_profiles(_write(tmp_path))
    c = course_rows(df, 1)
    chi = c.course_chi.to_numpy()
    assert np.all(np.diff(chi) > 0)
    # the join cell takes path 1's chi there (0.5 * 4 = 2.0); path 0 continues from it
    assert chi[4] == pytest.approx(2.0)
    assert chi[5] == pytest.approx(3.0)


def test_select_by_name_starts_at_the_most_upstream_named_cell(tmp_path):
    df = read_burned_profiles(_write(tmp_path))
    c, label = select_course(df, name="Brook")
    assert label == "Brook"
    assert c.path_id.iloc[0] == 1 and c.point_index.iloc[0] == 0
    c0, _ = select_course(df)          # default: rank 0
    assert set(c0.path_id) == {0}


def test_one_patch_per_pit(tmp_path):
    df = read_burned_profiles(_write(tmp_path))
    c, _ = select_course(df, path_rank=0)
    fig, ax = plt.subplots()
    draw_burned_profile(ax, df, c)
    pit_rgba = mcolors.to_rgba(PIT_COLOUR)
    n = 0
    for col in ax.collections:
        fc = col.get_facecolor()
        if len(fc) and np.allclose(fc[0][:3], pit_rgba[:3]):
            n += len(col.get_paths())
    plt.close(fig)
    assert n == 1


def test_compare_gives_one_panel_per_run(tmp_path):
    a, b = _write(tmp_path, stem="a"), _write(tmp_path, stem="b")
    out = tmp_path / "cmp.png"
    plot_burned_profiles([a, b], out, labels=["A", "B"], name="Brook", mark=(40.0, 0.0, "join"))
    assert out.exists() and out.stat().st_size > 10_000


def test_chi_axis_needs_a_chi_run(tmp_path):
    df = read_burned_profiles(_write(tmp_path, chi=False))
    c, _ = select_course(df)
    fig, ax = plt.subplots()
    with pytest.raises(ValueError, match="chi"):
        draw_burned_profile(ax, df, c, x="chi")
    plt.close(fig)


def test_old_format_is_refused(tmp_path):
    p = tmp_path / "old_burned_profiles.csv"
    pd.DataFrame({"reach_id": [0], "z_original": [1.0], "z_filled": [1.0], "z_burned": [1.0]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="not a burn_to_corridor"):
        read_burned_profiles(p)


def test_cli(tmp_path):
    a = _write(tmp_path, stem="a")
    out = tmp_path / "cli.png"
    assert cli.main(["--csv", str(a), "--name", "Brook", "--x", "chi", "--out", str(out)]) == 0
    assert out.exists()
    with pytest.raises(SystemExit):
        cli.main(["--csv", str(a), "--csv", str(a), "--path-id", "0", "--out", str(out)])
