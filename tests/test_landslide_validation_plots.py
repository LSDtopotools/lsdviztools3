"""lsdviz3-landslide-validation: table loading, map categories and the CLI."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest

from lsdviztools3.charts import landslide_validation as ch
from lsdviztools3.cli.plot_landslide_validation import main


def _tables(d):
    names = ["Hazard area (a)", "Local slope", "Upslope area (MFD)", "scaLr: lower F", "DoD: |change|"]
    pd.DataFrame({"predictor": names, "auc": [.7, .65, .55, .66, .6], "ci95_low": [.68, .63, .53, .64, .58],
                  "ci95_high": [.72, .67, .57, .68, .62]}).to_csv(d / "validation_auc.csv", index=False)
    pd.DataFrame([{"predictor": n, "false_positive_rate": x, "true_positive_rate": x ** .7} for n in names
                  for x in (0, .5, 1)]).to_csv(d / "validation_roc.csv", index=False)
    pd.DataFrame([{"predictor": n, "x_median": x, "frequency_ratio": x + .5, "random_low": .95, "random_high": 1.05}
                  for n in ("Hazard area", "Local slope") for x in (0., .5, 1.)]).to_csv(
        d / "validation_frequency_ratio.csv", index=False)
    pd.DataFrame({"detector": ["scaLr: F <= 1", "DoD: loss + gain"], "precision": [.2, .3], "recall": [.5, .1],
                  "lift": [1.6, 2.4], "slope_matched_lift": [1.1, 1.8], "fraction_of_domain": [.3, .06]}).to_csv(
        d / "validation_binary.csv", index=False)
    pd.DataFrame({"area_low_m2": [0, 1e4], "area_high_m2": [1e4, 1e9], "n_polygons": [10, 3],
                  "any_cell_detected": [.3, .9], "half_detected": [.1, .2]}).to_csv(
        d / "validation_dod_polygon_recall.csv", index=False)
    pd.DataFrame({"change": ["loss", "gain"], "area_low_m2": [0, 0], "area_high_m2": [1e9, 1e9], "n_patches": [4, 5],
                  "majority_inside": [.3, .2]}).to_csv(d / "validation_dod_patch_precision.csv", index=False)
    pd.DataFrame({"initiation_slope": [.7, .7, .839, .839], "stopping_slope": [.1, .176, .1, .176],
                  "auc": [.6, .61, .62, .63]}).to_csv(d / "validation_hazard_sweep.csv", index=False)


def test_agreement_classes():
    p = np.array([[1, 1, 0, 0, 1]], bool)
    m = np.array([[1, 0, 1, 0, 1]], bool)
    v = np.array([[1, 1, 1, 1, 0]], bool)
    out = ch.agreement_classes(p, m, v)
    assert out[0, :4].tolist() == [1, 2, 3, 0]
    assert np.isnan(out[0, 4])


def test_missing_column_is_refused(tmp_path):
    pd.DataFrame({"predictor": ["a"], "auc": [.5]}).to_csv(tmp_path / "validation_auc.csv", index=False)
    with pytest.raises(ValueError, match="missing columns"):
        ch.load_validation_tables(tmp_path)


def test_cli_writes_every_chart(tmp_path):
    _tables(tmp_path)
    out = tmp_path / "fig"
    assert main(["--tables-dir", str(tmp_path), "--out-dir", str(out), "--prefix", "t", "--dpi", "50"]) == 0
    for name in ("skill", "frequency_ratio", "detectors", "dem_change_by_size", "hazard_sweep"):
        assert (out / f"t_{name}.png").stat().st_size > 0
