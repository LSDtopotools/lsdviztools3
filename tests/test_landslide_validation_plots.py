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
    pd.DataFrame({"is_zero": [True, False, False], "ah_low": [0., 1e-30, 10.], "ah_high": [0., 10., 1000.],
                  "ah_centre": [0., 1e-20, 100.], "n_cells": [20, 40, 40], "fraction_of_domain": [.2, .4, .4],
                  "ratio": [.2, .9, 2.], "ip_low": [.8, .9, .95], "ip_high": [1.2, 1.1, 1.05],
                  "significant": [True, False, True]}).to_csv(d / "validation_hazard_conditional.csv", index=False)
    pd.DataFrame({"catalogue": ["Inventory: polygons", "scaLr: amalgamated", "DoD: loss", "Hazard area >= 60"],
                  "n_patches": [100, 90, 80, 70], "area_min_m2": [5e4] * 4, "eta": [-2.3, -2.2, -2.5, -1.9],
                  "eta_low": [-2.5, -2.4, -2.8, -2.1], "eta_high": [-2.1, -2.0, -2.2, -1.7],
                  "n_tail": [40, 30, 20, 25]}).to_csv(d / "validation_frequency_area.csv", index=False)
    pd.DataFrame([{"catalogue": n, "area_centre_m2": a, "probability_density_per_m2": 1 / a ** 2, "n_landslides": 5}
                  for n in ("Inventory: polygons", "scaLr: amalgamated", "DoD: loss", "Hazard area >= 60")
                  for a in (1e4, 1e5)]).to_csv(d / "validation_frequency_area_density.csv", index=False)


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
    for name in ("skill", "frequency_ratio", "detectors", "dem_change_by_size", "hazard_sweep",
                 "hazard_conditional", "frequency_area"):
        assert (out / f"t_{name}.png").stat().st_size > 0


def test_hazard_conditional_needs_nonzero_bins():
    table = pd.DataFrame({"is_zero": [True], "ah_low": [0.], "ah_high": [0.], "ah_centre": [0.], "n_cells": [5],
                          "fraction_of_domain": [1.], "ratio": [1.], "ip_low": [.9], "ip_high": [1.1],
                          "significant": [False]})
    with pytest.raises(ValueError, match="no bins"):
        ch.plot_hazard_conditional(table)


def test_frequency_area_rejects_an_unknown_common_area(tmp_path):
    _tables(tmp_path)
    t = ch.load_validation_tables(tmp_path)
    with pytest.raises(ValueError, match="no fits"):
        ch.plot_frequency_area_comparison(t["frequency_fit"], t["frequency_density"], area_min_m2=123.0)


def test_frequency_area_curves_must_exist(tmp_path):
    _tables(tmp_path)
    t = ch.load_validation_tables(tmp_path)
    with pytest.raises(ValueError, match="not in the table"):
        ch.plot_frequency_area_comparison(t["frequency_fit"], t["frequency_density"], curves=["missing"])
