"""Floodplain polygons map and HAND QQ chart (lsdtt3 fluvial-surfaces outputs)."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

from lsdviztools3.charts.floodplain_qq import (
    plot_floodplain_hand_qq,
    read_floodplain_thresholds,
)


def _write_qq(path, n=200, upper=False):
    z = np.linspace(-3, 3, n)
    ref = 5.0 + 2.0 * z
    dq = ref.copy()
    if upper:
        dq[z > 1.0] += 10.0 * (z[z > 1.0] - 1.0) ** 2   # departs upward
    else:
        dq[z < -1.0] = 0.0                               # flat low tail
    pd.DataFrame({"standard_normal_variate": z, "data_quantile": dq,
                  "normal_reference": ref}).to_csv(path, index=False)


def _write_thresholds(tmp_path):
    _write_qq(tmp_path / "x_floodplain_hand_qq.csv")
    _write_qq(tmp_path / "x_floodplain_hand_qq_second_pass.csv", upper=True)
    t = tmp_path / "x_floodplain_thresholds.csv"
    pd.DataFrame({
        "pass": [1, 1, 2],
        "quantity": ["hand", "gradient", "hand"],
        "algorithm": ["qq_gaussian", "qq_gaussian", "qq_upper_departure"],
        "threshold": [3.0, 0.05, 9.0],
        "qq_index": [40, 40, 130],
        "n_values": [1000, 1000, 400],
        "table_fname": ["x_floodplain_hand_qq.csv", "x_floodplain_gradient_qq.csv",
                        "x_floodplain_hand_qq_second_pass.csv"],
    }).to_csv(t, index=False)
    return t


def test_read_thresholds_resolves_table_paths(tmp_path):
    t = read_floodplain_thresholds(_write_thresholds(tmp_path))
    assert list(t["pass"]) == [1, 1, 2]
    assert t.loc[2, "table_path"] == tmp_path / "x_floodplain_hand_qq_second_pass.csv"


def test_hand_qq_one_panel_per_hand_pass(tmp_path):
    out = tmp_path / "qq.png"
    fig = plot_floodplain_hand_qq(_write_thresholds(tmp_path), str(out))
    assert len(fig.axes) == 2                 # the gradient row is not drawn
    assert "pass 2" in fig.axes[1].get_title()
    assert out.exists() and out.stat().st_size > 0


def test_hand_qq_without_hand_rows_raises(tmp_path):
    t = tmp_path / "t.csv"
    pd.DataFrame({"pass": [1], "quantity": ["gradient"], "algorithm": ["qq_gaussian"],
                  "threshold": [0.1], "qq_index": [3], "n_values": [10],
                  "table_fname": ["g.csv"]}).to_csv(t, index=False)
    with pytest.raises(ValueError, match="No QQ-chosen HAND"):
        plot_floodplain_hand_qq(t)


def test_render_floodplain_polygons_writes_png(tmp_path):
    pytest.importorskip("pygmt")
    from lsdviztools3.render.floodplain import render_floodplain_polygons

    crs = "EPSG:32632"
    hs = tmp_path / "hs.tif"
    data = (np.indices((60, 60)).sum(axis=0) % 255).astype("float32")
    with rasterio.open(hs, "w", driver="GTiff", height=60, width=60, count=1,
                       dtype="float32", crs=crs,
                       transform=from_origin(440000, 5350600, 10, 10)) as ds:
        ds.write(data, 1)
    big = gpd.GeoDataFrame({"label_id": [1]}, geometry=[box(440100, 5350100, 440400, 5350400)], crs=crs)
    small = gpd.GeoDataFrame({"label_id": [1]}, geometry=[box(440150, 5350150, 440300, 5350300)], crs=crs)
    fig = render_floodplain_polygons(hs, [big, small], labels=["first pass", "second pass"])
    out = tmp_path / "map.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0

    with pytest.raises(ValueError, match="labels"):
        render_floodplain_polygons(hs, [big, small], labels=["only one"])
