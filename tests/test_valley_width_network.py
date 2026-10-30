"""Valley width network map and profile (lsdtt3 run_valley_width_network)."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import LineString

from lsdviztools3.charts.valley_width_network import (
    confluence_distances,
    plot_channel_and_centreline,
    largest_trunk,
    plot_valley_width_map,
    plot_valley_widths_along_channel,
)
from lsdviztools3.cli import plot_valley_widths

CRS = "EPSG:32632"


def _write_run(tmp_path):
    """Two trunks (0 small, 1 large) and a tributary (2) of trunk 1."""
    stem = tmp_path / "w"
    pd.DataFrame({
        "segment_id": [0, 1, 2], "parent_segment_id": [-1, -1, 1],
        "nesting_depth": [0, 0, 1], "n_channel_pixels": [500, 300, 50],
        "n_centreline_nodes": [10, 11, 5],
        "head_x": [0.0, 500.0, 700.0], "head_y": [0.0, 1000.0, 600.0],
        "end_x": [0.0, 500.0, 520.0], "end_y": [-100.0, 0.0, 600.0],
        "drains_to_edge": [1, 1, 1], "catchment_touches_edge": [0, 0, 0],
        "valley_width_estimate_m": [50.0, 100.0, 30.0],
        "medial_path_failed": [0, 0, 0], "n_samples_snapped_to_mask_edge": [0, 0, 0],
    }).to_csv(f"{stem}_valley_segments.csv", index=False)

    # Trunk 1 runs south along x = 500 from y = 1000; nodes every 100 m.
    nodes = [(1, i, 500.0, 1000.0 - 100.0 * i, float(100 * i)) for i in range(11)]
    nodes += [(2, i, 700.0 - 45.0 * i, 600.0, float(45 * i)) for i in range(5)]
    nodes += [(0, i, 0.0, -10.0 * i, float(10 * i)) for i in range(10)]
    pd.DataFrame(nodes, columns=["source_id", "node_seq", "x", "y", "distance_m"]).assign(
        elevation=0.0, trough_elevation=0.0, baseline_elevation=0.0, is_pin=0
    ).to_csv(f"{stem}_valley_centreline_nodes.csv", index=False)

    reasons1 = ["not_in_valley", "none", "none", "tributary_crossing", "none",
                "side_not_found", "none"]
    rows = []
    for k, r in enumerate(reasons1):
        rows.append((1, -1, 0, k, 50.0 + 150.0 * k, 500.0, 950.0 - 150.0 * k,
                     1e6 * (k + 1), 0.0 if r == "not_in_valley" else 100.0 + 10 * k, r))
    rows.append((0, -1, 0, 0, 5.0, 0.0, -5.0, 2e5, 40.0, "none"))
    rows.append((2, 1, 1, 0, 20.0, 680.0, 600.0, 3e5, 30.0, "none"))
    rows.append((2, 1, 1, 1, 60.0, 640.0, 600.0, 3e5, 30.0, "near_parent_channel"))
    w = pd.DataFrame(rows, columns=["segment_id", "parent_segment_id", "nesting_depth",
                                    "station_seq", "along_distance_m", "x", "y",
                                    "drainage_area_m2", "valley_width_m", "excluded_reason"])
    w["left_half_width_m"] = w["valley_width_m"] / 2
    w["right_half_width_m"] = w["valley_width_m"] / 2
    w["channel_elevation"] = 0.0
    w["drains_to_edge"] = 1
    w["catchment_touches_edge"] = 0
    w.to_csv(f"{stem}_valley_widths.csv", index=False)

    geoms, attrs = [], []
    for r in w.itertuples():
        half = r.valley_width_m / 2
        geoms.append(LineString([(r.x - half, r.y), (r.x + half, r.y)]))
        attrs.append({"source_id": r.segment_id, "station_seq": r.station_seq,
                      "segment_id": r.segment_id, "excluded_reason": r.excluded_reason,
                      "valley_width_m": r.valley_width_m,
                      "left_x": r.x - half, "left_y": r.y,
                      "right_x": r.x + half, "right_y": r.y})
    gpd.GeoDataFrame(attrs, geometry=geoms, crs=CRS).to_file(
        f"{stem}_valley_cross_sections.fgb", driver="FlatGeobuf")
    gpd.GeoDataFrame({"segment_id": [0, 1, 2], "parent_segment_id": [-1, -1, 1]},
                     geometry=[LineString([(0, 0), (5, -45), (0, -90)]),
                               LineString([(500, 1000), (530, 500), (500, 0)]),
                               LineString([(700, 600), (600, 610), (520, 600)])],
                     crs=CRS).to_file(f"{stem}_valley_segment_channels.fgb",
                                      driver="FlatGeobuf")
    gpd.GeoDataFrame({"source_id": [0, 1, 2]}, geometry=[
        LineString([(0, 0), (0, -90)]), LineString([(500, 1000), (500, 0)]),
        LineString([(700, 600), (520, 600)])], crs=CRS).to_file(
        f"{stem}_valley_centreline.fgb", driver="FlatGeobuf")

    transform = from_origin(-200.0, 1100.0, 10.0, 10.0)
    for name, data in (("hs.tif", np.full((130, 110), 180, dtype=np.float32)),
                       ("mask.tif", np.ones((130, 110), dtype=np.uint8))):
        with rasterio.open(tmp_path / name, "w", driver="GTiff", height=130, width=110,
                           count=1, dtype=data.dtype, crs=CRS,
                           transform=transform) as dst:
            dst.write(data, 1)
    return stem


def test_largest_trunk_is_the_one_reaching_the_largest_area(tmp_path):
    stem = _write_run(tmp_path)
    w = pd.read_csv(f"{stem}_valley_widths.csv")
    s = pd.read_csv(f"{stem}_valley_segments.csv")
    # Trunk 0 has more channel pixels; trunk 1 reaches the larger area.
    assert largest_trunk(w, s) == 1


def test_confluence_projects_onto_the_parent_centreline(tmp_path):
    stem = _write_run(tmp_path)
    c = confluence_distances(stem, 1)
    assert list(c["segment_id"]) == [2]
    # Tributary ends at (520, 600): nearest trunk node is (500, 600), 400 m down.
    assert c["distance_m"].iloc[0] == 400.0
    assert confluence_distances(stem, 2).empty


def test_profile_and_map_write(tmp_path):
    stem = _write_run(tmp_path)
    prof = plot_valley_widths_along_channel(stem, tmp_path / "p.png")
    assert prof.exists() and prof.stat().st_size > 0
    m = plot_valley_width_map(stem, tmp_path / "hs.tif", tmp_path / "m.png",
                              valley_mask=tmp_path / "mask.tif")
    assert m.exists() and m.stat().st_size > 0
    # A subset of segments and only the kept sections.
    m2 = plot_valley_width_map(stem, tmp_path / "hs.tif", tmp_path / "m2.png",
                               segment_ids=[2], show_excluded=False)
    assert m2.exists()


def test_channel_and_centreline_map(tmp_path):
    stem = _write_run(tmp_path)
    out = plot_channel_and_centreline(stem, tmp_path / "hs.tif", tmp_path / "c.png",
                                      valley_mask=tmp_path / "mask.tif")
    assert out.exists() and out.stat().st_size > 0
    out2 = plot_channel_and_centreline(stem, tmp_path / "hs.tif", tmp_path / "c2.png",
                                       segment_ids=[1, 2], bbox=(400, 0, 800, 1000))
    assert out2.exists()


def test_cli(tmp_path):
    stem = _write_run(tmp_path)
    rc = plot_valley_widths.main([
        "--stem", str(stem), "--base", str(tmp_path / "hs.tif"),
        "--out-map", str(tmp_path / "cm.png"), "--out-profile", str(tmp_path / "cp.png"),
        "--out-channels", str(tmp_path / "cc.png"), "--segment-id", "1"])
    assert rc == 0
    for f in ("cm.png", "cp.png", "cc.png"):
        assert (tmp_path / f).exists()
