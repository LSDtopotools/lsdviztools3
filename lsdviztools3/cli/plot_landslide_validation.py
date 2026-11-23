"""CLI: lsdviz3-landslide-validation - slope-stability outputs against an inventory.

Reads tables of scores (``validation_*.csv``) from ``--tables-dir`` and writes:

  <prefix>_skill.png             ROC curves and AUC of each predictor
  <prefix>_frequency_ratio.png   landslide probability along each predictor
  <prefix>_detectors.png         precision, recall and lift of each detector
  <prefix>_dem_change_by_size.png  DoD against inventory landslide and patch size
  <prefix>_hazard_sweep.png      hazard-area AUC over its two thresholds

With ``--dem`` and ``--inventory-raster`` (1 on landslide cells), and any of
``--rupture-unstable``, ``--hazard-area`` and ``--dem-change-class``:

  <prefix>_agreement_map.png     hits, false alarms and missed cells per method

Scores are read from the tables. Only the map categories are derived from rasters.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-landslide-validation``."""
    p = argparse.ArgumentParser(prog="lsdviz3-landslide-validation",
                                description="Figures comparing slope-stability outputs with a landslide inventory.")
    p.add_argument("--tables-dir", required=True, help="Directory holding validation_*.csv.")
    p.add_argument("--out-dir", required=True, help="Directory for the figures.")
    p.add_argument("--prefix", default="validation", help="Prefix of the output names (default validation).")
    p.add_argument("--base-rate", type=float, default=None,
                   help="Share of the domain in the inventory, marked on the precision axis.")
    p.add_argument("--dem", default=None, help="DEM for the map hillshade.")
    p.add_argument("--inventory-raster", default=None, help="Raster with 1 on inventory cells, on the DEM's grid.")
    p.add_argument("--rupture-unstable", default=None, help="<run>_rupture_unstable.tif.")
    p.add_argument("--hazard-area", default=None, help="<run>_hazard_area.tif.")
    p.add_argument("--hazard-threshold", type=float, default=60.0,
                   help="Hazard area (m2/m) drawn as a detection on the map (default 60).")
    p.add_argument("--dem-change-class", default=None, help="<run>_dem_change_class.tif.")
    p.add_argument("--dpi", type=int, default=300, help="DPI of the figures (default 300).")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from ..charts import landslide_validation as ch

    tables_dir, out_dir = Path(args.tables_dir), Path(args.out_dir)
    if not tables_dir.is_dir():
        raise SystemExit(f"{tables_dir} is not a directory.")
    t = ch.load_validation_tables(tables_dir)
    if not t:
        raise SystemExit(f"No validation_*.csv tables in {tables_dir}.")
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    def save(fig, name):
        path = out_dir / f"{args.prefix}_{name}.png"
        fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(path)

    if {"auc", "roc"} <= set(t):
        save(ch.plot_skill(t["auc"], t["roc"])[0], "skill")
    if "frequency_ratio" in t:
        save(ch.plot_frequency_ratio(t["frequency_ratio"])[0], "frequency_ratio")
    if "binary" in t:
        save(ch.plot_detector_scores(t["binary"], base_rate=args.base_rate)[0], "detectors")
    if {"polygon_recall", "patch_precision"} <= set(t):
        save(ch.plot_dem_change_by_size(t["polygon_recall"], t["patch_precision"])[0], "dem_change_by_size")
    if "sweep" in t:
        save(ch.plot_hazard_sweep(t["sweep"])[0], "hazard_sweep")

    if args.dem and args.inventory_raster:
        from ..io.slope_stability import read_band

        dem = read_band(args.dem)
        inventory = np.nan_to_num(read_band(args.inventory_raster)) == 1
        detections = {}
        scores = {}
        binary = t.get("binary")

        def score(key, name):
            if binary is not None and (binary.detector == key).any():
                row = binary[binary.detector == key].iloc[0]
                scores[name] = (float(row.precision), float(row.recall))

        if args.rupture_unstable:
            detections["scaLr: F ≤ 1"] = np.nan_to_num(read_band(args.rupture_unstable)) == 1
            score("scaLr: F <= 1", "scaLr: F ≤ 1")
        if args.hazard_area:
            name = f"Hazard area ≥ {args.hazard_threshold:g} m²/m"
            detections[name] = np.nan_to_num(read_band(args.hazard_area)) >= args.hazard_threshold
            score(f"Hazard area >= {args.hazard_threshold:g} m2/m", name)
        if args.dem_change_class:
            detections["DoD: loss + gain patches"] = np.isin(np.nan_to_num(read_band(args.dem_change_class)), [-1, 1])
            score("DoD: loss + gain", "DoD: loss + gain patches")
        if detections:
            with __import__("rasterio").open(args.dem) as d:
                cell = abs(d.transform.a)
            save(ch.plot_agreement_maps(dem, inventory, detections, scores, cell_size=cell)[0], "agreement_map")

    for path in written:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
