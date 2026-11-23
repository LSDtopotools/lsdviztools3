"""CLI: lsdviz3-slope-stability - the outputs of an lsdtt-slope-stability run.

Reads the run's rasters and tables (``<prefix>_<suffix>``; see
:mod:`lsdviztools3.io.slope_stability`) and writes, for each output present:

  <out-prefix>_sinmap_slope_area.png  SINMAP's slope - specific catchment area
      plot: cells coloured by stability class, with the class boundaries (one
      figure per region of the boundaries table)
  <out-prefix>_landslide_area_volume.png, _landslide_area_depth.png,
  _landslide_frequency_area.png, _landslide_length_width.png
      the rupture-plane catalogues (by daylight point and amalgamated); the
      frequency - area plot draws lsdtt3's binned p(A) and its power-law tail
      fit (eta, its standard error and A_min in the legend). Without the
      frequency - area tables (an older lsdtt3) _landslide_area_distribution.png,
      a histogram, is drawn instead
  <out-prefix>_slip_body_area_volume.png, ..._area_depth.png,
  ..._frequency_area.png, ..._length_width.png
      the slip-body catalogues (by toe and amalgamated)
  <out-prefix>_depth_comparison.png  per crown: the slip surface's deepest body
      thickness against the rupture-plane depth
  with --dem, maps on its hillshade:
  <out-prefix>_sinmap_class_map.png, _rupture_depth_map.png,
  _rupture_failure_probability_map.png, _slip_safety_factor_map.png (the
  screened-out crowns grey), _slip_failure_probability_map.png,
  _slip_body_depth_map.png, _slip_cover_safety_factor_map.png,
  _slip_cover_failure_probability_map.png, _water_table_depth_map.png
  from a DEM of difference (run_dem_change_landslides):
  <out-prefix>_dem_change_uncertainty.png  error of the difference: its
      distribution, NMAD against slope, and the spatial variogram
  <out-prefix>_dem_change_coverage.png  area in the study domain, and how it
      depends on the detection threshold
  <out-prefix>_dem_change_frequency_area.png  p(A) with the power-law tails
  <out-prefix>_dem_change_false_detection.png  patches found against those from
      noise alone, by size
  with --dem (the newer DEM): <out-prefix>_dem_change_map.png, the changed
  pixels, and _dem_difference_map.png, the bias-corrected difference

Everything drawn is read from the files; nothing is recomputed.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-slope-stability``."""
    p = argparse.ArgumentParser(
        prog="lsdviz3-slope-stability",
        description="Plots of an lsdtt-slope-stability run: SINMAP's slope-area plot and class map, "
                    "the landslide and slip-body catalogues (area-volume, area-depth, area and "
                    "length/width distributions), rupture against slip-surface depth, maps of "
                    "depth, factor of safety and probability of failure, and potential landslides "
                    "from a DEM of difference with their uncertainty.")
    p.add_argument("--data-dir", required=True, help="Directory holding the lsdtt-slope-stability outputs.")
    p.add_argument("--prefix", required=True, help="The run's write_fname.")
    p.add_argument("--out-dir", required=True, help="Directory for the outputs.")
    p.add_argument("--out-prefix", default=None, help="Prefix of the output names (default --prefix).")
    p.add_argument("--dem", default=None, help="DEM for the maps' hillshade (the run's input DEM). "
                                               "Without it no maps are drawn.")
    p.add_argument("--safety-factor-range", type=float, nargs=2, default=(0.5, 1.5), metavar=("LOW", "HIGH"),
                   help="Colour range of the factor-of-safety maps, centred on 1 by default (0.5 1.5).")
    p.add_argument("--max-points", type=int, default=200_000,
                   help="At most this many cells in the SINMAP and depth scatter plots (every k-th "
                        "cell; default 200000).")
    p.add_argument("--dem-difference-range", type=float, default=None, metavar="METRES",
                   help="Colour range of the DEM of difference map, +/- this (default 5 x NMAD).")
    p.add_argument("--coordinates", default="local", help="Map coordinates: local/projected or geographic.")
    p.add_argument("--dpi", type=int, default=300, help="DPI of the charts (default 300).")
    return p


def _save(fig, path, dpi, written):
    import matplotlib.pyplot as plt

    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    written.append(path)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    data_dir, out_dir = Path(args.data_dir), Path(args.out_dir)
    stem = args.out_prefix or args.prefix
    if not data_dir.is_dir():
        raise SystemExit(f"{data_dir} is not a directory.")

    from ..io.slope_stability import (load_catalogue, load_frequency_area, load_sinmap_class_areas,
                                      load_sinmap_curves, present, read_band, slope_stability_files)

    f = present(slope_stability_files(data_dir, args.prefix))
    if not f:
        raise SystemExit(f"No lsdtt-slope-stability outputs for prefix {args.prefix!r} in {data_dir}.")
    out_dir.mkdir(parents=True, exist_ok=True)

    from ..charts import slope_stability as ch

    written = []

    def note(msg):
        print(msg)

    # SINMAP slope-area plot.
    need = ("sinmap_class", "sinmap_curves", "dinf_slope", "specific_catchment_area")
    if all(k in f for k in need):
        curves = load_sinmap_curves(f["sinmap_curves"])
        slope, sca, cls = (read_band(f[k]) for k in ("dinf_slope", "specific_catchment_area", "sinmap_class"))
        regions = sorted(int(r) for r in curves["region_id"].unique())
        for r in regions:
            fig, _ = ch.plot_sinmap_slope_area(slope, sca, cls, curves, region_id=r, max_points=args.max_points,
                                               title=None if len(regions) == 1 else
                                               f"SINMAP slope - a: region {r} boundaries over all cells")
            name = f"{stem}_sinmap_slope_area.png" if len(regions) == 1 else f"{stem}_sinmap_slope_area_region_{r}.png"
            _save(fig, out_dir / name, args.dpi, written)
    else:
        note("no SINMAP slope-area plot: it needs " + ", ".join(
            f"{args.prefix}{'_' + k if k != 'sinmap_curves' else '_sinmap_slope_area_curves.csv'}" for k in need if k not in f)
             + " (run_sinmap, write_sinmap_rasters and write_sinmap_slope_area_curves_table)")

    # Catalogues.
    def catalogue_charts(keys, labels, out_word, title_word, fa_keys):
        cats = [(lab, load_catalogue(f[k])) for k, lab in zip(keys, labels) if k in f]
        cats = [(lab, df) for lab, df in cats if not df.empty]
        if not cats:
            note(f"no {title_word} catalogue charts: no non-empty {'/'.join(keys)} tables")
            return
        charts = [("area_volume", ch.plot_catalogue_area_volume, "area - volume"),
                  ("area_depth", ch.plot_catalogue_area_depth, "area - mean depth"),
                  ("length_width", ch.plot_catalogue_length_width, "length / width")]
        if not all(k in f for k in fa_keys):
            note(f"no {title_word} frequency - area plot ({args.prefix}_{out_word}_frequency_area*.csv not found, "
                 "an older lsdtt3); drawing the area histogram instead")
            charts.append(("area_distribution", ch.plot_catalogue_area_distribution, "area distribution"))
        for suffix, fn, word in charts:
            fig, _ = fn(cats, title=f"{title_word}: {word}")
            _save(fig, out_dir / f"{stem}_{out_word}_{suffix}.png", args.dpi, written)
        if all(k in f for k in fa_keys):
            density, fit = load_frequency_area(f[fa_keys[0]], f[fa_keys[1]])
            fig, _ = ch.plot_frequency_area(density, fit, title=f"{title_word}: frequency - area")
            _save(fig, out_dir / f"{stem}_{out_word}_frequency_area.png", args.dpi, written)

    catalogue_charts(("landslides", "landslides_amalgamated"), ("by daylight point", "amalgamated"),
                     "landslide", "rupture-plane landslides",
                     ("landslide_frequency_area", "landslide_frequency_area_fit"))
    catalogue_charts(("slip_bodies", "slip_bodies_amalgamated"), ("by toe", "amalgamated"),
                     "slip_body", "slip bodies",
                     ("slip_body_frequency_area", "slip_body_frequency_area_fit"))

    # Depth comparison.
    if "rupture_depth" in f and "slip_max_depth" in f:
        try:
            fig, _ = ch.plot_depth_comparison(read_band(f["rupture_depth"]), read_band(f["slip_max_depth"]),
                                              max_points=args.max_points)
            _save(fig, out_dir / f"{stem}_depth_comparison.png", args.dpi, written)
        except ValueError as e:
            note(f"no depth comparison: {e}")

    # DEM of difference.
    from ..charts import dem_change as dc
    written += dc.write_dem_change_charts(f, out_dir, stem, dpi=args.dpi, note=note)

    # Maps.
    if not args.dem:
        note("no maps: give --dem for the hillshade")
    else:
        dem = args.dem
        lo, hi = args.safety_factor_range
        maps = []
        if "sinmap_class" in f:
            areas = load_sinmap_class_areas(f["sinmap_class_areas"]) if "sinmap_class_areas" in f else None
            written.append(ch.render_sinmap_class_map(dem, f["sinmap_class"], out_dir / f"{stem}_sinmap_class_map.png",
                                                      class_areas=areas, coordinates=args.coordinates))
        maps += [
            ("rupture_depth", "_rupture_depth_map.png", ch.DEPTH_CMAP, "rupture depth h (m)", 0.0, None,
             "rupture planes: depth", None),
            ("rupture_failure_probability", "_rupture_failure_probability_map.png", ch.PROBABILITY_CMAP,
             "probability of failure", 0.0, 1.0, "rupture planes: P(F <= 1)", None),
            ("slip_safety_factor", "_slip_safety_factor_map.png", ch.SAFETY_FACTOR_CMAP, "factor of safety", lo, hi,
             "slip surfaces: factor of safety", "slip_evaluated"),
            ("slip_failure_probability", "_slip_failure_probability_map.png", ch.PROBABILITY_CMAP,
             "probability of failure", 0.0, 1.0, "slip surfaces: P(F <= 1), fixed surface", "slip_evaluated"),
            ("slip_body_depth", "_slip_body_depth_map.png", ch.DEPTH_CMAP, "body depth (m)", 0.0, None,
             "unstable slip bodies: depth", None),
            ("slip_cover_safety_factor", "_slip_cover_safety_factor_map.png", ch.SAFETY_FACTOR_CMAP,
             "factor of safety", lo, hi, "ellipsoids: lowest F covering each cell", None),
            ("slip_cover_failure_probability", "_slip_cover_failure_probability_map.png", ch.PROBABILITY_CMAP,
             "probability of failure", 0.0, 1.0, "ellipsoids: highest P(F <= 1) covering each cell", None),
            ("water_table_depth", "_water_table_depth_map.png", "Blues_r", "water-table depth (m)", 0.0, None,
             "groundwater: water-table depth", None),
        ]
        for key, suffix, cmap, label, vmin, vmax, title, overlay in maps:
            if key not in f:
                continue
            mask = None
            if overlay and overlay in f:
                mask = ch.screened_mask(f[overlay])
            written.append(ch.render_value_map(dem, f[key], out_dir / f"{stem}{suffix}", cmap=cmap, label=label,
                                               vmin=vmin, vmax=vmax, title=title, overlay_mask=mask,
                                               overlay_label="screened out, not evaluated",
                                               coordinates=args.coordinates))
    if args.dem:
        written += dc.write_dem_change_maps(f, args.dem, out_dir, stem, coordinates=args.coordinates,
                                            value_range=args.dem_difference_range)
    for w in written:
        print(w)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
