# Slope stability

## Potential landslides (`lsdviz3-slope-stability`)

Draws the outputs of one
[`lsdtt-slope-stability`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_slope_stability.html)
run. That CLI is still experimental. Everything is read from the run's files;
nothing is recomputed.

```bash
pixi run lsdviz3-slope-stability --data-dir runs/site --prefix site \
  --dem runs/site/site_DEM.tif --out-dir figures
```

`--prefix` is the run's `write_fname`. Each figure is drawn when the outputs
it needs are there, and the command says which ones it skipped and why.

**Charts:**

- `<prefix>_sinmap_slope_area.png`: SINMAP's slope against specific catchment
  area plot, the cells coloured by stability class, with the class boundaries
  from `_sinmap_slope_area_curves.csv`.
  - With a region table, there is one figure per region of boundaries (all
    cells are drawn in each).
  - At most `--max-points` cells are drawn (every k-th cell).
- From the rupture-plane landslide catalogues (by daylight point, and
  amalgamated):
  - `<prefix>_landslide_area_volume.png`
  - `<prefix>_landslide_area_depth.png` (mean depth)
  - `<prefix>_landslide_frequency_area.png`: the probability density p(A) in
    logarithmic bins, with the power-law tail lsdtt3 fitted, p(A) ~ A^η. The
    line starts at A_min (marked by a dotted line), and the legend gives η,
    its standard error and the number of landslides in the tail. With output
    from an older lsdtt3 that has no frequency – area tables, a plain area
    histogram, `_landslide_area_distribution.png`, is drawn instead.
  - `<prefix>_landslide_length_width.png`
- The same four for the slip-body catalogues (by toe, and amalgamated), named
  `<prefix>_slip_body_...`.
- `<prefix>_depth_comparison.png`: per crown, the slip surface's deepest body
  thickness against the rupture-plane depth, with the 1:1 line.

**Maps, with `--dem`** (the run's input DEM, for the hillshade):

- `<prefix>_sinmap_class_map.png`: the stability classes, with each class's
  share of the area.
- `<prefix>_rupture_depth_map.png` and
  `<prefix>_rupture_failure_probability_map.png`.
- `<prefix>_slip_safety_factor_map.png` and
  `<prefix>_slip_failure_probability_map.png`. The crowns screened out
  before the slip-surface search are shown in grey, since they were not
  evaluated and are not known to be stable.
- `<prefix>_slip_body_depth_map.png`: the deepest unstable body over each
  cell.
- For the 3-D ellipsoids, `<prefix>_slip_cover_safety_factor_map.png` and
  `<prefix>_slip_cover_failure_probability_map.png`.
- With a water table, `<prefix>_water_table_depth_map.png`.

**Colours:**

- The SINMAP classes run from blue (stable) to red (defended).
- Factors of safety use a diverging scale centred on 1, from red (unstable)
  to blue. `--safety-factor-range` sets its ends, 0.5 and 1.5 by default.
- Probabilities and depths use sequential scales.

The catalogues are statistics of *potential* landslides, not maps of
individual ones.

## Potential landslides from a DEM of difference

The same command draws the outputs of `run_dem_change_landslides`. Give the newer
DEM as `--dem`.

```bash
pixi run lsdviz3-slope-stability --data-dir runs/site --prefix site \
  --dem runs/site/newer_dem.tif --out-dir figures
```

- `<prefix>_dem_change_map.png`: the changed pixels on a hillshade of the newer DEM.
- `<prefix>_dem_difference_map.png`: the bias-corrected difference. Its colour range
  is five times the NMAD, or `--dem-difference-range`.
- `<prefix>_dem_change_uncertainty.png`: the distribution of the difference, the
  error against slope, and the spatial variogram.
- `<prefix>_dem_change_coverage.png`: the share of the study domain in loss and gain,
  the area expected from noise alone, and the sensitivity to the threshold.
- `<prefix>_dem_change_frequency_area.png`: p(A) with the power-law tails. Sizes where
  noise gives more than 10 % of the patches are shaded, and a fit that starts inside
  the shading is not reliable.
- `<prefix>_dem_change_false_detection.png`: patches found against those from noise
  alone, by size.

The patches are candidates. Vegetation, snow, water, radar penetration and DEM
artefacts also change elevation between two products.

## Comparing outputs with a landslide inventory

`lsdviz3-landslide-validation` draws tables of scores against a mapped inventory.
It reads `validation_*.csv` from `--tables-dir` and draws whichever are present.

```bash
pixi run lsdviz3-landslide-validation --tables-dir tables --out-dir figures --prefix site \
  --dem pre_event_dem.tif --inventory-raster inventory.tif \
  --rupture-unstable run_rupture_unstable.tif --hazard-area run_hazard_area.tif \
  --dem-change-class run_dem_change_class.tif
```

- `<prefix>_skill.png`: ROC curves, and each predictor's AUC with its interval.
- `<prefix>_frequency_ratio.png`: landslide probability along each predictor,
  relative to the study-area average.
- `<prefix>_detectors.png`: precision, recall and lift, with and without matching
  the slope distribution.
- `<prefix>_dem_change_by_size.png`: inventory landslides found and patches matched,
  by size.
- `<prefix>_hazard_sweep.png`: hazard-area AUC over its two thresholds.
- `<prefix>_hazard_conditional.png`: P(L | hazard area) / P(L) after Milledge et al.
  (2019), with a separate point for zero hazard area and the share of area per bin.
- `<prefix>_frequency_area.png`: frequency-area curves and each catalogue's tail
  exponent above one common area (`--frequency-area-min`). `--frequency-area-curves`
  chooses the catalogues on the density panel.
- `<prefix>_agreement_map.png`: cells in both, detected only, and mapped only.
