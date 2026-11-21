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
  - `<prefix>_landslide_area_distribution.png` (counts in logarithmic bins)
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
