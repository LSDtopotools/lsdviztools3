# Hillslopes and ridgelines

Hillslope process diagnostics from real landscapes (slope–area, gradient
against curvature, E*–R*, hilltop–channel coupling) and ridgeline profiles from
`lsdtt-ridgeline-analysis`. For the solver tests see
[Landscape evolution](landscape_evolution.md#hillslope-solver-comparisons).

## Hillslope charts

Each module docstring explains what the figure tests and how to read it; read
it before quoting a result.

| Function | Module | What it draws |
|---|---|---|
| `plot_slope_area` | `charts.slope_area` | The slope–area "boomerang" from slope + drainage-area rasters, with the binned median and the hillslope/channel rollover interval (`binned_median`, `rollover_area`). |
| `plot_gradient_curvature` | `charts.gradient_curvature` | Convexity vs gradient with the linear and Roering nonlinear flux-law predictions overlaid (after Roering, 2008). |
| `plot_estar_rstar` | `charts.estar_rstar` | Measured hilltop E*–R* against the Roering et al. (2007) steady relation, plus the R* residual panel; `summarise` gives per-group residual statistics. |

```bash
# E*-R* figure + per-group residual summary CSV
pixi run lsdviz3-estar-rstar --metrics hilltop_metrics.csv --outdir demo_output/estar
```

## Hilltop–channel coupling (`lsdviz3-hilltop-coupling`)

Hilltop/channel coupling (`lsdtt-hillslope-channel-coupling` output, loaded with
`io.hilltop_coupling.load_hilltop_outputs`) has its own chart set in
`charts.hilltop_coupling`: `plot_traces_on_hillshade`, `plot_trace_outcomes`,
`plot_ridgelines_and_sides`, `plot_length_distributions`,
`plot_ridgeline_profiles`, `plot_looped_trace_profiles`.

```bash
# hilltop coupling figure set from one output directory or <write_path>/<write_fname> stem
pixi run lsdviz3-hilltop-coupling --prefix ../data/run/run_stem --outdir demo_output/hilltop
```

Flags: `--prefix --outdir [--dem --channels --max-traces --dpi]`.

## Ridgeline profiles

Four commands plot ridgelines from `lsdtt-ridgeline-analysis` (its
`*_with_polyfit.csv` profile tables). Each takes its settings either as flags
or from a driver file: `--driver file.txt`, with one `key: value` per line,
each key setting the flag of the same name (`write_path` sets
`--write-path`); flags on the command line win. Each writes PNG, SVG and PDF
into `png/`, `svg/` and `pdf/` under `--write-path`, named `--write-name`.
`--ridgeline-number` picks the ridgeline (default 0), `--reverse-ridgeline`
flips its direction, and `--use-distance-range` with `--minimum-distance` /
`--maximum-distance` clips it.

| Command | What it draws | Main inputs |
|---|---|---|
| `lsdviz3-ridgeline-curvature` | Ridge-crest elevation, moving-window detrending, maximal curvature and the Laplacian along one ridgeline | `--ridgeline-profile-csv`; `--moving-window-sizes-pixels` and `--curvature-lambdas-m` (four values each, default 25,49,99,149); `--site-display-name` |
| `lsdviz3-ridgeline-summary` | Triple divide junctions, a ridgeline series with its curvature barcode, and a map of the ridgeline | `--ridgeline-profile-csv`; `--hillshade-raster-file` (map); up to three `--native-tdj-csv-N` and `--topological-tdj-csv-N` (threshold read from the file name); `--channel-network-vector-file`, `--paired-flowpaths-vector-file`; `--ridgeline-series`, `--highlighted-tdj-type`, `--highlighted-tdj-threshold` |
| `lsdviz3-ridgeline-spacing` | Maximal-curvature zero-crossing spacing against polyfit scale beside slope–area, one row per site | `--site NAME PROFILE_CSV SLOPE_AREA_CSV` (repeat per site) or `site_<N>_*` driver entries with `--site-count`; `--log-spacing`, `--spacing-y-max-m` |
| `lsdviz3-mdd-ridgeline` | A main-drainage-divide ridgeline profile above the model's drainage-area image | `--ridgeline-profile-csv` (range-defined, with a maximum-curvature metric); `--drainage-area-map-png` (`lsdtt-lem` log₁₀ drainage area of the same grid); `--model-name`, `--guide-interval` |
