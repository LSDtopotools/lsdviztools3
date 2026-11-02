# Landscape evolution

Animations of `lsdtt-lem` runs, and the figures used to test the hillslope
solvers against analytical and benchmark solutions.

## Animating a transient run (`lsdviz3-animate`)

```bash
pixi run lsdviz3-animate --in runs/lem --out figures/lem.gif --fps 8
```

Turns the `_model_iteration_<N>` rasters of an `lsdtt-lem` run (written every
`frequency_of_output_writing` steps) into a GIF, annotated with model time from
`_model_elevation_stats.csv`. Flags: `--in --out --glob --stats-csv --fps --loop
--no-annotate --keep-frames`, plus the shared style flags.

## Incision pulses (`lsdviz3-lem-pulse`)

```bash
pixi run lsdviz3-lem-pulse --in runs/pulse --prefix pulse --out figures/pulse.gif
```

Each frame is a hillshade with the erosion rate draped over it, the channel
network, and a panel of the run's history. Needs a run with
`write_erosion_rate_raster` and `write_channel_mask_raster`. Flags: `--in --out
--prefix --ext --stats-csv --erosion-cmap --erosion-vmin --erosion-vmax
--erosion-percentile --erosion-transparency --channel-color --no-timeseries
--no-annotate --fps --loop --keep-frames`, plus the shared style flags.

## Channelization-index ladders (`lsdviz3-lem-ci`)

```bash
pixi run lsdviz3-lem-ci --in runs/ci --out figures/ci.gif --view oblique
```

Animates a channelization-index ladder run: the CI history panel beside a map
or an oblique view of the landscape. Flags: `--in --out --view --prefix --ext
--stats-csv --no-initial --map-width --azimuth --elevation
--vertical-exaggeration --landscape-cmap --fps --frame-dpi --loop
--keep-frames`, plus the shared style flags.

Python: `animate_transient_rasters`, `animate_lem_pulse` and `animate_ci_run`,
with the frame finders `find_transient_frames`, `find_pulse_frames` and
`find_ci_frames`.

## Hillslope solver comparisons

Figures comparing lsdtt3's hillslope solvers (Newton and cellular) with each
other and with analytical solutions. Each module docstring explains what the
figure tests and how to read it; read it before quoting a result.

| Function | Module | What it draws |
|---|---|---|
| `plot_estar_rstar_model_scan` | `charts.estar_rstar` | Analytical, Newton and cellular strip results in E*–R* space. |
| `plot_hillslope_profiles` / `plot_profile_residual_panels` | `charts.hillslope_profiles` | Elevation, gradient, curvature and residual panels for 1-D hillslope comparisons. |
| `plot_residual_density` / `plot_residual_raster` | `charts.hillslope_residuals` | Residual histograms (optionally hillslope-masked) and residual maps centred on zero; `summarise_residuals` tabulates RMS/L∞/percentiles. |
| `plot_transient_residuals` / `plot_strip_transient_metrics` | `charts.hillslope_transient_residuals` | RMS and max residuals through time; relief, max gradient and divide Laplacian changes per model. |

```bash
# paired Newton vs cellular hillslope solver animation
pixi run lsdviz3-hillslope-comparison --newton runs/newton --cellular runs/cellular \
  --out demo_output/hillslope_comparison.gif --layout 2d
```

The animation is built by
`render.hillslope_comparison_animation.render_hillslope_comparison_animation`;
`examples/hillslope_comparison.py` drives the whole comparison workflow.

The Perron and Ren benchmark figures are made by `examples/perron_ren_p1_p2.py`
… `perron_ren_p6.py` (P1–P2 outputs, P3 timestep independence, P4 transient
convergence and runtimes, P5 stability envelope, P6 spatial convergence and
anisotropy).
