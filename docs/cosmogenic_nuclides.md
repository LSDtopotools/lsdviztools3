# Cosmogenic nuclides

Figures from lsdtt3's cosmogenic tools: comparisons against CRONUS and
published data, transient (step-change) and landslide experiments, and basin
erosion-rate maps.

## Chi and cosmogenic basins

**Basins coloured by erosion rate (`lsdviz3-cosmo-basins`).** Spatially joins
point erosion rates from a cosmogenic CSV to basin polygons and draws the
basins, coloured by erosion rate, over the standard hillshade, with the
scaling, muon and atmosphere schemes in the caption:

```bash
pixi run lsdviz3-cosmo-basins --base runs/site_DEM.tif \
  --basins runs/site_junction_basin_polygons.fgb --cosmo-csv runs/site_erates.csv \
  --out figures/site_erates.png --colorbar
```

Flags: `--base --basins --cosmo-csv --out [--hillshade --joined-basins-out
--colorbar --scaling-scheme --muon-scheme --atmosphere --erate-input-units
--allow-incomplete-basins]`. Python: `join_basins_with_cosmo_erates`,
`render_cosmo_basin_erate_map`.

**Chi with erosion rates (`lsdviz3-chi-cosmo`).** From a chi run's fit tables,
writes the stacked and per-basin chi plots, a per-basin k_sn summary and map,
and, given `--cosmo-csv`, k_sn against erosion rate with error bars:

```bash
pixi run lsdviz3-chi-cosmo --fit-nodes-csv runs/site_chi_fit_nodes_mn_0p450.csv \
  --fit-segments-csv runs/site_chi_fit_segments_mn_0p450.csv --outdir figures \
  --base runs/site_DEM.tif --basins runs/site_chi_basin_polygons.fgb \
  --cosmo-csv runs/site_erates.csv --theta 0.45
```

Flags: `--fit-nodes-csv --fit-segments-csv --outdir [--base --hillshade
--basins --basin-id-col --cosmo-csv --erate-col --lon-col --lat-col
--sample-col --theta --scaling-scheme --atmosphere --muon-scheme
--erate-input-units --allow-incomplete-basins --stem --dpi --cmap]`. Python:
`render_chi_and_cosmo_basin_workflow`, `summarize_basin_ksn`,
`render_chi_basin_ksn_map`, `charts.plot_ksn_vs_erate`.

## Cosmo comparison workflows

For cosmogenic workflows, `lsdviztools3` can render standard comparison figures
from lsdtt3 CSV artefacts.

Expected primary input:

- `single_point_scan_comparison.csv`, written by
  `tests/scripts/run_cosmo_single_point_scaling_scan.py` in a clone of lsdtt3
  into the folder its `--outdir` names (it is generated, not tracked), for
  example:

  ```bash
  cd ../lsdtt3
  python3 tests/scripts/run_cosmo_single_point_scaling_scan.py \
      --pred-exe build/lsdtt-cosmo-prediction --tool-exe build/lsdtt-cosmo-tool \
      --with-cronus-v3 --sample-age-yr 0 \
      --outdir ../lsdviztools3/scan_out
  ```

  `--with-cronus-v3` queries the online CRONUS v3 calculator for the reference
  values that `reference="cronus_v3"` below compares against (without it those
  columns are empty), and `--sample-age-yr 0` matches how CRONUS treats the
  samples.

The published-summary inputs below are tracked in lsdtt3, under
`tests/expected/cosmo/published/`; the commands assume a clone of lsdtt3 next to
this one.

Optional secondary input:

- sample-level comparison tables (for example OCTOPUS/riversand summary style
  CSVs containing a percent-difference column).

Python API:

```python
from lsdviztools3.charts import load_cosmo_table, plot_cosmo_sweep_panels
from lsdviztools3.workflows import render_cosmo_comparison_figures

scan = load_cosmo_table("scan_out/single_point_scan_comparison.csv")
fig, axes = plot_cosmo_sweep_panels(scan, reference="cronus_v3", representative_pressure_hpa=900)

outputs = render_cosmo_comparison_figures(
    scan_comparison_csv="scan_out/single_point_scan_comparison.csv",
    sample_comparison_csv="../lsdtt3/tests/expected/cosmo/published/dethier2014_riversand_summary.csv",
    out_dir="demo_output/cosmo",
)
print(outputs)
```

Scripted example:

```bash
python examples/cosmo_comparison_workflow.py \
  --scan-csv scan_out/single_point_scan_comparison.csv \
  --outdir demo_output/cosmo

# Published-summary residuals (stone/lsdn/riversand)
pixi run lsdviz3-cosmo-published \
  --published-dir ../lsdtt3/tests/expected/cosmo/published \
  --datasets dethier2014,clapp2002,dibiase2010 \
  --mode auto \
  --outdir demo_output/cosmo
```

Main chart helpers:

- `plot_cosmo_sweep_panels`: nuclide x scaling panels of residual vs latitude,
  including separate muon-scheme traces.
- `plot_cosmo_sample_residuals`: grouped boxplots for sample-level residual
  tables, optionally split by `scaling_scheme`.

Main workflow helpers:

- `render_cosmo_comparison_figures`: sweep residual panels and optional
  sample-level residual figure.
- `render_cosmo_published_residuals`: dataset-level residual figure directly
  from `<dataset>_summary.csv`, `<dataset>_lsdn_summary.csv`, and
  `<dataset>_riversand_summary.csv` style files.
## Transient cosmogenic and landslide charts

Charts over the [`lsdtt-cosmo-column-tool`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_cosmo_column_tool.html) `step_change_sweep` CSVs
(`cosmo_transient_sweep_<label>.csv`) and the per-sample
`<stem>_joined_samples.csv` from a routed step-change landscape. Every quantity
drawn is a column lsdtt3 already wrote: nothing is normalised or inverted here.

| Function | What it draws |
|---|---|
| `plot_cosmo_transient_concentration` | Normalised concentration vs time since the step change, one series per nuclide (column prefix `approach_frac_`, `norm_conc_` or `norm_full_`). |
| `plot_cosmo_transient_erate_ratio` | Apparent erosion-rate ratios vs time since the step change. |
| `plot_nuclide_removal_response` | Three panels explaining why stripping a surface moves the apparent ¹⁴C/¹⁰Be ratio. |
| `plot_ratio_vs_steptime` | Per-sample apparent-rate ratio vs time since the knickpoint passed, coloured by true basin-mean erosion rate. |
| `plot_ratio_vs_mean_erate` | The same ratio vs true basin-mean erosion rate, coloured by time since the step. |
| `plot_ratio_overlay_1d` | `plot_ratio_vs_steptime` with the equivalent 1-D column sweep drawn over it. |
| `plot_scenario_delta` | Paired change in the apparent-rate ratio, one series per scenario. |
| `plot_landslide_area_volume_larsen` | Realised landslide area–volume pairs (`*_LS_INVENTORY.csv`) against a Larsen-style power law. |
| `plot_sample_effect_apparent_vs_true` / `plot_sample_effect_ratio` | Landslide-vs-control apparent rates from `*_LS_SAMPLE_EFFECT.csv`. |

Loaders: `load_cosmo_transient_table`, `load_joined_samples`,
`load_landslide_inventory`. The sweep figures also have a CLI:

```bash
pixi run lsdviz3-cosmo-transient \
  --sweep-csv cosmo_transient_sweep_fast.csv cosmo_transient_sweep_slow.csv \
  --outdir demo_output/cosmo_transient \
  --families concentration erate_ratio
```
