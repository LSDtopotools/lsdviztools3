# Valleys and cross sections

Floodplains, terraces, valley floors, valley widths and cross sections from
[`lsdtt-fluvial-surfaces`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_fluvial_surfaces.html), [`lsdtt-flat-surface-finder`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_flat_surface_finder.html),
[`lsdtt-valley-centreline`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_valley_centreline.html) and [`lsdtt-channel-extraction`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_channel_extraction.html)
(`run_perimeter_sources`). Most figures need a hillshade on the run's grid,
from [`lsdtt-raster-preprocessing`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_raster_preprocessing.html) (`write_hillshade_raster`) or
`lsdtt-valley-centreline` (`write_hillshade_raster`).

## Floodplain maps (`lsdviz3-floodplain`)

```bash
pixi run lsdviz3-floodplain --stem runs/site --base runs/site_hillshade.tif \
  --drape runs/site_reference_field.tif \
  --out-map figures/site_floodplain.png --out-qq figures/site_hand_qq.png
```

For an `lsdtt-fluvial-surfaces` run (`--stem` = path + `write_fname`) with
`write_floodplain_polygons: true`. `--out-map` draws the floodplain polygons
over the hillshade: the first pass in blue and, where a second pass
(`floodplain_second_pass_algorithm_choice`) changed it, the second in orange.
`--drape` colours a raster underneath, normally `<stem>_reference_field.tif`
(HAND). `--out-qq` draws the HAND quantile–quantile plots the floodplain bounds
were read from (`clubb2017` / `qq_gaussian` runs).

- `--levels` maps `<stem>_surface_patch_polygons` by ladder level instead of
  the floodplain passes; `--colours` sets the colours (two passes, or one per
  level).
- `--out-training` draws HAND box-and-whiskers per training polygon
  (`<stem>_training_polygons.csv`); `--training-log-y` puts HAND on a log axis.
- `--base-is-dem` illuminates `--base` instead of drawing it flat;
  `--fill-transparency`, `--outline-width`, `--outline-transparency`,
  `--drape-label`, `--drape-vmin`, `--drape-vmax`, `--drape-cmap`,
  `--figure-size` and `--dpi` set the look.

## Comparing floodplain and terrace methods (`lsdviz3-surface-methods`)

```bash
pixi run lsdviz3-surface-methods --stem runs/site --base runs/site_hillshade.tif \
  --out-classes figures/site_classes.png --out-levels figures/site_levels.png
```

For an `lsdtt-fluvial-surfaces` run, including one that ran several methods at
once. `--out-classes` maps each method's floodplain/terrace classes side by
side; `--out-levels` plots patch median HAND per method (`--levels-log-y` for a
log axis). With mapped surfaces (`--mapped` polygons, outlined on the maps),
`--out-agreement` writes `<stem>_floodplain.png` and `<stem>_terrace.png`
showing where each method agrees with the mapping, and `--out-scores` plots the
scores against them. `--methods` restricts the figures to a subset.

## Fluvial-surface ladder and valley sections

`lsdtt-fluvial-surfaces` bands a height reference (normally HAND) per
along-channel window into a "ladder" of levels. These figures check whether the
levels it found are real.

| Function | CLI | Input | What it draws |
|---|---|---|---|
| `plot_ladder_windows` | `lsdviz3-ladder-windows` | `<stem>_ladder_windows.csv` | Break selection along each channel: breaks, level centres, and windows the DEM could not resolve. `summarise_ladder` gives one row per channel. |
| `plot_ladder_density` | `lsdviz3-ladder-density` | `<stem>_ladder_histograms.csv` (`write_ladder_histograms_csv: true`) + windows CSV | Each window's HAND density with the breaks and their separation ratios. |
| `plot_ladder_band_map` | `lsdviz3-ladder-bands` | `<stem>_band_index.tif` | The band index over a hillshade. Bands are per window, so read patches, not a continuous surface. |
| `plot_valley_sections` | `lsdviz3-valley-sections` | `<stem>_valley_cross_section_samples.csv` | One panel per station: ground profile, picked valley floor, low point, and an optional HAND field coloured by band. `summarise_sections` gives width and centreline-to-low-point offset. |
| `plot_valley_section_map` | `lsdviz3-valley-section-map` | `<stem>_valley_cross_sections.fgb` (+ centreline, channel, sides) | Section lines against the centreline and channel network, flagging sections that cross no channel cell. |
| `plot_valley_width_map`, `plot_valley_widths_along_channel`, `plot_channel_and_centreline`, `plot_centreline_profile`, `plot_terrace_profile` | `lsdviz3-valley-widths` | `<stem>_valley_cross_sections.fgb`, `_valley_widths.csv`, `_valley_segments.csv`, `_valley_segment_channels.fgb`, `_valley_centreline.fgb`, `_valley_centreline_nodes.csv` (`run_valley_width_network: true`) | Every segment's section width lines on the hillshade, coloured by `excluded_reason`; valley width down one segment (default the trunk reaching the largest drainage area) with tributary confluences marked and drainage area in a panel below; each segment's extracted channel against its valley centreline; and down one centreline the DEM under each node against the centreline baseline, with their difference (any lsdtt-valley-centreline run); and terrace patches binned along a centreline (`write_terrace_profile_csv`: `<stem>_terrace_profile.csv`, `_terrace_patches.csv`) as 25-75 percentile bars over the valley floor profile, numbered by patch and coloured by level. |

```bash
pixi run lsdviz3-ladder-windows --windows run_ladder_windows.csv --out demo_output/ladder_windows.png
pixi run lsdviz3-ladder-density --histograms run_ladder_histograms.csv \
  --windows run_ladder_windows.csv --out demo_output/ladder_density.png
pixi run lsdviz3-ladder-bands --hillshade run_hs.tif --bands run_band_index.tif \
  --channel run_channel_network_lines.fgb --out demo_output/ladder_bands.png
pixi run lsdviz3-valley-sections --samples run_valley_cross_section_samples.csv \
  --field run_hand.tif --bands run_band_index.tif --out demo_output/valley_sections.png
pixi run lsdviz3-valley-section-map --hillshade run_hs.tif \
  --sections run_valley_cross_sections.fgb --centreline run_valley_centreline.fgb \
  --channel run_channel_network_lines.fgb --out demo_output/valley_section_map.png
pixi run lsdviz3-valley-widths --stem run_W --base run_hs.tif \
  --valley-mask run_floodplain_mask.tif \
  --out-map demo_output/valley_widths_map.png --out-profile demo_output/valley_widths_profile.png \
  --out-channels demo_output/valley_widths_channels.png \
  --out-centreline-profile demo_output/centreline_profile.png \
  --out-terrace-profile demo_output/terrace_profile.png --max-relative-elevation 120
```

These sit alongside the matplotlib cross-section charts in
[Valley cross sections](#valley-cross-sections), which work from the same samples table.
## Valley floor / HAND charts

**These charts read files that current lsdtt3 no longer writes**: the HAND
profile switch (`write_hand_profile_csv`) was retired from
`lsdtt-flat-surface-finder`, and the band-break tables went with it. Use them
with output from older lsdtt3 runs; for current runs see the floodplain and
ladder figures above.

Charts over `lsdtt-flat-surface-finder` output: HAND (height above nearest
drainage) per pixel in `*_hand_profile.csv`, the reach-relative band breaks in
`*_hand_band_breaks.csv`, and the reach-id / valley-floor rasters. Surface
classes follow the lsdtt3 `_surface_class` convention (1 = floodplain,
2 = terrace).

| Function | What it draws |
|---|---|
| `plot_hand_distance` | HAND against distance along the channel, coloured by surface class, optionally with the band-break levels. |
| `plot_hand_window_populations` | Per-window HAND populations with the selected valley break. |
| `plot_valley_break_metric` | Reviewer figure: how the valley HAND break is chosen. |
| `plot_valley_break_fraction` | Per-reach relief-adaptive acceptance fraction along the channel. |
| `plot_valley_reach_map` | Map of the along-channel reaches over a hillshade. |
| `plot_valley_floor_map` | Map of the extracted valley floor over a hillshade. |

The map functions need `rasterio` in addition to matplotlib. Full parameter
lists are in the function docstrings (`lsdviztools3/charts/hand.py`).
## Valley cross sections

Charts over `lsdtt-valley-centreline` output: the per-sample transect table
`*_valley_cross_section_samples.csv` (`write_valley_cross_sections`) and the
feature table written by `write_cross_section_features`. Every transect on a
channel shares a length set by the widest valley on that channel, so panels are
directly comparable.

| Function | What it draws |
|---|---|
| `plot_valley_cross_sections` | Grid of cross sections orthogonal to the centreline, with the in-valley run shaded and the detected sides marked. |
| `plot_valley_cross_section_map` | Map of the numbered transects over a hillshade. |
| `plot_valley_curved_cross_section_map` | The same for the curved (boundary-fitted) sections. |
| `plot_valley_cross_section_features` | Individual sections annotated with detected features. |
| `plot_valley_cross_section_features_map` | Map of transects with valley-edge and terrace/floodplain features. |
| `plot_valley_cross_section_hydraulics` | Flow-width step surfaces and terrace/floodplain edges, per section. |

`cross_section_hydraulics`, `detect_hydraulic_flats`, `combine_width_steps` and
`compute_flat_surface_edges` are the numeric helpers behind the hydraulics
plot. They are **prototyping code**: lsdviztools3 visualises, lsdtt3 computes,
so treat their results as exploratory until the equivalent lives in C++.
## Hydraulic breaks

Cross-section flooding (`lsdtt-valley-centreline` with
`write_cross_section_features`) records the water-surface elevations separating
flat surfaces from the risers between them, written as local HAND values in
`*_valley_cross_section_hydraulic_breaks.csv`.

| Function | What it draws |
|---|---|
| `plot_hydraulic_breaks_downstream` | Break HAND of every section against distance downstream, one series per break rank. |
| `hydraulic_break_clusters` | Binned 1-D clustering test (default 10 sections per bin) asking whether break elevations fall into distinct levels. Numeric, not a plot — same prototyping caveat as [Valley cross sections](#valley-cross-sections). |
| `plot_hydraulic_break_clusters` | Break-HAND scatter down-valley with the per-bin cluster levels overlaid. |
| `plot_cross_section_informed_valley_mask` | The cross-section-informed valley floor over a hillshade, beside the histogram-derived mask. |

## Perimeter valley sources (`lsdviz3-perimeter-sources`)

```bash
pixi run lsdviz3-perimeter-sources --stem runs/site --base runs/site_hillshade.tif \
  --out-map figures/site_sources_map.png --out-profile figures/site_sources_profile.png
```

For an `lsdtt-channel-extraction` run with `run_perimeter_sources`, which finds
where valleys enter the DEM across its edge. `--out-map` draws the flat runs,
outlets and sources on the hillshade (or one computed from `--dem`), with the
`--n-unselected-outlets` largest unselected outlet runs (default 3);
`--out-profile` plots elevation and gradient round the perimeter, with a
threshold line at `--max-gradient` (default the run's
`perimeter_max_gradient`). `--bbox` zooms the map.
