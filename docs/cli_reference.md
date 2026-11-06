# Command reference

Every `lsdviz3-*` command, grouped as in this guide. Run one with
`pixi run lsdviz3-<name> ...` from the lsdviztools3 checkout, and add `--help`
for its full list of flags. Each group links to the page that explains its
figures.

## Basic plots

See [Basic plots](basic_plots.md).

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-hillshade` | shaded-relief PNG | `--in --out --input-is-hillshade` + style (incl. `--clean`) |
| `lsdviz3-drape` | gray hillshade + coloured semi-transparent drape | `--base --drape --out --shade-base --log --transparency --vmin --vmax --percentile-min --percentile-max --out-of-bounds --colorbar --colorbar-label --colorbar-position` + style (`--cmap` = drape; GMT or matplotlib maps, or a solid colour) |
| `lsdviz3-basins` | hillshade + basin polygons | `--base --basins --out --fill --transparency` + style |
| `lsdviz3-basemap` | global/aerial locator | `--in --out --mode --elevation-km --view-limit-deg --look-bearing --look-vertical-angle --grid-spacing-deg --no-north-arrow` + style |
| `lsdviz3-sti` | hillshade + STI rasters (single file or directory batch) | single mode: `--base --sti --out`; directory mode: `--sti-dir [--base]` (auto-plots all discovered STI rasters); options: `--transparency --colorbar` + style |

## Data and lsdtt3

See [Data and lsdtt3](getting_started.md).

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-grabdem` | download a DEM | `--bbox --source --dataset --out --dir --ot-api-key-file` |
| `lsdviz3-lsdtt` | run an lsdtt3 tool | `tool --params --data-path --lsdtt-dir --set-lsdtt-dir --cwd --dir` |

## Chi analysis

See [Chi analysis](chi_analysis.md).

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-chi` | an lsdtt-chi-analysis run (`fit_chi_networks: true`): k_sn map of the fitted segments over a hillshade, stacked chi-elevation plot, optional per-basin profiles | `--data-dir --prefix --m-over-n --out-dir [--dem --out-prefix --per-basin --no-stacked --cmap --ksn-vmin --ksn-vmax --remove-ksn-above --min-tributary-area-m2 --min-tributary-area-fraction --top-tributaries --keep-descendants --show-removed --coordinates --dpi]` |

## Channel networks

See [Channel networks](channel_networks.md).

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-channels` | hillshade + channel lines | `--base --lines --out` + style |
| `lsdviz3-channel-bearings` | channel-link bearing roses, junction-angle summaries and comparison tables by group, from lsdtt-channel-network-analysis | `--data-dir --prefix --out-dir [--group-key --stratum-key --plots --pdf-weight --measurands --compare-prefix --compare-label --map-dem --regions --region-field --region-title --label --title --dpi]` |
| `lsdviz3-burned-profiles` | long profiles of channels burned by lsdtt-stream-burning `burn_to_corridor`, one panel per run: input DEM, burned channel, minima, filled pits, cuts and raises | `--csv (repeat) --out [--labels --name \| --path-rank \| --path-id --x --mark --max-join-labels --title --dpi]` |
| `lsdviz3-divide-migration` | `lsdtt-divide-migration` figure set: hillshade + divide polylines + channel heads coloured by a metric, the side-vs-side pair chart (Δ = A − B, A the lower-numbered unit, error bar only when the run's uncertainty method is a dispersion), and a per-head scatter against distance to divide | `--prefix --outdir [--base --base-is-hillshade --metric --point-size --divide-color --alpha --stem --chart-dpi]` + style |

## Lithology

See [Lithology](lithology.md).

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-lithology-map` | geologic unit map with legend, credit, scale bar and north arrow, optionally with channels or regions | `--codes --out [--dem --names --crosswalk --channels --network-status --channel-color-by --regions --unit-colours --hide-unmapped --no-hillshade --legend-title --credit --no-credit]` + style |
| `lsdviz3-ksn-by-lithology` | k_sn histograms, box plots and unit-by-basin dots per lithologic unit | `--data-dir --prefix --out-dir [--m-over-n --out-prefix --order --layout --raw --box-source --unit-colours --title --dpi]` |
| `lsdviz3-erodibility-calibration` | erosion rate against k_sn and the sample-basin maps of an erodibility calibration | `--data-dir --prefix --out-dir [--m-over-n --out-prefix --plots --map-colour-by --dem --display-downsample --unit-colours --credit --title --dpi]` |

## Valleys and cross sections

See [Valleys and cross sections](valleys_and_cross_sections.md).

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-floodplain` | floodplain polygons (first and second pass) over a hillshade, and the HAND QQ plots | `--stem [--base --base-is-dem --out-map --levels --out-qq --out-training --training-log-y --drape --drape-label --drape-vmin --drape-vmax --drape-cmap --colours --fill-transparency --outline-width --outline-transparency --figure-size --dpi]` |
| `lsdviz3-surface-methods` | compares the floodplain/terrace methods of one lsdtt-fluvial-surfaces run | `--stem [--base --mapped --methods --out-classes --out-agreement --out-levels --levels-log-y --out-scores --panel-width --dpi]` |
| `lsdviz3-valley-widths` | valley-width network: section map by excluded_reason, width profile down the largest trunk, channel vs centreline map | `--stem [--base --valley-mask --out-map --out-profile --out-channels --out-centreline-profile --out-terrace-profile --terrace-field --include-floodplain --min-bin-pixels --max-relative-elevation --segment-id --map-segments --bbox --hide-excluded --title --dpi]` |
| `lsdviz3-valley-sections` | valley cross-section panels with picked floor and optional HAND field | `--samples --out [--source-id --field --bands --field-label --max-panels --n-cols --summary --title --dpi]` |
| `lsdviz3-valley-section-map` | cross-section lines vs centreline and channel network | `--hillshade --sections --out [--centreline --channel --samples --sides --source-id --label-every --title --dpi]` |
| `lsdviz3-ladder-windows` | fluvial-surface ladder break selection along each channel | `--windows --out [--summary --title --dpi]` |
| `lsdviz3-ladder-density` | per-window HAND density with ladder breaks | `--histograms --windows --out [--source-id --max-panels --title --dpi]` |
| `lsdviz3-ladder-bands` | ladder band index over a hillshade | `--hillshade --bands --out [--channel --centreline --max-band --title --dpi]` |
| `lsdviz3-perimeter-sources` | map and perimeter profile of the valley sources from run_perimeter_sources | `--stem [--base --dem --out-map --out-profile --n-unselected-outlets --max-gradient --bbox --title --dpi]` |

## Hillslopes and ridgelines

See [Hillslopes and ridgelines](hillslopes_and_ridgelines.md).

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-hilltop-coupling` | hilltop/channel coupling figure set | `--prefix --outdir [--dem --channels --max-traces --dpi]` |
| `lsdviz3-estar-rstar` | measured E*–R* vs Roering et al. (2007) theory + residual summary CSV | `--metrics --outdir [--group-col --title --dpi]` |
| `lsdviz3-ridgeline-curvature` | ridge-crest elevation, detrending, maximal curvature and Laplacian along one ridgeline | `[--driver] --ridgeline-profile-csv [--ridgeline-number --reverse-ridgeline --use-distance-range --minimum-distance --maximum-distance --moving-window-sizes-pixels --curvature-lambdas-m --site-display-name] --write-path --write-name` |
| `lsdviz3-ridgeline-summary` | triple divide junctions, a ridgeline series with curvature barcode, and a ridgeline map | `[--driver] --ridgeline-profile-csv [--hillshade-raster-file --native-tdj-csv-1..3 --topological-tdj-csv-1..3 --channel-network-vector-file --paired-flowpaths-vector-file --ridgeline-series --highlighted-tdj-type --highlighted-tdj-threshold ...] --write-path --write-name` |
| `lsdviz3-ridgeline-spacing` | maximal-curvature zero-crossing spacing against polyfit scale beside slope–area, per site | `[--driver] --site NAME PROFILE_CSV SLOPE_AREA_CSV (repeat) [--site-count --log-spacing --spacing-y-max-m ...] --write-path --write-name` |
| `lsdviz3-mdd-ridgeline` | main-drainage-divide ridgeline profile above the model's drainage-area image | `[--driver] --ridgeline-profile-csv --drainage-area-map-png [--model-name --ridgeline-number --guide-interval] --write-path --write-name` |

## Landscape evolution

See [Landscape evolution](landscape_evolution.md).

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-animate` | animate a transient `lsdtt-lem` run as a GIF | `--in --out --glob --stats-csv --fps --loop --no-annotate --keep-frames` + style |
| `lsdviz3-lem-pulse` | animate an `lsdtt-lem` incision-pulse run (hillshade + erosion drape + channels + history panel) | `--in --out --prefix --ext --stats-csv --erosion-cmap --erosion-vmin --erosion-vmax --erosion-percentile --erosion-transparency --channel-color --no-timeseries --no-annotate --fps --loop --keep-frames` + style |
| `lsdviz3-lem-ci` | animate an `lsdtt-lem` channelization-index ladder run (CI history panel + map or oblique landscape) | `--in --out --view --prefix --ext --stats-csv --no-initial --map-width --azimuth --elevation --vertical-exaggeration --landscape-cmap --fps --frame-dpi --loop --keep-frames` + style |
| `lsdviz3-hillslope-comparison` | animate paired Newton vs cellular hillslope rasters | `--newton --cellular --out [--glob --stats-csv --fps --loop --keep-frames --layout]` |

## Cosmogenic nuclides

See [Cosmogenic nuclides](cosmogenic_nuclides.md).

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-cosmo-basins` | standard hillshade + basins colored by inferred erosion rate from cosmo CSV | `--base --basins --cosmo-csv --out [--hillshade --joined-basins-out --colorbar --scaling-scheme --muon-scheme --atmosphere --erate-input-units --allow-incomplete-basins]` |
| `lsdviz3-chi-cosmo` | chi basin workflow bundle: stacked chi + per-basin chi PNGs + basin k_sn summary/map + optional k_sn-vs-erate | `--fit-nodes-csv --fit-segments-csv --outdir [--base --hillshade --basins --cosmo-csv --theta --scaling-scheme --atmosphere --muon-scheme --erate-input-units --allow-incomplete-basins --stem --dpi]` |
| `lsdviz3-cosmo-compare` | cosmogenic comparison figures from standard lsdtt3 CSVs | `--scan-csv --outdir [--sample-csv --reference --pressure --erate --stem --dpi]` |
| `lsdviz3-cosmo-published` | residual figure(s) from lsdtt3 published summary tables | `--published-dir --datasets [--mode --outdir --stem --dpi]` |
| `lsdviz3-cosmo-transient` | transient (step-change) cosmogenic figure families | `--sweep-csv --outdir [--value-prefix --families --concentration-xmax --erate-ratio-xmax --no-title --dpi]` |

Shared **style flags** (all PyGMT map/animation CLIs — the matplotlib chart
CLIs take only the flags listed in their row): `--coordinates --cmap --projection
--figure-size --line-width --line-color --font-annotation --font-label
--font-title --dpi --high-res --opaque --clean`.

The PyGMT map and animation commands, `lsdviz3-grabdem` and `lsdviz3-lsdtt`
also take `--dir DIR` — a target directory created if absent, into which a
relative `--out` (and, for `lsdviz3-lsdtt`, the lsdtt3 read/write/cwd) is
placed. See [Target directory](getting_started.md#target-directory---dir).
