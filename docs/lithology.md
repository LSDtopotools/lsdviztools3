# Lithology

Geologic maps and channel statistics by rock type. The statistics come from
lsdtt3; these commands draw them:

| lsdtt3 run | What it writes | Command |
|---|---|---|
| [`lsdtt-fetch-lithology`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_fetch_lithology.html) | `<prefix>_macrostrat_units.fgb` (Macrostrat geologic map) | — (input to the runs below) |
| [`lsdtt-chi-analysis`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_chi_analysis.html) with `lithology_fname` | k_sn by lithologic unit, `<prefix>_lithology_code.tif` (`write_lithology_code_raster`) | `lsdviz3-ksn-by-lithology`, `lsdviz3-lithology-map` |
| `lsdtt-chi-analysis` with `erosion_rate_fname` | erodibility calibration (erosion rate against k_sn) | `lsdviz3-erodibility-calibration` |
| [`lsdtt-channel-network-analysis`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_channel_network_analysis.html) with `lithology_fname` and `channel_bearing_group_by: lithology` | bearing and junction-angle statistics by unit | `lsdviz3-channel-bearings --group-key lithology` ([Channel networks](channel_networks.md#bearing-and-junction-angle-statistics-lsdviz3-channel-bearings)) |

All three lithology commands accept `--unit-colours`: a `name,colour` CSV, or a
JSON `{name: colour}`, fixing each unit's colour by name. Use the same file for
every figure of a study so a unit has one colour throughout.

## Lithology map (`lsdviz3-lithology-map`)

```bash
pixi run lsdviz3-lithology-map --codes runs/site_lithology_code.tif \
  --dem runs/site_DEM.tif --out figures/site_lithology_map.png
```

Draws the unit code raster with a legend of unit names, a credit line for the
geologic map and the DEM, a scale bar and a north arrow. By default the unit
colours are semi-transparent over a hillshade of `--dem`; `--no-hillshade`
draws flat opaque colours.

- **Unit names** are found next to the code raster
  (`<prefix>_lithology_units.csv`, then the chi k_sn tables), or given with
  `--names` (a `code,name` CSV) or `--crosswalk`. A unit with no name stops the
  plot with a message: lsdtt3 writes an empty name for map polygons whose
  lithology field is empty, so give such a unit a name with `--names`.
- **Channels:** `--channels` overlays a line or point vector, or the chi
  lithology node CSV (`*_chi_ksn_lithology_nodes_mn_*.csv`), in one colour, by
  stream order (`--channel-order-column`), or by k_sn
  (`--channel-color-by k_sn`, which mutes the geology). `--network-status`
  draws a channel-network-analysis run's whole network, analysed links thick.
- **Regions:** `--regions` outlines region polygons (e.g. physiographic
  provinces), named by `--region-field`.
- `--hide-unmapped` leaves unmapped cells out; `--legend-title`,
  `--legend-columns`, `--transparency`, `--relief` and the shared style flags
  control the look; `--credit`, `--attribution`, `--dem-source` or
  `--no-credit` set the credit line.

## k_sn by lithology (`lsdviz3-ksn-by-lithology`)

```bash
pixi run lsdviz3-ksn-by-lithology --data-dir runs --prefix site_chi \
  --m-over-n 0.45 --out-dir figures
```

Reads the lithology tables of a chi run and writes, prefixed `<out-prefix>_`:

| File | Contents |
|---|---|
| `ksn_hist_length.png` | Length-weighted k_sn histograms per unit |
| `ksn_hist_count.png` | Node-count histograms per unit |
| `ksn_box.png` | Per-unit box plots with sample sizes |
| `reach_ksn_box.png` | Per-reach k_sn, when the reach table exists |
| `ksn_unit_basin.png` | Unit-by-basin dots, when that table exists |

`--order` sets the unit order (`auto`, `median`, `reach_median`, `code`,
`name`, `length`); `--layout` overlays the histograms (up to 6 units by
default) or draws them as facets; `--raw` plots raw length or count instead of
each unit's fraction; `--box-source` takes the boxes from the per-node values
or the summary quartiles.

Before reading a difference between units as an erodibility contrast, check the
caveats lsdtt3 prints with these tables: a knickzone crossing a unit, or
strike-parallel drainage that puts a unit mostly at channel heads, moves its
k_sn for reasons other than rock strength.

## Erodibility calibration (`lsdviz3-erodibility-calibration`)

```bash
pixi run lsdviz3-erodibility-calibration --data-dir runs --prefix site_cal \
  --m-over-n 0.45 --out-dir figures
```

For a chi run given basin-averaged erosion rates (`erosion_rate_fname`), this
draws `<out-prefix>_erosion_rate_vs_ksn.png` (erosion rate against k_sn on log
axes, with lsdtt3's fit) and the sample-basin maps
`<out-prefix>_basin_map_status.png` and `_basin_map_unit.png` (coloured by
whether each basin was used, and by its dominant unit). `--plots` picks
`scatter` and/or `map`; `--map-colour-by` picks `status`, `unit` or both; the
map's hillshade comes from `--dem` or the run's own DEM.
