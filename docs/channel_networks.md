# Channel networks

Channel network maps, channel-link bearings and junction angles, network
geometry across synthetic networks, divide migration, and the stream-burning
quality-control maps and burned long profiles. For k_sn on the network see
[Chi analysis](chi_analysis.md); for statistics by rock type see
[Lithology](lithology.md).

## Channel network map (`lsdviz3-channels`)

```bash
pixi run lsdviz3-channels --base work/eigg/eigg_hillshade.tif \
  --lines work/eigg/eigg_channel_network_lines.fgb --out demo_output/channels.png
```

Draws the channel lines from `lsdtt-channel-extraction`
(`write_channel_network_lines`) over the base raster, with the shared style
flags. In Python, `render_channels(lines, color_by="stream_order")` colours
them by any numeric attribute; see
[Chi analysis](chi_analysis.md#channel-network-coloured-by-log10-k_sn) for the
full signature.

## Bearing and junction-angle statistics (`lsdviz3-channel-bearings`)

`lsdtt-channel-network-analysis` with `write_channel_bearing_statistics` and
`write_junction_angle_statistics` groups channel links (by Strahler order,
lithology, region, ...) and compares the groups with circular statistics.
`lsdviz3-channel-bearings` draws that run:

```bash
pixi run lsdviz3-channel-bearings --data-dir runs/site --prefix site \
  --out-dir figures --group-key strahler_order --stratum-key none
```

| Flag | Meaning |
|---|---|
| `--data-dir`, `--prefix` | The run's directory and `write_fname`. |
| `--group-key` | `group_key`(s) to plot: `strahler_order`, `lithology`, `lithology_interior`, ... |
| `--stratum-key` | The run's `channel_bearing_stratify_by` (e.g. `strahler_order`), or `none`. **It must match the run**: with the wrong stratum no roses are drawn. |
| `--plots` | `rose`, `junction`, `tables`, `blocks`, `cdf`, `pdf` (default the first four). |
| `--pdf-weight` | Fraction of links (`count`) or of channel length (`length`). |
| `--measurands` | Junction measurands: `alpha`, `beta_min`, `beta_max`. |
| `--compare-prefix`, `--compare-label` | A second run drawn beside this one (e.g. an isotropic null). |
| `--map-dem`, `--regions`, `--region-field`, `--region-title` | A channel map coloured by stream order, optionally with region outlines. |

Outputs, prefixed `<prefix>_` (plus the group key when it is not the default):
roses per group (`rose_<axial|directional>_<count|length>[_<stratum>].png`),
junction-angle summaries and comparisons, the comparison tables, and the
block-order counts. The Python functions behind the roses and comparisons are
below.
## Rose diagrams (bearings)

A polar **rose / wind-rose** of bearings in degrees, drawn North-up and
clockwise (geographic convention). Pass either a flat `bearings` array for a
single rose, or `groups={label: bearings}` to compare/stack sub-populations —
e.g. channel-link bearings split by Strahler stream order.

```python
import pandas as pd
from lsdviztools3.charts import rose_diagram

links = pd.read_csv("work/site/site_channel_link_bearings.csv")

# one rose per stream order, stacked
groups = {f"order {o}": links.loc[links.stream_order == o, "bearing_deg"].to_numpy()
          for o in sorted(links.stream_order.unique())}
fig, ax = rose_diagram(groups=groups, bin_width_deg=10, mode="directional",
                       stacked=True, title="Channel bearings by order")
fig.savefig("demo_output/rose_bearings.png", dpi=300, bbox_inches="tight")

# a single rose from a flat array
fig, ax = rose_diagram(links["bearing_deg"].to_numpy(), bin_width_deg=15)
```

`rose_diagram(bearings=None, *, groups=None, bin_width_deg=15.0,
mode="directional", density=False, stacked=True, cmap="viridis", colors=None,
label_fmt="{}", title=None, ax=None)`. Returns `(fig, ax)` with `ax` a North-up,
clockwise `PolarAxes`.

- `mode="directional"` (default) keeps the full 360° (e.g. upstream→downstream
  channel links, outlet→centroid bearings); `mode="axial"` folds bearings to
  180° and mirrors each petal (orientation only, e.g. lineaments).
- `stacked=True` stacks grouped petals (low order at the base); `stacked=False`
  draws overlaid outline-only petals.
- `density=True` normalises each group to fractions so shapes can be compared
  regardless of sample size. `bin_width_deg` sets the angular bin width.
- `cmap` samples per-group colours; pass `colors={label: colour}` to override.
## Circular / directional comparison

Compare two angular samples (e.g. channel bearings from two DEMs) with circular
statistics and an overlay plot. These mirror the authoritative C++ circular
tests in lsdtt3's `lsdtt_statistical_tests`; for publication numbers use
lsdtt3's own comparison (`lsdtt-channel-network-analysis` with
`write_channel_bearing_statistics`, which writes
`<prefix>_channel_bearing_comparisons.csv`) and feed its result into the plot.

```python
import pandas as pd
from lsdviztools3.charts import (
    circular_wasserstein, kuiper_two_sample, plot_circular_comparison,
)

a = pd.read_csv("work/site_a/site_a_channel_link_bearings.csv")["bearing_deg"].to_numpy()
b = pd.read_csv("work/site_b/site_b_channel_link_bearings.csv")["bearing_deg"].to_numpy()

dist = circular_wasserstein(a, b)        # degrees; a true metric, handles wraparound
v, p = kuiper_two_sample(a, b)           # rotation-invariant circular KS analogue

fig, (ax_pdf, ax_cdf) = plot_circular_comparison(
    a, b, labels=("dendritic", "trellis"), bin_width_deg=10, kind="both")
fig.savefig("demo_output/bearing_comparison.png", dpi=300, bbox_inches="tight")
```

- `circular_wasserstein(a, b, *, period=360.0, rotation_invariant=False,
  rotation_step=1.0)` — circular 1-Wasserstein (earth-mover) distance in the
  input units (degrees). A true metric, handles wraparound (359° ≈ 1°); 0 for
  identical distributions, bounded by `period / 2`. No p-value. Use
  `period=180` for axial data; `rotation_invariant=True` also minimises over a
  rigid rotation of `b` (shape only, ignoring absolute orientation).
- `kuiper_two_sample(a, b, *, period=360.0) -> (V, p_value)` — two-sample
  Kuiper test, the rotation-invariant circular analogue of the
  Kolmogorov–Smirnov test, giving significance (`V = D⁺ + D⁻`, small `p` ⇒ the
  distributions differ).
- `plot_circular_comparison(a, b, *, labels=("A","B"), bin_width_deg=15.0,
  kind="both", period=360.0, annotate=True, annotation=None,
  rotation_invariant=False, styles=None, title=None, ax=None)` — overlays the
  two samples' PDFs and/or CDFs and annotates the distance. `kind="both"`
  returns `(fig, (ax_pdf, ax_cdf))`; `kind="pdf"` or `"cdf"` returns
  `(fig, ax)`. Pass `annotation=<str>` to inject externally-computed stats
  (e.g. from lsdtt3's `_channel_bearing_comparisons.csv`)
  instead of the internally-computed Wasserstein/Kuiper note.

These bearing analyses need **lsdtt3**'s `lsdtt-channel-network-analysis`
(see [Running lsdtt3](getting_started.md#running-lsdtt3)): `only_complete_basins`
and `min_basin_pixels` restrict the bearing and junction outputs to complete
basins, and `write_channel_bearing_statistics` computes the circular-Wasserstein
distance and Kuiper test between groups. (These comparisons were once a separate
`lsdtt-bearing-comparison` CLI, which lsdtt3 has removed.)

## Distributions

`plot_distributions(series, ...)` overlays step-histograms of several labelled
value series — used in the worked examples for junction-angle populations and
for fault-relative deflection split by side of the fault
([Python workflows](getting_started.md#python-workflows)).

```python
from lsdviztools3.charts import plot_distributions
fig, ax = plot_distributions({"dendritic": angles_a, "trellis": angles_b})
```

## Network geometry and synthetic networks

| Function | Module | What it draws |
|---|---|---|
| `plot_junction_angle_cdf`, `plot_delta_L_cdf`, `plot_edge_basin_cdf`, `plot_metric_panels`, `plot_metric_vs_temperature` | `charts.network_geometry` | Network-geometry metrics across a sweep (see the rows in [lsdtt3 outputs](lsdtt3_outputs.md)). |
| `plot_hillshade`, `plot_network`, `plot_network_on_hillshade` | `charts.network_maps` | Matching map panels per synthetic-network case: relief, the network alone (line width ~ log area), and both together. |
| `plot_metric_by_temperature`, `plot_branching_range`, `plot_flow_change_fraction` | `charts.snapping` | The OCN snapping experiment: does annealing temperature produce a range of branching, and does snapping reduce it? |

`examples/synthetic_networks_demo.py` drives `lsdtt-synthetic-networks` and
uses these charts end to end.
## Divide migration charts

`lsdviz3-divide-migration` ([Command reference](cli_reference.md)) renders the full figure set
for an `lsdtt-divide-migration` run. The pieces are available individually:

| Function | What it draws |
|---|---|
| `render.divide_migration.render_divide_migration` | PyGMT map: grey hillshade + divide polylines + channel heads coloured by one metric. |
| `plot_pair_comparison` | Side A and side B means for one metric, with Δ = A − B and its uncertainty. Side A is always the lower-numbered unit. |
| `plot_head_metric_vs_distance` | Per-head metric against distance to the nearest divide, split by unit. |
| `aggressor_annotation` | Text describing which side Δ points at, and whether the uncertainty supports it. |

Loaders live in `io.divide_migration` (`find_outputs`, `load_divides`,
`load_channel_heads`, `load_pair_summary`).

## Stream-burning quality control

These presets build on
[`hillshade_with_vectors`](basic_plots.md#overlaying-several-vector-layers-on-a-hillshade).
**Stream-burning preset.** `compare_stream_networks(base_raster, networks,
png_path, **kwargs)` is a thin wrapper over `hillshade_with_vectors` for
stream-burn QC — overlaying a flow-routed network against the vector layers that
conditioned the DEM (it just defaults `clip_to_valid=True`).

**Segment-classification preset.** `plot_segment_classification(base_raster,
segments, png_path=None, *, status_field="status", followed_value=0,
burned_gaps=None, flow_traces=None, extra_layers=None, clip_to_valid=True,
**kwargs)` maps the `lsdtt-stream-burning` trace-divergence burn
(`snap_to_accumulation` mode) segment by segment. It splits the `_segments.fgb`
diagnostic layer (written when `write_burn_diagnostics: true`) on its `status`
field and draws each fixed-length vector segment by how the burn treated it —
**followed** (`status == 0`; the D8 trace already carried the channel) vs
**unresolved** (still diverging after every pass, the cases to investigate) —
optionally with the burned trenches (`_burned_gaps.fgb`) and sampled D8 traces
(`_flow_traces.fgb`) on top. Forwards `region=` etc. to `hillshade_with_vectors`,
so the same call zooms to a road crossing. This is the most useful QC view for
the accurate-vector (VectorMap Local) burn.

```python
from lsdviztools3.workflows import plot_segment_classification

plot_segment_classification(
    "sd_acc_burned.tif", "sd_acc_segments.fgb", "qc_segments.png",
    burned_gaps="sd_acc_burned_gaps.fgb",
    flow_traces="sd_acc_flow_traces.fgb",
    region=[171600, 172400, 634200, 634800])  # native-CRS zoom
```

### Burned long profiles (`lsdviz3-burned-profiles`)

For `lsdtt-stream-burning` runs with `stream_burn_mode: burn_to_corridor`
(a vector network taken as the truth, e.g. OS Open Rivers prepared by
`lsdtt-fetch-river-lines`). In that mode lsdtt3 builds each burned channel's
elevations from the minima of the input DEM along the channel. Small pits
are filled first (`burn_profile_max_pit_pixels`, default 16). The minima are
joined by straight lines in distance, or in chi
(`burn_profile_fit_algorithm_choice: chi`). The rest of the DEM is then
carved and filled to drain to the channel. lsdtt3 writes
`<stem>_burned_profiles.csv` (switch `write_burned_profiles_csv`, on by
default), one row per channel cell. This command draws it.

```bash
pixi run lsdviz3-burned-profiles \
  --csv runs/site_x0_burned_profiles.csv \
  --csv runs/site_x60_burned_profiles.csv \
  --labels "X = 0 m" "X = 60 m" \
  --name "River Name" --out figures/site_burned_profiles.png
```

One panel per `--csv`, all following the same course on a shared x axis. A
course is a channel path from its head down to the network outlet, following
the paths it joins. Choose it with one of:

- `--name`: from the most upstream cell carrying that river name. The run
  needs `vector_name_field`; lines prepared by `lsdtt-fetch-river-lines`
  carry the name in `name`.
- `--path-rank`: 0, the default, is the first path lsdtt3 fitted, the
  longest.
- `--path-id`: one `--csv` only, since path ids differ between runs.

Each panel draws:

| Mark | Meaning |
|---|---|
| grey line | the input DEM under the channel |
| blue line | the burned channel in the output DEM |
| dots / crosses | minima the profile joins / minima rejected by `burn_profile_max_cut_metres` |
| orange band | cut below the input DEM |
| green band | a small pit filled before the profile was sampled, where the channel sits in it |
| purple band | raised above the input DEM for another reason: a tributary raised to meet the channel it joins, or a rejected minimum |
| grey ticks | where other paths join the course, the longest named |

The bands split each cell's change against the input DEM without overlap.
Where a filled pit was then cut through, only the cut shows. The panel title
gives the counts: pits filled, minima used and rejected, and how many steps
of the burned channel rise downstream (lsdtt3 also checks this in its log;
it should be 0).

Comparison panels are the main use. Give one `--csv` per run to compare
`vector_corridor_half_width_metres` values, `burn_profile_max_pit_pixels`
values, or a linear against a chi fit. Other options:

- `--x chi` plots against chi; it needs a chi-fit run. Each path's chi in the
  CSV counts from that path's head, and the plot shifts it to be continuous
  where paths join, so a chi fit draws as straight lines between minima.
- `--mark x,y,label` draws a vertical line at the course cell nearest a
  point, in the run's CRS.
- `--max-join-labels` limits how many joining rivers are named.

In Python:

```python
from lsdviztools3.charts.burned_profiles import (
    plot_burned_profiles, read_burned_profiles, select_course)

plot_burned_profiles(["site_x0_burned_profiles.csv", "site_x60_burned_profiles.csv"],
                     "profiles.png", labels=["X = 0 m", "X = 60 m"], name="River Name")
course, label = select_course(read_burned_profiles("site_x0_burned_profiles.csv"),
                              path_rank=0)   # the rows, downstream, with course_distance_m
```

Files from older `burn_to_corridor` runs, with `z_filled` and `z_burned`
columns, are refused: that design was replaced.

## Smoothed channel profiles and slopes (`lsdviz3-channel-slopes`)

`lsdtt-channel-extraction` with `write_channel_slope_points: true` writes
`<write_fname>_channel_slope_points.<csv|fgb|gpkg>` (format from
`channel_slope_points_write_vector_format`, else `write_vector_format`). It
has one point per channel cell, with these columns:

| Column | Meaning |
|---|---|
| `drainage_area_m2`, `flow_distance_m` | drainage area; flow distance to the outlet |
| `elevation_m` | the pit-filled DEM |
| `slope_m_per_m` | least-squares slope of elevation against distance along the channel, over ± `slope_window_half_length_metres` |
| `smoothed_elevation_m` | that fitted line's elevation at the cell |
| `slope_valid`, `slope_window_length_m` | whether the fit was usable; the along-channel span it used |

Upstream, the window follows the donor with the largest drainage area. A
fit is unusable where the window covers fewer than 3 cells or spans less
than the half-length, which happens near sources. Unusable fits are NaN and
are left out of the slope plots.

```bash
pixi run lsdviz3-channel-slopes \
  --slope-points runs/site_channel_slope_points.csv \
  --dem runs/site.tif --min-area-m2 1e6 --out-dir figures
```

This writes:

- `<stem>_channel_profiles.png`: elevation against flow distance to the
  outlet. The pit-filled DEM is grey; the smoothed profile is blue.
- `<stem>_slope_distance.png`: smoothed slope against flow distance, with a
  log slope axis.
- `<stem>_slope_area.png`: smoothed slope against drainage area, log-log.
- `<stem>_channel_slope_map.png` (with `--dem`): channel cells coloured by
  log₁₀ slope over a hillshade.

The colour range is set by `--slope-vmin` / `--slope-vmax`, in m/m. Values
outside it show the end colours, with arrows on the colourbar.

A flat in the DEM appears as a band of slopes at the gradient that pit
filling imposes on flats (`fill_min_slope`, 0.0001 by default). Pit filling makes those
slopes, not the landscape. The slope measures elevation that pit filling has
raised in places, so treat low-gradient reaches with that in mind.

`--min-area-m2` keeps only cells that drain at least that area.

`--field-points` overlays your own points. The file can be csv with
`latitude`/`longitude`, fgb or gpkg. It needs `drainage_area_m2` and
`slope_m_per_m` columns. For the profiles it also needs `flow_distance_m`,
plus `smoothed_elevation_m` or `elevation_m`.

`--field-color-by <column>` colours the field points by one of their
attributes, such as a grain size, and adds a colourbar. Without it they are
a single colour. On the map they are open circles.

In Python:

```python
from lsdviztools3.charts.channel_slopes import (
    load_channel_points, plot_channel_profiles, plot_slope_area,
    plot_slope_distance, render_channel_slope_map)

pts = load_channel_points("site_channel_slope_points.csv")
fig, ax = plot_slope_area(pts, min_area_m2=1e6)
render_channel_slope_map("site.tif", pts, "slope_map.png", slope_vmin=1e-3, slope_vmax=0.3)
```

## Worked examples
### Channel-bearing comparison (dendritic vs trellis)

`examples/channel_bearing_comparison.py` is an end-to-end script that compares
drainage-network *orientation* between two sites: an Allegheny Plateau tile
(dendritic) and a Valley and Ridge tile (trellis, strongly structurally
aligned). For each site it downloads a COP30 tile, fills and hillshades it
(`lsdtt-raster-preprocessing` with `write_fill_raster` / `write_hillshade_raster`),
runs `lsdtt-channel-network-analysis` with `only_complete_basins: true` for
channel-link bearings, junction angles and outlet→centroid bearings in complete
basins, and `lsdtt-channel-extraction` for the plottable network lines. It
renders a hillshade + channel-network map per site and draws:

- a rose of channel bearings by Strahler stream order ([Rose diagrams (bearings)](#rose-diagrams-bearings)),
- a rose of outlet→centroid bearings,
- a junction-angle distribution comparison (`plot_distributions`),
- a circular comparison of the two sites' channel bearings
  ([Circular / directional comparison](#circular--directional-comparison)).

**The statistics of the comparison come from lsdtt3.** The script downloads one
more DEM covering both sites, writes each site's box as a region polygon
(`site_regions.gpkg`), and runs `lsdtt-channel-network-analysis` once with
`write_channel_bearing_statistics: true`, `channel_bearing_group_by: region`
and `channel_region_fname` pointing at those polygons. lsdtt3 writes the
circular Wasserstein distance, the Kuiper test and a block-permutation p-value
between the two regions to `both_sites_COP30_channel_bearing_comparisons.csv`;
the plot draws the links of that same run
(`_channel_bearing_group_membership.csv`) and annotates those numbers. Links
outside both boxes belong to no region and take no part.

Read the block counts with the p-values: links in one valley are not
independent, so the iid Kuiper p-value is far smaller than the block-permutation
one when there are few independence blocks (`n_blocks_a`, `n_blocks_b`).

Outputs land under `../data/channel_bearing_comparison` (one subdirectory per
site, plus a `comparison/` dir; set `CHANNEL_BEARING_DIR` to write elsewhere).

```bash
pixi run python examples/channel_bearing_comparison.py
pixi run python examples/channel_bearing_comparison.py --skip-compute  # replot only
```

### Basin geometry vs a fault (Parkfield, San Andreas)

`examples/fault_basin_deformation.py` compares basin geometry to a baseline
structure — here the San Andreas Fault. It downloads a COP30 tile around
Parkfield, fills + hillshades it, clips the global GEM active-fault database to
the footprint, then extracts the basins whose channels cross the fault using
`lsdtt-basin-extraction`'s `select_basins_from_line_crossings` mode. **All the
analysis is done in lsdtt3 (C++)** — the tool reprojects the fault internally and
writes `..._line_crossing_basins.csv` with each basin's outlet, centre of mass,
area, outlet→COM bearing, local fault azimuth, side of the fault, and bearing
relative to the fault-orthogonal. The script is purely orchestration + plotting:
it reads that CSV (and the basin polygons) and draws

- the hillshade with the fault trace, with basins, and with outlet→COM vectors;
- a rose of outlet→COM bearings ([Rose diagrams (bearings)](#rose-diagrams-bearings)) with the
  fault direction marked;
- the deflection-from-fault-orthogonal distribution split by side of the fault
  ([Circular / directional comparison](#circular--directional-comparison) / `plot_distributions`), to test
  for shear-related bending.

This illustrates the intended division of labour: **lsdtt3 computes, lsdviztools3
plots**. See the lsdtt3 basin-extraction reference for the
`select_basins_from_line_crossings` parameters and the CSV columns. Outputs land
under `../data/parkfield`.

```bash
pixi run python examples/fault_basin_deformation.py
pixi run python examples/fault_basin_deformation.py --skip-compute  # replot only
```

### Basins draining into a single channel (Nujiang)

`examples/single_channel_basins.py` demonstrates lsdtt3's
`extract_basins_from_single_channel` mode. From an upstream point (and an optional
downstream point) lsdtt-basin-extraction builds **one main channel** by following
D8 — the upstream point is identified by D8 reachability and the channel is
truncated at the common-downstream point where the downstream point's flow joins
— then extracts the side-tributary basins draining into that channel (complete +
minimum-size filtered). All analysis is in lsdtt3; the script preprocesses the
DEM, runs the mode, and plots the single channel (blue) with the side basins
(orange), outlets, and the input points. Test site: a stretch of the Nujiang
(Nu/Salween) gorge; outputs under `../data/nujang`. See the lsdtt3 basin-extraction
reference for the mode's parameters.

```bash
pixi run python examples/single_channel_basins.py
pixi run python examples/single_channel_basins.py --skip-compute  # replot only
```
