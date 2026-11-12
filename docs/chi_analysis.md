# Chi analysis

Maps and profiles of an [`lsdtt-chi-analysis`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_chi_analysis.html) run: channel steepness (k_sn) on
a map, stacked chi–elevation plots of many basins, and single-basin profiles.
Nothing is computed here: chi, the piecewise fit and k_sn all come from
lsdtt3.

## The `lsdviz3-chi` command

For a run made with `fit_chi_networks: true` (and `write_fit_tables`,
`write_chi_fit_segment_lines`), one command draws the k_sn map, the stacked
chi plot, and optionally one profile per basin:

```bash
pixi run lsdviz3-chi --data-dir work/eigg --prefix eigg --m-over-n 0.45 \
  --dem work/eigg/eigg_remove_seas.tif --out-dir demo_output/chi --per-basin
```

| Flag | Meaning |
|---|---|
| `--data-dir`, `--prefix` | The run's directory and `write_fname`. |
| `--m-over-n` | m/n of the run, `0.45` or the file suffix `0p450`. |
| `--dem` | Projected DEM for the map's hillshade. Without it no map is drawn. |
| `--out-dir`, `--out-prefix` | Where the PNGs go, and their prefix (default `--prefix`). |
| `--per-basin` | Also write one chi–elevation profile per basin. |
| `--no-stacked` | Skip the stacked plot. |
| `--cmap`, `--coordinates`, `--dpi` | k_sn colour map (default `viridis`), map coordinates, profile DPI. |
| `--ksn-vmin`, `--ksn-vmax` | The k_sn colour range, **in k_sn units** (see below). Default: the data's range. |
| `--ksn-classes`, `--n-classes`, `--class-colors`, `--class-color-mode` | Colour the map's k_sn in classes, with a class legend instead of the colourbar (see [Colouring k_sn in classes](#colouring-k_sn-in-classes)). |

It writes `<prefix>_ksn_map.png`, `<prefix>_chi_stacked.png` and
`<prefix>_chi_profile_basin<N>.png`. The sections below describe the plots
and the Python functions behind them.

### Setting the k_sn colour range

k_sn is coloured on a log₁₀ scale, so a few very steep links (a channel
routed over a dam, a DEM artefact) stretch the scale and wash out the rest of
the network. `--ksn-vmin` and `--ksn-vmax` set the range in **k_sn**, not in
log₁₀ k_sn. `--ksn-vmax 400` puts the top of the scale at log₁₀ 400 = 2.6.

- Links outside the range are drawn in the end colours, not dropped.
- The colourbar shows an arrow (triangle) at each end that data run past.
- The same range is used for the map, the stacked plot and every per-basin
  profile, so per-basin profiles share one colour scale and can be compared.
- Both limits must be above 0.
- Check the spread of the run's k_sn first. The fit-segments table
  `*_chi_fit_segments_mn_*.csv` has one `k_sn` per segment; its 95th or
  99th percentile is a reasonable `--ksn-vmax`.

```bash
pixi run lsdviz3-chi --data-dir work/eigg --prefix eigg --m-over-n 0.45 \
  --dem work/eigg/eigg_remove_seas.tif --out-dir demo_output/chi \
  --ksn-vmin 10 --ksn-vmax 400
```

In Python the same limits are `vmin=` / `vmax=` on `plot_chi_profile`,
`plot_chi_stacked` and `render_channels`. On the workflows
`render_chi_ksn_map`, `render_chi_ksn_relief3d` and `render_chi_profiles`
they are `ksn_vmin=` / `ksn_vmax=`. All of them take k_sn units.

### Colouring k_sn in classes

`--ksn-classes` colours the k_sn map in classes, with a legend of the classes
instead of the colourbar: give the breaks in k_sn units (`--ksn-classes
20,50,100,200` makes five classes), or a scheme, `equal_interval`, `quantile`
or `natural_breaks` (`jenks`), with `--n-classes`. The schemes work in log10
k_sn, so `equal_interval` gives classes of equal ratio. `--ksn-vmin` and
`--ksn-vmax` set the ends of the colour gradient (each class takes the colour
at its centre) and the range the schemes use. `--class-colors` gives one
colour per class, and a qualitative `--cmap` (`tab10`, `Set1`, ...) gives
distinct colours. The same options work on the map of `lsdviz3-knickpoints`.
The details are in [Colouring in classes](basic_plots.md#colouring-in-classes).

```bash
pixi run lsdviz3-chi --data-dir . --prefix salmon_ksn --m-over-n 0.45 \
  --dem salmon_DEM.tif --out-dir . --ksn-classes natural_breaks --n-classes 5
```

### Removing tributaries

A run made with a current lsdtt3 also writes `<prefix>_chi_tributaries_mn_*.csv`
(`write_chi_tributary_table`, on by default). It has one row per tributary of
the fit: a branch from its head down to where it joins its parent. Each row
gives:
- its `source_key` (the chi data maps' channel key, also on every fit
  segment);
- the drainage area just before it joins (`junction_drainage_area_m2`), and
  that area as a fraction of the basin;
- its rank by that area within the basin;
- its parent, and how many tributaries drain into it;
- its steepest segment's k_sn (`ksn_max`).

`lsdviz3-chi` uses it to leave tributaries out of the plots. Main stems are
never removed:

| Flag | Removes a tributary when |
|---|---|
| `--remove-ksn-above VALUE` | any of its segments is steeper than VALUE (k_sn). `q3_plus_3iqr` uses the outlier value lsdtt3 suggests in `<prefix>_chi_knickpoint_thresholds_mn_*.csv`. For links routed over a dam onto a steep face. |
| `--min-tributary-area-m2 A` | its drainage area where it joins is below A. |
| `--min-tributary-area-fraction F` | that area is below F of its basin's. |
| `--top-tributaries N` | it is not among the N largest tributaries of its basin. |
| `--keep-descendants` | (turns off) removing the tributaries that drain into a removed one, which is the default. |
| `--show-removed ghost\|hide` | how removed tributaries appear on the profiles: thin and light grey (default), or not at all. |

```bash
pixi run lsdviz3-chi --data-dir work/eigg --prefix eigg --m-over-n 0.45 \
  --dem work/eigg/eigg_remove_seas.tif --out-dir demo_output/chi \
  --remove-ksn-above q3_plus_3iqr --min-tributary-area-fraction 0.01
```

With any of these flags the command also writes:
- `<prefix>_removed_tributaries.csv`: the removed rows of the tributary
  table, with a `reason` column (`ksn_outlier`, `small_area`,
  `small_fraction`, `rank` or `drains_to_removed`; the first rule that
  applied, in that order);
- with `--dem`, `<prefix>_removed_tributaries_map.png`: the kept network
  coloured by k_sn, each removed tributary dashed in its reason's style,
  junctions of removed tributaries as open circles, and a legend with counts.

The k_sn map and the profiles then show only the kept tributaries. Their
colour range is taken from the kept ones, so removing a few dam-routed links
also stops them stretching the scale. A steep step on a main stem (a dam on the
trunk) is not removed: use `--ksn-vmax` for that.

In Python:

```python
from lsdviztools3.io.chi import load_chi_fit_segments, load_chi_tributaries, select_tributaries
from lsdviztools3.charts import plot_chi_stacked

segments = load_chi_fit_segments("eigg_chi_fit_segments_mn_0p450.csv")
kept, removed = select_tributaries(load_chi_tributaries("eigg_chi_tributaries_mn_0p450.csv"),
                                   ksn_max=800, min_area_fraction=0.01)
fig, ax = plot_chi_stacked(segments, removed=removed, show_removed="ghost")
```

`removed=` is accepted by `plot_chi_profile`, `plot_chi_stacked`,
`render_chi_ksn_map` and `render_chi_profiles`, and `render_chi_removed_map`
draws the map. Selection only chooses what is drawn: the chi fit is not
redone without the removed tributaries.
## Knickpoints (`lsdviz3-knickpoints`)

A current lsdtt3 chi run also writes `<prefix>_chi_knickpoints_mn_*.csv`
(`write_chi_knickpoint_table`, on by default). It has one row per node where
two fitted k_sn segments meet, with:
- `delta_ksn`: upstream minus downstream k_sn, positive where the channel
  steepens upstream;
- `log_ratio_ksn`: log₁₀ of their ratio;
- `at_junction`: 1 where a tributary's first segment meets its parent, at
  its mouth.

It also writes `<prefix>_chi_knickpoint_thresholds_mn_*.csv`, suggested
thresholds for a "major" knickpoint over the within-tributary rows: for
|Δk_sn| and |log ratio|, the 90th and 95th percentiles, the median plus three
times the MAD (scaled by 1.4826) and Otsu's two-class split. lsdtt3 computes
all of these; this command draws them and applies the threshold you choose.

```bash
pixi run lsdviz3-knickpoints --data-dir work/eigg --prefix eigg --m-over-n 0.45 \
  --out-dir demo_output/knickpoints --dem work/eigg/eigg_remove_seas.tif \
  --threshold otsu
```

| Flag | Meaning |
|---|---|
| `--measure delta\|log_ratio` | The change in k_sn: the difference (default) or the log ratio. |
| `--threshold` | A value of \|measure\|, or a suggestion from the thresholds table: `p90`, `p95` (default), `median_plus_3mad`, `otsu`. |
| `--direction both\|steepening\|shallowing` | Count changes in both directions (default), or only where the channel steepens, or shallows, upstream. |
| `--include-junctions` | Also count changes at tributary mouths (left out by default). |
| `--along chi\|chi_normalised\|elevation` | x axis of the along-chi plot. |
| `--basin N` | One basin's profile instead of the stacked plot. |
| `--ksn-vmin`, `--ksn-vmax`, `--cmap` | The k_sn colours of the profile and map, as `lsdviz3-chi`. |
| `--ksn-classes`, `--n-classes`, `--class-colors`, `--class-color-mode` | The map's k_sn network in classes, as `lsdviz3-chi`. |
| tributary flags | `--remove-ksn-above`, `--min-tributary-area-m2`, `--min-tributary-area-fraction`, `--top-tributaries`, `--keep-descendants`, `--show-removed` as in [Removing tributaries](#removing-tributaries). The knickpoints of removed tributaries are dropped. |

It writes:

- `<prefix>_knickpoint_distribution.png`. Left: a histogram of the signed
  change, with counts on a log axis so the tails show. Right: the absolute
  changes ranked largest first, where a break in slope can be picked by eye.
  Both carry every suggested threshold (with how many knickpoints lie
  beyond it) and the chosen one in black. To choose by hand, read a value
  off the ranked curve and pass it as `--threshold`.
- `<prefix>_knickpoint_profile.png`: the chi profile (stacked, or one
  basin) with the major knickpoints. An open upward triangle marks where the
  channel steepens upstream, a filled downward one where it shallows; the
  size scales with the change.
- `<prefix>_knickpoints_along_chi.png`: a histogram of the major
  knickpoints along chi over that of all segment boundaries, and one row per
  basin. Knickpoints at the same chi across tributaries or basins are what a
  wave of incision moving upstream would leave.
- `<prefix>_knickpoint_map.png` (with `--dem`): the same symbols on the
  hillshade over the k_sn network.
- `<prefix>_major_knickpoints.csv`: the selected rows, with the measure,
  threshold and direction used.

A steep segment between two gentle ones gives a pair: a steepening knickpoint
at its foot and a shallowing one at its top. A staircase of such pairs along
one channel often marks a single knickzone fitted in pieces.

In Python: `select_major_knickpoints`, `plot_knickpoint_distribution`,
`plot_knickpoints_on_profile`, `plot_knickpoints_along_chi` and
`render_knickpoint_map` in `lsdviztools3.charts.knickpoints`;
`load_chi_knickpoints` and `load_chi_knickpoint_thresholds` in
`lsdviztools3.io.chi`.

## Channel network coloured by log₁₀ k_sn

The chi network with per-segment `k_sn` is in
`*_chi_fit_segment_lines_mn_*.fgb` (from `lsdtt-chi-analysis`; see
[Python workflows](getting_started.md#python-workflows)). `render_channels(color_by="k_sn")` draws each
segment through a CPT (log₁₀ by default) and can add a colourbar legend.

The easiest route is the workflow helper, which composes hillshade + network +
colourbar + scale bar:

```python
from lsdviztools3.workflows import lsdtt3_chi_analysis, render_chi_ksn_map

out = lsdtt3_chi_analysis("demo_output/eigg_glo30.tif", "work/eigg", stem="eigg",
                          source_threshold=300, n_largest_basins=5, m_over_n=0.45)
render_chi_ksn_map("work/eigg/eigg_remove_seas.tif", out["segment_lines"],
                   "demo_output/eigg_ksn_map.png", cmap="viridis")
```

Or compose it by hand to control each layer:

```python
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.render.channels import render_channels
from lsdviztools3.render.decorations import add_scalebar
from lsdviztools3.render.base import save_figure, grid_to_display
from lsdviztools3.io.raster import load_raster, to_coordinates, to_gmt_grid, get_crs
from lsdviztools3.render.style import MapStyle

dem = "work/eigg/eigg_remove_seas.tif"
da = to_coordinates(load_raster(dem), "local"); crs = get_crs(da)
_, extent = grid_to_display(to_gmt_grid(da), crs)   # km extent for the scale bar

fig = render_hillshade(dem, coordinates="local", shade=True,
                       style=MapStyle(cmap="gray", figure_size="17c", transparent=False))
fig = render_channels(
    out["segment_lines"],
    color_by="k_sn", log_color=True, cmap="viridis",
    colorbar=True, colorbar_label="log@-10@-k@-sn@-",
    coordinates="local",
    style=MapStyle(line_width="1.2p", cmap="viridis", figure_size="17c"),
    fig=fig,
)
add_scalebar(fig, extent, crs)
save_figure(fig, "demo_output/eigg_ksn_map.png", style=MapStyle(transparent=False))
```

`render_channels(lines, *, color_by=None, log_color=True, cmap=None,
colorbar=False, colorbar_label=None, coordinates="local", style=None, fig=None)`.
With `color_by=None` it draws plain lines in `style.line_color`. You can colour
by any numeric attribute (e.g. `stream_order` on
`*_channel_network_lines.fgb`).
## Chi profiles

Two views: a single basin's chi–elevation profile
([Single-basin profiles (`plot_chi_profile`)](#single-basin-profiles-plot_chi_profile)), and the lsdviztools
**stacked chi plot** — every basin normalised to χ=0, then offset along the χ
axis so the profiles sit side-by-side on one pair of axes, coloured by
log₁₀ k_sn ([Stacked chi plots (`plot_chi_stacked`)](#stacked-chi-plots-plot_chi_stacked)).

Quick start:

```python
from lsdviztools3.io.chi import load_chi_fit
from lsdviztools3.charts import plot_chi_profile, plot_chi_stacked

# Attach per-node k_sn from the segment table (one step):
nodes = load_chi_fit("work/eigg/eigg_chi_fit_nodes_mn_0p450.csv",
                     "work/eigg/eigg_chi_fit_segments_mn_0p450.csv")

fig, ax = plot_chi_stacked(nodes, color_by="k_sn", log_color=True, cmap="viridis")
fig.savefig("demo_output/chi_stacked.png", dpi=300, bbox_inches="tight")

fig, ax = plot_chi_profile(nodes, basin=0)        # one basin
fig.savefig("demo_output/chi_basin0.png", dpi=300, bbox_inches="tight")
```

Or run the analysis and write a stacked plot plus one profile per basin in one
call:

```python
from lsdviztools3.workflows import lsdtt3_chi_analysis, render_chi_profiles
out = lsdtt3_chi_analysis("demo_output/eigg_glo30.tif", "work/eigg", stem="eigg",
                          source_threshold=300, n_largest_basins=5, m_over_n=0.45)
render_chi_profiles(out, "demo_output/chi", stem="eigg")
# -> demo_output/chi/eigg_stacked.png + eigg_profile_basin<N>.png
```

### Inputs: which lsdtt3 call produces which file

Every chi chart consumes output from **`lsdtt-chi-analysis`**. Nothing is
computed here: chi, the piecewise-linear fit and k_sn all come from the C++
side. The minimal parameter set — exactly what
`workflows.lsdtt3_chi_analysis` sends — is:

```python
from lsdviztools3.driver import run_lsdtt

run_lsdtt("lsdtt-chi-analysis", {
    "read_path": "work/eigg", "read_fname": "eigg_remove_seas.tif",
    "write_path": "work/eigg", "write_fname": "eigg",
    "threshold_contributing_pixels": 300,
    "n_largest_basins": 5,            # keep >0: see the warning in 5.1.4
    "write_chi_data_maps": "true",    # -> observed points
    "chi_data_maps_write_vector_format": "flatgeobuf",
    "fit_chi_networks": "true",       # -> the fit itself
    "write_fit_tables": "true",       # -> nodes + segments CSVs
    "write_chi_fit_segment_lines": "true",
    "write_vector_format": "flatgeobuf",
    "mn_min": 0.45, "mn_max": 0.45, "mn_step": 0.02,
}, data_path="work/eigg")
```

`mn_*` (or `m_over_n`) sets the concavity; the value appears in the filename as
`mn_0p450`. The four products and their columns:

| File | Written when | Loader | Columns |
|---|---|---|---|
| `*_chi_fit_nodes_mn_*.csv` | `write_fit_tables` | `load_chi_fit_nodes` | `basin_index`, `outlet_flat`, `flat`, `receiver_flat`, `chi`, `observed_elevation`, `fitted_elevation`, `residual_elevation`, `abs_residual_elevation`, `branch_id`, `is_main_stem`, `is_source`, `is_junction`, `is_split_point` — **no `k_sn`** |
| `*_chi_fit_segments_mn_*.csv` | `write_fit_tables` | `load_chi_fit_segments` | `basin_index`, `outlet_flat`, `segment_id`, `parent_segment_id`, `branch_id`, `downstream_flat`, `upstream_flat`, `downstream_chi`, `upstream_chi`, `downstream_elevation`, `upstream_elevation`, **`k_sn`**, `max_abs_residual`, `n_nodes` |
| `*_chi_fit_segment_lines_mn_*.fgb` | `write_chi_fit_segment_lines` | `load_chi_segment_lines` | the segment columns *minus* `parent_segment_id` and `max_abs_residual`, plus map-space geometry (dropped on load) |
| `*_chi_data_maps.fgb` / `.csv` | `write_chi_data_maps` | `load_chi_points` | `chi`, `drainage_area`, `flow_distance`, `elevation`, `source_key`, **`basin_key`** |

`k_sn` lives on the **segments**, not the nodes. `load_chi_fit(nodes, segments)`
(or `assign_node_ksn`) grafts it on by finding, within each
`(basin_index, branch_id)`, the segment whose `[downstream_chi, upstream_chi]`
interval contains each node — nodes just outside every interval take the
nearest segment in the branch.

Two more `lsdtt-chi-analysis` products are used elsewhere in this guide, not by
the chi charts: `*_chi_ksn_basin_summary_mn_*.csv` (per-basin k_sn statistics)
and the segment-lines fgb again, as a **map** layer in
[Channel network coloured by log₁₀ k_sn](#channel-network-coloured-by-log10-k_sn).

### Stacked chi plots (`plot_chi_stacked`)

```python
plot_chi_stacked(data, *, basins=None, x_offset=None, color_by="k_sn",
                 log_color=True, cmap="viridis", point_size=3.0, show_fit=True,
                 observed=None, observed_color="0.6", observed_size=1.5,
                 label_basins=True, ax=None, colorbar=True)  ->  (fig, ax)
```

| Parameter | Default | Meaning |
|---|---|---|
| `data` | — | Nodes **or** segments table; see [Nodes table vs segments table](#nodes-table-vs-segments-table). Must carry `basin_index`. |
| `basins` | all present | Ordered `basin_index` values to stack — also the only way to reorder the stack. |
| `x_offset` | 20% of the widest basin's χ span | Horizontal gap between successive basins, in χ units. |
| `color_by` | `"k_sn"` | Column used for colour. |
| `log_color` | `True` | Colour on a log₁₀ scale. See the flooring rule below. |
| `cmap` | `"viridis"` | Any matplotlib colourmap. |
| `point_size` | `3.0` | Marker size — nodes table only. |
| `show_fit` | `True` | Overlay the fitted profile per `branch_id` — nodes table only (a segments table *is* the fit). |
| `observed` | `None` | Chi data-map frame overlaid as neutral points; see [Overlaying observed points](#overlaying-observed-points). |
| `observed_color` / `observed_size` | `"0.6"` / `1.5` | Style of that overlay. |
| `label_basins` | `True` | Annotate each stack with its `basin_index`. |
| `ax` | `None` | Draw into an existing `Axes` instead of a new figure. |
| `colorbar` | `True` | Attach a colourbar. |

Conventions worth knowing before you read a stacked plot:

- **Stacking.** Each basin's χ is shifted so the basin starts at 0, then
  translated by a running cursor: `cursor += basin_chi_span + x_offset`. The x
  axis is therefore *offset χ*, not comparable between stacks — only the shape
  and the slope (k_sn) are.
- **Order.** The default is ascending `basin_index`, i.e. the order lsdtt3
  retained the basins in, **not** by drainage area or k_sn. Pass `basins=[…]`
  to impose your own order.
- **Colour normalisation is global** across the selected basins, so colours are
  comparable between stacks in one figure but not between figures.
- **log₁₀ k_sn flooring.** With `log_color=True`, values `>= 0.1` are
  log₁₀-scaled and anything below (including 0 and negatives) is mapped to
  `0.0` rather than dropped — mirroring lsdviztools. On low-relief basins this
  piles genuinely small k_sn onto the same colour as k_sn = 1. Set
  `log_color=False` if that matters for your site.

### Nodes table vs segments table

Both plotters detect the input type: a frame carrying `downstream_chi` **and**
`upstream_chi` is treated as *segments*, otherwise as *nodes*.

| Input | What is drawn |
|---|---|
| **Nodes** (`load_chi_fit` / `load_chi_fit_nodes`) | Scatter of `observed_elevation` vs `chi`, coloured by the grafted-on `k_sn`, with `fitted_elevation` overlaid as a line per `branch_id` (`show_fit`). |
| **Segments** (`load_chi_fit_segments` or `load_chi_segment_lines`) | A `LineCollection` of the fitted segments, each coloured by its own `k_sn`. No observed scatter (use `observed=`), and no node→segment join needed. |

The segments route lets you plot straight from the FlatGeobuf, with no CSVs:

```python
from lsdviztools3.io.chi import load_chi_segment_lines, load_chi_points
segs = load_chi_segment_lines("work/eigg/eigg_chi_fit_segment_lines_mn_0p450.fgb")
obs  = load_chi_points("work/eigg/eigg_chi_data_maps.fgb")  # observed pts (no k_sn)
fig, ax = plot_chi_stacked(segs, observed=obs)
```

> **`data=` must carry `basin_index`.** Both plotters group and label by
> `basin_index`, so a chi **data-map** frame (which has `basin_key` instead)
> works only as `observed=`, never as `data=`. Passing one as `data=` raises a
> `KeyError` on `basin_index`.

### Overlaying observed points

`observed=` takes a chi data-map frame and scatters its `chi`/`elevation` under
the fit, using the same per-basin offset. Those points carry no `k_sn`, so they
are drawn in a neutral colour (`observed_color`).

The overlay is joined **by basin number**: `basin_index` on the fit side,
`basin_key` on the data-map side. These are two independent counters in lsdtt3
— `basin_key` is assigned in `lsdtt_chi_data_maps_cpu.cpp` over the *complete*
basins, `basin_index` is the position in the fitted-basin list from
`lsdtt_chi_network_analysis_cpu.cpp`.

> **Keep `n_largest_basins > 0`** (the default is 5) whenever you overlay
> observed points. With `n_largest_basins > 0`, lsdtt-chi-analysis filters the
> basin list to complete basins before both products are written and the two
> numberings agree. With `n_largest_basins: 0` the fit keeps incomplete basins
> while the data maps still drop them, the numbering shifts, and the observed
> points are silently attached to the **wrong stack**. If you must use `0`,
> join on `outlet_flat` yourself instead of relying on the basin number.

### Single-basin profiles (`plot_chi_profile`)

```python
plot_chi_profile(data, *, basin=None, color_by="k_sn", log_color=True,
                 cmap="viridis", point_size=4.0, show_fit=True, observed=None,
                 observed_color="0.6", observed_size=2.0, ax=None,
                 colorbar=True)  ->  (fig, ax)
```

Same input rules, same colour conventions; `basin` selects the `basin_index` to
draw (default: the first one present) and the title is set to `Basin <n>`.
There is no offset — the x axis is true χ. Raises `ValueError` if the requested
basin is absent.

## Chi network on a 3D surface
coloured by log₁₀ k_sn (with a colourbar), over a hillshaded 3D relief:

```python
from lsdviztools3.workflows import lsdtt3_chi_analysis, render_chi_ksn_relief3d

out = lsdtt3_chi_analysis("demo_output/eigg_glo30.tif", "work/eigg", stem="eigg",
                          source_threshold=300, n_largest_basins=5, m_over_n=0.45)
render_chi_ksn_relief3d("work/eigg/eigg_remove_seas.tif", out["segment_lines"],
                        "demo_output/eigg_ksn_relief3d.png",
                        azimuth=135, elevation=30, vertical_exaggeration=2.0)
```

`render_chi_ksn_relief3d(dem_path, segment_lines, png_path, *,
coordinates="local", azimuth=135.0, elevation=30.0, vertical_exaggeration=2.0,
cmap="viridis", width=MapStyle.MAX_WIDTH, high_res=False)`. `segment_lines` is
the `*_chi_fit_segment_lines_mn_*.fgb` from `lsdtt3_chi_analysis`.

## Basin k_sn maps

`summarize_basin_ksn` reduces the fit nodes to one k_sn per basin (median with
25th/75th percentiles), and `render_chi_basin_ksn_map` colours basin polygons
by it. `lsdviz3-chi-cosmo` bundles these with the stacked and per-basin chi
plots and, given cosmogenic erosion rates, a k_sn-versus-erosion-rate plot; see
[Cosmogenic nuclides](cosmogenic_nuclides.md#chi-and-cosmogenic-basins).
