# Getting started

This explains how lsdviztools3 draws maps and charts (coordinates, styles, saving), how to
get a DEM, and how to run lsdtt3 from Python. Every plot page assumes these
conventions.

Run Python via the project environment (`pixi run python …`) so PyGMT/GDAL are
on the path. In a notebook started from the pixi env, use `fig.show()` instead
of `save_figure(...)`.
## Conventions

### Coordinate modes (`coordinates=`)

Every map function and CLI takes a coordinate mode:

- `"local"` (default) — reproject to the estimated UTM zone and draw with a
  linear projection; axes are shown in **km** and a `CRS is EPSG:XXXXX` title is
  added. Correct scale bars, undistorted local maps.
- `"geographic"` — EPSG:4326, Mercator; axes in degrees.
- an explicit CRS string (e.g. `"EPSG:3857"`).

Vectors are reprojected to match the base raster automatically.

### `MapStyle` — the styling options shared by all maps

```python
from lsdviztools3.render.style import MapStyle

style = MapStyle(
    cmap="gray",            # raster colour map (or a channel/point CPT)
    figure_size="12c",      # width; MapStyle.SINGLE_COLUMN=8.5c, DOUBLE_COLUMN/MAX_WIDTH=17c
    line_width="1p",        # vector line width
    line_color="black",     # vector colour when not colour-by
    point_size="0.2c",      # default point size when not size-by
    projection=None,        # None => auto projection for the CRS
    frame=None,             # None => auto publication frame
    font_annotation="9p",   # tick labels (>= 8pt)
    font_label="11p",       # axis labels
    font_title="11p",       # CRS title
    dpi=300, high_res_dpi=600,
    transparent=True,       # transparent background + nodata; set False for opaque
)
```

The CLIs expose these as flags: `--cmap`, `--figure-size`, `--line-width`,
`--line-color`, `--projection`, `--font-annotation`, `--font-label`,
`--font-title`, `--dpi`, `--high-res`, `--opaque`, `--coordinates`.

Colourmaps can use both `matplotlib` and `pyGMT` colourmaps.  

### Saving

```python
from lsdviztools3.render.base import save_figure
save_figure(fig, "out.png", style=style, high_res=False)  # 300 dpi (600 if high_res)
```

### Target directory (`--dir`)

The map, animation and download commands take an optional `--dir DIR`: a single target directory for that run's
outputs, **created (with parents) if it does not exist**. A relative `--out` is
written inside it; an absolute `--out` is left where it points (but `DIR` is still
created). For `lsdviz3-lsdtt`, `--dir` is also the default `--data-path` and
`--cwd`, so lsdtt3 reads and writes there.

```bash
# download, run lsdtt3, and plot all into ./tucuy (made on first use)
pixi run python -m lsdviztools3.cli.grab_dem --bbox -71.5 -13.6 -71.2 -13.3 \
  --out rt.tif --dir tucuy
pixi run python -m lsdviztools3.cli.plot_hillshade --in rt.tif --out hs.png --dir tucuy
```

In Python the same is `download_dem(..., directory="tucuy")`,
`run_lsdtt(..., directory="tucuy")`, or `lsdviztools3.paths.resolve_in_dir`.

## Getting data: download a DEM

**NOTE** The `lsdtt3` command line tools have a number of download options using [`lsdtt-fetch-raster`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_fetch_raster.html) and the preferred way to get toopographic data is to use that. 

Keyless Copernicus GLO-30 via the Microsoft Planetary Computer (no credentials).

**CLI**

```bash
pixi run python -m lsdviztools3.cli.grab_dem \
  --bbox -6.20 56.85 -6.06 56.95 \
  --source planetary_computer \
  --out demo_output/eigg_glo30.tif
```

`--bbox WEST SOUTH EAST NORTH` (WGS84 degrees), `--source
{planetary_computer,opentopography}`, `--dataset` (optional), `--out`, `--dir`
(target directory, see [Conventions](#conventions)), `--ot-api-key-file` (OpenTopography
only — the key is never a direct argument).

**Python**

```python
from lsdviztools3.download.ot_scraper import download_dem

download_dem((-6.20, 56.85, -6.06, 56.95), source="planetary_computer",
             out_file="demo_output/eigg_glo30.tif")
```

Named demo regions are available too:

```python
from lsdviztools3.regions import REGIONS
download_dem(REGIONS["guadalcanal"].bbox, source="planetary_computer",
             out_file="guadalcanal.tif")
```

To choose your own box with Google Maps, right-click the lower left
(south-west) corner of the area and click the numbers at the top of the menu
to copy them. Do the same for the upper right (north-east) corner, and paste
both:

```python
from lsdviztools3.regions import box_from_corners, parse_lat_lon

box = box_from_corners("55.90, -3.25", "55.99, -3.10")   # prints the box and its size
west, east, south, north = box                           # or box.west, box.north, ...
download_dem(box.bbox, source="planetary_computer", out_file="edinburgh.tif")

latitude, longitude = parse_lat_lon("55.9476273103734, -3.158220317286066")
```

Either pair of opposite corners gives the same box. `box.bbox` is the
`(west, south, east, north)` order that `download_dem` and `Region.bbox` use.
Each corner must be decimal degrees with the latitude first, which is what
Google Maps copies. Degrees-minutes-seconds text is rejected with a message
saying how to copy the decimal form.

## Running lsdtt3

The map-from-raster and download paths need no lsdtt3. Channel/chi/swath
products come from lsdtt3 tools, driven via `lsdviztools3.driver` /
`lsdviztools3.workflows`. Point lsdviztools3 at your lsdtt3 `bin`:

```bash
export LSDTT3_DIR=/path/to/lsdtt3/bin
# or persist it:
pixi run python -m lsdviztools3.cli.run_lsdtt --set-lsdtt-dir /path/to/lsdtt3/bin
```

Run any tool with a parameter dict:

```python
from lsdviztools3.driver import run_lsdtt

run_lsdtt("lsdtt-channel-extraction", {
    "read_path": "demo_output", "read_fname": "eigg_glo30.tif",
    "write_path": "demo_output", "write_fname": "eigg",
    "threshold_contributing_pixels": 300,
    "write_channel_network_lines": "true",
    "write_channel_network_points": "true",   # stream_order + drainage_area
})
```

See [Python workflows](#python-workflows) for the high-level workflows that chain
preprocessing → analysis → plot.
## Python workflows

`lsdviztools3.workflows` chains lsdtt3 + rendering. Each returns the output
paths (or a Figure/PNG path).

| Workflow | What it does |
|----------|--------------|
| `lsdtt3_hillshade_channels_basin(dem, workdir, stem=…, source_threshold=…)` | preprocess → channel + largest-basin extraction; returns hillshade/channels/basin paths |
| `render_hillshade_channels_basin(outputs, png)` | composite hillshade + channels + largest basin |
| `lsdtt3_chi_analysis(dem, workdir, stem=…, source_threshold=…, n_largest_basins=5, m_over_n=0.45)` | preprocess → chi-analysis (fit), N largest complete basins; returns fit CSV/fgb paths |
| `render_chi_profiles(outputs, out_dir)` | stacked chi plot + one profile per basin |
| `render_chi_ksn_map(dem, segment_lines, png)` | hillshade + log₁₀ k_sn network + colourbar + scale bar |
| `render_chi_ksn_relief3d(dem, segment_lines, png, *, vertical_exaggeration=2.0)` | 3D hillshaded relief with the chi network draped on the surface, coloured by log₁₀ k_sn + colourbar |
| `lsdtt3_swath_profile(dem, workdir, stem=…, half_width_metres=…, bin_width_metres=…)` | preprocess → main-channel baseline → swath; returns swath CSV + rasters |
| `render_swath_profile(outputs, png)` | swath envelope plot |
| `lsdtt3_stream_order_basins(dem, workdir, stem=…, stream_order=…)` | extract all junction basins of one Strahler order; returns polygon paths |
| `render_hillshade_stream_order_basins(outputs, png)` | hillshade + channels + every selected-order basin polygon |
| `render_cosmo_comparison_figures(scan_csv, out_dir, *, sample_comparison_csv=None)` | sweep residual panels (+ optional sample residual boxplots) from standard lsdtt3 cosmo CSVs |
| `render_cosmo_published_residuals(published_dir, out_dir, *, datasets, mode=...)` | residual boxplots directly from lsdtt3 published summary tables |
| `join_basins_with_cosmo_erates(basins, cosmo_csv, ...)` | spatially join point erosion rates to basin polygons and aggregate to per-basin values |
| `render_cosmo_basin_erate_map(base, basins, cosmo_csv, out_png, ...)` | canonical hillshade + semi-transparent basin overlay colored by inferred erosion rate with scaling/muon/atmosphere metadata |
| `summarize_basin_ksn(fit_nodes_csv, fit_segments_csv, ...)` | per-basin `k_sn` summary from chi nodes: median + 25th/75th percentiles + asymmetric IQR errors |
| `render_chi_basin_ksn_map(base, basins, fit_nodes_csv, fit_segments_csv, out_png, ...)` | hillshade + basin polygons colored by median `k_sn` |
| `render_chi_and_cosmo_basin_workflow(outputs, out_dir, *, basin_polygons=..., base_raster=..., cosmo_csv=...)` | writes one chi profile PNG per basin, basin `k_sn` summary CSV/map, and optional `k_sn` vs inferred erosion-rate errorbar plot (`charts.plot_ksn_vs_erate`) |
| `render_hand_profile(…)` / `render_hand_window_populations(…)` | the [Valley floor / HAND charts](valleys_and_cross_sections.md#valley-floor--hand-charts) HAND charts, straight from lsdtt3 output paths |
| `hillshade_with_vectors(dem, layers, png)` | hillshade with any number of vector layers overlaid ([Overlaying several vector layers on a hillshade](basic_plots.md#overlaying-several-vector-layers-on-a-hillshade)) |
| `compare_stream_networks(…)` | stream-burning preset over `hillshade_with_vectors` |
| `plot_segment_classification(…)` / `plot_algoviz(…)` | stream-burning diagnostics: trace-divergence per segment; replay of the `snap_to_accumulation` source-path burn |
| `sample_dem_along_line(dem, line)` / `trace_d8_flowpath(dem, start)` | sample elevation along a polyline; trace a D8 steepest-descent path |
| `plot_long_profile(profiles, …)` | elevation-vs-distance long profiles on shared axes |

`lsdtt3_chi_analysis` selects the **N largest complete basins** and moves
interior (coastal) outlets upstream so island basins pass the completeness test
— see the switches `n_largest_basins` and
`n_pixels_move_outlet_upstream` in the
[lsdtt3 chi analysis reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_chi_analysis.html).

## Example scripts
| Script | What it does |
|---|---|
| `examples/island_demo.py` | Keyless DEM download + hillshade (`pixi run demo-eigg`, `demo-rum`, `demo-fiji`, `demo-solomon`). |
| `examples/lsdtt3_workflow.py` | Download → lsdtt3 hillshade/channels/largest basin → composite map (`pixi run demo-lsdtt3`). |
| `examples/cosmo_comparison_workflow.py` | Cosmogenic sweep comparison figures ([Cosmo comparison workflows](cosmogenic_nuclides.md#cosmo-comparison-workflows)). |
| `examples/synthetic_networks_demo.py` | Synthetic OCN networks: generate, render and compare junction geometry ([Channel networks](channel_networks.md#network-geometry-and-synthetic-networks)). |
| `examples/hillslope_comparison.py` | Newton vs cellular hillslope solver comparison ([Landscape evolution](landscape_evolution.md#hillslope-solver-comparisons)). |
| `examples/perron_ren_p1_p2.py` … `perron_ren_p6.py` | Perron/Ren hillslope benchmark figures: P1–P2 outputs, P3 timestep independence, P4 transient convergence and runtimes, P5 stability envelope, P6 spatial convergence and anisotropy. |
