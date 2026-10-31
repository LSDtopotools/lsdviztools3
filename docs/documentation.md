# lsdviztools3 — Usage documentation

Welcome to lsdviztools3 documention. `lsdviztools3` is an update to `lsdviztools`, and is designed to produce plots from the outputs of `lsdtt3`, a topographic analysis software package. 

This documentation includes worked examples for different plot types, with the options each one accepts, shown both as **command-line** invocations (headless PNG output) and as **Python /
notebook** code. For installation and a one-line quick start, see the
[home page](index.md) or the [README on GitHub](https://github.com/LSDtopotools/lsdviztools3#readme).

Contents:

1. [Conventions: coordinates, MapStyle, saving](#1-conventions)
2. [Getting data: download a DEM](#2-getting-data-download-a-dem)
3. [Running lsdtt3 (driver + workflows)](#3-running-lsdtt3)
4. [Maps](#4-maps)
   - [Hillshade](#41-hillshade)
   - [Channel network coloured by log₁₀ k_sn](#42-channel-network-coloured-by-log10-k_sn)
   - [Points sized by drainage area / coloured by stream order](#43-points-sized-by-drainage-area)
   - [Basins](#44-basins)
   - [Locator / orthographic basemap](#45-locator--orthographic-basemap)
   - [Decorations: scale bar and colourbar](#46-decorations-scale-bar-and-colourbar)
   - [3D perspective relief (rayshader-style)](#47-3d-perspective-relief-rayshader-style)
   - [Draping a data raster over a precomputed hillshade](#48-draping-a-data-raster-over-a-precomputed-hillshade)
   - [Overlaying several vector layers on a hillshade](#49-overlaying-several-vector-layers-on-a-hillshade)
5. [Charts (x–y plots)](#5-charts-xy-plots)
   - [Chi profiles (individual + stacked)](#51-chi-profiles)
     - [Inputs: which lsdtt3 call produces which file](#511-inputs-which-lsdtt3-call-produces-which-file)
     - [Stacked chi plots](#512-stacked-chi-plots-plot_chi_stacked)
     - [Nodes table vs segments table](#513-nodes-table-vs-segments-table)
     - [Overlaying observed points](#514-overlaying-observed-points)
     - [Single-basin profiles](#515-single-basin-profiles-plot_chi_profile)
   - [Swath profiles](#52-swath-profiles)
   - [Rose diagrams (bearings)](#53-rose-diagrams-bearings)
   - [Circular / directional comparison](#54-circular--directional-comparison)
   - [Distributions](#55-distributions)
   - [Valley floor / HAND charts](#56-valley-floor--hand-charts)
   - [Valley cross sections](#57-valley-cross-sections)
   - [Hydraulic breaks](#58-hydraulic-breaks)
   - [Cosmo comparison workflows](#59-cosmo-comparison-workflows)
   - [Transient cosmogenic and landslide charts](#510-transient-cosmogenic-and-landslide-charts)
   - [Hillslope charts](#511-hillslope-charts)
   - [Fluvial-surface ladder and valley sections](#512-fluvial-surface-ladder-and-valley-sections)
   - [Divide migration charts](#513-divide-migration-charts)
   - [Network geometry, synthetic networks and valley spacing](#514-network-geometry-synthetic-networks-and-valley-spacing)
6. [End-to-end workflows](#6-end-to-end-workflows)
7. [CLI reference](#7-cli-reference)
8. [lsdtt3 output → plot function map](#8-lsdtt3-output--plot-function-map)

Run Python via the project environment (`pixi run python …`) so PyGMT/GDAL are
on the path. In a notebook started from the pixi env, use `fig.show()` instead
of `save_figure(...)`.

---

## 1. Conventions

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

Every CLI takes an optional `--dir DIR`: a single target directory for that run's
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

---

## 2. Getting data: download a DEM

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
(target directory, see [§1](#1-conventions)), `--ot-api-key-file` (OpenTopography
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

---

## 3. Running lsdtt3

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
    "source_pixels_threshold": 300,
    "write_channel_network_lines": "true",
    "write_channel_network_points": "true",   # stream_order + drainage_area
})
```

See [§6](#6-end-to-end-workflows) for the high-level workflows that chain
preprocessing → analysis → plot.

---

## 4. Maps

All map functions return a `pygmt.Figure` and accept an existing `fig=` so you
can layer them (hillshade → channels/points → basins → decorations).

### 4.1 Hillshade

**CLI**

```bash
pixi run python -m lsdviztools3.cli.plot_hillshade \
  --in demo_output/eigg_glo30.tif --out demo_output/hillshade.png \
  --coordinates local --cmap gray --opaque
```

Flags: `--in`, `--out`, `--input-is-hillshade` (treat input as a precomputed
hillshade), `--clean` (see below), plus the shared style flags.

**Clean maps (`--clean` / `MapStyle(clean=True)`).** Strips all map furniture —
frame ticks, tick marks, tick labels, axis labels, the CRS title, and the
colourbar — leaving just the map inside a thin border, with a small **solid black
rectangle scale bar** in the bottom-left. The bar length is the largest of a set
of round numbers (0.1, 0.5, 1, 2, 3, 4, 5, 8, 10, 15, 20, 25, 40, 50, 75, 100,
150, 200 km) that fits within 1/10 of the map width, labelled e.g. `8 km`. Works
on `render_hillshade` and `render_drape` (and any map sharing `MapStyle`).

**Python**

```python
from lsdviztools3.io.raster import load_raster
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.style import MapStyle

style = MapStyle(cmap="gray", transparent=False)
fig = render_hillshade(load_raster("demo_output/eigg_glo30.tif"),
                       coordinates="local", shade=True, style=style)
save_figure(fig, "demo_output/hillshade.png", style=style)
```

`render_hillshade(raster, *, coordinates="local", shade=True, style=None, fig=None)`.
Set `shade=False` if the raster is already a hillshade.

### 4.2 Channel network coloured by log₁₀ k_sn

The chi network with per-segment `k_sn` is in
`*_chi_fit_segment_lines_mn_*.fgb` (from `lsdtt-chi-analysis`; see
[§6](#6-end-to-end-workflows)). `render_channels(color_by="k_sn")` draws each
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

### 4.3 Points sized by drainage area

Points carrying `drainage_area` (channel-extraction points, or chi data-map
points) can be size-scaled by log₁₀(drainage_area) so big channels render
thicker, and coloured by another column (e.g. Strahler `stream_order`).

```python
from lsdviztools3.render.points import render_points
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.style import MapStyle

fig = render_hillshade("work/eigg/eigg_remove_seas.tif", coordinates="local",
                       shade=True, style=MapStyle(cmap="gray", transparent=False))
fig = render_points(
    "work/eigg/eigg_channel_network_points.fgb",
    size_by="drainage_area",        # log10 -> size_range
    size_range=(1.0, 3.0),          # smallest..largest point, in pt
    color_by="stream_order", log_color=False,
    cmap="batlow", colorbar=True, colorbar_label="Strahler order",
    coordinates="local", style=MapStyle(cmap="batlow"), fig=fig,
)
save_figure(fig, "demo_output/eigg_points.png", style=MapStyle(transparent=False))
```

`render_points(df, *, color_by=None, size_by=None, log_color=True,
size_range=(1,3), cmap=None, colorbar=False, colorbar_label=None,
coordinates="local", style=None, fig=None)`. Input may be a GeoDataFrame or a
path to a point vector. Drop `size_by` for a constant `style.point_size`; drop
`color_by` for a solid colour.

### 4.4 Basins

**CLI**

```bash
pixi run python -m lsdviztools3.cli.plot_basins \
  --base demo_output/eigg_glo30.tif \
  --basins work/eigg/eigg_junction_basin_polygons.fgb \
  --out demo_output/basins.png \
  --coordinates local --fill skyblue --transparency 70 --opaque
```

Flags: `--base`, `--basins`, `--out`, `--fill`, `--transparency`, shared style.

**Python**

```python
from lsdviztools3.render.basins import render_basins
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.io.vector import load_vector
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.style import MapStyle

style = MapStyle(transparent=False)
fig = render_hillshade("demo_output/eigg_glo30.tif", coordinates="local",
                       style=MapStyle(cmap="gray", transparent=False))
fig = render_basins(load_vector("work/eigg/eigg_junction_basin_polygons.fgb"),
                    coordinates="local", fill="skyblue", transparency=70,
                    style=MapStyle(line_color="navy", transparent=False), fig=fig)
save_figure(fig, "demo_output/basins.png", style=style)
```

### 4.5 Locator / orthographic basemap

A globe / aerial locator outlining the DEM footprint.

**CLI**

```bash
pixi run python -m lsdviztools3.cli.plot_basemap \
  --in demo_output/eigg_glo30.tif --out demo_output/basemap.png \
  --mode aerial --elevation-km 50 --look-bearing 0 --look-vertical-angle 0 --opaque
```

Flags: `--in`, `--out`, `--mode {auto,global,aerial}`, `--elevation-km`,
`--view-limit-deg`, `--look-bearing`, `--look-vertical-angle`,
`--grid-spacing-deg`, `--no-north-arrow`, shared style.

**Python**

```python
from lsdviztools3.render.basemap import render_orthographic_basemap
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.style import MapStyle

style = MapStyle(figure_size="10c", transparent=False)
fig = render_orthographic_basemap("demo_output/eigg_glo30.tif", mode="global", style=style)
save_figure(fig, "demo_output/basemap_global.png", style=style)

fig = render_orthographic_basemap("demo_output/eigg_glo30.tif", mode="aerial",
                                  elevation_km=50, look_bearing=0,
                                  look_vertical_angle=0, style=style)
save_figure(fig, "demo_output/basemap_aerial.png", style=style)
```

### 4.6 Decorations: scale bar and colourbar

```python
from lsdviztools3.render.decorations import add_scalebar, add_colourbar
```

- `add_scalebar(fig, extent, crs, *, length_km=None, anchor="BL")` — branches on
  the CRS: a hand-drawn km ruler for projected (linear-km) maps, native GMT
  `map_scale` for geographic ones. `length_km=None` picks a "nice" ~1/5 of the
  map width; `anchor` is `"BL"|"BR"|"TL"|"TR"`. `extent` is the display extent
  `[w, e, s, n]` (km for projected) — obtain it via `grid_to_display` (see
  [§4.2](#42-channel-network-coloured-by-log10-k_sn)).
- `add_colourbar(fig, *, label=None, placement="bottom", position=None, frame=None)`
  — colourbar for the active CPT. `placement` is `"bottom"` (default), `"top"`,
  `"left"`, or `"right"` (offsets are chosen to clear the axis labels); pass a
  raw GMT `-D` string as `position` to override it. `render_channels`/
  `render_points`/`render_drape` call this for you when `colorbar=True`; call it
  directly if you build a CPT yourself.

### 4.7 3D perspective relief (rayshader-style)

A 3D perspective view of the DEM via PyGMT `grdview`: an illuminated terrain
surface seen from a viewing azimuth/elevation, coloured by height (or by an
optional drape raster). Projected DEMs get km x/y axes and an Elevation (m)
z-axis.

```python
from lsdviztools3.render.relief3d import render_relief3d
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.style import MapStyle

style = MapStyle(cmap="geo", figure_size="17c", transparent=False)
fig = render_relief3d(
    "demo_output/eigg_glo30.tif",
    coordinates="local",
    azimuth=135, elevation=30,       # camera bearing + angle above horizon
    vertical_exaggeration=2.0,       # true vertical:horizontal scale ratio
    shade=True,                      # illuminate from the topography
    style=style,
)
save_figure(fig, "demo_output/eigg_relief3d.png", style=style)
```

`render_relief3d(raster, *, coordinates="local", azimuth=135.0, elevation=30.0,
vertical_exaggeration=None, zsize="4c", shade=True, drape=None, cmap=None,
style=None, fig=None)`. `vertical_exaggeration` (projected CRS only) sets the
true vertical:horizontal scale ratio and overrides `zsize`; pass `zsize` (the
raw GMT plotted height of the elevation range) directly for explicit control.

**Drape option.** Pass `drape=<raster>` to colour the surface by another raster
(reprojected and resampled onto the DEM grid) instead of by elevation — e.g. a
k_sn, gradient, or land-cover grid. `cmap` defaults to `style.cmap`.

**Chi network on the 3D surface.** The workflow helper drapes the chi network,
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

### 4.8 Draping a data raster over a precomputed hillshade

A common 2D figure is a continuous data raster — slope, curvature, a k_sn grid,
relief, land cover — coloured by a CPT and laid *semi-transparently* over a
grayscale hillshade, so the shaded-relief texture shows through the colour. This
is the 2D analogue of the 3D `drape=` option ([§4.7](#47-3d-perspective-relief-rayshader-style)).

`render_drape` builds the two layers for you: a **gray** hillshade base, then the
data raster on top through your chosen `cmap` at a given `transparency`. The base
is always gray (it provides texture, not colour); the drape carries the colour.

Both rasters here are **precomputed and live in a different directory** from your
working dir (e.g. an lsdtt3 output folder on another disk) — just point at them
by path. They must share the same CRS and grid (lsdtt3 outputs for one DEM
already do); reproject/resample the drape onto the hillshade grid first if not.

**CLI**

```bash
pixi run python -m lsdviztools3.cli.plot_drape \
  --base /data/lsdtt3_runs/eigg/eigg_hs.tif \
  --drape /data/lsdtt3_runs/eigg/eigg_curvature.tif \
  --out demo_output/eigg_curvature_drape.png \
  --coordinates local --cmap roma --transparency 40 \
  --colorbar --colorbar-label "Curvature" --opaque
```

Flags: `--base`, `--drape`, `--out`, `--shade-base` (treat `--base` as elevation
to illuminate rather than as a precomputed hillshade), `--log` (colour by log10,
clamping values <= 0 to transparent), `--transparency` (0 opaque … 100 invisible,
default 40), `--vmin`/`--vmax` (fixed colour-scale bounds), `--percentile-min`/
`--percentile-max` (bounds from data percentiles), `--out-of-bounds`
(`clamp`/`clip`), `--colorbar`, `--colorbar-label`, `--colorbar-position`
(`bottom`/`top`/`left`/`right`, default `bottom`), `--clean` (strip all furniture,
draw just a solid scale bar — see [§4.1](#41-hillshade)), plus the shared style
flags (`--cmap` sets the drape colour map; the base stays gray).

**Colour maps.** `--cmap` accepts GMT's own maps (`geo`, `roma`, `batlow`,
`viridis`, `seis`, …), **matplotlib colour maps** (`cividis`, `Spectral`, `RdBu`,
`twilight`, … — converted to a GMT CPT on the fly; needs the `charts`/matplotlib
extra), a single named colour (`black`, `red`) for a solid mask, a `colour,colour`
list, or a `.cpt` path. Plain GMT cannot read matplotlib maps; this conversion is
done by `render.colormaps.resolve_cmap`, which all drape calls route through.

**Bounds and thresholding.** By default the colour scale spans the drape's data
range. Narrow it with fixed `vmin`/`vmax`, or derive it from the data with
`percentile_min`/`percentile_max` (0–100). `out_of_bounds` decides what happens
to values outside the range: `"clamp"` (default) paints them the nearest end
colour; `"clip"` leaves them **unrendered** (transparent — the hillshade shows
through). Combining `percentile_min` + `clip` + a solid `cmap` gives a threshold
mask, e.g. "paint everything above the 90th percentile black, show bare hillshade
elsewhere":

```bash
pixi run python -m lsdviztools3.cli.plot_drape \
  --base /data/rio_tucuy/rt_hillshade.tif \
  --drape /data/rio_tucuy/rt_surf_direct_polyfit_curvedness.tif \
  --out figures/curvedness_p90_black.png \
  --coordinates local --cmap black \
  --percentile-min 90 --out-of-bounds clip --transparency 0 --opaque
```

**Python**

```python
from pathlib import Path

from lsdviztools3.render.channels import render_drape
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.style import MapStyle

data_dir = Path("/data/lsdtt3_runs/eigg")     # precomputed rasters elsewhere

style = MapStyle(cmap="roma", figure_size="17c", transparent=False)
fig = render_drape(
    str(data_dir / "eigg_hs.tif"),            # gray hillshade base
    str(data_dir / "eigg_curvature.tif"),     # coloured drape
    coordinates="local",
    base_is_hillshade=True,                   # base file is already a hillshade
    transparency=40,                          # 0 opaque … 100 invisible
    cmap="roma",                              # drape colour map
    colorbar=True, colorbar_label="Curvature",
    colorbar_position="bottom",               # or "top" / "left" / "right"
    style=style,
)
save_figure(fig, "demo_output/eigg_curvature_drape.png", style=style)
```

`render_drape(base_raster, drape_raster, *, coordinates="local",
base_is_hillshade=True, log_color=False, transparency=40.0, cmap=None,
vmin=None, vmax=None, percentile_min=None, percentile_max=None,
out_of_bounds="clamp", colorbar=False, colorbar_label=None,
colorbar_position="bottom", style=None, fig=None)`.

Notes:

- `transparency` is the GMT blend of the colour layer over the hillshade. Lower
  it for more colour, raise it for more terrain texture.
- `log_color=True` (CLI `--log`) colours by `log10` of the drape values — useful
  for k_sn or drainage-area grids — and **clamps values <= 0 to NaN** (drawn
  transparent), so zeros/negatives never break the log. The colourbar then reads
  in log10 units; label it accordingly (e.g. `"log@-10@-k@-sn@-"`).
- For an *elevation* drape (colour by height with hillshade illumination from the
  same DEM) you don't need this — a single
  `render_hillshade(dem, shade=True, style=MapStyle(cmap="geo"))` colours by
  height *and* illuminates. `render_drape` is for draping a **different** raster
  than the one providing the shading.

### 4.9 Overlaying several vector layers on a hillshade

`hillshade_with_vectors` is a general, dataset-agnostic map builder: give it a
base raster and a list of vector layers (any mix of **lines, points and
polygons**) and it draws a hillshade, overlays each layer in its own style,
reprojects everything to a common CRS so the layers line up, and adds a legend.
Use it for any overlay job — channel networks, faults, sample sites, roads,
basins — not just the stream-burning case it grew out of. It dispatches each
layer to `render_channels` / `render_points` / `render_basins` by geometry type,
so plain GeoDataFrames draw reliably (a bare `fig.plot(data=gdf)` does **not**).

Each layer is a dict. Only `path` is required; the rest fall back to sensible
defaults:

| key | applies to | default | meaning |
|-----|-----------|---------|---------|
| `path` | all | — | vector file path, or a `GeoDataFrame` |
| `layer` | all | `None` | sub-layer name for multi-layer files (GeoPackage) |
| `kind` | all | `"auto"` | `"auto"` (detect) / `"line"` / `"point"` / `"polygon"` |
| `label` | all | — | legend entry (omit to leave the layer out of the legend) |
| `color` | all | `"black"` | line / outline / point colour |
| `width` | line, polygon outline | `"1p"` | pen width |
| `size` | point | `"0.15c"` | symbol size |
| `symbol` | point | `"c"` | GMT symbol code (`c` circle, `s` square, `t` triangle…) |
| `fill` | point, polygon | `= color` | fill colour |
| `transparency` | polygon | `60` | fill transparency (0 opaque … 100 invisible) |

Layers draw in list order (first = bottom).

**Python**

```python
from lsdviztools3.workflows import hillshade_with_vectors

hillshade_with_vectors(
    "dem.tif",
    [
        {"path": "basins.fgb", "kind": "polygon", "fill": "green",
         "transparency": 60, "label": "basins"},
        {"path": "channels.fgb", "kind": "line", "color": "red",
         "width": "0.8p", "label": "channels"},
        {"path": "sites.gpkg", "layer": "samples", "kind": "point",
         "color": "yellow", "symbol": "c", "size": "0.18c", "label": "samples"},
    ],
    "overlay.png",
    region=[xmin, xmax, ymin, ymax],   # optional zoom, in the raster's native CRS
)
```

`hillshade_with_vectors(base_raster, layers, savepath=None, *,
coordinates="local", shade=True, base_cmap=None, width=MapStyle.DOUBLE_COLUMN,
region=None, clip_to_valid=False, legend=True, legend_position="JTR+jTR+o0.2c",
high_res=False, fig=None)`.

Notes:

- `shade=True` treats `base_raster` as elevation (shaded relief); `shade=False`
  draws a pre-computed hillshade. `base_cmap` defaults to `"gray"` for a neutral
  backdrop that makes coloured overlays pop; pass `"geo"` for a coloured
  elevation base.
- `region=[xmin, xmax, ymin, ymax]` (in the raster's **native** CRS) zooms by
  cropping the base raster to a temp GeoTIFF, so the hillshade extent and the
  overlays all follow the zoom.
- `clip_to_valid=True` clips overlays to the raster's finite-data footprint —
  handy when a vector tile is larger than the DEM.
- `savepath=None` returns the `pygmt.Figure` without writing (e.g. for
  `fig.show()` in a notebook, or further layering via `fig=`).

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

---

## 5. Charts (x–y plots)

These use **matplotlib** (the optional `charts` extra) and return
`(fig, ax)`. They consume `lsdtt-chi-analysis` / `lsdtt-swath-profile` output,
the bearing/junction outputs of `lsdtt-channel-network-analysis` (rose and
circular-comparison charts, [§5.3](#53-rose-diagrams-bearings)–[§5.4](#54-circular--directional-comparison)),
and the valley-floor products of `lsdtt-flat-surface-finder` /
`lsdtt-valley-centreline` ([§5.6](#56-valley-floor--hand-charts)–[§5.8](#58-hydraulic-breaks)).
[§8](#8-lsdtt3-output--plot-function-map) maps every lsdtt3 output file to the
function that plots it.

### 5.1 Chi profiles

Two views: a single basin's chi–elevation profile
([§5.1.5](#515-single-basin-profiles-plot_chi_profile)), and the lsdviztools
**stacked chi plot** — every basin normalised to χ=0, then offset along the χ
axis so the profiles sit side-by-side on one pair of axes, coloured by
log₁₀ k_sn ([§5.1.2](#512-stacked-chi-plots-plot_chi_stacked)).

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

#### 5.1.1 Inputs: which lsdtt3 call produces which file

Every chi chart consumes output from **`lsdtt-chi-analysis`**. Nothing is
computed here: chi, the piecewise-linear fit and k_sn all come from the C++
side. The minimal parameter set — exactly what
`workflows.lsdtt3_chi_analysis` sends — is:

```python
from lsdviztools3.driver import run_lsdtt

run_lsdtt("lsdtt-chi-analysis", {
    "read_path": "work/eigg", "read_fname": "eigg_remove_seas.tif",
    "write_path": "work/eigg", "write_fname": "eigg",
    "source_pixels_threshold": 300,
    "n_largest_basins": 5,            # keep >0: see the warning in 5.1.4
    "print_chi_data_maps": "true",    # -> observed points
    "chi_data_maps_output_format": "flatgeobuf",
    "fit_chi_networks": "true",       # -> the fit itself
    "write_fit_tables": "true",       # -> nodes + segments CSVs
    "write_chi_fit_segment_lines": "true",
    "chi_fit_segment_lines_output_format": "flatgeobuf",
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
| `*_chi_data_maps.fgb` / `.csv` | `print_chi_data_maps` | `load_chi_points` | `chi`, `drainage_area`, `flow_distance`, `elevation`, `source_key`, **`basin_key`** |

`k_sn` lives on the **segments**, not the nodes. `load_chi_fit(nodes, segments)`
(or `assign_node_ksn`) grafts it on by finding, within each
`(basin_index, branch_id)`, the segment whose `[downstream_chi, upstream_chi]`
interval contains each node — nodes just outside every interval take the
nearest segment in the branch.

Two more `lsdtt-chi-analysis` products are used elsewhere in this guide, not by
the chi charts: `*_chi_ksn_basin_summary_mn_*.csv` (per-basin k_sn statistics)
and the segment-lines fgb again, as a **map** layer in
[§4.2](#42-channel-network-coloured-by-log10-k_sn).

#### 5.1.2 Stacked chi plots (`plot_chi_stacked`)

```python
plot_chi_stacked(data, *, basins=None, x_offset=None, color_by="k_sn",
                 log_color=True, cmap="viridis", point_size=3.0, show_fit=True,
                 observed=None, observed_color="0.6", observed_size=1.5,
                 label_basins=True, ax=None, colorbar=True)  ->  (fig, ax)
```

| Parameter | Default | Meaning |
|---|---|---|
| `data` | — | Nodes **or** segments table; see [§5.1.3](#513-nodes-table-vs-segments-table). Must carry `basin_index`. |
| `basins` | all present | Ordered `basin_index` values to stack — also the only way to reorder the stack. |
| `x_offset` | 20% of the widest basin's χ span | Horizontal gap between successive basins, in χ units. |
| `color_by` | `"k_sn"` | Column used for colour. |
| `log_color` | `True` | Colour on a log₁₀ scale. See the flooring rule below. |
| `cmap` | `"viridis"` | Any matplotlib colourmap. |
| `point_size` | `3.0` | Marker size — nodes table only. |
| `show_fit` | `True` | Overlay the fitted profile per `branch_id` — nodes table only (a segments table *is* the fit). |
| `observed` | `None` | Chi data-map frame overlaid as neutral points; see [§5.1.4](#514-overlaying-observed-points). |
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

#### 5.1.3 Nodes table vs segments table

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

#### 5.1.4 Overlaying observed points

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

#### 5.1.5 Single-basin profiles (`plot_chi_profile`)

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

### 5.2 Swath profiles

Distance-along-baseline vs value, with a min/max envelope, an inter-quartile
band, and median + mean.

```python
from lsdviztools3.workflows import lsdtt3_swath_profile, render_swath_profile
out = lsdtt3_swath_profile("demo_output/eigg_glo30.tif", "work/eigg_swath",
                           stem="eigg", source_threshold=300,
                           half_width_metres=1000, bin_width_metres=500)
render_swath_profile(out, "demo_output/eigg_swath.png")
```

Plot a swath table yourself:

```python
from lsdviztools3.io.swath import load_swath
from lsdviztools3.charts import plot_swath
df = load_swath("work/eigg_swath/eigg_swath_profile.csv")
fig, ax = plot_swath(df, bands=((0, 100), (25, 75)), value_label="Elevation (m)")
```

`plot_swath(df, *, bands=((0,100),(25,75)), show_median=True, show_mean=True,
distance_units="km", value_label="Elevation (m)", color="tab:blue", ax=None)`.

### 5.3 Rose diagrams (bearings)

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

### 5.4 Circular / directional comparison

Compare two angular samples (e.g. channel bearings from two DEMs) with circular
statistics and an overlay plot. These mirror the authoritative C++ circular
tests in lsdtt3's `lsdtt_statistical_tests`; for publication numbers, run
lsdtt3's `lsdtt-bearing-comparison` CLI (see below) and feed its result into the
plot.

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
  (e.g. from lsdtt3's `lsdtt-bearing-comparison` CLI, called via the driver)
  instead of the internally-computed Wasserstein/Kuiper note.

These bearing analyses depend on a recent **lsdtt3** build (see
[§3](#3-running-lsdtt3)): `lsdtt-channel-network-analysis` gained
`only_complete_basins` / `minimum_basin_pixels` to restrict bearing/junction
outputs to complete basins, and there is a new `lsdtt-bearing-comparison` CLI
that computes the authoritative circular-Wasserstein distance and Kuiper test
between two bearing files. See the worked example in
[§6](#6-end-to-end-workflows).

### 5.5 Distributions

`plot_distributions(series, ...)` overlays step-histograms of several labelled
value series — used in the worked examples for junction-angle populations and
for fault-relative deflection split by side of the fault
([§6](#6-end-to-end-workflows)).

```python
from lsdviztools3.charts import plot_distributions
fig, ax = plot_distributions({"dendritic": angles_a, "trellis": angles_b})
```

### 5.6 Valley floor / HAND charts

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

### 5.7 Valley cross sections

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

### 5.8 Hydraulic breaks

Cross-section flooding (`lsdtt-valley-centreline` with
`write_cross_section_features`) records the water-surface elevations separating
flat surfaces from the risers between them, written as local HAND values in
`*_valley_cross_section_hydraulic_breaks.csv`.

| Function | What it draws |
|---|---|
| `plot_hydraulic_breaks_downstream` | Break HAND of every section against distance downstream, one series per break rank. |
| `hydraulic_break_clusters` | Binned 1-D clustering test (default 10 sections per bin) asking whether break elevations fall into distinct levels. Numeric, not a plot — same prototyping caveat as §5.7. |
| `plot_hydraulic_break_clusters` | Break-HAND scatter down-valley with the per-bin cluster levels overlaid. |
| `plot_cross_section_informed_valley_mask` | The cross-section-informed valley floor over a hillshade, beside the histogram-derived mask. |

### 5.9 Cosmo comparison workflows

For cosmogenic workflows, `lsdviztools3` can render standard comparison figures
from lsdtt3 CSV artefacts.

Expected primary input:

- `single_point_scan_comparison.csv` from
  `tests/scripts/run_cosmo_single_point_scaling_scan.py` in lsdtt3.

Optional secondary input:

- sample-level comparison tables (for example OCTOPUS/riversand summary style
  CSVs containing a percent-difference column).

Python API:

```python
from lsdviztools3.charts import load_cosmo_table, plot_cosmo_sweep_panels
from lsdviztools3.workflows import render_cosmo_comparison_figures

scan = load_cosmo_table("../lsdtt3/tests/expected/cosmo/generated/single_point_scan_comparison.csv")
fig, axes = plot_cosmo_sweep_panels(scan, reference="cronus_v3", representative_pressure_hpa=900)

outputs = render_cosmo_comparison_figures(
    scan_comparison_csv="../lsdtt3/tests/expected/cosmo/generated/single_point_scan_comparison.csv",
    sample_comparison_csv="../lsdtt3/tests/expected/cosmo/published/dethier2014_riversand_summary.csv",
    out_dir="demo_output/cosmo",
)
print(outputs)
```

Scripted example:

```bash
python examples/cosmo_comparison_workflow.py \
  --scan-csv ../lsdtt3/tests/expected/cosmo/generated/single_point_scan_comparison.csv \
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

### 5.10 Transient cosmogenic and landslide charts

Charts over the `lsdtt-cosmo-column-tool` `step_change_sweep` CSVs
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

### 5.11 Hillslope charts

Process diagnostics for hillslopes. Each module docstring explains what the
figure tests and how to read it; read it before quoting a result.

| Function | Module | What it draws |
|---|---|---|
| `plot_slope_area` | `charts.slope_area` | The slope–area "boomerang" from slope + drainage-area rasters, with the binned median and the hillslope/channel rollover interval (`binned_median`, `rollover_area`). |
| `plot_gradient_curvature` | `charts.gradient_curvature` | Convexity vs gradient with the linear and Roering nonlinear flux-law predictions overlaid (after Roering, 2008). |
| `plot_estar_rstar` | `charts.estar_rstar` | Measured hilltop E*–R* against the Roering et al. (2007) steady relation, plus the R* residual panel; `summarise` gives per-group residual statistics. |
| `plot_estar_rstar_model_scan` | `charts.estar_rstar` | Analytical, Newton and cellular strip results in E*–R* space. |
| `plot_hillslope_profiles` / `plot_profile_residual_panels` | `charts.hillslope_profiles` | Elevation, gradient, curvature and residual panels for 1-D hillslope comparisons. |
| `plot_residual_density` / `plot_residual_raster` | `charts.hillslope_residuals` | Residual histograms (optionally hillslope-masked) and residual maps centred on zero; `summarise_residuals` tabulates RMS/L∞/percentiles. |
| `plot_transient_residuals` / `plot_strip_transient_metrics` | `charts.hillslope_transient_residuals` | RMS and max residuals through time; relief, max gradient and divide Laplacian changes per model. |

Hilltop/channel coupling (`lsdtt-hillslope-channel-coupling` output, loaded with
`io.hilltop_coupling.load_hilltop_outputs`) has its own chart set in
`charts.hilltop_coupling`: `plot_traces_on_hillshade`, `plot_trace_outcomes`,
`plot_ridgelines_and_sides`, `plot_length_distributions`,
`plot_ridgeline_profiles`, `plot_looped_trace_profiles`.

```bash
# hilltop coupling figure set from one output directory or <write_path>/<write_fname> stem
pixi run lsdviz3-hilltop-coupling --prefix ../data/run/run_stem --outdir demo_output/hilltop

# E*-R* figure + per-group residual summary CSV
pixi run lsdviz3-estar-rstar --metrics hilltop_metrics.csv --outdir demo_output/estar

# paired Newton vs cellular hillslope solver animation
pixi run lsdviz3-hillslope-comparison --newton runs/newton --cellular runs/cellular \
  --out demo_output/hillslope_comparison.gif --layout 2d
```

The animation is built by
`render.hillslope_comparison_animation.render_hillslope_comparison_animation`;
`examples/hillslope_comparison.py` drives the whole comparison workflow.

### 5.12 Fluvial-surface ladder and valley sections

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
[§5.7](#57-valley-cross-sections), which work from the same samples table.

### 5.13 Divide migration charts

`lsdviz3-divide-migration` ([§7](#7-cli-reference)) renders the full figure set
for an `lsdtt-divide-migration` run. The pieces are available individually:

| Function | What it draws |
|---|---|
| `render.divide_migration.render_divide_migration` | PyGMT map: grey hillshade + divide polylines + channel heads coloured by one metric. |
| `plot_pair_comparison` | Side A and side B means for one metric, with Δ = A − B and its uncertainty. Side A is always the lower-numbered unit. |
| `plot_head_metric_vs_distance` | Per-head metric against distance to the nearest divide, split by unit. |
| `aggressor_annotation` | Text describing which side Δ points at, and whether the uncertainty supports it. |

Loaders live in `io.divide_migration` (`find_outputs`, `load_divides`,
`load_channel_heads`, `load_pair_summary`).

### 5.14 Network geometry, synthetic networks and valley spacing

| Function | Module | What it draws |
|---|---|---|
| `plot_junction_angle_cdf`, `plot_delta_L_cdf`, `plot_edge_basin_cdf`, `plot_metric_panels`, `plot_metric_vs_temperature` | `charts.network_geometry` | Network-geometry metrics across a sweep (see the rows in [§8](#8-lsdtt3-output--plot-function-map)). |
| `plot_hillshade`, `plot_network`, `plot_network_on_hillshade` | `charts.network_maps` | Matching map panels per synthetic-network case: relief, the network alone (line width ~ log area), and both together. |
| `plot_metric_by_temperature`, `plot_branching_range`, `plot_flow_change_fraction` | `charts.snapping` | The OCN snapping experiment: does annealing temperature produce a range of branching, and does snapping reduce it? |
| `plot_valley_spectrum`, `plot_spacing_vs_ci`, `plot_landscape_panel` | `charts.valley_spacing` | Valley spacing from `lsdtt-spectral-analysis` (`_valley_spacing.csv`, `_valley_spacing_spectrum.csv`): one landscape's spectrum with both spacing measures, spacing vs channelization index across a sweep, and a side-by-side hillshade panel. |
| `plot_significant_frequencies_polar` | `charts.spectral` | Significant 2D spectrum cells from `lsdtt-spectral-analysis` (`_significant_frequencies_2d.csv`, or the legacy `_significant_freq.csv`) on a North-up polar axis: bearing (wave vector, or crest with `orientation="crest"`) against frequency, rings labelled in wavelength, coloured by normalised PSD. A run with nothing significant draws empty axes with a note. |
| `plot_spectral_tile_summary`, `plot_normalised_spectrum_2d`, `plot_radial_spectrum` | `charts.spectral` | One `lsdtt-spectral-analysis` run on a 2 × 2 figure: the input raster (found from `_ingested_params.param`), the normalised 2D spectrum north-up with significant cells outlined, the binned radial spectrum with background, normalised PSD and rollover, and the polar plot. The two spectrum panels are also available on their own. Load the run with `io.spectral.load_spectral_analysis`; panels whose outputs were switched off show a note. |

`examples/synthetic_networks_demo.py` drives `lsdtt-synthetic-networks` and
uses these charts end to end.

---

## 6. End-to-end workflows

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
| `render_hand_profile(…)` / `render_hand_window_populations(…)` | the [§5.6](#56-valley-floor--hand-charts) HAND charts, straight from lsdtt3 output paths |
| `hillshade_with_vectors(dem, layers, png)` | hillshade with any number of vector layers overlaid ([§4.9](#49-overlaying-several-vector-layers-on-a-hillshade)) |
| `compare_stream_networks(…)` | stream-burning preset over `hillshade_with_vectors` |
| `plot_segment_classification(…)` / `plot_algoviz(…)` | stream-burning diagnostics: trace-divergence per segment; replay of the `snap_to_accumulation` source-path burn |
| `sample_dem_along_line(dem, line)` / `trace_d8_flowpath(dem, start)` | sample elevation along a polyline; trace a D8 steepest-descent path |
| `plot_long_profile(profiles, …)` | elevation-vs-distance long profiles on shared axes |

`lsdtt3_chi_analysis` selects the **N largest complete basins** and moves
interior (coastal) outlets upstream so island basins pass the completeness test
— see the lsdtt3 switches `n_largest_basins` and
`n_pixels_move_outlet_upstream`.

### Worked example: channel-bearing comparison (dendritic vs trellis)

`examples/channel_bearing_comparison.py` is an end-to-end script that compares
drainage-network *orientation* between two sites: an Allegheny Plateau tile
(dendritic) and a Valley and Ridge tile (trellis, strongly structurally
aligned). For each site it downloads a COP30 tile, fills and hillshades it
(`lsdtt-raster-preprocessing` with `write_fill_raster` / `write_hillshade_raster`),
then runs `lsdtt-channel-network-analysis` with `only_complete_basins=true` to
get channel-link bearings, junction angles, and outlet→centroid bearings
restricted to complete basins, plus `lsdtt-channel-extraction` for the plottable
network lines. It then renders a clean hillshade + channel-network overlay per
site and produces:

- a rose of channel bearings by Strahler stream order ([§5.3](#53-rose-diagrams-bearings)),
- a rose of outlet→centroid bearings,
- a junction-angle distribution comparison (`plot_distributions`),
- a circular channel-bearing comparison ([§5.4](#54-circular--directional-comparison))
  annotated with lsdtt3's authoritative `lsdtt-bearing-comparison` result.

Outputs land under `../data/channel_bearing_comparison` (one subdirectory per
site, plus a `comparison/` dir). A recent lsdtt3 build is required for the
`only_complete_basins` switch and the `lsdtt-bearing-comparison` CLI; the script
falls back gracefully (Python circular stats) if either is missing.

```bash
pixi run python examples/channel_bearing_comparison.py
pixi run python examples/channel_bearing_comparison.py --skip-compute  # replot only
```

### Worked example: basin geometry vs a fault (Parkfield, San Andreas)

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
- a rose of outlet→COM bearings ([§5.3](#53-rose-diagrams-bearings)) with the
  fault direction marked;
- the deflection-from-fault-orthogonal distribution split by side of the fault
  ([§5.4](#54-circular--directional-comparison) / `plot_distributions`), to test
  for shear-related bending.

This illustrates the intended division of labour: **lsdtt3 computes, lsdviztools3
plots**. See the lsdtt3 basin-extraction reference for the
`select_basins_from_line_crossings` parameters and the CSV columns. Outputs land
under `../data/parkfield`.

```bash
pixi run python examples/fault_basin_deformation.py
pixi run python examples/fault_basin_deformation.py --skip-compute  # replot only
```

### Worked example: basins draining into a single channel (Nujiang)

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

### Other example scripts

| Script | What it does |
|---|---|
| `examples/island_demo.py` | Keyless DEM download + hillshade (`pixi run demo-eigg`, `demo-rum`, `demo-fiji`, `demo-solomon`). |
| `examples/lsdtt3_workflow.py` | Download → lsdtt3 hillshade/channels/largest basin → composite map (`pixi run demo-lsdtt3`). |
| `examples/cosmo_comparison_workflow.py` | Cosmogenic sweep comparison figures ([§5.9](#59-cosmo-comparison-workflows)). |
| `examples/synthetic_networks_demo.py` | Synthetic OCN networks: generate, render and compare junction geometry ([§5.14](#514-network-geometry-synthetic-networks-and-valley-spacing)). |
| `examples/hillslope_comparison.py` | Newton vs cellular hillslope solver comparison ([§5.11](#511-hillslope-charts)). |
| `examples/perron_ren_p1_p2.py` … `perron_ren_p6.py` | Perron/Ren hillslope benchmark figures: P1–P2 outputs, P3 timestep independence, P4 transient convergence and runtimes, P5 stability envelope, P6 spatial convergence and anisotropy. |

---

## 7. CLI reference

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-hillshade` | shaded-relief PNG | `--in --out --input-is-hillshade` + style (incl. `--clean`) |
| `lsdviz3-drape` | gray hillshade + coloured semi-transparent drape | `--base --drape --out --shade-base --log --transparency --vmin --vmax --percentile-min --percentile-max --out-of-bounds --colorbar --colorbar-label --colorbar-position` + style (`--cmap` = drape; GMT or matplotlib maps, or a solid colour) |
| `lsdviz3-sti` | hillshade + STI rasters (single file or directory batch) | single mode: `--base --sti --out`; directory mode: `--sti-dir [--base]` (auto-plots all discovered STI rasters); options: `--transparency --colorbar` + style |
| `lsdviz3-channels` | hillshade + channel lines | `--base --lines --out` + style |
| `lsdviz3-basins` | hillshade + basin polygons | `--base --basins --out --fill --transparency` + style |
| `lsdviz3-divide-migration` | `lsdtt-divide-migration` figure set: hillshade + divide polylines + channel heads coloured by a metric, the side-vs-side pair chart (Δ = A − B, A the lower-numbered unit, error bar only when the run's uncertainty method is a dispersion), and a per-head scatter against distance to divide | `--prefix --outdir [--base --base-is-hillshade --metric --point-size --divide-color --alpha --stem --chart-dpi]` + style |
| `lsdviz3-cosmo-basins` | standard hillshade + basins colored by inferred erosion rate from cosmo CSV | `--base --basins --cosmo-csv --out [--hillshade --joined-basins-out --colorbar --scaling-scheme --muon-scheme --atmosphere --erate-input-units --allow-incomplete-basins]` |
| `lsdviz3-chi-cosmo` | chi basin workflow bundle: stacked chi + per-basin chi PNGs + basin k_sn summary/map + optional k_sn-vs-erate | `--fit-nodes-csv --fit-segments-csv --outdir [--base --hillshade --basins --cosmo-csv --theta --scaling-scheme --atmosphere --muon-scheme --erate-input-units --allow-incomplete-basins --stem --dpi]` |
| `lsdviz3-basemap` | global/aerial locator | `--in --out --mode --elevation-km --view-limit-deg --look-bearing --look-vertical-angle --grid-spacing-deg --no-north-arrow` + style |
| `lsdviz3-animate` | animate a transient `lsdtt-lem` run as a GIF | `--in --out --glob --stats-csv --fps --loop --no-annotate --keep-frames` + style |
| `lsdviz3-lem-pulse` | animate an `lsdtt-lem` incision-pulse run (hillshade + erosion drape + channels + history panel) | `--in --out --prefix --ext --stats-csv --erosion-cmap --erosion-vmin --erosion-vmax --erosion-percentile --erosion-transparency --channel-color --no-timeseries --no-annotate --fps --loop --keep-frames` + style |
| `lsdviz3-lem-ci` | animate an `lsdtt-lem` channelization-index ladder run (CI history panel + map or oblique landscape) | `--in --out --view --prefix --ext --stats-csv --no-initial --map-width --azimuth --elevation --vertical-exaggeration --landscape-cmap --fps --frame-dpi --loop --keep-frames` + style |
| `lsdviz3-hillslope-comparison` | animate paired Newton vs cellular hillslope rasters ([§5.11](#511-hillslope-charts)) | `--newton --cellular --out [--glob --stats-csv --fps --loop --keep-frames --layout]` |
| `lsdviz3-hilltop-coupling` | hilltop/channel coupling figure set ([§5.11](#511-hillslope-charts)) | `--prefix --outdir [--dem --channels --max-traces --dpi]` |
| `lsdviz3-estar-rstar` | measured E*–R* vs Roering et al. (2007) theory + residual summary CSV | `--metrics --outdir [--group-col --title --dpi]` |
| `lsdviz3-ladder-windows` | fluvial-surface ladder break selection along each channel ([§5.12](#512-fluvial-surface-ladder-and-valley-sections)) | `--windows --out [--summary --title --dpi]` |
| `lsdviz3-ladder-density` | per-window HAND density with ladder breaks | `--histograms --windows --out [--source-id --max-panels --title --dpi]` |
| `lsdviz3-ladder-bands` | ladder band index over a hillshade | `--hillshade --bands --out [--channel --centreline --max-band --title --dpi]` |
| `lsdviz3-valley-sections` | valley cross-section panels with picked floor and optional HAND field | `--samples --out [--source-id --field --bands --field-label --max-panels --n-cols --summary --title --dpi]` |
| `lsdviz3-valley-section-map` | cross-section lines vs centreline and channel network | `--hillshade --sections --out [--centreline --channel --samples --sides --source-id --label-every --title --dpi]` |
| `lsdviz3-valley-widths` | valley-width network: section map by excluded_reason, width profile down the largest trunk, channel vs centreline map | `--stem [--base --valley-mask --out-map --out-profile --out-channels --out-centreline-profile --out-terrace-profile --terrace-field --include-floodplain --min-bin-pixels --max-relative-elevation --segment-id --map-segments --bbox --hide-excluded --title --dpi]` |
| `lsdviz3-cosmo-compare` | cosmogenic comparison figures from standard lsdtt3 CSVs | `--scan-csv --outdir [--sample-csv --reference --pressure --erate --stem --dpi]` |
| `lsdviz3-cosmo-published` | residual figure(s) from lsdtt3 published summary tables | `--published-dir --datasets [--mode --outdir --stem --dpi]` |
| `lsdviz3-cosmo-transient` | transient (step-change) cosmogenic figure families ([§5.10](#510-transient-cosmogenic-and-landslide-charts)) | `--sweep-csv --outdir [--value-prefix --families --concentration-xmax --erate-ratio-xmax --no-title --dpi]` |
| `lsdviz3-grabdem` | download a DEM | `--bbox --source --dataset --out --dir --ot-api-key-file` |
| `lsdviz3-lsdtt` | run an lsdtt3 tool | `tool --params --data-path --lsdtt-dir --set-lsdtt-dir --cwd --dir` |

Shared **style flags** (all PyGMT map/animation CLIs — the matplotlib chart
CLIs take only the flags listed in their row): `--coordinates --cmap --projection
--figure-size --line-width --line-color --font-annotation --font-label
--font-title --dpi --high-res --opaque --clean`.

Every command also takes `--dir DIR` — a target directory created if absent, into
which a relative `--out` (and, for `lsdviz3-lsdtt`, the lsdtt3 read/write/cwd)
is placed. See [§1](#1-conventions).

> Chi, swath, channel-by-k_sn, point, and 3D-relief plots currently have **no
> dedicated CLI** — use the Python API / workflows above. The `lsdviz3-lsdtt`
> command runs the underlying lsdtt3 tools that produce their input.

---

## 8. lsdtt3 output → plot function map

lsdviztools3 never computes topographic metrics: every chart and map is a view
of a file **lsdtt3** wrote. This table is the linkage — which CLI, which switch,
which file, which columns, which loader, which plotter. Filenames are shown
with the `<stem>` = lsdtt3 `write_fname` prefix stripped.

| lsdtt3 CLI | Switch that writes it | Output file | Key columns / bands | Loader | Plot function |
|---|---|---|---|---|---|
| `lsdtt-raster-preprocessing` | `write_hillshade_raster` | `_hillshade.tif` | single band | (path) | `render_hillshade` ([§4.1](#41-hillshade)) |
| `lsdtt-raster-preprocessing` | `write_fill_raster` / default | `_fill.tif`, `_remove_seas.tif` | elevation | (path) | base layer for every map; `render_relief3d` ([§4.7](#47-3d-perspective-relief-rayshader-style)) |
| `lsdtt-channel-extraction` | `write_channel_network_lines` | `_channel_network_lines.fgb` | `stream_order` | `load_vector` | `render_channels` ([§4.2](#42-channel-network-coloured-by-log10-k_sn)) |
| `lsdtt-channel-extraction` | `write_channel_network_points` | `_channel_network_points.fgb` | `flat`, `row`, `col`, `stream_order`, `drainage_area` — the **`.csv` form omits `drainage_area`**, so ask for `flatgeobuf`/`gpkg` if you want to size points by it | `load_vector` | `render_points` ([§4.3](#43-points-sized-by-drainage-area)) |
| `lsdtt-basin-extraction` | `write_basin_outlines` (+ `_as_polygon`) | `_junction_basin_polygons.fgb` | basin id | `load_vector` | `render_basins` ([§4.4](#44-basins)) |
| `lsdtt-basin-extraction` | `select_basins_from_line_crossings` | `_line_crossing_basins.csv` | outlet, centre of mass, area, outlet→COM bearing, fault azimuth, side, relative bearing | `load_table` | `rose_diagram`, `plot_distributions` ([§6](#6-end-to-end-workflows)) |
| `lsdtt-chi-analysis` | `write_fit_tables` | `_chi_fit_nodes_mn_*.csv` | `basin_index`, `chi`, `observed_elevation`, `fitted_elevation`, `branch_id` (no `k_sn`) | `load_chi_fit_nodes` / `load_chi_fit` | `plot_chi_stacked`, `plot_chi_profile` ([§5.1](#51-chi-profiles)) |
| `lsdtt-chi-analysis` | `write_fit_tables` | `_chi_fit_segments_mn_*.csv` | `basin_index`, `branch_id`, `downstream_chi`, `upstream_chi`, `k_sn` | `load_chi_fit_segments` | `plot_chi_stacked`, `plot_chi_profile` |
| `lsdtt-chi-analysis` | `write_chi_fit_segment_lines` | `_chi_fit_segment_lines_mn_*.fgb` | as above + geometry | `load_chi_segment_lines` | `plot_chi_*`; `render_chi_ksn_map` / `render_chi_ksn_relief3d` ([§4.2](#42-channel-network-coloured-by-log10-k_sn)) |
| `lsdtt-chi-analysis` | `print_chi_data_maps` | `_chi_data_maps.fgb` | `chi`, `elevation`, `drainage_area`, `flow_distance`, `basin_key`, `source_key` | `load_chi_points` | `observed=` overlay ([§5.1.4](#514-overlaying-observed-points)); `render_points` |
| `lsdtt-chi-analysis` | `write_ksn_basin_summary` | `_chi_ksn_basin_summary_mn_*.csv` | per-basin k_sn statistics | `load_table` | — (tabular) |
| `lsdtt-swath-profile` | default | `_swath_profile.csv` | `along_axis_centre`, `n_pixels`, `mean`, `p0`…`p100` (21 levels, step 5) | `load_swath` | `plot_swath` ([§5.2](#52-swath-profiles)) |
| `lsdtt-channel-network-analysis` | `write_channel_link_bearings` | `_channel_link_bearings.csv/.fgb` | `bearing_deg`, `bearing_rad`, `stream_order`, `drainage_area_m2`, `total_length_metres`, `fit_linearity`, … | `load_table` / `load_vector` | `rose_diagram`, `plot_circular_comparison` ([§5.3](#53-rose-diagrams-bearings)–[§5.4](#54-circular--directional-comparison)) |
| `lsdtt-channel-network-analysis` | `write_junction_angles` | `_junction_angles.csv/.fgb` | `donor_opening_angle_rad`, `donor_*_receiver_angle_rad`, per-limb `*_bearing_deg`, `*_order`, `*_area_m2` | `load_table`, `load_junction_angle_sets` | `plot_distributions` ([§5.5](#55-distributions)), `plot_junction_angle_cdf` |
| `lsdtt-channel-network-analysis` | `write_area_scaling_summary` | `_area_scaling_summary.csv`, `_area_scaling_curve.csv` | `beta` (MLE), `beta_stderr`, `beta_ols`, `n_used`, **`n_decades`** — check this before quoting `beta` | `load_area_scaling`, `load_area_scaling_curve` | `plot_metric_panels`, `plot_metric_vs_temperature` |
| `lsdtt-ridgeline-analysis` | `print_paired_flowpath_statistics` | `_paired_flowpaths.csv`, `_paired_flowpaths_summary.csv` | `delta_L`, `delta_chi`, `path_length_a/b`; summary adds `n_pairs` and the **rejection counts** — read those before the median | `load_paired_flowpaths`, `load_paired_flowpath_pairs` | `plot_delta_L_cdf`, `plot_metric_panels` |
| `lsdtt-channel-network-analysis` | `write_junction_centre_of_mass_bearings` | `_junction_com_bearings.csv` | outlet→centre-of-mass bearing | `load_table` | `rose_diagram` |
| `lsdtt-bearing-comparison` | default | stdout / result file | circular Wasserstein + Kuiper | — | `annotation=` for `plot_circular_comparison` ([§5.4](#54-circular--directional-comparison)) |
| `lsdtt-flat-surface-finder` | `write_hand_profile_csv` | `_hand_profile.csv` | HAND, distance, `surface_class` (1 = floodplain, 2 = terrace) | `load_table` | `plot_hand_distance` and the rest of [§5.6](#56-valley-floor--hand-charts) |
| `lsdtt-valley-centreline` | `write_valley_cross_sections` | `_valley_cross_section_samples.csv` | per-sample offset, elevation, section id | `load_table` | [§5.7](#57-valley-cross-sections) |
| `lsdtt-valley-centreline` | `write_cross_section_features` | `_valley_cross_section_hydraulic_breaks.csv` | break HAND per section and rank | `load_table` | [§5.8](#58-hydraulic-breaks) |
| `lsdtt-lem` | `frequency_of_output_writing`, `write_erosion_rate_raster`, `write_channel_mask_raster` | `_model_iteration_<N>` rasters, `_model_elevation_stats.csv` | elevation, erosion rate, channel mask, model time | `find_pulse_frames` / `find_transient_frames` | `animate_lem_pulse`, `animate_transient_rasters` (`lsdviz3-lem-pulse`, `lsdviz3-animate`) |
| `lsdtt-lem` | channelization-index ladder run | `_model_iteration_<N>` rasters, `_initial_topography`, `_model_elevation_stats.csv` | elevation, CI history | `find_ci_frames` | `animate_ci_run` (`lsdviz3-lem-ci`) |
| `lsdtt-fluvial-surfaces` | default / `write_ladder_histograms_csv` | `_ladder_windows.csv`, `_ladder_histograms.csv`, `_band_index.tif` | per-window breaks, level centres, status; HAND densities; band index | (path) | `plot_ladder_windows`, `plot_ladder_density`, `plot_ladder_band_map` ([§5.12](#512-fluvial-surface-ladder-and-valley-sections)) |
| `lsdtt-valley-centreline` | `write_valley_cross_sections` | `_valley_cross_sections.fgb`, `_valley_centreline.fgb`, `_valley_cross_section_sides.fgb` | section geometry, centreline, margins | (path) | `plot_valley_sections`, `plot_valley_section_map` ([§5.12](#512-fluvial-surface-ladder-and-valley-sections)) |
| `lsdtt-valley-centreline` | `run_valley_width_network` | `_valley_widths.csv`, `_valley_segments.csv`, `_valley_cross_sections.fgb` (with `segment_id`, `excluded_reason`) | width per station, channel tree, exclusions | (path) | `plot_valley_width_map`, `plot_valley_widths_along_channel`, `plot_channel_and_centreline` (`lsdviz3-valley-widths`, [§5.12](#512-fluvial-surface-ladder-and-valley-sections)) |
| `lsdtt-hillslope-channel-coupling` | default | trace / ridgeline / hilltop outputs | trace outcome, hillslope length, relief, E*, R* | `load_hilltop_outputs` | `charts.hilltop_coupling`, `plot_estar_rstar` (`lsdviz3-hilltop-coupling`, `lsdviz3-estar-rstar`) |
| `lsdtt-cosmo-column-tool` | `step_change_sweep` | `cosmo_transient_sweep_<label>.csv`, `_joined_samples.csv` | normalised concentrations, apparent-rate ratios | `load_cosmo_transient_table`, `load_joined_samples` | [§5.10](#510-transient-cosmogenic-and-landslide-charts) (`lsdviz3-cosmo-transient`) |
| `lsdtt-spectral-analysis` | valley spacing | `_valley_spacing.csv`, `_valley_spacing_spectrum.csv` | spectrum, peaks, spacing measures | `load_valley_spacing`, `load_valley_spacing_sweep` | `plot_valley_spectrum`, `plot_spacing_vs_ci` ([§5.14](#514-network-geometry-synthetic-networks-and-valley-spacing)) |
| `lsdtt-spectral-analysis` | `write_significant_frequencies_2d` (default on) | `_significant_frequencies_2d.csv` (legacy: `_significant_freq.csv`) | bearing, frequency, wavelength, normalised PSD per significant 2D cell | `load_significant_frequencies_2d`, `load_significant_frequencies_2d_sweep`, `summarise_significant_frequencies_2d` (`io.spectral`) | `plot_significant_frequencies_polar` ([§5.14](#514-network-geometry-synthetic-networks-and-valley-spacing)) |
| `lsdtt-spectral-analysis` | `write_radial_spectrum: true`, `write_2d_spectrum_rasters` (default on) | `_radial_spectrum.csv`, `_significant_frequencies.csv`, `_psd_2d.tif`, `_background_psd_2d.tif`, `_significant_frequencies_2d.csv` | radial and 2D spectra, rollover, significant cells | `load_spectral_analysis` (`io.spectral`) | `plot_spectral_tile_summary`, `plot_normalised_spectrum_2d`, `plot_radial_spectrum` ([§5.14](#514-network-geometry-synthetic-networks-and-valley-spacing)) |
| `lsdtt-divide-migration` | default | divides, channel heads, pair summary | per-head metrics, Δ = A − B | `io.divide_migration` loaders | [§5.13](#513-divide-migration-charts) (`lsdviz3-divide-migration`) |

Two rules follow from this table:

1. **If a column is missing, fix the lsdtt3 call, not the plot.** The switch in
   column 2 is what produces it; `lsdviz3-lsdtt` / `run_lsdtt` runs the tool.
2. **Basin numbering is per-product.** `basin_index` (fit tables) and
   `basin_key` (data maps) are separate counters that only coincide when
   lsdtt-chi-analysis filtered to complete basins — see the warning in
   [§5.1.4](#514-overlaying-observed-points).
