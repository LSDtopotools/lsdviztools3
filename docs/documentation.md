# lsdviztools3 — Usage documentation

Worked examples for every plot type, with the options each one accepts, shown
both as **command-line** invocations (headless PNG output) and as **Python /
notebook** code. For installation and a one-line quick start, see the
[README](../README.md). For the architecture, see [DESIGN.md](../DESIGN.md).

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
5. [Statistics (x–y) plots](#5-statistics-xy-plots)
   - [Chi profiles (individual + stacked)](#51-chi-profiles)
   - [Swath profiles](#52-swath-profiles)
6. [End-to-end workflows](#6-end-to-end-workflows)
7. [CLI reference](#7-cli-reference)

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

### Saving

```python
from lsdviztools3.render.base import save_figure
save_figure(fig, "out.png", style=style, high_res=False)  # 300 dpi (600 if high_res)
```

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
{planetary_computer,opentopography}`, `--dataset` (optional), `--out`,
`--ot-api-key-file` (OpenTopography only — the key is never a direct argument).

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
`twilight`, … — converted to a GMT CPT on the fly; needs the `stats`/matplotlib
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

---

## 5. Statistics (x–y) plots

These use **matplotlib** (the optional `stats` extra) and return
`(fig, ax)`. They consume `lsdtt-chi-analysis` / `lsdtt-swath-profile` output.

### 5.1 Chi profiles

Two views: a single basin's chi–elevation profile, and the lsdviztools
"stacked" plot (basins normalised to χ=0 then offset), coloured by log₁₀ k_sn.

```python
from lsdviztools3.io.chi import load_chi_fit
from lsdviztools3.stats import plot_chi_profile, plot_chi_stacked

# Attach per-node k_sn from the segment table (one step):
nodes = load_chi_fit("work/eigg/eigg_chi_fit_nodes_mn_0p450.csv",
                     "work/eigg/eigg_chi_fit_segments_mn_0p450.csv")

fig, ax = plot_chi_stacked(nodes, color_by="k_sn", log_color=True, cmap="viridis")
fig.savefig("demo_output/chi_stacked.png", dpi=300, bbox_inches="tight")

fig, ax = plot_chi_profile(nodes, basin=0)        # one basin
fig.savefig("demo_output/chi_basin0.png", dpi=300, bbox_inches="tight")
```

`plot_chi_stacked(data, *, basins=None, x_offset=None, color_by="k_sn",
log_color=True, cmap="viridis", show_fit=True, observed=None, label_basins=True, ...)`
and `plot_chi_profile(data, *, basin=None, color_by="k_sn", show_fit=True,
observed=None, ...)`.

**FlatGeobuf input.** Both plotters also accept a *segments* table — the fitted
segment-lines fgb — so you can plot without the CSVs, and overlay the observed
points from the data-map fgb:

```python
from lsdviztools3.io.chi import load_chi_segment_lines, load_chi_points
segs = load_chi_segment_lines("work/eigg/eigg_chi_fit_segment_lines_mn_0p450.fgb")
obs  = load_chi_points("work/eigg/eigg_chi_data_maps.fgb")  # observed pts (no k_sn)
fig, ax = plot_chi_stacked(segs, observed=obs)
```

Or run analysis + write a stacked plot and one profile per basin in one call:

```python
from lsdviztools3.workflows import lsdtt3_chi_analysis, render_chi_profiles
out = lsdtt3_chi_analysis("demo_output/eigg_glo30.tif", "work/eigg", stem="eigg",
                          source_threshold=300, n_largest_basins=5)
render_chi_profiles(out, "demo_output/chi", stem="eigg")
```

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
from lsdviztools3.stats import plot_swath
df = load_swath("work/eigg_swath/eigg_swath_profile.csv")
fig, ax = plot_swath(df, bands=((0, 100), (25, 75)), value_label="Elevation (m)")
```

`plot_swath(df, *, bands=((0,100),(25,75)), show_median=True, show_mean=True,
distance_units="km", value_label="Elevation (m)", color="tab:blue", ax=None)`.

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

`lsdtt3_chi_analysis` selects the **N largest complete basins** and moves
interior (coastal) outlets upstream so island basins pass the completeness test
— see the lsdtt3 switches `n_largest_basins` and
`n_pixels_move_outlet_upstream`.

---

## 7. CLI reference

| Command | Purpose | Key flags |
|---------|---------|-----------|
| `lsdviz3-hillshade` | shaded-relief PNG | `--in --out --input-is-hillshade` + style (incl. `--clean`) |
| `lsdviz3-drape` | gray hillshade + coloured semi-transparent drape | `--base --drape --out --shade-base --log --transparency --vmin --vmax --percentile-min --percentile-max --out-of-bounds --colorbar --colorbar-label --colorbar-position` + style (`--cmap` = drape; GMT or matplotlib maps, or a solid colour) |
| `lsdviz3-channels` | hillshade + channel lines | `--base --lines --out` + style |
| `lsdviz3-basins` | hillshade + basin polygons | `--base --basins --out --fill --transparency` + style |
| `lsdviz3-basemap` | global/aerial locator | `--in --out --mode --elevation-km --view-limit-deg --look-bearing --look-vertical-angle --grid-spacing-deg --no-north-arrow` + style |
| `lsdviz3-grabdem` | download a DEM | `--bbox --source --dataset --out --ot-api-key-file` |
| `lsdviz3-lsdtt` | run an lsdtt3 tool | `tool --params --data-path --lsdtt-dir --set-lsdtt-dir --cwd` |

Shared **style flags** (all map CLIs): `--coordinates --cmap --projection
--figure-size --line-width --line-color --font-annotation --font-label
--font-title --dpi --high-res --opaque`.

> Chi, swath, channel-by-k_sn, point, and 3D-relief plots currently have **no
> dedicated CLI** — use the Python API / workflows above. The `lsdviz3-lsdtt`
> command runs the underlying lsdtt3 tools that produce their input.
