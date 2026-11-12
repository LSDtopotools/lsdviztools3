# Basic plots

Maps that work for any DEM or lsdtt3 raster: hillshades, drapes, basins,
points, locator maps, 3D relief, vector overlays, and the swath profile.
All map functions return a `pygmt.Figure` and accept an existing `fig=` so you
can layer them (hillshade → channels/points → basins → decorations). The
conventions they share (coordinate modes, `MapStyle`, saving) are in
[Getting started](getting_started.md).
## Hillshade

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
## Draping a data raster over a precomputed hillshade

A common 2D figure is a continuous data raster — slope, curvature, a k_sn grid,
relief, land cover — coloured by a CPT and laid *semi-transparently* over a
grayscale hillshade, so the shaded-relief texture shows through the colour. This
is the 2D analogue of the 3D `drape=` option ([3D perspective relief (rayshader-style)](#3d-perspective-relief-rayshader-style)).

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
draw just a solid scale bar — see [Hillshade](#hillshade)), plus the shared style
flags (`--cmap` sets the drape colour map; the base stays gray).

**Colour maps.** `--cmap` accepts GMT's own maps (`geo`, `roma`, `batlow`,
`viridis`, `seis`, …), **matplotlib colour maps** (`cividis`, `Spectral`, `RdBu`,
`twilight`, … — converted to a GMT CPT on the fly), a single named colour (`black`, `red`) for a solid mask, a `colour,colour`
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

## Colouring in classes

Instead of a continuous colour scale, a map can be coloured in a few classes,
as in a GIS: every value between two breaks gets the same colour, and a legend
listing the classes replaces the colourbar. It works the same way for rasters
(`render_drape`, `lsdviz3-drape`), lines (`render_channels`,
`lsdviz3-channels`, and the k_sn maps of [`lsdviz3-chi` and
`lsdviz3-knickpoints`](chi_analysis.md#colouring-k_sn-in-classes)) and points
(`render_points`).

**Choosing the classes.** Give the breaks yourself, or name a scheme that finds
them in the data:

| `classes` (CLI `--classes`) | Classes |
|---|---|
| a list, tuple or NumPy array of breaks, e.g. `[50, 100, 200]` (CLI `50,100,200`) | `n` breaks make `n + 1` classes: `below 50`, `50 to 100`, `100 to 200`, `200 and above` |
| `"equal_interval"` | `n_classes` equal steps from the lower to the upper limit |
| `"quantile"` | the same number of values in every class (breaks at data values) |
| `"natural_breaks"` (or `"jenks"`) | the Fisher-Jenks optimum: the classes that are each as uniform as possible (least sum of squared deviations from the class means) |

`n_classes` (CLI `--n-classes`, default 5) sets the number of classes for a
scheme. A value equal to a break goes in the class above it; nodata gets no
class and stays transparent. Natural breaks on more than 2000 values use an
even sample of the sorted data, as GIS packages do.

**Limits.** `vmin` and `vmax` (CLI `--vmin`, `--vmax`; `--ksn-vmin`,
`--ksn-vmax` on the chi maps) are the ends of the colour gradient, and bound the
schemes: `equal_interval` runs from `vmin` to `vmax`, and `quantile` and
`natural_breaks` use only the values between them. Values beyond the limits
still fall in the end classes. Without limits the data's minimum and maximum
are used.

**Log scale.** With `log_color=True` (CLI `--log`; always on for k_sn), the
schemes and the colours work in log10: equal intervals become equal ratios.
Breaks you give are in data units. Values <= 0 get no class.

**Colours.** By default each class takes the colour of the colour map at the
**centre of the class**, on the gradient from `vmin` to `vmax`: the colour a
continuous map with the same limits would give it there. Unevenly spaced
breaks therefore give unevenly spaced colours, and a classed map matches its
continuous version. Two other choices:

- `class_color_mode="distinct"` (CLI `--class-color-mode distinct`) takes the
  colours of a qualitative map in order, one per class, for classes that should
  look different rather than ordered. This is automatic for matplotlib's
  qualitative maps (`tab10`, `tab20`, `Set1`, `Set2`, `Set3`, `Dark2`,
  `Paired`, `Accent`, `Pastel1`, `Pastel2`); a map with fewer colours than
  classes is refused.
- `class_colors` (CLI `--class-colors red,orange,#2c7fb8`) gives one colour per
  class yourself.

The colour map can be a matplotlib one (`viridis`, `magma_r`, ...) or a GMT one
(`geo`, `roma`, ...).

**Legend.** The class legend sits outside the right edge of the map, so it
covers no data, with its width set from the longest label and title so the text
stays inside the box. Its title is `legend_title` (CLI `--legend-title` on
`lsdviz3-channels`, `--colorbar-label` on `lsdviz3-drape`), in GMT text:
`"k@-sn@-"` gives a subscript. `legend=False` leaves it off.

```bash
# Gradient (m/m) in 6 natural-breaks classes
pixi run lsdviz3-drape --base eigg_hillshade.tif --drape eigg_gradient.tif \
  --out eigg_gradient_classes.png --classes natural_breaks --n-classes 6 \
  --cmap viridis --colorbar-label "Gradient (m/m)"

# Channels by k_sn, breaks given
pixi run lsdviz3-channels --base eigg_hillshade.tif --lines eigg_chi_fit_segment_lines_mn_0p450.fgb \
  --out eigg_ksn_classes.png --color-by k_sn --classes 20,50,100,200 \
  --cmap magma --legend-title "k@-sn@-"
```

```python
from lsdviztools3.render.channels import render_drape, render_channels

fig = render_drape("eigg_hillshade.tif", "eigg_gradient.tif",
                   classes="quantile", n_classes=5, vmin=0, vmax=1.2,
                   cmap="viridis", legend_title="Gradient (m/m)")
fig = render_channels(segments, color_by="k_sn", classes=[20, 50, 100, 200],
                      cmap="magma", legend_title="k@-sn@-", fig=fig)
```

The classification itself is in `lsdviztools3.colour_classes`:
`make_colour_classes(values, classes, ...)` returns the breaks, a colour and
a label per class, and `classify(values)` gives each value's class, for use in
your own plots.

## Basins

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
## Points sized by drainage area

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
## Locator / orthographic basemap

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
## Decorations: scale bar and colourbar

```python
from lsdviztools3.render.decorations import add_scalebar, add_colourbar
```

- `add_scalebar(fig, extent, crs, *, length_km=None, anchor="BL")` — branches on
  the CRS: a hand-drawn km ruler for projected (linear-km) maps, native GMT
  `map_scale` for geographic ones. `length_km=None` picks a "nice" ~1/5 of the
  map width; `anchor` is `"BL"|"BR"|"TL"|"TR"`. `extent` is the display extent
  `[w, e, s, n]` (km for projected) — obtain it via `grid_to_display` (see
  [Channel network coloured by log₁₀ k_sn](chi_analysis.md#channel-network-coloured-by-log10-k_sn)).
- `add_colourbar(fig, *, label=None, placement="bottom", position=None, frame=None)`
  — colourbar for the active CPT. `placement` is `"bottom"` (default), `"top"`,
  `"left"`, or `"right"` (offsets are chosen to clear the axis labels); pass a
  raw GMT `-D` string as `position` to override it. `render_channels`/
  `render_points`/`render_drape` call this for you when `colorbar=True`; call it
  directly if you build a CPT yourself.
## 3D perspective relief (rayshader-style)

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

## Overlaying several vector layers on a hillshade

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

## Swath profiles

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
