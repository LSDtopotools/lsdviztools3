# lsdviztools3 — Design Document

Status (2026-11-12): maps, x–y charts and animations are implemented for most
lsdtt3 outputs: 41 `lsdviz3-*` commands (see the
[Command reference](cli_reference.md)) over the IO, render and charts layers
described here. Three stubs remain deferred (§8). This document describes the
architecture and the rules the code follows; the plot pages describe the
figures, and [lsdtt3 outputs](lsdtt3_outputs.md) maps every lsdtt3 file to the
function that plots it.

## 1. Rationale

`lsdviztools` (v1/v2) grew around matplotlib + a GDAL/OGR IO layer
(`lsdmap_gdalio`) that mixed plotting and IO concerns and predates the modern
geospatial Python stack. `lsdtt3` now emits cloud-native, OGC-friendly formats
(COG/GeoTIFF rasters, FlatGeobuf/GeoPackage vectors, EPSG:4326 CSV points), and
the project standardises on **PyGMT** as the rendering engine for
publication-quality maps.

`lsdviztools3` is therefore a clean rewrite with three firm principles:

1. **Backend-agnostic IO/data core.** Everything in `io/` uses
   rioxarray/xarray, geopandas, and pandas. It contains **no PyGMT imports**
   (the one nuance: `io/grids.py` may call `pygmt.grdgradient` for hillshade
   prep, which is data prep rather than rendering — flagged in its docstring).
   This keeps the data layer testable and reusable without GMT installed.
2. **Thin PyGMT render layer on top.** Everything in `render/` consumes the
   core's xarray/geopandas/pandas objects and emits a `pygmt.Figure`.
3. **lsdtt3 computes, lsdviztools3 draws.** Every map and chart is a view of a
   file lsdtt3 wrote; topographic metrics, fits and statistics belong in lsdtt3
   (C++). The few numeric helpers in `charts/` (the circular statistics in
   `charts/circular.py`, the cross-section hydraulics and break clustering in
   `charts/cross_sections.py` and `charts/hydraulic_breaks.py`) are mirrors or
   prototypes: a figure that needs authoritative numbers annotates lsdtt3's
   (e.g. `_channel_bearing_comparisons.csv`), and a prototype whose result
   matters is moved into lsdtt3. If a plot lacks a column, fix the lsdtt3 call,
   not the plot.

The layering holds in the code: no `io/` or `charts/` module imports PyGMT.
`render/` uses matplotlib only to turn matplotlib colour-map names into GMT
CPTs (`render/colormaps.py`) and for the Newton-vs-cellular animation.

### Figure model

Two engines, chosen by the kind of figure:

- **Maps: PyGMT.** Every `render_*` function returns a `pygmt.Figure`
  (`.show()` in a notebook, `.savefig(path)` to a file) and optionally accepts
  an existing `fig=` so layers (hillshade → drape → channels → basins → points
  → decorations) can be composed onto one figure.
- **x–y charts: matplotlib.** Every `plot_*` function in `charts/` returns
  `(fig, ax)` and accepts an existing `ax=`. matplotlib is a core
  dependency (the `charts` extra is kept, empty, for old install commands).
- **Animations** (`render/animation.py`, `lem_pulse.py`, `lem_ci_animation.py`,
  `hillslope_comparison_animation.py`) render frames and assemble a GIF with
  Pillow.
- **Headless CLIs** build the figure and save it directly — suitable for WSL2 /
  headless servers where no X display is available.

### Publication appearance (defaults in `MapStyle` + `render/base`)

- **Page width:** journal-aware — `MapStyle.SINGLE_COLUMN` = 8.5 cm,
  `DOUBLE_COLUMN`/`MAX_WIDTH` = 17 cm (the hard cap); the bare default is 12 cm.
- **Fonts (8 pt floor):** tick labels `font_annotation` = 9 pt, axis labels
  `font_label` = 11 pt, title `font_title` = 11 pt. Applied via `pygmt.config`
  in `new_figure`.
- **Resolution:** `dpi` = 300 (low-res), `high_res_dpi` = 600; chosen at save
  time by `save_figure(..., high_res=...)`.
- **Projected CRS:** axes are shown in **km** (coords / 1000, via
  `grid_to_display`/`gdf_to_display`) to avoid 6-7 digit easting/northing, with
  `Easting (km)`/`Northing (km)` labels and a `CRS is EPSG:XXXXX` title above the
  map. Geographic CRS keeps degrees with Longitude/Latitude labels.
- **Nodata is transparent**, not grey: `grdimage(nan_transparent=True)` plus a
  transparent background from `save_figure` (`transparent=True`). Note GLO-30
  encodes *sea* as valid 0 m (not nodata) — mask it (lsdtt3 `min_elevation`)
  to get transparent sea.

## 2. Package tree and responsibilities

```
lsdviztools3/
├── docs/                      This guide (mkdocs); design.md is this document.
├── README.md                  Quick start and console-script table.
├── pyproject.toml             PEP 621 metadata, deps, console_scripts.
├── pixi.toml                  conda-forge env (includes gmt, which pip cannot install) + tasks.
├── lsdviztools3/
│   ├── __init__.py            Package version + top-level re-exports.
│   ├── frames.py              Oblique analysis frames (lsdtt3 rotated feature frames), grid convergence.
│   ├── paths.py               Target-directory helpers (--dir).
│   ├── regions.py             Named demo regions; box_from_corners / parse_lat_lon.
│   ├── workflows.py           Composite lsdtt3 -> lsdviztools3 workflows.
│   ├── io/                    Backend-agnostic IO/data core (NO pygmt), one module per lsdtt3 product:
│   │   ├── raster.py, vector.py, table.py   Grids, vectors, CSV points; coordinate handling.
│   │   ├── grids.py           Derived grids (hillshade prep; deferred, §8).
│   │   ├── chi.py, swath.py   Chi fit tables / data maps (node->segment k_sn join); swath percentiles.
│   │   ├── channel_bearings.py, network_geometry.py, snapping.py, coverage.py   Channel networks.
│   │   ├── lithology.py, erodibility_calibration.py                             Lithology.
│   │   ├── divide_migration.py, hilltop_coupling.py, ridgeline.py               Divides, hillslopes, ridges.
│   │   └── spectral.py, valley_spacing.py                                       Spectral analysis.
│   ├── render/                Thin PyGMT render layer; each returns pygmt.Figure:
│   │   ├── base.py, style.py, colormaps.py, decorations.py, graticule.py   Lifecycle, MapStyle, CPTs, furniture.
│   │   ├── hillshade.py, channels.py (render_drape, render_channels), basins.py, points.py,
│   │   │   categorical.py, basemap.py, relief3d.py, sti.py                General maps.
│   │   ├── lithology_map.py, lithology_colours.py, erodibility_calibration_map.py   Lithology.
│   │   ├── floodplain.py, surface_methods.py   Floodplains and terraces.
│   │   ├── divide_migration.py, ridgeline.py   Divides and ridges.
│   │   └── animation.py, lem_pulse.py, lem_ci_animation.py, hillslope_comparison_animation.py   GIFs.
│   ├── charts/                x-y charts (matplotlib), grouped as in the guide:
│   │   ├── chi.py, swath.py, distributions.py
│   │   ├── rose.py, circular.py, channel_bearings.py, network_geometry.py, network_maps.py, snapping.py
│   │   ├── lithology_ksn.py, erodibility_calibration.py
│   │   ├── floodplain_qq.py, surface_methods.py, training_polygons.py, ladder_windows.py,
│   │   │   ladder_density.py, ladder_band_map.py, valley_sections.py, valley_section_map.py,
│   │   │   valley_width_network.py, cross_sections.py, hydraulic_breaks.py, hand.py, perimeter_sources.py
│   │   ├── slope_area.py, gradient_curvature.py, estar_rstar.py, hilltop_coupling.py, ridgeline.py
│   │   ├── hillslope_profiles.py, hillslope_residuals.py, hillslope_transient_residuals.py
│   │   ├── cosmo_comparison.py, cosmo_transient.py, cosmo_transient_spatial.py, landslides.py
│   │   └── spectral.py, valley_spacing.py, divide_migration.py
│   ├── download/              DEM acquisition: apikey.py (secure key resolution), ot_scraper.py.
│   ├── driver/                lsdtt3 process orchestration: config.py (locate the CLIs), lsdtt_driver.py (run_lsdtt).
│   └── cli/                   argparse entry points, one module per lsdviz3-* command;
│                              _common.py (style flags, --dir, IO, save), _ridgeline.py (driver files).
├── examples/                  End-to-end scripts (see Getting started).
└── tests/                     pytest suite (unit, smoke, and lsdtt3 integration tests).
```

## 3. Public API: the core

The topic pages document the chart and map functions with their signatures.
This section covers the core every figure is built on.

### io/raster.py
- `load_raster(path) -> xr.DataArray` — rioxarray; COG/GeoTIFF + legacy ENVI .bil (read-only).
- `to_gmt_grid(da)` — hand a `DataArray` to PyGMT (grdimage-ready grid).
- `get_extent(da) -> [w, e, s, n]`.
- `get_crs(da)`.
- `reproject_raster(da, crs)`.
- `resolve_target_crs(da, coordinates="local") -> CRS` — pick the CRS for the
  selected `coordinates` mode (`"local"` UTM / `"geographic"` EPSG:4326 / explicit).
- `to_coordinates(da, coordinates="local") -> xr.DataArray` — reproject the
  raster into the resolved coordinate system (decision #6).

### io/vector.py
- `load_vector(path, layer=None) -> gpd.GeoDataFrame` — `.fgb` / `.gpkg` / `.shp`.
- `to_xy(gdf)` — extract x/y (and attributes) for PyGMT `plot`.
- `to_coordinates(gdf, coordinates="local") -> gpd.GeoDataFrame` — reproject a
  vector to match the raster's coordinate system.

### io/table.py
- `load_table(path) -> pd.DataFrame` — CSV.
- `latlon_columns(df) -> (lat_col, lon_col)` — detect lat/long columns.

### io/grids.py
- `make_hillshade(da, *, azimuth=315, altitude=45, ...) -> xr.DataArray` —
  gradient-based hillshade prep (may call `pygmt.grdgradient`).

### render/style.py
- `MapStyle` dataclass: `projection`, `frame`, `cmap`, line/colour widths,
  scalebar options, figure size.

### render/base.py
- `new_figure(style) -> pygmt.Figure`.
- `auto_projection(extent, crs, width="12c") -> str`.
- `auto_frame(extent) -> str | list`.
- `is_projected(crs) -> bool`.
- `grid_to_display(grid, crs)` / `gdf_to_display(gdf, crs)` — convert projected
  coords to km for display (decision in §1 appearance).
- `display_frame(crs, style) -> str | list` — Easting/Northing (km) vs Lon/Lat frame.
- `add_crs_label(fig, extent, crs)` — `CRS is EPSG:XXXXX` title above the map.
- `set_grid_gtype(grid, crs)` — set the GMT grid type for the CRS.
- `save_figure(fig, path, *, style=None, high_res=False) -> None`.
- `finalize(fig, savepath=None, *, style=None) -> pygmt.Figure`.


### render: the general map functions
- `render_hillshade(raster, *, coordinates="local", shade=True, style=None, fig=None)`.
- `render_drape(base_raster, drape_raster, *, coordinates="local",
  base_is_hillshade=True, log_color=False, transparency=40.0, cmap=None,
  vmin=None, vmax=None, percentile_min=None, percentile_max=None,
  out_of_bounds="clamp", colorbar=False, ...)` — a coloured, semi-transparent
  data raster over a gray hillshade.
- `render_channels(lines, *, color_by=None, log_color=True, cmap=None,
  colorbar=False, colorbar_label=None, style, fig=None)` — channel lines, with
  `color_by` (e.g. `k_sn`) drawn through a CPT per segment (log10 by default) and
  an optional colourbar legend.
- `render_basins(basins_gdf, *, label_col=None, style, fig=None)`.
- `render_categorical_drape(...)`, `render_binary_drape(...)`.
- `render_points(df, *, color_by=None, size_by=None, log_color=True,
  size_range=(1,3), cmap=None, colorbar=False, style, fig=None)` — points
  coloured by a column and/or size-scaled by log10 of a column (e.g.
  `size_by="drainage_area"` fattens big channels from 1→3 pt).

- Colour maps go through `render.colormaps.resolve_cmap`: GMT names, matplotlib
  names (converted to a CPT), a single colour, a `colour,colour` list or a `.cpt`.
- Lithology unit colours are assigned **by unit name** (`render.lithology_colours`),
  so a unit keeps its colour across the map, the k_sn charts and the
  erodibility figures; `--unit-colours` fixes them from a file.

### render/basemap.py
- `render_orthographic_basemap(raster, *, style=None, fig=None, mode="auto", ...)`
  — a locator map outlining the raster footprint. `mode="global"` draws an
  orthographic hemisphere (never with a north arrow: GMT's rose on a
  whole-globe map breaks the PNG export), `"regional"` calls
  `render_regional_basemap`, `"aerial"` a perspective view from a camera
  `elevation_km` above Earth, and `"auto"` chooses global for a footprint of
  25 square degrees or more and regional otherwise. Aerial controls include
  `look_bearing`, `look_vertical_angle`, `view_twist` and graticule spacing.
- `render_regional_basemap(raster, *, extent_km=None, relief=False, ...)` — a
  square Mercator map around the footprint (ten times its diagonal, at least
  250 km): coastlines, country and state borders, rivers, lakes, optional
  `earth_relief` shading (downloaded by GMT), scale bar, north arrow.
- `add_locator_inset(fig, raster, *, mode="global"|"regional", position="TR",
  width="4c", ...)` — either map in a corner of a map already drawn, via
  `fig.inset`; also `render.decorations.add_inset`. The map commands expose it
  as `--locator` (`cli/_common.add_locator_args`).

### render/relief3d.py
- `render_relief3d(raster, *, coordinates="local", azimuth=135, elevation=30,
  vertical_exaggeration=None, zsize="4c", shade=True, drape=None, cmap=None, ...)`
  — rayshader-style 3D perspective relief via PyGMT `grdview`: an illuminated
  terrain surface coloured by height (or by an optional resampled `drape`
  raster). `vertical_exaggeration` (projected CRS) sets the true vertical:
  horizontal scale ratio and overrides `zsize`.

### render/decorations.py
- `add_colourbar(fig, *, label=None, position=..., frame=None)` — colourbar for
  the active CPT (e.g. the log10 k_sn palette).
- `add_scalebar(fig, extent, crs, *, length_km=None, anchor="BL")` — scale bar
  that branches on CRS: a hand-drawn km ruler on the linear (projected) maps
  (where GMT `map_scale` doesn't apply) and native `basemap(map_scale=...)` on
  geographic (Mercator) maps. Auto length = a 'nice' ~1/5 of the map width.

### download/apikey.py
- `resolve_ot_api_key(key_file=None) -> str`.
- `mask_key_in_url(url, key) -> str`.

### download/ot_scraper.py
- `download_dem(bbox, *, source, dataset, out_file, api_key_file=None)`.

### driver/config.py
Locating the lsdtt3 executables (they are often not on `PATH`, e.g.
`~/.local/lsdtt3/bin`). Resolution order, first found wins:
1. explicit `lsdtt_dir` argument → 2. env var `LSDTT3_DIR` →
3. `[lsdtt] dir` in a config file (`$LSDVIZTOOLS3_CONFIG`, `./lsdviztools3.toml`,
   `~/.config/lsdviztools3/config.toml`) → 4. `PATH` (`shutil.which`).
- `resolve_lsdtt_dir(lsdtt_dir=None) -> Path | None`.
- `resolve_lsdtt_tool(tool, lsdtt_dir=None) -> str` (raises `FileNotFoundError`
  with config guidance if the tool can't be found).
- `write_user_config(lsdtt_dir) -> Path` (persists `[lsdtt] dir`).

### driver/lsdtt_driver.py
- `run_lsdtt(tool, parameters, *, data_path=None, lsdtt_dir=None, extra_args=None, cwd=None, check=True) -> subprocess.CompletedProcess`.
  Resolves the executable via `driver.config`, invokes `<exe> <path> <parameter_file>`
  (lsdtt3's convention; `<path>` defaults to the params file dir), and serialises
  a `dict` of parameters to a temporary `.driver` file. The `lsdviz3-lsdtt` CLI
  exposes `--lsdtt-dir` and `--set-lsdtt-dir DIR`. `clean_gdal_env=True` (default)
  strips conda/pixi `GDAL_DRIVER_PATH`/`GDAL_DATA`/`PROJ_LIB`/`PROJ_DATA` from the
  subprocess env so the external lsdtt3 binary uses its own GDAL (avoids plugin/ABI
  clashes when run inside the pixi environment).

## 4. Commands (`cli/`)

One module per command, registered in `pyproject.toml`; the
[Command reference](cli_reference.md) lists all of them with their flags.
Conventions:

- **PyGMT map and animation commands** share the style flags
  (`_common.add_style_args`: `--coordinates --cmap --projection --figure-size
  --line-width --line-color --font-* --dpi --high-res --opaque --clean`) and
  `--dir` (`_common.add_target_dir_arg`: a target directory, created if
  absent).
- **matplotlib chart commands** take only their own flags, usually `--dpi` and
  an output path or directory.
- **Commands that read one lsdtt3 run** take its directory and `write_fname`
  (`--data-dir` + `--prefix`, or `--stem` = path + `write_fname`) and find the
  files from lsdtt3's naming, rather than one flag per file. A missing input is
  an error naming the lsdtt3 switch that writes it.
- **The ridgeline commands** also read a `key: value` driver file
  (`--driver`), with flags on the command line winning, and write PNG, SVG and
  PDF into `png/`, `svg/` and `pdf/`.

### Demo data (examples & tests)
No raster/vector data is committed (`.gitignore` excludes `*.tif`/`*.bil`/`*.fgb`/…).
Examples and tests **download** small DEMs on demand from the keyless Planetary
Computer (Copernicus GLO-30). Named regions live in `lsdviztools3/regions.py`:
`eigg` (tiny, single tile — the fast test fixture), `rum`, `viti_levu` (Fiji),
`guadalcanal` (Solomon Is.). The larger islands span several GLO-30 tiles and
are merged on download. See `examples/island_demo.py` and the `demo-*` pixi tasks.
Network tests are marked `network` and skip cleanly when offline.

## 5. Data ingestion

[lsdtt3 outputs](lsdtt3_outputs.md) is the full mapping from each lsdtt3
switch and file to its loader and plot function. By format:

| lsdtt3 output                | Format               | Core reader (io/)                  | Drawn with            |
|------------------------------|----------------------|------------------------------------|-----------------------|
| DEM / hillshade / drape grid | COG / GeoTIFF        | `raster.load_raster` (rioxarray)   | `Figure.grdimage`     |
| Legacy raster (read-only)    | ENVI `.bil` + `.hdr` | `raster.load_raster` (rioxarray)   | `Figure.grdimage`     |
| Channel network / lines      | FlatGeobuf `.fgb`    | `vector.load_vector` (geopandas)   | `Figure.plot` (lines) |
| Basins / polygons            | GeoPackage `.gpkg` / `.fgb` | `vector.load_vector` (geopandas) | `Figure.plot` (polygons) |
| Point data (e.g. samples)    | CSV (EPSG:4326 lat/long) | `table.load_table` (pandas)    | `Figure.plot` (points) |
| Tables (fits, statistics, profiles) | CSV           | the topic's `io/` module (pandas)  | matplotlib (`charts/`) |

## 6. Secure key-handling spec

Mirrors the `download_dem_tiles.py` pattern. `resolve_ot_api_key` resolves the
OpenTopography key in this strict order (first found wins):

1. **Environment variable** `OT_API_KEY`.
2. **Key file** (`key_file` argument): the first non-empty, non-comment
   (`#`-prefixed) line. Recommend `chmod 600` and never committing it.

The API key is **never** accepted as a direct CLI argument (to keep it out of
shell history and process listings). Any URL that embeds the key is passed
through `mask_key_in_url(url, key)` before logging/printing so the secret never
appears in console output. `.gitignore` excludes common key-file names.

Planetary Computer is **keyless** (STAC + `planetary_computer.sign`), so the
default DEM path requires no secret at all.

## 7. Charts (x-y) plotting

Every chart module lives in `charts/` and returns `(fig, ax)`. Rules the
modules follow:

- **Draw what lsdtt3 wrote.** Columns are read by name from lsdtt3's tables;
  nothing is normalised, inverted or refitted unless the module says so in its
  docstring.
- **NaN is never turned into 0.** Missing values are dropped, and the newer
  modules (divide migration, lithology, valley widths and others) report how
  many; carry that into any new module.
- **Uncertainty columns are drawn as what they are.** For example, the
  divide-migration pair chart draws `uncertainty_*` as an error bar only when
  the run's `pair_uncertainty_method` is a dispersion (`std_err`, `bootstrap`);
  under `ttest` it is a p-value and is shown as a significance annotation.
- **Each module docstring says what the figure tests and how to read it.**
- **Comparable colours:** colour normalisation is global within a figure, and
  lithology units are coloured by name across figures.

## 8. Deferred

These stubs raise `NotImplementedError`:

| Function | Module | Note |
|---|---|---|
| `plot_profile` | `charts/__init__.py` | Longitudinal channel profiles. |
| `make_hillshade` | `io/grids.py` | Hillshades come from lsdtt3 (`write_hillshade_raster`) or `render_hillshade(shade=True)`. |
| `to_xy` | `io/vector.py` | Vectors are plotted from GeoDataFrames by the render layer instead. |

WCS DEM download is also deferred (decision 7).

## 9. Resolved decisions

> Decisions 1–7 were resolved on 2026-06-17; 8–11 record later decisions.

1. **Stats-plot engine: matplotlib** (not PyGMT-only). Deferred x–y/statistical
   plots (chi, profiles, slope-area, ksn, joyplots) will use matplotlib behind an
   *optional* dependency, keeping the core maps install lean. PyGMT for maps,
   matplotlib for x–y.
2. **Colourmaps: GMT defaults + bundled Scientific Colour Maps** (Crameri:
   `batlow`, `grayC`, `oleron`, …). No bespoke house CPTs in v1; `MapStyle.cmap`
   selects by name. Grey for hillshade, a scientific map for drapes. Revisit only
   if a signature look is wanted.
3. **Default DEM source: keyless Planetary Computer** (Copernicus GLO-30).
   Zero credentials for new users. OpenTopography is opt-in via the secure key
   path (§6).
4. **Basin labels: read the basin/junction id attribute** from the `.fgb`/`.gpkg`
   feature; place the label at the polygon's `representative_point()` (guaranteed
   inside concave basins). Fall back to the junctions table only if the outline
   layer carries no id attribute.
5. **Legacy ENVI `.bil`: read-only — confirmed.** Drop all `.bil` write paths.
   New raster output (e.g. on-the-fly hillshade) is GeoTIFF.
6. **Coordinates: user-selectable**, default **local** (native projected/UTM).
   Render functions take `coordinates="local"` (UTM estimated from extent, GMT
   linear projection — correct scalebars, undistorted local maps) or
   `"geographic"` (EPSG:4326, Mercator), or an explicit CRS. `io.raster.to_coordinates`
   does the reprojection; the render layer sets the GMT grid type accordingly.
   Vectors are reprojected to match the raster.
7. **WCS: deferred** to a later release. v1 `grabdem` = Planetary Computer
   (default) + OpenTopography. The generic `apikey` infra leaves room to add it.
8. **lsdtt3 computes, lsdviztools3 draws** (principle 3). Python numeric
   routines in `charts/` were prototyping; the computation belongs in lsdtt3.
9. **External lsdtt3 runs get a clean GDAL environment.** `run_lsdtt` strips
   the pixi/conda `GDAL_DATA`, `GDAL_DRIVER_PATH`, `PROJ_LIB` and `PROJ_DATA`
   (`clean_gdal_env=True`), so the lsdtt3 binary uses its own GDAL.
10. **Current lsdtt3 parameter names.** The driver and workflows send lsdtt3's
    current names (`threshold_contributing_pixels`, `min_elevation`,
    `write_chi_data_maps`, one `write_vector_format`, ...); lsdtt3 still accepts
    the old ones, with a deprecation warning.
11. **Channel-bearing comparisons come from lsdtt3.** The separate
    `lsdtt-bearing-comparison` CLI is gone; [`lsdtt-channel-network-analysis`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_channel_network_analysis.html)
    (`write_channel_bearing_statistics`) compares groups (Strahler orders,
    basins, regions, lithologic units) with block-permutation p-values, and
    lsdviztools3 draws and annotates its tables.
