# lsdviztools3

Welcome to `lsdviztools3`, the visualisation package for [`lsdtts3`](https://github.com/LSDtopotools/lsdtt3). `lsdtt3` is a software package for topographic analysis, landscape evolution modelling, and hydrology, which runs in a command line environment. Outputs from `lsdtt3` can be visualised in any GIS software package, but `lsdviztools3` has been built to 

* Enable production of publication-ready figures from `lsdtt3` outputs. It can also work with other standard geographic outputs generated from alternative topographic analysis packages. 
* Provide python based wrappers around `lsdtt3` command line programs so analysis can be run in python notebooks. 

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`lsdviztools3` is a ground-up rewrite of the `lsdviztools` plotting stack around
[PyGMT](https://www.pygmt.org/) (maps) and matplotlib (x–y charts). It covers
hillshades, channel networks coloured by log₁₀ k_sn, points size-scaled by
drainage area, basins, locator basemaps, general multi-layer vector overlays on
a hillshade (lines/points/polygons), chi profiles (individual + stacked), swath
profiles, rose diagrams and circular comparison of bearings, and many other
options.

**Full usage guide with examples for every plot (CLI and notebook):
[`docs/documentation.md`](docs/documentation.md).** This README covers install
and the most basic usage only.

## Design at a glance

- **Backend-agnostic IO/data core** (`lsdviztools3.io`): rioxarray/xarray,
  geopandas, pandas. No PyGMT imports here.
- **Thin PyGMT render layer** (`lsdviztools3.render`): every render function
  returns a `pygmt.Figure` (use `.show()` / `.savefig()`).
- **Headless CLIs** (`lsdviztools3.cli`): render straight to PNG, suitable for
  WSL2 / headless servers.

See [`DESIGN.md`](DESIGN.md) for the full design, the data-ingestion mapping,
the secure key-handling spec, and the open questions for review.

## Install (pixi)

We use the `pixi` package manager (it allows mixing of `conda/mamba` and `pip` python packages) for the installation environment. Testing of this package (i.e. pixi and the pixi environment) has been done in a WSL2 instance running Ubuntu 22.04. 

### Installing pixi itself

[pixi](https://pixi.sh) is a fast conda/PyPI package manager that handles the
whole `lsdviztools3` environment for you. Install it once:

```bash
# Linux / macOS / WSL2
curl -fsSL https://pixi.sh/install.sh | bash
```

```powershell
# Windows (PowerShell)
powershell -ExecutionPolicy ByPass -c "irm -useb https://pixi.sh/install.ps1 | iex"
```

Restart your shell (or `source ~/.bashrc`) and confirm it is on your `PATH`:

```bash
pixi --version
```

If `pixi: command not found`, see the [WSL2 notes](#wsl2-notes-ubuntu-2204-default)
below for fixing your `PATH`.

### Project setup

Once you have cloned `lsdviztools3`, navigate into that directory and run:

```bash
# one-time install of project env + editable package
pixi install
pixi run install

# run tests
pixi run test

# run a small keyless demo
pixi run demo-eigg
```

GMT and PyGMT are conda-only; use the provided `pixi.toml` (conda-forge).

Users will be pleased to know that this will work on a server even if you do not have install permissions: `pixi` only requires you to have a writable file system. 

## lsdtt3 requirement for driver/workflows

The `lsdviztools3` render/download modules run inside this package, but the
driver/workflow path requires a separate installation of **lsdtt3** executables
(`lsdtt-raster-preprocessing`, `lsdtt-channel-extraction`, etc.).

- If you only use plotting/downloading, you do not need lsdtt3.
- If you use `lsdviz3-lsdtt` or `lsdviztools3.workflows`, install lsdtt3 first.

Point `lsdviztools3` at your lsdtt3 `bin` directory using one of:

```bash
# one-off per shell
export LSDTT3_DIR=/path/to/lsdtt3/bin

# or persist in lsdviztools3 config
pixi run lsdviz3-lsdtt --set-lsdtt-dir /path/to/lsdtt3/bin
```

## WSL2 notes (Ubuntu 22.04 default)

`lsdviztools3` is designed to run headless in WSL2 Linux environments. The
default/tested setup is **Ubuntu 22.04** using bash.

If `pixi` is installed but `pixi` is not found in your shell, add pixi to
`PATH`:

```bash
export PATH="$HOME/.pixi/bin:$PATH"
pixi --version
```

To make that permanent in Ubuntu 22.04 bash:

```bash
echo 'export PATH="$HOME/.pixi/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

Then use pixi-run commands from the repo root, for example:

```bash
pixi run test
pixi run python -m lsdviztools3.cli.plot_hillshade --help
```

If you need to call pixi without relying on shell startup files, use the
absolute path:

```bash
"$HOME/.pixi/bin/pixi" run test
```

## Troubleshooting

- `pixi: command not found`:
  - Add `~/.pixi/bin` to `PATH` (see WSL2 section), or run
    `"$HOME/.pixi/bin/pixi" ...` directly.
- `ModuleNotFoundError` for `rioxarray`, `pygmt`, etc.:
  - You are likely running outside the pixi environment. Use `pixi run ...`.
- `FileNotFoundError` for `lsdtt-*` tools:
  - lsdtt3 is not installed/configured. Set `LSDTT3_DIR` or run
    `pixi run lsdviz3-lsdtt --set-lsdtt-dir /path/to/lsdtt3/bin`.
- Sea appears coloured (for example red) instead of transparent:
  - Copernicus GLO-30 encodes sea as valid elevation near 0 m, not nodata.
    Use the lsdtt3 preprocessing workflow (`minimum_elevation=0`) or mask sea
    values before rendering.
- Transparent output looks wrong in a viewer/printer:
  - Save opaque output instead (white background): use CLI flag `--opaque`.
    Example: `pixi run python -m lsdviztools3.cli.plot_hillshade --in dem.tif --out map.png --opaque`.
- Network tests skip/fail while offline:
  - Some tests download DEM data (Planetary Computer). Re-run when online, or
    filter tests to local-only subsets.
- Running in WSL2 with GUI concerns:
  - CLI outputs are headless PNG writes; no X display is required.

## Console scripts

| Command            | Module                              |
|--------------------|-------------------------------------|
| `lsdviz3-hillshade`| `lsdviztools3.cli.plot_hillshade`   |
| `lsdviz3-drape`    | `lsdviztools3.cli.plot_drape`       |
| `lsdviz3-channels` | `lsdviztools3.cli.plot_channels`    |
| `lsdviz3-basins`   | `lsdviztools3.cli.plot_basins`      |
| `lsdviz3-basemap`  | `lsdviztools3.cli.plot_basemap`     |
| `lsdviz3-cosmo-compare` | `lsdviztools3.cli.plot_cosmo_compare` |
| `lsdviz3-cosmo-published` | `lsdviztools3.cli.plot_cosmo_published` |
| `lsdviz3-grabdem`  | `lsdviztools3.cli.grab_dem`         |
| `lsdviz3-lsdtt`    | `lsdviztools3.cli.run_lsdtt`        |

## Basic usage

Download a small DEM and render a hillshade — CLI:

```bash
# 1) grab a keyless GLO-30 DEM
pixi run python -m lsdviztools3.cli.grab_dem \
  --bbox -6.20 56.85 -6.06 56.95 --out demo_output/eigg_glo30.tif

# 2) shaded-relief PNG
pixi run python -m lsdviztools3.cli.plot_hillshade \
  --in demo_output/eigg_glo30.tif --out demo_output/eigg_hillshade.png \
  --coordinates local --opaque
```

…or the same in Python / a notebook:

```python
from lsdviztools3.io.raster import load_raster
from lsdviztools3.render.hillshade import render_hillshade
from lsdviztools3.render.base import save_figure
from lsdviztools3.render.style import MapStyle

style = MapStyle(cmap="gray", transparent=False)
fig = render_hillshade(load_raster("demo_output/eigg_glo30.tif"),
                       coordinates="local", style=style)
save_figure(fig, "demo_output/eigg_hillshade.png", style=style)  # or fig.show() in a notebook
```

For every other plot — channels coloured by log₁₀ k_sn, points sized by drainage
area, basins, locator basemaps, chi profiles, swath profiles — and the full list
of options for each (CLI flags and Python keyword arguments), see
**[`docs/documentation.md`](docs/documentation.md)**.

Cosmogenic comparison CLIs:

```bash
# Sweep-based lsdtt3 vs CRONUS/cronuscalc figures
pixi run lsdviz3-cosmo-compare \
  --scan-csv ../lsdtt3/tests/expected/cosmo/generated/single_point_scan_comparison.csv \
  --outdir demo_output/cosmo

# Dataset summary residual figures from published summaries
pixi run lsdviz3-cosmo-published \
  --published-dir ../lsdtt3/tests/expected/cosmo/published \
  --datasets dethier2014,clapp2002,dibiase2010 \
  --mode auto \
  --outdir demo_output/cosmo
```

## What's implemented

Maps (`render.*`, PyGMT): hillshade, channels (incl. colour-by log₁₀ k_sn),
points (incl. size-by drainage area), basins, orthographic/aerial basemap,
3D perspective relief (`render.relief3d`, `grdview`), scale bar + colourbar.
Charts (`charts.*`, matplotlib): chi profiles
(individual + stacked), swath profiles, rose diagrams of bearings (North-up,
grouped/stacked by e.g. stream order), and circular/directional comparison of
two bearing samples (`circular_wasserstein`, `kuiper_two_sample`,
`plot_circular_comparison`), plus cosmogenic comparison charts for lsdtt3-vs-
reference CSVs (`plot_cosmo_sweep_panels`, `plot_cosmo_sample_residuals`).
Plus IO (`io.*`), DEM download
(`download.*`), the lsdtt3 driver (`driver.*`), and chained `workflows`.
Deferred (raise `NotImplementedError`): `charts.plot_profile`,
`charts.plot_slope_area`, `render.decorations.add_inset`. The matplotlib chart
plots need the optional `charts` extra (`pip install lsdviztools3[charts]`); it is
already in the pixi env.

For a quick cosmogenic end-to-end example, see
`examples/cosmo_comparison_workflow.py` and docs section
`5.9 Cosmo comparison workflows` in `docs/documentation.md`.
