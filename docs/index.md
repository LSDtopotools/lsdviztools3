# lsdviztools3

Welcome to `lsdviztools3`, the visualisation package for
[`lsdtt3`](https://github.com/LSDtopotools/lsdtt3) (its documentation is at
<https://lsdtopotools.github.io/lsdtt3/>). It is a ground-up rewrite of
the `lsdviztools` plotting stack around [PyGMT](https://www.pygmt.org/) (maps)
and matplotlib (x–y charts). lsdviztools3 only draws: every map and chart is a
view of a file lsdtt3 wrote.

Installation (pixi) and a quick start are in the
[README](https://github.com/LSDtopotools/lsdviztools3#readme). Each page below
shows its plots both as **command-line** invocations (headless PNG output) and
as **Python / notebook** code.

## Guide

| Page | What it covers |
|---|---|
| [Getting started](getting_started.md) | Coordinate modes, `MapStyle`, saving, target directories; downloading a DEM; running lsdtt3 from Python; the workflow helpers and example scripts |
| [Basic plots](basic_plots.md) | Hillshades, raster drapes, basins, points, locator maps, decorations, 3D relief, vector overlays, swath profiles |
| [Chi analysis](chi_analysis.md) | k_sn maps, stacked chi plots and single-basin profiles (`lsdviz3-chi`) |
| [Channel networks](channel_networks.md) | Channel maps, bearing roses and circular comparisons, junction angles (`lsdviz3-channel-bearings`), network geometry, divide migration, stream-burning QC, worked examples |
| [Lithology](lithology.md) | Geologic unit maps, k_sn by lithologic unit, erodibility calibration |
| [Valleys and cross sections](valleys_and_cross_sections.md) | Floodplains and terraces, method comparisons, the fluvial-surface ladder, valley widths, cross sections, hydraulic breaks, perimeter valley sources |
| [Hillslopes and ridgelines](hillslopes_and_ridgelines.md) | Slope–area, gradient–curvature, E*–R*, hilltop–channel coupling, ridgeline profiles |
| [Landscape evolution](landscape_evolution.md) | Animations of [`lsdtt-lem`](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_lem.html) runs; hillslope solver comparisons |
| [Cosmogenic nuclides](cosmogenic_nuclides.md) | Basin erosion-rate maps, chi with erosion rates, CRONUS and published-data comparisons, transient and landslide experiments |
| [Spectral analysis](spectral_analysis.md) | 2D and radial spectra, significant frequencies, valley spacing |
| [Command reference](cli_reference.md) | Every `lsdviz3-*` command and its flags, by topic |
| [lsdtt3 outputs](lsdtt3_outputs.md) | Which lsdtt3 switch writes each file, and which function plots it |
| [Design](design.md) | For developers: the package's architecture and design decisions |

## Layout

```
docs/                   This guide
examples/               End-to-end scripts (channel bearings, faults, cosmo, ...)
lsdviztools3/           io (data) · render (PyGMT maps) · charts (matplotlib) · cli (headless) · workflows
```

## Quick start

```bash
# one-time install of project env + editable package
pixi install
pixi run install

# render a keyless hillshade demo
pixi run demo-eigg
```
