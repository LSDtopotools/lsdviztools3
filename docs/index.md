# lsdviztools3

Welcome to `lsdviztools3`, the visualisation package for
[`lsdtt3`](https://github.com/LSDtopotools/lsdtt3). It is a ground-up rewrite of
the `lsdviztools` plotting stack around [PyGMT](https://www.pygmt.org/) (maps)
and matplotlib (x–y charts): hillshades, channel networks coloured by log₁₀
k_sn, basins, locator basemaps, chi/swath profiles, rose diagrams, and more.

- **Quick start and installation (pixi):** see the [README](https://github.com/LSDtopotools/lsdviztools3#readme).
- **Full usage guide with worked examples for every plot, in CLI and Python:**
  [Usage guide](documentation.md).

## Layout

```
docs/documentation.md   Full usage guide (the main reference)
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

For a single-file introduction to the CLI and Python APIs, open the
[Usage guide](documentation.md) and follow §1 (conventions) then §4 (maps).
