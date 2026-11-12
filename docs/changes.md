# What's new

## Unreleased

- Fixed: in the vertical box plots of `lsdviz3-ksn-by-lithology` (node and
  reach k<sub>sn</sub>), the "(low sample: few qualifying reaches)" note under
  each low-sample unit ran into its neighbours' labels. The labels now say
  "(low sample)", and the caption explains it once.

## 0.4.2

**Colouring in classes.** Rasters, lines and points can be coloured in a few
classes instead of a continuous scale, with a class legend instead of the
colourbar: give the breaks (a list or array), or choose `equal_interval`,
`quantile` or `natural_breaks` (Jenks). `vmin`/`vmax` set the ends of the
colour gradient, each class takes the colour at its centre, and qualitative
colour maps (`tab10`, `Set1`, ...) give distinct colours. In Python:
`classes=` on `render_drape`, `render_channels` and `render_points`; on the
command line: `--classes` (`lsdviz3-drape`, `lsdviz3-channels`) and
`--ksn-classes` (`lsdviz3-chi`, `lsdviz3-knickpoints`). See
[Colouring in classes](basic_plots.md#colouring-in-classes).

**Locator maps.**

- A `regional` mode (`render_regional_basemap`, `lsdviz3-basemap --mode
  regional`): the region around the raster with coastlines, country and state
  borders, rivers, lakes, optional shaded relief, a scale bar and the
  footprint. `auto` now chooses it for small rasters, instead of the aerial
  close-up.
- Insets: `add_locator_inset`, and `--locator global|regional` on
  `lsdviz3-hillshade`, `-drape`, `-basins`, `-channels`, `-chi` and
  `-knickpoints`.
- Fixed: `global` mode failed to save a PNG (Ghostscript error) because of
  its north arrow, which a globe no longer has.

See [Locator maps](basic_plots.md#locator-maps-where-is-this-raster).

**Other changes**

- The `lsdviz3-*` commands always use matplotlib's `Agg` backend. Run from a
  notebook (`!lsdviz3-...` in Colab) they used to inherit the notebook's
  `MPLBACKEND` and fail to start.
- `lsdviz3-basins` no longer writes `scaling=stone, muon=BraucherBorchers,
  atmosphere=ERA40` across every basin map. The cosmogenic labels appear only
  when given (`--scaling-scheme`, `--muon-scheme`, `--atmosphere`), and with
  `--plot-title` in the map's top-left corner.
- `lsdviz3-channels` takes `--vmin`, `--vmax`, `--colorbar`, `--legend-title`
  and `--linear` (colour by the value rather than its log10).
- `render_points` takes `vmin` and `vmax`.
- `render.decorations.add_inset` is implemented (it raised
  `NotImplementedError`).

## 0.4.1

The first release on [PyPI](https://pypi.org/project/lsdviztools3/):
`pip install lsdviztools3`, for Python 3.12 or later. Without GMT, the map
commands stop with a message saying how to install it, and the matplotlib
charts still work.
