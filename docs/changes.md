# What's new

## Unreleased

- New: `lsdviz3-slope-stability` draws an lsdtt-slope-stability run: SINMAP's
  slope - specific catchment area plot and class map, the landslide and
  slip-body catalogues (area-volume, area-depth, area and length/width
  distributions), rupture against slip-surface depth, and maps of depth,
  factor of safety and probability of failure. The frequency – area plot
  draws lsdtt3's power-law tail fit and its exponent η. See
  [Slope stability](slope_stability.md).

## 0.4.3

**Knickpoints (`lsdviz3-knickpoints`).** See [Knickpoints](chi_analysis.md#knickpoints-lsdviz3-knickpoints).

- **Changed defaults:** major knickpoints are ranked by lsdtt3's split
  hierarchy (`--measure split_normalised`), and the threshold is lsdtt3's own
  choice (`--threshold major`). A table from an older lsdtt3 without the
  split columns still uses `delta` and `p95`. The knickpoints are always the
  stage 1 ones; `--stepaware-ksn` colours the k_sn of the step-aware fit
  under them.
- **New plots:**
  - `_knickpoints_by_elevation.png`;
  - `_knickpoint_order_kde.png`, the 1st to 4th strongest knickpoint of each
    channel, shallowing and steepening (`--max-order`, `--order-majors-only`);
  - `_step_long_profiles.png`, steps on the longitudinal profiles, green at a
    basin outlet or tributary mouth.
- **New overlays:** knickzones on the map (`--no-map-knickzones` to leave
  them off, `--profile-knickzones` to add them to the chi profile,
  `--no-knickzones`). Steps can be filtered with `--min-step-drop` and
  `--hide-flagged-steps`, or left off with `--no-steps`.
- `--basins` limits every figure and CSV to those basins.

**Chi plots (`lsdviz3-chi`, and the profiles of `lsdviz3-knickpoints`).**

- **Choosing what to plot:**
  - `--basins` picks the basins and their order;
  - `--basins-per-row` wraps the stacked plot onto rows;
  - `--ksn-gradient-key` adds a key of chi–elevation gradients at the plot's
    scale;
  - `--longest-channels N` keeps the main stem and the N-1 longest
    tributaries of each basin.

  See [Stacked plot layout](chi_analysis.md#stacked-plot-layout) and
  [Removing tributaries](chi_analysis.md#removing-tributaries).
- **New plot:** `_ksn_distribution.png`, the k_sn distribution per channel
  pixel and per fitted segment (`--no-ksn-distribution` to skip it).
  `--ksn-class-basis` chooses which one a class scheme is found in, and
  `--ksn-class-min` / `--ksn-class-max` bound the scheme (defaults 5 and 250).
  See [k_sn distributions](chi_analysis.md#k_sn-distributions).
- **Fit choice:** a run with step detection is drawn from its step-aware
  (stage 2) fit; `--stage1` draws the ordinary fit.

- New: `lsdviz3-slope-area` draws an lsdtt-slope-area run: binned gradient
  against area per polyfit lambda or D8 window, a comparison of two runs, and
  the low-concavity theta histogram, slope-area cloud and segments map. See
  [Slope-area and low-concavity reaches](channel_networks.md#slope-area-and-low-concavity-reaches-lsdviz3-slope-area).
  `load_slope_area` moved to `lsdviztools3.io.slope_area` (still importable
  from `lsdviztools3.io.ridgeline`).
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
