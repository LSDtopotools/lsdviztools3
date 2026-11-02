# Spectral analysis and valley spacing

Figures from `lsdtt-spectral-analysis`: the 2D and radial spectra of a tile,
significant frequencies on a polar plot, and valley spacing.

| Function | Module | What it draws |
|---|---|---|
| `plot_valley_spectrum`, `plot_spacing_vs_ci`, `plot_landscape_panel` | `charts.valley_spacing` | Valley spacing from `lsdtt-spectral-analysis` (`_valley_spacing.csv`, `_valley_spacing_spectrum.csv`): one landscape's spectrum with both spacing measures, spacing vs channelization index across a sweep, and a side-by-side hillshade panel. |
| `plot_significant_frequencies_polar` | `charts.spectral` | Significant 2D spectrum cells from `lsdtt-spectral-analysis` (`_significant_frequencies_2d.csv`, or the legacy `_significant_freq.csv`) on a North-up polar axis: bearing (wave vector, or crest with `orientation="crest"`) against frequency, rings labelled in wavelength, coloured by normalised PSD. A run with nothing significant draws empty axes with a note. |
| `plot_spectral_tile_summary`, `plot_normalised_spectrum_2d`, `plot_radial_spectrum` | `charts.spectral` | One `lsdtt-spectral-analysis` run on a 2 × 2 figure: the input raster (found from `_ingested_params.param`), the normalised 2D spectrum north-up with significant cells outlined, the binned radial spectrum with background, normalised PSD and rollover, and the polar plot. The two spectrum panels are also available on their own. Load the run with `io.spectral.load_spectral_analysis`; panels whose outputs were switched off show a note. |
