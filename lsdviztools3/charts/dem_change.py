"""Figures of potential landslides from a DEM of difference.

Draws the outputs of ``run_dem_change_landslides`` in ``lsdtt-slope-stability``.
Every statistic is read from the run's files; nothing is recomputed. Histogram
bins and subsampling are drawing choices.

* The error of the difference: its distribution, NMAD against slope, and the
  spatial variogram with lsdtt3's fit.
* Area coverage in the study domain, and its sensitivity to the threshold.
* Frequency - area of the patches.
* Patches found against those produced by noise alone, by size.
* Maps on a hillshade of the newer DEM: changed pixels and the difference.

The patches are candidates. Vegetation, snow, water, radar penetration and DEM
artefacts also change elevation between two products.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from ..io.slope_stability import load_frequency_area, load_quantity_table, load_table, read_band
from .slope_stability import _import_mpl, _map_style, _finish_map, plot_frequency_area, render_value_map

#: Loss and gain: the ends of ColorBrewer RdBu, safe for red-green colour blindness.
LOSS_COLOUR = "#b2182b"
GAIN_COLOUR = "#2166ac"
#: Results from noise alone: a neutral, so it never reads as a third kind of change.
NOISE_COLOUR = "#7f7f7f"
INK = "#333333"
GRID = "#e3e3e3"
DIFFERENCE_CMAP = "RdBu"


def _tidy(ax) -> None:
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _change_colour(name: str) -> str:
    return LOSS_COLOUR if name == "loss" else GAIN_COLOUR


def plot_dem_change_uncertainty(corrected: np.ndarray, summary: Dict[str, float], by_slope: pd.DataFrame,
                                variogram: pd.DataFrame, *, title: Optional[str] = None):
    """Three panels: the distribution of the difference, NMAD by slope, the variogram."""
    plt = _import_mpl()
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)

    ax = axes[0]
    values = corrected[np.isfinite(corrected)]
    nmad = summary["nmad_metres"]
    limit = 6.0 * nmad if nmad > 0 else max(abs(np.percentile(values, [1, 99])))
    # Values beyond the axis are left out of the bars and counted in the text.
    ax.hist(values, bins=121, range=(-limit, limit), color="#6d8fb3", linewidth=0)
    outside = 100.0 * np.mean(np.abs(values) > limit)
    for x, style, label in ((-1.96 * nmad, "--", "±1.96 NMAD"), (1.96 * nmad, "--", None)):
        ax.axvline(x, color=INK, lw=1.1, ls=style, label=label)
    ax.set_xlim(-limit, limit)
    ax.set_xlabel("difference, bias removed (m)")
    ax.set_ylabel("cells")
    source = "stable terrain" if summary.get("used_stable_mask", 0) else "all terrain, clipped at 4 NMAD"
    ax.set_title("error of the difference", fontsize=10.5)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.42)   # room for the text above the peak
    text = (f"error from {source}\nbias {summary['bias_metres']:+.2f} m, NMAD {nmad:.2f} m, "
            f"RMSE {summary['rmse_metres']:.2f} m\nn {int(summary['n_stable_cells']):,}; "
            f"{outside:.1f} % of cells beyond the axis")
    ax.text(0.02, 0.98, text, transform=ax.transAxes, va="top", fontsize=7.5, color=INK, zorder=5,
            bbox=dict(fc="white", ec="none", pad=1.5))
    ax.legend(fontsize=8, loc="center right", frameon=False)
    _tidy(ax)

    ax = axes[1]
    ax.plot(by_slope["slope_median"], by_slope["nmad_metres"], "-o", color=INK, lw=1.8, ms=4.5)
    ax.set_xlabel("slope (m/m)")
    ax.set_ylabel("NMAD (m)")
    ax.set_ylim(bottom=0)
    ax.set_title("error against slope", fontsize=10.5)
    _tidy(ax)

    ax = axes[2]
    if len(variogram):
        ax.scatter(variogram["lag_metres"], variogram["semivariance"], s=22, color=INK, linewidths=0,
                   label="empirical")
        if "fitted_semivariance" in variogram and variogram["fitted_semivariance"].notna().any():
            ax.plot(variogram["lag_metres"], variogram["fitted_semivariance"], color=LOSS_COLOUR, lw=1.8,
                    label="nugget + spherical")
            r = summary.get("correlation_range_metres", np.nan)
            if np.isfinite(r) and r > 0:
                ax.axvline(r, color=NOISE_COLOUR, lw=1.0, ls=":",
                           label=f"range {r:,.0f} m, nugget {summary.get('nugget_fraction', np.nan):.2f}")
        ax.set_xscale("log")
        ax.legend(fontsize=8, loc="lower right", frameon=False)
    ax.set_xlabel("distance (m)")
    ax.set_ylabel("standardised semivariance")
    ax.set_ylim(bottom=0)
    ax.set_title("spatial correlation of the error", fontsize=10.5)
    _tidy(ax)
    if title:
        fig.suptitle(title, fontsize=11.5)
    return fig, axes


def plot_dem_change_coverage(coverage: pd.DataFrame, sensitivity: pd.DataFrame,
                             false_detection: Optional[Dict[str, float]] = None,
                             *, chosen_sigma: Optional[float] = None, title: Optional[str] = None):
    """Share of the study domain that changed, and its dependence on the threshold."""
    plt = _import_mpl()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), constrained_layout=True)
    assessed = float(coverage["assessed_area_m2"].iloc[0])

    ax = axes[0]
    rows = coverage[coverage["n_patches"] > 0].reset_index(drop=True)
    names = rows["change"].tolist()
    area = rows["area_m2"].to_numpy(dtype=float)
    error = rows["area_uncertainty_m2"].to_numpy(dtype=float)
    counts = rows["n_patches"].to_numpy(dtype=int)
    colours = [_change_colour(c) for c in names]
    if len(rows) == 2:
        # The noise simulation counts every patch it finds, so it is compared
        # with the total, never with one sign.
        names = names + ["loss + gain"]
        area = np.append(area, area.sum())
        error = np.append(error, np.hypot(error[0], error[1]))
        counts = np.append(counts, counts.sum())
        colours = colours + ["#8a8a94"]
    x = np.arange(len(names))
    percent = 100.0 * area / assessed
    percent_error = 100.0 * error / assessed
    ax.bar(x, percent, width=0.58, yerr=percent_error, capsize=4, ecolor=INK, color=colours, edgecolor="white",
           linewidth=2)
    for k in range(len(names)):
        ax.annotate(f"{area[k] / 1e6:,.1f} km\u00b2\n{counts[k]:,} patches", (k, percent[k] + percent_error[k]),
                    xytext=(0, 4), textcoords="offset points", ha="center", fontsize=8, color=INK)
    if len(names) and false_detection is not None and np.isfinite(false_detection.get("mean_false_area_m2", np.nan)):
        noise = 100.0 * false_detection["mean_false_area_m2"] / assessed
        ax.hlines(noise, x[-1] - 0.36, x[-1] + 0.36, color="black", lw=2.0, ls="--",
                  label=f"from noise alone: {noise:.2f} %")
        ax.legend(fontsize=8, loc="upper left", frameon=False)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel(f"% of assessed area ({assessed / 1e6:,.0f} km\u00b2)")
    ax.set_ylim(0, (percent + percent_error).max() * 1.32 if len(names) else 1)
    ax.set_title("area coverage", fontsize=10.5)
    _tidy(ax)
    ax.grid(False, axis="x")

    ax = axes[1]
    s = sensitivity.sort_values("detection_sigma")
    ax.plot(s["detection_sigma"], 100.0 * s["fraction_of_assessed_area"], "-o", color=INK, lw=1.8, ms=5)
    if chosen_sigma is not None:
        ax.axvline(chosen_sigma, color=NOISE_COLOUR, lw=1.1, ls=":", label=f"used: {chosen_sigma:.2f}")
        ax.legend(fontsize=8, frameon=False)
    ax.set_yscale("log")
    ax.set_xlabel("detection threshold (multiples of the error)")
    ax.set_ylabel("% of assessed area detected")
    ax.set_title("sensitivity to the threshold", fontsize=10.5)
    _tidy(ax)
    if title:
        fig.suptitle(title, fontsize=11.5)
    return fig, axes


def plot_false_detection_by_area(by_area: pd.DataFrame, false_detection: Optional[Dict[str, float]] = None,
                                 *, title: Optional[str] = None, ax=None):
    """Patches found in each size class, against the mean number from noise alone."""
    plt = _import_mpl()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    else:
        fig = ax.figure
    rows = by_area[(by_area["n_observed"] > 0) | (by_area["mean_n_simulated"] > 0)]
    centre = np.sqrt(rows["area_lower_m2"] * rows["area_upper_m2"])
    ax.plot(centre, rows["n_observed"].where(rows["n_observed"] > 0), "-o", color=LOSS_COLOUR, lw=1.8, ms=4.5,
            label="patches found")
    ax.plot(centre, rows["mean_n_simulated"].where(rows["mean_n_simulated"] > 0), "--s", color=NOISE_COLOUR, lw=1.6,
            ms=4, label="mean from noise alone")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("patch area (m²)")
    ax.set_ylabel("patches per size class")
    ax.legend(fontsize=8, frameon=False)
    if false_detection is not None and np.isfinite(false_detection.get("area_false_discovery_ratio", np.nan)):
        ax.text(0.97, 0.80, f"noise alone could explain\n{100 * false_detection['area_false_discovery_ratio']:.0f} % "
                            "of the detected area\n(a lower bound)", transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color=INK)
    ax.set_title(title or "patches found against noise alone", fontsize=10.5)
    _tidy(ax)
    return fig, ax


def noise_limited_area(by_area: pd.DataFrame, max_ratio: float = 0.1) -> float:
    """Upper edge of the largest size class where noise gives over ``max_ratio`` of the patches found.

    Read from lsdtt3's table. NaN if no size class exceeds it.
    """
    rows = by_area[by_area["n_observed"] > 0]
    over = rows[pd.to_numeric(rows["false_discovery_ratio"], errors="coerce") > max_ratio]
    return float(over["area_upper_m2"].max()) if len(over) else float("nan")


def render_dem_change_map(dem_path, class_raster, png_path, *, coverage: Optional[pd.DataFrame] = None,
                          coordinates: str = "local", high_res: bool = False):
    """Changed pixels on a hillshade of the newer DEM."""
    from ..io.raster import load_raster
    from ..render.categorical import render_categorical_drape

    da = load_raster(class_raster)
    values = np.asarray(da.values, dtype=float)
    drape = da.copy(data=np.where(values < -0.5, 1.0, np.where(values > 0.5, 2.0, np.nan)))
    label = {"loss": "elevation loss", "gain": "elevation gain"}
    shares = {}
    if coverage is not None:
        shares = {str(r["change"]): 100.0 * float(r["fraction_of_assessed_area"]) for _, r in coverage.iterrows()}
    categories = {}
    for key, name, colour in ((1, "loss", LOSS_COLOUR), (2, "gain", GAIN_COLOUR)):
        if (drape.values == key).any():
            suffix = f" ({shares[name]:.1f} %)" if name in shares else ""
            categories[key] = (label[name] + suffix, colour)
    if not categories:
        return None
    fig = render_categorical_drape(dem_path, drape, categories=categories, coordinates=coordinates, style=_map_style(),
                                   transparency=5.0, legend_header="changed pixels", require_labels=False,
                                   decorate=False)
    return _finish_map(fig, dem_path, png_path, coordinates, None, high_res=high_res)


def write_dem_change_charts(f: Dict[str, Path], out_dir, stem: str, *, dpi: int = 300,
                            note: Callable[[str], None] = print) -> List[Path]:
    """Every chart whose inputs are in ``f`` (see :func:`lsdviztools3.io.slope_stability.present`)."""
    written: List[Path] = []
    if not any(k.startswith("dem_change") for k in f):
        return written
    plt = _import_mpl()
    out_dir = Path(out_dir)

    def save(fig, name):
        path = out_dir / f"{stem}_{name}.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(path)

    false_detection = load_quantity_table(f["dem_change_false_detection"]) if "dem_change_false_detection" in f else None

    need = ("dem_change_corrected", "dem_change_uncertainty", "dem_change_uncertainty_by_slope", "dem_change_variogram")
    if all(k in f for k in need):
        fig, _ = plot_dem_change_uncertainty(
            read_band(f["dem_change_corrected"]), load_quantity_table(f["dem_change_uncertainty"]),
            load_table(f["dem_change_uncertainty_by_slope"], ("slope_median", "nmad_metres")),
            load_table(f["dem_change_variogram"], ("lag_metres", "semivariance")),
            title="uncertainty of the DEM of difference")
        save(fig, "dem_change_uncertainty")
    else:
        note("no DEM of difference uncertainty chart: missing " + ", ".join(k for k in need if k not in f))

    if "dem_change_coverage" in f and "dem_change_sensitivity" in f:
        fig, _ = plot_dem_change_coverage(
            load_table(f["dem_change_coverage"], ("change", "n_patches", "area_m2", "area_uncertainty_m2",
                                                  "fraction_of_assessed_area", "assessed_area_m2")),
            load_table(f["dem_change_sensitivity"], ("detection_sigma", "fraction_of_assessed_area")), false_detection,
            title="potential landslides: area coverage")
        save(fig, "dem_change_coverage")

    if "dem_change_frequency_area" in f and "dem_change_frequency_area_fit" in f:
        density, fit = load_frequency_area(f["dem_change_frequency_area"], f["dem_change_frequency_area_fit"])
        fig, ax = plot_frequency_area(density, fit, title="potential landslides: frequency - area",
                                      colours={"loss": LOSS_COLOUR, "gain": GAIN_COLOUR})
        if "dem_change_false_detection_by_area" in f:
            limit = noise_limited_area(load_table(f["dem_change_false_detection_by_area"],
                                                  ("area_upper_m2", "n_observed", "false_discovery_ratio")))
            if np.isfinite(limit):
                left = ax.get_xlim()[0]
                ax.axvspan(left, limit, color=NOISE_COLOUR, alpha=0.14, lw=0, zorder=0,
                           label=f"noise > 10 % of patches (A < {limit:,.0f} m\u00b2)")
                ax.set_xlim(left=left)
                ax.legend(fontsize=7.5, loc="lower left")
        save(fig, "dem_change_frequency_area")
    else:
        note("no DEM of difference frequency - area chart: no patches, or missing tables")

    if "dem_change_false_detection_by_area" in f:
        by_area = load_table(f["dem_change_false_detection_by_area"], ("area_lower_m2", "area_upper_m2", "n_observed",
                                                                       "mean_n_simulated"))
        fig, _ = plot_false_detection_by_area(by_area, false_detection)
        save(fig, "dem_change_false_detection")
    return written


def write_dem_change_maps(f: Dict[str, Path], dem, out_dir, stem: str, *, coordinates: str = "local",
                          value_range: Optional[float] = None) -> List[Path]:
    """The maps on a hillshade of ``dem``, the newer DEM."""
    written: List[Path] = []
    out_dir = Path(out_dir)
    if "dem_change_class" in f:
        coverage = load_table(f["dem_change_coverage"], ("change", "fraction_of_assessed_area")) \
            if "dem_change_coverage" in f else None
        path = render_dem_change_map(dem, f["dem_change_class"], out_dir / f"{stem}_dem_change_map.png",
                                     coverage=coverage, coordinates=coordinates)
        if path is not None:
            written.append(path)
    if "dem_change_corrected" in f:
        limit = value_range
        if limit is None and "dem_change_uncertainty" in f:
            limit = 5.0 * load_quantity_table(f["dem_change_uncertainty"])["nmad_metres"]
        if limit is None or not limit > 0:
            values = read_band(f["dem_change_corrected"])
            limit = float(np.nanpercentile(np.abs(values), 99))
        written.append(render_value_map(dem, f["dem_change_corrected"], out_dir / f"{stem}_dem_difference_map.png",
                                        cmap=DIFFERENCE_CMAP, label="elevation change, bias removed (m)", vmin=-limit,
                                        vmax=limit, title="DEM of difference (newer - older)", coordinates=coordinates))
    return written
