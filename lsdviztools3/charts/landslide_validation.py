"""Figures comparing slope-stability outputs with a mapped landslide inventory.

Draws tables of scores (ROC curves, AUC, frequency ratio, precision and recall,
and recall by landslide size) and maps where each method agrees with the
inventory. Every score is read from the tables; only the map categories are
derived, by comparing the rasters cell by cell.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

#: Okabe-Ito colours, safe for colour-blind readers. A method keeps its colour in every figure.
METHOD_COLOURS: Dict[str, str] = {
    "hazard": "#0072B2",
    "scalr": "#D55E00",
    "dod": "#009E73",
    "baseline": "#7f7f7f",
}
HIT_COLOUR = "#0072B2"
FALSE_COLOUR = "#E69F00"
MISS_COLOUR = "#CC79A7"
INK = "#333333"
GRID = "#e3e3e3"

AUC_COLUMNS = ("predictor", "auc", "ci95_low", "ci95_high")
ROC_COLUMNS = ("predictor", "false_positive_rate", "true_positive_rate")
RATIO_COLUMNS = ("predictor", "x_median", "frequency_ratio", "random_low", "random_high")
BINARY_COLUMNS = ("detector", "precision", "recall", "lift", "slope_matched_lift", "fraction_of_domain")
POLYGON_COLUMNS = ("area_low_m2", "area_high_m2", "n_polygons", "any_cell_detected", "half_detected")
PATCH_COLUMNS = ("change", "area_low_m2", "area_high_m2", "n_patches", "majority_inside")
SWEEP_COLUMNS = ("initiation_slope", "stopping_slope", "auc")
CONDITIONAL_COLUMNS = ("is_zero", "ah_low", "ah_high", "ah_centre", "n_cells", "fraction_of_domain", "ratio",
                       "ip_low", "ip_high", "significant")
FREQUENCY_FIT_COLUMNS = ("catalogue", "n_patches", "area_min_m2", "eta", "eta_low", "eta_high", "n_tail")
FREQUENCY_DENSITY_COLUMNS = ("catalogue", "area_centre_m2", "probability_density_per_m2", "n_landslides")

TABLES: Dict[str, Tuple[str, Sequence[str]]] = {
    "auc": ("validation_auc.csv", AUC_COLUMNS),
    "roc": ("validation_roc.csv", ROC_COLUMNS),
    "frequency_ratio": ("validation_frequency_ratio.csv", RATIO_COLUMNS),
    "binary": ("validation_binary.csv", BINARY_COLUMNS),
    "polygon_recall": ("validation_dod_polygon_recall.csv", POLYGON_COLUMNS),
    "patch_precision": ("validation_dod_patch_precision.csv", PATCH_COLUMNS),
    "sweep": ("validation_hazard_sweep.csv", SWEEP_COLUMNS),
    "hazard_conditional": ("validation_hazard_conditional.csv", CONDITIONAL_COLUMNS),
    "frequency_fit": ("validation_frequency_area.csv", FREQUENCY_FIT_COLUMNS),
    "frequency_density": ("validation_frequency_area_density.csv", FREQUENCY_DENSITY_COLUMNS),
}


def load_validation_tables(directory) -> Dict[str, pd.DataFrame]:
    """Read the tables present in ``directory``; a table missing a column is an error."""
    out: Dict[str, pd.DataFrame] = {}
    for key, (name, columns) in TABLES.items():
        path = Path(directory) / name
        if not path.exists():
            continue
        table = pd.read_csv(path)
        missing = [c for c in columns if c not in table.columns]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        out[key] = table
    return out


def method_of(label: str) -> str:
    """The method a predictor or detector label belongs to."""
    low = label.lower()
    if low.startswith("hazard"):
        return "hazard"
    if low.startswith("scalr"):
        return "scalr"
    if low.startswith("dod"):
        return "dod"
    return "baseline"


def _import_mpl():
    import matplotlib.pyplot as plt

    return plt


def _tidy(ax) -> None:
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _size_labels(table: pd.DataFrame) -> list:
    def text(v):
        if v >= 1e8:
            return "∞"
        return f"{v / 1e4:g}"

    return [f"{text(a)}–{text(b)}" for a, b in zip(table.area_low_m2, table.area_high_m2)]


def plot_skill(auc: pd.DataFrame, roc: pd.DataFrame, *, title: Optional[str] = None):
    """ROC curves, and each predictor's AUC with its 95 % interval."""
    plt = _import_mpl()
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.4, 5.0), constrained_layout=True,
                                      gridspec_kw={"width_ratios": [1.0, 1.15]})
    dashes = {"hazard": "-", "scalr": "-", "dod": "-", "baseline": "--"}
    shown = [p for p in auc.predictor if p in set(roc.predictor) and "|z|" not in p and "elevation loss" not in p]
    for name in shown:
        part = roc[roc.predictor == name]
        value = float(auc.loc[auc.predictor == name, "auc"].iloc[0])
        method = method_of(name)
        left.plot(part.false_positive_rate, part.true_positive_rate, lw=2, ls=dashes[method],
                  color="#b0b0b0" if name.lower().startswith("upslope") else METHOD_COLOURS[method],
                  label=f"{name} ({value:.2f})")
    left.plot([0, 1], [0, 1], color=INK, lw=0.8, ls=":", label="no skill (0.50)")
    left.set(xlim=(0, 1), ylim=(0, 1), xlabel="False positive rate", ylabel="True positive rate", aspect="equal",
             title="ROC curves")
    left.legend(frameon=False, fontsize=7.5, loc="lower right")
    _tidy(left)

    ordered = auc.sort_values("auc").reset_index(drop=True)
    y = np.arange(len(ordered))
    colours = [METHOD_COLOURS[method_of(p)] for p in ordered.predictor]
    right.hlines(y, ordered.ci95_low, ordered.ci95_high, colors=colours, lw=2)
    right.scatter(ordered.auc, y, s=46, c=colours, edgecolors="white", linewidths=1.2, zorder=3)
    for k, row in ordered.iterrows():
        right.text(row.ci95_high + 0.006, k, f"{row.auc:.3f}", va="center", fontsize=8, color=INK)
    right.axvline(0.5, color=INK, lw=0.8, ls=":")
    right.set_yticks(y, ordered.predictor, fontsize=8)
    right.set_xlim(0.48, max(0.8, float(ordered.ci95_high.max()) + 0.05))
    right.set(xlabel="AUC (95 % interval, 3 km block bootstrap)", title="Skill at ranking landslide cells")
    _tidy(right)
    if title:
        fig.suptitle(title)
    return fig, (left, right)


def plot_frequency_ratio(table: pd.DataFrame, *, title: Optional[str] = None):
    """Landslide probability relative to the study-area average, along each predictor."""
    plt = _import_mpl()
    names = list(dict.fromkeys(table.predictor))
    fig, axes = plt.subplots(1, len(names), figsize=(3.3 * len(names), 3.7), sharey=True, constrained_layout=True,
                             squeeze=False)
    labels = {"Hazard area": "Hazard area (m²/m)", "Local slope": "Gradient (m/m)",
              "scaLr factor of safety": "Factor of safety", "DoD |change|": "|Elevation change| (m)"}
    for ax, name in zip(axes[0], names):
        part = table[table.predictor == name].sort_values("x_median")
        method = "hazard" if name == "Hazard area" else "scalr" if name.startswith("scaLr") else \
            "dod" if name.startswith("DoD") else "baseline"
        x = part.x_median.to_numpy(float)
        if name == "Hazard area":
            x = np.maximum(x, 0.01)  # Ah = 0 is drawn at the left edge of the log axis
            ax.set_xscale("log")
        ax.fill_between(x, part.random_low, part.random_high, color="#d9d9d9", lw=0, label="expected by chance")
        ax.axhline(1.0, color=INK, lw=0.8)
        ax.plot(x, part.frequency_ratio, "-o", ms=4, lw=2, color=METHOD_COLOURS[method], markeredgecolor="white",
                markeredgewidth=0.6)
        ax.set(xlabel=labels.get(name, name) + (", 0 drawn at 0.01" if name == "Hazard area" else ""), title=name)
        _tidy(ax)
    axes[0][0].set_ylabel("P(landslide | x) / P(landslide)")
    axes[0][0].legend(frameon=False, fontsize=8, loc="upper left")
    if title:
        fig.suptitle(title)
    return fig, axes


def plot_hazard_conditional(table: pd.DataFrame, *, reference_area: Optional[float] = 60.0,
                            title: Optional[str] = None):
    """Normalised conditional probability against hazard area, after Milledge et al. (2019, Fig. 6).

    The upper panel is P(L | Ah) / P(L) on a log axis, with a separate point for
    Ah = 0 and an asterisk where the value lies outside the interval of their
    Eq. 5. The lower panel is the share of the study area in each bin.
    """
    plt = _import_mpl()
    zero = table[table.is_zero.astype(bool)]
    bins = table[~table.is_zero.astype(bool)].sort_values("ah_centre")
    if bins.empty:
        raise ValueError("the conditional-probability table has no bins with hazard area above zero")
    # The lowest bin can start at a vanishingly small area: draw it one decade wide.
    low = bins.ah_low.to_numpy(float, copy=True)
    high = bins.ah_high.to_numpy(float)
    low[0] = max(low[0], high[0] / 10.0)
    x = np.clip(bins.ah_centre.to_numpy(float), low, high)
    x_zero = low[0] / 4.0
    fig, (ax, hist) = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True, height_ratios=(3.2, 1), constrained_layout=True)
    colour = METHOD_COLOURS["hazard"]
    ax.fill_between(x, bins.ip_low, bins.ip_high, color="#d9d9d9", lw=0, label="Interval for random landslide cells (Eq. 5)")
    ax.axhline(1.0, color=INK, lw=0.8)
    ax.plot(x, bins.ratio, "-o", color=colour, lw=2, ms=5, markeredgecolor="white", markeredgewidth=0.6,
            label="Hazard area > 0")
    flagged = bins.significant.astype(bool).to_numpy()
    ax.plot(x[flagged], bins.ratio.to_numpy(float)[flagged] * 1.13, ls="none", marker="*", ms=5, color=INK,
            label="Outside the interval")
    if not zero.empty:
        row = zero.iloc[0]
        ax.errorbar([x_zero], [row.ratio], yerr=[[max(row.ratio - row.ratio * row.ip_low, 0)],
                                                 [row.ratio * row.ip_high - row.ratio]],
                    fmt="s", ms=7, color=colour, markeredgecolor="white", label="Hazard area = 0")
        if bool(row.significant):
            ax.plot([x_zero], [row.ratio * 1.18], marker="*", ms=5, color=INK, ls="none")
    if reference_area is not None:
        for a in (ax, hist):
            a.axvline(reference_area, color=INK, lw=0.8, ls=":")
        ax.text(reference_area * 1.08, 0.03, f"{reference_area:g} m²/m", transform=ax.get_xaxis_transform(), va="bottom",
                fontsize=8, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("P(L | hazard area) / P(L)")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    _tidy(ax)
    edges_low, edges_high = low, high
    hist.bar(edges_low, bins.fraction_of_domain, width=edges_high - edges_low, align="edge", color="#b2182b",
             edgecolor="white", lw=0.5)
    if not zero.empty:
        hist.bar([x_zero], zero.fraction_of_domain, width=x_zero * 0.7, color="#b2182b", edgecolor="white", lw=0.5)
    hist.set(xlabel="Hazard area (m²/m)", ylabel="Share of area")
    _tidy(hist)
    if not zero.empty:
        ticks = [x_zero] + [t for t in (0.1, 1, 10, 100, 1000, 10000) if x_zero < t <= edges_high[-1]]
        hist.set_xticks(ticks)
        hist.set_xticklabels(["0"] + [f"{t:g}" for t in ticks[1:]])
        hist.minorticks_off()
    hist.set_xlim(x_zero * 0.5, edges_high[-1] * 1.2)
    if title:
        fig.suptitle(title)
    return fig, (ax, hist)


def plot_frequency_area_comparison(fits: pd.DataFrame, density: pd.DataFrame, *, area_min_m2: Optional[float] = None,
                                   curves: Optional[Sequence[str]] = None, title: Optional[str] = None):
    """Frequency-area density of each catalogue, and its tail exponent above one common area.

    ``curves`` names the catalogues drawn on the density panel (default: all).
    """
    plt = _import_mpl()
    if area_min_m2 is None:
        area_min_m2 = float(fits.area_min_m2.max())
    part = fits[np.isclose(fits.area_min_m2, area_min_m2)].reset_index(drop=True)
    if part.empty:
        raise ValueError(f"no fits at area_min_m2 = {area_min_m2:g}")
    names = list(part.catalogue) if curves is None else list(curves)
    unknown = [n for n in names if n not in set(part.catalogue)]
    if unknown:
        raise ValueError(f"catalogues not in the table: {unknown}")
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.4, max(4.6, 0.32 * len(part) + 1.6)), constrained_layout=True,
                                      width_ratios=(1, 1.05))
    styles = ("-", "--", ":", "-.")
    used: Dict[str, int] = {}
    for name in names:
        method = method_of(name)
        d = density[(density.catalogue == name) & (density.n_landslides > 0)].sort_values("area_centre_m2")
        if d.empty:
            continue
        k = used.get(method, 0)
        used[method] = k + 1
        left.plot(d.area_centre_m2, d.probability_density_per_m2, styles[k % len(styles)], lw=1.8,
                  color=METHOD_COLOURS[method], label=name)
    left.axvline(area_min_m2, color=INK, lw=0.8, ls=":")
    left.set(xscale="log", yscale="log", xlabel="Patch area (m²)", ylabel="Probability density (per m²)")
    left.legend(frameon=False, fontsize=7)
    _tidy(left)
    y = np.arange(len(part))[::-1]
    for yy, row in zip(y, part.itertuples()):
        colour = METHOD_COLOURS[method_of(row.catalogue)]
        right.plot([row.eta_low, row.eta_high], [yy, yy], color=colour, lw=2)
        right.plot([row.eta], [yy], "o", ms=6, color=colour, markeredgecolor="white")
        right.text(1.01, yy, f"n = {int(row.n_tail)}", transform=right.get_yaxis_transform(), va="center", fontsize=7,
                   color=INK)
    right.set_yticks(y)
    right.set_yticklabels(part.catalogue, fontsize=8)
    right.set_xlabel(f"Exponent η of p(A) ~ A^η for A ≥ {area_min_m2 / 1e4:g} ha (95% bootstrap)")
    _tidy(right)
    if title:
        fig.suptitle(title)
    return fig, (left, right)


def plot_detector_scores(binary: pd.DataFrame, *, base_rate: Optional[float] = None, title: Optional[str] = None):
    """Precision, recall and lift of each binary detector."""
    plt = _import_mpl()
    table = binary.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(table))
    colours = [METHOD_COLOURS[method_of(d)] for d in table.detector]
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 0.42 * len(table) + 1.7), sharey=True, constrained_layout=True)
    panels = (("precision", "Precision: detected area in the inventory", (0, 1)),
              ("recall", "Recall: inventory area detected", (0, 1)))
    for ax, (column, label, limits) in zip(axes[:2], panels):
        ax.barh(y, table[column], color=colours, height=0.62)
        for k, v in enumerate(table[column]):
            ax.text(v + 0.012, k, f"{v:.2f}", va="center", fontsize=8, color=INK)
        ax.set(xlim=limits, xlabel=label)
        _tidy(ax)
    if base_rate is not None:
        axes[0].axvline(base_rate, color=INK, ls=":", lw=1)
        axes[0].text(base_rate + 0.01, len(table) - 0.45, f"inventory share {base_rate:.2f}", fontsize=7.5, color=INK)
    ax = axes[2]
    ax.hlines(y, table.slope_matched_lift, table.lift, color="#c8c8c8", lw=2)
    ax.scatter(table.lift, y, s=42, c=colours, label="lift", zorder=3)
    ax.scatter(table.slope_matched_lift, y, s=42, facecolors="white", edgecolors=colours, linewidths=1.5,
               label="lift after matching slope", zorder=3)
    ax.axvline(1.0, color=INK, lw=0.8, ls=":")
    ax.set(xlabel="Precision ÷ chance")
    ax.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    _tidy(ax)
    axes[0].set_yticks(y, table.detector, fontsize=8)
    if title:
        fig.suptitle(title)
    return fig, axes


def plot_dem_change_by_size(polygon_recall: pd.DataFrame, patch_precision: pd.DataFrame,
                            *, title: Optional[str] = None):
    """How much of each inventory size class the DoD finds, and how many patches are in the inventory."""
    plt = _import_mpl()
    from .dem_change import GAIN_COLOUR, LOSS_COLOUR

    fig, (left, right) = plt.subplots(1, 2, figsize=(11.6, 4.3), constrained_layout=True)
    x = np.arange(len(polygon_recall))
    left.bar(x - 0.2, polygon_recall.any_cell_detected, width=0.38, color=METHOD_COLOURS["dod"], label="any cell changed")
    left.bar(x + 0.2, polygon_recall.half_detected, width=0.38, color="#9ad5c1", label="at least half changed")
    for k, n in enumerate(polygon_recall.n_polygons):
        left.text(k, 1.015, f"n={n}", ha="center", fontsize=7, color=INK)
    left.set_xticks(x, _size_labels(polygon_recall), fontsize=8)
    left.set(ylim=(0, 1.08), xlabel="Mapped landslide area (ha)", ylabel="Fraction of mapped landslides",
             title="Inventory landslides found by the DoD")
    left.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0, 0.95))
    _tidy(left)

    bins = patch_precision[["area_low_m2", "area_high_m2"]].drop_duplicates().sort_values("area_low_m2")
    labels = _size_labels(bins)
    xr = np.arange(len(bins))
    for change, colour, offset in (("loss", LOSS_COLOUR, -0.2), ("gain", GAIN_COLOUR, 0.2)):
        part = patch_precision[patch_precision.change == change].set_index("area_low_m2")
        values = [float(part.majority_inside.get(v, np.nan)) for v in bins.area_low_m2]
        right.bar(xr + offset, values, width=0.38, color=colour, label=change)
    right.set_xticks(xr, labels, fontsize=8)
    right.set(ylim=(0, 1), xlabel="DoD patch area (ha)", ylabel="Fraction of patches mostly inside the inventory",
              title="DoD patches that match the inventory")
    right.legend(frameon=False, fontsize=8)
    _tidy(right)
    if title:
        fig.suptitle(title)
    return fig, (left, right)


def plot_hazard_sweep(sweep: pd.DataFrame, *, published: Tuple[float, float] = (0.839, 0.176),
                      title: Optional[str] = None):
    """AUC of the hazard area over the two thresholds."""
    plt = _import_mpl()
    grid = sweep.pivot(index="initiation_slope", columns="stopping_slope", values="auc").sort_index()
    fig, ax = plt.subplots(figsize=(5.8, 4.4), constrained_layout=True)
    image = ax.imshow(grid.to_numpy(), origin="lower", cmap="Blues", aspect="auto")
    ax.set_xticks(range(grid.shape[1]), [f"{v:g}" for v in grid.columns])
    ax.set_yticks(range(grid.shape[0]), [f"{v:g}" for v in grid.index])
    middle = float(np.nanmean([np.nanmin(grid.to_numpy()), np.nanmax(grid.to_numpy())]))
    for r in range(grid.shape[0]):
        for c in range(grid.shape[1]):
            v = grid.iat[r, c]
            ax.text(c, r, f"{v:.3f}", ha="center", va="center", fontsize=8, color="white" if v > middle else INK)
    if published[0] in grid.index and published[1] in grid.columns:
        ax.add_patch(plt.Rectangle((list(grid.columns).index(published[1]) - 0.5,
                                    list(grid.index).index(published[0]) - 0.5), 1, 1, fill=False, lw=2, ec=INK))
    ax.set(xlabel="Stopping gradient (m/m)", ylabel="Initiation gradient (m/m)",
           title=title or "Hazard area AUC (box: published values)")
    fig.colorbar(image, ax=ax, label="AUC", shrink=0.85)
    return fig, ax


def agreement_classes(predicted: np.ndarray, inventory: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """0 neither, 1 hit, 2 false alarm, 3 missed; NaN outside ``valid``."""
    out = np.full(predicted.shape, np.nan)
    p, m = predicted.astype(bool), inventory.astype(bool)
    out[valid] = 0
    out[valid & p & m] = 1
    out[valid & p & ~m] = 2
    out[valid & ~p & m] = 3
    return out


def plot_agreement_maps(dem: np.ndarray, inventory: np.ndarray, detections: Mapping[str, np.ndarray],
                        scores: Optional[Mapping[str, Tuple[float, float]]] = None, *, cell_size: float = 30.0,
                        title: Optional[str] = None):
    """A hillshade per method with hits, false alarms and missed inventory cells."""
    plt = _import_mpl()
    from matplotlib.colors import LightSource, ListedColormap
    from matplotlib.patches import Patch

    n = len(detections)
    aspect = dem.shape[0] / dem.shape[1]
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 4.6 * aspect + 0.9), constrained_layout=True, squeeze=False)
    filled = np.where(np.isfinite(dem), dem, np.nanmin(dem))
    shade = LightSource(azdeg=315, altdeg=45).hillshade(filled, dx=cell_size, dy=cell_size)
    valid = np.isfinite(dem)
    colours = ListedColormap([HIT_COLOUR, FALSE_COLOUR, MISS_COLOUR])
    extent = (0, dem.shape[1] * cell_size / 1000, 0, dem.shape[0] * cell_size / 1000)
    for ax, (name, predicted) in zip(axes[0], detections.items()):
        cls = agreement_classes(predicted, inventory, valid)
        ax.imshow(np.ma.masked_where(~valid, shade), cmap="gray", vmin=0, vmax=1, extent=extent, interpolation="bilinear")
        overlay = np.ma.masked_where(~np.isin(cls, [1, 2, 3]), cls)
        ax.imshow(overlay, cmap=colours, vmin=0.5, vmax=3.5, extent=extent, interpolation="nearest", alpha=0.9)
        heading = name
        if scores and name in scores:
            heading += f"\nprecision {scores[name][0]:.2f}, recall {scores[name][1]:.2f}"
        ax.set_title(heading, fontsize=9)
        ax.set_xlabel("Easting from window edge (km)", fontsize=8)
        ax.tick_params(labelsize=7)
    axes[0][0].set_ylabel("Northing from window edge (km)", fontsize=8)
    handles = [Patch(color=HIT_COLOUR, label="in both"), Patch(color=FALSE_COLOUR, label="detected, not in inventory"),
               Patch(color=MISS_COLOUR, label="in inventory, not detected")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=9)
    if title:
        fig.suptitle(title)
    return fig, axes
