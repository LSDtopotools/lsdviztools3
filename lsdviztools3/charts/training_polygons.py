"""HAND distribution of each training polygon (lsdtt3 fluvial-surfaces, plan 7i).

lsdtt3 ``lsdtt-fluvial-surfaces`` with ``training_surfaces_fname`` writes
``<stem>_training_polygons.csv`` (polygon_id, label, class, terrace_number,
n_cells, n_with_field, median_m, p10_m, p25_m, p75_m, p90_m, iqr_m, flagged)
and ``<stem>_training_fit.csv``. This draws one box per polygon:

  box       p25 .. p75 - the core lsdtt3 chooses thresholds from;
  line      the median;
  whiskers  p10 .. p90;
  colour    floodplain or terrace; a hatched box is a polygon flagged as not a
            single level (left out of the terrace breaks).

The thresholds fitted over all training (the ``global`` rows of the fit CSV)
are drawn as horizontal lines. Visualisation only: every number is read from
lsdtt3's CSVs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd

_COLOURS = {"floodplain": "#3b7dd8", "terrace": "#d8833b"}


def plot_training_polygon_hand(polygons_csv, out_png: Optional[str] = None,
                               fit_csv=None, *, order: str = "median",
                               log_y: bool = False, figsize=None):
    """Box-and-whisker of HAND per training polygon.

    Args:
        polygons_csv: lsdtt3 ``<stem>_training_polygons.csv``.
        out_png: Optional output path; the figure is also returned.
        fit_csv: Optional ``<stem>_training_fit.csv``; its global thresholds
            are drawn as horizontal lines. Default: beside polygons_csv.
        order: ``median`` (ascending HAND) or ``id`` (file order).
        log_y: Log HAND axis from 0.1 m, for tiles where terraces sit hundreds
            of metres up.
        figsize: Inches; default widens with the polygon count.

    Returns:
        The matplotlib ``Figure``.
    """
    polygons_csv = Path(polygons_csv)
    p = pd.read_csv(polygons_csv).dropna(subset=["median_m"])
    if order == "median":
        p = p.sort_values("median_m")
    if fit_csv is None:
        cand = Path(str(polygons_csv).replace("_training_polygons.csv",
                                              "_training_fit.csv"))
        fit_csv = cand if cand.exists() else None

    n = len(p)
    fig, ax = plt.subplots(figsize=figsize or (max(5.0, 0.45 * n + 2.0), 4.5))
    for x, (_, r) in enumerate(p.iterrows()):
        c = _COLOURS.get(r["class"], "grey")
        ax.plot([x, x], [r["p10_m"], r["p25_m"]], color=c, lw=1)
        ax.plot([x, x], [r["p75_m"], r["p90_m"]], color=c, lw=1)
        ax.add_patch(Rectangle((x - 0.32, r["p25_m"]), 0.64,
                               r["p75_m"] - r["p25_m"], facecolor=c,
                               alpha=0.35 if r["flagged"] else 0.6,
                               edgecolor=c, hatch="///" if r["flagged"] else None))
        ax.plot([x - 0.32, x + 0.32], [r["median_m"]] * 2, color="black", lw=1.2)

    handles = [Patch(facecolor=_COLOURS["floodplain"], alpha=0.6, label="floodplain"),
               Patch(facecolor=_COLOURS["terrace"], alpha=0.6, label="terrace"),
               Patch(facecolor="white", edgecolor="grey", hatch="///",
                     label="flagged: not one level")]
    if fit_csv is not None:
        f = pd.read_csv(fit_csv)
        g = f[(f["scope"] == "global") & f["threshold_m"].notna()]
        g = g.drop_duplicates("threshold")
        # Many thresholds: lines only, labels would pile up.
        label_lines = len(g) <= 6
        for _, r in g.iterrows():
            style = "-" if r["threshold"] == "floodplain_bound" else "--"
            ax.axhline(r["threshold_m"], color="0.3", lw=0.8, ls=style)
            if label_lines:
                ax.annotate(f"{r['threshold']} {r['threshold_m']:.1f} m",
                            (n - 0.5, r["threshold_m"]), fontsize=7, ha="right",
                            va="bottom", color="0.3")
        handles += [plt.Line2D([], [], color="0.3", lw=0.8, ls="-",
                               label="floodplain / terrace bound"),
                    plt.Line2D([], [], color="0.3", lw=0.8, ls="--",
                               label="terrace breaks (all training)")]

    ax.set_xticks(range(n))
    ax.set_xticklabels([f"{int(i)}: {l}" for i, l in zip(p["polygon_id"], p["label"])],
                       rotation=60, ha="right", fontsize=8)
    ax.set_xlim(-0.6, n - 0.4)
    if log_y:
        # HAND 0 (on the channel) sits on the lower edge.
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.1)
    ax.set_ylabel("HAND (m)")
    ax.set_title("Training polygons: box p25-p75, whiskers p10-p90", fontsize=10)
    ax.legend(handles=handles, fontsize=8, loc="upper left")
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=200)
    return fig
