"""Per-window HAND densities with the ladder's breaks drawn on them.

This is the figure that says whether a break is real.

A break value on its own tells you nothing: 2 m and 20 m look equally like
answers in a table. What distinguishes them is the shape of the distribution the
break was cut from. Two surfaces separated by a riser give two modes with a deep
dip between them, and the break sits in the dip. One broad surface gives one
mode, and a break through it is an arbitrary line that will move from window to
window as the noise moves.

The resolvability floor (`required_gap_m` = kappa * vertical accuracy) does NOT
answer this. It is an instrument floor: it says two levels cannot be claimed
closer together than the DEM can measure. On a lidar DEM two genuinely distinct
terraces may sit tens of centimetres apart and clear a floor set to that DEM's
accuracy, while a metre-scale break through a single mode on a coarse DEM clears
a much larger floor and means nothing. The floor is necessary and it is not
sufficient.

`separation_ratio` is the sufficient part, and it is annotated on every panel:
the smoothed density AT the break over the smaller flanking peak.

    <0.2   clean - the break sits in a real dip
     0.6   the antimode finder's default acceptance limit
    ~1.0   no dip at all; the break cuts through one population
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _breaks_for(windows: pd.DataFrame, sid, centre) -> list[float]:
    row = windows[(windows["source_id"] == sid) &
                  (np.isclose(windows["window_centre_m"], centre))]
    if row.empty:
        return []
    raw = row.iloc[0].get("breaks_m")
    if pd.isna(raw) or not str(raw).strip():
        return []
    return [float(v) for v in str(raw).split(";") if v.strip()]


def _seps_for(windows: pd.DataFrame, sid, centre) -> list[float]:
    row = windows[(windows["source_id"] == sid) &
                  (np.isclose(windows["window_centre_m"], centre))]
    if row.empty or "separation_ratio" not in row.columns:
        return []
    raw = row.iloc[0].get("separation_ratio")
    if pd.isna(raw) or not str(raw).strip():
        return []
    return [float(v) for v in str(raw).split(";") if v.strip()]


def plot_ladder_density(
    hist_csv: str | Path,
    windows_csv: str | Path,
    out_png: str | Path,
    source_id: int | None = None,
    max_panels: int = 12,
    title: str | None = None,
    dpi: int = 160,
) -> Path:
    """Plot each window's reference-field density with its breaks marked.

    Args:
        hist_csv: A `<stem>_ladder_histograms.csv`.
        windows_csv: The matching `<stem>_ladder_windows.csv`.
        out_png: Output PNG path.
        source_id: Restrict to one channel; None plots every channel found.
        max_panels: Cap on panels per channel, sampled evenly downstream.
        title: Optional figure title.
        dpi: Output resolution.

    Returns:
        The path written.
    """
    hist = pd.read_csv(hist_csv)
    wins = pd.read_csv(windows_csv)
    if hist.empty:
        raise ValueError(f"{hist_csv} has no histogram rows - was "
                         "write_ladder_histograms_csv set?")

    sources = ([source_id] if source_id is not None
               else sorted(hist["source_id"].unique()))

    panels = []
    for sid in sources:
        sub = hist[hist["source_id"] == sid]
        centres = sorted(sub["window_centre_m"].unique())
        if len(centres) > max_panels:
            # Sample evenly along the channel rather than taking the first N,
            # so the figure spans the reach instead of its top end.
            idx = np.linspace(0, len(centres) - 1, max_panels).astype(int)
            centres = [centres[i] for i in idx]
        for c in centres:
            panels.append((sid, c))

    n = len(panels)
    ncol = min(4, max(1, n))
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.7 * nrow),
                             squeeze=False)

    for k, (sid, centre) in enumerate(panels):
        ax = axes[k // ncol][k % ncol]
        sub = hist[(hist["source_id"] == sid) &
                   (np.isclose(hist["window_centre_m"], centre))]
        mid = 0.5 * (sub["bin_lo_m"].to_numpy() + sub["bin_hi_m"].to_numpy())
        cnt = sub["count"].to_numpy().astype(float)
        width = (sub["bin_hi_m"] - sub["bin_lo_m"]).mean()

        ax.bar(mid, cnt, width=width, color="#9ecae1", edgecolor="none")

        brks = _breaks_for(wins, sid, centre)
        seps = _seps_for(wins, sid, centre)
        for j, b in enumerate(brks):
            s = seps[j] if j < len(seps) else float("nan")
            # Colour the break by how well it actually separates, so a bad break
            # is visible without reading the annotation.
            colour = ("#238b45" if s < 0.2 else
                      "#fe9929" if s < 0.6 else "#cb181d")
            ax.axvline(b, color=colour, lw=1.8)
            ax.annotate(f"{b:.1f} m\nsep {s:.2f}", xy=(b, ax.get_ylim()[1]),
                        xytext=(3, -10), textcoords="offset points",
                        fontsize=7, color=colour, va="top")

        ax.set_title(f"ch {sid}  @ {centre/1000.0:.1f} km", fontsize=9)
        ax.set_xlabel("height above channel (m)", fontsize=8)
        ax.set_ylabel("count", fontsize=8)
        ax.tick_params(labelsize=7)

    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png
