"""The reference-field ladder's level selection, per channel and along it.

`lsdtt-fluvial-surfaces` bands a height reference (normally height above the
nearest drainage) and labels connected components within each band. It does that
per along-channel WINDOW, so the break between one level and the next is free to
change downstream as incision changes - and whether it succeeded at all changes
downstream too.

`<stem>_ladder_windows.csv` is the record of that, one row per window. This
module turns it into the figure the CSV is standing in for: where the breaks sat
along each channel, how far apart the levels ended up, and which windows the DEM
could not actually resolve.

WHY THE STATUS PANEL IS NOT DECORATION. A window reporting `unresolved` did not
find levels the DEM can separate, and the ladder emits ONE band there rather
than inventing a number. A reader who sees only the break trace will read a gap
in it as missing data instead of as a statement. The bottom panel is where that
statement is made, and it is why the two panels share an x axis.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Status -> (colour, human label). Resolved and forced are both "the ladder
# produced the level count asked of it"; they differ only in who chose it.
_STATUS_STYLE = {
    "resolved": ("#2c7fb8", "resolved"),
    "forced": ("#41b6c4", "forced (explicit level count)"),
    "unresolved": ("#d95f0e", "UNRESOLVED - one band emitted"),
    "degenerate": ("#999999", "degenerate - too few pixels"),
}


def _split_breaks(value) -> list[float]:
    """Parse a semicolon-separated break list; empty and NaN give []."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    out = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        try:
            out.append(float(part))
        except ValueError:
            continue
    return out


def plot_ladder_windows(
    windows_csv: str | Path,
    out_png: str | Path,
    title: str | None = None,
    dpi: int = 160,
) -> Path:
    """Plot the ladder's break selection along each channel.

    Args:
        windows_csv: A `<stem>_ladder_windows.csv` from lsdtt-fluvial-surfaces.
        out_png: Output PNG path.
        title: Optional figure title.
        dpi: Output resolution.

    Returns:
        The path written.
    """
    df = pd.read_csv(windows_csv)
    if df.empty:
        raise ValueError(f"{windows_csv} has no windows to plot")

    sources = sorted(df["source_id"].unique())
    n_src = len(sources)

    fig, axes = plt.subplots(
        2, n_src,
        figsize=(1.0 + 5.0 * n_src, 7.0),
        sharex="col",
        gridspec_kw={"height_ratios": [3, 1]},
        squeeze=False,
    )

    for col, sid in enumerate(sources):
        sub = df[df["source_id"] == sid].sort_values("window_centre_m")
        x = sub["window_centre_m"].to_numpy() / 1000.0     # km along channel
        ax = axes[0][col]

        # Every break in every window. A k-level ladder has k-1 of them, so a
        # window can contribute more than one point to this trace.
        bx, by = [], []
        for xi, raw in zip(x, sub["breaks_m"]):
            for b in _split_breaks(raw):
                bx.append(xi)
                by.append(b)
        if bx:
            ax.plot(bx, by, "o-", color="#08519c", ms=4, lw=1.2,
                    label="break between levels", zorder=3)

        # Band centres: where each level's surface actually SITS, which is what
        # the resolvability floor is measured between - not the breaks.
        cx, cy = [], []
        for xi, raw in zip(x, sub["band_centres_m"]):
            for c in _split_breaks(raw):
                cx.append(xi)
                cy.append(c)
        if cx:
            ax.scatter(cx, cy, s=18, facecolor="none", edgecolor="#756bb1",
                       lw=1.0, label="band centre (median of the level)",
                       zorder=2)

        # Shade windows the DEM could not resolve, so a gap in the break trace
        # reads as a statement rather than as missing data.
        for xi, st in zip(x, sub["status"]):
            if st in ("unresolved", "degenerate"):
                ax.axvspan(xi - 0.05, xi + 0.05,
                           color=_STATUS_STYLE[st][0], alpha=0.18, zorder=0)

        ax.set_ylabel("height above channel (m)")
        ax.set_title(f"channel {sid}")
        ax.grid(alpha=0.3)
        if col == 0:
            ax.legend(fontsize=8, loc="upper left")

        # Status panel.
        axs = axes[1][col]
        for st, (colour, label) in _STATUS_STYLE.items():
            m = sub["status"] == st
            if m.any():
                axs.scatter(x[m.to_numpy()], sub.loc[m, "k_selected"],
                            color=colour, s=28, label=label, zorder=3)
        # The gap the ladder had to clear. Constant per run, but drawn per panel
        # so neither panel has to be read against the other.
        gap = sub["required_gap_m"].max()
        axs.plot(x, sub["min_band_centre_gap_m"], "-", color="#525252", lw=1.0,
                 alpha=0.7, label="min band-centre gap (m)")
        if gap > 0:
            axs.axhline(gap, ls="--", color="#cb181d", lw=1.0,
                        label=f"required gap {gap:.2f} m")
        axs.set_xlabel("distance along channel (km)")
        axs.set_ylabel("levels / gap (m)")
        axs.grid(alpha=0.3)
        if col == 0:
            axs.legend(fontsize=7, loc="upper left", ncol=1)

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png


def summarise_ladder(windows_csv: str | Path) -> pd.DataFrame:
    """One row per channel: how the ladder fared along it.

    Args:
        windows_csv: A `<stem>_ladder_windows.csv`.

    Returns:
        A per-source summary frame.
    """
    df = pd.read_csv(windows_csv)
    rows = []
    for sid, sub in df.groupby("source_id"):
        breaks = [b for raw in sub["breaks_m"] for b in _split_breaks(raw)]
        rows.append({
            "source_id": sid,
            "n_windows": len(sub),
            "n_resolved": int((sub["status"].isin(["resolved", "forced"])).sum()),
            "n_unresolved": int((sub["status"] == "unresolved").sum()),
            "n_degenerate": int((sub["status"] == "degenerate").sum()),
            "break_min_m": min(breaks) if breaks else float("nan"),
            "break_median_m": float(pd.Series(breaks).median()) if breaks else float("nan"),
            "break_max_m": max(breaks) if breaks else float("nan"),
            "required_gap_m": float(sub["required_gap_m"].max()),
        })
    return pd.DataFrame(rows)
