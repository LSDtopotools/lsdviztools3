"""QQ plots of the HAND floodplain thresholds, first and second pass.

lsdtt3 ``lsdtt-fluvial-surfaces`` (``floodplain_threshold_algorithm_choice
qq_gaussian``) writes the QQ table each HAND bound was read from and a summary,
``<stem>_floodplain_thresholds.csv`` (pass, quantity, algorithm, threshold,
qq_index, n_values, table_fname). This draws one panel per HAND pass:

  pass 1, qq_gaussian (Clubb et al. 2017): the LOWEST quantile from which the
      data agree with the Gaussian line - the floodplain bound;
  pass 2, qq_upper_departure: over the pass-1 floodplain cells, the HIGHEST
      quantile at which they still agree - above it HAND departs upward
      (fans, upland flats).

Visualisation only: the thresholds and their quantiles are read from lsdtt3's
CSV, never recomputed here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

_PASS_TEXT = {
    "qq_gaussian": "lowest quantile agreeing with the Gaussian line",
    "qq_upper_departure": "highest quantile agreeing with the Gaussian line",
}


def read_floodplain_thresholds(path) -> pd.DataFrame:
    """Read ``<stem>_floodplain_thresholds.csv``; table paths made absolute."""
    path = Path(path)
    t = pd.read_csv(path, dtype={"table_fname": str})
    t["table_path"] = [path.parent / f if isinstance(f, str) and f else None
                       for f in t["table_fname"]]
    return t


def plot_floodplain_hand_qq(thresholds_csv, out_png: Optional[str] = None,
                            *, figsize=(6.0, 4.2)):
    """One QQ panel per HAND pass in ``thresholds_csv``.

    Each panel: the data quantiles against the standard normal variate, the
    Gaussian reference line, and the threshold - a horizontal line at its value
    and a marker at the quantile it was read at.

    Args:
        thresholds_csv: lsdtt3 ``<stem>_floodplain_thresholds.csv``.
        out_png: Optional output path; the figure is also returned.
        figsize: Size of ONE panel (inches); panels sit side by side.

    Returns:
        The matplotlib ``Figure``.
    """
    t = read_floodplain_thresholds(thresholds_csv)
    hand = t[(t["quantity"] == "hand") & t["table_path"].notna()]
    hand = hand[hand["algorithm"].isin(_PASS_TEXT)].sort_values("pass")
    if hand.empty:
        raise ValueError(f"No QQ-chosen HAND threshold in {thresholds_csv}.")

    n = len(hand)
    fig, axes = plt.subplots(1, n, figsize=(figsize[0] * n, figsize[1]),
                             squeeze=False)
    for ax, (_, row) in zip(axes[0], hand.iterrows()):
        qq = pd.read_csv(row["table_path"])
        z, dq, ref = (qq["standard_normal_variate"], qq["data_quantile"],
                      qq["normal_reference"])
        i = int(row["qq_index"])
        thr = float(row["threshold"])
        ax.plot(z, dq, ".", ms=1.5, color="0.25", label="data quantiles")
        ax.plot(z, ref, "-", lw=1.0, color="tab:red", label="Gaussian line (25-75 %)")
        ax.axhline(thr, ls="--", lw=1.0, color="tab:blue")
        ax.plot([z.iloc[i]], [dq.iloc[i]], "o", ms=6, mfc="none", mec="tab:blue",
                label=f"threshold {thr:.3g} m")
        lo, hi = float(dq.min()), float(dq.max())
        pad = 0.05 * (hi - lo)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlabel("standard normal variate")
        ax.set_ylabel("HAND (m)")
        ax.set_title(f"pass {int(row['pass'])}: {row['algorithm']}\n"
                     f"{_PASS_TEXT[row['algorithm']]}; "
                     f"n = {int(row['n_values']):,}", fontsize=9)
        ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=200)
    return fig
