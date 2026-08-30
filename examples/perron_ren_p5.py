"""Render P5 stability and robustness envelope diagnostics."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


COLORS = {
    "perron_qimp": "#0072B2", "ren_modular": "#D55E00",
    "newton_fixed": "#009E73", "newton_adaptive": "#56B4E9",
    "cellular_c99": "#CC79A7", "cellular_c85": "#7A3E9D",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    frame = pd.read_csv(args.csv)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex="col")
    for row, test in enumerate(("P1", "P2")):
        subset = frame[frame.test == test]
        for (model, initial), group in subset.groupby(["model", "initial"], sort=False):
            color = COLORS[model]
            style = "-" if initial == "flat" else "--"
            label = f"{model.replace('_', ' ')} ({initial})"
            axes[row, 0].plot(group.requested_dt_yr, group.max_slope_ratio,
                              style, color=color, label=label)
            continued = group.continuation_completed.astype(bool)
            axes[row, 0].scatter(group.loc[continued, "requested_dt_yr"],
                                 group.loc[continued, "max_slope_ratio"], color=color, s=28)
            axes[row, 0].scatter(group.loc[~continued, "requested_dt_yr"],
                                 group.loc[~continued, "max_slope_ratio"], color=color,
                                 marker="x", s=60)
            axes[row, 1].plot(group.requested_dt_yr, group.elapsed_s,
                              style, color=color, label=label)
        axes[row, 0].axhline(1.0, color="firebrick", ls=":", label="Sc threshold")
        axes[row, 0].set_title(f"{test}: slope ratio and continuation status")
        axes[row, 0].set_ylabel("Maximum slope / Sc")
        axes[row, 1].set_title(f"{test}: one requested interval wall time")
        axes[row, 1].set_ylabel("Wall time (s)")
        axes[row, 0].legend(fontsize=6, ncol=2)
        axes[row, 1].legend(fontsize=6, ncol=2)
    for ax in axes.ravel():
        ax.set_xscale("log"); ax.grid(alpha=0.25)
    axes[0, 1].set_yscale("log"); axes[1, 1].set_yscale("log")
    axes[1, 0].set_xlabel("Requested interval (yr)")
    axes[1, 1].set_xlabel("Requested interval (yr)")
    fig.suptitle("P5 stability envelope: circles continue; x marks next-step rejection")
    fig.tight_layout(); args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
