"""Plot P3 equilibrium timestep-independence results."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--p1-summary", required=True, type=Path)
    parser.add_argument("--p2-summary", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    frame = pd.read_csv(args.csv)
    p1 = pd.read_csv(args.p1_summary).set_index("model")
    p2 = pd.read_csv(args.p2_summary).set_index("model")

    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    groups = (
        ("P1", "perron_qimp", "P1 Perron Q-imp", "#0072B2"),
        ("P2", "perron_qimp", "P2 Perron Q-imp", "#009E73"),
        ("P2-ren-archive", "ren_modular", "P2 Ren archive convention", "#D55E00"),
    )
    for test, model, label, color in groups:
        subset = frame[(frame.test == test) & (frame.model == model)].sort_values("dt_yr")
        row = 0 if test == "P1" else 1
        completed = subset.completed.astype(bool)
        axes[row, 0].plot(subset.loc[completed, "dt_yr"],
                          subset.loc[completed, "rmse_over_relief"], "o-",
                          color=color, label=label)
        axes[row, 0].scatter(subset.loc[~completed, "dt_yr"],
                             subset.loc[~completed, "rmse_over_relief"], marker="x",
                             s=80, color=color, label=f"{label} rejected")
        axes[row, 1].plot(subset.dt_yr, subset.elapsed_s, "o-", color=color,
                          label=label)

    for model, style in (("newton", "--"), ("cellular_c99", ":"),
                         ("cellular_c85", "-.")):
        axes[0, 0].axhline(p1.loc[model, "rmse_over_relief"], ls=style,
                           color="0.35", label=f"P1 {model} reference")
    for model, style in (("newton", "--"), ("cellular_c99", ":")):
        axes[1, 0].axhline(p2.loc[model, "rmse_over_relief"], ls=style,
                           color="0.35", label=f"P2 {model} reference")

    for row, summary, title in ((0, p1, "P1 all-model, 3 Myr wall time"),
                                (1, p2, "P2 all-model, 3 Myr wall time")):
        ordered = summary.sort_values("elapsed_s")
        axes[row, 2].bar(ordered.index, ordered.elapsed_s)
        axes[row, 2].set_yscale("log")
        axes[row, 2].set_ylabel("Wall time (s)")
        axes[row, 2].set_title(title)
        axes[row, 2].tick_params(axis="x", rotation=35)
        axes[row, 2].grid(axis="y", alpha=0.25)
    axes[0, 0].set_title("P1 equilibrium accuracy")
    axes[1, 0].set_title("P2 equilibrium accuracy")
    axes[0, 1].set_title("P1 Q-imp timestep runtime")
    axes[1, 1].set_title("P2 Q-imp/Ren timestep runtime")
    for ax in axes[:, 0]:
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylabel("RMSE / relief")
        ax.grid(alpha=0.25); ax.legend(fontsize=7)
    for ax in axes[:, 1]:
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylabel("Wall time (s)")
        ax.grid(alpha=0.25); ax.legend(fontsize=7)
    axes[1, 0].set_xlabel("Timestep (yr)"); axes[1, 1].set_xlabel("Timestep (yr)")
    fig.suptitle("P3: equilibrium timestep independence and stability")
    fig.tight_layout(); args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
