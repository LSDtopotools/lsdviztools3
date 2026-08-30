"""Plot P4 transient timestep convergence and matched-duration runtimes."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    frame = pd.read_csv(args.csv)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for test, ax in (("P1", axes[0, 0]), ("P2", axes[1, 0])):
        subset = frame[(frame.test == test) & (frame.model == "perron_qimp")]
        ax.plot(subset.dt_yr, subset.rms_over_relief, "o-", label="Perron Q-imp")
        ax.set(xscale="log", yscale="log", xlabel="Timestep (yr)",
               ylabel="RMS transient error / relief", title=f"{test} transient convergence")
        ax.grid(alpha=0.25); ax.legend()
    ren = frame[frame.test == "P2-ren-archive"]
    axes[1, 0].plot(ren.dt_yr, ren.rms_over_relief, "s--", color="#D55E00",
                    label="Ren modular (archive convention)")
    axes[1, 0].legend(fontsize=8)

    for test, ax in (("P1", axes[0, 1]), ("P2", axes[1, 1])):
        subset = frame[frame.test == test].copy()
        subset.loc[subset.dt_yr == 0, "dt_yr"] = np.nan
        for model, group in subset.groupby("model", sort=False):
            x = group.dt_yr.to_numpy(float)
            label = model.replace("_", " ")
            if np.all(~np.isfinite(x)):
                ax.scatter([1], group.elapsed_s, marker="D", label=f"{label} (adaptive/CFL)")
            else:
                ax.plot(x, group.elapsed_s, "o-", label=label)
        ax.set(xscale="log", yscale="log", xlabel="Requested timestep (yr)",
               ylabel="Wall time (s)", title=f"{test}: 250 kyr matched-duration wall time")
        ax.grid(alpha=0.25); ax.legend(fontsize=7)
    fig.suptitle("P4 transient timestep error and wall time")
    fig.tight_layout(); args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
