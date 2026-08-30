"""Render P6 spatial convergence, runtime, and P2 anisotropy diagnostics."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def observed_order(frame):
    rows = []
    for (test, model), group in frame.groupby(["test", "model"]):
        group = group[(group.completed == 1) & (group.regularization_active == 0)]
        if len(group) < 3:
            continue
        slope = np.polyfit(np.log(group.dx_m), np.log(group.rmse_m), 1)[0]
        rows.append({"test": test, "model": model, "n_points": len(group),
                     "observed_spatial_order": slope})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    frame = pd.read_csv(args.csv)
    orders = observed_order(frame); orders.to_csv(args.summary, index=False)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for row, test in enumerate(("P1", "P2")):
        subset = frame[frame.test == test] if test == "P1" else frame[frame.test.isin(("P2", "P2-ren-archive"))]
        for model, group in subset.groupby("model", sort=False):
            regularized = bool(group.regularization_active.any())
            label = model.replace("_", " ") + (" (regularized)" if regularized else "")
            axes[row, 0].plot(group.dx_m, group.rmse_over_relief, "o-", label=label)
            axes[row, 1].plot(group.dx_m, group.elapsed_s, "o-", label=label)
        axes[row, 0].set(xscale="log", yscale="log", ylabel="RMSE / relief",
                          title=f"{test} spatial convergence")
        axes[row, 1].set(xscale="log", yscale="log", ylabel="Wall time (s)",
                          title=f"{test} 3 Myr runtime")
        axes[row, 0].legend(fontsize=7); axes[row, 1].legend(fontsize=7)
        axes[row, 0].grid(alpha=0.25); axes[row, 1].grid(alpha=0.25)
    p2 = frame[frame.test.isin(("P2", "P2-ren-archive"))]
    for model, group in p2.groupby("model", sort=False):
        axes[0, 2].plot(group.dx_m, group.radial_azimuthal_rmse_m, "o-",
                         label=model.replace("_", " "))
        axes[1, 2].plot(group.dx_m, group.maximum_absolute_error_m, "o-",
                         label=model.replace("_", " "))
    axes[0, 2].set(xscale="log", yscale="log", ylabel="Radial-bin azimuthal RMSE (m)",
                    title="P2 axisymmetric anisotropy")
    axes[1, 2].set(xscale="log", yscale="log", ylabel="Maximum absolute error (m)",
                    title="P2 maximum error")
    for ax in axes[:, 2]: ax.legend(fontsize=7); ax.grid(alpha=0.25)
    for ax in axes[-1]: ax.set_xlabel("Grid spacing dx (m)")
    fig.suptitle("P6: fixed-domain spatial refinement")
    fig.tight_layout(); args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
