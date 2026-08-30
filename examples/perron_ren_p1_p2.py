"""Render Perron/Ren P1-P2 C++ benchmark outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

from lsdviztools3.charts.hillslope_profiles import profile_derivatives


LABELS = {
    "analytical_m": "Analytical Roering",
    "perron_qimp_m": "Perron Q-imp",
    "newton_m": "Newton",
    "cellular_c99_m": "Cellular C99",
    "cellular_c85_m": "Cellular C85",
}


def read_raster(path):
    with rasterio.open(path) as src:
        return src.read(1).astype(float), src.res[0]


def save_p1(root, out):
    frame = pd.read_csv(root / "p1_profile.csv")
    x = frame.x_m.to_numpy()
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for column, label in LABELS.items():
        elevation = frame[column].to_numpy()
        gradient, laplacian = profile_derivatives(x, elevation)
        style = {"lw": 2.0} if column == "analytical_m" else {"lw": 1.2}
        axes[0].plot(x, elevation, label=label, **style)
        axes[1].plot(x, gradient, label=label, **style)
        axes[2].plot(x, laplacian, label=label, **style)
    axes[0].set_ylabel("Elevation (m)")
    axes[1].set_ylabel("Gradient dz/dx")
    axes[2].set_ylabel("Laplacian (m$^{-1}$)")
    axes[2].set_xlabel("Distance from divide (m)")
    for ax in axes: ax.grid(alpha=0.25)
    axes[0].legend(ncol=3, fontsize=8)
    fig.suptitle("P1: Perron one-dimensional equilibrium")
    fig.tight_layout(); fig.savefig(out / "p1_profiles.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    analytical = frame.analytical_m.to_numpy()
    for column, label in list(LABELS.items())[1:]:
        residual = frame[column].to_numpy() - analytical
        axes[0].plot(x, residual, label=label)
        axes[1].plot(x[2:-2], residual[2:-2], label=label)
    axes[0].set_ylabel("Residual (m)"); axes[0].set_title("All profile cells")
    axes[1].set_ylabel("Residual (m)"); axes[1].set_title("Two-cell boundary buffer excluded")
    axes[1].set_xlabel("Distance from divide (m)")
    for ax in axes: ax.axhline(0, color="0.3", lw=0.7); ax.grid(alpha=0.25)
    axes[0].legend(ncol=2, fontsize=8)
    fig.tight_layout(); fig.savefig(out / "p1_residuals.png", dpi=180); plt.close(fig)

    save_summary(root / "p1_summary.csv", out / "p1_accuracy.png", "P1 accuracy")


def save_summary(csv_path, output, title):
    frame = pd.read_csv(csv_path)
    colors = ["#4D4D4D", "#0072B2", "#009E73", "#D55E00", "#CC79A7"][:len(frame)]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].bar(frame.model, frame.rmse_over_relief, color=colors)
    axes[0].set_yscale("log"); axes[0].set_ylabel("RMSE / relief")
    axes[1].bar(frame.model, frame.elapsed_s, color=colors)
    axes[1].set_yscale("log"); axes[1].set_ylabel("Wall time (s)")
    for ax in axes:
        ax.tick_params(axis="x", rotation=35); ax.grid(axis="y", alpha=0.25)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(output, dpi=180); plt.close(fig)


def save_p2(root, out):
    names = ("analytical", "perron_qimp", "ren_modular", "newton",
             "cellular_c99", "cellular_c85")
    arrays = {name: read_raster(root / f"p2_{name}.tif")[0] for name in names}
    analytical = arrays["analytical"]
    extent = (0, analytical.shape[1] * 2, 0, analytical.shape[0] * 2)
    common_names = tuple(name for name in names if name != "ren_modular")
    vmin = min(float(np.nanmin(arrays[name])) for name in common_names)
    vmax = max(float(np.nanmax(arrays[name])) for name in common_names)
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    for ax, name in zip(axes.ravel(), names):
        image = ax.imshow(arrays[name], extent=extent, origin="lower", cmap="terrain",
                          vmin=vmin, vmax=vmax)
        ax.set_title(name.replace("_", " ").title()); ax.set(xlabel="x (m)", ylabel="y (m)")
    fig.colorbar(image, ax=axes, label="Elevation (m)", shrink=0.8)
    fig.suptitle("P2: analytical and numerical axisymmetric equilibria")
    fig.savefig(out / "p2_elevation.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for ax, name in zip(axes.ravel(), names[1:]):
        residual = arrays[name] - analytical
        limit = float(np.nanmax(np.abs(residual))) or 1.0
        image = ax.imshow(residual, extent=extent, origin="lower", cmap="RdBu_r",
                          vmin=-limit, vmax=limit)
        ax.set_title(f"{name.replace('_', ' ').title()} - analytical\nmax |error|={limit:.4g} m")
        ax.set(xlabel="x (m)", ylabel="y (m)")
        fig.colorbar(image, ax=ax, label="Residual (m)", shrink=0.75)
    axes.ravel()[-1].axis("off")
    fig.suptitle("P2 elevation residuals; each panel has its own symmetric scale")
    fig.savefig(out / "p2_residuals.png", dpi=180); plt.close(fig)
    save_summary(root / "p2_summary.csv", out / "p2_accuracy.png", "P2 accuracy and runtime")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    save_p1(args.root, args.out); save_p2(args.root, args.out)


if __name__ == "__main__":
    main()
