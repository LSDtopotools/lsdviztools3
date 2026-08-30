"""Render the lsdtt3 hillslope comparison workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

from lsdviztools3.charts.hillslope_profiles import (
    load_profile_raster,
    plot_hillslope_profiles,
    plot_profile_residual_panels,
    profile_derivatives,
)
from lsdviztools3.charts.hillslope_residuals import (
    plot_residual_density,
    plot_residual_raster,
    residual_statistics,
    summarise_residuals,
)
from lsdviztools3.render.hillslope_comparison_animation import (
    discover_frame_pairs,
    render_hillslope_comparison_animation,
)


U = 4.0e-4
D = 0.1
SC = 0.8


def raster(path):
    with rasterio.open(path) as src:
        values = src.read(1, masked=True).astype(float).filled(np.nan)
        return values, src.profile


def render_comparison_hillshade(path):
    """Render a local-coordinate hillshade without requiring the optional GMT stack."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource

    elevation, profile = raster(path)
    dx = abs(float(profile["transform"].a))
    shade = LightSource(azdeg=315, altdeg=45).hillshade(
        elevation, vert_exag=1.0, dx=dx, dy=dx,
    )
    fig, ax = plt.subplots(figsize=(7, 6))
    extent = (0, elevation.shape[1] * dx, 0, elevation.shape[0] * dx)
    ax.imshow(shade, cmap="gray", origin="lower", extent=extent)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.tight_layout()
    return fig


def laplacian_2d(values, profile):
    transform = profile["transform"]
    dx, dy = abs(float(transform.a)), abs(float(transform.e))
    dz_dy, dz_dx = np.gradient(values, dy, dx)
    return np.gradient(dz_dx, dx, axis=1) + np.gradient(dz_dy, dy, axis=0)


def strip_metrics(values, profile):
    dx = abs(float(profile["transform"].a))
    elevation = values[values.shape[0] // 2]
    x = (np.arange(elevation.size) - (elevation.size - 1) / 2.0) * dx
    gradient = np.gradient(elevation, x)
    centre = elevation.size // 2
    window = slice(centre - 8, centre + 9)
    divide_laplacian = 2.0 * np.polyfit(x[window], elevation[window], 2)[0]
    return {
        "relief_m": float(elevation[centre] - 0.5 * (elevation[0] + elevation[-1])),
        "max_abs_gradient": float(np.max(np.abs(gradient[2:-2]))),
        "divide_laplacian_m-1": float(divide_laplacian),
    }


def analytic_strip(n=257, dx=5.0):
    x = (np.arange(n) - (n - 1) / 2.0) * dx
    L = abs(x[-1])
    z_linear = U * (L * L - x * x) / (2.0 * D)
    # Integrate the Roering steady relation by solving for S at each x.
    q = U * np.abs(x)
    slope = np.zeros_like(x)
    lo, hi = np.zeros_like(x), np.full_like(x, SC * (1.0 - 1e-10))
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = D * mid / (1.0 - (mid / SC) ** 2) - q
        lo = np.where(f < 0.0, mid, lo)
        hi = np.where(f < 0.0, hi, mid)
    slope = 0.5 * (lo + hi)
    z_roering = np.zeros_like(x)
    for i in range(n // 2 + 1, n):
        z_roering[i] = z_roering[i - 1] - 0.5 * (slope[i] + slope[i - 1]) * dx
    for i in range(n // 2 - 1, -1, -1):
        z_roering[i] = z_roering[i + 1] - 0.5 * (slope[i] + slope[i + 1]) * dx
    z_roering -= z_roering[0]
    curves = {}
    for label, elevation in (("linear", z_linear), ("roering", z_roering)):
        gradient, curvature = profile_derivatives(x, elevation)
        curves[label] = {
            "elevation": elevation,
            "gradient": gradient,
            "curvature": curvature,
        }
    return x, curves


def save_strip(root: Path):
    out = root / "plots" / "strip"
    out.mkdir(parents=True, exist_ok=True)
    analytic_x, analytic_z = analytic_strip()
    newton_x, newton_z = load_profile_raster(
        root / "strip_1d" / "newton_steady" / "strip_newton_hillslope_on_loaded_dem.tif"
    )
    cell_x, cell_z = load_profile_raster(
        root / "strip_1d" / "cellular_steady" / "strip_cellular_modelled_topography.tif"
    )
    newton_x = newton_x - newton_x[len(newton_x) // 2]
    cell_x = cell_x - cell_x[len(cell_x) // 2]
    if not (np.allclose(newton_x, analytic_x) and np.allclose(cell_x, analytic_x)):
        raise ValueError("Strip raster coordinates do not match the analytic grid")
    profiles = {"Newton": newton_z, "cellular": cell_z}
    fig, _ = plot_hillslope_profiles(
        analytic_x, profiles, analytic={"roering": analytic_z["roering"]},
        channel_indices=(0, len(analytic_x) - 1),
        title=("5 x 257 near-critical hillslope; cellular differences above "
               "r_match are constitutive"),
        include_residuals=False,
    )
    fig.savefig(out / "strip_profiles.png", dpi=160)
    import matplotlib.pyplot as plt
    plt.close(fig)

    residuals = {
        "Newton - Roering": newton_z - analytic_z["roering"]["elevation"],
        "cellular - Roering": cell_z - analytic_z["roering"]["elevation"],
        "cellular-Newton": cell_z - newton_z,
    }
    table = pd.DataFrame({"distance_m": analytic_x, **residuals})
    table.to_csv(out / "strip_residuals.csv", index=False)
    fig, _ = plot_residual_density(residuals, title="Strip residual density")
    fig.savefig(out / "strip_residual_density.png", dpi=160)
    plt.close(fig)
    hillslope = np.ones(len(analytic_x), dtype=bool)
    hillslope[[0, -1]] = False
    summary = summarise_residuals(residuals, hillslope_mask=hillslope)
    summary.to_csv(out / "strip_residual_summary.csv", index=False)

    fig, _ = plot_profile_residual_panels(
        analytic_x, residuals,
        title="Steady near-critical strip residuals (cellular continuation active)"
    )
    fig.savefig(out / "strip_residual_panels.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for label, elevation in (("Newton", newton_z), ("cellular", cell_z)):
        _, curvature = profile_derivatives(analytic_x, elevation)
        ax.plot(analytic_x, curvature, label=f"{label} curvature")
        ax.plot(analytic_x, -curvature, "--", label=f"{label} convexity")
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set(xlabel="Distance from divide (m)", ylabel="m$^{-1}$",
           title="Signed curvature and convexity")
    ax.grid(alpha=0.25); ax.legend(ncol=2, fontsize=8)
    fig.tight_layout(); fig.savefig(out / "strip_curvature_convexity.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    phase_profiles = {
        "Analytical Roering": analytic_z["roering"],
    }
    for label, elevation in (("Newton", newton_z), ("Cellular", cell_z)):
        gradient, laplacian = profile_derivatives(analytic_x, elevation)
        phase_profiles[label] = {"gradient": gradient, "curvature": laplacian}
    for label, quantities in phase_profiles.items():
        ax.scatter(quantities["gradient"][2:-2], quantities["curvature"][2:-2], s=13,
                   alpha=0.7, label=label)
    ax.set(xlabel="Gradient dz/dx", ylabel="Laplacian (m$^{-1}$)",
           title="Roering gradient-Laplacian phase space")
    ax.grid(alpha=0.25); ax.legend(); fig.tight_layout()
    fig.savefig(out / "strip_gradient_laplacian.png", dpi=160)
    plt.close(fig)

    diagnostics = []
    for label, elevation in (("Newton", newton_z), ("cellular", cell_z)):
        gradient, curvature = profile_derivatives(analytic_x, elevation)
        diagnostics.append({
            "model": label,
            "maximum_slope_ratio": float(np.max(np.abs(gradient)) / SC),
            "ridgetop_curvature_m-1": float(curvature[len(curvature) // 2]),
        })
    pd.DataFrame(diagnostics).to_csv(out / "strip_diagnostics.csv", index=False)


def save_2d(root: Path):
    out = root / "plots" / "model_2d"
    out.mkdir(parents=True, exist_ok=True)
    newton = root / "model_2d" / "newton_steady" / "two_d_newton_hillslope_on_loaded_dem.tif"
    cellular = root / "model_2d" / "cellular_steady" / "two_d_cellular_modelled_topography.tif"
    a, profile = raster(newton)
    b, _ = raster(cellular)
    residual = b - a
    hillslope = np.isfinite(residual)
    mask_path = root / "inputs" / "imposed_channel_101x101" / "two_d_channel_mask.asc"
    with rasterio.open(mask_path) as src:
        channel = src.read(1) > 0.5
    hillslope &= ~channel
    rows = []
    for scope, values in (("full", residual), ("hillslope", residual[hillslope]),
                          ("channel", residual[channel])):
        rows.append({"scope": scope, "comparison": "cellular-Newton",
                     **residual_statistics(values)})
    pd.DataFrame(rows).to_csv(out / "two_d_residual_summary.csv", index=False)
    with rasterio.open(out / "two_d_residual.tif", "w", **profile) as dst:
        dst.write(residual.astype("float32"), 1)
    fig, _ = plot_residual_density(
        {"full": residual, "hillslope": residual[hillslope]},
        title="101 x 101 cellular - Newton residual density",
    )
    fig.savefig(out / "two_d_residual_density.png", dpi=160)
    import matplotlib.pyplot as plt
    plt.close(fig)
    extent = (0, a.shape[1] * abs(profile["transform"].a),
              0, a.shape[0] * abs(profile["transform"].e))
    limit = float(np.nanmax(np.abs(residual)))
    fig, ax = plt.subplots(figsize=(7, 6))
    image = ax.imshow(residual, cmap="RdBu_r", vmin=-limit, vmax=limit,
                      extent=extent, origin="lower")
    ax.set(xlabel="x (m)", ylabel="y (m)", title="Elevation residual: cellular - Newton")
    fig.colorbar(image, ax=ax, label="residual (m)")
    fig.tight_layout()
    fig.savefig(out / "two_d_residual_map.png", dpi=160)
    plt.close(fig)
    lap_residual = laplacian_2d(b, profile) - laplacian_2d(a, profile)
    lap_limit = float(np.nanmax(np.abs(lap_residual)))
    fig, ax = plt.subplots(figsize=(7, 6))
    image = ax.imshow(lap_residual, cmap="PuOr_r", vmin=-lap_limit, vmax=lap_limit,
                      extent=extent, origin="lower")
    ax.set(xlabel="x (m)", ylabel="y (m)",
           title="Laplacian residual: cellular - Newton")
    fig.colorbar(image, ax=ax, label="Laplacian residual (m$^{-1}$)")
    fig.tight_layout(); fig.savefig(out / "two_d_laplacian_residual_map.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 6))
    dx = abs(float(profile["transform"].a))
    for values, label in ((a, "Newton"), (b, "Cellular")):
        gy, gx = np.gradient(values, dx, dx)
        gradient = np.hypot(gx, gy)
        laplacian = laplacian_2d(values, profile)
        interior = np.s_[2:-2, 2:-2]
        ax.scatter(gradient[interior].ravel(), laplacian[interior].ravel(),
                   s=4, alpha=0.2, label=label)
    ax.set(xlabel="Gradient magnitude", ylabel="Laplacian (m$^{-1}$)",
           title="2-D gradient-Laplacian phase space")
    ax.grid(alpha=0.25); ax.legend(); fig.tight_layout()
    fig.savefig(out / "two_d_gradient_laplacian.png", dpi=160)
    plt.close(fig)
    finite = np.concatenate((a[np.isfinite(a)], b[np.isfinite(b)]))
    vmin, vmax = float(finite.min()), float(finite.max())
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    image = None
    for ax, values, label in ((axes[0], a, "Newton"), (axes[1], b, "Cellular")):
        image = ax.imshow(values, cmap="terrain", vmin=vmin, vmax=vmax,
                          extent=extent, origin="lower")
        ax.set_title(f"{label} steady elevation")
        ax.set(xlabel="x (m)", ylabel="y (m)")
    fig.colorbar(image, ax=axes, label="Elevation (m)", shrink=0.85)
    fig.savefig(out / "two_d_elevation.png", dpi=160)
    plt.close(fig)
    for label, path in (("newton", newton), ("cellular", cellular)):
        fig = render_comparison_hillshade(path)
        fig.savefig(out / f"two_d_{label}_hillshade.png", dpi=160)
        plt.close(fig)
    (out / "two_d_manifest.txt").write_text(
        "Newton: " + str(newton) + "\nCellular: " + str(cellular)
        + "\nResidual: cellular - Newton\nChannel cells: 40\n"
    )


def save_scan(root: Path):
    from lsdviztools3.charts.estar_rstar import plot_estar_rstar_model_scan
    import matplotlib.pyplot as plt

    summary = root / "strip_scan" / "strip_scan_summary.csv"
    if not summary.is_file():
        raise FileNotFoundError(f"Run run_strip_scan.py before plotting: {summary}")
    out = root / "plots" / "strip_scan"
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(summary)
    fig, _ = plot_estar_rstar_model_scan(frame)
    fig.savefig(out / "strip_scan_estar_rstar.png", dpi=180)
    plt.close(fig)
    frame["R_star_theory"] = frame["E_star"].map(
        lambda e: (np.sqrt(1 + e * e) - np.log((1 + np.sqrt(1 + e * e)) / 2) - 1) / e
    )
    frame["R_star_residual"] = frame["R_star"] - frame["R_star_theory"]
    frame.to_csv(out / "strip_scan_estar_rstar.csv", index=False)


def save_transient(root: Path):
    out = root / "plots" / "animations"
    out.mkdir(parents=True, exist_ok=True)
    for stem, newton, cellular, mask_path in (
        ("strip", root / "strip_1d" / "newton_transient",
         root / "strip_1d" / "cellular_transient",
         root / "inputs" / "strip_5x257" / "strip_hillslope_strip_channel_mask.tif"),
        ("two_d", root / "model_2d" / "newton_transient",
         root / "model_2d" / "cellular_transient",
         root / "inputs" / "imposed_channel_101x101" / "two_d_channel_mask.asc"),
    ):
        pairs = discover_frame_pairs(newton, cellular, pattern="*_model_iteration_*.tif")
        nstats = next(newton.glob("*_model_elevation_stats.csv"))
        na = pd.read_csv(nstats)
        times = dict(zip(na["iteration"].astype(int), na["time_yr"].astype(float)))
        channel, _ = raster(mask_path)
        hillslope = np.isfinite(channel) & (channel <= 0.5)
        rows = []
        metric_rows = []
        for n, c, iteration in pairs:
            an, frame_profile = raster(n); ac, _ = raster(c)
            valid = np.isfinite(ac) & np.isfinite(an)
            difference = ac - an
            full = difference[valid]
            hill = difference[valid & hillslope]
            stats = residual_statistics(full)
            rows.append({"iteration": int(iteration), "time_yr": times[int(iteration)],
                         **{key: stats[key] for key in ("mean", "rms", "linf", "p05", "median", "p95")},
                         "hillslope_mean": float(hill.mean()),
                         "hillslope_rms": float(np.sqrt(np.mean(hill * hill)))})
            if stem == "strip":
                for model, values in (("Newton", an), ("Cellular", ac)):
                    metric_rows.append({"iteration": int(iteration),
                                        "time_yr": times[int(iteration)],
                                        "model": model,
                                        **strip_metrics(values, frame_profile)})
        stats = pd.DataFrame(rows)
        stats_path = root / ("strip_1d" if stem == "strip" else "model_2d") \
            / f"{stem}_transient_residuals.csv"
        stats.to_csv(stats_path, index=False)
        from lsdviztools3.charts.hillslope_transient_residuals import plot_transient_residuals
        fig, _ = plot_transient_residuals(stats, title=f"{stem} transient residuals")
        fig.savefig(out / f"{stem}_transient_residuals.png", dpi=160)
        import matplotlib.pyplot as plt
        plt.close(fig)
        if stem == "strip":
            from lsdviztools3.charts.hillslope_transient_residuals import (
                plot_strip_transient_metrics,
            )
            metrics = pd.DataFrame(metric_rows)
            metrics_path = root / "strip_1d" / "strip_transient_metrics.csv"
            metrics.to_csv(metrics_path, index=False)
            fig, _ = plot_strip_transient_metrics(
                metrics, title="Strip response to doubled uplift at 5,000 yr"
            )
            fig.savefig(out / "strip_transient_metrics.png", dpi=160)
            plt.close(fig)
        gif = render_hillslope_comparison_animation(
            newton, cellular, out / f"{stem}_newton_cellular.gif",
            pattern="*_model_iteration_*.tif", stats_csv=stats_path,
            fps=4, keep_frames=False, layout="strip" if stem == "strip" else "2d",
        )
        print(f"Wrote {gif}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--static-only", action="store_true")
    parser.add_argument("--transient-only", action="store_true")
    args = parser.parse_args()
    if args.static_only and args.transient_only:
        parser.error("--static-only and --transient-only are mutually exclusive")
    if not args.transient_only:
        save_strip(args.root)
        save_2d(args.root)
        save_scan(args.root)
    if not args.static_only:
        save_transient(args.root)


if __name__ == "__main__":
    main()
