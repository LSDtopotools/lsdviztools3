"""Matplotlib/Pillow animation for paired Newton and cellular rasters."""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import numpy as np

_TIME_RE = re.compile(r"(?:iteration|iter|time)[^0-9]*([0-9]+(?:\.[0-9]+)?)", re.I)


def _raster(path):
    try:
        import rasterio
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Hillslope animation raster paths need rasterio.") from exc
    with rasterio.open(path) as src:
        return np.asarray(src.read(1, masked=True).filled(np.nan), float), src.transform


def _validate_pair(a, transform_a, b, transform_b):
    if a.shape != b.shape:
        raise ValueError("Newton and cellular frame shapes differ")
    if not transform_a.almost_equals(transform_b):
        raise ValueError("Newton and cellular frame transforms differ")


def _cell_sizes(transform):
    if transform.b != 0 or transform.d != 0:
        raise ValueError("Rotated rasters are not supported")
    return abs(float(transform.a)), abs(float(transform.e))


def _laplacian(values, transform):
    dx, dy = _cell_sizes(transform)
    dz_dy, dz_dx = np.gradient(values, dy, dx)
    return np.gradient(dz_dx, dx, axis=1) + np.gradient(dz_dy, dy, axis=0)


def _strip_quantities(values, transform):
    dx, _ = _cell_sizes(transform)
    elevation = np.asarray(values[values.shape[0] // 2], float)
    distance = (np.arange(elevation.size) - (elevation.size - 1) / 2.0) * dx
    gradient = np.gradient(elevation, distance)
    laplacian = np.gradient(gradient, distance)
    return distance, elevation, gradient, laplacian


def _key(path: Path):
    match = _TIME_RE.search(path.stem)
    return (float(match.group(1)) if match else None)


def discover_frame_pairs(newton_dir, cellular_dir, *, pattern="*.tif"):
    """Pair frames by an iteration/time token or exact filename stem.

    A missing partner is an error rather than silently shortening the movie.
    ``.tif`` is the default because it is the supported analysis output in the
    approved plan; callers may explicitly pass another rasterio-readable glob.
    """
    def collect(directory):
        paths = sorted(Path(directory).glob(pattern))
        if not paths:
            raise FileNotFoundError(f"No frames matching {pattern!r} in {directory}")
        keyed = {}
        for path in paths:
            key = _key(path)
            key = ("name", path.stem) if key is None else ("number", key)
            if key in keyed:
                raise ValueError(f"Duplicate frame key {key!r} in {directory}")
            keyed[key] = path
        return keyed
    left, right = collect(newton_dir), collect(cellular_dir)
    if set(left) != set(right):
        raise ValueError("Newton and cellular frame sets do not match: "
                         f"Newton-only={sorted(set(left)-set(right))}, "
                         f"cellular-only={sorted(set(right)-set(left))}")
    keys = sorted(left, key=lambda k: (k[0], k[1]))
    return [(left[k], right[k], k[1]) for k in keys]


def _limits(pairs):
    zmin, zmax, rmax, lap_rmax = np.inf, -np.inf, 0.0, 0.0
    for newton, cellular, _ in pairs:
        a, transform_a = _raster(newton)
        b, transform_b = _raster(cellular)
        _validate_pair(a, transform_a, b, transform_b)
        values = np.concatenate((a[np.isfinite(a)], b[np.isfinite(b)]))
        if values.size:
            zmin, zmax = min(zmin, values.min()), max(zmax, values.max())
        residual = b - a
        finite = residual[np.isfinite(residual)]
        if finite.size:
            rmax = max(rmax, float(np.max(np.abs(finite))))
        lap_residual = _laplacian(b, transform_b) - _laplacian(a, transform_a)
        finite_lap = lap_residual[np.isfinite(lap_residual)]
        if finite_lap.size:
            lap_rmax = max(lap_rmax, float(np.max(np.abs(finite_lap))))
    if not np.isfinite(zmin):
        raise ValueError("No finite values in paired frames")
    if zmin == zmax:
        zmax = zmin + 1.0
    if rmax == 0:
        rmax = 1.0
    if lap_rmax == 0:
        lap_rmax = 1.0
    return float(zmin), float(zmax), float(rmax), float(lap_rmax)


def _strip_limits(pairs, reference_pair):
    limits = [[np.inf, -np.inf] for _ in range(3)]
    anomaly_limits = [[np.inf, -np.inf] for _ in range(3)]
    difference_limits = [[np.inf, -np.inf] for _ in range(3)]
    references = []
    for path in reference_pair[:2]:
        values, transform = _raster(path)
        references.append(_strip_quantities(values, transform)[1:])
    for newton, cellular, _ in pairs:
        pair_quantities = []
        for model_index, path in enumerate((newton, cellular)):
            values, transform = _raster(path)
            quantities = _strip_quantities(values, transform)[1:]
            pair_quantities.append(quantities)
            for index, (bound, quantity) in enumerate(zip(limits, quantities)):
                finite = quantity[np.isfinite(quantity)]
                bound[0] = min(bound[0], float(finite.min()))
                bound[1] = max(bound[1], float(finite.max()))
                anomaly = quantity - references[model_index][index]
                finite_anomaly = anomaly[np.isfinite(anomaly)]
                anomaly_limits[index][0] = min(
                    anomaly_limits[index][0], float(finite_anomaly.min()))
                anomaly_limits[index][1] = max(
                    anomaly_limits[index][1], float(finite_anomaly.max()))
        for index in range(3):
            response_difference = ((pair_quantities[1][index] - references[1][index])
                                   - (pair_quantities[0][index] - references[0][index]))
            finite = response_difference[np.isfinite(response_difference)]
            difference_limits[index][0] = min(difference_limits[index][0], float(finite.min()))
            difference_limits[index][1] = max(difference_limits[index][1], float(finite.max()))
    return limits, anomaly_limits, difference_limits, references


def render_strip_comparison_frame(newton, cellular, *, limits, anomaly_limits,
                                  difference_limits, references, time,
                                  figsize=(18, 8)):
    """Render absolute profiles and changes from the pre-step reference."""
    import matplotlib.pyplot as plt
    a, transform_a = _raster(newton)
    b, transform_b = _raster(cellular)
    _validate_pair(a, transform_a, b, transform_b)
    qa = _strip_quantities(a, transform_a)
    qb = _strip_quantities(b, transform_b)
    fig, axes = plt.subplots(3, 3, figsize=figsize, sharex="col")
    labels = ("Elevation (m)", "Gradient dz/dx", "Laplacian (m$^{-1}$)")
    for index, ylabel in enumerate(labels, start=1):
        ax = axes[index - 1, 0]
        anomaly_ax = axes[index - 1, 1]
        difference_ax = axes[index - 1, 2]
        ax.plot(qa[0], qa[index], label="Newton", lw=1.6)
        ax.plot(qb[0], qb[index], label="Cellular", lw=1.3, ls="--")
        lo, hi = limits[index - 1]
        pad = 0.04 * max(hi - lo, 1e-12)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
        anomaly_ax.plot(qa[0], qa[index] - references[0][index - 1],
                        label="Newton", lw=1.6)
        anomaly_ax.plot(qb[0], qb[index] - references[1][index - 1],
                        label="Cellular", lw=1.3, ls="--")
        anomaly_ax.axhline(0, color="0.3", lw=0.7)
        alo, ahi = anomaly_limits[index - 1]
        apad = 0.06 * max(ahi - alo, 1e-12)
        anomaly_ax.set_ylim(alo - apad, ahi + apad)
        anomaly_ax.set_ylabel(f"Change in {ylabel.lower()}")
        anomaly_ax.grid(alpha=0.25)
        response_difference = ((qb[index] - references[1][index - 1])
                               - (qa[index] - references[0][index - 1]))
        difference_ax.plot(qa[0], response_difference, color="#7A3E9D", lw=1.4)
        difference_ax.axhline(0, color="0.3", lw=0.7)
        dlo, dhi = difference_limits[index - 1]
        dpad = 0.06 * max(dhi - dlo, 1e-12)
        difference_ax.set_ylim(dlo - dpad, dhi + dpad)
        difference_ax.set_ylabel(f"Cellular - Newton\nresponse difference")
        difference_ax.grid(alpha=0.25)
    axes[0, 0].legend(); axes[0, 1].legend()
    axes[0, 0].set_title("Absolute profiles")
    axes[0, 1].set_title("Change from t = 5,000 yr pre-step state")
    axes[0, 2].set_title("Difference between model responses")
    axes[-1, 0].set_xlabel("Distance from divide (m)")
    axes[-1, 1].set_xlabel("Distance from divide (m)")
    axes[-1, 2].set_xlabel("Distance from divide (m)")
    fig.suptitle(f"Roering strip comparison, t = {time:g} yr")
    fig.tight_layout()
    return fig


def render_hillslope_comparison_frame(newton, cellular, *, elevation_limits,
                                      residual_limit, laplacian_residual_limit,
                                      time, history=None, figsize=(15, 8)):
    """Render one paired comparison frame and return ``(figure, residual)``."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Hillslope animation needs matplotlib.") from exc
    a, transform = _raster(newton)
    b, transform_b = _raster(cellular)
    _validate_pair(a, transform, b, transform_b)
    residual = b - a
    laplacian_residual = _laplacian(b, transform_b) - _laplacian(a, transform)
    dx, dy = _cell_sizes(transform)
    extent = (0.0, dx * a.shape[1], 0.0, dy * a.shape[0])
    fig = plt.figure(figsize=figsize)
    grid = fig.add_gridspec(2, 3)
    axes = np.empty((2, 3), dtype=object)
    axes[0, 0] = fig.add_subplot(grid[0, 0])
    axes[0, 1] = fig.add_subplot(grid[0, 1])
    axes[1, 0] = fig.add_subplot(grid[1, 0])
    axes[1, 1] = fig.add_subplot(grid[1, 1])
    history_ax = fig.add_subplot(grid[:, 2])
    axes[0, 2] = history_ax
    axes[1, 2] = history_ax
    for ax, image, label in ((axes[0, 0], a, "Newton elevation"),
                             (axes[0, 1], b, "Cellular elevation")):
        ax.imshow(image, cmap="terrain", vmin=elevation_limits[0],
                  vmax=elevation_limits[1], extent=extent, origin="lower")
        ax.set_title(label)
    image = axes[1, 0].imshow(residual, cmap="RdBu_r", vmin=-residual_limit,
                               vmax=residual_limit, extent=extent, origin="lower")
    axes[1, 0].set_title("Elevation residual: cellular - Newton")
    fig.colorbar(image, ax=axes[1, 0], shrink=0.8, label="residual (m)")
    lap_image = axes[1, 1].imshow(
        laplacian_residual, cmap="PuOr_r", vmin=-laplacian_residual_limit,
        vmax=laplacian_residual_limit, extent=extent, origin="lower")
    axes[1, 1].set_title("Laplacian residual: cellular - Newton")
    fig.colorbar(lap_image, ax=axes[1, 1], shrink=0.8,
                 label="Laplacian residual (m$^{-1}$)")
    if history:
        history_ax.plot(history[0], history[1], label="RMS residual")
        history_ax.plot(history[0], history[2],
                        label="maximum absolute residual")
        history_ax.axvline(time, color="0.4", ls=":")
        history_ax.legend()
    history_ax.set_title(f"Elevation residual through time: {time:g} yr")
    history_ax.set_xlabel("Time (yr)")
    history_ax.set_ylabel("residual (m)")
    for ax in (axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_aspect("equal")
    for ax in set(axes.ravel()):
        ax.grid(alpha=0.2)
    fig.suptitle(f"Hillslope model comparison, t = {time:g} yr")
    fig.tight_layout()
    return fig, residual


def render_hillslope_comparison_animation(newton_dir, cellular_dir, out_gif,
                                          *, pattern="*.tif", fps=4, loop=0,
                                          stats_csv=None, keep_frames=False,
                                          layout="2d"):
    """Write a fixed-scale four-panel Newton/cellular comparison GIF.

    If ``stats_csv`` is supplied, its ``time_yr``, ``rms`` and ``linf`` columns
    are used for the time-series panel. Otherwise these metrics are calculated
    from the paired rasters. Frame keys are used as time labels when no CSV is
    supplied.
    """
    pairs = discover_frame_pairs(newton_dir, cellular_dir, pattern=pattern)
    if layout not in {"2d", "strip"}:
        raise ValueError("layout must be '2d' or 'strip'")
    zlim = _limits(pairs) if layout == "2d" else None
    strip_reference = min(pairs, key=lambda pair: abs(float(pair[2]) - 50.0))
    strip_limits = (_strip_limits(pairs, strip_reference)
                    if layout == "strip" else None)
    if stats_csv is not None:
        import pandas as pd
        table = pd.read_csv(stats_csv)
        required = {"time_yr", "rms", "linf"}
        if not required.issubset(table.columns):
            raise KeyError(f"stats CSV must contain {sorted(required)}")
        if "iteration" in table.columns:
            indexed = table.set_index(table["iteration"].astype(int), drop=False)
            requested = [int(pair[2]) for pair in pairs]
            missing = sorted(set(requested) - set(indexed.index))
            if missing:
                raise ValueError(f"stats CSV is missing frame iterations: {missing}")
            table = indexed.loc[requested].reset_index(drop=True)
        times = table.time_yr.to_numpy(float)
        rms = table.rms.to_numpy(float)
        linf = table.linf.to_numpy(float)
        if len(times) != len(pairs):
            raise ValueError("stats CSV and frame count do not match")
    else:
        times, rms, linf = [], [], []
        for left, right, key in pairs:
            a, _ = _raster(left); b, _ = _raster(right)
            d = (b - a); d = d[np.isfinite(d)]
            times.append(float(key)); rms.append(float(np.sqrt(np.mean(d*d))))
            linf.append(float(np.max(np.abs(d))))
        times, rms, linf = np.asarray(times), np.asarray(rms), np.asarray(linf)
    out_gif = Path(out_gif)
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    if keep_frames:
        frame_dir = out_gif.parent / f"{out_gif.stem}_frames"
        frame_dir.mkdir(parents=True, exist_ok=True)
        context = None
    else:
        context = tempfile.TemporaryDirectory()
        frame_dir = Path(context.name)
    try:
        pngs = []
        history = (times, rms, linf)
        for i, (left, right, _) in enumerate(pairs):
            if layout == "strip":
                fig = render_strip_comparison_frame(
                    left, right, limits=strip_limits[0],
                    anomaly_limits=strip_limits[1],
                    difference_limits=strip_limits[2], references=strip_limits[3],
                    time=float(times[i]))
            else:
                fig, _ = render_hillslope_comparison_frame(
                    left, right, elevation_limits=zlim[:2], residual_limit=zlim[2],
                    laplacian_residual_limit=zlim[3], time=float(times[i]),
                    history=history)
            png = frame_dir / f"frame_{i:04d}.png"
            fig.savefig(png, dpi=120)
            pngs.append(png)
            import matplotlib.pyplot as plt
            plt.close(fig)
        _write_gif(pngs, out_gif, fps=fps, loop=loop)
    finally:
        if context is not None:
            context.cleanup()
    return out_gif


def _write_gif(frame_pngs, out_gif, *, fps=4, loop=0):
    """Write a GIF without importing the PyGMT/rioxarray animation module."""
    try:
        from PIL import Image
    except ImportError as exc:  # pragma: no cover
        raise ImportError("GIF animation needs Pillow.") from exc
    images = [Image.open(path).convert("RGB") for path in frame_pngs]
    if not images:
        raise ValueError("No frames to write.")
    images[0].save(str(out_gif), save_all=True, append_images=images[1:],
                   duration=int(round(1000.0 / fps)), loop=loop, optimize=True)


animate_hillslope_comparison = render_hillslope_comparison_animation
