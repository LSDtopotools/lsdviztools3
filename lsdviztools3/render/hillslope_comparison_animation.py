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
    zmin, zmax, rmax = np.inf, -np.inf, 0.0
    for newton, cellular, _ in pairs:
        a, _ = _raster(newton)
        b, _ = _raster(cellular)
        if a.shape != b.shape:
            raise ValueError(f"Paired frame shapes differ: {newton} and {cellular}")
        values = np.concatenate((a[np.isfinite(a)], b[np.isfinite(b)]))
        if values.size:
            zmin, zmax = min(zmin, values.min()), max(zmax, values.max())
        residual = b - a
        finite = residual[np.isfinite(residual)]
        if finite.size:
            rmax = max(rmax, float(np.max(np.abs(finite))))
    if not np.isfinite(zmin):
        raise ValueError("No finite values in paired frames")
    if zmin == zmax:
        zmax = zmin + 1.0
    if rmax == 0:
        rmax = 1.0
    return float(zmin), float(zmax), float(rmax)


def render_hillslope_comparison_frame(newton, cellular, *, elevation_limits,
                                      residual_limit, time, history=None,
                                      figsize=(12, 8)):
    """Render one paired comparison frame and return ``(figure, residual)``."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Hillslope animation needs matplotlib.") from exc
    a, transform = _raster(newton)
    b, _ = _raster(cellular)
    if a.shape != b.shape:
        raise ValueError("Newton and cellular frame shapes differ")
    residual = b - a
    extent = (transform.c, transform.c + transform.a * a.shape[1],
              transform.f + transform.e * a.shape[0], transform.f)
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    for ax, image, label in ((axes[0, 0], a, "Newton elevation"),
                             (axes[0, 1], b, "Cellular elevation")):
        ax.imshow(image, cmap="terrain", vmin=elevation_limits[0],
                  vmax=elevation_limits[1], extent=extent)
        ax.set_title(label)
    image = axes[1, 0].imshow(residual, cmap="RdBu_r", vmin=-residual_limit,
                              vmax=residual_limit, extent=extent)
    axes[1, 0].set_title("Cellular - Newton residual")
    fig.colorbar(image, ax=axes[1, 0], shrink=0.8, label="m")
    if history:
        axes[1, 1].plot(history[0], history[1], label="RMS")
        axes[1, 1].plot(history[0], history[2], label="Linf")
        axes[1, 1].axvline(time, color="0.4", ls=":")
        axes[1, 1].legend()
    axes[1, 1].set_title(f"Residual through time: {time:g} yr")
    axes[1, 1].set_xlabel("Time (yr)")
    axes[1, 1].set_ylabel("m")
    for ax in axes.ravel():
        ax.set_aspect("equal" if ax is not axes[1, 1] else "auto")
        ax.grid(alpha=0.2)
    fig.suptitle(f"Hillslope model comparison, t = {time:g} yr")
    fig.tight_layout()
    return fig, residual


def render_hillslope_comparison_animation(newton_dir, cellular_dir, out_gif,
                                          *, pattern="*.tif", fps=4, loop=0,
                                          stats_csv=None, keep_frames=False):
    """Write a fixed-scale four-panel Newton/cellular comparison GIF.

    If ``stats_csv`` is supplied, its ``time_yr``, ``rms`` and ``linf`` columns
    are used for the time-series panel. Otherwise these metrics are calculated
    from the paired rasters. Frame keys are used as time labels when no CSV is
    supplied.
    """
    pairs = discover_frame_pairs(newton_dir, cellular_dir, pattern=pattern)
    zlim = _limits(pairs)
    if stats_csv is not None:
        import pandas as pd
        table = pd.read_csv(stats_csv)
        required = {"time_yr", "rms", "linf"}
        if not required.issubset(table.columns):
            raise KeyError(f"stats CSV must contain {sorted(required)}")
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
            fig, _ = render_hillslope_comparison_frame(
                left, right, elevation_limits=zlim[:2], residual_limit=zlim[2],
                time=float(times[i]), history=history)
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
