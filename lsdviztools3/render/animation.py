"""Animate a transient model run as a GIF.

A transient ``lsdtt-lem`` run writes a sequence of intermediate elevation
rasters named ``<prefix>_model_iteration_<N>.{tif,bil}`` (with ``N`` the LEM
iteration number), alongside a ``<prefix>_model_elevation_stats.csv`` whose
``iteration`` / ``time_yr`` columns give the model time of each frame. This
module turns that sequence into an animated GIF: each frame is rendered with the
shared PyGMT hillshade layer on a single, consistent colour scale (so the colour
of a given elevation does not flicker between frames), optionally annotated with
the iteration number and model time, and the frames are assembled into a GIF
with Pillow.

The per-frame renderer (:func:`render_transient_frame`) returns a
``pygmt.Figure`` like the other ``render_*`` functions; the orchestrator
(:func:`animate_transient_rasters`) writes PNG frames to a temporary directory
and stitches them into the output GIF.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np

from ..io.raster import load_raster, to_coordinates, to_gmt_grid
from .hillshade import render_hillshade
from .style import MapStyle

PathLike = Union[str, Path]

#: Default glob for the intermediate rasters lsdtt-lem writes each ``write_freq``.
DEFAULT_FRAME_GLOB = "*_model_iteration_*.tif"

#: Captures the integer iteration index from a ``..._model_iteration_<N>.<ext>`` name.
_ITER_RE = re.compile(r"_model_iteration_(\d+)\b")


def _iteration_index(path: PathLike) -> int:
    """Return the LEM iteration number embedded in a frame filename.

    Falls back to ``-1`` when the ``_model_iteration_<N>`` token is absent, so
    that an arbitrary file list still sorts deterministically (by name).
    """
    m = _ITER_RE.search(Path(path).name)
    return int(m.group(1)) if m else -1


def find_transient_frames(
    directory: PathLike,
    *,
    glob: str = DEFAULT_FRAME_GLOB,
) -> list[Path]:
    """Find and order the transient frame rasters in a directory.

    Args:
        directory: Folder holding the ``lsdtt-lem`` output (its ``output/`` dir).
        glob: Filename pattern for the frames. The default matches the
            ``_model_iteration_<N>`` rasters; pass e.g. ``"*_model_iteration_*.bil"``
            for ENVI output.

    Returns:
        Frame paths sorted by iteration number (then by name).
    """
    frames = sorted(Path(directory).glob(glob), key=lambda p: (_iteration_index(p), p.name))
    if not frames:
        raise FileNotFoundError(
            f"No transient frames matching {glob!r} in {directory!r}."
        )
    return frames


def _time_lookup(stats_csv: Optional[PathLike]) -> dict[int, float]:
    """Map iteration number -> model time (yr) from a stats CSV, or ``{}``."""
    if stats_csv is None:
        return {}
    path = Path(stats_csv)
    if not path.is_file():
        return {}
    import pandas as pd

    df = pd.read_csv(path)
    if "iteration" not in df.columns or "time_yr" not in df.columns:
        return {}
    return {int(it): float(t) for it, t in zip(df["iteration"], df["time_yr"])}


def _global_clim(frames: Sequence[PathLike], coordinates: str) -> tuple:
    """Return ``(zmin, zmax)`` over every frame, for a shared colour scale."""
    zmin = np.inf
    zmax = -np.inf
    for f in frames:
        grid = to_gmt_grid(to_coordinates(load_raster(f), coordinates))
        values = np.asarray(grid.values, dtype=float)
        if np.isfinite(values).any():
            zmin = min(zmin, float(np.nanmin(values)))
            zmax = max(zmax, float(np.nanmax(values)))
    if not np.isfinite(zmin) or not np.isfinite(zmax):
        raise ValueError("None of the transient frames have finite values.")
    if zmin == zmax:
        zmax = zmin + 1.0
    return zmin, zmax


def _frame_label(path: PathLike, times: dict[int, float]) -> Optional[str]:
    """Build the corner annotation for a frame (iteration and time, if known)."""
    it = _iteration_index(path)
    if it < 0:
        return None
    if it in times:
        t = times[it]
        t_str = f"{t/1000:.1f} kyr" if abs(t) >= 1000 else f"{t:.0f} yr"
        return f"iter {it}  |  t = {t_str}"
    return f"iter {it}"


def write_gif(
    frame_pngs: Sequence[PathLike],
    out_gif: PathLike,
    *,
    fps: float = 4.0,
    loop: int = 0,
    background: str = "white",
) -> Path:
    """Stitch rendered PNG frames into an animated GIF.

    The frames are composited onto an opaque background first: the render layer
    writes transparent PNGs, and GIF has no alpha channel, so without this the
    transparent nodata areas come out as whatever the palette happens to put in
    slot 0.

    Args:
        frame_pngs: Ordered PNG frame paths.
        out_gif: Output ``.gif`` path.
        fps: Frames per second.
        loop: GIF loop count (``0`` = loop forever).
        background: Solid colour composited under each frame.

    Returns:
        The path to the written GIF.
    """
    from PIL import Image

    frame_pngs = [Path(p) for p in frame_pngs]
    if not frame_pngs:
        raise ValueError("No frames to write.")

    images = []
    for png in frame_pngs:
        rgba = Image.open(png).convert("RGBA")
        bg = Image.new("RGBA", rgba.size, background)
        images.append(
            Image.alpha_composite(bg, rgba).convert("P", palette=Image.ADAPTIVE)
        )

    out_gif = Path(out_gif)
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    images[0].save(
        str(out_gif),
        save_all=True,
        append_images=images[1:],
        duration=int(round(1000.0 / fps)),
        loop=loop,
        disposal=2,
        optimize=True,
    )
    return out_gif


def render_transient_frame(
    raster: PathLike,
    *,
    coordinates: str = "local",
    clim: Optional[tuple] = None,
    style: Optional[MapStyle] = None,
    label: Optional[str] = None,
):
    """Render one transient frame as a hillshade with an optional time label.

    A thin wrapper over :func:`lsdviztools3.render.hillshade.render_hillshade`
    that pins the colour scale (``clim``) and stamps a top-left annotation.

    Args:
        raster: Path to a single timestep elevation raster.
        coordinates: Coordinate system, see ``io.raster.to_coordinates``.
        clim: ``(zmin, zmax)`` shared colour limits across the animation.
        style: Styling options (colour map, figure size, fonts).
        label: Corner annotation (e.g. ``"iter 7 | t = 0.4 kyr"``). ``None``
            draws no label.

    Returns:
        The ``pygmt.Figure``.
    """
    style = style or MapStyle()
    fig = render_hillshade(
        raster,
        coordinates=coordinates,
        shade=True,
        style=style,
        clim=clim,
    )
    if label:
        # Place the label in the top-left of the current map frame. The region
        # and projection set by render_hillshade are still active, so a
        # reference-point position ("TL") with a small inset offset works.
        fig.text(
            position="TL",
            text=label,
            offset="0.2c/-0.2c",
            justify="TL",
            no_clip=True,
            font="10p,Helvetica-Bold,black",
            fill="white@20",
            pen="0.5p,black",
            clearance="0.1c/0.1c+tO",
        )
    return fig


def animate_transient_rasters(
    frames: Union[PathLike, Sequence[PathLike]],
    out_gif: PathLike,
    *,
    glob: str = DEFAULT_FRAME_GLOB,
    stats_csv: Optional[PathLike] = None,
    coordinates: str = "local",
    style: Optional[MapStyle] = None,
    clim: Optional[tuple] = None,
    annotate: bool = True,
    fps: float = 4.0,
    loop: int = 0,
    background: str = "white",
    keep_frames: bool = False,
) -> Path:
    """Assemble a transient raster sequence into an animated GIF.

    Each frame is rendered with :func:`render_transient_frame` on one shared
    colour scale, then the PNGs are stitched into a GIF with Pillow.

    Args:
        frames: A directory of frame rasters (globbed with ``glob`` and ordered
            by iteration number), or an explicit ordered sequence of raster paths.
        out_gif: Output ``.gif`` path.
        glob: Filename pattern used when ``frames`` is a directory.
        stats_csv: ``_model_elevation_stats.csv`` for the run; its
            ``iteration``/``time_yr`` columns drive the time annotation. When
            ``None`` and ``frames`` is a directory, a sibling ``*_model_elevation_stats.csv``
            is auto-detected.
        coordinates: Coordinate system, see ``io.raster.to_coordinates``.
        style: Styling options shared by every frame.
        clim: ``(zmin, zmax)`` colour limits. ``None`` computes a single range
            spanning all frames so the colour scale is consistent throughout.
        annotate: Stamp each frame with its iteration number and model time.
        fps: Frames per second of the GIF.
        loop: GIF loop count (``0`` = loop forever).
        background: Solid colour composited under each (transparent) frame.
        keep_frames: Keep the intermediate PNG frames next to ``out_gif`` instead
            of using a temporary directory.

    Returns:
        The path to the written GIF.
    """
    style = style or MapStyle()

    # Resolve the ordered frame list.
    if isinstance(frames, (str, Path)):
        directory = Path(frames)
        frame_paths = find_transient_frames(directory, glob=glob)
        if stats_csv is None:
            cand = sorted(directory.glob("*_model_elevation_stats.csv"))
            stats_csv = cand[0] if cand else None
    else:
        frame_paths = [Path(f) for f in frames]
        if not frame_paths:
            raise ValueError("No frames supplied.")

    times = _time_lookup(stats_csv) if annotate else {}

    if clim is None:
        clim = _global_clim(frame_paths, coordinates)

    out_gif = Path(out_gif)
    out_gif.parent.mkdir(parents=True, exist_ok=True)

    if keep_frames:
        frame_dir_ctx = None
        frame_dir = out_gif.parent / (out_gif.stem + "_frames")
        frame_dir.mkdir(parents=True, exist_ok=True)
    else:
        frame_dir_ctx = tempfile.TemporaryDirectory()
        frame_dir = Path(frame_dir_ctx.name)

    try:
        pngs = []
        for i, fpath in enumerate(frame_paths):
            label = _frame_label(fpath, times) if annotate else None
            fig = render_transient_frame(
                fpath,
                coordinates=coordinates,
                clim=clim,
                style=style,
                label=label,
            )
            png = frame_dir / f"frame_{i:04d}.png"
            fig.savefig(str(png), dpi=style.dpi, transparent=True)
            pngs.append(png)

        write_gif(pngs, out_gif, fps=fps, loop=loop, background=background)
    finally:
        if frame_dir_ctx is not None:
            frame_dir_ctx.cleanup()

    return out_gif
