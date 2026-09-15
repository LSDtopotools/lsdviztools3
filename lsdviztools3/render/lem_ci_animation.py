"""Animate an LEM run whose channelization index is walked through time.

A coupled ``lsdtt-lem`` run (fluvial and hillslope processes on every pixel) has
its drainage density set by the channelization index of Bonetti et al. (2020),
``CI = K l^(m+n) / (D^n U^(1-n))``. Stepping ``K`` or ``D`` with a forcing
history walks CI upward, and the network responds by dissecting further. This
module renders that as an animation with two registered panels per frame:

* a **narrow CI time-series panel on the left**, time on the vertical axis
  running downward, drawing the WHOLE run's CI history with a dot marking the
  current frame, so the viewer can see both where the run is going and where
  this frame sits in it; and
* a **large landscape panel on the right**, either

  - ``view="map"``: a shaded-relief map of the frame's elevation, or
  - ``view="oblique"``: a 3D perspective (``grdview``) of the same surface, so
    valleys and hillslopes are seen in relief.

Every frame of an animation shares one elevation colour range and, for the
oblique view, one vertical scale, taken over all frames. Without that each frame
would be rescaled to itself and relief change would be invisible.

The rasters and the CI series all come from lsdtt3. ``lsdtt-lem`` writes the CI it
applied at every step (``channelization_index`` in
``<prefix>_model_elevation_stats.csv``); this module only draws it.

Inputs, per run directory::

    <prefix>_initial_topography.tif      frame at t = 0            (optional)
    <prefix>_model_iteration_<N>.tif     elevation frames          (required)
    <prefix>_model_elevation_stats.csv   time + channelization_index (required)
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np

from ..io.raster import load_raster
from .animation import write_gif
from .base import new_figure
from .hillshade import render_hillshade
from .lem_pulse import PulseFrame, _read_stats, find_pulse_frames
from .relief3d import render_relief3d
from .style import MapStyle

PathLike = Union[str, Path]

#: Stats CSV column holding the channelization index lsdtt-lem applied.
CI_COLUMN = "channelization_index"


def find_ci_frames(
    directory: PathLike,
    *,
    prefix: Optional[str] = None,
    ext: str = ".tif",
    stats_csv: Optional[PathLike] = None,
    include_initial: bool = True,
) -> list[PulseFrame]:
    """Collect a run's frames in time order, optionally led by the t = 0 surface.

    Args:
        directory: The run's output directory.
        prefix: The run's ``write_fname``; ``None`` accepts any prefix.
        ext: Raster extension.
        stats_csv: The run's stats CSV (auto-detected when ``None``).
        include_initial: Prepend ``<prefix>_initial_topography`` as the t = 0
            frame when it exists, so the animation starts from the start.

    Returns:
        Frames ordered by model time.
    """
    directory = Path(directory)
    frames = find_pulse_frames(directory, prefix=prefix, ext=ext, stats_csv=stats_csv)
    if include_initial:
        stem = prefix or "*"
        initial = sorted(directory.glob(f"{stem}_initial_topography{ext}"))
        if initial:
            frames.insert(0, PulseFrame(iteration=0, elevation=initial[0], time_yr=0.0))
    return frames


def _stats_for(directory: Path, prefix: Optional[str], stats_csv: Optional[PathLike]):
    if stats_csv is None:
        found = sorted(directory.glob(f"{prefix or '*'}_model_elevation_stats.csv"))
        stats_csv = found[0] if found else None
    stats = _read_stats(stats_csv)
    if stats is None or CI_COLUMN not in stats.columns:
        raise ValueError(
            f"No {CI_COLUMN!r} column in the stats CSV ({stats_csv}). It is written "
            "by lsdtt-lem builds that support K/D forcing; re-run the model with a "
            "current lsdtt-lem."
        )
    return stats


def elevation_limits(frames: Sequence[PulseFrame]) -> tuple:
    """``(zmin, zmax)`` over every frame, for a colour scale shared by all frames."""
    lo, hi = np.inf, -np.inf
    for fr in frames:
        values = np.asarray(load_raster(fr.elevation).values, dtype=float)
        finite = values[np.isfinite(values)]
        if finite.size:
            lo, hi = min(lo, float(finite.min())), max(hi, float(finite.max()))
    if not np.isfinite(lo) or not np.isfinite(hi):
        raise ValueError("None of the frames has finite elevations.")
    return (lo, hi if hi > lo else lo + 1.0)


def _time_label(t_yr: float) -> str:
    return f"t = {t_yr / 1e6:.2f} Myr" if abs(t_yr) >= 1e5 else f"t = {t_yr / 1e3:.0f} kyr"


def add_ci_timeseries_panel(
    fig,
    stats,
    *,
    current_time: Optional[float],
    width: str = "3.2c",
    height: str = "10c",
    line_pen: str = "1.5p,black",
    marker_fill: str = "firebrick",
):
    """Draw the CI history as a narrow vertical panel at the current origin.

    Time runs DOWN the vertical axis and CI across the horizontal one. The full
    series is always drawn; ``current_time`` adds a dot on it and a faint
    horizontal guide, so the frame's position in the run is readable at a glance.

    Args:
        fig: ``pygmt.Figure`` to draw onto (at its current origin).
        stats: The run's stats DataFrame with ``time_yr`` and ``channelization_index``.
        current_time: Model time of this frame in years, or ``None`` for no marker.
        width, height: Panel size (GMT lengths).
        line_pen: Pen for the CI series.
        marker_fill: Fill of the current-time dot.

    Returns:
        The figure.
    """
    t_myr = stats["time_yr"].to_numpy(dtype=float) / 1e6
    ci = stats[CI_COLUMN].to_numpy(dtype=float)
    ok = np.isfinite(t_myr) & np.isfinite(ci)
    if not ok.any():
        raise ValueError("The CI series has no finite values.")
    t_myr, ci = t_myr[ok], ci[ok]

    # The stats row is logged at the END of each step; prepend the start so the
    # series spans the run from t = 0.
    if t_myr[0] > 0.0:
        t_myr = np.concatenate(([0.0], t_myr))
        ci = np.concatenate(([ci[0]], ci))

    t_end = float(t_myr.max())
    ci_lo, ci_hi = float(ci.min()), float(ci.max())
    pad = 0.08 * max(ci_hi - ci_lo, 1.0)

    # Negative height flips the y axis so time increases downward.
    fig.basemap(
        region=[ci_lo - pad, ci_hi + pad, 0.0, t_end],
        projection=f"X{width}/-{height}",
        frame=["WSne", "xa+lCI", "yaf+ltime (Myr)"],
    )
    fig.plot(x=ci, y=t_myr, pen=line_pen)

    if current_time is not None:
        t_now = float(current_time) / 1e6
        ci_now = float(np.interp(t_now, t_myr, ci))
        fig.plot(x=[ci_lo - pad, ci_hi + pad], y=[t_now, t_now], pen="0.5p,gray60,-")
        fig.plot(x=[ci_now], y=[t_now], style="c0.28c", fill=marker_fill, pen="0.6p,black")
    return fig


def render_ci_frame(
    frame: PulseFrame,
    stats,
    *,
    view: str = "map",
    zlimits: Optional[tuple] = None,
    coordinates: str = "local",
    map_width: str = "15c",
    panel_width: str = "3.2c",
    panel_gap: str = "2.2c",
    oblique_azimuth: float = 200.0,
    oblique_elevation: float = 35.0,
    vertical_exaggeration: float = 2.0,
    cmap: str = "gray",
    style: Optional[MapStyle] = None,
    annotate: bool = True,
):
    """Render one frame: CI panel on the left, landscape on the right.

    Args:
        frame: The frame (see :func:`find_ci_frames`).
        stats: The run's stats DataFrame (must hold ``channelization_index``).
        view: ``"map"`` (shaded relief) or ``"oblique"`` (3D perspective).
        zlimits: ``(zmin, zmax)`` shared by all frames (see :func:`elevation_limits`).
            ``None`` scales to this frame alone.
        coordinates: Coordinate mode for the landscape panel.
        map_width: Width of the landscape panel.
        panel_width: Width of the CI panel.
        panel_gap: Gap between the CI panel and the landscape.
        oblique_azimuth, oblique_elevation: Camera bearing and elevation angle for
            ``view="oblique"``. The defaults look from the south-south-west, low
            enough to see valley walls.
        vertical_exaggeration: True vertical:horizontal ratio for the oblique view.
        cmap: Colour map for the landscape.
        style: Base style; ``figure_size`` is overridden by ``map_width``.
        annotate: Stamp the model time on the landscape panel.

    Returns:
        The ``pygmt.Figure``.
    """
    if view not in ("map", "oblique"):
        raise ValueError(f"view must be 'map' or 'oblique', not {view!r}")

    base = style or MapStyle()
    # Opaque: GIF has no alpha, and a transparent style makes the grey colour map
    # take its anti-NaN nudge, which tints the lowest elevations blue.
    map_style = MapStyle(**{**base.__dict__, "figure_size": map_width, "cmap": cmap,
                            "transparent": False})

    # The CI panel matches the landscape panel's height: the raster aspect for
    # the map, a fixed fraction of the width for the tilted oblique view.
    da = load_raster(frame.elevation)
    ny, nx = int(da.shape[-2]), int(da.shape[-1])
    map_w_cm = float(str(map_width).rstrip("c"))
    panel_h_cm = map_w_cm * ny / nx if view == "map" else 0.62 * map_w_cm

    fig = new_figure(map_style)
    add_ci_timeseries_panel(
        fig, stats, current_time=frame.time_yr,
        width=panel_width, height=f"{panel_h_cm:.3g}c",
    )

    gap_cm = float(str(panel_gap).rstrip("c"))
    fig.shift_origin(xshift=f"{float(str(panel_width).rstrip('c')) + gap_cm:.3g}c")

    if view == "map":
        render_hillshade(frame.elevation, coordinates=coordinates, shade=True,
                         style=map_style, fig=fig, clim=zlimits)
    else:
        render_relief3d(frame.elevation, coordinates=coordinates,
                        azimuth=oblique_azimuth, elevation=oblique_elevation,
                        vertical_exaggeration=vertical_exaggeration,
                        shade=True, cmap=cmap, style=map_style, fig=fig,
                        zlimits=zlimits)

    if annotate and frame.time_yr is not None:
        fig.text(
            position="TL", text=_time_label(frame.time_yr),
            offset="0.2c/-0.2c", justify="TL", no_clip=True,
            font="12p,Helvetica-Bold,black", fill="white@20",
            pen="0.5p,black", clearance="0.1c/0.1c+tO",
        )
    return fig


def animate_ci_run(
    directory: PathLike,
    out_gif: PathLike,
    *,
    view: str = "map",
    prefix: Optional[str] = None,
    ext: str = ".tif",
    stats_csv: Optional[PathLike] = None,
    include_initial: bool = True,
    coordinates: str = "local",
    map_width: str = "15c",
    oblique_azimuth: float = 200.0,
    oblique_elevation: float = 35.0,
    vertical_exaggeration: float = 2.0,
    cmap: str = "gray",
    fps: float = 10.0,
    dpi: int = 110,
    loop: int = 0,
    background: str = "white",
    keep_frames: bool = False,
    style: Optional[MapStyle] = None,
) -> Path:
    """Animate a CI-ladder ``lsdtt-lem`` run as a GIF (see module docstring).

    Args:
        directory: Run output directory.
        out_gif: Output ``.gif`` path.
        view: ``"map"`` or ``"oblique"``.
        prefix: The run's ``write_fname`` (``None`` = any).
        ext: Raster extension.
        stats_csv: Stats CSV; auto-detected when ``None``.
        include_initial: Start from the t = 0 surface when it exists.
        coordinates: Coordinate mode.
        map_width: Landscape panel width.
        oblique_azimuth, oblique_elevation, vertical_exaggeration: Oblique camera.
        cmap: Landscape colour map.
        fps: Frames per second (duration = n_frames / fps).
        dpi: Frame resolution. GIFs grow fast with it; ~110 keeps a 100-frame
            animation to a sensible size.
        loop: GIF loop count (0 = forever).
        background: Colour composited under each frame.
        keep_frames: Keep the PNG frames beside the GIF.
        style: Base style.

    Returns:
        The path to the written GIF.
    """
    directory = Path(directory)
    frames = find_ci_frames(directory, prefix=prefix, ext=ext,
                            stats_csv=stats_csv, include_initial=include_initial)
    stats = _stats_for(directory, prefix, stats_csv)
    zlimits = elevation_limits(frames)

    out_gif = Path(out_gif)
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    if keep_frames:
        ctx = None
        frame_dir = out_gif.parent / f"{out_gif.stem}_frames"
        frame_dir.mkdir(parents=True, exist_ok=True)
    else:
        ctx = tempfile.TemporaryDirectory()
        frame_dir = Path(ctx.name)

    try:
        pngs = []
        for i, frame in enumerate(frames):
            fig = render_ci_frame(
                frame, stats, view=view, zlimits=zlimits, coordinates=coordinates,
                map_width=map_width, oblique_azimuth=oblique_azimuth,
                oblique_elevation=oblique_elevation,
                vertical_exaggeration=vertical_exaggeration, cmap=cmap, style=style,
            )
            png = frame_dir / f"frame_{i:04d}.png"
            fig.savefig(str(png), dpi=dpi, transparent=False)
            pngs.append(png)
        write_gif(pngs, out_gif, fps=fps, loop=loop, background=background)
    finally:
        if ctx is not None:
            ctx.cleanup()
    return out_gif
