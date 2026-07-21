"""Animate a coupled LEM run responding to a pulse of incision.

This is the transient-LEM animation for the case the plain elevation animation in
:mod:`lsdviztools3.render.animation` cannot tell: a landscape with Roering
nonlinear hillslopes and a channel network being driven by a *pulse* of incision
-- a window of rapid base-level fall (or elevated uplift) that sends a wave of
knickpoints up the network, which the hillslopes then relax toward.

The pulse is nearly invisible in elevation: a knickpoint is a few metres of
relief on a hillside that is hundreds of metres tall. It is obvious in the
**erosion rate**. So each frame draws three registered layers over one map, plus
a chart:

1. the elevation raster as an illuminated grey hillshade (the landscape),
2. the per-pixel erosion rate draped on top through a colour scale pinned across
   every frame (the pulse -- knickpoints are a bright band migrating upstream),
3. the channel network of that timestep as a solid overlay (the structure the
   pulse travels along), and
4. underneath, a time series of mean erosion rate with the pulse window shaded
   and a marker on the current frame, so a viewer can see where in the forcing
   history the frame sits.

The rasters all come from one ``lsdtt-lem`` run with ``frequency_of_output_writing``,
``write_erosion_rate_raster`` and ``write_channel_mask_raster`` set, driven by a
``forcing_history_fname`` CSV. This module only *renders* them -- every number it
draws was computed by lsdtt3.

Frame files, per iteration ``N``::

    <prefix>_model_iteration_<N>.tif                 elevation      (required)
    <prefix>_model_iteration_<N>_erosion_rate.tif    erosion rate   (optional)
    <prefix>_model_iteration_<N>_channel_mask.tif    channel mask   (optional)
    <prefix>_model_elevation_stats.csv               time + forcing (optional)

Missing optional layers degrade gracefully: without the erosion raster the frame
is a plain hillshade, without the channel mask there is no network overlay, and
without the stats CSV there is no chart panel.
"""

from __future__ import annotations

import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np

from ..io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
from .animation import write_gif
from .base import auto_projection, grid_to_display, set_grid_gtype
from .channels import render_drape
from .hillshade import render_hillshade
from .style import MapStyle

PathLike = Union[str, Path]

#: Captures the iteration index from a ``..._model_iteration_<N>...`` filename.
_ITER_RE = re.compile(r"_model_iteration_(\d+)")

#: Suffixes lsdtt-lem gives the two optional per-frame diagnostic rasters.
_EROSION_SUFFIX = "_erosion_rate"
_CHANNEL_SUFFIX = "_channel_mask"


@dataclass
class PulseFrame:
    """The rasters and model time belonging to one animation frame."""

    iteration: int
    elevation: Path
    erosion_rate: Optional[Path] = None
    channel_mask: Optional[Path] = None
    time_yr: Optional[float] = None

    def label(self) -> str:
        """Corner annotation: model time if known, else the iteration number."""
        if self.time_yr is None:
            return f"iter {self.iteration}"
        t = self.time_yr
        return f"t = {t/1000:.0f} kyr" if abs(t) >= 1000 else f"t = {t:.0f} yr"


def _read_stats(stats_csv: Optional[PathLike]):
    """Load the LEM stats CSV, or ``None`` if it is absent/unreadable."""
    if stats_csv is None:
        return None
    path = Path(stats_csv)
    if not path.is_file():
        return None

    import pandas as pd

    df = pd.read_csv(path)
    if "iteration" not in df.columns or "time_yr" not in df.columns:
        return None
    return df


def find_pulse_frames(
    directory: PathLike,
    *,
    prefix: Optional[str] = None,
    ext: str = ".tif",
    stats_csv: Optional[PathLike] = None,
) -> list[PulseFrame]:
    """Collect the per-iteration frames of an ``lsdtt-lem`` run.

    Args:
        directory: The run's output directory.
        prefix: ``write_fname`` of the run. ``None`` (default) accepts any prefix,
            which is what you want unless the directory holds several runs.
        ext: Raster extension (``".tif"``; pass ``".bil"`` for ENVI output).
        stats_csv: The run's ``_model_elevation_stats.csv``. ``None`` auto-detects
            a matching one in ``directory``; it supplies each frame's model time.

    Returns:
        Frames ordered by iteration number.

    Raises:
        FileNotFoundError: If no elevation frames are found.
    """
    directory = Path(directory)
    stem = prefix or "*"

    # The elevation frame is the one WITHOUT a diagnostic suffix, so match the
    # extension straight after the iteration number.
    elevations = [
        p
        for p in directory.glob(f"{stem}_model_iteration_*{ext}")
        if _ITER_RE.search(p.name)
        and p.name.endswith(f"{_ITER_RE.search(p.name).group(0)}{ext}")
    ]
    if not elevations:
        raise FileNotFoundError(
            f"No elevation frames matching {stem}_model_iteration_*{ext} in "
            f"{directory!r}. Run lsdtt-lem with frequency_of_output_writing > 0."
        )

    if stats_csv is None:
        found = sorted(directory.glob(f"{stem}_model_elevation_stats.csv"))
        stats_csv = found[0] if found else None
    stats = _read_stats(stats_csv)
    times = (
        {int(i): float(t) for i, t in zip(stats["iteration"], stats["time_yr"])}
        if stats is not None
        else {}
    )

    frames: list[PulseFrame] = []
    for elev in elevations:
        iteration = int(_ITER_RE.search(elev.name).group(1))
        base = elev.name[: -len(ext)]

        erosion = directory / f"{base}{_EROSION_SUFFIX}{ext}"
        channel = directory / f"{base}{_CHANNEL_SUFFIX}{ext}"

        frames.append(
            PulseFrame(
                iteration=iteration,
                elevation=elev,
                erosion_rate=erosion if erosion.is_file() else None,
                channel_mask=channel if channel.is_file() else None,
                time_yr=times.get(iteration),
            )
        )

    frames.sort(key=lambda f: f.iteration)
    return frames


def erosion_colour_limits(
    frames: Sequence[PulseFrame],
    *,
    percentile: float = 99.0,
) -> tuple:
    """Choose one erosion-rate colour range spanning every frame.

    A shared range is essential: the whole point of the animation is to compare
    erosion rates *between* frames, which a per-frame rescale would destroy (the
    pulse would look identical to the quiescent background).

    The top of the range is a high percentile rather than the maximum, because a
    single cell at a knickpoint lip can be an order of magnitude above the rest
    of the network and would flatten the colour scale onto nothing. The bottom is
    pinned at 0 -- deposition (negative rates) is not the subject here and is
    left to clamp onto the low end colour.

    Args:
        frames: The frames to scan (those without an erosion raster are skipped).
        percentile: Upper percentile of the pooled positive rates (default 99).

    Returns:
        ``(0.0, vmax)``.

    Raises:
        ValueError: If no frame carries an erosion raster with finite values.
    """
    pooled = []
    for frame in frames:
        if frame.erosion_rate is None:
            continue
        values = np.asarray(load_raster(frame.erosion_rate).values, dtype=float)
        values = values[np.isfinite(values)]
        pooled.append(values[values > 0])

    if not pooled or not any(p.size for p in pooled):
        raise ValueError("No frame has an erosion-rate raster with positive values.")

    vmax = float(np.percentile(np.concatenate(pooled), percentile))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1e-4
    return 0.0, vmax


def _add_channel_overlay(fig, mask_raster, coordinates, style, color):
    """Paint the channel pixels of a mask raster as a solid overlay."""
    import pygmt

    from .colormaps import _solid_cpt

    grid = to_gmt_grid(to_coordinates(load_raster(mask_raster), coordinates))
    crs = get_crs(grid)
    grid, extent = grid_to_display(grid, crs)
    set_grid_gtype(grid, crs)

    # Non-channel cells become NaN so nan_transparent leaves the map beneath them
    # untouched; the channels themselves all carry the same colour.
    values = np.asarray(grid.values, dtype=float)
    grid = grid.copy(data=np.where(values > 0.5, 1.0, np.nan))
    if not np.isfinite(grid.values).any():
        return fig  # no channels this timestep — nothing to draw

    projection = style.projection or auto_projection(extent, crs, style.figure_size)
    fig.grdimage(
        grid=grid,
        region=extent,
        projection=projection,
        cmap=_solid_cpt(color),
        nan_transparent=True,
    )
    return fig


def _true_runs(flags) -> list[tuple]:
    """Index ranges ``(first, last)`` of each contiguous run of True in ``flags``.

    Both ends are inclusive, so they index directly into the time array. Used to
    shade every interval in which the forcing was on -- a history may contain
    more than one pulse.
    """
    flags = np.asarray(flags, dtype=bool)
    if not flags.any():
        return []

    padded = np.concatenate(([False], flags, [False]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    return [
        (int(start), int(end - 1)) for start, end in zip(edges[0::2], edges[1::2])
    ]


def _add_timeseries_panel(
    fig,
    stats,
    *,
    current_time: Optional[float],
    style: MapStyle,
    height: str = "3c",
    gap: str = "1.4c",
):
    """Draw the mean-erosion-rate history under the map, marking this frame.

    The pulse window (the steps whose ``baselevel_fall_rate`` is non-zero) is
    shaded, so the viewer can see the forcing that produced what the map shows.
    """
    time_kyr = stats["time_yr"].to_numpy(dtype=float) / 1000.0
    erosion = stats["mean_erosion_rate"].to_numpy(dtype=float) * 1000.0  # mm/yr

    finite = np.isfinite(erosion)
    if not finite.any():
        return fig

    t0, t1 = float(time_kyr.min()), float(time_kyr.max())
    if t1 <= t0:
        t1 = t0 + 1.0
    ymax = float(np.nanmax(erosion[finite])) * 1.15
    if not np.isfinite(ymax) or ymax <= 0:
        ymax = 1.0

    width = style.figure_size
    fig.shift_origin(yshift=f"-{gap}")
    fig.basemap(
        region=[t0, t1, 0.0, ymax],
        projection=f"X{width}/{height}",
        frame=["WSne", "xaf+lmodel time (kyr)", "yaf+lmean erosion (mm/yr)"],
    )

    # Shade the interval(s) of base-level fall — the pulse itself.
    if "baselevel_fall_rate" in stats.columns:
        fall = stats["baselevel_fall_rate"].to_numpy(dtype=float)
        pulse = np.isfinite(fall) & (fall != 0)
        for start, end in _true_runs(pulse):
            fig.plot(
                x=[time_kyr[start], time_kyr[end], time_kyr[end],
                   time_kyr[start], time_kyr[start]],
                y=[0.0, 0.0, ymax, ymax, 0.0],
                fill="220/120/120",
                transparency=70,
            )

    # The uplift rate is the erosion rate the landscape relaxes back toward.
    if "uplift_rate" in stats.columns:
        uplift = stats["uplift_rate"].to_numpy(dtype=float) * 1000.0
        if np.isfinite(uplift).any():
            fig.plot(x=time_kyr, y=uplift, pen="0.75p,gray50,--")

    fig.plot(x=time_kyr[finite], y=erosion[finite], pen="1.25p,black")

    if current_time is not None:
        t_now = float(current_time) / 1000.0
        fig.plot(x=[t_now, t_now], y=[0.0, ymax], pen="0.9p,firebrick")
        idx = int(np.argmin(np.abs(time_kyr - t_now)))
        if finite[idx]:
            fig.plot(
                x=[t_now],
                y=[float(erosion[idx])],
                style="c0.18c",
                fill="firebrick",
                pen="0.5p,black",
            )
    return fig


def render_pulse_frame(
    frame: PulseFrame,
    *,
    coordinates: str = "local",
    erosion_clim: Optional[tuple] = None,
    erosion_cmap: str = "inferno",
    erosion_transparency: float = 30.0,
    channel_color: str = "dodgerblue",
    stats=None,
    style: Optional[MapStyle] = None,
    annotate: bool = True,
    colorbar: bool = True,
):
    """Render one frame: hillshade + erosion drape + channels + history panel.

    Args:
        frame: The frame's rasters and model time (see :func:`find_pulse_frames`).
        coordinates: Coordinate mode for every layer. See ``io.raster.to_coordinates``.
        erosion_clim: ``(vmin, vmax)`` erosion-rate colour limits, shared across
            the animation (see :func:`erosion_colour_limits`). ``None`` scales to
            this frame alone, which is almost never what you want in an animation.
        erosion_cmap: Colour map for the erosion-rate drape.
        erosion_transparency: Drape transparency, so relief shows through.
        channel_color: Colour of the channel-network overlay.
        stats: The run's stats DataFrame; supplies the time-series panel. ``None``
            draws no panel.
        style: Styling options shared by every frame.
        annotate: Stamp the model time on the map.
        colorbar: Draw the erosion-rate colour bar.

    Returns:
        The ``pygmt.Figure``.
    """
    style = style or MapStyle()

    if frame.erosion_rate is not None:
        vmin, vmax = erosion_clim if erosion_clim is not None else (None, None)
        fig = render_drape(
            frame.elevation,
            frame.erosion_rate,
            coordinates=coordinates,
            base_is_hillshade=False,  # the base is elevation: illuminate it
            cmap=erosion_cmap,
            vmin=vmin,
            vmax=vmax,
            out_of_bounds="clamp",
            transparency=erosion_transparency,
            colorbar=colorbar,
            colorbar_label="erosion rate (m/yr)",
            style=style,
        )
    else:
        fig = render_hillshade(
            frame.elevation, coordinates=coordinates, shade=True, style=style
        )

    if frame.channel_mask is not None:
        _add_channel_overlay(fig, frame.channel_mask, coordinates, style, channel_color)

    if annotate:
        fig.text(
            position="TL",
            text=frame.label(),
            offset="0.2c/-0.2c",
            justify="TL",
            no_clip=True,
            font="11p,Helvetica-Bold,black",
            fill="white@20",
            pen="0.5p,black",
            clearance="0.1c/0.1c+tO",
        )

    if stats is not None:
        _add_timeseries_panel(
            fig, stats, current_time=frame.time_yr, style=style
        )

    return fig


def animate_lem_pulse(
    directory: PathLike,
    out_gif: PathLike,
    *,
    prefix: Optional[str] = None,
    ext: str = ".tif",
    stats_csv: Optional[PathLike] = None,
    coordinates: str = "local",
    erosion_clim: Optional[tuple] = None,
    erosion_percentile: float = 99.0,
    erosion_cmap: str = "inferno",
    erosion_transparency: float = 30.0,
    channel_color: str = "dodgerblue",
    style: Optional[MapStyle] = None,
    show_timeseries: bool = True,
    annotate: bool = True,
    fps: float = 5.0,
    loop: int = 0,
    background: str = "white",
    keep_frames: bool = False,
) -> Path:
    """Animate an ``lsdtt-lem`` incision-pulse run as a GIF.

    Args:
        directory: The run's output directory (holds the ``_model_iteration_*``
            rasters and the stats CSV).
        out_gif: Output ``.gif`` path.
        prefix: The run's ``write_fname``. ``None`` accepts any prefix -- set it
            when the directory holds more than one run.
        ext: Raster extension (``".tif"``, or ``".bil"`` for ENVI output).
        stats_csv: The run's ``_model_elevation_stats.csv``; auto-detected when
            ``None``. Supplies frame times and the time-series panel.
        coordinates: Coordinate mode for every layer.
        erosion_clim: ``(vmin, vmax)`` erosion-rate colour limits shared by all
            frames. ``None`` derives them from every frame (see
            :func:`erosion_colour_limits`), which is what keeps the colour of a
            given erosion rate constant through the animation.
        erosion_percentile: Percentile setting the top of the derived colour
            range (ignored when ``erosion_clim`` is given).
        erosion_cmap: Colour map for the erosion-rate drape.
        erosion_transparency: Drape transparency, so relief shows through.
        channel_color: Colour of the channel-network overlay.
        style: Styling shared by every frame.
        show_timeseries: Draw the mean-erosion-rate panel under each map.
        annotate: Stamp the model time on each map.
        fps: Frames per second.
        loop: GIF loop count (``0`` = forever).
        background: Solid colour composited under each frame.
        keep_frames: Keep the PNG frames next to ``out_gif`` instead of using a
            temporary directory.

    Returns:
        The path to the written GIF.
    """
    style = style or MapStyle()
    directory = Path(directory)

    frames = find_pulse_frames(
        directory, prefix=prefix, ext=ext, stats_csv=stats_csv
    )

    if stats_csv is None:
        stem = prefix or "*"
        found = sorted(directory.glob(f"{stem}_model_elevation_stats.csv"))
        stats_csv = found[0] if found else None
    stats = _read_stats(stats_csv) if show_timeseries else None

    have_erosion = any(f.erosion_rate is not None for f in frames)
    if erosion_clim is None and have_erosion:
        erosion_clim = erosion_colour_limits(frames, percentile=erosion_percentile)

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
        for i, frame in enumerate(frames):
            fig = render_pulse_frame(
                frame,
                coordinates=coordinates,
                erosion_clim=erosion_clim,
                erosion_cmap=erosion_cmap,
                erosion_transparency=erosion_transparency,
                channel_color=channel_color,
                stats=stats,
                style=style,
                annotate=annotate,
            )
            png = frame_dir / f"frame_{i:04d}.png"
            fig.savefig(str(png), dpi=style.dpi, transparent=True)
            pngs.append(png)

        write_gif(pngs, out_gif, fps=fps, loop=loop, background=background)
    finally:
        if frame_dir_ctx is not None:
            frame_dir_ctx.cleanup()

    return out_gif
