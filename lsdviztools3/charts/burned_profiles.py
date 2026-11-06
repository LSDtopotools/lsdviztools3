"""Long profiles of channels burned by lsdtt3 ``lsdtt-stream-burning``.

Reads ``<write_fname>_burned_profiles.csv`` from a run with
``stream_burn_mode: burn_to_corridor``: one row
per channel cell, downstream order within each path (a chain of reaches from
its head to where it joins another path or leaves the network), with

  path_id, path_rank, reach_id, name, point_index, distance_m, flow_distance_m,
  x, y, z_original, z_cleaned, z_profile, z_final, pit_filled, pit_id,
  is_record_low, is_rejected_low, drainage_area_m2, chi, joins_path_id,
  join_distance_m, starts_on_path_id

Everything here is drawing: which cells form a course (a path and the paths it
joins, down to the outlet) is read from ``joins_path_id`` / ``join_distance_m``;
nothing is recomputed. The one derived quantity is the chi axis of a course,
where each joined path's chi (counted from that path's head in the CSV) is
shifted so the course's chi is continuous at the joins.

``course_rows``       the rows of a course, ordered downstream, with
                      ``course_distance_m`` and ``course_chi``
``select_course``     a course by path id, path rank or river name
``draw_burned_profile`` one course on one axis
``plot_burned_profiles`` one panel per CSV (e.g. X = 0 / 60 / 150 m, or a
                      linear against a chi fit), shared x axis
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Fixed roles (validated together, light surface; green is below 3:1 contrast,
# so every band carries a legend label).
FINAL_COLOUR = "#2a78d6"       # the burned channel in the output DEM
CUT_COLOUR = "#eb6834"         # cut below the cleaned DEM
PIT_COLOUR = "#1baf7a"         # small pits filled before the profile was sampled
RAISE_COLOUR = "#9a5bd6"       # raised above the cleaned DEM (junction pins, rejected minima)
ORIGINAL_COLOUR = "#8a8a85"    # the input DEM under the channel
INK = "#3d3d3a"

REQUIRED = ["path_id", "path_rank", "name", "point_index", "distance_m", "flow_distance_m",
            "x", "y", "z_original", "z_cleaned", "z_profile", "z_final", "pit_filled",
            "pit_id", "is_record_low", "is_rejected_low", "chi", "joins_path_id",
            "join_distance_m"]


def read_burned_profiles(csv) -> pd.DataFrame:
    """Read a ``_burned_profiles.csv`` written by burn_to_corridor."""
    p = Path(csv)
    if not p.exists():
        raise FileNotFoundError(f"{p} (run lsdtt-stream-burning with stream_burn_mode: "
                                "burn_to_corridor and write_burned_profiles_csv: true).")
    df = pd.read_csv(p, keep_default_na=False, na_values={"chi": [""], "drainage_area_m2": [""],
                                                         "join_distance_m": [""]})
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{p}: not a burn_to_corridor profiles file (missing {missing}); "
                         "older runs (z_filled / z_burned) are not supported.")
    df["name"] = df["name"].astype(str)
    return df


def course_rows(df: pd.DataFrame, path_id: int, start_index: int = 0) -> pd.DataFrame:
    """Rows of the course from ``path_id`` (from its cell ``start_index``) down
    to the outlet, following ``joins_path_id``. Adds ``course_distance_m``
    (downstream from the course's first cell) and ``course_chi`` (continuous
    at the joins; NaN where the run has no chi)."""
    parts = []
    pid, lo = int(path_id), int(start_index)
    offset = 0.0          # added to the chi of the current path
    seen = set()
    while pid >= 0 and pid not in seen:
        seen.add(pid)
        g = df[df.path_id == pid].sort_values("point_index")
        if g.empty:
            raise ValueError(f"no path {pid} in the profiles file")
        g = g[g.point_index >= lo].copy()
        g["course_chi"] = g.chi.astype(float) + offset
        nxt = int(g.joins_path_id.iloc[0])
        if nxt >= 0:
            # the joined path continues from the join cell (this path's last
            # cell); its chi is shifted to match this path's value there
            h = df[df.path_id == nxt].sort_values("point_index")
            jd = float(g.join_distance_m.iloc[0])
            at = h[np.isclose(h.distance_m, jd, atol=1e-6)]
            if at.empty:
                raise ValueError(f"path {pid} joins path {nxt} at {jd} m, which has no such cell")
            lo = int(at.point_index.min())
            offset = float(g.course_chi.iloc[-1]) - float(at.chi.iloc[0])
            g = g.iloc[:-1]
        parts.append(g)
        pid = nxt
    c = pd.concat(parts, ignore_index=True)
    c["course_distance_m"] = c.flow_distance_m.iloc[0] - c.flow_distance_m
    return c


def select_course(df: pd.DataFrame, path_id=None, path_rank=None, name=None):
    """A course by path id, path rank (0 = the first, longest path) or river
    name (from the most upstream cell with that name). Returns (rows, label)."""
    given = [v is not None for v in (path_id, path_rank, name)]
    if sum(given) > 1:
        raise ValueError("give one of path_id, path_rank, name")
    if name is not None:
        hit = df[df.name == name]
        if hit.empty:
            raise ValueError(f"no cells named {name!r} (was vector_name_field set?)")
        top = hit.loc[hit.flow_distance_m.idxmax()]
        return course_rows(df, int(top.path_id), int(top.point_index)), name
    if path_rank is None and path_id is None:
        path_rank = 0
    if path_rank is not None:
        hit = df[df.path_rank == int(path_rank)]
        if hit.empty:
            raise ValueError(f"no path of rank {path_rank}")
        path_id = int(hit.path_id.iloc[0])
    rows = course_rows(df, int(path_id))
    names = list(dict.fromkeys(n for n in rows.name if n))
    label = (" / ".join(names[:3]) + (" ..." if len(names) > 3 else "")) if names else f"path {path_id}"
    return rows, label


def _joins_on_course(df, course, max_labels):
    """(course x position, name, length) of paths joining the course."""
    on = {(round(a, 3), round(b, 3)): i for i, (a, b) in enumerate(zip(course.x, course.y))}
    course_paths = set(course.path_id)
    out = []
    for pid, g in df.groupby("path_id"):
        if pid in course_paths:
            continue
        g = g.sort_values("point_index")
        last = g.iloc[-1]
        k = on.get((round(last.x, 3), round(last.y, 3)))
        if k is None or int(last.joins_path_id) < 0:
            continue
        nm = next((n for n in g.name[::-1] if n), "")
        out.append((k, nm, float(g.distance_m.max())))
    out.sort(key=lambda t: -t[2])
    return out, out[:max_labels]


def draw_burned_profile(ax, df: pd.DataFrame, course: pd.DataFrame, x: str = "distance",
                        mark=None, max_join_labels: int = 6, legend: bool = True):
    """Draw one course on ``ax``. ``x``: "distance" (km downstream) or "chi"."""
    if x == "chi":
        if not np.isfinite(course.course_chi).all():
            raise ValueError("this run has no chi (burn_profile_fit_algorithm_choice: chi)")
        xs = course.course_chi.to_numpy()
        ax.set_xlabel("chi (m)")
    else:
        xs = course.course_distance_m.to_numpy() / 1000.0
        ax.set_xlabel("distance downstream (km)")
    zo = course.z_original.to_numpy()
    zc = course.z_cleaned.to_numpy()
    zf = course.z_final.to_numpy()
    zp = course.z_profile.to_numpy()

    # Bands first, lines over them. Each cell's change against the input DEM is
    # split without overlap: cut (final below the DEM); pit fill (from the DEM
    # up to the lower of the cleaned and final values: the part of a small pit
    # the channel still sits in); other raise (above both: junction pins,
    # rejected minima).
    pit = zc > zo
    n_pits = int(course.loc[pit, "pit_id"].nunique())

    def widen(m):
        # each run widened by one cell either side (the band is zero there),
        # so a one-cell run still draws
        w = m.copy()
        w[1:] |= m[:-1]
        w[:-1] |= m[1:]
        return w

    cut = zf < zo
    pit_top = np.minimum(zc, zf)
    pit_band = pit & (zf > zo)
    raise_lo = np.maximum(zo, zc)
    raise_band = zf > raise_lo + 1e-6
    ax.fill_between(xs, zf, zo, where=widen(cut), color=CUT_COLOUR, alpha=0.35, lw=0)
    ax.fill_between(xs, zo, np.where(pit_band, pit_top, zo), where=widen(pit_band),
                    color=PIT_COLOUR, alpha=0.6, lw=0)
    ax.fill_between(xs, raise_lo, np.where(raise_band, zf, raise_lo), where=widen(raise_band),
                    color=RAISE_COLOUR, alpha=0.5, lw=0)
    ax.plot(xs, zo, color=ORIGINAL_COLOUR, lw=0.8, zorder=3)
    if np.nanmax(np.abs(zp - zf)) > 0.01:
        ax.plot(xs, zp, color=INK, lw=0.8, ls="--", zorder=4)
    ax.plot(xs, zf, color=FINAL_COLOUR, lw=1.6, zorder=5)
    lows = course.is_record_low.to_numpy() == 1
    ax.plot(xs[lows], zf[lows], "o", ms=3.0, color=INK, mec="white", mew=0.5, zorder=6, ls="none")
    rej = course.is_rejected_low.to_numpy() == 1
    if rej.any():
        ax.plot(xs[rej], zo[rej], "x", ms=4.5, color=INK, zorder=6, ls="none")

    joins, labelled = _joins_on_course(df, course, max_join_labels)
    for k, nm, _ in joins:
        ax.axvline(xs[k], color=ORIGINAL_COLOUR, lw=0.5, alpha=0.6, zorder=1)
    ytop = ax.get_ylim()[1]
    span = float(np.nanmax(xs) - np.nanmin(xs)) or 1.0
    placed = []
    for k, nm, _ in labelled:
        if nm and all(abs(xs[k] - q) > 0.03 * span for q in placed):
            ax.text(xs[k], ytop, " " + nm, rotation=90, va="top", ha="right", fontsize=7, color=INK)
            placed.append(xs[k])
    if mark is not None:
        mx, my, mlabel = mark
        k = int(np.argmin(np.hypot(course.x - mx, course.y - my)))
        ax.axvline(xs[k], color=INK, lw=0.8, ls=":", zorder=2)
        ax.text(xs[k], ytop, " " + mlabel, va="top", ha="left", fontsize=8, color=INK)
    ax.set_ylabel("elevation (m)")
    ax.grid(alpha=0.25, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    rises = int((np.diff(zf) > 0).sum())
    stats = (f"{n_pits} pits filled ({int(pit.sum())} cells), {int(lows.sum())} minima used, "
             f"{int(rej.sum())} rejected; final rises at {rises} of {len(zf) - 1} steps")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Line2D([], [], color=ORIGINAL_COLOUR, lw=0.8, label="input DEM under the channel"),
               Line2D([], [], color=FINAL_COLOUR, lw=1.6, label="burned channel (output DEM)"),
               Line2D([], [], color=INK, marker="o", ms=3, mec="white", ls="none",
                      label="minima the profile joins"),
               Patch(color=CUT_COLOUR, alpha=0.35, label="cut below the input DEM"),
               Patch(color=PIT_COLOUR, alpha=0.6, label="small pit filled before sampling (channel in it)"),
               Patch(color=RAISE_COLOUR, alpha=0.5, label="raised (junction, rejected minimum)")]
    if rej.any():
        handles.append(Line2D([], [], color=INK, marker="x", ls="none", label="minima rejected by the cut cap"))
    if np.nanmax(np.abs(zp - zf)) > 0.01:
        handles.append(Line2D([], [], color=INK, lw=0.8, ls="--", label="fitted profile"))
    if legend:
        ax.legend(handles=handles, fontsize=7, loc="upper right", frameon=False)
    draw_burned_profile.last_handles = handles
    return stats


def plot_burned_profiles(csvs, out, labels=None, path_id=None, path_rank=None, name=None,
                         x: str = "distance", mark=None, title=None, dpi: int = 200,
                         panel_height: float = 3.4, max_join_labels: int = 6):
    """One panel per CSV for the same course (chosen in the first CSV, by path
    id, path rank or name; by name or rank in the others, whose path ids may
    differ), shared x axis. Writes ``out`` (PNG)."""
    csvs = [csvs] if isinstance(csvs, (str, Path)) else list(csvs)
    labels = list(labels) if labels else [Path(c).name.replace("_burned_profiles.csv", "") for c in csvs]
    if len(labels) != len(csvs):
        raise ValueError("give one label per CSV")
    fig, axes = plt.subplots(len(csvs), 1, figsize=(11, panel_height * len(csvs)), sharex=True, squeeze=False)
    course_label = None
    for i, (csv, lab) in enumerate(zip(csvs, labels)):
        ax = axes[i, 0]
        df = read_burned_profiles(csv)
        if i > 0 and path_id is not None and name is None and path_rank is None:
            raise ValueError("compare several runs by name or path rank, not path id")
        course, course_label = select_course(df, path_id=path_id, path_rank=path_rank, name=name)
        stats = draw_burned_profile(ax, df, course, x=x, mark=mark, max_join_labels=max_join_labels,
                                    legend=False)
        ax.set_title(f"{lab}: {stats}", fontsize=9, loc="left")
        if i < len(csvs) - 1:
            ax.set_xlabel("")
    fig.suptitle(title or course_label, fontsize=11)
    handles = draw_burned_profile.last_handles
    fig.legend(handles=handles, loc="lower center", ncol=min(4, len(handles)), fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.06 if len(csvs) > 1 else 0.14, 1, 1))
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out
