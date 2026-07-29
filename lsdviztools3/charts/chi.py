"""Chi-profile plotting (matplotlib).

Two views over ``lsdtt-chi-analysis`` fit output:

* :func:`plot_chi_profile` — a single basin's chi-elevation profile.
* :func:`plot_chi_stacked` — many basins on one axis, each normalised to start
  at chi = 0 and shifted right by a cumulative offset (the "stacked chi plot"
  of lsdviztools), coloured globally by ``k_sn``.

Both accept **either** representation of the fit, so plotting works whether the
analysis was written as CSV or FlatGeobuf:

* a **nodes** table — the CSV ``*_chi_fit_nodes_mn_*.csv`` with per-node
  ``observed_elevation``/``fitted_elevation`` and a ``k_sn`` column attached by
  :func:`lsdviztools3.io.chi.assign_node_ksn`. Drawn as observed points coloured
  by k_sn with the fitted profile overlaid.
* a **segments** table — the CSV ``*_chi_fit_segments_mn_*.csv`` or the fgb
  ``*_chi_fit_segment_lines_mn_*.fgb`` (via
  :func:`lsdviztools3.io.chi.load_chi_segment_lines`). Each fitted segment is
  drawn as a line in (chi, elevation) coloured by its own ``k_sn`` — no observed
  scatter, and no node->segment join needed.

The plotting engine is matplotlib (DESIGN decision #1), imported lazily so the
core maps install stays lean.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

_MPL_MSG = (
    "Chi plotting needs matplotlib. Install the optional 'charts' extra, "
    "e.g. `pip install lsdviztools3[charts]` or `pip install matplotlib`."
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
        from matplotlib import colors  # noqa: F401
        from matplotlib.collections import LineCollection  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without mpl
        raise ImportError(_MPL_MSG) from exc
    return plt, colors, LineCollection


def _is_segments(data: pd.DataFrame) -> bool:
    """A segments table has chi/elevation segment endpoints; nodes do not."""
    return {"downstream_chi", "upstream_chi"}.issubset(data.columns)


def _colour_values(values: np.ndarray, log_scale: bool) -> np.ndarray:
    """Return colour values, optionally log10-scaled (floored like lsdviztools)."""
    v = np.asarray(values, dtype=float)
    if not log_scale:
        return v
    # Mirror lsdviztools: floor tiny/zero k_sn to avoid -inf in log space.
    out = np.full_like(v, np.nan)
    positive = v >= 0.1
    out[positive] = np.log10(v[positive])
    out[(~positive) & np.isfinite(v)] = 0.0
    return out


def _colour_label(color_by: str, log_scale: bool) -> str:
    if color_by == "k_sn":
        return r"$\log_{10}\,k_{sn}$" if log_scale else r"$k_{sn}$"
    return f"log10 {color_by}" if log_scale else color_by


def _basin_chi_bounds(data: pd.DataFrame, segments: bool):
    """Return (min_chi, max_chi) for a single-basin frame."""
    if segments:
        return float(data["downstream_chi"].min()), float(data["upstream_chi"].max())
    return float(data["chi"].min()), float(data["chi"].max())


def _basin_column(df: pd.DataFrame) -> str:
    """Basin grouping column: ``basin_index`` (fit tables) or ``basin_key`` (data maps)."""
    for name in ("basin_index", "basin_key"):
        if name in df.columns:
            return name
    raise KeyError(
        "Frame has no 'basin_index' or 'basin_key' column; "
        f"columns are {list(df.columns)}."
    )


def _scatter_observed(ax, observed, basin, shift, *, color, size):
    """Overlay observed chi-elevation points for one basin at the given x-shift.

    ``observed`` is a chi data-map frame (e.g. from
    :func:`lsdviztools3.io.chi.load_chi_points`) with ``chi``/``elevation`` and a
    basin column. Coloured neutrally because the data-map points carry no k_sn.
    """
    col = _basin_column(observed)
    pts = observed[observed[col] == basin]
    if pts.empty:
        return
    ax.scatter(
        pts["chi"].to_numpy(dtype=float) + shift, pts["elevation"],
        s=size, c=color, edgecolors="none", zorder=0,
    )


def _segment_lines(data: pd.DataFrame, x_shift: float):
    """Build matplotlib line segments [(chi0,z0),(chi1,z1)] in (chi, elevation)."""
    d_chi = data["downstream_chi"].to_numpy(dtype=float) + x_shift
    u_chi = data["upstream_chi"].to_numpy(dtype=float) + x_shift
    d_z = data["downstream_elevation"].to_numpy(dtype=float)
    u_z = data["upstream_elevation"].to_numpy(dtype=float)
    return [
        [(d_chi[i], d_z[i]), (u_chi[i], u_z[i])] for i in range(len(data))
    ]


def plot_chi_profile(
    data: pd.DataFrame,
    *,
    basin: Optional[int] = None,
    color_by: str = "k_sn",
    log_color: bool = True,
    cmap: str = "viridis",
    point_size: float = 4.0,
    show_fit: bool = True,
    observed: Optional[pd.DataFrame] = None,
    observed_color: str = "0.6",
    observed_size: float = 2.0,
    ax=None,
    colorbar: bool = True,
):
    """Plot a single basin's chi-elevation profile.

    Accepts a **nodes** table (observed points coloured by k_sn + fitted
    overlay) or a **segments** table (fitted segments coloured by their own
    k_sn). See the module docstring.

    Args:
        data: Nodes or segments table for one or more basins.
        basin: ``basin_index`` to plot. If ``None``, the first basin present.
        color_by: Colour column (default ``k_sn``).
        log_color: Colour on a log10 scale (default True; sensible for k_sn).
        cmap: matplotlib colourmap name.
        point_size: Observed-point marker size (nodes table only).
        show_fit: Draw the fitted profile (nodes table only; always drawn for
            a segments table).
        observed: Optional chi data-map frame (e.g. from
            :func:`lsdviztools3.io.chi.load_chi_points` on ``*_chi_data_maps.fgb``)
            whose observed chi-elevation points are overlaid as neutral markers.
            Useful with a segments table to show the data under the fit.
        observed_color: Marker colour for the observed overlay.
        observed_size: Marker size for the observed overlay.
        ax: Existing matplotlib ``Axes`` or ``None`` for a new one.
        colorbar: Attach a colourbar.

    Returns:
        ``(fig, ax)``.
    """
    plt, colors, LineCollection = _import_mpl()
    segments = _is_segments(data)

    if basin is None:
        basin = int(data["basin_index"].iloc[0])
    rows = data[data["basin_index"] == basin]
    if rows.empty:
        raise ValueError(f"No data found for basin_index={basin}.")

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = ax.figure

    if observed is not None:
        _scatter_observed(ax, observed, basin, 0.0,
                          color=observed_color, size=observed_size)

    mappable = None
    if segments:
        cvals = _colour_values(rows[color_by].to_numpy(), log_color)
        finite = cvals[np.isfinite(cvals)]
        norm = colors.Normalize(*(float(finite.min()), float(finite.max()))
                                if finite.size else (0.0, 1.0))
        lc = LineCollection(_segment_lines(rows, 0.0), cmap=cmap, norm=norm,
                            linewidths=1.4)
        lc.set_array(cvals)
        ax.add_collection(lc)
        ax.autoscale()
        mappable = lc
    else:
        cvals = _colour_values(rows[color_by].to_numpy(), log_color)
        mappable = ax.scatter(
            rows["chi"], rows["observed_elevation"],
            c=cvals, cmap=cmap, s=point_size, edgecolors="none", zorder=2,
        )
        if show_fit and "fitted_elevation" in rows.columns:
            for _, branch in rows.groupby("branch_id"):
                branch = branch.sort_values("chi")
                ax.plot(branch["chi"], branch["fitted_elevation"],
                        color="0.25", linewidth=0.8, zorder=1)

    ax.set_xlabel(r"$\chi$ (m)")
    ax.set_ylabel("Elevation (m)")
    ax.set_title(f"Basin {basin}")
    if colorbar and mappable is not None:
        fig.colorbar(mappable, ax=ax, label=_colour_label(color_by, log_color))
    return fig, ax


def plot_chi_stacked(
    data: pd.DataFrame,
    *,
    basins: Optional[Sequence[int]] = None,
    x_offset: Optional[float] = None,
    color_by: str = "k_sn",
    log_color: bool = True,
    cmap: str = "viridis",
    point_size: float = 3.0,
    show_fit: bool = True,
    observed: Optional[pd.DataFrame] = None,
    observed_color: str = "0.6",
    observed_size: float = 1.5,
    label_basins: bool = True,
    ax=None,
    colorbar: bool = True,
):
    """Stacked chi-elevation plot: basins normalised to chi=0 and offset right.

    For each basin the chi axis is shifted to start at zero and then translated
    by a cumulative offset so the profiles sit side-by-side on one axis, as in
    lsdviztools' ``StackedProfilesGradient``. Points/segments are coloured
    globally by ``color_by`` (typically log10 ``k_sn``).

    Accepts a **nodes** or **segments** table (module docstring).

    Args:
        data: Nodes or segments table.
        basins: Ordered ``basin_index`` values to stack. Default: all present.
        x_offset: Horizontal gap (chi units) between successive basins.
            Default: 20% of the widest single-basin chi span.
        color_by: Colour column (default ``k_sn``).
        log_color: Colour on a log10 scale (default True).
        cmap: matplotlib colourmap name.
        point_size: Marker size (nodes table only).
        show_fit: Overlay the fitted profile (nodes table only).
        observed: Optional chi data-map frame (e.g. from
            :func:`lsdviztools3.io.chi.load_chi_points`) whose observed points are
            overlaid as neutral markers, using the same per-basin offset.
        observed_color: Marker colour for the observed overlay.
        observed_size: Marker size for the observed overlay.
        label_basins: Annotate each stack with its ``basin_index``.
        ax: Existing matplotlib ``Axes`` or ``None``.
        colorbar: Attach a colourbar.

    Returns:
        ``(fig, ax)``.
    """
    plt, colors, LineCollection = _import_mpl()
    segments = _is_segments(data)

    if basins is None:
        basins = sorted(data["basin_index"].unique())
    if len(basins) == 0:
        raise ValueError("No basins to plot.")

    # Global colour normalisation across all selected basins.
    sel = data[data["basin_index"].isin(basins)]
    all_c = _colour_values(sel[color_by].to_numpy(), log_color)
    finite = all_c[np.isfinite(all_c)]
    vmin, vmax = (float(finite.min()), float(finite.max())) if finite.size else (0.0, 1.0)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    # Default offset: a fraction of the widest basin's chi span.
    if x_offset is None:
        widths = []
        for b in basins:
            rows = data[data["basin_index"] == b]
            if rows.empty:
                continue
            lo, hi = _basin_chi_bounds(rows, segments)
            widths.append(hi - lo)
        x_offset = 0.2 * max(widths) if widths else 0.0

    if ax is None:
        fig, ax = plt.subplots(figsize=(max(6, 1.6 * len(basins)), 4))
    else:
        fig = ax.figure

    mappable = None
    all_lines, all_cvals = [], []
    cursor = 0.0
    for b in basins:
        rows = data[data["basin_index"] == b]
        if rows.empty:
            continue
        chi0, chi_max = _basin_chi_bounds(rows, segments)
        width = chi_max - chi0
        shift = -chi0 + cursor

        if observed is not None:
            _scatter_observed(ax, observed, b, shift,
                              color=observed_color, size=observed_size)

        if segments:
            all_lines.extend(_segment_lines(rows, shift))
            all_cvals.append(_colour_values(rows[color_by].to_numpy(), log_color))
            zmin = float(rows["downstream_elevation"].min())
        else:
            cvals = _colour_values(rows[color_by].to_numpy(), log_color)
            mappable = ax.scatter(
                rows["chi"].to_numpy() + shift, rows["observed_elevation"],
                c=cvals, cmap=cmap, norm=norm, s=point_size,
                edgecolors="none", zorder=2,
            )
            if show_fit and "fitted_elevation" in rows.columns:
                for _, branch in rows.groupby("branch_id"):
                    branch = branch.sort_values("chi")
                    ax.plot(branch["chi"].to_numpy() + shift,
                            branch["fitted_elevation"],
                            color="0.25", linewidth=0.6, zorder=1)
            zmin = float(rows["observed_elevation"].min())

        if label_basins:
            ax.text(cursor + 0.05 * (width or 1.0), zmin, str(int(b)),
                    fontsize=8, va="bottom", ha="left")

        cursor += width + x_offset

    if segments and all_lines:
        lc = LineCollection(all_lines, cmap=cmap, norm=norm, linewidths=1.2)
        lc.set_array(np.concatenate(all_cvals))
        ax.add_collection(lc)
        ax.autoscale()
        mappable = lc

    ax.set_xlabel(r"$\chi$ (m, offset per basin)")
    ax.set_ylabel("Elevation (m)")
    if colorbar and mappable is not None:
        fig.colorbar(mappable, ax=ax, label=_colour_label(color_by, log_color))
    return fig, ax


def plot_ksn_vs_erate(
    data: pd.DataFrame,
    *,
    erate_col: str = "inferred_erate_mm_yr",
    ksn_col: str = "ksn_median",
    erate_err_low_col: str = "erate_err_low_mm_yr",
    erate_err_high_col: str = "erate_err_high_mm_yr",
    ksn_q25_col: str = "ksn_q25",
    ksn_q75_col: str = "ksn_q75",
    basin_label_col: Optional[str] = "basin_index",
    log_axes: bool = True,
    theta: Optional[float] = None,
    ax=None,
):
    """Plot basin-median k_sn versus inferred erosion rate with uncertainties.

    Y-errors use ``ksn_median`` against ``ksn_q25``/``ksn_q75`` (IQR). X-errors
    use the provided erosion-rate uncertainty columns.
    """
    plt, _, _ = _import_mpl()

    cols = [erate_col, ksn_col, ksn_q25_col, ksn_q75_col]
    for col in cols:
        if col not in data.columns:
            raise KeyError(f"Required column missing for k_sn vs erate plot: {col}")

    work = data.copy()
    for col in [erate_col, ksn_col, ksn_q25_col, ksn_q75_col, erate_err_low_col, erate_err_high_col]:
        if col in work.columns:
            work[col] = pd.to_numeric(work[col], errors="coerce")

    work = work[np.isfinite(work[erate_col]) & np.isfinite(work[ksn_col])]
    if work.empty:
        raise ValueError("No finite k_sn/erosion-rate rows available for plotting.")

    yerr_low = (work[ksn_col] - work[ksn_q25_col]).clip(lower=0).to_numpy(dtype=float)
    yerr_high = (work[ksn_q75_col] - work[ksn_col]).clip(lower=0).to_numpy(dtype=float)

    if erate_err_low_col in work.columns:
        xerr_low = work[erate_err_low_col].fillna(0).clip(lower=0).to_numpy(dtype=float)
    else:
        xerr_low = np.zeros(len(work), dtype=float)
    if erate_err_high_col in work.columns:
        xerr_high = work[erate_err_high_col].fillna(0).clip(lower=0).to_numpy(dtype=float)
    else:
        xerr_high = np.zeros(len(work), dtype=float)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.5, 4.0))
    else:
        fig = ax.figure

    ax.errorbar(
        work[erate_col].to_numpy(dtype=float),
        work[ksn_col].to_numpy(dtype=float),
        xerr=np.vstack([xerr_low, xerr_high]),
        yerr=np.vstack([yerr_low, yerr_high]),
        fmt="o",
        ms=4,
        color="0.2",
        ecolor="0.45",
        elinewidth=0.8,
        capsize=2,
    )

    if basin_label_col is not None and basin_label_col in work.columns:
        for _, row in work.iterrows():
            ax.text(row[erate_col], row[ksn_col], str(row[basin_label_col]), fontsize=7, ha="left", va="bottom")

    if log_axes:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel("Inferred erosion rate E (mm/yr)")
    ax.set_ylabel("Median ksn")
    if theta is None:
        ax.set_title("Basin median ksn vs inferred erosion rate")
    else:
        ax.set_title(rf"Basin median ksn vs inferred erosion rate (theta={theta:g})")
    return fig, ax
