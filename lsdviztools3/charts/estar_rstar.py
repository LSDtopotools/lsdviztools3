"""E*-R* charts: measured hilltop metrics against nonlinear hillslope theory.

For a one-dimensional steady hillslope with the Roering nonlinear flux law
``q = D S / (1 - (S/S_c)^2)``, the non-dimensional erosion rate and relief are
related by a single curve:

    R* = (1/E*) * [ sqrt(1 + E*^2) - ln( (1 + sqrt(1 + E*^2)) / 2 ) - 1 ]

Roering, J. J., Perron, J. T. and Kirchner, J. W. (2007). Functional
relationships between denudation, hillslope form and relief. *Earth and
Planetary Science Letters* 264:245-258.
Hurst, M. D., Mudd, S. M., Walcott, R., Attal, M. and Yoo, K. (2012). Using
hilltop curvature to derive the spatial distribution of erosion rates. *JGR
Earth Surface* 117:F02017.
Grieve, S. W. D., Mudd, S. M., Hurst, M. D. and Milodowski, D. T. (2016). A
nondimensional framework for exploring the relief structure of landscapes.
*Earth Surface Dynamics* 4:309-325.

The curve is 1-D theory that knows nothing about how E* and R* were measured, so
plotting measured values against it tests the tracer, the ridge extraction, the
curvature convention and the physics at once. Two limits are worth checking by
eye because they are easy to get wrong:

* ``E* -> 0``   gives ``R* -> E*/4``  (linear/diffusive limit)
* ``E* -> inf`` gives ``R* -> 1``     (threshold limit; relief saturates at S_c*Lh)

**R* > 1 is not possible.** Systematic exceedance means Lh, R or S_c is wrong,
or the landscape is not at hillslope steady state at all -- it is the first
thing to look at, before any curve fitting.

Nothing here computes a landscape quantity: E*, R* and the trace status all come
from lsdtt3's own output.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

_MPL_MSG = ("The E*-R* charts need matplotlib. "
            "Install it e.g. `pip install matplotlib`.")


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def rstar_theory(e_star: np.ndarray | float) -> np.ndarray:
    """The Roering et al. (2007) steady nonlinear-hillslope relation.

    Args:
        e_star: Non-dimensional erosion rate, ``2 |Cht| Lh / S_c``. Must be > 0;
            zero and negative values return NaN rather than raising, so a whole
            column can be passed without pre-filtering.

    Returns:
        ``R*`` for each ``E*``.
    """
    e = np.asarray(e_star, dtype=float)
    out = np.full(e.shape, np.nan)
    good = e > 0
    eg = e[good]
    s = np.sqrt(1.0 + eg * eg)
    out[good] = (s - np.log((1.0 + s) / 2.0) - 1.0) / eg
    return out


def plot_estar_rstar_model_scan(data, *, figsize=(11, 5), title=None):
    """Compare analytical, Newton and cellular strip results in E*-R* space."""
    plt = _import_mpl()
    frame = data.copy() if isinstance(data, pd.DataFrame) else pd.read_csv(data)
    required = {"case", "model", "E_star", "R_star"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise KeyError(f"Missing E*-R* scan columns: {missing}")
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    e_min = max(0.05, float(frame.E_star.min()) * 0.8)
    e_max = float(frame.E_star.max()) * 1.25
    theory_e = np.logspace(np.log10(e_min), np.log10(e_max), 400)
    axes[0].plot(theory_e, rstar_theory(theory_e), "k-", lw=2,
                 label="analytical Roering")
    styles = {"analytical": ("o", "black"), "newton": ("s", "#0072B2"),
              "cellular": ("^", "#D55E00")}
    for model in ("analytical", "newton", "cellular"):
        subset = frame[frame.model == model]
        marker, color = styles[model]
        axes[0].scatter(subset.E_star, subset.R_star, marker=marker, color=color,
                        s=45, label=model.capitalize(), zorder=4)
        residual = subset.R_star.to_numpy() - rstar_theory(subset.E_star.to_numpy())
        axes[1].scatter(subset.E_star, residual, marker=marker, color=color,
                        s=45, label=model.capitalize())
    for case, subset in frame.groupby("case"):
        axes[0].plot(subset.E_star, subset.R_star, color="0.75", lw=0.8, zorder=1)
        row = subset[subset.model == "analytical"].iloc[0]
        axes[0].annotate(str(case), (row.E_star, row.R_star), xytext=(3, 3),
                         textcoords="offset points", fontsize=7)
    for ax in axes:
        ax.set_xscale("log"); ax.set_xlabel("E* = 2 |Laplacian at divide| Lh / Sc")
        ax.grid(alpha=0.25); ax.legend(fontsize=8)
    axes[0].set_ylabel("R* = relief / (Lh Sc)")
    axes[0].set_title("Dimensionless erosion-relief relation")
    axes[1].axhline(0, color="k", lw=1)
    axes[1].set_ylabel("R* residual from analytical Roering")
    axes[1].set_title("Numerical model residual")
    fig.suptitle(title or "Roering strip scan across D, Sc and Lh")
    fig.tight_layout()
    return fig, axes


def usable_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Rows whose E* and R* are meaningful.

    Only traces that reached a channel have a defined hillslope length, and E*
    is undefined where the curvature raster was nodata. Everything else must be
    excluded from the comparison and counted, not quietly dropped.
    """
    if df.empty:
        return df
    ok = (df["trace_status"] == 0) & df["E_star"].notna() & df["R_star"].notna()
    return df[ok]


def plot_estar_rstar(
    df: pd.DataFrame,
    *,
    group_col: Optional[str] = None,
    title: Optional[str] = None,
    e_range: Sequence[float] = (0.05, 200.0),
    figsize=(11, 5),
):
    """Measured E*-R* against theory, plus the residual in R*.

    Args:
        df: Pooled hilltop metrics. Needs ``E_star``, ``R_star``,
            ``trace_status``; ``group_col`` if supplied colours the points.
        group_col: Column to colour by, e.g. the ``U/D`` label of each landscape.
            Landscapes with the same U/D should land in the same place, which is
            itself a prediction worth being able to see.
        title: Figure title.
        e_range: E* range for the theoretical curve.
        figsize: Figure size in inches.

    Returns:
        The figure.
    """
    plt = _import_mpl()

    good = usable_rows(df)
    if good.empty:
        raise ValueError("plot_estar_rstar: no rows with a defined E* and R*")

    e = good["E_star"].to_numpy()
    r = good["R_star"].to_numpy()
    resid = r - rstar_theory(e)

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # ---- panel 1: the cloud against the curve
    ax = axes[0]
    ec = np.logspace(np.log10(e_range[0]), np.log10(e_range[1]), 400)
    ax.plot(ec, rstar_theory(ec), "k-", lw=2, zorder=5,
            label="Roering et al. (2007) theory")
    ax.plot(ec, ec / 4.0, "k--", lw=1, alpha=0.6, zorder=4,
            label="linear limit  R* = E*/4")
    ax.axhline(1.0, color="firebrick", ls=":", lw=1.5, zorder=4,
               label="R* = 1 (hard bound)")

    if group_col is not None and group_col in good.columns:
        cmap = plt.get_cmap("viridis")
        keys = sorted(good[group_col].unique())
        for i, key in enumerate(keys):
            sub = good[good[group_col] == key]
            ax.scatter(sub["E_star"], sub["R_star"], s=5, alpha=0.45,
                       color=cmap(i / max(len(keys) - 1, 1)), label=str(key))
    else:
        ax.scatter(e, r, s=5, alpha=0.35, color="#1f77b4", label="measured")

    ax.set_xscale("log")
    ax.set_xlabel("E*  =  2 |Cht| Lh / S_c")
    ax.set_ylabel("R*  =  R / (Lh S_c)")
    ax.set_ylim(0, max(1.4, float(np.nanpercentile(r, 99)) * 1.1))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="lower right")
    ax.set_title("Measured hilltop metrics against nonlinear theory")

    # ---- panel 2: the residual, which says WHICH term is wrong
    ax = axes[1]
    ax.axhline(0.0, color="k", lw=1.5, zorder=5)
    ax.scatter(e, resid, s=5, alpha=0.35, color="#8c564b")
    ax.set_xscale("log")
    ax.set_xlabel("E*")
    ax.set_ylabel("R* measured  -  R* theory")
    ax.grid(alpha=0.25)
    ax.set_title(
        f"Residual: median {np.nanmedian(resid):+.3f}, "
        f"IQR {np.nanpercentile(resid, 75) - np.nanpercentile(resid, 25):.3f}\n"
        "a vertical offset implicates Cht or S_c; a tilt implicates Lh",
        fontsize=9)

    n_all = int(len(df))
    n_used = int(len(good))
    n_over = int((r > 1.0).sum())
    fig.suptitle(
        (title or "E*-R* validation") +
        f"   |   {n_used} of {n_all} rows usable, "
        f"{n_over} with R* > 1 ({100.0 * n_over / max(n_used, 1):.1f} %)",
        fontsize=11, y=1.02)

    fig.tight_layout()
    return fig


def summarise(df: pd.DataFrame, group_col: Optional[str] = None) -> pd.DataFrame:
    """Per-group residual statistics against the theoretical curve.

    Returns:
        One row per group (or one overall) with the counts, the median and IQR
        of the R* residual, and the fraction violating ``R* < 1``. This is the
        table the acceptance criteria in the plan are read off.
    """
    def _stats(sub: pd.DataFrame, label) -> dict:
        good = usable_rows(sub)
        if good.empty:
            return {"group": label, "n_rows": len(sub), "n_usable": 0}
        e = good["E_star"].to_numpy()
        r = good["R_star"].to_numpy()
        resid = r - rstar_theory(e)
        return {
            "group": label,
            "n_rows": len(sub),
            "n_usable": len(good),
            "frac_usable": len(good) / max(len(sub), 1),
            "E_star_median": float(np.nanmedian(e)),
            "R_star_median": float(np.nanmedian(r)),
            "resid_median": float(np.nanmedian(resid)),
            "resid_iqr": float(np.nanpercentile(resid, 75)
                               - np.nanpercentile(resid, 25)),
            "frac_Rstar_gt_1": float((r > 1.0).mean()),
        }

    if group_col is None or group_col not in df.columns:
        return pd.DataFrame([_stats(df, "all")])
    return pd.DataFrame([_stats(df[df[group_col] == k], k)
                         for k in sorted(df[group_col].unique())])
