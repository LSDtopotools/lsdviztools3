"""Erosion rate against k_sn for lsdtt3's erodibility calibration (matplotlib).

:func:`plot_erosion_rate_vs_ksn` draws, log-log:

* the basins lsdtt3 used (``status == used``), coloured by their dominant
  lithologic unit with the shared name -> colour palette, with error bars from
  lsdtt3's ``k_sn_ci_low/high`` (95 % bootstrap) and
  ``erosion_rate_uncertainty_m_per_yr`` (1 sigma, as given);
* every other sample that has a k_sn as a hollow grey marker labelled with its
  status; samples without a k_sn are listed in the caption;
* the fit line(s) lsdtt3 wrote in the ``_fit_`` file: E = K A0^m k_sn^n
  (lsdtt3's formula, m = (m/n) n), evaluated with the written K, m and n across
  the plotted k_sn range — the fitted (York) n and each fixed n.

No confidence band is drawn: lsdtt3 writes separate n and K intervals, not
their joint distribution, so a band cannot be drawn from the files. The
intervals are printed as text instead. Every number on the figure comes from
lsdtt3's CSVs.
"""

from __future__ import annotations

from typing import Iterable, Mapping, Optional

import numpy as np
import pandas as pd

from ..io.erodibility_calibration import FIT_COLUMNS, SAMPLE_COLUMNS, USED_STATUS

_MPL_MSG = (
    "Calibration charts need matplotlib. Install the optional 'charts' extra, e.g. "
    "`pip install lsdviztools3[charts]` or `pip install matplotlib`."
)

#: Caption: why there is no band.
NO_BAND_NOTE = ("No confidence band: lsdtt3 writes separate n and K intervals, not "
                "their joint (n, K) distribution, so none is drawn.")

#: Line styles per fit row.
_FIT_STYLES = {"fitted_n": ("-", 2.0), 1.0: ("--", 1.3), 2.0: (":", 1.6)}


def _import_mpl():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_MPL_MSG) from exc
    return plt


def _num(v) -> Optional[float]:
    v = pd.to_numeric(pd.Series([v]), errors="coerce").iloc[0]
    return float(v) if np.isfinite(v) else None


def _g(v, fmt=".3g") -> str:
    v = _num(v)
    return "nan" if v is None else format(v, fmt)


def fit_line(row, k_sn) -> np.ndarray:
    """E = K A0^m k_sn^n for one lsdtt3 fit row (display of the written fit)."""
    c = FIT_COLUMNS
    K, n, m = _num(row[c["K"]]), _num(row[c["n"]]), _num(row.get(c["m"], np.nan))
    a0 = _num(row.get(c["A_0"], 1.0)) or 1.0
    if K is None or n is None:
        return np.full(np.shape(k_sn), np.nan)
    scale = a0 ** m if (m is not None and a0 != 1.0) else 1.0
    return K * scale * np.asarray(k_sn, dtype=float) ** n


def fit_label(row) -> str:
    """One-line description of a fit row from its own columns."""
    c = FIT_COLUMNS
    kind = str(row[c["fit_type"]])
    n = _g(row[c["n"]])
    head = f"fitted n = {n}" if kind == "fitted_n" else f"fixed n = {n}"
    if kind == "fitted_n":
        head += f" [{_g(row.get(c['n_ci_low']))}, {_g(row.get(c['n_ci_high']))}]"
    head += (f", K = {_g(row[c['K']])} [{_g(row.get(c['K_ci_low']))}, "
             f"{_g(row.get(c['K_ci_high']))}] {row.get(c['K_units'], '')}")
    if _num(row.get(c["mswd"])) is not None:
        head += f", MSWD {_g(row.get(c['mswd']))}"
    return head


def plot_erosion_rate_vs_ksn(
    samples: pd.DataFrame,
    fit: Optional[pd.DataFrame] = None,
    *,
    unit_names: Optional[Iterable[str]] = None,
    unit_colours: Optional[Mapping[str, str]] = None,
    show_fixed_n: Iterable[float] = (1.0, 2.0),
    title: Optional[str] = None,
    ax=None,
):
    """Log-log E against k_sn with lsdtt3's fit line(s); see the module docstring.

    Args:
        samples: ``*_erodibility_calibration_samples_mn_*.csv`` table.
        fit: ``*_erodibility_calibration_fit_mn_*.csv`` table (optional).
        unit_names: All unit names of the run (e.g. the run's
            ``_lithology_units.csv``), so colours match the run's other
            figures. Defaults to the samples' dominant units.
        unit_colours: ``{name: colour}`` overrides.
        show_fixed_n: Fixed-n rows to overlay (if lsdtt3 wrote them with
            status ``ok``).
        title: Title prefix.
        ax: Existing axes or ``None``.

    Returns:
        ``(fig, ax)``.
    """
    from ..render.lithology_colours import lithology_palette
    from .lithology_ksn import ksn_label

    plt = _import_mpl()
    c = SAMPLE_COLUMNS
    s = samples.copy()
    st = s[c["status"]].astype(str)
    k = pd.to_numeric(s[c["k_sn"]], errors="coerce")
    e = pd.to_numeric(s[c["E"]], errors="coerce")
    used = (st == USED_STATUS) & (k > 0) & (e > 0)
    other = ~used & (k > 0) & (e > 0)
    unplotted = s[~used & ~other]

    names = s[c["dominant_lith_name"]].dropna().astype(str)
    pal = lithology_palette(list(unit_names) if unit_names is not None else names,
                            unit_colours)
    for n in names:  # a unit missing from unit_names still gets a colour
        if n not in pal:
            pal = lithology_palette(list(pal) + [n], unit_colours)

    if ax is None:
        fig, ax = plt.subplots(figsize=(11.5, 6.2))
    else:
        fig = ax.figure

    u = s[used]
    lo = pd.to_numeric(u[c["k_sn_ci_low"]], errors="coerce")
    hi = pd.to_numeric(u[c["k_sn_ci_high"]], errors="coerce")
    sig = pd.to_numeric(u[c["E_sigma"]], errors="coerce")
    seen = set()
    clipped = 0
    for (_, r), kl, kh, sg in zip(u.iterrows(), lo, hi, sig):
        name = r.get(c["dominant_lith_name"])
        name = "no dominant unit" if pd.isna(name) else str(name)
        col = pal.get(name, "#555557")
        kk, ee = float(r[c["k_sn"]]), float(r[c["E"]])
        xerr = [[kk - kl], [kh - kk]] if np.isfinite(kl) and np.isfinite(kh) else None
        yerr = None
        if np.isfinite(sg):
            # A lower bar reaching E <= 0 cannot be drawn on a log axis.
            clipped += int(sg >= ee)
            yerr = [[min(sg, ee * 0.999)], [sg]]
        ax.errorbar(kk, ee, xerr=xerr, yerr=yerr, fmt="o", color=col, ecolor=col,
                    elinewidth=0.9, capsize=2, markersize=6, markeredgecolor="black",
                    markeredgewidth=0.5, label=None if name in seen else name, zorder=3)
        seen.add(name)

    o = s[other]
    if len(o):
        ax.scatter(pd.to_numeric(o[c["k_sn"]]), pd.to_numeric(o[c["E"]]), s=36,
                   facecolors="none", edgecolors="0.45", linewidths=1.0, zorder=2,
                   label="not used (status as labelled)")
        for _, r in o.iterrows():
            ax.annotate(str(r[c["status"]]), (float(r[c["k_sn"]]), float(r[c["E"]])),
                        xytext=(4, -3), textcoords="offset points", fontsize=6, color="0.35")

    lines = []
    if fit is not None and len(u):
        f = FIT_COLUMNS
        kgrid = np.geomspace(float(k[used | other].min()), float(k[used | other].max()), 200)
        for _, r in fit.iterrows():
            if str(r[f["status"]]) != "ok":
                lines.append(f"{r[f['fit_type']]} (n={_g(r[f['n']])}): status "
                             f"{r[f['status']]}, not drawn")
                continue
            kind = str(r[f["fit_type"]])
            nval = _num(r[f["n"]])
            if kind == "fitted_n":
                ls, lw = _FIT_STYLES["fitted_n"]
            elif nval is not None and any(abs(nval - x) < 1e-9 for x in show_fixed_n):
                ls, lw = _FIT_STYLES.get(nval, ("-.", 1.2))
            else:
                continue
            ax.plot(kgrid, fit_line(r, kgrid), ls, color="black", linewidth=lw,
                    label=fit_label(r), zorder=4)
            if kind == "fitted_n":
                lines.append(
                    f"fitted n: N_basins {_g(r[f['N_basins']], '.0f')}, N_groups "
                    f"{_g(r[f['N_groups']], '.0f')}, corr(n, ln K) {_g(r.get(f['corr_n_lnK']))}, "
                    f"{_g(r.get(f['n_bootstrap']), '.0f')} bootstrap replicates over "
                    f"{r.get(f['bootstrap_cluster'], '')}")
        r0 = fit.iloc[0]
        lines.insert(0, f"m/n = {_g(r0.get(f['m_over_n']))}, A0 = {_g(r0.get(f['A_0']))} m^2; "
                        "lines: E = K A0^m k_sn^n (lsdtt3), m = (m/n) n")

    ax.set_xscale("log")
    ax.set_yscale("log")
    mn = _num(s[c["m_over_n"]].iloc[0]) if c["m_over_n"] in s.columns and len(s) else None
    ax.set_xlabel("basin " + ksn_label(mn) + "  (bars: 95 % bootstrap, lsdtt3)")
    ax.set_ylabel("erosion rate E (m yr$^{-1}$)  (bars: 1$\\sigma$ as given)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=6.5, loc="upper left", bbox_to_anchor=(1.02, 1.0), framealpha=0.9,
              borderaxespad=0.0)

    notes = list(lines) + [NO_BAND_NOTE]
    if clipped:
        notes.append(f"{clipped} lower E bar(s) reach E <= 0 and are drawn cut at the "
                     "axis (log scale); the sigma values are as written.")
    if len(unplotted):
        notes.append("Not plotted (no k_sn): " + ", ".join(
            f"{r[c['sample_id']]} ({r[c['status']]})" for _, r in unplotted.iterrows()))
    import textwrap

    ax.set_title((title + "\n" if title else "") + "\n".join(
        textwrap.fill(n, 120) for n in notes), fontsize=6.5, loc="left")
    fig.tight_layout()
    return fig, ax
