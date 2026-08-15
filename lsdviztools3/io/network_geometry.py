"""Readers for the network-geometry sweep outputs.

The sweep compares drainage networks produced by different mechanisms -- an
Optimal Channel Network annealed across temperature, against stochastic
surfaces (diamond-square, red noise, white noise) whose networks are simply
whatever D8 finds. The question is which measurable properties actually
separate an optimised network from an arbitrary one.

This module reads the lsdtt3 outputs; it computes nothing. Two metrics here are
newer than the ones in :mod:`lsdviztools3.io.snapping` and are the reason this
module exists:

``beta``
    Exponent of the drainage-area exceedance distribution, P(A > a) ~ a^-beta,
    from ``lsdtt-channel-network-analysis``. Fitted two-sided by MLE, because
    the finite grid rolls the distribution over at the top and a tail-only OLS
    fit reads the rollover rather than the scaling.

``delta_L``
    Lengthwise asymmetry of paired flow paths (Goren & Shelef, 2024), from
    ``lsdtt-ridgeline-analysis``. For the two flow paths that diverge across a
    divide and rejoin at a common junction,
    ``delta_L = 2|L_i - L_j| / (L_i + L_j)``. Published medians run 0.19 for
    sub-parallel trunk streams to 1.00 for tortuous networks.

Every loader returns NaNs for a missing file rather than raising, so one failed
cell of a sweep does not take the whole table down. The counterpart of that
choice is that you must check ``n_pairs`` / ``n_used`` before trusting a
summary statistic -- see :func:`load_paired_flowpaths`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# Individual metric files
# --------------------------------------------------------------------------

def load_area_scaling(directory, stem: str) -> dict:
    """Drainage-area exceedance exponent beta for one DEM.

    Reads ``<stem>_area_scaling_summary.csv``.

    ``n_decades`` is the one to check before quoting ``beta``: a power law
    fitted over less than about 1.5 decades is not constrained, and lsdtt3
    reports the fitted range precisely so this can be judged rather than
    assumed.
    """
    path = Path(directory) / f"{stem}_area_scaling_summary.csv"
    out = {
        "beta": float("nan"),
        "beta_stderr": float("nan"),
        "beta_ols": float("nan"),
        "beta_r2": float("nan"),
        "beta_n_used": float("nan"),
        "beta_n_decades": float("nan"),
    }
    if not path.exists():
        return out

    df = pd.read_csv(path)
    if df.empty:
        return out

    r = df.iloc[0]
    for src, dst in (
        ("beta", "beta"),
        ("beta_stderr", "beta_stderr"),
        ("beta_ols", "beta_ols"),
        ("r_squared", "beta_r2"),
        ("n_used", "beta_n_used"),
        ("n_decades", "beta_n_decades"),
    ):
        if src in df.columns:
            out[dst] = float(r[src])
    return out


def load_area_scaling_curve(directory, stem: str) -> Optional[pd.DataFrame]:
    """The exceedance curve itself, for plotting. ``None`` if absent."""
    path = Path(directory) / f"{stem}_area_scaling_curve.csv"
    return pd.read_csv(path) if path.exists() else None


def _paired_flowpath_path(directory, stem: str, suffix: str) -> Path:
    """Locate a paired-flowpath product for ``stem``.

    The ΔL step runs `lsdtt-ridgeline-analysis` under a ``pf_`` write prefix, so
    its many other outputs (READMEs, ridgeline CSVs) cannot collide with the
    products the main sweep already wrote for the same DEM. That leaves the
    paired-flowpath files named ``pf_<stem>_...`` rather than ``<stem>_...``.

    Both spellings are accepted, unprefixed first. Doing this here rather than
    copying files into place matters for reproducibility: a copy goes stale
    silently when the ΔL step is rerun, and the plot would then show the old
    numbers with no indication anything was wrong.
    """
    directory = Path(directory)
    direct = directory / f"{stem}{suffix}"
    if direct.exists():
        return direct
    return directory / f"pf_{stem}{suffix}"


def load_paired_flowpaths(directory, stem: str) -> dict:
    """Lengthwise-asymmetry summary for one DEM.

    Reads ``<stem>_paired_flowpaths_summary.csv``.

    The rejection counts are returned alongside the medians on purpose. A
    paired-flow-path run can retain very few pairs -- divide edges are rejected
    when the two paths never rejoin, which happens whenever the flanking basins
    sit in different parents -- and a median over a few hundred pairs is not
    comparable with one over tens of thousands. Read ``n_pairs`` before reading
    ``median_delta_L``; the counts are how you tell "no complexity" from "no
    surviving pairs".
    """
    path = _paired_flowpath_path(directory, stem, "_paired_flowpaths_summary.csv")
    out = {
        "median_delta_L": float("nan"),
        "mean_delta_L": float("nan"),
        "median_delta_chi": float("nan"),
        "n_pairs": float("nan"),
        "n_edges_examined": float("nan"),
        "n_rejected_no_junction": float("nan"),
        "n_rejected_distance": float("nan"),
    }
    if not path.exists():
        return out

    df = pd.read_csv(path)
    if df.empty:
        return out

    r = df.iloc[0]
    for key in out:
        if key in df.columns:
            out[key] = float(r[key])

    # Fraction of examined divide edges whose two flanking paths actually
    # rejoin. This starts life as a sampling diagnostic -- a low value means the
    # median delta_L rests on few pairs -- but it is also a structural property
    # in its own right: it is the fraction of the divide network that is
    # INTERNAL to a shared parent basin, as opposed to separating basins that
    # never meet. A landscape of long parallel strips draining to opposite
    # boundaries has almost no internal divides; a deeply nested, branching one
    # has many. Reported separately from delta_L because it survives even when
    # delta_L cannot be computed at all.
    edges = out.get("n_edges_examined", float("nan"))
    pairs = out.get("n_pairs", float("nan"))
    out["internal_divide_fraction"] = (
        pairs / edges if np.isfinite(edges) and edges > 0 else float("nan")
    )
    return out


def load_paired_flowpath_pairs(directory, stem: str) -> Optional[pd.DataFrame]:
    """Per-pair delta_L table, for distribution plots. ``None`` if absent."""
    path = _paired_flowpath_path(directory, stem, "_paired_flowpaths.csv")
    return pd.read_csv(path) if path.exists() else None


def load_edge_basin_sizes(directory, stem: str, prefix: str = "eb_"):
    """Sizes of the edge-draining (top-level) basins, as domain fractions.

    Reads ``<prefix><stem>_selected_basin_outlets.csv`` from
    ``lsdtt-basin-extraction`` run with ``select_all_basins`` and
    ``filter_nested_basins``, which selects exactly the non-nested basins. With
    ``outlet_mode: north_south`` every basin drains to a boundary, so those
    basins partition the tile and their fractions sum to ~1 — worth checking,
    because a sum well below 1 means the selection was not a partition and the
    shares are not interpretable as such.

    Returns the ``basin_fraction_of_domain`` array, or ``None`` if absent.
    """
    path = Path(directory) / f"{prefix}{stem}_selected_basin_outlets.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    if "basin_fraction_of_domain" not in df.columns or df.empty:
        return None
    values = pd.to_numeric(
        df["basin_fraction_of_domain"], errors="coerce").to_numpy(float)
    values = values[np.isfinite(values)]
    return values if values.size else None


def summarise_edge_basins(fractions) -> dict:
    """Summary of an edge-draining basin size distribution.

    The plain median is close to useless here and is reported only so that is
    visible: most top-level basins are single boundary cells with negligible
    area, so the median sits near zero for every landscape regardless of how
    the tile is really divided up.

    ``area_weighted_median`` is the statistic that answers "do a few basins take
    most of the tile?" — it is the basin size such that half the total drained
    area lies in basins at least that big. Together with ``max_fraction`` and
    ``top5_share`` it distinguishes a landscape of many equal strips from one
    dominated by a couple of large catchments.
    """
    f = np.asarray(fractions, dtype=float)
    f = f[np.isfinite(f)]
    if f.size == 0:
        nan = float("nan")
        return {"n_edge_basins": 0, "median_fraction": nan,
                "area_weighted_median_fraction": nan, "max_fraction": nan,
                "top5_share": nan, "fraction_sum": nan}

    desc = np.sort(f)[::-1]
    total = desc.sum()
    cumulative = np.cumsum(desc) / total if total > 0 else np.zeros_like(desc)
    idx = int(np.searchsorted(cumulative, 0.5))
    idx = min(idx, desc.size - 1)

    return {
        "n_edge_basins": int(f.size),
        "median_fraction": float(np.median(f)),
        "area_weighted_median_fraction": float(desc[idx]),
        "max_fraction": float(desc[0]),
        "top5_share": float(desc[:5].sum()),
        "fraction_sum": float(total),
    }


def load_hack(directory, stem: str) -> dict:
    """Hack exponent from the channel-extraction summary, if present."""
    path = Path(directory) / f"{stem}_hack_summary.csv"
    out = {"hack_h": float("nan"), "hack_r2": float("nan")}
    if not path.exists():
        return out
    df = pd.read_csv(path)
    if df.empty:
        return out
    r = df.iloc[0]
    if "h" in df.columns:
        out["hack_h"] = float(r["h"])
    if "r_squared" in df.columns:
        out["hack_r2"] = float(r["r_squared"])
    return out


# --------------------------------------------------------------------------
# Family labelling
# --------------------------------------------------------------------------

#: Ordered so that plots put the optimised networks first and the
#: zero-correlation end member last.
FAMILY_ORDER = (
    "Optimal Channel Network",
    "Diamond-square fractal surface",
    "Red-noise surface",
    "White-noise surface",
)

FAMILY_COLOURS = {
    "Optimal Channel Network": "#1f77b4",
    "Diamond-square fractal surface": "#d62728",
    "Red-noise surface": "#2ca02c",
    "White-noise surface": "#7f7f7f",
}


def family_of(stem: str) -> str:
    """Generation family from a DEM stem.

    Deliberately explicit rather than clever: the stems are set by the sweep
    script, and an unrecognised one should be visible in the legend as
    ``unknown`` rather than silently folded into a neighbouring family.
    """
    if stem.startswith("ocn"):
        return "Optimal Channel Network"
    if stem.startswith("ds"):
        return "Diamond-square fractal surface"
    if stem.startswith("rn"):
        return "Red-noise surface"
    if stem.startswith("wn"):
        return "White-noise surface"
    return "unknown"


def temperature_of(stem: str, directory=None) -> float:
    """Annealing temperature for an OCN stem, else NaN.

    THE TAG ALONE IS AMBIGUOUS AND MUST NOT BE TRUSTED. The sweep names cases by
    writing the temperature with the decimal point deleted, which is not
    invertible: ``ocn_T025`` is 0.25 if it came from "0.25" but 2.5 if it came
    from "02.5", and ``ocn_T20`` is 2.0 rather than the 20.0 a naive reading
    gives. That misreading silently rescales the x-axis of every
    temperature plot.

    So when ``directory`` is given, the value is read from the generating
    parameter file (``<tag>.param``, key ``temperature_scale``), which is
    authoritative. The tag heuristic is used only as a fallback, and only for
    the ``0.xx`` forms where it happens to be unambiguous.
    """
    if not stem.startswith("ocn_T"):
        return float("nan")

    tag = "ocn_T" + stem[len("ocn_T"):].split("_")[0]

    if directory is not None:
        path = Path(directory) / f"{tag}.param"
        if path.exists():
            for line in path.read_text().splitlines():
                if line.strip().startswith("temperature_scale:"):
                    try:
                        return float(line.split(":", 1)[1].strip())
                    except ValueError:
                        break

    digits = tag[len("ocn_T"):]
    if not digits.isdigit():
        return float("nan")
    # Fallback: only the leading-zero forms are unambiguous ("00" -> 0.0,
    # "001" -> 0.01, "025" -> 0.25, "05" -> 0.5). Anything else needs the
    # parameter file, so say so with NaN rather than guess wrong.
    if digits.startswith("0"):
        return float("0." + digits[1:]) if len(digits) > 1 else 0.0
    return float("nan")


# --------------------------------------------------------------------------
# Sweep collection
# --------------------------------------------------------------------------

def collect_network_sweep(
    directory,
    stems: Sequence[str],
    *,
    include_junction_angles: bool = True,
) -> pd.DataFrame:
    """One tidy row per DEM, combining every metric the sweep produces.

    Args:
        directory: Directory holding the lsdtt3 outputs.
        stems: DEM stems, e.g. ``["ocn_T00_ocn_dem", "wn_uniform_random"]``.
        include_junction_angles: Also summarise the junction-angle CSV. This is
            the slowest part, since it reads a row per junction.

    Returns:
        DataFrame with one row per stem, carrying ``stem``, ``family``,
        ``temperature`` and every metric column. Missing metrics are NaN.
    """
    from .snapping import (
        load_branching_metrics,
        load_junction_angles,
        summarise_junction_angles,
    )

    directory = Path(directory)
    rows = []

    for stem in stems:
        row: dict = {
            "stem": stem,
            "family": family_of(stem),
            "temperature": temperature_of(stem, directory),
        }
        row.update(load_branching_metrics(directory, stem))
        row.update(load_area_scaling(directory, stem))
        row.update(load_paired_flowpaths(directory, stem))
        row.update(load_hack(directory, stem))

        edge = load_edge_basin_sizes(directory, stem)
        if edge is not None:
            row.update(summarise_edge_basins(edge))

        if include_junction_angles:
            ja = directory / f"{stem}_junction_angles.csv"
            if ja.exists():
                try:
                    row.update(summarise_junction_angles(load_junction_angles(ja)))
                except Exception:
                    # A malformed angle table should cost this row its angle
                    # columns, not the whole sweep.
                    pass

        rows.append(row)

    table = pd.DataFrame(rows)
    if "family" in table.columns:
        table["family"] = pd.Categorical(
            table["family"],
            categories=list(FAMILY_ORDER) + ["unknown"],
            ordered=True,
        )
    return table


def load_junction_angle_sets(
    directory,
    stems: Sequence[str],
    *,
    angle_column: Optional[str] = None,
) -> "dict[str, np.ndarray]":
    """Junction angles per stem, as arrays ready for a CDF.

    Returns a mapping ``stem -> angles``. Stems whose CSV is missing or empty
    are omitted entirely rather than mapped to an empty array, so a caller
    plotting the result cannot draw an invisible line and think it plotted
    something.
    """
    from .snapping import _angle_column, load_junction_angles

    directory = Path(directory)
    out: "dict[str, np.ndarray]" = {}

    for stem in stems:
        path = directory / f"{stem}_junction_angles.csv"
        if not path.exists():
            continue
        df = load_junction_angles(path)
        if df.empty:
            continue
        col = angle_column or _angle_column(df)
        values = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if values.size:
            out[stem] = values

    return out
