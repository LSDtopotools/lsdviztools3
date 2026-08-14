"""Data layer for the synthetic-network snapping experiment.

The experiment asks two questions about landscapes generated as Optimal Channel
Networks (OCNs) and then *snapped* to a self-consistent fluvial + hillslope
form:

  (a) Does varying the OCN annealing temperature produce a real range of
      branching geometry?
  (b) Does snapping *reduce* that range -- i.e. in lithologically homogeneous
      material, is the achievable set of network geometries constrained by the
      requirement that gradients follow fluvial and hillslope form with no
      elevation discontinuity across divides?

Nothing here computes geomorphology. Flow directions and junction angles are
produced by lsdtt3 (``lsdtt-channel-extraction`` and
``lsdtt-channel-network-analysis``); this module loads those products and
compares them. The one derived quantity, the flow-direction change mask, is a
cell-by-cell comparison of two rasters lsdtt3 already wrote -- not a
reimplementation of D8 routing.

No PyGMT imports here, per the package's IO/render split.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import xarray as xr

from .raster import load_raster

# Suffix written by lsdtt-channel-extraction with
# write_arcmap_d8_flowdirection_raster: true
FLOWDIR_SUFFIX = "_arcmap_d8_flowdirection.tif"

# The experiment's variant codes, in the order they should appear in figures:
# the un-snapped surface first, then the three hillslope treatments.
VARIANTS = ("pre", "cs", "roe", "lin")

VARIANT_LABELS = {
    "pre": "no snapping",
    "cs": "critical slope",
    "roe": "Roering",
    "lin": "linear",
}

TEMPERATURE_LABELS = {
    "T00": "T = 0.0",
    "T05": "T = 0.5",
    "T10": "T = 1.0",
}


def dem_stem(tag: str, variant: str) -> str:
    """Filename stem lsdtt3 wrote for one (temperature, variant) cell.

    The un-snapped surface is the generator's own DEM; the snapped ones are the
    equilibrated DEM. Both come from the same invocation, so the pre-snap
    surface is shared by every variant at a given temperature.
    """
    if variant == "pre":
        return f"{tag}_pre_ocn_dem"
    return f"{tag}_{variant}_ocn_equilibrated_dem"


def load_flow_direction(path) -> xr.DataArray:
    """Load an lsdtt3 ArcMap-convention D8 flow-direction raster."""
    return load_raster(path)


def flow_direction_change_mask(
    pre_path,
    post_path,
    *,
    return_stats: bool = False,
):
    """Cells whose D8 flow direction differs between two surfaces.

    Both inputs must be flow-direction rasters lsdtt3 produced from the *same*
    grid with the *same* boundary conditions -- otherwise the comparison mixes a
    change in the landscape with a change in how it was routed.

    Cells that are nodata in either raster are excluded rather than counted as
    changed: a cell that gained or lost validity has not "changed direction" in
    any meaningful sense, and counting it would inflate the result.

    Args:
        pre_path: Flow direction before snapping.
        post_path: Flow direction after snapping.
        return_stats: Also return a dict of counts and the changed fraction.

    Returns:
        Boolean ``xarray.DataArray`` (True = direction changed), or
        ``(mask, stats)`` when ``return_stats`` is True.
    """
    pre = load_raster(pre_path)
    post = load_raster(post_path)

    if pre.shape != post.shape:
        raise ValueError(
            f"flow direction rasters differ in shape: {pre.shape} vs {post.shape}"
        )

    a = pre.values
    b = post.values

    valid = np.isfinite(a) & np.isfinite(b)
    # lsdtt3 writes ArcMap D8 codes (powers of two); 0 is also used for
    # unresolved/outlet cells in some paths, so treat it as not-a-direction.
    valid &= (a != 0) & (b != 0)

    changed = valid & (a != b)

    mask = xr.zeros_like(pre, dtype=bool)
    mask.values = changed
    mask = mask.where(valid)          # nodata outside the comparable set
    mask.name = "flow_direction_changed"

    if not return_stats:
        return mask

    n_valid = int(valid.sum())
    n_changed = int(changed.sum())
    stats = {
        "n_comparable": n_valid,
        "n_changed": n_changed,
        "fraction_changed": (n_changed / n_valid) if n_valid else float("nan"),
    }
    return mask, stats


def load_junction_angles(path) -> pd.DataFrame:
    """Load the junction-angle CSV from lsdtt-channel-network-analysis."""
    return pd.read_csv(path)


def _angle_column(df: pd.DataFrame) -> str:
    """Find the junction-angle column without hard-coding one spelling.

    The CSV schema has changed before; failing loudly with the available names
    beats silently picking a wrong column.
    """
    # donor_opening_angle_deg is the branching angle proper: the angle between
    # the two donor limbs at the junction. alpha_deg is its counterpart measured
    # against the receiver, and is kept as a fallback rather than a preference
    # because the two answer different questions.
    preferred = [
        "donor_opening_angle_deg",
        "junction_angle_deg",
        "junction_angle",
        "alpha_deg",
        "angle_deg",
        "angle",
    ]
    for name in preferred:
        if name in df.columns:
            return name
    candidates = [c for c in df.columns if "angle" in c.lower()]
    if len(candidates) == 1:
        return candidates[0]
    raise KeyError(
        "could not identify the junction-angle column; "
        f"looked for {preferred}, saw {list(df.columns)}"
    )


def summarise_junction_angles(df: pd.DataFrame) -> dict:
    """Summary statistics for one junction-angle table.

    The spread statistics are the point of the experiment: question (b) is
    whether snapping narrows the *distribution* of junction angles, so the IQR
    and standard deviation matter more than the mean.
    """
    col = _angle_column(df)
    a = pd.to_numeric(df[col], errors="coerce").dropna()
    if a.empty:
        return {
            "n_junctions": 0,
            "mean_deg": float("nan"),
            "median_deg": float("nan"),
            "std_deg": float("nan"),
            "iqr_deg": float("nan"),
            "p10_deg": float("nan"),
            "p90_deg": float("nan"),
        }
    q1, q3 = np.percentile(a, [25, 75])

    # D8 quantization diagnostic. Bearings are fitted over a sampled reach, so
    # the angles are mostly continuous -- but where a link is short and straight
    # the fit returns a bearing that is exactly a D8 direction, and the opening
    # angle lands exactly on a multiple of 45 degrees. Those spikes can carry
    # enough mass to pin the MEDIAN onto 45 or 90, which then looks like a
    # physical shift when it is really a grid artefact.
    #
    # Report the spike mass so any comparison of medians can be judged against
    # it rather than taken at face value.
    on_grid = np.isclose(np.mod(a.to_numpy(), 45.0), 0.0, atol=1e-6) | np.isclose(
        np.mod(a.to_numpy(), 45.0), 45.0, atol=1e-6
    )
    frac_on_d8 = float(on_grid.mean())

    off = a.to_numpy()[~on_grid]

    return {
        "n_junctions": int(a.size),
        "mean_deg": float(a.mean()),
        "median_deg": float(a.median()),
        "std_deg": float(a.std(ddof=1)) if a.size > 1 else float("nan"),
        "iqr_deg": float(q3 - q1),
        "p10_deg": float(np.percentile(a, 10)),
        "p90_deg": float(np.percentile(a, 90)),
        # Quantization diagnostics, and the same statistics with the exact-D8
        # angles removed, so a claim can be checked for robustness to them.
        "frac_on_d8_multiple": frac_on_d8,
        "n_off_grid": int(off.size),
        "median_deg_off_grid": float(np.median(off)) if off.size else float("nan"),
        "iqr_deg_off_grid": (
            float(np.subtract(*np.percentile(off, [75, 25]))) if off.size > 1
            else float("nan")
        ),
    }


def collect_experiment(
    directory,
    *,
    tags=("T00", "T05", "T10"),
    variants=VARIANTS,
) -> pd.DataFrame:
    """Assemble one tidy row per (temperature, variant) cell.

    Missing cells are reported with NaN rather than dropped, so a run that
    failed for one variant is visible in the table instead of quietly absent.
    """
    directory = Path(directory)
    rows = []

    for tag in tags:
        pre_stem = dem_stem(tag, "pre")
        pre_fd = directory / f"{pre_stem}{FLOWDIR_SUFFIX}"

        for variant in variants:
            stem = dem_stem(tag, variant)
            row = {
                "tag": tag,
                "temperature": TEMPERATURE_LABELS.get(tag, tag),
                "variant": variant,
                "variant_label": VARIANT_LABELS.get(variant, variant),
                "dem": str(directory / f"{stem}.tif"),
                "dem_exists": (directory / f"{stem}.tif").exists(),
            }

            ja_path = directory / f"{stem}_junction_angles.csv"
            if ja_path.exists():
                row.update(summarise_junction_angles(load_junction_angles(ja_path)))
                row["junction_angles_csv"] = str(ja_path)
            else:
                row["junction_angles_csv"] = None

            # Flow-direction change is only defined against the pre-snap surface.
            post_fd = directory / f"{stem}{FLOWDIR_SUFFIX}"
            if variant != "pre" and pre_fd.exists() and post_fd.exists():
                _, stats = flow_direction_change_mask(
                    pre_fd, post_fd, return_stats=True
                )
                row.update(stats)

            rows.append(row)

    return pd.DataFrame(rows)


def load_branching_metrics(directory, stem: str) -> dict:
    """Tokunaga fit and junction area-ratio summary for one DEM.

    Both are order-topological or ratio-based rather than bearing-based, so
    neither inherits the 45-degree grid lineation that contaminates junction
    angles on OCN surfaces -- which is why they are the right metrics here.

    The two are not interchangeable. The Tokunaga OLS fit is typically
    constrained by only a handful of order offsets on a domain this size, so its
    ``a`` can swing wildly; the area-ratio summary uses every interior
    confluence and carries bootstrap confidence intervals. Prefer
    ``c_equivalent`` and read ``tok_n_offsets_used`` before trusting ``tok_a``.

    Missing files give NaNs rather than raising, so one failed cell does not
    take down a sweep.
    """
    directory = Path(directory)
    out: dict = {}

    fit = directory / f"{stem}_tokunaga_summary_fit.csv"
    if fit.exists():
        df = pd.read_csv(fit)
        ols = df[df["method"] == "ols"] if "method" in df.columns else df
        if not ols.empty:
            r = ols.iloc[0]
            out.update(
                {
                    "tok_a": float(r.get("a", float("nan"))),
                    "tok_c": float(r.get("c", float("nan"))),
                    "tok_r2": float(r.get("r_squared", float("nan"))),
                    "tok_n_offsets_used": float(r.get("n_offsets_used", float("nan"))),
                }
            )

    ratio = directory / f"{stem}_junction_area_ratio_summary.csv"
    if ratio.exists():
        df = pd.read_csv(ratio)
        if not df.empty:
            r = df.iloc[0]
            for src, dst in (
                ("n_junctions", "ar_n_junctions"),
                ("median_ratio", "ar_median_ratio"),
                ("lambda_tail", "ar_lambda_tail"),
                ("r_a", "ar_r_a"),
                ("c_equivalent", "ar_c_equivalent"),
                ("c_equivalent_ci_low", "ar_c_ci_low"),
                ("c_equivalent_ci_high", "ar_c_ci_high"),
                ("n_branches", "ar_n_branches"),
            ):
                if src in df.columns:
                    out[dst] = float(r[src])

    return out


def collect_gamma_sweep(
    directory,
    *,
    gammas=("0.20", "0.35", "0.50", "0.65", "0.80"),
    seeds=(1, 2, 3),
    variants=(("pre", "no snapping"), ("roe", "Roering")),
) -> pd.DataFrame:
    """One tidy row per (gamma, seed, variant) of the gamma sweep."""
    directory = Path(directory)
    rows = []
    for g in gammas:
        gt = "g" + g.replace(".", "")
        for s in seeds:
            for code, label in variants:
                stem = (
                    f"{gt}_s{s}_pre_ocn_dem"
                    if code == "pre"
                    else f"{gt}_s{s}_{code}_ocn_equilibrated_dem"
                )
                row = {
                    "gamma": float(g),
                    "seed": s,
                    "variant": code,
                    "variant_label": label,
                    "stem": stem,
                    "exists": (directory / f"{stem}.tif").exists(),
                }
                row.update(load_branching_metrics(directory, stem))
                rows.append(row)
    return pd.DataFrame(rows)


def gamma_gate(table: pd.DataFrame, metric: str = "ar_c_equivalent") -> pd.DataFrame:
    """Does gamma separate the metric by more than the seed-to-seed spread?

    This is the gate, stated the way the OCN plan requires it: the signal is the
    spread of the PER-GAMMA MEANS across gamma, and the noise is the
    within-gamma spread across SEEDS -- never the within-network scatter, which
    is much smaller and would pass a gate that should fail.

    Returns one row per variant with the signal, the noise, and their ratio. A
    ratio comfortably above 1 means gamma is a usable lever for that variant.
    """
    out = []
    for (variant, label), grp in table.groupby(
        ["variant", "variant_label"], sort=False
    ):
        vals = grp.dropna(subset=[metric])
        if vals.empty:
            continue
        per_gamma = vals.groupby("gamma")[metric].mean()
        # Noise: mean within-gamma standard deviation across seeds.
        within = vals.groupby("gamma")[metric].std(ddof=1).mean()
        signal = float(per_gamma.max() - per_gamma.min())
        out.append(
            {
                "variant": variant,
                "variant_label": label,
                "metric": metric,
                "n_gamma": int(per_gamma.size),
                "signal_range_across_gamma": signal,
                "noise_mean_seed_sd": float(within) if pd.notna(within) else float("nan"),
                "signal_to_noise": (
                    signal / float(within)
                    if pd.notna(within) and within > 0
                    else float("nan")
                ),
                "min": float(per_gamma.min()),
                "max": float(per_gamma.max()),
            }
        )
    return pd.DataFrame(out)


def branching_range(table: pd.DataFrame, metric: str = "median_deg") -> pd.DataFrame:
    """Spread of a junction-angle metric ACROSS temperatures, per variant.

    This is the number question (b) turns on. For each snapping variant it
    measures how much the chosen metric still varies across the three annealing
    temperatures. A variant whose range collapses relative to ``pre`` is one
    where snapping has removed geometric freedom that the OCN annealing had put
    in.
    """
    out = []
    for variant, grp in table.groupby("variant", sort=False):
        vals = pd.to_numeric(grp[metric], errors="coerce").dropna()
        if vals.empty:
            continue
        out.append(
            {
                "variant": variant,
                "variant_label": VARIANT_LABELS.get(variant, variant),
                "n_temperatures": int(vals.size),
                "min": float(vals.min()),
                "max": float(vals.max()),
                "range": float(vals.max() - vals.min()),
                "std": float(vals.std(ddof=1)) if vals.size > 1 else float("nan"),
            }
        )
    df = pd.DataFrame(out)
    if not df.empty and "pre" in set(df["variant"]):
        base = float(df.loc[df["variant"] == "pre", "range"].iloc[0])
        df["range_relative_to_pre"] = df["range"] / base if base else float("nan")
    return df
