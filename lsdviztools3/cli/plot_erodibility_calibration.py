"""CLI: lsdviz3-erodibility-calibration — plots of lsdtt3's K calibration.

Reads the erodibility-calibration outputs of ``lsdtt-chi-analysis``
(``<prefix>_erodibility_calibration_{samples,fit}_mn_<mn>.csv`` and
``<prefix>_erodibility_calibration_basins.fgb``) and writes:

* ``<out-prefix>_erosion_rate_vs_ksn.png`` — log-log E against k_sn with the
  fit line(s) lsdtt3 wrote;
* ``<out-prefix>_basin_map_status.png`` / ``_basin_map_unit.png`` — the sample
  basins over a hillshade, by final status / dominant unit.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-erodibility-calibration``."""
    p = argparse.ArgumentParser(
        prog="lsdviz3-erodibility-calibration",
        description="Plot lsdtt3's erodibility calibration (E vs k_sn, basin map).")
    p.add_argument("--data-dir", required=True, help="Directory with the lsdtt3 outputs.")
    p.add_argument("--prefix", required=True, help="Chi run write_fname (e.g. wvcal_chi).")
    p.add_argument("--m-over-n", default=None,
                   help="m/n (float or suffix like 0p450); default: the only one present.")
    p.add_argument("--out-dir", required=True, help="Directory for the PNGs.")
    p.add_argument("--out-prefix", default=None, help="Output file prefix (default --prefix).")
    p.add_argument("--plots", default="scatter,map",
                   help="Comma list of scatter, map (default both).")
    p.add_argument("--map-colour-by", default="status,unit",
                   help="Map colouring(s): status, unit, or both (default).")
    p.add_argument("--dem", default=None,
                   help="DEM for the map hillshade (default <data-dir>/<read_fname> from the "
                   "run's ingested params, if it exists). Without one the map has a plain "
                   "background.")
    p.add_argument("--display-downsample", type=int, default=0,
                   help="Block-mean factor for drawing the DEM (display only); 0 = auto "
                   "(about 2000 pixels across).")
    p.add_argument("--unit-colours", default=None,
                   help="name,colour CSV or JSON {name: colour} (shared with the other "
                   "lithology plots).")
    p.add_argument("--credit", default=None, help="Explicit map credit (overrides discovery).")
    p.add_argument("--title", default=None, help="Title prefix.")
    p.add_argument("--dpi", type=int, default=200, help="Output DPI (default 200).")
    return p


def resolve_credit(directory, prefix, dem_path=None) -> str:
    """Geology (Macrostrat) / erosion rates (OCTOPUS) / DEM credits, one per line."""
    from ..io import lithology as L
    from ..io.erodibility_calibration import (
        discover_erosion_rate_attribution, erosion_rate_credit,
    )

    lines = []
    geo = L.discover_attribution(directory, prefix)
    if geo is not None:
        lines.append(L.condense_attribution(L.read_attribution(geo)))
    er = discover_erosion_rate_attribution(directory, prefix)
    if er is not None:
        lines.append(erosion_rate_credit(er))
    dem_params = L.discover_dem_params(dem_path) if dem_path else None
    dem_src = L.dem_source_from_params(dem_params) if dem_params else None
    if dem_src is None and dem_path:
        # Fetch params named <stem>_fetch_dem.param / <stem>_fetch.param.
        stem = Path(dem_path).stem.replace("_DEM", "")
        for cand in (f"{stem}_fetch_dem.param", f"{stem}_fetch.param",
                     f"{stem}_ingested_params.param"):
            q = Path(directory) / cand
            if q.exists() and L.dem_source_from_params(q):
                dem_src = L.dem_source_from_params(q)
                break
    if dem_src:
        lines.append(f"DEM: {dem_src}")
    return "\n".join(ln for ln in lines if ln)


def main(argv=None) -> int:
    """Entry point."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ..charts.erodibility_calibration import plot_erosion_rate_vs_ksn
    from ..io import lithology as L
    from ..io.erodibility_calibration import load_calibration_outputs

    args = build_parser().parse_args(argv)
    d = Path(args.data_dir)
    mn = args.m_over_n
    if mn is not None and "p" not in mn:
        mn = float(mn)
    data = load_calibration_outputs(d, args.prefix, mn)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = args.out_prefix or args.prefix
    plots = {x.strip() for x in args.plots.split(",") if x.strip()}
    uc = L.read_unit_colours(args.unit_colours) if args.unit_colours else None
    units_csv = d / f"{args.prefix}{L.SUFFIXES['units']}"
    unit_names = list(L.read_code_names(units_csv).values()) if units_csv.exists() else None
    written = []

    if "scatter" in plots:
        fig, _ = plot_erosion_rate_vs_ksn(data["samples"], data["fit"], unit_names=unit_names,
                                          unit_colours=uc, title=args.title)
        path = out / f"{stem}_erosion_rate_vs_ksn.png"
        fig.savefig(path, dpi=args.dpi)
        plt.close(fig)
        written.append(path)

    if "map" in plots:
        from ..render.base import save_figure
        from ..render.erodibility_calibration_map import render_calibration_basin_map
        from ..render.style import MapStyle

        if data["basins"] is None:
            raise SystemExit("No basin polygons (*_erodibility_calibration_basins.fgb).")
        dem = args.dem
        if dem is None:
            params = d / f"{args.prefix}{L.SUFFIXES['ingested_params']}"
            rf = L.read_param_file(params).get("read_fname") if params.exists() else None
            if rf and (d / rf).exists():
                dem = str(d / rf)
        factor = args.display_downsample
        if dem is not None and factor == 0:
            import rasterio

            with rasterio.open(dem) as src:
                factor = max(1, int(round(max(src.width, src.height) / 2000)))
        credit = args.credit if args.credit is not None else resolve_credit(d, args.prefix, dem)
        style = MapStyle(clean=True, transparent=False, figure_size="14c")
        for cb in [x.strip() for x in args.map_colour_by.split(",") if x.strip()]:
            fig = render_calibration_basin_map(
                data["samples"], data["basins"], dem=dem, colour_by=cb,
                unit_names=unit_names, unit_colours=uc, credit=credit,
                display_downsample=max(factor, 1), style=style)
            path = out / f"{stem}_basin_map_{cb}.png"
            save_figure(fig, path, style=style)
            written.append(path)

    for p in written:
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
