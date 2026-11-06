"""k_sn colour-scale limits (vmin / vmax in k_sn units on a log10 scale)."""

from __future__ import annotations

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import LineCollection
from shapely.geometry import LineString

from lsdviztools3.colour_limits import colour_limits, gmt_extend


def _segments():
    """Two basins, k_sn 5 / 50 / 500 / 2000."""
    rows = []
    for b, ksns in ((0, (5.0, 50.0)), (1, (500.0, 2000.0))):
        z, chi = 100.0, 0.0
        for i, k in enumerate(ksns):
            rows.append(dict(basin_index=b, outlet_flat=0, segment_id=i, parent_segment_id=i - 1,
                             branch_id=0, downstream_flat=i, upstream_flat=i + 1,
                             downstream_chi=chi, upstream_chi=chi + 1.0,
                             downstream_elevation=z, upstream_elevation=z + k * 1.0,
                             k_sn=k, max_abs_residual=0.0, n_nodes=10))
            chi += 1.0
            z += k
    return pd.DataFrame(rows)


def test_colour_limits_in_data_units_on_a_log_scale():
    c = np.log10([5.0, 50.0, 500.0, 2000.0])
    lo, hi, ext = colour_limits(c, 10.0, 1000.0, log_scale=True)
    assert lo == pytest.approx(1.0) and hi == pytest.approx(3.0) and ext == "both"
    lo, hi, ext = colour_limits(c, None, 1000.0, log_scale=True)
    assert lo == pytest.approx(np.log10(5.0)) and ext == "max"
    lo, hi, ext = colour_limits(c, None, None, log_scale=True)
    assert ext == "neither"
    lo, hi, ext = colour_limits([1.0, 9.0], 2.0, None, log_scale=False)
    assert lo == 2.0 and hi == 9.0 and ext == "min"


def test_colour_limits_refuse_bad_values():
    with pytest.raises(ValueError, match="> 0"):
        colour_limits(np.log10([5.0]), 0.0, 10.0, log_scale=True)
    with pytest.raises(ValueError, match="below"):
        colour_limits(np.log10([5.0]), 100.0, 10.0, log_scale=True)


def test_gmt_extend_modifiers():
    assert gmt_extend("max") == "+ef" and gmt_extend("min") == "+eb"
    assert gmt_extend("both") == "+e" and gmt_extend("neither") == ""


def test_stacked_vmax_sets_the_norm_and_extend():
    from lsdviztools3.charts import plot_chi_stacked
    fig, ax = plot_chi_stacked(_segments(), vmax=100.0)
    lc = next(c for c in ax.collections if isinstance(c, LineCollection))
    assert lc.norm.vmax == pytest.approx(2.0)
    assert lc.norm.vmin == pytest.approx(np.log10(5.0))
    # segments above vmax take the top colour
    top = lc.cmap(1.0)
    colours = lc.cmap(lc.norm(lc.get_array()))
    assert np.allclose(colours[-1], top) and np.allclose(colours[-2], top)
    cbar = fig.axes[-1]
    assert cbar is not ax
    assert fig.axes[-1].get_label() == "<colorbar>"
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_profile_vmin_vmax_and_colourbar_extend():
    from lsdviztools3.charts import plot_chi_profile
    fig, ax = plot_chi_profile(_segments(), basin=1, vmin=100.0, vmax=1000.0)
    lc = next(c for c in ax.collections if isinstance(c, LineCollection))
    assert (lc.norm.vmin, lc.norm.vmax) == pytest.approx((2.0, 3.0))
    cb = lc.colorbar
    assert cb is not None and cb.extend == "max"     # 2000 > 1000, 500 > 100
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_profile_without_limits_is_unchanged():
    from lsdviztools3.charts import plot_chi_profile
    fig, ax = plot_chi_profile(_segments(), basin=0)
    lc = next(c for c in ax.collections if isinstance(c, LineCollection))
    assert (lc.norm.vmin, lc.norm.vmax) == pytest.approx((np.log10(5.0), np.log10(50.0)))
    assert lc.colorbar.extend == "neither"
    import matplotlib.pyplot as plt
    plt.close(fig)


def _lines():
    return gpd.GeoDataFrame(
        {"k_sn": [5.0, 50.0, 500.0]},
        geometry=[LineString([(0, i), (1, i + 0.5)]) for i in range(3)], crs="EPSG:4326")


def test_render_channels_cpt_range_and_extend(monkeypatch):
    import pygmt
    from lsdviztools3.render import channels, decorations
    from lsdviztools3.render.style import MapStyle

    calls = {}
    real = pygmt.makecpt

    def spy(**kw):
        calls.update(kw)
        return real(**kw)

    monkeypatch.setattr(pygmt, "makecpt", spy)
    seen = {}
    real_cb = decorations.add_colourbar

    def spy_cb(fig, **kw):
        seen.update(kw)
        return real_cb(fig, **kw)

    monkeypatch.setattr(channels, "add_colourbar", spy_cb)
    channels.render_channels(_lines(), color_by="k_sn", coordinates="geographic",
                             style=MapStyle(), colorbar=True, vmin=10.0, vmax=100.0)
    assert calls["series"] == pytest.approx([1.0, 2.0])
    assert calls["background"] is True
    assert seen["extend"] == "both"


def test_render_channels_writes_png_with_limits(tmp_path):
    from lsdviztools3.render import channels
    from lsdviztools3.render.style import MapStyle
    fig = channels.render_channels(_lines(), color_by="k_sn", coordinates="geographic",
                                   style=MapStyle(), colorbar=True, vmax=100.0)
    out = tmp_path / "ksn.png"
    fig.savefig(str(out))
    assert out.exists() and out.stat().st_size > 1000


def test_cli_refuses_non_positive_limits(tmp_path):
    from lsdviztools3.cli import plot_chi
    with pytest.raises(SystemExit):
        plot_chi.main(["--data-dir", str(tmp_path), "--prefix", "x", "--m-over-n", "0.45",
                       "--out-dir", str(tmp_path), "--ksn-vmax", "0"])
    with pytest.raises(SystemExit):
        plot_chi.main(["--data-dir", str(tmp_path), "--prefix", "x", "--m-over-n", "0.45",
                       "--out-dir", str(tmp_path), "--ksn-vmin", "100", "--ksn-vmax", "10"])
