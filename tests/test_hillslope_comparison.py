import numpy as np
import pytest

from lsdviztools3.charts.hillslope_profiles import (
    plot_hillslope_profiles,
    plot_profile_residual_panels,
    profile_derivatives,
)
from lsdviztools3.charts.hillslope_residuals import (
    residual_statistics,
    summarise_residuals,
)
from lsdviztools3.charts.hillslope_transient_residuals import (
    load_transient_residuals,
)
from lsdviztools3.render.hillslope_comparison_animation import discover_frame_pairs
from lsdviztools3.render.hillslope_comparison_animation import (
    render_hillslope_comparison_animation,
)


def test_profile_derivatives_and_plot_use_array_only():
    pytest.importorskip("matplotlib")
    x = np.arange(5.0)
    z = x * x
    gradient, curvature = profile_derivatives(x, z)
    np.testing.assert_allclose(gradient[1:-1], [2, 4, 6])
    np.testing.assert_allclose(curvature[1:-1], [1.5, 2.0, 1.5])
    fig, axes = plot_hillslope_profiles(
        x, {"Newton": z, "cellular": z + 0.1},
        analytic={"linear": {"elevation": z, "gradient": 2 * x,
                               "curvature": np.full_like(x, 2.0)}},
    )
    assert len(axes) == 4
    fig.clf()
    fig, axes = plot_profile_residual_panels(x, {"a": z, "b": -z})
    assert len(axes) == 2
    fig.clf()


def test_residual_statistics_and_hillslope_scope():
    values = np.array([-2.0, 1.0, np.nan, 3.0])
    stats = residual_statistics(values)
    assert stats["n"] == 3
    assert stats["linf"] == 3.0
    table = summarise_residuals({"cellular - Newton": values},
                                hillslope_mask=[True, False, True, False])
    assert set(table["scope"]) == {"full", "hillslope"}
    assert table.loc[table.scope == "hillslope", "n"].item() == 1


def test_transient_csv_requires_plan_columns(tmp_path):
    path = tmp_path / "residuals.csv"
    path.write_text("time_yr,mean,rms,linf,p05,median,p95,hillslope_mean,hillslope_rms\n"
                    "10,0,1,2,-1,0,1,0,1\n"
                    "0,0,0.5,1,-.5,0,.5,0,.5\n")
    frame = load_transient_residuals(path)
    assert frame.time_yr.tolist() == [0, 10]


def test_frame_pairing_rejects_missing_partner(tmp_path):
    newton = tmp_path / "newton"
    cellular = tmp_path / "cellular"
    newton.mkdir()
    cellular.mkdir()
    (newton / "surface_iteration_0.tif").touch()
    (cellular / "surface_iteration_1.tif").touch()
    with pytest.raises(ValueError, match="frame sets do not match"):
        discover_frame_pairs(newton, cellular)


def test_animation_smoke_with_tiny_geotiffs(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    pytest.importorskip("matplotlib")
    pytest.importorskip("PIL")
    from rasterio.transform import from_origin

    newton = tmp_path / "newton"
    cellular = tmp_path / "cellular"
    newton.mkdir()
    cellular.mkdir()
    transform = from_origin(0, 2, 1, 1)
    for i in (0, 1):
        for directory, offset in ((newton, 0.0), (cellular, 0.25)):
            with rasterio.open(directory / f"surface_iteration_{i}.tif", "w",
                               driver="GTiff", height=2, width=2, count=1,
                               dtype="float32", crs="EPSG:32630",
                               transform=transform) as dst:
                dst.write(np.full((2, 2), i + offset, dtype="float32"), 1)
    output = render_hillslope_comparison_animation(newton, cellular,
                                                    tmp_path / "comparison.gif")
    assert output.is_file()
