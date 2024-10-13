from .test_build import run_test_build, workdir
import os
import pytest
import pyvista as pv


@pytest.fixture(scope="session")
def load_geometry(run_test_build):
    """Fixture to run step s2. Checks if output exists to avoid re-running."""
    return pv.read(os.path.join(workdir, "test_blade_joined.vtu"))


def test_geometry_bounding_box(load_geometry):
    """Test if the geometry bounding box is correct."""

    assert load_geometry.bounds == (
        -3.0369677543640137,
        2.499872922897339,
        -3.3991336822509766,
        3.0669267177581787,
        0.0,
        100.0,
    )
