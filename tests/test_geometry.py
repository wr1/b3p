from .conftest import built_blade
import os
import pytest
import pyvista as pv
import pandas as pd
from io import StringIO
import glob


@pytest.fixture(scope="session")
def load_geometry(built_blade):
    """Fixture to load the joined geometry from the build."""
    workdir = built_blade["workdir"]
    joined_vtu = glob.glob(f"{workdir}/drape/test_blade_joined.vtu")[0]
    return pv.read(joined_vtu)


def test_planform(built_blade):
    """Test if the planform CSV matches expected values."""
    workdir = built_blade["workdir"]
    temp_dir = built_blade["temp_dir"]  # Access the parent temp dir
    pref = glob.glob(f"{workdir}/mesh/*50.csv")[0]
    planform = pd.read_csv(pref, sep=";")

    # Load the reference CSV from the copied data directory
    csv_path = temp_dir / "data" / "test_blade_csv.csv"
    if not csv_path.exists():
        pytest.skip("Reference CSV 'test_blade_csv.csv' not found in tests/data/")
    expected_df = pd.read_csv(csv_path, sep=";")

    print(expected_df)
    pd.testing.assert_frame_equal(expected_df, planform)


def test_geometry_bounding_box(load_geometry):
    """Test if the geometry bounding box is correct."""
    assert load_geometry.bounds == pytest.approx(
        (
            -3.0369677543640137,
            2.499872922897339,
            -3.3991336822509766,
            3.0669267177581787,
            0.0,
            100.0,
        ),
        rel=1e-5,
    ), "Bounding box does not match expected values"
