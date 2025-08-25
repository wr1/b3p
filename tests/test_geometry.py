import pytest
import pyvista as pv
import pandas as pd
import logging

logger = logging.getLogger(__name__)


@pytest.fixture(scope="session")
def load_geometry(built_blade):
    """Fixture to load the joined geometry from the build."""
    workdir = built_blade["workdir"]
    joined_vtu = workdir / "drape" / "test_blade_joined.vtu"
    return pv.read(joined_vtu)


def test_planform(built_blade):
    """Test if the planform CSV matches expected values."""
    workdir = built_blade["workdir"]
    temp_dir = built_blade["temp_dir"]
    planform = workdir / "mesh" / "test_blade_sca_50.csv"
    planform = pd.read_csv(planform, sep=";")

    csv_path = temp_dir / "data" / "test_blade_csv.csv"
    if not csv_path.exists():
        pytest.skip("Reference CSV 'test_blade_csv.csv' not found in tests/data/")
    expected_df = pd.read_csv(csv_path, sep=";")

    logger.info(f"Expected DataFrame:\n{expected_df}")
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


def test_geometry_n_points(load_geometry):
    """Test if the geometry has the expected number of points."""
    assert load_geometry.n_points == 10200, "Number of points does not match expected value"


def test_geometry_n_cells(load_geometry):
    """Test if the geometry has the expected number of cells."""
    assert load_geometry.n_cells == 5049, "Number of cells does not match expected value"


def test_geometry_cell_types(load_geometry):
    """Test if the geometry contains only quad cells."""
    assert (
        load_geometry.celltypes == [9]
    ), "Geometry should only contain quad cells (pyvista.VTK_QUAD)"

def test_geometry_cell_data(load_geometry):
    """Test if the geometry has the expected cell data arrays."""
    expected_arrays = [
        "Cell Normal",
        "ply_angle",
        "n_plies",
        "is_le",
        "is_te",
        "is_web",
        "is_inner",
        "is_outer",
        "is_bondline",
        "is_profile",
        "is_fill",
        "is_core",
        "id",
        "Material",
        "rel_thick",
        "abs_thick",
        "area",
    ]
    for array in expected_arrays:
        assert array in load_geometry.cell_data, f"Expected cell data array '{array}' not found"

    # Check for layer arrays
    layer_arrays = [name for name in load_geometry.cell_data if name.startswith("layer_")]
    assert len(layer_arrays) > 0, "No layer arrays found in cell data"
