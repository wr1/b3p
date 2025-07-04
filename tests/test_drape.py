import pyvista as pv
import numpy as np
import json
import pytest
import logging

logger = logging.getLogger(__name__)


def test_laminate_number_of_plies(built_blade):
    """Test if the number of plies is correct."""
    workdir = built_blade["workdir"]
    temp_dir = built_blade["temp_dir"]
    joined = workdir / "drape" / "test_blade_joined.vtu"
    # joined = glob.glob(f"{workdir}/drape/*_joined.vtu")[-1]
    vtu = pv.read(joined)
    cell_arrays = [i for i in vtu.cell_data.keys() if "ply" in i]
    assert len(cell_arrays) == 348, f"Expected 348 plies, got {len(cell_arrays)}"

    sums = [float(vtu.cell_data[i][:, 1].sum()) for i in cell_arrays]
    logger.info(f"Computed ply sums: {sums}")

    json_path = temp_dir / "data" / "laminate_sums.json"
    if not json_path.exists():
        pytest.skip("Reference JSON 'laminate_sums.json' not found in tests/data/")
    with open(json_path, "r") as f:
        expected_sums = json.load(f)

    np.testing.assert_allclose(
        sums, expected_sums, atol=1e-4, err_msg="Ply sums do not match expected values"
    )
