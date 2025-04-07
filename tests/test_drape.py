import glob
import pyvista as pv
import numpy as np
from .conftest import run_build  # Explicitly import run_build
import json


def test_laminate_number_of_plies(run_build):
    """Test if the number of plies is correct."""
    workdir = run_build["workdir"]
    temp_dir = run_build["temp_dir"]  # Access the parent temp dir
    joined = glob.glob(f"{workdir}/*_joined.vtu")[-1]
    vtu = pv.read(joined)
    cell_arrays = [i for i in vtu.cell_data.keys() if "ply" in i]
    assert len(cell_arrays) == 348, f"Expected 348 plies, got {len(cell_arrays)}"

    sums = [float(vtu.cell_data[i][:, 1].sum()) for i in cell_arrays]
    print(sums)

    # Load expected sums from JSON file in the temporary data directory
    json_path = temp_dir / "data" / "laminate_sums.json"
    if not json_path.exists():
        pytest.skip("Reference JSON 'laminate_sums.json' not found in tests/data/")
    with open(json_path, "r") as f:
        expected_sums = json.load(f)

    np.testing.assert_allclose(
        sums, expected_sums, atol=1e-4, err_msg="Ply sums do not match expected values"
    )
