from .test_build import run_build

import pyvista as pv


def test_bondline_model(run_build):
    """Test if the bondline model is created."""
    bondline_vtu = run_build["workdir"].glob("*bondline.vtu")
    print(bondline_vtu)

    # assert True
    # assert build_output_exists()
