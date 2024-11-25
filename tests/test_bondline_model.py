from .test_build import run_test_build, workdir

import pyvista as pv


def test_bondline_model(run_test_build):
    """Test if the bondline model is created."""
    bondline_vtu = workdir.glob("*bondline.vtu")
    print(bondline_vtu)

    # assert True
    # assert build_output_exists()
