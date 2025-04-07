import pyvista as pv
from .test_build import run_build


def test_bondline_model(run_build):
    """Test if the bondline model is created."""
    workdir = run_build["workdir"]
    bondline_vtu = list(workdir.glob("*bondline.vtu"))
    assert bondline_vtu, "Bondline VTU file should exist"
    vtu = pv.read(bondline_vtu[0])
    assert vtu.n_points > 0, "Bondline VTU should contain points"
