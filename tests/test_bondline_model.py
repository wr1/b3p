import pyvista as pv
from .test_build import built_blade


def test_bondline_model(built_blade):
    """Test if the bondline model is created."""
    workdir = built_blade["workdir"]
    bondline_vtu = list(workdir.glob("drape/*bondline.vtu"))
    assert bondline_vtu, "Bondline VTU file should exist"
    vtu = pv.read(bondline_vtu[0])
    assert vtu.n_points > 0, "Bondline VTU should contain points"
