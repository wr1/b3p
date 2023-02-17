import pytest
from tests.test_geometry import get_yaml, build_blade
from tests.test_structure import build_structure
from tests.test_plybook import plybook
import pyvista as pv
from b3p import drape_mesh
import numpy as np
import os

skeys = "slab_thickness_leading_edge_core,slab_thickness_shell_triax,slab_thickness_sparcap,slab_thickness_trailing_edge_core,slab_thickness_trailing_edge_ud,slab_thickness_transition_triax".split(
    ","
)


@pytest.fixture(scope="session")
def drape(get_yaml, build_blade, build_structure, plybook):
    """Build a plybook"""
    prefix = get_yaml["general"]["prefix"]
    wdir = get_yaml["general"]["workdir"]
    pref = os.path.join(wdir, prefix)
    return [
        drape_mesh.drape_mesh(
            f"{pref}_{i if i != 'shell' else 'web'}.vtp",
            plybook,
            i,
            pref + f"_{i}.vtu",
        )
        for i in ["shell", "w0"]
    ]


def test_drape_mesh_basic(drape):
    """test some basic properties of the draped mesh, if the arrays are there and point 3000 is in the same place

    Args:
        drape ([UnstructuredGrid]): grid with ply data in cell_data"""
    assert drape[0].n_points == 4000
    assert drape[0].n_cells == 3960
    assert drape[0].points[3000] == pytest.approx([-1.0641935, 1.800456, 75.757576])
    assert list(
        drape[0].cell_data
    ) == "Area,Length,Normals,Volume,chord_length,d_abs_dist_from_bte,d_abs_dist_from_te,d_chord,d_le,d_le_r,d_lechord,d_lechord_abs,d_miny,d_rel_dist_from_te,d_sl,d_sla,d_te,d_te_offset,d_te_r,d_w0,d_w0_r,d_w1,d_w1_r,d_w2,d_w2_r,d_w3,d_w3_r,d_w4,d_w4_r,d_x,d_y,is_web,n_plies,ply_00000000_shell_triax,ply_00000001_shell_triax,ply_00000100_sparcap,ply_00000101_sparcap,ply_00000102_sparcap,ply_00000103_sparcap,ply_00000104_sparcap,ply_00000105_sparcap,ply_00000106_sparcap,ply_00000107_sparcap,ply_00000108_sparcap,ply_00000109_sparcap,ply_00000110_sparcap,ply_00000111_sparcap,ply_00000112_sparcap,ply_00000113_sparcap,ply_00000114_sparcap,ply_00000115_sparcap,ply_00000116_sparcap,ply_00000117_sparcap,ply_00000118_sparcap,ply_00000119_sparcap,ply_00000120_sparcap,ply_00000121_sparcap,ply_00000122_sparcap,ply_00000123_sparcap,ply_00000124_sparcap,ply_00000125_sparcap,ply_00000126_sparcap,ply_00000127_sparcap,ply_00000128_sparcap,ply_00000129_sparcap,ply_00000130_sparcap,ply_00000131_sparcap,ply_00000132_sparcap,ply_00000133_sparcap,ply_00000134_sparcap,ply_00000135_sparcap,ply_00000136_sparcap,ply_00000137_sparcap,ply_00000138_sparcap,ply_00000139_sparcap,ply_00000140_sparcap,ply_00000141_sparcap,ply_00000142_sparcap,ply_00000143_sparcap,ply_00000144_sparcap,ply_00000145_sparcap,ply_00000146_sparcap,ply_00000147_sparcap,ply_00000148_sparcap,ply_00000149_sparcap,ply_00000150_sparcap,ply_00000151_sparcap,ply_00000152_sparcap,ply_00000153_sparcap,ply_00000154_sparcap,ply_00000155_sparcap,ply_00000156_sparcap,ply_00000157_sparcap,ply_00000300_trailing_edge_ud,ply_00000301_trailing_edge_ud,ply_00000302_trailing_edge_ud,ply_00000303_trailing_edge_ud,ply_00000304_trailing_edge_ud,ply_00000305_trailing_edge_ud,ply_00000306_trailing_edge_ud,ply_00000307_trailing_edge_ud,ply_00000308_trailing_edge_ud,ply_00000309_trailing_edge_ud,ply_00000310_trailing_edge_ud,ply_00000311_trailing_edge_ud,ply_00000312_trailing_edge_ud,ply_00000313_trailing_edge_ud,ply_00000314_trailing_edge_ud,ply_00000315_trailing_edge_ud,ply_00000316_trailing_edge_ud,ply_00000317_trailing_edge_ud,ply_00000318_trailing_edge_ud,ply_00000319_trailing_edge_ud,ply_00000600_transition_triax,ply_00000601_transition_triax,ply_00000602_transition_triax,ply_00000603_transition_triax,ply_00000604_transition_triax,ply_00000605_transition_triax,ply_00000606_transition_triax,ply_00000607_transition_triax,ply_00000608_transition_triax,ply_00000609_transition_triax,ply_00000610_transition_triax,ply_00000611_transition_triax,ply_00000612_transition_triax,ply_00000613_transition_triax,ply_00000614_transition_triax,ply_00000615_transition_triax,ply_00000616_transition_triax,ply_00000617_transition_triax,ply_00000618_transition_triax,ply_00000619_transition_triax,ply_00000620_transition_triax,ply_00000621_transition_triax,ply_00000622_transition_triax,ply_00000623_transition_triax,ply_00000624_transition_triax,ply_00000625_transition_triax,ply_00000626_transition_triax,ply_00000627_transition_triax,ply_00000628_transition_triax,ply_00000629_transition_triax,ply_00000630_transition_triax,ply_00000631_transition_triax,ply_00000632_transition_triax,ply_00000633_transition_triax,ply_00000634_transition_triax,ply_00000635_transition_triax,ply_00000636_transition_triax,ply_00000637_transition_triax,ply_00000638_transition_triax,ply_00000639_transition_triax,ply_00000640_transition_triax,ply_00000641_transition_triax,ply_00000642_transition_triax,ply_00000643_transition_triax,ply_00000644_transition_triax,ply_00000645_transition_triax,ply_00000646_transition_triax,ply_00000647_transition_triax,ply_00000648_transition_triax,ply_00000649_transition_triax,ply_00000650_transition_triax,ply_00000651_transition_triax,ply_00000652_transition_triax,ply_00000653_transition_triax,ply_00000654_transition_triax,ply_00000655_transition_triax,ply_00000656_transition_triax,ply_00000657_transition_triax,ply_00000658_transition_triax,ply_00000659_transition_triax,ply_00000660_transition_triax,ply_00000661_transition_triax,ply_00000662_transition_triax,ply_00000663_transition_triax,ply_00000664_transition_triax,ply_00000665_transition_triax,ply_00000666_transition_triax,ply_00000667_transition_triax,ply_00000668_transition_triax,ply_00000669_transition_triax,ply_00000670_transition_triax,ply_00000671_transition_triax,ply_00000672_transition_triax,ply_00000673_transition_triax,ply_00000674_transition_triax,ply_00000675_transition_triax,ply_00000676_transition_triax,ply_00000677_transition_triax,ply_00000678_transition_triax,ply_00000679_transition_triax,ply_00000680_transition_triax,ply_00000681_transition_triax,ply_00000682_transition_triax,ply_00000683_transition_triax,ply_00000684_transition_triax,ply_00000685_transition_triax,ply_00000686_transition_triax,ply_00000687_transition_triax,ply_00000688_transition_triax,ply_00000900_trailing_edge_core,ply_00000901_trailing_edge_core,ply_00000902_trailing_edge_core,ply_00000903_trailing_edge_core,ply_00001100_leading_edge_core,ply_00001101_leading_edge_core,ply_00001102_leading_edge_core,ply_00001103_leading_edge_core,ply_00007999_shell_triax,ply_00008000_shell_triax,radius,slab_thickness_leading_edge_core,slab_thickness_shell_triax,slab_thickness_sparcap,slab_thickness_trailing_edge_core,slab_thickness_trailing_edge_ud,slab_thickness_transition_triax,thickness,x_dir,y_dir,zone_ps,zone_ss".split(
        ","
    )


def test_drape_slab_thicknesses(drape):
    """Check the max thicknesses of the slabs"""
    max_thicknesses = [drape[0].cell_data[k].max() for k in skeys]
    assert max_thicknesses == pytest.approx([30.0, 5.0, 42.34, 50.0, 13.799999, 111.25])


def test_drape_slab_volume(drape):
    """Check the volumes of the slabs"""
    volumes = [
        (drape[0].cell_data[k] * drape[0].cell_data["Area"]).sum() for k in skeys
    ]
    assert volumes == pytest.approx(
        [
            5624.37386,
            5071.47777,
            7652.12303,
            15432.86749,
            823.525190,
            15681.79849,
        ]
    )
