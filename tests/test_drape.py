# import pytest
# from tests.test_geometry import get_yaml, build_blade
# from tests.test_structure import build_structure
# from tests.test_plybook import plybook

# # import pyvista as pv
# from b3p import drape_mesh

# # import numpy as np
# import os

# skeys = "slab_thickness_leading_edge_core,slab_thickness_shell_triax,slab_thickness_sparcap,slab_thickness_trailing_edge_core,slab_thickness_trailing_edge_ud,slab_thickness_transition_triax".split(
#     ","
# )


# @pytest.fixture(scope="session")
# def drape(get_yaml, build_blade, build_structure, plybook):
#     """Build a plybook"""
#     prefix = get_yaml["general"]["prefix"]
#     wdir = get_yaml["general"]["workdir"]
#     pref = os.path.join(wdir, prefix)
#     return [
#         drape_mesh.drape_mesh(
#             f"{pref}_{i if i != 'shell' else 'web'}.vtp",
#             plybook,
#             i,
#             f"{pref}_{i}.vtu",
#         )
#         for i in ["shell", "w0"]
#     ]


# # @pytest.mark.skip
# def test_drape_mesh_basic(drape):
#     """test some basic properties of the draped mesh, if the arrays are there and point 3000 is in the same place

#     Args:
#         drape ([UnstructuredGrid]): grid with ply data in cell_data"""
#     assert drape[0].n_points == 4000
#     assert drape[0].n_cells == 3960
#     assert drape[0].points[3000] == pytest.approx([-1.0641935, 1.800456, 75.757576])
#     assert list(
#         drape[0].cell_data
#     ) == "Area,Length,Normals,Volume,chord_length,d_abs_dist_from_bte,d_abs_dist_from_te,d_chord,d_le,d_le_r,d_lechord,d_lechord_abs,d_miny,d_rel_dist_from_te,d_sl,d_sla,d_te,d_te_offset,d_te_r,d_w0,d_w0_r,d_w1,d_w1_r,d_w2,d_w2_r,d_w3,d_w3_r,d_w4,d_w4_r,d_x,d_y,is_web,n_plies,ply_00000000_shell_triax,ply_00000001_shell_triax,ply_00000050_transition_triax,ply_00000051_transition_triax,ply_00000052_transition_triax,ply_00000053_transition_triax,ply_00000054_transition_triax,ply_00000055_transition_triax,ply_00000056_transition_triax,ply_00000057_transition_triax,ply_00000058_transition_triax,ply_00000059_transition_triax,ply_00000060_transition_triax,ply_00000061_transition_triax,ply_00000062_transition_triax,ply_00000063_transition_triax,ply_00000064_transition_triax,ply_00000065_transition_triax,ply_00000066_transition_triax,ply_00000067_transition_triax,ply_00000068_transition_triax,ply_00000069_transition_triax,ply_00000070_transition_triax,ply_00000071_transition_triax,ply_00000072_transition_triax,ply_00000073_transition_triax,ply_00000074_transition_triax,ply_00000075_transition_triax,ply_00000076_transition_triax,ply_00000077_transition_triax,ply_00000078_transition_triax,ply_00000079_transition_triax,ply_00000080_transition_triax,ply_00000081_transition_triax,ply_00000082_transition_triax,ply_00000083_transition_triax,ply_00000084_transition_triax,ply_00000085_transition_triax,ply_00000086_transition_triax,ply_00000087_transition_triax,ply_00000088_transition_triax,ply_00000089_transition_triax,ply_00000090_transition_triax,ply_00000091_transition_triax,ply_00000092_transition_triax,ply_00000093_transition_triax,ply_00000100_sparcap,ply_00000101_sparcap,ply_00000102_sparcap,ply_00000103_sparcap,ply_00000104_sparcap,ply_00000105_sparcap,ply_00000106_sparcap,ply_00000107_sparcap,ply_00000108_sparcap,ply_00000109_sparcap,ply_00000110_sparcap,ply_00000111_sparcap,ply_00000112_sparcap,ply_00000113_sparcap,ply_00000114_sparcap,ply_00000115_sparcap,ply_00000116_sparcap,ply_00000117_sparcap,ply_00000118_sparcap,ply_00000119_sparcap,ply_00000120_sparcap,ply_00000121_sparcap,ply_00000122_sparcap,ply_00000123_sparcap,ply_00000124_sparcap,ply_00000125_sparcap,ply_00000126_sparcap,ply_00000127_sparcap,ply_00000128_sparcap,ply_00000129_sparcap,ply_00000130_sparcap,ply_00000131_sparcap,ply_00000132_sparcap,ply_00000133_sparcap,ply_00000134_sparcap,ply_00000135_sparcap,ply_00000136_sparcap,ply_00000137_sparcap,ply_00000138_sparcap,ply_00000139_sparcap,ply_00000140_sparcap,ply_00000141_sparcap,ply_00000142_sparcap,ply_00000143_sparcap,ply_00000144_sparcap,ply_00000145_sparcap,ply_00000146_sparcap,ply_00000147_sparcap,ply_00000148_sparcap,ply_00000149_sparcap,ply_00000150_sparcap,ply_00000151_sparcap,ply_00000152_sparcap,ply_00000153_sparcap,ply_00000154_sparcap,ply_00000155_sparcap,ply_00000156_sparcap,ply_00000157_sparcap,ply_00000300_trailing_edge_ud,ply_00000301_trailing_edge_ud,ply_00000302_trailing_edge_ud,ply_00000303_trailing_edge_ud,ply_00000304_trailing_edge_ud,ply_00000305_trailing_edge_ud,ply_00000306_trailing_edge_ud,ply_00000307_trailing_edge_ud,ply_00000308_trailing_edge_ud,ply_00000309_trailing_edge_ud,ply_00000310_trailing_edge_ud,ply_00000311_trailing_edge_ud,ply_00000312_trailing_edge_ud,ply_00000313_trailing_edge_ud,ply_00000314_trailing_edge_ud,ply_00000315_trailing_edge_ud,ply_00000316_trailing_edge_ud,ply_00000317_trailing_edge_ud,ply_00000318_trailing_edge_ud,ply_00000319_trailing_edge_ud,ply_00000900_trailing_edge_core,ply_00000901_trailing_edge_core,ply_00000902_trailing_edge_core,ply_00000903_trailing_edge_core,ply_00000904_trailing_edge_core,ply_00000905_trailing_edge_core,ply_00001100_leading_edge_core,ply_00001101_leading_edge_core,ply_00001102_leading_edge_core,ply_00001103_leading_edge_core,ply_00006957_transition_triax,ply_00006958_transition_triax,ply_00006959_transition_triax,ply_00006960_transition_triax,ply_00006961_transition_triax,ply_00006962_transition_triax,ply_00006963_transition_triax,ply_00006964_transition_triax,ply_00006965_transition_triax,ply_00006966_transition_triax,ply_00006967_transition_triax,ply_00006968_transition_triax,ply_00006969_transition_triax,ply_00006970_transition_triax,ply_00006971_transition_triax,ply_00006972_transition_triax,ply_00006973_transition_triax,ply_00006974_transition_triax,ply_00006975_transition_triax,ply_00006976_transition_triax,ply_00006977_transition_triax,ply_00006978_transition_triax,ply_00006979_transition_triax,ply_00006980_transition_triax,ply_00006981_transition_triax,ply_00006982_transition_triax,ply_00006983_transition_triax,ply_00006984_transition_triax,ply_00006985_transition_triax,ply_00006986_transition_triax,ply_00006987_transition_triax,ply_00006988_transition_triax,ply_00006989_transition_triax,ply_00006990_transition_triax,ply_00006991_transition_triax,ply_00006992_transition_triax,ply_00006993_transition_triax,ply_00006994_transition_triax,ply_00006995_transition_triax,ply_00006996_transition_triax,ply_00006997_transition_triax,ply_00006998_transition_triax,ply_00006999_transition_triax,ply_00007000_transition_triax,ply_00007999_shell_triax,ply_00008000_shell_triax,radius,slab_thickness_leading_edge_core,slab_thickness_shell_triax,slab_thickness_sparcap,slab_thickness_trailing_edge_core,slab_thickness_trailing_edge_ud,slab_thickness_transition_triax,thickness,x_dir,y_dir,zone_ps,zone_ss".split(
#         ","
#     )


# def test_drape_slab_thicknesses(drape):
#     """Check the max thicknesses of the slabs"""
#     max_thicknesses = [drape[0].cell_data[k].max() for k in skeys]
#     assert max_thicknesses == pytest.approx([30.0, 5.0, 42.34, 50.0, 13.799999, 111.25])


# def test_drape_slab_volume(drape):
#     """Check the volumes of the slabs"""
#     volumes = [
#         (drape[0].cell_data[k] * drape[0].cell_data["Area"]).sum() for k in skeys
#     ]
#     assert volumes == pytest.approx(
#         [
#             5624.37386,
#             5071.47777,
#             7548.25637363090,
#             14595.484817536844,
#             823.525190,
#             15681.79849,
#         ]
#     )
