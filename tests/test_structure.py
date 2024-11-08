# import pytest
# from tests.test_geometry import get_yaml
# import pyvista as pv
# from b3p import build_blade_structure
# import numpy as np


# @pytest.fixture(scope="session")
# def build_structure(get_yaml):
#     """Build a blade structure, using the yml input file in examples/"""
#     return pv.PolyData(build_blade_structure.build_blade_structure(get_yaml).poly)


# def test_blade_structure_point_on_mesh(build_structure):
#     """Check a point somewhere in the mesh, not the root or tip"""
#     assert build_structure.points[900] == pytest.approx(
#         [-0.8162439, -2.6600213, 22.222221]
#     )


# def test_blade_structure_bounding_box(build_structure):
#     """Check that the blade sits in the exact same bounding box"""
#     assert build_structure.bounds == pytest.approx(
#         (
#             -3.052877902984619,
#             2.5136566162109375,
#             -2.783158302307129,
#             3.6869149208068848,
#             0.0,
#             100.0,
#         )
#     )


# def test_blade_structure_has_coordinate_systems(build_structure):
#     """Check if all point coordinate systems are present"""
#     assert build_structure.point_data.keys() == [
#         "d_te",
#         "radius",
#         "d_rel_dist_from_te",
#         "d_abs_dist_from_te",
#         "d_abs_dist_from_bte",
#         "chord_length",
#         "d_miny",
#         "d_le",
#         "zone_ss",
#         "zone_ps",
#         "d_te_r",
#         "d_le_r",
#         "d_chord",
#         "d_lechord",
#         "d_lechord_abs",
#         "d_x",
#         "d_y",
#         "d_sl",
#         "d_sla",
#         "is_web",
#         "d_w0",
#         "d_w0_r",
#         "d_w1",
#         "d_w1_r",
#         "d_w2",
#         "d_w2_r",
#         "d_w3",
#         "d_w3_r",
#         "d_w4",
#         "d_w4_r",
#         "d_te_offset",
#     ]


# def test_blade_structure_check_connectivity_size(build_structure):
#     assert build_structure.faces.shape[0] == 19800


# def test_blade_structure_check_connectivity_mid(build_structure):
#     assert np.all(
#         build_structure.faces[4000:4012]
#         == np.array([4, 800, 840, 841, 801, 4, 801, 841, 842, 802, 4, 802])
#     )


# def test_blade_structure_check_connectivity_end(build_structure):
#     assert np.all(
#         build_structure.faces[-12:]
#         == np.array([3998, 3958, 4, 3958, 3998, 3999, 3959, 4, 3959, 3999, 3960, 3920])
#     )


# def test_blade_structure_check_connectivity_start(build_structure):
#     assert np.all(
#         build_structure.faces[:12] == np.array([4, 0, 40, 41, 1, 4, 1, 41, 42, 2, 4, 2])
#     )
