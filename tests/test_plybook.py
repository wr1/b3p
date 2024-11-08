# import pytest
# from tests.test_geometry import get_yaml
# import pyvista as pv
# from b3p import build_plybook
# import numpy as np


# @pytest.fixture(scope="session")
# def plybook(get_yaml):
#     """Build a plybook"""
#     return build_plybook.lamplan2plies(get_yaml)


# def test_build_plybook_nslabs(plybook):
#     """Check a point somewhere in the mesh, not the root or tip"""
#     assert len(plybook) == 10


# def test_plybook_numbering(plybook):
#     assert plybook[0]["numbering"][:12] == [
#         100,
#         101,
#         102,
#         103,
#         104,
#         105,
#         106,
#         107,
#         108,
#         109,
#         110,
#         111,
#     ]
#     assert plybook[5]["numbering"] == [1100, 1101, 1102, 1103]


# def test_plybook_names(plybook):
#     names = [i["name"] for i in plybook]
#     assert names == [
#         "sparcap",
#         "trailing_edge_ud",
#         "shell_triax",
#         "transition_triax",
#         "trailing_edge_core",
#         "leading_edge_core",
#         "web3_biax",
#         "web3_core",
#         "web4_biax",
#         "web4_core",
#     ]


# def test_plybook_thickness(plybook):
#     assert plybook[6]["stack"] == [
#         [7, 1.2, 0.0, 123.0],
#         [7, 1.2, 0.0, 123.0],
#         [7, 1.2, 0.0, 123.0],
#         [7, 1.2, 0.0, 123.0],
#     ]
