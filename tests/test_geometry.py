from b3p import build_blade_geometry, yml_portable
import os
from ruamel import yaml
import pytest
import pyvista as pv
import numpy as np

absolute_path = os.path.dirname(__file__)


@pytest.fixture(scope="session")
def get_yaml():
    fn = os.path.join(absolute_path, "../examples/blade_test.yml")
    return yml_portable.yaml_make_portable(fn)


@pytest.fixture(scope="session")
def build_blade(get_yaml):
    """Build a blade, using the yml input file in examples/"""
    return build_blade_geometry.build_blade_geometry(get_yaml)


def test_yaml_loadable(get_yaml):
    """check that the yaml file contains the right sections"""
    sections = ["general", "aero", "materials", "mesh", "mesh2d", "loads", "laminates"]
    for s in sections:
        assert s in get_yaml


def test_build_blade_planform_interpolate(build_blade):
    """Checks some values of the interpolated planform, in order to check consistency"""
    assert build_blade.dx[1][12] == -0.011400483062328838
    assert build_blade.dy[1][25] == 0.2716041985416647
    assert build_blade.chord[1][4] == 5.234238195149214


def test_build_blade_geometry_point_on_mesh(build_blade):
    """Check a point somewhere in the mesh, not the root or tip"""
    assert build_blade.poly.points[50, 1] == -2.4598988864366738


def test_build_blade_mesh_generation_bounding_box(build_blade):
    """Check that the blade sits in the exact same bounding box"""
    assert build_blade.poly.points.max(axis=0).tolist() == [
        2.499435016496742,
        3.0673998220471317,
        126.0,
    ]


def test_build_blade_mesh_size(build_blade):
    """Check that the blade has the right number of points and cells"""
    assert build_blade.poly.n_points == 7000
    assert build_blade.poly.n_cells == 6900
