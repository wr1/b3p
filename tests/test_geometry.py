from b3p import build_blade_geometry
import os
from ruamel import yaml
import pytest
import pyvista as pv
import numpy as np

absolute_path = os.path.dirname(__file__)


@pytest.fixture(scope="session")
def get_yaml():
    fn = os.path.join(absolute_path, "../examples/blade_test_portable.yml")
    yml = yaml.YAML(typ="rt")
    return yml.load(open(fn, "rb"))


@pytest.fixture(scope="session")
def build_blade(get_yaml):
    """Build a blade, using the yml input file in examples/"""
    return build_blade_geometry.build_blade_geometry(get_yaml)


def test_build_blade_geometry_key_parameters(build_blade):
    """Check some of the global parameters defining the geometry"""
    assert build_blade.np_chordwise == 100


def test_build_blade_geometry_planform(build_blade):
    """"""
    assert build_blade.dx[1][12] == -0.18536853892047267


def test_build_blade_geometry_mesh(build_blade):
    assert build_blade.poly.points[50, 1] == -2.340645738717508
    assert build_blade.poly.points.max(axis=0).tolist() == [
        2.514195211238277,
        3.6879477184670346,
        126.0,
    ]
