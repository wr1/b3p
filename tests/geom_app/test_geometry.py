#!/usr/bin/env python3
"""Tests for geometry.py in geom_app."""

import pytest
import numpy as np
from pathlib import Path
import tempfile
import json

# Import the functions to test
from src.b3p.geom_app.geometry import optspace, build_blade_geometry, blade


def test_optspace():
    """Test the optspace function."""
    n_points = 10
    result = optspace(n_points)
    assert len(result) == n_points
    assert result[0] == 0.0
    assert result[-1] == 1.0
    assert all(result[i] <= result[i + 1] for i in range(len(result) - 1))


def test_build_blade_geometry():
    """Test build_blade_geometry with sample config."""
    sample_config = {
        "planform": {
            "chord": [[0, 1], [1, 0.5]],
            "thickness": [[0, 0.2], [1, 0.1]],
            "twist": [[0, 0], [1, 10]],
            "dx": [[0, 0], [1, 0]],
            "dy": [[0, 0], [1, 0]],
            "z": [[0, 0], [1, 10]],
            "npchord": 50,
            "npspan": 20,
        },
        "aero": {
            "airfoils": {
                0.2: {"xy": [[0, 0], [0.5, 0.1], [1, 0]]},
                0.1: {"xy": [[0, 0], [0.5, 0.05], [1, 0]]},
            }
        },
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)
        blade_obj = build_blade_geometry(sample_config, workdir)
        assert isinstance(blade_obj, blade)
        assert (workdir / "blade_geometry.vtp").exists()
        assert (workdir / "blade_geometry.pck").exists()
        assert (workdir / "blade_geometry_variables.json").exists()


def test_blade_initialization():
    """Test blade class initialization."""
    chord = [[0, 1], [1, 0.5]]
    thickness = [[0, 0.2], [1, 0.1]]
    twist = [[0, 0], [1, 10]]
    dx = [[0, 0], [1, 0]]
    dy = [[0, 0], [1, 0]]
    z = [[0, 0], [1, 10]]
    airfoils = {0.2: {"xy": [[0, 0], [0.5, 0.1], [1, 0]]}}
    chordwise_sampling = np.linspace(0, 1, 50)

    blade_obj = blade(
        chord, thickness, twist, dx, dy, z, airfoils, chordwise_sampling, np_spanwise=20
    )
    assert hasattr(blade_obj, "sections")
    assert len(blade_obj.sections) > 0
