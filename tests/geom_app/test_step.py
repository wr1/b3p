#!/usr/bin/env python3
"""Tests for step.py in geom_app."""

import pytest
import tempfile
from pathlib import Path
import yaml

# Import the step
from src.b3p.geom_app.step import GeometryStep


def test_geometry_step():
    """Test GeometryStep execution."""
    sample_config = {
        "general": {"workdir": "test_workdir"},
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
            }
        },
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        config_path = Path(tmpdir) / "test_config.yml"
        with open(config_path, "w") as f:
            yaml.dump(sample_config, f)

        step = GeometryStep(str(config_path))
        step.run()

        workdir = Path(tmpdir) / "test_workdir"
        assert (workdir / "blade_geometry.vtp").exists()
        assert (workdir / "blade_geometry.pck").exists()
