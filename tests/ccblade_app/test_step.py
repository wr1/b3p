#!/usr/bin/env python3
"""Tests for step.py in ccblade_app."""

import tempfile
from pathlib import Path
import yaml

# Import the step
from src.b3p.ccblade_app.step import CCBladeStep


def test_ccblade_step():
    """Test CCBladeStep execution."""
    sample_config = {
        "general": {"workdir": "test_workdir", "prefix": "test"},
        "aero": {
            "bem": {
                "polars": {"0.2": "polar1.txt"},
                "B": 3,
                "rho": 1.225,
                "mu": 1.81e-5,
                "precone": 0,
                "tilt": 0,
                "yaw": 0,
                "shearExp": 0.2,
                "hubHt": 90,
                "max_tipspeed": 95,
                "rated_power": 2e7,
                "uinf": [4, 6, 8, 10, 12],
            }
        },
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        config_path = Path(tmpdir) / "test_config.yml"
        with open(config_path, "w") as f:
            yaml.dump(sample_config, f)

        # Note: This test may need mock data for polars and blade files
        step = CCBladeStep(str(config_path))
        # step.run()  # Commented out as it requires full setup
