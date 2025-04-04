import json
import os
import subprocess
from pathlib import Path

import pytest
import yaml
from b3p.cli2 import AppState, TwoDApp

# from b3p.utils import copy_example_files
from utils import prep_temp_dir


@pytest.fixture(scope="session")
def tmp_dir(tmp_path_factory):
    """Create a temporary directory with example files for the test session."""
    tmp_dir = tmp_path_factory.mktemp("test_2d")
    prep_temp_dir(tmp_dir)  # Assuming this correctly sets up the environment
    return tmp_dir


@pytest.fixture
def two_d_app():
    """Provide a fresh TwoDApp instance for each test."""
    state = AppState.get_instance()
    return TwoDApp(state)


def test_mesh2d_output(two_d_app, tmp_dir):
    """Test that mesh2d generates a non-None result for test.yml."""
    yml_path = tmp_dir / "test.yml"
    print(f"Testing mesh2d with {yml_path}")
    result = two_d_app.mesh2d(yml_path)
    print(f"mesh2d result: {result}")
    assert result is not None, "mesh2d should return a non-None result"


def test_anba4_output_json(two_d_app, tmp_dir):
    """Test that ANBA4 produces at least one JSON output file with valid content."""
    yml_path = tmp_dir / "test.yml"

    # Load YAML to get workdir
    with open(yml_path, "r") as f:
        config = yaml.safe_load(f)
    workdir = config["general"]["workdir"]

    # Find JSON output files
    output_jsons = list(Path(tmp_dir / workdir).glob("msec_*.json"))
    assert output_jsons, "ANBA4 should create at least one output JSON file"

    # Load and inspect the first JSON
    with open(output_jsons[0], "r") as f:
        json_data = json.load(f)
    print(f"Found JSON files: {output_jsons}")
    print(f"First JSON content: {json_data}")


def test_anba4_output_file(tmp_dir):
    """Test that ANBA4 creates the expected XDMF output file."""
    yml_path = tmp_dir / "test.yml"

    # Load YAML to get workdir
    with open(yml_path, "r") as f:
        config = yaml.safe_load(f)
    workdir = config["general"]["workdir"]

    # Check for output file
    output_file = tmp_dir / workdir / "msec_100000_results.xdmf"
    assert output_file.exists(), "ANBA4 should create msec_100000_results.xdmf"
