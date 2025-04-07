import json
import os
from pathlib import Path

import pytest
import yaml
from b3p.cli2 import AppState, TwoDApp
import shutil


def _prep_temp_dir(temp_dir: Path):
    """Helper function to populate a directory with example files."""
    example_dir = Path(__file__).parent.parent / "examples"
    for item in os.listdir(example_dir):
        s = example_dir / item
        d = temp_dir / item
        if s.is_file():
            shutil.copy2(s, d)
        else:
            shutil.copytree(s, d)


@pytest.fixture(scope="function")
def prep_temp_dir():
    """Provide the helper function for setting up temp dirs."""
    return _prep_temp_dir


@pytest.fixture(scope="function")
def tmp_dir(tmp_path_factory, prep_temp_dir):
    """Create a temporary directory with example files for the test session."""
    tmp_dir = tmp_path_factory.mktemp("test_2d")
    prep_temp_dir(tmp_dir)
    return tmp_dir


@pytest.fixture
def two_d_app():
    """Provide a fresh TwoDApp instance for each test."""
    state = AppState.get_instance()
    return TwoDApp(state)


def test_mesh2d_output(two_d_app, tmp_dir):
    """Test that mesh2d generates a non-None result for test.yml."""
    yml_path = tmp_dir / "blade_test.yml"
    print(f"Testing mesh2d with {yml_path}")
    result = two_d_app.mesh2d(yml_path)
    print(f"mesh2d result: {result}")
    assert result is not None, "mesh2d should return a non-None result"


def test_anba4_output_json(two_d_app, tmp_dir):
    """Test that ANBA4 produces at least one JSON output file with valid content."""
    yml_path = tmp_dir / "blade_test.yml"
    original_dir = os.getcwd()
    os.chdir(tmp_dir)
    try:
        two_d_app.mesh2d(yml_path)  # Generate meshes
        # Assuming TwoDApp has an anba4 method; adjust if named differently
        two_d_app.run_anba4(yml_path)  # Run ANBA4 analysis
        with open(yml_path, "r") as f:
            config = yaml.safe_load(f)
        workdir = config["general"]["workdir"] + "_portable"

        output_jsons = list(Path(tmp_dir / workdir).glob("msec_*.json"))
        assert output_jsons, "ANBA4 should create at least one output JSON file"

        with open(output_jsons[0], "r") as f:
            json_data = json.load(f)
        print(f"Found JSON files: {output_jsons}")
        print(f"First JSON content: {json_data}")
    finally:
        os.chdir(original_dir)


def test_anba4_output_file(two_d_app, tmp_dir):
    """Test that ANBA4 creates the expected XDMF output file."""
    yml_path = tmp_dir / "blade_test.yml"
    original_dir = os.getcwd()
    os.chdir(tmp_dir)
    try:
        two_d_app.mesh2d(yml_path)  # Generate meshes
        # Assuming TwoDApp has an anba4 method; adjust if named differently
        two_d_app.run_anba4(yml_path)  # Run ANBA4 analysis
        with open(yml_path, "r") as f:
            config = yaml.safe_load(f)
        workdir = config["general"]["workdir"] + "_portable"

        output_file = tmp_dir / workdir / "msec_100000.xdmf"
        assert output_file.exists(), "ANBA4 should create msec_100000.xdmf"
    finally:
        os.chdir(original_dir)
