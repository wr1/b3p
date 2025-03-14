import pytest
from pathlib import Path
import shutil
from b3p.cli2 import AppState, TwoDApp
import subprocess
import os
import yaml
import glob
import json

# from b3p.two_d import TwoD
from b3p.utils import copy_example_files


@pytest.fixture(scope="session", autouse=True)
def run_anba4_once(tmp_path_factory):
    """
    Runs run_anba4 once for the entire test session.
    """
    tmp_path = tmp_path_factory.mktemp("anba4_run")
    copy_example_files(tmp_path)
    yml_path = tmp_path / "test.yml"

    # Create an instance of TwoD
    two_d_app = TwoD()

    # Run anba4
    print("Running anba4 once for the test session", tmp_path)
    two_d_app.run_anba4(yml_path)

    # Store tmp_path for later use in tests
    return tmp_path


@pytest.fixture
def two_d_app():
    return TwoDApp()


@pytest.fixture
def tmp_path(run_anba4_once):
    return run_anba4_once


def test_mesh2d_blade_test_yaml(two_d_app, tmp_path):
    print(tmp_path)
    yml_path = tmp_path / "test.yml"
    result = two_d_app.mesh2d(yml_path)
    print(result)
    assert result is not None, "mesh2d should return a non-None result"


def test_anba4_output_json_values(two_d_app, tmp_path):
    yml_path = tmp_path / "test.yml"

    # Load the YAML file to get the workdir
    with open(yml_path, "r") as f:
        config = yaml.safe_load(f)
    workdir = config["general"]["workdir"]

    # Construct the path to the output JSON file
    output_jsons = glob.glob(os.path.join(tmp_path, workdir, "msec_*.json"))

    jsons = [json.load(open(i, "r")) for i in output_jsons]

    print(output_jsons)

    print(jsons[0])

    assert len(output_jsons) > 0, "ANBA4 simulation should create an output JSON file"


def test_run_anba4_valid_config(tmp_path):
    yml_path = tmp_path / "test.yml"
    # Load the YAML file to get the workdir
    with open(yml_path, "r") as f:
        config = yaml.safe_load(f)
    workdir = config["general"]["workdir"]
    # Construct the path to the output file
    output_file = os.path.join(tmp_path, workdir, "msec_100000_results.xdmf")
    # Run anba4 and check if it creates the output file
    assert os.path.exists(output_file), "ANBA4 simulation should create an output file"


def test_mesh2d_blade_test_yaml(two_d_app, tmp_path):
    print(tmp_path)
    yml_path = tmp_path / "test.yml"
    result = two_d_app.mesh2d(yml_path)
    print(result)
    assert result is not None, "mesh2d should return a non-None result"
