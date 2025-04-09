import json
import os
from pathlib import Path

import pytest
import yaml
import meshio
import numpy as np
from b3p.cli.app_state import AppState
from b3p.cli.two_d_app import TwoDApp
import shutil


def _prep_temp_dir(temp_dir: Path):
    """Helper function to populate a directory with example files and test data."""
    base_dir = Path(__file__).parent.parent
    example_dir = base_dir / "examples"
    for item in os.listdir(example_dir):
        s = example_dir / item
        d = temp_dir / item
        if s.is_file():
            shutil.copy2(s, d)
        else:
            shutil.copytree(s, d)
    # Copy reference data directly into temp_dir
    data_dir = base_dir / "tests" / "data"
    if data_dir.exists():
        shutil.copytree(data_dir, temp_dir / "data", dirs_exist_ok=True)
        print(f"Copied test data from {data_dir} to {temp_dir / 'data'}")
    else:
        print(f"Warning: Test data directory {data_dir} not found")


@pytest.fixture(scope="function")
def prep_temp_dir():
    """Provide the helper function for setting up temp dirs."""
    return _prep_temp_dir


@pytest.fixture(scope="function")
def tmp_dir(tmp_path_factory, prep_temp_dir):
    """Create a temporary directory with example files and test data for the test."""
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
    result = two_d_app.mesh2d(yml_path, parallel=False)  # Generate meshes
    print(f"mesh2d result: {result}")
    assert result is not None, "mesh2d should return a non-None result"


def test_anba4_output_json(two_d_app, tmp_dir):
    """Test that ANBA4 produces at least one JSON output file with valid content."""
    yml_path = tmp_dir / "blade_test.yml"
    original_dir = os.getcwd()
    os.chdir(tmp_dir)
    try:
        two_d_app.mesh2d(yml_path, parallel=False)  # Generate meshes
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
        two_d_app.mesh2d(yml_path, parallel=False)  # Generate meshes
        two_d_app.run_anba4(yml_path)  # Run ANBA4 analysis
        with open(yml_path, "r") as f:
            config = yaml.safe_load(f)
        workdir = config["general"]["workdir"] + "_portable"

        output_file = tmp_dir / workdir / "msec_100000.xdmf"
        assert output_file.exists(), "ANBA4 should create msec_100000.xdmf"
    finally:
        os.chdir(original_dir)


def test_mesh2d_compare_reference(two_d_app, tmp_dir):
    """Test that 2D meshes match reference meshes in tests/data/reference_2d_meshes/."""
    yml_path = tmp_dir / "blade_test.yml"
    original_dir = os.getcwd()
    os.chdir(tmp_dir)
    try:
        # Run mesh2d to generate meshes
        two_d_app.mesh2d(yml_path, parallel=False)
        with open(yml_path, "r") as f:
            config = yaml.safe_load(f)
        workdir = config["general"]["workdir"] + "_portable"

        # Get generated mesh files
        generated_meshes = list(Path(tmp_dir / workdir).glob("msec_*0.xdmf"))
        assert generated_meshes, "mesh2d should generate at least one .xdmf file"

        # Reference directory (now directly in tmp_dir/data/)
        ref_dir = tmp_dir / "data" / "reference_2d_meshes"
        assert ref_dir.exists(), f"Reference directory {ref_dir} not found"

        # Compare each generated mesh with its reference
        for gen_mesh in generated_meshes:
            ref_mesh_path = ref_dir / gen_mesh.name
            assert ref_mesh_path.exists(), f"Reference mesh {ref_mesh_path} not found"

            # Load meshes with meshio
            gen_mesh_data = meshio.read(gen_mesh)
            ref_mesh_data = meshio.read(ref_mesh_path)

            # Compare points and cells
            assert np.allclose(gen_mesh_data.points, ref_mesh_data.points), (
                f"Points in {gen_mesh.name} do not match reference"
            )
            assert len(gen_mesh_data.cells) == len(ref_mesh_data.cells), (
                f"Number of cells in {gen_mesh.name} do not match reference"
            )
            for gen_cells, ref_cells in zip(gen_mesh_data.cells, ref_mesh_data.cells):
                assert gen_cells.type == ref_cells.type, (
                    f"Cell type in {gen_mesh.name} does not match reference"
                )
                print(
                    "test",
                    gen_mesh.name,
                    np.allclose(gen_cells.data, ref_cells.data, rtol=1e-5, atol=1e-3),
                    gen_cells.data,
                    ref_cells.data,
                )
                print("test1", (gen_cells.data - ref_cells.data).sum())
                assert np.allclose(gen_cells.data, ref_cells.data), (
                    f"Cells in {gen_mesh.name} do not match reference"
                )
    finally:
        os.chdir(original_dir)
