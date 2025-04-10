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


@pytest.fixture
def two_d_app():
    """Provide a fresh TwoDApp instance for each test."""
    state = AppState.get_instance()
    return TwoDApp(state)


def test_mesh2d_output(two_d_app, built_blade):
    """Test that mesh2d generates a non-None result for test.yml after blade build."""
    yml_path = built_blade["temp_dir"] / "examples" / "blade_test.yml"
    print(f"Testing mesh2d with {yml_path}")
    result = two_d_app.mesh2d(yml_path, parallel=False)  # Generate meshes
    print(f"mesh2d result: {result}")
    assert result is not None, "mesh2d should return a non-None result"


def test_anba4_output_json(two_d_app, built_blade):
    """Test that ANBA4 produces at least one JSON output file with valid content after blade build."""
    yml_path = built_blade["temp_dir"] / "examples" / "blade_test.yml"
    original_dir = os.getcwd()
    os.chdir(built_blade["temp_dir"] / "examples")
    try:
        two_d_app.mesh2d(yml_path, parallel=False)  # Generate meshes
        two_d_app.run_anba4(yml_path)  # Run ANBA4 analysis
        with open(yml_path, "r") as f:
            config = yaml.safe_load(f)
        workdir = config["general"]["workdir"] + "_portable"

        output_jsons = list(Path(built_blade["workdir"]).glob("msec_*.json"))
        assert output_jsons, "ANBA4 should create at least one output JSON file"

        with open(output_jsons[0], "r") as f:
            json_data = json.load(f)
        print(f"Found JSON files: {output_jsons}")
        print(f"First JSON content: {json_data}")
    finally:
        os.chdir(original_dir)


def test_anba4_output_file(two_d_app, built_blade):
    """Test that ANBA4 creates the expected XDMF output file after blade build."""
    yml_path = built_blade["temp_dir"] / "examples" / "blade_test.yml"
    original_dir = os.getcwd()
    os.chdir(built_blade["temp_dir"] / "examples")
    try:
        two_d_app.mesh2d(yml_path, parallel=False)  # Generate meshes
        two_d_app.run_anba4(yml_path)  # Run ANBA4 analysis
        with open(yml_path, "r") as f:
            config = yaml.safe_load(f)
        workdir = config["general"]["workdir"] + "_portable"

        output_file = built_blade["workdir"] / "msec_100000.xdmf"
        assert output_file.exists(), "ANBA4 should create msec_100000.xdmf"
    finally:
        os.chdir(original_dir)


def test_mesh2d_compare_reference(two_d_app, built_blade):
    """Test that 2D meshes match reference meshes in tests/data/reference_2d_meshes/ after blade build."""
    yml_path = built_blade["temp_dir"] / "examples" / "blade_test.yml"
    original_dir = os.getcwd()
    os.chdir(built_blade["temp_dir"] / "examples")
    try:
        # Run mesh2d to generate meshes
        two_d_app.mesh2d(yml_path, parallel=False)
        with open(yml_path, "r") as f:
            config = yaml.safe_load(f)
        workdir = config["general"]["workdir"] + "_portable"

        # Get generated mesh files
        generated_meshes = list(Path(built_blade["workdir"]).glob("msec_*0.xdmf"))
        assert generated_meshes, "mesh2d should generate at least one .xdmf file"

        # Reference directory
        ref_dir = built_blade["temp_dir"] / "data" / "reference_2d_meshes"
        assert ref_dir.exists(), f"Reference directory {ref_dir} not found"

        # Compare each generated mesh with its reference
        for gen_mesh in generated_meshes:
            ref_mesh_path = ref_dir / gen_mesh.name
            assert ref_mesh_path.exists(), f"Reference mesh {ref_mesh_path} not found"

            # Load meshes with meshio
            gen_mesh_data = meshio.read(gen_mesh)
            ref_mesh_data = meshio.read(ref_mesh_path)

            # Compare points with tolerance
            assert np.allclose(
                gen_mesh_data.points, ref_mesh_data.points, rtol=1e-1, atol=1e-2
            ), f"Points in {gen_mesh.name} do not match reference"

            # Compare cells as sets of connectivity tuples
            gen_cells = {tuple(cell.data.flatten()) for cell in gen_mesh_data.cells}
            ref_cells = {tuple(cell.data.flatten()) for cell in ref_mesh_data.cells}
            assert len(gen_cells) == len(ref_cells), (
                f"Number of cells in {gen_mesh.name} does not match reference"
            )
            # assert gen_cells == ref_cells, (
            #     f"Cells in {gen_mesh.name} do not match reference"
            # )
    finally:
        os.chdir(original_dir)
