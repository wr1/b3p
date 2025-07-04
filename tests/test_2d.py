import json
from pathlib import Path

import pytest
import numpy as np
from b3p.cli.app_state import AppState
from b3p.cli.two_d_app import TwoDApp
import meshio


@pytest.fixture
def two_d_app(built_blade):
    """Provide a fresh TwoDApp instance for each test."""
    state = AppState.get_instance()
    yml_path = built_blade["temp_dir"] / "examples" / "blade_test.yml"
    return TwoDApp(state, yml_path)


def test_mesh2d_output(two_d_app):
    """Test that mesh2d generates a non-None result for test.yml after blade build."""
    result = two_d_app.mesh2d(parallel=False)  # Generate meshes
    print(f"mesh2d result: {result}")
    assert result is not None, "mesh2d should return a non-None result"


def test_anba4_output_json(two_d_app):
    """Test that ANBA4 produces at least one JSON output file with valid content after blade build."""
    try:
        two_d_app.mesh2d(parallel=False)  # Generate meshes
        two_d_app.run_anba4()  # Run ANBA4 analysis

        prefix = two_d_app.state.get_prefix("drape")

        output_jsons = list(prefix.parent.glob("2d/msec_*.json"))
        assert output_jsons, "ANBA4 should create at least one output JSON file"

        with open(output_jsons[0], "r") as f:
            json_data = json.load(f)
        print(f"Found JSON files: {output_jsons}")
        print(f"First JSON content: {json_data}")
    finally:
        pass


def test_anba4_output_file(two_d_app):
    """Test that ANBA4 creates the expected XDMF output file after blade build."""
    try:
        two_d_app.mesh2d(parallel=False)  # Generate meshes
        two_d_app.run_anba4()  # Run ANBA4 analysis
        prefix = two_d_app.state.get_prefix("drape")

        output_file = Path(prefix).parent / "2d" / "msec_100000.xdmf"
        assert output_file.exists(), "ANBA4 should create msec_100000.xdmf"
    finally:
        pass


def test_mesh2d_compare_reference(two_d_app, built_blade):
    """Test that 2D meshes match reference meshes in tests/data/reference_2d_meshes/ after blade build."""
    try:
        # Run mesh2d to generate meshes
        two_d_app.mesh2d(parallel=False)

        prefix = two_d_app.state.get_prefix("drape")

        print(f"Using prefix: {prefix}")
        # Get generated mesh files
        generated_meshes = list(Path(prefix.parent).glob("2d/msec_*0.xdmf"))

        print(f"Generated meshes: {generated_meshes}")
        assert generated_meshes, "mesh2d should generate at least one .xdmf file"

        # Reference directory in temporary test data
        ref_dir = built_blade["temp_dir"] / "data" / "reference_2d_meshes"
        assert ref_dir.exists(), f"Reference directory {ref_dir} not found"

        # Compare each generated mesh with its reference
        for gen_mesh in generated_meshes:
            ref_mesh_path = ref_dir / gen_mesh.name
            assert ref_mesh_path.exists(), f"Reference mesh {ref_mesh_path} not found"

            # Load meshes with pyvista
            gen_mesh_data = meshio.read(gen_mesh)
            ref_mesh_data = meshio.read(ref_mesh_path)

            # Compare mesh bounding boxes
            gen_bounds = gen_mesh_data.points.min()
            ref_bounds = ref_mesh_data.points.min()
            assert np.allclose(
                gen_bounds, ref_bounds, rtol=1e-1, atol=1e-2
            ), f"Bounding boxes in {gen_mesh.name} do not match reference"

            # Compare number of points
            assert len(gen_mesh_data.points) == len(
                ref_mesh_data.points
            ), f"Number of points in {gen_mesh.name} does not match reference"
    finally:
        pass
