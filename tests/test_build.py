import pytest
import os
import shutil
import glob
from pathlib import Path
import tempfile
import subprocess
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp


@pytest.fixture(scope="session")
def temp_example_dir():
    """Fixture to create a temporary copy of the examples directory."""
    example_dir = Path("examples")
    if not example_dir.exists():
        pytest.skip("Examples directory not found; skipping build tests.")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_example_dir = Path(tmpdir) / "examples"
        shutil.copytree(example_dir, tmp_example_dir)
        yield tmp_example_dir


@pytest.fixture(scope="session")
def run_build(temp_example_dir):
    """Fixture to run the build process and capture output."""
    original_dir = os.getcwd()
    os.chdir(temp_example_dir)
    try:
        # Run the build using BuildApp directly
        state = AppState()
        build_app = BuildApp(state)
        build_app.build(Path("blade_test.yml"))

        # Capture stdout by running the CLI command
        result = subprocess.run(
            ["python", "-m", "b3p.cli2", "build", "blade_test.yml"],
            capture_output=True,
            text=True,
        )
        yield {
            "workdir": temp_example_dir / "temp_blade_portable",
            "stdout": result.stdout,
            "stderr": result.stderr,
            "returncode": result.returncode,
        }
    finally:
        os.chdir(original_dir)


def test_build_success(run_build):
    """Test if the build completes successfully."""
    assert run_build["returncode"] == 0, (
        f"Build failed with stderr: {run_build['stderr']}"
    )
    assert run_build["workdir"].exists(), "Working directory was not created."


def test_build_output_files(run_build):
    """Test if all expected output files are generated."""
    workdir = run_build["workdir"]
    output_files = {os.path.basename(f) for f in glob.glob(str(workdir / "*"))}

    expected_files = {
        "test_blade_base.vtp",
        "test_blade.png",
        "test_blade_shell.vtp",
        "test_blade_w0.vtp",
        "test_blade_w1.vtp",
        "test_blade_w2.vtp",
        "test_blade_w3.vtp",
        "test_blade_w4.vtp",
        "test_blade_w3_dr.vtu",
        "test_blade_shell_dr.vtu",
        "test_blade_w4_dr.vtu",
        "__matdb.yml",
        "test_blade_joined.vtu",
        "test_blade_portable.yml",
        "test_blade.var",
        "material_map.json",
        "test_blade_loads.png",
        "__plybook.pck",
        "test_blade.vtp",
        "test_blade_joined_bondline.vtu",
        "test_blade_mass.csv",
        "test_blade_mass.tex",
    }

    missing_files = expected_files - output_files
    assert not missing_files, f"Missing expected files: {missing_files}"


def test_build_stdout(run_build):
    """Test if key stdout messages are present."""
    stdout = run_build["stdout"]

    print(stdout)
    expected_messages = [
        "** Loading yaml file blade_test.yml and loading linked files",
        "saving to temp_blade_portable/test_blade.vtp",
        "written to: temp_blade_portable/test_blade_portable.yml",
        "** wrote mesh to temp_blade_portable/test_blade_base.vtp",
        "** wrote mesh to temp_blade_portable/test_blade_shell.vtp",
        "written material map to temp_blade_portable/material_map.json",
        "written to temp_blade_portable/__plybook.pck",
        "written mesh to temp_blade_portable/test_blade_joined.vtu",
        "Saved temp_blade_portable/test_blade_joined_bondline.vtu",
        "Loaded material database from temp_blade_portable/__matdb.yml",
        "Total Volume and Mass:",
        "Mass table per material",
        "** written load plot to temp_blade_portable/test_blade_loads.png",
    ]

    for message in expected_messages:
        assert message in stdout, f"Expected message not found in stdout: {message}"


def test_build_mass_table(run_build):
    """Test if the mass table in stdout matches expected values."""
    stdout = run_build["stdout"]

    mass_table_lines = [
        "11457.585028",
        "937.880619",
        "10047.614323",
        "31119.325388",
        "1125.32413",
        "783.68881",
        "3789.491169",
        "1570.496484",
        "60831.405951",
    ]

    for line in mass_table_lines:
        assert line in stdout, f"Expected mass table entry not found: {line}"


def test_build_geometry_only(temp_example_dir):
    """Test the geometry subcommand independently."""
    original_dir = os.getcwd()
    os.chdir(temp_example_dir)
    try:
        state = AppState()
        build_app = BuildApp(state)
        build_app.geometry(Path("blade_test.yml"))

        workdir = temp_example_dir / "temp_blade_portable"
        assert workdir.exists(), "Working directory not created."
        assert (workdir / "test_blade.vtp").exists(), "Geometry file not created."
        assert (workdir / "test_blade_portable.yml").exists(), (
            "Portable YAML not created."
        )
    finally:
        os.chdir(original_dir)
