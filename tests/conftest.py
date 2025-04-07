import pytest
import os
import shutil
from pathlib import Path
import subprocess
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp


@pytest.fixture(scope="session")
def temp_example_dir(tmp_path_factory):
    """Fixture to create a temporary copy of the examples directory."""
    example_dir = Path("examples")
    if not example_dir.exists():
        pytest.skip("Examples directory not found; skipping build tests.")
    tmp_dir = tmp_path_factory.mktemp("build_examples")
    shutil.copytree(example_dir, tmp_dir / "examples")

    # Copy tests/data/ directory if it exists
    data_dir = Path("tests/data")
    if data_dir.exists():
        shutil.copytree(data_dir, tmp_dir / "data")

    return tmp_dir / "examples"


@pytest.fixture(scope="session")
def run_build(temp_example_dir):
    """Fixture to run the build process and capture output."""
    original_dir = os.getcwd()
    os.chdir(temp_example_dir)
    try:
        state = AppState()
        build_app = BuildApp(state)
        build_app.build(Path("blade_test.yml"))
        result = subprocess.run(
            ["python", "-m", "b3p.__main__", "build", "blade_test.yml"],
            capture_output=True,
            text=True,
        )
        yield {
            "workdir": temp_example_dir / "temp_blade_portable",
            "stdout": result.stdout,
            "stderr": result.stderr,
            "returncode": result.returncode,
            "temp_dir": temp_example_dir.parent,  # Expose parent temp dir for data access
        }
    finally:
        os.chdir(original_dir)
