import pytest
import os
import shutil
from pathlib import Path
import subprocess
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp


@pytest.fixture(scope="session")
def temp_example_dir(tmp_path_factory):
    """Fixture to create a temporary copy of the examples directory and test data."""
    # Use path relative to this conftest file
    base_dir = Path(__file__).parent.parent  # /home/wr1/projects/b3p/
    example_dir = base_dir / "examples"  # /home/wr1/projects/b3p/examples/
    if not example_dir.exists():
        raise FileNotFoundError(
            f"Examples directory not found at {example_dir}. "
            "Ensure 'examples/' exists in /home/wr1/projects/b3p/."
        )

    # Create a base temporary directory
    tmp_base = tmp_path_factory.mktemp("build_examples")
    tmp_dir = tmp_base / "examples"
    shutil.copytree(example_dir, tmp_dir)

    # Copy tests/data/ directory to tmp_base/data/ to match test expectations
    data_dir = base_dir / "tests" / "data"
    if data_dir.exists():
        shutil.copytree(data_dir, tmp_base / "data", dirs_exist_ok=True)
        print(f"Copied test data from {data_dir} to {tmp_base / 'data'}")
    else:
        print(f"Warning: Test data directory {data_dir} not found")

    return tmp_dir


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
