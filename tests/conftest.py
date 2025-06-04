import pytest
import os
import shutil
from pathlib import Path
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp


@pytest.fixture(scope="session")
def temp_example_dir(tmp_path_factory):
    """Fixture to create a temporary copy of the examples directory, excluding temp_* workdirs."""
    base_dir = Path(__file__).parent.parent  # /home/wr1/projects/b3p/
    example_dir = base_dir / "examples"  # /home/wr1/projects/b3p/examples/
    if not example_dir.exists():
        raise FileNotFoundError(f"Examples directory not found at {example_dir}")

    # Create a base temporary directory
    tmp_base = tmp_path_factory.mktemp("build_examples")
    tmp_dir = tmp_base / "examples"

    # Copy examples dir, excluding temp_* directories
    ignore_patterns = shutil.ignore_patterns("temp_*")
    shutil.copytree(example_dir, tmp_dir, ignore=ignore_patterns)

    # Copy test data to tmp_base/data/
    data_dir = base_dir / "tests" / "data"
    if data_dir.exists():
        shutil.copytree(data_dir, tmp_base / "data", dirs_exist_ok=True)
        print(f"Copied test data from {data_dir} to {tmp_base / 'data'}")
    else:
        print(f"Warning: Test data directory {data_dir} not found")

    return tmp_dir


@pytest.fixture(scope="session")
def built_blade(temp_example_dir):
    """Fixture to build the blade once and provide the resulting workdir."""
    original_dir = os.getcwd()
    os.chdir(temp_example_dir)
    try:
        state = AppState()
        yml_path = Path("blade_test.yml")

        build_app = BuildApp(state, yml_path)
        build_app.build()
        workdir = temp_example_dir / "temp_blade"  # New default workdir
        assert workdir.exists(), f"Blade build failed: workdir {workdir} not created"
        yield {
            "workdir": workdir,
            "temp_dir": temp_example_dir.parent,  # Expose parent temp dir for data access
        }
    finally:
        os.chdir(original_dir)
