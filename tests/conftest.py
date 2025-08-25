import pytest
import os
import shutil
from pathlib import Path
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp
from b3p.cli.two_d_app import TwoDApp
from b3p.cli.ccx_app import CcxApp


@pytest.fixture(scope="session")
def temp_example_dir(tmp_path_factory):
    """Fixture to create a temporary copy of the examples directory, excluding temp_* workdirs."""
    base_dir = Path(__file__).parent.parent
    example_dir = base_dir / "examples"
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


@pytest.fixture(scope="function")
def app_state():
    """Fixture for a clean AppState singleton."""
    state = AppState.get_instance()
    state.reset()
    yield state
    state.reset()


def _run_build(temp_dir):
    original_dir = os.getcwd()
    os.chdir(temp_dir)
    try:
        state = AppState()
        yml_path = Path("blade_test.yml")

        build_app = BuildApp(state, yml_path)
        build_app.build()
        workdir = temp_dir / "temp_blade"  # New default workdir
        assert workdir.exists(), f"Blade build failed: workdir {workdir} not created"
        return {
            "workdir": workdir,
            "temp_dir": temp_dir.parent,  # Expose parent temp dir for data access
        }
    finally:
        os.chdir(original_dir)


@pytest.fixture(scope="session")
def built_blade(temp_example_dir):
    """Fixture to build the blade once and provide the resulting workdir."""
    return _run_build(temp_example_dir)


@pytest.fixture(scope="session")
def meshed_blade(built_blade):
    """Fixture to run the 2D meshing process."""
    original_dir = os.getcwd()
    os.chdir(built_blade["workdir"].parent)
    try:
        state = AppState()
        yml_path = Path("blade_test.yml")

        two_d_app = TwoDApp(state, yml_path)
        two_d_app.mesh2d(rotz=0.0, parallel=False)
        return built_blade
    finally:
        os.chdir(original_dir)


@pytest.fixture(scope="session")
def ccx_analyzed_blade(built_blade):
    """Fixture to run the CCX analysis."""
    original_dir = os.getcwd()
    os.chdir(built_blade["workdir"].parent)
    try:
        state = AppState()
        yml_path = Path("blade_test.yml")

        ccx_app = CcxApp(state, yml_path)
        ccx_app.prep(bondline=False)
        return built_blade
    finally:
        os.chdir(original_dir)
