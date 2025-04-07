# tests/conftest.py
import pytest
import os
import shutil
from pathlib import Path
import subprocess
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp


@pytest.fixture(scope="session")
def temp_example_dir(tmp_path_factory):
    example_dir = Path("examples")
    if not example_dir.exists():
        pytest.skip("Examples directory not found; skipping build tests.")
    tmp_dir = tmp_path_factory.mktemp("build_examples")
    shutil.copytree(example_dir, tmp_dir / "examples")
    return tmp_dir / "examples"


@pytest.fixture(scope="session")
def run_build(temp_example_dir):
    original_dir = os.getcwd()
    os.chdir(temp_example_dir)
    try:
        state = AppState()
        build_app = BuildApp(state)
        build_app.build(Path("blade_test.yml"))
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
