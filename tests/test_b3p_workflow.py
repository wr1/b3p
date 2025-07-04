import pytest
import logging
from pathlib import Path
import shutil
import os
from unittest.mock import patch
from b3p.cli.app_state import AppState
from b3p.cli.clean_app import CleanApp
from b3p.cli.build_app import BuildApp
from b3p.cli.two_d_app import TwoDApp
from b3p.cli.ccx_app import CcxApp
from b3p.cli.utils import check_existing_outputs

# Configure logging for tests
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@pytest.fixture
def temp_workdir(tmp_path):
    """Set up a temporary working directory with blade_test.yml."""
    files = [
        Path("examples/blade_test.yml"),
        Path("examples/materials_iso.yml"),
        Path("examples/loads_v0.yml"),
        Path("examples/laminates_v0.yml"),
    ]
    airfoils = Path("examples/airfoils")

    for src in files:
        if not src.exists():
            pytest.skip(f"{src} not found in examples/")
        dest = tmp_path / src.name
        shutil.copy(src, dest)
    shutil.copytree(airfoils, tmp_path / "airfoils", dirs_exist_ok=True)

    # Set CWD to tmp_path
    os.chdir(tmp_path)

    yield tmp_path, Path("blade_test.yml")
    # Cleanup
    state = AppState.get_instance()
    state.reset()


@pytest.fixture
def app_state():
    """Fixture for AppState singleton."""
    state = AppState.get_instance()
    yield state
    state.reset()


def test_clean(temp_workdir, app_state, caplog):
    """Test cleaning the working directory."""
    tmp_path, yml_path = temp_workdir
    caplog.set_level(logging.INFO)

    # Create a dummy workdir
    workdir = tmp_path / "temp_blade"
    workdir.mkdir()
    (workdir / "dummy.txt").write_text("test")

    clean_app = CleanApp(app_state, yml_path)
    clean_app.clean()

    assert not workdir.exists()
    assert "Removed workdir" in caplog.text


def test_build(temp_workdir, app_state, caplog):
    """Test the build process, skipping if outputs exist."""
    tmp_path, yml_path = temp_workdir
    caplog.set_level(logging.INFO)

    build_app = BuildApp(app_state, yml_path)

    # Check if key outputs already exist
    prefix = app_state.load_yaml(yml_path)["general"].get("prefix", "b3p")
    workdir = app_state.get_workdir()
    expected_files = [
        workdir / "drape" / f"{prefix}_joined.vtu",
        workdir / "drape" / f"{prefix}_mass.csv",
        workdir / "drape" / f"{prefix}_loads.png",
    ]

    if check_existing_outputs(expected_files):
        logger.info("Build outputs already exist, skipping build")
    else:
        build_app.build(bondline=True)

    for f in expected_files:
        assert f.exists(), f"Expected file {f} not found"
    assert "Mass table per material" in caplog.text


@pytest.mark.skip("seems a double test, output already checked in test_anba")
def test_2d_analysis(temp_workdir, app_state, caplog):
    """Test 2D analysis, mocking ANBA4 subprocess."""
    tmp_path, yml_path = temp_workdir
    caplog.set_level(logging.INFO)

    two_d_app = TwoDApp(app_state, yml_path)

    # Ensure build outputs exist
    build_app = BuildApp(app_state, yml_path)
    prefix = app_state.load_yaml(yml_path)["general"].get("prefix", "b3p")
    workdir = app_state.get_workdir()
    joined_vtu = workdir / "drape" / f"{prefix}_joined.vtu"

    if not joined_vtu.exists():
        build_app.build(bondline=True)

    # Mock subprocess.run to avoid running ANBA4
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        two_d_app.mesh2d(rotz=0.0, parallel=False)
        two_d_app.run_anba4(anba_env="anba4-env")

    expected_files = [
        workdir / "2d" / "msec_1000.xdmf",
        workdir / "2d" / "msec_80000.xdmf",
    ]
    for f in expected_files:
        assert f.exists(), f"Expected file {f} not found"
    assert "ANBA4 script completed successfully" in caplog.text


def test_ccx_analysis(temp_workdir, app_state, caplog):
    """Test CCX analysis, mocking subprocess."""
    tmp_path, yml_path = temp_workdir
    caplog.set_level(logging.INFO)

    ccx_app = CcxApp(app_state, yml_path)

    # Ensure build outputs exist
    build_app = BuildApp(app_state, yml_path)
    prefix = app_state.load_yaml(yml_path)["general"].get("prefix", "b3p")
    workdir = app_state.get_workdir()
    joined_vtu = workdir / "drape" / f"{prefix}_joined.vtu"

    if not joined_vtu.exists():
        build_app.build(bondline=True)

    # Mock subprocess.run to avoid running CCX
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        ccx_app.ccx(bondline=True)

    expected_files = [
        workdir / "fea" / f"{prefix}_ccx_lc_forward_flap.inp",
    ]
    for f in expected_files:
        assert f.exists(), f"Expected file {f} not found"
    assert f"Written: {workdir}/fea" in caplog.text


def test_full_workflow(temp_workdir, app_state, caplog):
    """Test the full workflow in sequence."""
    tmp_path, yml_path = temp_workdir
    caplog.set_level(logging.INFO)

    # Clean
    clean_app = CleanApp(app_state, yml_path)
    clean_app.clean()

    # Build
    build_app = BuildApp(app_state, yml_path)
    prefix = app_state.load_yaml(yml_path)["general"].get("prefix", "b3p")
    workdir = app_state.get_workdir()
    build_app.build(bondline=True)

    # 2D Analysis
    two_d_app = TwoDApp(app_state, yml_path)
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        two_d_app.mesh2d(rotz=0.0, parallel=False)
        two_d_app.run_anba4(anba_env="anba4-env")

    # CCX Analysis
    ccx_app = CcxApp(app_state, yml_path)
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        ccx_app.ccx(bondline=False)

    # Verify key outputs
    expected_files = [
        workdir / "drape" / f"{prefix}_joined.vtu",
        workdir / "2d" / "msec_1000.xdmf",
        workdir / "fea" / f"{prefix}_ccx_lc_forward_flap.inp",
    ]
    for f in expected_files:
        assert f.exists(), f"Expected file {f} not found"
