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


def test_clean(temp_example_dir, app_state, caplog):
    """Test cleaning the working directory."""
    yml_path = temp_example_dir / "blade_test.yml"
    caplog.set_level(logging.INFO)

    # Create a dummy workdir
    workdir = temp_example_dir / "temp_blade"
    workdir.mkdir(exist_ok=True)
    (workdir / "dummy.txt").write_text("test")

    clean_app = CleanApp(app_state, yml_path)
    clean_app.clean()

    assert not workdir.exists()
    assert "Removed workdir" in caplog.text


def test_build(built_blade, app_state, caplog):
    """Test the build process, skipping if outputs exist."""
    workdir = built_blade["workdir"]
    caplog.set_level(logging.INFO)

    # Check if key outputs already exist
    app_state.load_yaml(workdir.parent / "blade_test.yml")
    prefix = app_state.config.general.prefix
    expected_files = [
        workdir / "drape" / f"{prefix}_joined.vtu",
        workdir / "drape" / f"{prefix}_mass.csv",
        workdir / "drape" / f"{prefix}_loads.png",
    ]

    if check_existing_outputs(expected_files):
        logger.info("Build outputs already exist, skipping build")
    else:
        build_app = BuildApp(app_state, workdir.parent / "blade_test.yml")
        build_app.build(bondline=True)

    for f in expected_files:
        assert f.exists(), f"Expected file {f} not found"
    assert "Mass table per material" in caplog.text


@pytest.mark.skip("seems a double test, output already checked in test_anba")
def test_2d_analysis(meshed_blade, app_state, caplog):
    """Test 2D analysis, mocking ANBA4 subprocess."""
    workdir = meshed_blade["workdir"]
    yml_path = workdir.parent / "blade_test.yml"
    caplog.set_level(logging.INFO)

    two_d_app = TwoDApp(app_state, yml_path)

    # Mock subprocess.run to avoid running ANBA4
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        two_d_app.run_anba4(anba_env="anba4-env")

    app_state.load_yaml(yml_path)
    prefix = app_state.config.general.prefix
    expected_files = [
        workdir / "2d" / "msec_1000.xdmf",
        workdir / "2d" / "msec_80000.xdmf",
    ]
    for f in expected_files:
        assert f.exists(), f"Expected file {f} not found"
    assert "ANBA4 script completed successfully" in caplog.text


def test_ccx_analysis(ccx_analyzed_blade, app_state, caplog):
    """Test CCX analysis, mocking subprocess."""
    workdir = ccx_analyzed_blade["workdir"]
    yml_path = workdir.parent / "blade_test.yml"
    caplog.set_level(logging.INFO)

    # Mock subprocess.run to avoid running CCX
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        ccx_app = CcxApp(app_state, yml_path)
        ccx_app.ccx(bondline=True)

    app_state.load_yaml(yml_path)
    prefix = app_state.config.general.prefix
    expected_files = [
        workdir / "fea" / f"{prefix}_ccx_lc_forward_flap.inp",
    ]
    for f in expected_files:
        assert f.exists(), f"Expected file {f} not found"
    assert f"{workdir}/fea/{prefix}_ccx_lc_forward_flap.inp" in caplog.text


def test_full_workflow(temp_example_dir, app_state, caplog):
    """Test the full workflow in sequence."""
    yml_path = temp_example_dir / "blade_test.yml"
    caplog.set_level(logging.INFO)

    # Clean
    clean_app = CleanApp(app_state, yml_path)
    clean_app.clean()

    # Build
    build_app = BuildApp(app_state, yml_path)
    app_state.load_yaml(yml_path)
    prefix = app_state.config.general.prefix
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
    workdir / "drape" / "2d" / "msec_1000.xdmf",
    for f in expected_files:
        assert f.exists(), f"Expected file {f} not found"
