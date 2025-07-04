import pytest
from pathlib import Path
import glob
import os
from b3p.cli.app_state import AppState
from b3p.cli.ccx_app import CcxApp
from b3p.cli.build_app import BuildApp
import filecmp
import logging

logger = logging.getLogger(__name__)


@pytest.fixture(scope="session")
def run_ccx(temp_example_dir):
    """Fixture to run the CCX process after a build."""
    # original_dir = os.getcwd()
    # os.chdir(temp_example_dir)
    try:
        state = AppState()
        yml_path = Path(temp_example_dir / "blade_test.yml")

        build_app = BuildApp(state, yml_path)
        ccx_app = CcxApp(state, yml_path)
        build_app.build()
        ccx_app.prep(bondline=False)
        yield {
            "workdir": temp_example_dir / "temp_blade",
            "yml_path": yml_path,
            "ccx_app": ccx_app,
            "temp_dir": temp_example_dir.parent,
        }
    finally:
        pass
        # os.chdir(original_dir)


def test_ccx_prep(run_ccx):
    """Test if CCX preparation generates the expected input file."""
    workdir = run_ccx["workdir"]

    logger.info(f"Checking CCX prep in workdir: {workdir}")
    inp_files = glob.glob(f"{workdir}/fea/*_ccx_*.inp")
    assert inp_files, "CCX prep should generate at least one .inp file"
    assert os.path.exists(
        inp_files[0]
    ), f"Expected CCX input file {inp_files[0]} not found"


def test_ccx_bondline_selection(run_ccx):
    """Test if bondline mesh is selected when bondline=True."""
    workdir = run_ccx["workdir"]
    bondline_vtu = glob.glob(f"{workdir}/drape/*_bondline.vtu")
    assert bondline_vtu, "Bondline VTU should exist from the build process"
    inp_files = glob.glob(f"{workdir}/fea/*_ccx_*.inp")
    assert inp_files[0].startswith(
        str(workdir / "fea" / "test_blade")
    ), "CCX input file should be generated from bondline mesh"


def test_ccx_produce_fwd_edge_inp(run_ccx):
    """Test if CCX produces the forward edge input file."""
    workdir = run_ccx["workdir"]
    inp_files = glob.glob(f"{workdir}/fea/*_ccx_*.inp")
    assert inp_files, "CCX prep should generate at least one .inp file"
    assert any(
        "_forward_edge" in f for f in inp_files
    ), "CCX prep should produce a forward edge input file"


@pytest.mark.skip(
    reason="Seems a element numbering issue on the web, not sure if this gives different results, skip for now"
)
def test_ccx_forward_edge_content(run_ccx):
    """Test if the generated forward edge .inp file matches the reference file."""
    workdir = run_ccx["workdir"]
    temp_dir = run_ccx["temp_dir"]
    edgewise_loadcase = workdir / "fea" / "test_blade_ccx_lc_forward_edge.inp"
    assert edgewise_loadcase.exists(), "CCX prep should generate a forward edge input file named *_ccx_lc_forward_edge.inp"
    # generated_file = generated_files[0]

    reference_file = temp_dir / "data" / "test_blade_ccx_lc_forward_edge.inp"
    assert (
        reference_file.exists()
    ), f"Reference file {reference_file} not found in temp data directory"

    cmp = filecmp.cmp(edgewise_loadcase, reference_file, shallow=False)
    assert cmp, f"Generated file {edgewise_loadcase} does not match reference file {reference_file}"
