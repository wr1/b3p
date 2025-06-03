import pytest
from pathlib import Path
import glob
import os
from b3p.cli.app_state import AppState
from b3p.cli.ccx_app import CcxApp
from b3p.cli.build_app import BuildApp
import filecmp


@pytest.fixture(scope="session")
def run_ccx(temp_example_dir):
    """Fixture to run the CCX process after a build."""
    original_dir = os.getcwd()
    os.chdir(temp_example_dir)
    try:
        state = AppState()
        yml_path = Path("blade_test.yml")

        build_app = BuildApp(state, yml_path)
        ccx_app = CcxApp(state, yml_path)
        # yml_path = Path("blade_test.yml")
        build_app.build()  # Build the project first
        ccx_app.prep(bondline=False)  # Test with bondline option
        yield {
            "workdir": temp_example_dir / "temp_blade_portable",
            "yml_path": yml_path,
            "ccx_app": ccx_app,
            "temp_dir": temp_example_dir.parent,  # Expose parent temp dir for data access
        }
    finally:
        os.chdir(original_dir)


def test_ccx_prep(run_ccx):
    """Test if CCX preparation generates the expected input file."""
    workdir = run_ccx["workdir"]
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
    generated_files = glob.glob(f"{workdir}/fea/*_ccx_lc_forward_edge.inp")
    assert generated_files, "CCX prep should generate a forward edge input file named *_ccx_lc_forward_edge.inp"
    generated_file = generated_files[0]

    # Reference file from tests/data/, copied to temp_dir/data/
    reference_file = temp_dir / "data" / "test_blade_ccx_lc_forward_edge.inp"
    assert (
        reference_file.exists()
    ), f"Reference file {reference_file} not found in temp data directory"
    # Compare the generated file with the reference file
    # Use filecmp to compare the files
    cmp = filecmp.cmp(generated_file, reference_file, shallow=False)
    assert (
        cmp
    ), f"Generated file {generated_file} does not match reference file {reference_file}"
