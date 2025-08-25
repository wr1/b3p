import pytest
import glob
import os
import filecmp
import logging

logger = logging.getLogger(__name__)


def test_ccx_prep(ccx_analyzed_blade):
    """Test if CCX preparation generates the expected input file."""
    workdir = ccx_analyzed_blade["workdir"]

    logger.info(f"Checking CCX prep in workdir: {workdir}")
    inp_files = glob.glob(f"{workdir}/fea/*_ccx_*.inp")
    assert inp_files, "CCX prep should generate at least one .inp file"
    assert os.path.exists(inp_files[0]), (
        f"Expected CCX input file {inp_files[0]} not found"
    )


def test_ccx_bondline_selection(ccx_analyzed_blade):
    """Test if bondline mesh is selected when bondline=True."""
    workdir = ccx_analyzed_blade["workdir"]
    bondline_vtu = glob.glob(f"{workdir}/drape/*_bondline.vtu")
    assert bondline_vtu, "Bondline VTU should exist from the build process"
    inp_files = glob.glob(f"{workdir}/fea/*_ccx_*.inp")
    assert inp_files[0].startswith(str(workdir / "fea" / "test_blade")), (
        "CCX input file should be generated from bondline mesh"
    )


def test_ccx_produce_fwd_edge_inp(ccx_analyzed_blade):
    """Test if CCX produces the forward edge input file."""
    workdir = ccx_analyzed_blade["workdir"]
    inp_files = glob.glob(f"{workdir}/fea/*_ccx_*.inp")
    assert inp_files, "CCX prep should generate at least one .inp file"
    assert any("_forward_edge" in f for f in inp_files), (
        "CCX prep should produce a forward edge input file"
    )


@pytest.mark.skip(
    reason="Seems a element numbering issue on the web, not sure if this gives different results, skip for now"
)
def test_ccx_forward_edge_content(ccx_analyzed_blade):
    """Test if the generated forward edge .inp file matches the reference file."""
    workdir = ccx_analyzed_blade["workdir"]
    temp_dir = ccx_analyzed_blade["temp_dir"]
    edgewise_loadcase = workdir / "fea" / "test_blade_ccx_lc_forward_edge.inp"
    assert edgewise_loadcase.exists(), (
        "CCX prep should generate a forward edge input file named *_ccx_lc_forward_edge.inp"
    )
    # generated_file = generated_files[0]

    reference_file = temp_dir / "data" / "test_blade_ccx_lc_forward_edge.inp"
    assert reference_file.exists(), (
        f"Reference file {reference_file} not found in temp data directory"
    )

    cmp = filecmp.cmp(edgewise_loadcase, reference_file, shallow=False)
    assert cmp, (
        f"Generated file {edgewise_loadcase} does not match reference file {reference_file}"
    )
