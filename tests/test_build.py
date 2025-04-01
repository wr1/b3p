import pytest
import glob
import os
from b3p.cli2 import AppState, BuildApp
from pathlib import Path

example_dir = Path("examples")
workdir = example_dir / "temp_blade_portable"


def build_output_exists():
    """Check if the output of the build exists."""
    print(workdir)
    return os.path.exists(workdir)


def run_tbuild():
    """Run step s1."""
    if not build_output_exists():
        cd = os.getcwd()
        os.chdir(example_dir)
        BuildApp(AppState()).build("blade_test.yml")
        # b3p_cli.build("blade_test.yml")
        os.chdir(cd)
    else:
        print("Output already exists.")


@pytest.fixture(scope="session")
def run_test_build():
    """Fixture to run blade build."""
    run_tbuild()


def test_build_success(run_test_build):
    """Test if the build was successful."""
    assert build_output_exists()


def test_build_has_all_files(run_test_build):
    """Test if the build has all the expected files."""
    expected_files = """test_blade_base.vtp test_blade.png test_blade_shell.vtp test_blade_w1.vtp         test_blade_w3_dr.vtu test_blade_w4_dr.vtu __matdb.yml test_blade_joined.vtu test_blade_portable.yml  test_blade.var test_blade_w0.vtp material_map.json test_blade_loads.png test_blade_sca_50.csv test_blade.vtp __plybook.pck  test_blade.pck test_blade_shell_dr.vtu  test_blade_w2.vtp test_blade_w3.vtp test_blade_w4.vtp""".strip().split()

    output_files = [os.path.basename(i) for i in glob.glob(os.path.join(workdir, "*"))]

    for i in expected_files:
        assert i in output_files
