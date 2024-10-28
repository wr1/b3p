import pytest
import glob
import subprocess
import os
from b3p import b3p_cli


example_dir = "examples"
workdir = os.path.join(example_dir, "temp_blade_portable")


def run_tbuild():
    """Run step s1."""
    cd = os.getcwd()
    os.chdir(example_dir)
    cli = b3p_cli.cli("blade_test.yml")
    cli.build()
    cli.bondline()
    cli.ccxprep()
    # cli.ccxrun()
    cli.mesh2d(90.0, parallel=False)
    os.chdir(cd)


# def model_bondline():
# subprocess.run(
#     ["b3p", "--yml=blade_test.yml", "build"], check=True, cwd=example_dir
# )


def build_output_exists():
    """Check if the output of the build exists."""
    print(workdir)
    return os.path.exists(workdir)


@pytest.fixture(scope="session")
def run_test_build():
    """Fixture to run blade build."""
    if not build_output_exists():
        run_tbuild()
    else:
        print("build output exists.")


def test_build_success(run_test_build):
    """Test if the build was successful."""
    assert build_output_exists()


def test_build_has_all_files(run_test_build):
    """Test if the build has all the expected files."""
    expected_files = """          test_blade_base.vtp    test_blade.png           test_blade_shell.vtp      test_blade_w0.txt         test_blade_w1.vtp         test_blade_w3_dr.vtu      test_blade_w4_dr.vtu
__matdb.yml        test_blade_joined.vtu  test_blade_portable.yml  test_blade.var            test_blade_w0.vtp         test_blade_w2_points.txt  test_blade_w3_points.txt  test_blade_w4_points.txt
material_map.json  test_blade_loads.png   test_blade_sca_50.csv    test_blade.vtp            test_blade_w1_points.txt  test_blade_w2.txt         test_blade_w3.txt         test_blade_w4.txt
__plybook.pck      test_blade.pck         test_blade_shell_dr.vtu  test_blade_w0_points.txt  test_blade_w1.txt         test_blade_w2.vtp         test_blade_w3.vtp         test_blade_w4.vtp
    """.strip().split()

    output_files = [os.path.basename(i) for i in glob.glob(os.path.join(workdir, "*"))]

    for i in expected_files:
        assert i in output_files
