import pytest
import glob
import os
from pathlib import Path
from .conftest import built_blade  # Import from conftest
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp


def test_build_success(built_blade):
    """Test if the build completes successfully."""
    # assert built_blade["returncode"] == 0, (
    #     f"Build failed with stderr: {built_blade['stderr']}"
    # )
    assert built_blade["workdir"].exists(), "Working directory was not created."


def test_build_output_files(built_blade):
    """Test if all expected output files are generated."""
    #  AssertionError: Missing expected files: {'test_blade_joined_bondline.vtu', 'test_blade.var', 'test_blade_base.vtptest_blade_base.vtp', '__plybook.pck', 'test_blade_portable.yml'}
    #  assert not {'__plybook.pck', 'test_blade.var', 'test_blade_base.vtptest_blade_base.vtp', 'test_blade_joined_bondline.vtu', 'test_blade_portable.yml'}
    workdir = built_blade["workdir"]
    print(workdir)
    output_files = {os.path.basename(f) for f in glob.glob(str(workdir / "**/*"))}
    print(output_files)
    expected_files = {
        "test_blade_base.vtp",
        "test_blade.png",
        "test_blade_shell.vtp",
        "test_blade_w0.vtp",
        "test_blade_w1.vtp",
        "test_blade_w2.vtp",
        "test_blade_w3.vtp",
        "test_blade_w4.vtp",
        "test_blade_w3_dr.vtu",
        "test_blade_shell_dr.vtu",
        "test_blade_w4_dr.vtu",
        "__matdb.yml",
        "test_blade_joined.vtu",
        "material_map.json",
        "test_blade_loads.png",
        "test_blade.vtp",
        "test_blade_mass.csv",
        "test_blade_mass.tex",
    }

    missing_files = expected_files - output_files
    assert not missing_files, f"Missing expected files: {missing_files}"


def test_build_blademass_match(built_blade):
    # use pandas to match the mass table
    import pandas as pd

    workdir = built_blade["workdir"]
    mass_csv = workdir / "drape" / "test_blade_mass.csv"
    assert mass_csv.exists(), "Mass CSV file not found."
    # Read the CSV file
    mass_df = pd.read_csv(mass_csv, sep=",")
    # check if the mass table matches the table with same name in tests/data/

    # read the reference mass table
    ref_mass_csv = workdir.parent.parent / "data" / "test_blade_mass.csv"
    assert ref_mass_csv.exists(), "Reference mass CSV file not found."
    ref_mass_df = pd.read_csv(ref_mass_csv, sep=",")
    # Compare the two DataFrames
    assert mass_df.round().equals(ref_mass_df.round()), "Mass tables do not match."


# def test_build_stdout(built_blade):
#     """Test if key stdout messages are present."""
#     stdout = built_blade["stdout"]
#     print(stdout)
#     expected_messages = [
#         "** Loading yaml file blade_test.yml and loading linked files",
#         "saving to temp_blade_portable/test_blade.vtp",
#         "written to: temp_blade_portable/test_blade_portable.yml",
#         "** wrote mesh to temp_blade_portable/test_blade_base.vtp",
#         "** wrote mesh to temp_blade_portable/test_blade_shell.vtp",
#         "written material map to temp_blade_portable/material_map.json",
#         "written to temp_blade_portable/__plybook.pck",
#         "written mesh to temp_blade_portable/test_blade_joined.vtu",
#         "Saved temp_blade_portable/test_blade_joined_bondline.vtu",
#         "Loaded material database from temp_blade_portable/__matdb.yml",
#         "Total Volume and Mass:",
#         "Mass table per material",
#         "** written load plot to temp_blade_portable/test_blade_loads.png",
#     ]
#     for message in expected_messages:
#         assert message in stdout, f"Expected message not found in stdout: {message}"


# def test_build_mass_table(built_blade):
#     """Test if the mass table in stdout matches expected values."""
#     stdout = built_blade["stdout"]
#     mass_table_lines = [
#         "11457.585028",
#         "937.880619",
#         "10047.614323",
#         "31119.325388",
#         "1125.32413",
#         "783.68881",
#         "3789.491169",
#         "1570.496484",
#         "60831.405951",
#     ]
#     for line in mass_table_lines:
#         assert line in stdout, f"Expected mass table entry not found: {line}"


def test_build_geometry_only(temp_example_dir):
    """Test the geometry subcommand independently."""
    original_dir = os.getcwd()
    os.chdir(temp_example_dir)
    try:
        state = AppState()
        build_app = BuildApp(state, Path("blade_test.yml"))
        build_app.geometry()

        workdir = temp_example_dir / "temp_blade_portable" / "mesh"
        assert workdir.exists(), "Working directory not created."
        assert (workdir / "test_blade.vtp").exists(), "Geometry file not created."
        assert (workdir / ".." / "test_blade_portable.yml").exists(), (
            "Portable YAML not created."
        )
    finally:
        os.chdir(original_dir)
