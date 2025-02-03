import pytest
from pathlib import Path
import shutil
from b3p.cli2 import AppState, TwoDApp
import subprocess

@pytest.fixture
def app_state():
    return AppState.get_instance()

@pytest.fixture
def two_d_app(app_state):
    return TwoDApp(app_state)

def copy_example_files(tmp_path):
    example_dir = Path("examples/temp_blade_portable")
    out = "test.yml"
    for file_name in ["test_blade_portable.yml"]:  # Add all required subfiles here
        shutil.copy(example_dir / file_name, tmp_path / out)
    subprocess.run(['b3p', 'build', out], cwd=tmp_path)    

# def test_mesh2d_no_mesh2d_section(two_d_app, tmp_path):
#     yml_path = tmp_path / "config.yml"
#     yml_path.write_text("""
#     general:
#       workdir: /tmp/workdir
#       prefix: test_prefix
#     """)
#     result = two_d_app.mesh2d(yml_path)
#     assert result is None

# def test_mesh2d_no_sections(two_d_app, tmp_path):
#     yml_path = tmp_path / "config.yml"
#     yml_path.write_text("""
#     general:
#       workdir: /tmp/workdir
#       prefix: test_prefix
#     mesh2d:
#       some_other_key: value
#     """)
#     result = two_d_app.mesh2d(yml_path)
#     assert result is None

# def test_mesh2d_valid_config(two_d_app, tmp_path, mocker):
#     yml_path = tmp_path / "config.yml"
#     yml_path.write_text("""
#     general:
#       workdir: /tmp/workdir
#       prefix: test_prefix
#     mesh2d:
#       sections: [section1, section2]
#     """)
#     mocker.patch("b3p.cli2.mesh_2d.cut_blade_parallel", return_value=["mesh1", "mesh2"])
#     mocker.patch("b3p.cli2.anba4_prep.anba4_prep", return_value=["anba4_mesh1", "anba4_mesh2"])

#     result = two_d_app.mesh2d(yml_path)
#     assert result == ["anba4_mesh1", "anba4_mesh2"]

def test_mesh2d_blade_test_yaml(two_d_app, tmp_path):
    print(tmp_path)
    copy_example_files(tmp_path)
    yml_path = tmp_path / "test.yml"
    result = two_d_app.mesh2d(yml_path)
    print(result)
    assert result is not None, "mesh2d should return a non-None result"