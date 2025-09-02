#!/usr/bin/env python3
"""Test for geometry consistency using prepared test data."""

import pytest
import subprocess
from pathlib import Path

TEST_DIR = Path("tests/geom_app/test_data")
OUTPUT_DIR = TEST_DIR / "output"
REFERENCE_DIR = TEST_DIR / "reference"
TEST_YAML = TEST_DIR / "test_config.yml"


@pytest.fixture(scope="module")
def run_geometry():
    """Fixture to run the geometry app once."""
    if not TEST_YAML.exists():
        pytest.skip("Test data not prepared. Run test_geom_consistency.sh first.")
    result = subprocess.run(
        ["b3p", "geom", "-y", str(TEST_YAML), "run"], capture_output=True, text=True
    )
    assert result.returncode == 0, f"Geometry run failed: {result.stderr}"
    return result


@pytest.mark.skipif(
    not TEST_YAML.exists(),
    reason="Test data not prepared. Run test_geom_consistency.sh first.",
)
def test_geometry_run(run_geometry):
    """Test that the geometry run succeeds."""
    assert run_geometry.returncode == 0


@pytest.mark.skipif(
    not TEST_YAML.exists(),
    reason="Test data not prepared. Run test_geom_consistency.sh first.",
)
def test_blade_geometry_variables_json():
    """Test comparison of blade_geometry_variables.json."""
    output_json = OUTPUT_DIR / "blade_geometry_variables.json"
    ref_json = REFERENCE_DIR / "blade_geometry_variables.json"
    assert output_json.exists(), f"Output JSON not found: {output_json}"
    assert ref_json.exists(), f"Reference JSON not found: {ref_json}"
    with open(output_json) as f:
        output_data = f.read()
    with open(ref_json) as f:
        ref_data = f.read()
    assert output_data == ref_data, "Geometry variables differ"


@pytest.mark.skipif(
    not TEST_YAML.exists(),
    reason="Test data not prepared. Run test_geom_consistency.sh first.",
)
def test_blade_geometry_vtp():
    """Test comparison of blade_geometry.vtp file size."""
    output_file = OUTPUT_DIR / "blade_geometry.vtp"
    ref_file = REFERENCE_DIR / "blade_geometry.vtp"
    assert output_file.exists(), f"Output VTP not found: {output_file}"
    assert ref_file.exists(), f"Reference VTP not found: {ref_file}"
    output_size = output_file.stat().st_size
    ref_size = ref_file.stat().st_size
    assert output_size == ref_size, f"VTP sizes differ: {output_size} vs {ref_size}"


@pytest.mark.skipif(
    not TEST_YAML.exists(),
    reason="Test data not prepared. Run test_geom_consistency.sh first.",
)
def test_blade_geometry_pck():
    """Test comparison of blade_geometry.pck file size."""
    output_file = OUTPUT_DIR / "blade_geometry.pck"
    ref_file = REFERENCE_DIR / "blade_geometry.pck"
    assert output_file.exists(), f"Output PCK not found: {output_file}"
    assert ref_file.exists(), f"Reference PCK not found: {ref_file}"
    output_size = output_file.stat().st_size
    ref_size = ref_file.stat().st_size
    assert output_size == ref_size, f"PCK sizes differ: {output_size} vs {ref_size}"
