#!/usr/bin/env python3
"""Test for CCBlade consistency using prepared test data."""

import pytest
import subprocess
from pathlib import Path

TEST_DIR = Path("tests/ccblade_app/test_data")
OUTPUT_DIR = TEST_DIR / "output"
REFERENCE_DIR = TEST_DIR / "reference"
TEST_YAML = TEST_DIR / "test_config.yml"


@pytest.fixture(scope="module")
def run_ccblade():
    """Fixture to run the CCBlade app once."""
    if not TEST_YAML.exists():
        pytest.skip("Test data not prepared. Run test_ccblade_consistency.sh first.")
    result = subprocess.run(
        ["b3p", "ccblade", "-y", str(TEST_YAML), "run"], capture_output=True, text=True
    )
    assert result.returncode == 0, f"CCBlade run failed: {result.stderr}"
    return result


@pytest.mark.skipif(
    not TEST_YAML.exists(),
    reason="Test data not prepared. Run test_ccblade_consistency.sh first.",
)
def test_ccblade_run(run_ccblade):
    """Test that the CCBlade run succeeds."""
    assert run_ccblade.returncode == 0


@pytest.mark.skipif(
    not TEST_YAML.exists(),
    reason="Test data not prepared. Run test_ccblade_consistency.sh first.",
)
def test_ccblade_output_csv():
    """Test comparison of ccblade_output.csv."""
    output_csv = OUTPUT_DIR / "ccblade_output.csv"
    ref_csv = REFERENCE_DIR / "ccblade_output.csv"
    assert output_csv.exists(), f"Output CSV not found: {output_csv}"
    assert ref_csv.exists(), f"Reference CSV not found: {ref_csv}"
    with open(output_csv) as f:
        output_data = f.read()
    with open(ref_csv) as f:
        ref_data = f.read()
    assert output_data == ref_data, "CCBlade output differs"


@pytest.mark.skipif(
    not TEST_YAML.exists(),
    reason="Test data not prepared. Run test_ccblade_consistency.sh first.",
)
def test_ccblade_bladeloads_png():
    """Test comparison of ccblade_bladeloads.png file size."""
    output_file = OUTPUT_DIR / "ccblade_bladeloads.png"
    ref_file = REFERENCE_DIR / "ccblade_bladeloads.png"
    assert output_file.exists(), f"Output PNG not found: {output_file}"
    assert ref_file.exists(), f"Reference PNG not found: {ref_file}"
    output_size = output_file.stat().st_size
    ref_size = ref_file.stat().st_size
    assert output_size == ref_size, f"PNG sizes differ: {output_size} vs {ref_size}"
