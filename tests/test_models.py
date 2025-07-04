"""Tests for b3p models."""

import pytest
from pathlib import Path
from b3p.models.config import BladeConfig
from b3p.cli.yml_portable import yaml_make_portable
import logging

logger = logging.getLogger(__name__)
def test_default_yaml():
    """Test loading a default YAML configuration."""
    config = yaml_make_portable(Path("examples/blade_test.yml"))
    assert isinstance(config, BladeConfig)
    assert config.general.workdir == "temp_blade"
    assert config.general.prefix == "test_blade"
    # assert config.mesh.bondline["type"] == "default"

def test_invalid_yaml(tmp_path):
    """Test loading an invalid YAML configuration."""
    invalid_yaml = tmp_path / "invalid.yml"
    invalid_yaml.write_text("""
general:
  prefix: test
aero:
  airfoils:
    0.1:
      xy: [[1, 2, 3]]  # Invalid: should be lists of exactly two floats
""")
    with pytest.raises(ValueError, match="xy must be a list of \\[x, y\\] coordinates"):
        yaml_make_portable(invalid_yaml)

def test_airfoil_path_loading(tmp_path):
    """Test loading configuration with a valid airfoil path."""
    test_file = tmp_path / "airfoil_test.dat"
    test_file.write_text("0.0 0.0\n1.0 1.0")
    yaml_content = f"""
aero:
  airfoils:
    0.5:
      xy: "{test_file}"
"""
    yaml_file = tmp_path / "test.yml"
    yaml_file.write_text(yaml_content)


    config = yaml_make_portable(yaml_file)
    logger.debug(f"Loaded airfoil: {config.aero}")
    assert isinstance(config.aero.airfoils[0.5].xy, list)
    assert len(config.aero.airfoils[0.5].xy) > 0  # Ensure points are loaded
