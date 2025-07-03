"""CLI app for validating YAML configuration files."""

import logging
from pathlib import Path
from b3p.cli.app_state import AppState
from b3p.models.config import BladeConfig  # Updated import for fixed config

logger = logging.getLogger(__name__)


class ValidateApp:
    """Validate a YAML configuration file."""

    def __init__(self, yml: Path):
        self.yml = yml

    def validate(self):
        """Load and validate the YAML file."""
        state = AppState.get_instance()
        try:
            config = state.load_yaml(self.yml)  # Uses existing load_yaml to validate
            print(config.materials)
            logger.info(f"YAML file {self.yml} is valid.")
            return True
        except Exception as e:
            logger.error(f"YAML validation failed for {self.yml}: {e}")
            return False


# Note: To integrate this into the CLI, ensure the main CLI script (e.g., in src/b3p/cli/__init__.py or a main.py) adds a subcommand for 'validate'.
