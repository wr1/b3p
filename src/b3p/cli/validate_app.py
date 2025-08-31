"""CLI app for validating YAML configuration files."""

import logging
from pathlib import Path
from b3p.cli.app_state import AppState
from b3p.models.config import BladeConfig  # Updated import for fixed config
from treeparse import cli, command, option

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


def validate_callback(yml: Path):
    app = ValidateApp(yml)
    app.validate()

validate_cli = cli(
    name="validate",
    help="Validate YAML configuration",
    line_connect=True,
    show_types=True,
    show_defaults=True,
    options=[
        option(
            flags=["--yml", "-y"],
            arg_type=Path,
            required=True,
            help="Path to YAML config file",
        ),
    ],
)

validate_cli.commands.append(
    command(
        name="run",
        help="Validate the YAML file",
        callback=validate_callback,
        arguments=[],
    )
)
