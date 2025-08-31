import logging
from pathlib import Path
import os
import shutil
from ..models.config import BladeConfig  # Import for type hinting

from treeparse import cli, command, option

logger = logging.getLogger(__name__)


class CleanApp:
    """CLI app for cleaning up directories based on configuration."""

    def __init__(self, state, yml: Path):
        self.state = state
        self.yml = yml
        self.state.load_yaml(self.yml)

    def clean(self):
        """Remove the working directory specified in the configuration."""
        workdir_path = self.state.get_workdir()  # Path(config.general.workdir)
        if workdir_path.is_dir():
            shutil.rmtree(workdir_path)
            logger.info(f"Removed workdir {workdir_path}")
        else:
            logger.info(f"Workdir {workdir_path} does not exist")


def clean_callback(yml: Path):
    state = AppState.get_instance()
    app = CleanApp(state, yml)
    app.clean()

clean_cli = cli(
    name="clean",
    help="Clean working directory",
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

clean_cli.commands.append(
    command(
        name="run",
        help="Clean the working directory",
        callback=clean_callback,
        arguments=[],
    )
)
