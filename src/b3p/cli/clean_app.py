import logging
from pathlib import Path
import os
import shutil
from ..models.config import BladeConfig  # Import for type hinting

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
