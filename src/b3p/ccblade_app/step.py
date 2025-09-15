#!/usr/bin/env python3
"""Statesman step for running CCBlade analysis."""

import logging
from pathlib import Path
from .ccblade import ccblade_run
from ..cli import yml_portable
from statesman.core.base import Statesman

logger = logging.getLogger(__name__)


class CCBladeStep(Statesman):
    """Statesman step for running CCBlade analysis."""

    dependent_sections = ["general", "aero"]
    output_files = [
        "ccblade_output.csv",
        "ccblade_bladeloads.csv",
        "ccblade_moments.csv",
    ]

    def __init__(self, config_path, force=False):
        super().__init__(config_path)
        self.force = force

    def run(self):
        if self.force:
            self._execute()
        else:
            super().run()

    def _execute(self):
        """Execute the CCBlade analysis step."""
        # Load config using custom loader
        config_data = yml_portable.yaml_make_portable(Path(self.config_path))
        self.config = config_data.model_dump()

        # Set workdir relative to YAML file
        self.workdir = Path(self.config_path).parent / self.config["general"]["workdir"]
        self.workdir.mkdir(parents=True, exist_ok=True)  # Ensure directory exists

        # Run CCBlade analysis
        ccblade = ccblade_run(self.config, Path(self.config_path).parent)
        ccblade.run()
        logger.info(f"CCBlade analysis completed and saved to {self.workdir}")
