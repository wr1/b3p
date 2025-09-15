#!/usr/bin/env python3
"""Statesman step for building blade geometry."""

import logging
from pathlib import Path
from .geometry import build_blade_geometry
from ..cli import yml_portable
from statesman.core.base import Statesman

logger = logging.getLogger(__name__)


class GeometryStep(Statesman):
    """Statesman step for building blade geometry."""

    dependent_sections = ["general", "planform", "aero"]
    output_files = [
        "blade_geometry.vtp",
        "blade_geometry.pck",
        "blade_geometry_variables.json",
        "blade_geometry_portable.yml",
    ]
    workdir_key = "general.workdir"

    def __init__(self, config_path, force=False):
        super().__init__(config_path)
        self.force = force
        self.workdir.mkdir(parents=True, exist_ok=True)

    def run(self):
        if self.force:
            self._execute()
        else:
            super().run()

    def _execute(self):
        """Execute the geometry building step."""
        # Load config using custom loader
        config_data = yml_portable.yaml_make_portable(Path(self.config_path))
        self.config = config_data.model_dump()

        # Set workdir relative to YAML file
        self.workdir = Path(self.config_path).parent / self.config["general"]["workdir"]
        self.workdir.mkdir(parents=True, exist_ok=True)  # Ensure directory exists

        # Build geometry with fixed file names
        build_blade_geometry(self.config, self.workdir)
        yml_portable.save_yaml(
            self.workdir / "blade_geometry_portable.yml", config_data
        )
        logger.info(f"Geometry built and saved to {self.workdir}")
