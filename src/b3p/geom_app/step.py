#!/usr/bin/env python3
"""Statesman step for building blade geometry."""

import logging
from pathlib import Path
from ..models.config import BladeConfig
from .geometry import build_blade_geometry
from ..cli import yml_portable
from statesman.core.base import Statesman, ManagedFile

logger = logging.getLogger(__name__)


class GeometryStep(Statesman):
    """Statesman step for building blade geometry."""

    dependent_sections = ["general", "planform", "aero"]
    output_files = ["blade_geometry.vtu", "blade_geometry.pck", "blade_geometry_variables.json", "blade_geometry_portable.yml"]

    def _execute(self):
        """Execute the geometry building step."""
        # Load config using custom loader
        config_data = yml_portable.yaml_make_portable(Path(self.config_path))
        self.config = config_data

        # Build geometry with fixed file names
        build_blade_geometry(self.config.model_dump(), self.workdir)
        yml_portable.save_yaml(self.workdir / "blade_geometry_portable.yml", self.config)
        logger.info(f"Geometry built and saved to {self.workdir}")
