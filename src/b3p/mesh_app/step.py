#!/usr/bin/env python3
"""Statesman step for building blade mesh."""

import logging
from pathlib import Path
from ..models.config import BladeConfig
from .mesh import build_blade_mesh
from ..cli import yml_portable
from statesman.core.base import Statesman, ManagedFile

logger = logging.getLogger(__name__)


class MeshStep(Statesman):
    """Statesman step for building blade mesh."""

    dependent_sections = ["general", "planform", "mesh"]
    output_files = ["blade_mesh.vtp"]
    input_files = [
        ManagedFile(name="blade_geometry.vtp", non_empty=True),
        ManagedFile(name="blade_geometry.pck", non_empty=True),
    ]

    def _execute(self):
        """Execute the mesh building step."""
        # Load config using custom loader
        config_data = yml_portable.yaml_make_portable(Path(self.config_path))
        self.config = config_data.model_dump()

        # Set workdir relative to YAML file
        self.workdir = Path(self.config_path).parent / self.config["general"]["workdir"]
        self.workdir.mkdir(parents=True, exist_ok=True)  # Ensure directory exists

        # Build mesh with fixed file names
        build_blade_mesh(self.config, self.workdir)
        logger.info(f"Mesh built and saved to {self.workdir}")
