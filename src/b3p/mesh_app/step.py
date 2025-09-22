#!/usr/bin/env python3
"""Statesman step for building blade mesh using self-contained mesh_app logic."""

import logging
import shutil
from pathlib import Path
from b3p.mesh_app.mesh import build_blade_mesh  # Use self-contained mesh_app function
from ..cli import yml_portable
from statesman.core.base import Statesman, ManagedFile

logger = logging.getLogger(__name__)


class MeshStep(Statesman):
    """Statesman step for building blade mesh, self-contained in mesh_app."""

    dependent_sections = ["general", "planform", "mesh"]
    output_files = ["blade_mesh.vtp"]
    input_files = [
        ManagedFile(name="blade_geometry.vtp", non_empty=True),
        ManagedFile(name="blade_geometry.pck", non_empty=True),
    ]

    def __init__(self, config_path, force=False):
        super().__init__(config_path)
        self.force = force
        # Load config and set workdir to the config's workdir
        self.config = yml_portable.yaml_make_portable(Path(self.config_path))
        self.config = self.config.model_dump()  # Convert to dict for Statesman
        self.workdir = Path(self.config_path).parent / self.config["general"]["workdir"]
        self.workdir.mkdir(parents=True, exist_ok=True)

    def run(self):
        logger.info("Starting mesh step")
        if self.force:
            logger.info("Force mode: skipping dependency checks")
            self._execute()
        else:
            logger.info("Checking dependencies with Statesman")
            super().run()

    def _execute(self):
        """Execute the mesh building step using mesh_app's build_blade_mesh."""
        logger.info("Executing mesh building")
        # Config already loaded in __init__
        # Copy and rename input files to expected names in workdir
        prefix = self.config["general"]["prefix"]
        src_pck = self.workdir / "blade_geometry.pck"
        dst_pck = self.workdir / f"{prefix}.pck"
        if src_pck.exists():
            shutil.copy(src_pck, dst_pck)
            logger.info(f"Copied and renamed {src_pck} to {dst_pck}")

        src_vtp = self.workdir / "blade_geometry.vtp"
        dst_vtp = (
            self.workdir / f"{prefix}_base.vtp"
        )  # Match expectation in mesh_app/mesh.py
        if src_vtp.exists():
            shutil.copy(src_vtp, dst_vtp)
            logger.info(f"Copied and renamed {src_vtp} to {dst_vtp}")

        # Use self-contained build_blade_mesh from mesh_app
        build_blade_mesh(self.config, self.workdir)
        logger.info("Mesh built using mesh_app's build_blade_mesh")
        logger.info(f"Output files: blade_mesh.vtp in {self.workdir}")
