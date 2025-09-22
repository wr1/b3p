#!/usr/bin/env python3
"""CLI for Mesh App using treeparse, synced with build mesh."""

from pathlib import Path
from treeparse import cli, command, option
from b3p.mesh_app.step import MeshStep
import logging


# Treeparse CLI for mesh_app
def run_mesh_callback(yml: Path, force: bool):
    """Callback for running the mesh step, now synced."""
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.info("Starting mesh run")
    step = MeshStep(str(yml), force=force)
    step.run()


mesh_cli = cli(
    name="mesh",
    help="Build blade mesh (synced with build mesh)",
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
        option(
            flags=["--force", "-f"],
            arg_type=bool,
            default=False,
            help="Force rebuild the mesh",
        )
    ],
)

mesh_cli.commands.append(
    command(
        name="run",
        help="Run mesh building (synced)",
        callback=run_mesh_callback,
        arguments=[],
    )

)


if __name__=="__main__":
    mesh_cli.run()