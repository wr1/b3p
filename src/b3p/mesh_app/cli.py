#!/usr/bin/env python3
"""CLI for Mesh App using treeparse."""

from pathlib import Path
from treeparse import cli, command, option
from .step import MeshStep
from ..cli.yml_portable import yaml_make_portable


# Treeparse CLI for mesh_app
def run_mesh_callback(yml: Path):
    """Callback for running the mesh step."""
    step = MeshStep(str(yml))
    step.run()


mesh_cli = cli(
    name="mesh",
    help="Build blade mesh",
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

mesh_cli.commands.append(
    command(
        name="run",
        help="Run mesh building",
        callback=run_mesh_callback,
        arguments=[],
    )
)
