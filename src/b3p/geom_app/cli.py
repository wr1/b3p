#!/usr/bin/env python3
"""CLI for Geometry App using treeparse."""

from pathlib import Path
from treeparse import cli, command, option
from .step import GeometryStep
from ..cli.yml_portable import yaml_make_portable


# Treeparse CLI for geom_app
def run_geometry_callback(yml: Path):
    """Callback for running the geometry step."""
    step = GeometryStep(str(yml))
    step.run()

geom_cli = cli(
    name="geom",
    help="Build blade geometry",
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

geom_cli.commands.append(
    command(
        name="run",
        help="Run geometry building",
        callback=run_geometry_callback,
        arguments=[],
    )
)
