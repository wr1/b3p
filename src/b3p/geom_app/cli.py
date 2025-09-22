#!/usr/bin/env python3
"""CLI for Geometry App using treeparse."""

from pathlib import Path
from treeparse import cli, command, option
from b3p.geom_app.step import GeometryStep
import logging


# Treeparse CLI for geom_app
def run_geometry_callback(yml: Path, force: bool = True):
    """Callback for running the geometry step."""
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.info("Starting geometry run")
    step = GeometryStep(str(yml), force=force)
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
        option(
            flags=["--force", "-f"],
            arg_type=bool,
            default=False,
            help="Force execution even if dependencies haven't changed",
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


if __name__=='__main__':
    geom_cli.run()