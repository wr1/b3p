#!/usr/bin/env python3
"""CLI for CCBlade App using treeparse."""

from pathlib import Path
from treeparse import cli, command, option
from .step import CCBladeStep


# Treeparse CLI for ccblade_app
def run_ccblade_callback(yml: Path):
    """Callback for running the CCBlade step."""
    step = CCBladeStep(str(yml))
    step.run()


ccblade_cli = cli(
    name="ccblade",
    help="Run CCBlade analysis",
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

ccblade_cli.commands.append(
    command(
        name="run",
        help="Run CCBlade analysis",
        callback=run_ccblade_callback,
        arguments=[],
    )
)
