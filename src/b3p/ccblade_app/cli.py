#!/usr/bin/env python3
"""CLI for CCBlade App using treeparse."""

from pathlib import Path
from treeparse import cli, command, option
from .step import CCBladeStep


# Treeparse CLI for ccblade_app
def run_ccblade_callback(yml: Path, force: bool = False):
    """Callback for running the CCBlade step."""
    step = CCBladeStep(str(yml), force=force)
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
        option(
            flags=["--force", "-f"],
            arg_type=bool,
            default=False,
            help="Force run despite statesman",
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
