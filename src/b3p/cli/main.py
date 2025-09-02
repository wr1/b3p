"""Main entry point for the B3P CLI application."""

import logging
from rich.logging import RichHandler

handler = RichHandler(rich_tracebacks=True, show_time=False)
# handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
logging.basicConfig(handlers=[handler], level=logging.INFO)

logger = logging.getLogger(__name__)
# File handler for output.log
file_handler = logging.FileHandler("output.log")
file_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
logger.addHandler(file_handler)

from treeparse import cli, option
from .build_app import build_cli
from .two_d_app import twod_cli
from .ccx_app import ccx_cli
from .clean_app import clean_cli
from .validate_app import validate_cli
from ..geom_app.cli import geom_cli
from ..mesh_app.cli import mesh_cli  # Add the new mesh_cli
from .ccblade_app import ccblade_cli

from pathlib import Path

clean_cli.sort_key = 10
validate_cli.sort_key = 20
geom_cli.sort_key = 25
mesh_cli.sort_key = 30  # Add mesh_cli
build_cli.sort_key = 35
twod_cli.sort_key = 40
ccx_cli.sort_key = 50
ccblade_cli.sort_key = 60

app = cli(
    name="b3p",
    help="Blade Design CLI",
    line_connect=True,
    show_types=True,
    show_defaults=True,
    max_width=120,
    # Removed --yml from top level
    subgroups=[
        clean_cli,
        validate_cli,
        geom_cli,
        mesh_cli,  # Add mesh_cli
        build_cli,
        twod_cli,
        ccx_cli,
        ccblade_cli,
    ],
)


def main():
    app.run()


if __name__ == "__main__":
    main()
