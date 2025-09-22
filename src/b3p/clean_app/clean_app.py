import logging
from pathlib import Path
import shutil
from ruamel.yaml import YAML
from treeparse import cli, command, option

logger = logging.getLogger(__name__)


class CleanApp:
    """CLI app for cleaning up directories based on configuration."""

    def __init__(self, yml: Path):
        self.yml = yml
        self.workdir = self._get_workdir()

    def _get_workdir(self):
        yaml = YAML()
        with open(self.yml) as f:
            data = yaml.load(f)
        workdir = data.get('general', {}).get('workdir', 'output')
        if not Path(workdir).is_absolute():
            workdir = self.yml.parent / workdir
        return workdir.resolve()

    def clean(self):
        """Remove the working directory specified in the configuration."""
        if self.workdir.is_dir():
            shutil.rmtree(self.workdir)
            logger.info(f"Removed workdir {self.workdir}")
        else:
            logger.info(f"Workdir {self.workdir} does not exist")


def clean_callback(yml: Path):
    app = CleanApp(yml)
    app.clean()


clean_cli = cli(
    name="clean",
    help="Clean working directory",
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

clean_cli.commands.append(
    command(
        name="run",
        help="Clean the working directory",
        callback=clean_callback,
        arguments=[],
    )
)