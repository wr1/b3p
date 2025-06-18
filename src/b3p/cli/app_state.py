from pathlib import Path
import os
from b3p.cli import yml_portable
from b3p.laminates import build_plybook
import logging

logger = logging.getLogger(__name__)


class AppState:
    """Singleton to manage application state."""

    _state = None

    def __init__(self):
        self.dct = None
        self.yml_dir = None  # Store directory of YAML file

    @classmethod
    def get_instance(cls):
        if cls._state is None:
            cls._state = cls()
        return cls._state

    def load_yaml(self, yml: Path):
        if self.dct is None:
            self.yml_dir = Path(yml).parent  # Store YAML file's directory
            self.dct = yml_portable.yaml_make_portable(yml)
            self._set_workdir()
            self.make_workdir()
            self.expand_chamfered_cores()

        logger.info(f"Loaded YAML data from {self.dct['general']['workdir']}")

        return self.dct

    def _set_workdir(self):
        """Set workdir to dir(ymlfile)/output/ unless specified in YAML."""
        if self.dct and "general" in self.dct and "workdir" in self.dct["general"]:
            workdir = self.dct["general"]["workdir"]
            # If workdir is relative, resolve it relative to yml_dir
            if not os.path.isabs(workdir):
                workdir = self.yml_dir / workdir
            logger.info(f"Setting workdir to {workdir}")
            self.dct["general"]["workdir"] = str(workdir)
        else:
            # Default to dir(ymlfile)/output/
            self.dct["general"]["workdir"] = str(self.yml_dir / "output")

    def expand_chamfered_cores(self):
        if self.dct:
            self.dct = build_plybook.expand_chamfered_cores(self.dct)

    def make_workdir(self):
        if self.dct is None:
            raise ValueError("No YAML data loaded")
        workdir = Path(self.dct["general"]["workdir"])
        if not workdir.is_dir():
            os.makedirs(workdir, exist_ok=True)

    def get_prefix(self, subdir=None):
        if self.dct is None:
            return None
        wd = Path(self.dct["general"]["workdir"])
        prefix_name = self.dct["general"].get("prefix", "b3p")
        if subdir is None:
            prefix = wd / prefix_name
        else:
            prefix = wd / subdir / prefix_name
            if not (wd / subdir).is_dir():
                os.makedirs(wd / subdir, exist_ok=True)
        return prefix

    def get_workdir(self):
        if self.dct is None:
            return None

        return Path(self.dct["general"]["workdir"])

    def reset(self):
        self.dct = None
        self.yml_dir = None
