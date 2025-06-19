from pathlib import Path
import os
from b3p.cli import yml_portable
from b3p.laminates import build_plybook
import logging
from typing import Optional, Dict

logger = logging.getLogger(__name__)


class AppState:
    """Singleton to manage application state."""

    _state: Optional["AppState"] = None

    def __init__(self):
        self.dct: Optional[Dict] = None
        self.yml_dir: Optional[Path] = None
        self.workdir_path: Optional[Path] = None  # Stores the resolved workdir as Path

    @classmethod
    def get_instance(cls):
        if cls._state is None:
            cls._state = cls()
        return cls._state

    def load_yaml(self, yml: Path) -> Optional[Dict]:
        if self.dct is None:
            yml = yml.resolve()
            self.yml_dir = yml.parent
            self.dct = yml_portable.yaml_make_portable(yml)
            self._set_workdir()
            self.make_workdir()
            self.expand_chamfered_cores()

        logger.info(f"Loaded YAML data from {self.dct['general']['workdir']}")
        return self.dct

    def _set_workdir(self):
        if self.dct and "general" in self.dct and "workdir" in self.dct["general"]:
            workdir_str = self.dct["general"]["workdir"]
            workdir_path = Path(workdir_str)
            if not workdir_path.is_absolute():
                workdir_path = self.yml_dir / workdir_path
                workdir_path = workdir_path.resolve()
            self.workdir_path = workdir_path
            logger.info(f"Setting workdir to {self.workdir_path}")
        else:
            default_workdir = self.yml_dir / "output"
            self.workdir_path = default_workdir.resolve()

    def expand_chamfered_cores(self):
        if self.dct:
            self.dct = build_plybook.expand_chamfered_cores(
                self.dct,
                self.workdir_path / (self.dct["general"]["prefix"] + "_expanded.yml"),
            )

    def make_workdir(self):
        if self.dct is None:
            raise ValueError("No YAML data loaded")
        if self.workdir_path is None:
            raise ValueError("Workdir not set")
        if not self.workdir_path.is_dir():
            os.makedirs(self.workdir_path, exist_ok=True)

    def get_prefix(self, subdir: Optional[str] = None) -> Optional[Path]:
        if self.dct is None:
            logger.warning("No YAML data loaded, returning None")
            return None
        if self.workdir_path is None:
            logger.warning("Workdir not set, returning None")
            return None
        prefix_name = self.dct["general"].get("prefix", "b3p")
        if subdir is None:
            prefix = self.workdir_path / prefix_name
        else:
            prefix = self.workdir_path / subdir / prefix_name
            if not (self.workdir_path / subdir).is_dir():
                os.makedirs(self.workdir_path / subdir, exist_ok=True)
        return prefix.resolve()

    def get_workdir(self) -> Optional[Path]:
        if self.workdir_path is None:
            logger.warning("Workdir not set, returning None")
            return None
        return self.workdir_path

    def reset(self):
        self.dct = None
        self.yml_dir = None
        self.workdir_path = None
