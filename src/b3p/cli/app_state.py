from pathlib import Path
import os
from ..models.config import BladeConfig
from ..laminates import build_plybook
import logging
from typing import Optional

logger = logging.getLogger(__name__)

class AppState:
    """Singleton to manage application state."""

    _state: Optional["AppState"] = None

    def __init__(self):
        self.config: Optional[BladeConfig] = None
        self.yml_dir: Optional[Path] = None
        self.workdir_path: Optional[Path] = None

    @classmethod
    def get_instance(cls):
        if cls._state is None:
            cls._state = cls()
        return cls._state

    def load_yaml(self, yml: Path) -> BladeConfig:
        from .yml_portable import yaml_make_portable
        if self.config is None:
            yml = yml.resolve()
            self.yml_dir = yml.parent
            self.config = yaml_make_portable(yml)
            self._set_workdir()
            self.make_workdir()
            self.expand_chamfered_cores()
        logger.info(f"Loaded YAML data from {self.config.general.workdir}")
        return self.config

    def _set_workdir(self):
        if self.config and self.config.general.workdir:
            workdir_path = Path(self.config.general.workdir)
            if not workdir_path.is_absolute():
                workdir_path = self.yml_dir / workdir_path
                workdir_path = workdir_path.resolve()
            self.workdir_path = workdir_path
            logger.info(f"Setting workdir to {self.workdir_path}")
        else:
            default_workdir = self.yml_dir / "output"
            self.workdir_path = default_workdir.resolve()

    def expand_chamfered_cores(self):
        if self.config:
            self.config = BladeConfig(**build_plybook.expand_chamfered_cores(
                self.config.dict(),
                self.workdir_path / (self.config.general.prefix + "_expanded.yml"),
            ))

    def make_workdir(self):
        if self.config is None:
            raise ValueError("No YAML data loaded")
        if self.workdir_path is None:
            raise ValueError("Workdir not set")
        if not self.workdir_path.is_dir():
            os.makedirs(self.workdir_path, exist_ok=True)

    def get_prefix(self, subdir: Optional[str] = None) -> Optional[Path]:
        if self.config is None:
            logger.warning("No YAML data loaded, returning None")
            return None
        if self.workdir_path is None:
            logger.warning("Workdir not set, returning None")
            return None
        prefix_name = self.config.general.prefix
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
        self.config = None
        self.yml_dir = None
        self.workdir_path = None
