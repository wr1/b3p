from pathlib import Path
import os
from b3p import yml_portable
from b3p.laminates import build_plybook


class AppState:
    """Singleton to manage application state."""

    _state = None

    def __init__(self):
        self.dct = None

    @classmethod
    def get_instance(cls):
        if cls._state is None:
            cls._state = cls()
        return cls._state

    def load_yaml(self, yml: Path):
        if self.dct is None:
            self.dct = yml_portable.yaml_make_portable(yml)
            self.make_workdir(yml)
            self.expand_chamfered_cores()
        return self.dct

    def expand_chamfered_cores(self):
        if self.dct:
            self.dct = build_plybook.expand_chamfered_cores(self.dct)

    def make_workdir(self, yml: Path):
        if self.dct is None:
            self.load_yaml(yml)
        if not os.path.isdir(self.dct["general"]["workdir"]):
            os.makedirs(self.dct["general"]["workdir"])

    def get_prefix(self, subdir=None):
        if self.dct is None:
            return None
        wd = Path(self.dct["general"]["workdir"])
        if subdir is None:
            prefix = wd / self.dct["general"]["prefix"]
        else:
            prefix = wd / subdir / self.dct["general"]["prefix"]
            if not os.path.isdir(wd / subdir):
                os.makedirs(wd / subdir)
        return prefix

    def get_workdir(self):
        if self.dct is None:
            return None
        wd = Path(self.dct["general"]["workdir"])
        return wd

    def reset(self):
        self.dct = None
