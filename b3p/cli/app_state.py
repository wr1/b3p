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
        prefix = os.path.join(
            self.dct["general"]["workdir"], self.dct["general"]["prefix"]
        )
        if not os.path.isdir(prefix):
            os.makedirs(prefix)

    def get_prefix(self):
        if self.dct is None:
            return None
        wd = Path(self.dct["general"]["workdir"])
        prefix = wd / self.dct["general"]["prefix"]
        return prefix

    def reset(self):
        self.dct = None
