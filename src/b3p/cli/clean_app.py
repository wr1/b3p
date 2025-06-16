import logging
from pathlib import Path
import os
import shutil

logger = logging.getLogger(__name__)


class CleanApp:
    def __init__(self, state, yml: Path):
        self.state = state
        self.yml = yml

    def clean(self):
        dct = self.state.load_yaml(self.yml)
        prefix = dct["general"]["workdir"]
        if os.path.isdir(prefix):
            shutil.rmtree(prefix)
            logger.info(f"Removed workdir {prefix}")
        else:
            logger.info(f"Workdir {prefix} does not exist")
