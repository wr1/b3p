from pathlib import Path
import os
import shutil

class CleanApp:
    def __init__(self, state):
        self.state = state

    def clean(self, yml: Path):
        dct = self.state.load_yaml(yml)
        prefix = dct["general"]["workdir"]
        if os.path.isdir(prefix):
            shutil.rmtree(prefix)
            print(f"** Removing workdir {prefix}")
        else:
            print(f"** Workdir {prefix} does not exist")
