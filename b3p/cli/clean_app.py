from pathlib import Path
import os
import shutil

class CleanApp:
    def __init__(self, state, yml: Path):
        self.state = state
        self.yml = yml

    def clean(self):
        dct = self.state.load_yaml(self.yml)
        prefix = dct["general"]["workdir"]
        if os.path.isdir(prefix):
            shutil.rmtree(prefix)
            print(f"** Removing workdir {prefix}")
        else:
            print(f"** Workdir {prefix} does not exist")
