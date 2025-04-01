from pathlib import Path
import os
import subprocess
import shutil
from b3p.anba import anba4_prep
from b3p.anba import mesh_2d


class TwoDApp:
    def __init__(self, state):
        self.state = state

    def mesh2d(self, yml: Path, rotz=0.0, parallel=True):
        dct = self.state.load_yaml(yml)
        if "mesh2d" not in dct:
            print("** No mesh2d section in yml file")
            return
        if "sections" not in dct["mesh2d"]:
            print("** No sections in mesh2d section in yml file")
            return
        sections = dct["mesh2d"]["sections"]
        yml_dir = yml.parent
        prefix = os.path.join(
            yml_dir, dct["general"]["workdir"], dct["general"]["prefix"]
        )
        section_meshes = mesh_2d.cut_blade_parallel(
            f"{prefix}_joined.vtu",
            sections,
            if_bondline=False,
            rotz=rotz,
            var=f"{prefix}.var",
            parallel=parallel,
        )
        return anba4_prep.anba4_prep(section_meshes)

    def run_anba4(self, yml: Path, meshes: list = None, anba_env="anba4-env"):
        conda_path = os.environ.get("CONDA_EXE") or shutil.which("conda")
        if conda_path is None:
            print("** Conda not found - please install conda.")
            return
        result = subprocess.run(
            [conda_path, "env", "list"], capture_output=True, text=True
        )
        if result.returncode != 0 or anba_env not in result.stdout:
            print(f"** Conda environment {anba_env} not found - please create it")
            return
        else:
            print(f"** Using Conda environment for running anba4 {anba_env}")
        dct = self.state.load_yaml(yml)
        if meshes is None:
            meshes = self.mesh2d(yml)
        yml_dir = yml.parent
        material_map = os.path.join(
            yml_dir, dct["general"]["workdir"], "material_map.json"
        )
        script_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "anba", "anba4_solve.py")
        )
        print(f"** Running ANBA4 using {script_path}")
        return subprocess.run(
            [
                conda_path,
                "run",
                "-n",
                anba_env,
                "python",
                script_path,
                *meshes,
                material_map,
            ],
            env={
                **os.environ.copy(),
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "OMP_NUM_THREADS": "1",
                "CUDA_VISIBLE_DEVICES": "-1",
            },
        ).returncode

    def clean(self, yml: Path):
        dct = self.state.load_yaml(yml)
        workdir = Path(dct["general"]["workdir"])
        for msec_file in workdir.glob("msec*"):
            try:
                msec_file.unlink()
                print(f"Removed {msec_file}")
            except Exception as e:
                print(f"Failed to remove {msec_file}: {e}")
