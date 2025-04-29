from pathlib import Path
import os
import subprocess
import shutil
from b3p.anba import anba4_prep
from b3p.anba import mesh_2d
import glob


class TwoDApp:
    def __init__(self, state, yml: Path):
        self.state = state
        self.yml = yml

    def mesh2d(self, rotz=0.0, parallel=True):
        dct = self.state.load_yaml(self.yml)
        if "mesh2d" not in dct:
            print("** No mesh2d section in yml file")
            return
        if "sections" not in dct["mesh2d"]:
            print("** No sections in mesh2d section in yml file")
            return
        sections = dct["mesh2d"]["sections"]
        # yml_dir = self.yml.parent

        drape_prefix = self.state.get_prefix("drape")
        mesh_prefix = self.state.get_prefix("mesh")
        # prefix = self.state.get_prefix()
        # prefix = os.path.join(
        #     yml_dir, dct["general"]["workdir"], dct["general"]["prefix"]
        # )

        section_meshes = mesh_2d.cut_blade_parallel(
            f"{drape_prefix}_joined.vtu",
            sections,
            if_bondline=False,
            rotz=rotz,
            var=f"{mesh_prefix}_variables.json",
            parallel=parallel,
        )
        if not section_meshes or not all(
            Path(mesh).exists() for mesh in section_meshes
        ):
            print("** Error: Section meshes were not generated correctly.")
            return []
        return anba4_prep.anba4_prep(section_meshes, parallel=parallel)

    def run_anba4(self, anba_env="anba4-env"):
        prefix = self.state.get_prefix("drape")
        meshes = glob.glob(os.path.join(prefix.parent, "2d", "msec_*.xdmf"))
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
        dct = self.state.load_yaml(self.yml)
        if meshes is None:
            meshes = self.mesh2d()
        # yml_dir = self.yml.parent

        # prefix = self.state.get_prefix("drape")
        material_map = f"{prefix.parent}/material_map.json"
        script_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "anba", "anba4_solve.py")
        )
        print(f"** Running ANBA4 using {script_path} in env {anba_env}")
        result = subprocess.run(
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
            capture_output=True,
            text=True,
            env={
                **os.environ.copy(),
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "OMP_NUM_THREADS": "1",
                "CUDA_VISIBLE_DEVICES": "-1",
            },
        )
        if result.returncode != 0:
            print(f"** Error: ANBA4 script failed with return code {result.returncode}")
            print(f"** Stdout: {result.stdout}")
            print(f"** Stderr: {result.stderr}")
        else:
            print(f"** ANBA4 script completed successfully.")
            print(f"** Stdout: {result.stdout}")
        return result.returncode

    def clean(self):
        dct = self.state.load_yaml(self.yml)
        workdir = self.state.get_prefix("drape")
        # remove the 2d subdir in workdir
        workdir = workdir.parent / "2d"
        if not workdir.exists():
            print(f"** Workdir {workdir} does not exist - nothing to clean")
            return

        try:
            shutil.rmtree(workdir)
            print(f"Removed workdir {workdir}")
        except Exception as e:
            print(f"Failed to remove workdir {workdir}: {e}")
