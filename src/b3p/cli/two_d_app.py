import logging
from pathlib import Path
import os
import subprocess
import shutil
from b3p.anba import anba4_prep
from b3p.anba import mesh_2d
import glob

logger = logging.getLogger(__name__)


class TwoDApp:
    def __init__(self, state, yml: Path):
        """Initialize TwoDApp with state and YAML config file."""
        self.state = state
        self.yml = yml

    def mesh2d(self, rotz=0.0, parallel=True):
        """Create 2D meshes from blade sections."""
        dct = self.state.load_yaml(Path(self.yml))
        if "mesh2d" not in dct:
            logger.error("No mesh2d section in yml file")
            return
        if "sections" not in dct["mesh2d"]:
            logger.error("No sections in mesh2d section in yml file")
            return
        sections = dct["mesh2d"]["sections"]

        drape_prefix = self.state.get_prefix("drape")
        mesh_prefix = self.state.get_prefix("mesh")
        logger.info(f"drape and mesh prefixes: {drape_prefix}, {mesh_prefix}")
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
            logger.error("Section meshes were not generated correctly")
            return []
        return anba4_prep.anba4_prep(section_meshes, parallel=parallel)

    def run_anba4(self, anba_env="anba4-env"):
        """Run ANBA4 on 2D meshes."""
        yml = self.yml
        self.state.load_yaml(yml)
        prefix = self.state.get_prefix("drape")
        meshes = glob.glob(str(Path(prefix) / "2d" / "msec_*.xdmf"))
        conda_path = os.environ.get("CONDA_EXE") or shutil.which("conda")
        if conda_path is None:
            logger.error("Conda not found - please install conda")
            return
        result = subprocess.run(
            [conda_path, "env", "list"], capture_output=True, text=True
        )
        if result.returncode != 0 or anba_env not in result.stdout:
            logger.error(f"Conda environment {anba_env} not found - please create it")
            return

        logger.info(f"Using Conda environment for running anba4 {anba_env}")
        self.state.load_yaml(yml)
        if not meshes:
            meshes = self.mesh2d()

        material_map = str(Path(prefix).parent / "material_map.json")
        script_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "anba", "anba4_solve.py")
        )
        logger.info(f"Running ANBA4 using {script_path} in env {anba_env}")
        conda_command = [
            conda_path,
            "run",
            "-n",
            anba_env,
            "python",
            script_path,
            *meshes,
            material_map,
        ]

        logger.info(" ".join(conda_command))
        result = subprocess.run(
            conda_command,
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
            logger.error(f"ANBA4 script failed with return code {result.returncode}")
            logger.error(f"Stdout: {result.stdout}")
            logger.error(f"Stderr: {result.stderr}")
        else:
            logger.info("ANBA4 script completed successfully")
            logger.debug(f"Stdout: {result.stdout}")
        return result.returncode

    def clean(self):
        """Remove 2D working directory and its contents."""
        self.state.load_yaml(self.yml)
        workdir = Path(self.state.get_prefix("2d")).parent / "2d"
        if not workdir.exists():
            logger.info(f"Workdir {workdir} does not exist - nothing to clean")
            return

        try:
            shutil.rmtree(workdir)
            logger.info(f"Removed workdir {workdir}")
        except Exception as e:
            logger.error(f"Failed to remove workdir {workdir}: {e}")
