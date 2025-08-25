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
        self.config = self.state.load_yaml(yml)  # Load config once

    def mesh2d(self, rotz=0.0, parallel=True):
        """Create 2D meshes from blade sections."""
        if self.config.mesh2d is None:
            logger.error("No mesh2d section in config")
            return
        if "sections" not in self.config.mesh2d:
            logger.error("No sections in mesh2d section in config")
            return
        sections = self.config.mesh2d["sections"]

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
        if not meshes:
            meshes = self.mesh2d()  # This will use self.config

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

        # logger.info(f"Running command: {conda_command}")
        logger.info(" ".join(conda_command))
        result = subprocess.run(
            conda_command,
            capture_output=True,
            text=True,
            # shell=True,  # Use shell to ensure proper interpretation
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
        if self.config.general.workdir:
            workdir = Path(self.config.general.workdir) / "2d"
            if not workdir.exists():
                logger.info(f"Workdir {workdir} does not exist - nothing to clean")
                return

            try:
                shutil.rmtree(workdir)
                logger.info(f"Removed workdir {workdir}")
            except Exception as e:
                logger.error(f"Failed to remove workdir {workdir}: {e}")

from treeparse import cli, command, argument, option
from .app_state import AppState

def run_callback(yml: Path, rotz: float, parallel: bool, anba_env: str):
    state = AppState.get_instance()
    app = TwoDApp(state, yml)
    app.mesh2d(rotz=rotz, parallel=parallel)
    app.run_anba4(anba_env=anba_env)

def mesh2d_callback(yml: Path, rotz: float, parallel: bool):
    state = AppState.get_instance()
    app = TwoDApp(state, yml)
    app.mesh2d(rotz=rotz, parallel=parallel)

def run_anba4_callback(yml: Path, anba_env: str):
    state = AppState.get_instance()
    app = TwoDApp(state, yml)
    app.run_anba4(anba_env=anba_env)

def clean_callback(yml: Path):
    state = AppState.get_instance()
    app = TwoDApp(state, yml)
    app.clean()

twod_cli = cli(
    name="2d",
    help="2D mesh and ANBA4 operations",
    line_connect=True,
    show_types=True,
    show_defaults=True,
)

twod_cli.commands.append(
    command(
        name="run",
        help="Run full 2D process",
        callback=run_callback,
        arguments=[],
        options=[
            option(flags=["--rotz", "-r"], default=0.0, arg_type=float, help="Rotation around Z-axis (degrees)"),
            option(flags=["--parallel", "-P"], is_flag=True, default=True, help="Enable parallel processing"),
            option(flags=["--anba-env", "-e"], default="anba4-env", arg_type=str, help="Conda environment for ANBA4"),
        ],
    )
)

twod_cli.commands.append(
    command(
        name="mesh2d",
        help="Create 2D meshes",
        callback=mesh2d_callback,
        arguments=[],
        options=[
            option(flags=["--rotz", "-r"], default=0.0, arg_type=float, help="Rotation around Z-axis (degrees)"),
            option(flags=["--parallel", "-P"], is_flag=True, default=True, help="Enable parallel processing"),
        ],
    )
)

twod_cli.commands.append(
    command(
        name="run-anba4",
        help="Run ANBA4 on 2D meshes",
        callback=run_anba4_callback,
        arguments=[],
        options=[
            option(flags=["--anba-env", "-e"], default="anba4-env", arg_type=str, help="Conda environment for ANBA4"),
        ],
    )
)

twod_cli.commands.append(
    command(
        name="clean",
        help="Remove msec* files",
        callback=clean_callback,
        arguments=[],
    )
)
