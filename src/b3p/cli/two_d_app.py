import logging
from pathlib import Path
import os
import subprocess
import shutil
from ..anba import anba4_prep
from ..anba import mesh_2d
from statesman.core.base import Statesman, ManagedFile
from treeparse import cli, command, option

logger = logging.getLogger(__name__)


class Mesh2DStep(Statesman):
    """Step for creating 2D meshes."""

    dependent_sections = ["mesh2d"]
    output_files = ["2d_meshes.xdmf"]

    def _execute(self, rotz=0.0, parallel=True):
        if self.config.get("mesh2d") is None:
            logger.error("No mesh2d section in config")
            return
        sections = self.config["mesh2d"]["sections"]

        drape_prefix = self.workdir / self.config["general"]["prefix"]
        mesh_prefix = self.workdir / self.config["general"]["prefix"]
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


class RunAnba4Step(Statesman):
    """Step for running ANBA4."""

    input_files = [
        ManagedFile(name="2d_meshes.xdmf", non_empty=True),
    ]
    output_files = ["anba4_output.json"]
    dependent_sections = ["anba"]

    def _execute(self, anba_env="anba4-env"):
        meshes = [f for f in self.workdir.glob("2d/msec_*.xdmf")]
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

        material_map = str(self.workdir / "material_map.json")
        script_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "anba", "anba4_solve.py")
        )
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
        else:
            logger.info("ANBA4 script completed successfully")


# CLI code for 2d remains, but integrated with statesman
# ... (rest unchanged)

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
    options=[
        option(
            flags=["--yml", "-y"],
            arg_type=Path,
            required=True,
            help="Path to YAML config file",
        ),
    ],
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
