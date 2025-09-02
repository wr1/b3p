import logging
from pathlib import Path
import os
import pickle
from ..models.config import BladeConfig
from . import yml_portable
from ..geometry import build_blade_geometry
from ..mesh import (
    add_load_to_mesh,
    add_te_solids,
    build_blade_structure,
    combine_meshes,
)
from ..laminates import build_plybook, drape_mesh, drape_summary
from rich.logging import RichHandler
from statesman.core.base import Statesman, ManagedFile
from treeparse import cli, command, option
from .app_state import AppState

# logging.basicConfig(handlers=[RichHandler(rich_tracebacks=True)], level=logging.INFO)
logger = logging.getLogger(__name__)


class GeometryStep(Statesman):
    """Step for building blade geometry."""

    dependent_sections = ["general", "planform", "aero"]
    output_files = ["geometry_output.vtu", "geometry_output.pck"]

    def _execute(self):
        prefix = self.workdir / self.config["general"]["prefix"]
        build_blade_geometry.build_blade_geometry(self.config, prefix)
        yml_portable.save_yaml(f"{prefix}_portable.yml", self.config)


class MeshStep(Statesman):
    """Step for meshing blade structure."""

    input_files = [
        ManagedFile(name="geometry_output.vtu", non_empty=True),
        ManagedFile(name="geometry_output.pck", non_empty=True),
    ]
    output_files = ["mesh_output.vtp"]
    dependent_sections = ["mesh"]

    def _execute(self):
        prefix = self.workdir / self.config["general"]["prefix"]
        build_blade_structure.build_blade_structure(self.config, prefix)


class DrapeStep(Statesman):
    """Step for draping plies onto mesh."""

    input_files = [
        ManagedFile(name="mesh_output.vtp", non_empty=True),
    ]
    output_files = ["drape_output.vtu"]
    dependent_sections = ["laminates"]

    def _execute(self, bondline=True):
        plybookname = "_plybook.pck"
        prefix = self.workdir / self.config["general"]["prefix"]
        mesh_prefix = self.workdir / self.config["general"]["prefix"]
        pbookpath = str(prefix) + plybookname

        build_plybook.lamplan2plies(self.config, pbookpath)
        slb = self.config["laminates"]["slabs"]
        used_grids = {slb[i]["grid"] for i in slb}

        if os.path.exists(pbookpath):
            plybook = pickle.load(open(pbookpath, "rb"))
            meshes = []
            for grid in used_grids:
                out = f"{prefix}_{grid}_dr.vtu"
                drape_mesh.drape_mesh(f"{mesh_prefix}_{grid}.vtp", plybook, grid, out)
                meshes.append(out)
            grid = f"{prefix}_joined.vtu"
            combine_meshes.combine_meshes(meshes, grid)

            if bondline:
                add_te_solids.add_bondline(
                    grid,
                    self.workdir / "material_map.json",
                    bondline_config=self.config["mesh"]["bondline"],
                )


class MassStep(Statesman):
    """Step for calculating blade mass."""

    input_files = [
        ManagedFile(name="drape_output.vtu", non_empty=True),
    ]
    output_files = ["mass_output.csv"]
    dependent_sections = ["laminates"]

    def _execute(self):
        prefix = self.workdir / self.config["general"]["prefix"]
        mass_table = drape_summary.drape_summary(f"{prefix}_joined.vtu")
        mass_table.to_csv(f"{prefix}_mass.csv")
        mass_table.replace(to_replace="_", value="", regex=True).to_latex(
            f"{prefix}_mass.tex", index=False
        )


class ApplyLoadsStep(Statesman):
    """Step for applying loads to mesh."""

    input_files = [
        ManagedFile(name="drape_output.vtu", non_empty=True),
    ]
    output_files = ["loads_output.png"]
    dependent_sections = ["loads"]

    def _execute(self):
        prefix = self.workdir / self.config["general"]["prefix"]
        add_load_to_mesh.add_load_to_mesh(
            self.config,
            f"{prefix}_joined.vtu",
            f"{prefix}_loads.png",
        )


# CLI remains similar, but now uses statesman steps
# ... (rest of CLI code unchanged)


def run_callback(yml: Path, bondline: bool):
    state = AppState.get_instance()
    app = BuildApp(state, yml)
    app.build(bondline=bondline)


def geometry_callback(yml: Path):
    state = AppState.get_instance()
    app = BuildApp(state, yml)
    app.geometry()


def mesh_callback(yml: Path):
    state = AppState.get_instance()
    app = BuildApp(state, yml)
    app.mesh()


def drape_callback(yml: Path, bondline: bool):
    state = AppState.get_instance()
    app = BuildApp(state, yml)
    app.drape(bondline=bondline)


def mass_callback(yml: Path):
    state = AppState.get_instance()
    app = BuildApp(state, yml)
    app.mass()


def apply_loads_callback(yml: Path):
    state = AppState.get_instance()
    app = BuildApp(state, yml)
    app.apply_loads()


build_cli = cli(
    name="build",
    help="Build the blade model",
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

geometry_cmd = command(
    name="geometry",
    help="Build blade geometry",
    callback=geometry_callback,
    arguments=[],
)
geometry_cmd.sort_key = 10
build_cli.commands.append(geometry_cmd)

mesh_cmd = command(
    name="mesh",
    help="Mesh blade structure",
    callback=mesh_callback,
    arguments=[],
)
mesh_cmd.sort_key = 20
build_cli.commands.append(mesh_cmd)

drape_cmd = command(
    name="drape",
    help="Drape plies onto mesh",
    callback=drape_callback,
    arguments=[],
    options=[
        option(
            flags=["--bondline", "-b"],
            is_flag=True,
            default=True,
            help="Add bondline to mesh",
            arg_type=bool,
        ),
    ],
)
drape_cmd.sort_key = 30
build_cli.commands.append(drape_cmd)

apply_loads_cmd = command(
    name="apply-loads",
    help="Apply loads to mesh",
    callback=apply_loads_callback,
    arguments=[],
)
apply_loads_cmd.sort_key = 40
build_cli.commands.append(apply_loads_cmd)

mass_cmd = command(
    name="mass",
    help="Calculate blade mass",
    callback=mass_callback,
    arguments=[],
)
mass_cmd.sort_key = 50
build_cli.commands.append(mass_cmd)

run_cmd = command(
    name="run",
    help="Build the full blade model",
    callback=run_callback,
    arguments=[],
    options=[
        option(
            flags=["--bondline", "-b"],
            is_flag=True,
            default=True,
            help="Include bondline",
            arg_type=bool,
        ),
    ],
)
run_cmd.sort_key = 60
build_cli.commands.append(run_cmd)
