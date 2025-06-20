import logging
from pathlib import Path
import os
import pickle
from . import yml_portable
from b3p.geometry import build_blade_geometry
from b3p.mesh import (
    add_load_to_mesh,
    add_te_solids,
    build_blade_structure,
    combine_meshes,
)
from b3p.laminates import build_plybook, drape_mesh, drape_summary
from rich.logging import RichHandler  # Add rich log formatting


logging.basicConfig(handlers=[RichHandler(rich_tracebacks=True)], level=logging.INFO)
# Configure rich for logging
logger = logging.getLogger(__name__)


class BuildApp:
    def __init__(self, state, yml: Path):
        self.state = state
        self.dct = self.state.load_yaml(yml)

    def geometry(self):
        prefix = self.state.get_prefix("mesh")
        build_blade_geometry.build_blade_geometry(self.dct, prefix)
        prefix = self.state.get_prefix()
        yml_portable.save_yaml(f"{prefix}_portable.yml", self.dct)

    def mesh(self):
        prefix = self.state.get_prefix("mesh")
        build_blade_structure.build_blade_structure(self.dct, prefix)

    def drape(self, bondline: bool = True):
        plybookname = "_plybook.pck"
        workdir = self.state.get_workdir()
        prefix = self.state.get_prefix("drape")
        mesh_prefix = self.state.get_prefix("mesh")
        pbookpath = str(prefix) + plybookname

        logger.info(f"prefix and mesh prefix: {prefix}, {mesh_prefix}")

        build_plybook.lamplan2plies(self.dct, pbookpath)
        slb = self.dct["laminates"]["slabs"]
        used_grids = {slb[i]["grid"] for i in slb}

        if os.path.exists(pbookpath):
            plybook = pickle.load(open(pbookpath, "rb"))
            meshes = []
            for grid in used_grids:
                out = f"{prefix}_{grid}_dr.vtu"
                logger.info(f"Draping mesh for grid {grid} to {out}")
                drape_mesh.drape_mesh(f"{mesh_prefix}_{grid}.vtp", plybook, grid, out)
                meshes.append(out)
            grid = f"{prefix}_joined.vtu"
            combine_meshes.combine_meshes(meshes, grid)

            logger.info(f" running with bondline: {bondline}")
            if bondline:
                add_te_solids.add_bondline(
                    grid,
                    workdir / "drape" / "material_map.json",
                    bondline_config=self.dct["mesh"]["bondline"],
                )
                logger.info("Bondline added to mesh")
        else:
            logger.error(f"Plybook not found at {pbookpath}")
            raise FileNotFoundError(f"Plybook not found at {pbookpath}")

    def mass(self):
        wd = self.state.get_workdir()
        prefix = self.state.get_prefix("drape")
        if not wd.is_dir():
            logger.warning("No workdir found, building new")
            self.build(bondline=True)

        mass_table = drape_summary.drape_summary(f"{prefix}_joined.vtu")
        mass_table.to_csv(f"{prefix}_mass.csv")
        mass_table.replace(to_replace="_", value="", regex=True).to_latex(
            f"{prefix}_mass.tex", index=False
        )
        logger.info("Mass table per material:\n%s", mass_table)

    def apply_loads(self):
        prefix = self.state.get_prefix("drape")
        add_load_to_mesh.add_load_to_mesh(
            self.dct,
            f"{prefix}_joined.vtu",
            f"{prefix}_loads.png",
        )

    def build(self, bondline: bool = True):
        self.geometry()
        self.mesh()
        self.drape(bondline=bondline)
        self.mass()
        self.apply_loads()
