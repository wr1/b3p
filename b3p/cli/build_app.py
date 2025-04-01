from pathlib import Path
import os
import pickle
from b3p import (
    build_plybook,
    drape_mesh,
    combine_meshes,
    add_load_to_mesh,
    drape_summary,
    add_te_solids,
    yml_portable,
)
from b3p.core import build_blade_geometry, build_blade_structure


class BuildApp:
    def __init__(self, state):
        self.state = state

    def geometry(self, yml: Path):
        dct = self.state.load_yaml(yml)
        build_blade_geometry.build_blade_geometry(dct)
        prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
        yml_portable.save_yaml(f"{prefix}_portable.yml", dct)

    def mesh(self, yml: Path):
        dct = self.state.load_yaml(yml)
        build_blade_structure.build_blade_structure(dct)

    def drape(self, yml: Path, bondline: bool = False):
        dct = self.state.load_yaml(yml)
        plybookname = "__plybook.pck"
        prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
        build_plybook.lamplan2plies(dct, plybookname)
        slb = dct["laminates"]["slabs"]
        used_grids = {slb[i]["grid"] for i in slb}
        pbook = os.path.join(dct["general"]["workdir"], plybookname)
        if os.path.exists(pbook):
            plybook = pickle.load(open(pbook, "rb"))
            meshes = []
            for grid in used_grids:
                out = f"{prefix}_{grid}_dr.vtu"
                drape_mesh.drape_mesh(f"{prefix}_{grid}.vtp", plybook, grid, out)
                meshes.append(out)
            combine_meshes.combine_meshes(meshes, f"{prefix}_joined.vtu")
            if bondline:
                add_te_solids.add_bondline(dct)
        else:
            raise FileNotFoundError("Plybook not found")

    def mass(self, yml: Path):
        dct = self.state.load_yaml(yml)
        wd = Path(dct["general"]["workdir"])
        prefix = self.state.get_prefix()
        if not wd.is_dir():
            print("\n** No workdir found, building new\n")
            self.build(yml)
        mass_table = drape_summary.drape_summary(f"{prefix}_joined.vtu")
        mass_table.to_csv(f"{prefix}_mass.csv")
        mass_table.replace(to_replace="_", value="", regex=True).to_latex(
            f"{prefix}_mass.tex", index=False
        )
        print("Mass table per material")
        print(mass_table)

    def apply_loads(self, yml: Path):
        dct = self.state.load_yaml(yml)
        add_load_to_mesh.add_load_to_mesh(
            dct,
            f"{self.state.get_prefix()}_joined.vtu",
            f"{self.state.get_prefix()}_loads.png",
        )

    def build(self, yml: Path, bondline: bool = True):
        self.geometry(yml)
        self.mesh(yml)
        self.drape(yml, bondline=bondline)
        self.mass(yml)
        self.apply_loads(yml)
